//! Mesh canvas drawing and pointer interaction for the native GUI.
//!
//! The canvas works in screen coordinates, but all geometry starts in the
//! `Model2D` coordinate system. `CanvasProjector` performs that mapping and
//! keeps pan/zoom handling independent from the model itself.

use super::loaded_model::{LoadedModel, ModelBounds};
use super::results::{AnalysisResults, ResultField, ScalarRange};
use super::selection::{SelectedEntity, ViewOptions};
use super::theme::{ACCENT, BORDER, CANVAS_BACKGROUND, ERROR, MUTED_TEXT, TEXT};
use super::topology::{element_edge_segments, element_outline_node_ids, normalized_edge};
use super::workflow::{DrawTool, ElementDraftKind, WorkMode};
use eframe::egui::{self, Color32, FontId, Pos2, Rect, Sense, Stroke, Vec2};
use rusty_fem::elements::Element2D;
use rusty_fem::model::{Dof2D, ElementLoad2D, Model2D};
use rusty_fem::visualisation::scalar_color_rgb;
use std::collections::{HashMap, HashSet};

/// Current user-controlled canvas transform.
#[derive(Debug, Clone, Copy)]
pub(super) struct MeshViewState {
    zoom: f32,
    pan: Vec2,
    dragging_node: Option<usize>,
    drag_target: Option<[f64; 2]>,
}

impl MeshViewState {
    pub(super) fn fit(&mut self) {
        self.zoom = 1.0;
        self.pan = Vec2::ZERO;
        self.dragging_node = None;
        self.drag_target = None;
    }

    fn zoom_at(&mut self, factor: f32, anchor: Pos2, canvas_center: Pos2) {
        if !factor.is_finite() || factor <= 0.0 {
            return;
        }

        let previous_zoom = self.zoom;
        let next_zoom = (previous_zoom * factor).clamp(0.15, 20.0);
        let applied_factor = next_zoom / previous_zoom;

        if (applied_factor - 1.0).abs() <= f32::EPSILON {
            return;
        }

        let previous_center = canvas_center + self.pan;
        let next_center = anchor - (anchor - previous_center) * applied_factor;
        self.pan = next_center - canvas_center;
        self.zoom = next_zoom;
    }
}

impl Default for MeshViewState {
    fn default() -> Self {
        Self { zoom: 1.0, pan: Vec2::ZERO, dragging_node: None, drag_target: None }
    }
}

/// Canvas editing configuration owned by the application.
pub(super) struct CanvasEditor<'a> {
    pub(super) work_mode: WorkMode,
    pub(super) tool: DrawTool,
    pub(super) show_grid: bool,
    pub(super) snap_to_grid: bool,
    pub(super) grid_spacing: f64,
    pub(super) element_kind: ElementDraftKind,
    pub(super) element_draft_points: &'a [[f64; 2]],
}

/// Read-only view configuration passed to the canvas for one frame.
pub(super) struct CanvasView<'a> {
    pub(super) options: &'a ViewOptions,
    pub(super) editor: CanvasEditor<'a>,
    pub(super) results: Option<&'a AnalysisResults>,
    pub(super) result_field: ResultField,
    pub(super) scalar_range: Option<ScalarRange>,
    pub(super) deformation_scale: f64,
}

/// One model-changing gesture emitted by the canvas.
pub(super) enum CanvasAction {
    InsertNode { x: f64, y: f64 },
    PlaceElementPoint { x: f64, y: f64 },
    MoveNode { node_id: usize, x: f64, y: f64 },
    AddConstraint(usize),
    AddNodalLoad(usize),
    AddEdgeTraction { element_id: usize, node_ids: [usize; 2] },
}

/// Information returned after painting one canvas frame.
pub(super) struct CanvasOutput {
    pub(super) action: Option<CanvasAction>,
    pub(super) cursor_coordinates: Option<[f64; 2]>,
}

/// Draws the main model canvas and updates selection when the user clicks.
pub(super) fn draw_model_canvas(
    ui: &mut egui::Ui, loaded_model: Option<&LoadedModel>, view_state: &mut MeshViewState,
    selection: &mut Option<SelectedEntity>, view: &CanvasView<'_>,
) -> CanvasOutput {
    let available = ui.available_size().max(Vec2::new(200.0, 200.0));
    let (response, painter) = ui.allocate_painter(available, Sense::click_and_drag());
    let rect = response.rect;
    draw_canvas_background(&painter, rect);

    let Some(model) = loaded_model else {
        draw_canvas_message(&painter, rect, "No model loaded");
        return CanvasOutput { action: None, cursor_coordinates: None };
    };

    let bounds = model.bounds.unwrap_or_else(ModelBounds::drawing_default);
    let is_editing = view.editor.work_mode == WorkMode::DrawFem;
    let is_moving_node = is_editing && view.editor.tool == DrawTool::MoveNode;

    if response.dragged_by(egui::PointerButton::Middle)
        || (!is_editing && response.dragged_by(egui::PointerButton::Primary))
    {
        view_state.pan += response.drag_delta();
    }

    if response.hovered() {
        let (pinch_factor, scroll_delta) = ui.input(|input| (input.zoom_delta(), input.smooth_scroll_delta().y));
        let scroll_factor =
            if scroll_delta.abs() > f32::EPSILON { (1.0 + scroll_delta * 0.001).clamp(0.8, 1.2) } else { 1.0 };
        let zoom_factor = if (pinch_factor - 1.0).abs() > f32::EPSILON { pinch_factor } else { scroll_factor };

        if (zoom_factor - 1.0).abs() > f32::EPSILON {
            let anchor = response.hover_pos().unwrap_or(rect.center());
            view_state.zoom_at(zoom_factor, anchor, rect.shrink(28.0).center());
        }
    }

    let projector = CanvasProjector::new(bounds, rect.shrink(28.0), view_state.zoom, view_state.pan);
    let cursor_coordinates = response.hover_pos().map(|position| {
        snap_coordinates(projector.unproject(position), view.editor.snap_to_grid, view.editor.grid_spacing)
    });
    let node_lookup = node_lookup(&model.model);
    let boundary_edges = boundary_edge_set(&model.model);
    let loaded_edges = loaded_edge_set(&model.model);
    let constrained_nodes = constrained_node_set(&model.model);
    let picked =
        response.hover_pos().and_then(|hover_pos| picked_entity(&model.model, hover_pos, &node_lookup, &projector));

    if view.editor.show_grid {
        draw_grid(&painter, rect, &projector, view.editor.grid_spacing);
    }

    let mut action = None;

    if is_moving_node {
        handle_node_drag(&response, &model.model, &projector, &view.editor, view_state, selection, &mut action);
    } else if response.clicked_by(egui::PointerButton::Primary) && !response.dragged() {
        action = handle_canvas_click(&model.model, &view.editor, cursor_coordinates, picked.as_ref(), selection);
    }

    if view.result_field.is_scalar_contour() {
        draw_scalar_fill(&painter, &projector, view.results, view.result_field, view.scalar_range);
    }

    if view.options.show_mesh {
        draw_model_edges(&painter, &model.model, &node_lookup, &projector, &HashSet::new(), mesh_stroke());
    }

    if view.options.show_boundary_edges {
        draw_model_edges(&painter, &model.model, &node_lookup, &projector, &boundary_edges, boundary_stroke());
    }

    if view.options.show_loads {
        draw_model_edges(&painter, &model.model, &node_lookup, &projector, &loaded_edges, load_stroke());
    }

    draw_selected_entity(&painter, selection.as_ref(), &model.model, &node_lookup, &projector);

    if view.result_field == ResultField::Displacement {
        draw_deformed_mesh(&painter, &model.model, &projector, view.results, view.deformation_scale);
    }

    if view.editor.tool == DrawTool::InsertElement {
        draw_element_preview(
            &painter,
            view.editor.element_kind,
            view.editor.element_draft_points,
            cursor_coordinates,
            &projector,
        );
    }

    if view.options.show_element_ids {
        draw_element_labels(&painter, &model.model, &node_lookup, &projector);
    }

    draw_nodes(&painter, &model.model, &node_lookup, &projector, &constrained_nodes, view.options, selection.as_ref());

    if view.options.show_loads {
        draw_nodal_loads(&painter, &model.model, &node_lookup, &projector);
    }

    if view.options.show_node_ids {
        draw_node_labels(&painter, &model.model, &projector);
    }

    if let Some(picked) = picked {
        draw_hover_label(&painter, rect, &picked.label);
    } else if let Some([x, y]) = cursor_coordinates {
        draw_coordinate_label(&painter, rect, x, y);
    }

    if let Some(results) = view.results {
        draw_result_legend(&painter, rect, view.result_field, results, view.scalar_range, view.deformation_scale);
    }

    CanvasOutput { action, cursor_coordinates }
}

fn draw_canvas_background(painter: &egui::Painter, rect: Rect) {
    painter.rect_filled(rect, 6.0, CANVAS_BACKGROUND);
    painter.rect_stroke(rect, 6.0, Stroke::new(1.0, BORDER), egui::StrokeKind::Inside);
}

fn draw_canvas_message(painter: &egui::Painter, rect: Rect, text: &str) {
    painter.text(rect.center(), egui::Align2::CENTER_CENTER, text, FontId::proportional(18.0), MUTED_TEXT);
}

struct PickedEntity {
    entity: SelectedEntity,
    label: String,
}

struct CanvasProjector {
    bounds: ModelBounds,
    scale: f32,
    center: Pos2,
}

impl CanvasProjector {
    fn new(bounds: ModelBounds, rect: Rect, zoom: f32, pan: Vec2) -> Self {
        let mut width = (bounds.max_x - bounds.min_x).abs() as f32;
        let mut height = (bounds.max_y - bounds.min_y).abs() as f32;

        if width <= f32::EPSILON && height <= f32::EPSILON {
            width = 20.0;
            height = 20.0;
        } else {
            width = width.max(1e-9);
            height = height.max(1e-9);
        }
        let scale = (rect.width() / width).min(rect.height() / height) * zoom;
        let center = rect.center() + pan;

        Self { bounds, scale, center }
    }

    fn project(&self, x: f64, y: f64) -> Pos2 {
        let model_center_x = 0.5 * (self.bounds.min_x + self.bounds.max_x);
        let model_center_y = 0.5 * (self.bounds.min_y + self.bounds.max_y);
        let screen_x = self.center.x + ((x - model_center_x) as f32 * self.scale);
        let screen_y = self.center.y - ((y - model_center_y) as f32 * self.scale);

        Pos2::new(screen_x, screen_y)
    }

    fn unproject(&self, position: Pos2) -> [f64; 2] {
        let model_center_x = 0.5 * (self.bounds.min_x + self.bounds.max_x);
        let model_center_y = 0.5 * (self.bounds.min_y + self.bounds.max_y);

        [
            model_center_x + ((position.x - self.center.x) / self.scale) as f64,
            model_center_y - ((position.y - self.center.y) / self.scale) as f64,
        ]
    }
}

fn snap_coordinates(coordinates: [f64; 2], enabled: bool, spacing: f64) -> [f64; 2] {
    if !enabled || !spacing.is_finite() || spacing <= 0.0 {
        return coordinates;
    }

    [(coordinates[0] / spacing).round() * spacing, (coordinates[1] / spacing).round() * spacing]
}

fn draw_grid(painter: &egui::Painter, rect: Rect, projector: &CanvasProjector, requested_spacing: f64) {
    if !requested_spacing.is_finite() || requested_spacing <= 0.0 {
        return;
    }

    let [left, top] = projector.unproject(rect.left_top());
    let [right, bottom] = projector.unproject(rect.right_bottom());
    let width = (right - left).abs();
    let height = (top - bottom).abs();
    let mut spacing = requested_spacing;

    while width / spacing > 120.0 || height / spacing > 120.0 {
        spacing *= 2.0;
    }

    let stroke = Stroke::new(0.6, Color32::from_rgb(226, 232, 240));
    let first_x = (left.min(right) / spacing).floor() as i64;
    let last_x = (left.max(right) / spacing).ceil() as i64;
    let first_y = (bottom.min(top) / spacing).floor() as i64;
    let last_y = (bottom.max(top) / spacing).ceil() as i64;

    for index in first_x..=last_x {
        let x = index as f64 * spacing;
        painter.line_segment([projector.project(x, bottom), projector.project(x, top)], stroke);
    }

    for index in first_y..=last_y {
        let y = index as f64 * spacing;
        painter.line_segment([projector.project(left, y), projector.project(right, y)], stroke);
    }
}

fn handle_canvas_click(
    model: &Model2D, editor: &CanvasEditor<'_>, cursor_coordinates: Option<[f64; 2]>, picked: Option<&PickedEntity>,
    selection: &mut Option<SelectedEntity>,
) -> Option<CanvasAction> {
    match editor.tool {
        DrawTool::Select => {
            *selection = picked.map(|picked| picked.entity);
            None
        }
        DrawTool::InsertNode => {
            let [x, y] = cursor_coordinates?;
            Some(CanvasAction::InsertNode { x, y })
        }
        DrawTool::InsertElement => {
            let mut coordinates = cursor_coordinates?;

            if let Some(PickedEntity { entity: SelectedEntity::Node(node_id), .. }) = picked
                && let Some(node) = model.nodes().iter().find(|node| node.id() == *node_id)
            {
                coordinates = [node.x(), node.y()];
            }

            Some(CanvasAction::PlaceElementPoint { x: coordinates[0], y: coordinates[1] })
        }
        DrawTool::Constraint => match picked?.entity {
            SelectedEntity::Node(node_id) => {
                *selection = Some(SelectedEntity::Node(node_id));
                Some(CanvasAction::AddConstraint(node_id))
            }
            _ => None,
        },
        DrawTool::NodalLoad => match picked?.entity {
            SelectedEntity::Node(node_id) => {
                *selection = Some(SelectedEntity::Node(node_id));
                Some(CanvasAction::AddNodalLoad(node_id))
            }
            _ => None,
        },
        DrawTool::EdgeTraction => match picked?.entity {
            SelectedEntity::Edge { element_id, node_ids } => {
                Some(CanvasAction::AddEdgeTraction { element_id, node_ids })
            }
            _ => None,
        },
        DrawTool::MoveNode => None,
    }
}

fn handle_node_drag(
    response: &egui::Response, model: &Model2D, projector: &CanvasProjector, editor: &CanvasEditor<'_>,
    view_state: &mut MeshViewState, selection: &mut Option<SelectedEntity>, action: &mut Option<CanvasAction>,
) {
    if response.drag_started_by(egui::PointerButton::Primary)
        && let Some(position) = response.interact_pointer_pos()
        && let Some(picked) = hovered_node(model, position, projector)
        && let SelectedEntity::Node(node_id) = picked.entity
    {
        view_state.dragging_node = Some(node_id);
        *selection = Some(SelectedEntity::Node(node_id));
    }

    if view_state.dragging_node.is_some()
        && let Some(position) = response.interact_pointer_pos()
    {
        view_state.drag_target =
            Some(snap_coordinates(projector.unproject(position), editor.snap_to_grid, editor.grid_spacing));
    }

    if response.drag_stopped_by(egui::PointerButton::Primary) {
        if let (Some(node_id), Some([x, y])) = (view_state.dragging_node.take(), view_state.drag_target.take()) {
            *action = Some(CanvasAction::MoveNode { node_id, x, y });
        }
    } else if let (Some(_), Some([x, y])) = (view_state.dragging_node, view_state.drag_target) {
        let position = projector.project(x, y);
        response.ctx.set_cursor_icon(egui::CursorIcon::Grabbing);
        response.ctx.debug_painter().circle_filled(position, 6.0, ACCENT);
    }
}

fn draw_element_preview(
    painter: &egui::Painter, kind: ElementDraftKind, draft_points: &[[f64; 2]], cursor_coordinates: Option<[f64; 2]>,
    projector: &CanvasProjector,
) {
    let mut positions = draft_points.iter().map(|&[x, y]| projector.project(x, y)).collect::<Vec<_>>();

    if let Some([x, y]) = cursor_coordinates {
        positions.push(projector.project(x, y));
    }

    for points in positions.windows(2) {
        painter.line_segment([points[0], points[1]], Stroke::new(2.0, ACCENT));
    }

    for position in positions.iter().take(positions.len().saturating_sub(1)) {
        painter.circle_filled(*position, 5.0, ACCENT);
    }

    if kind.is_surface()
        && positions.len() == kind.placement_point_count()
        && let (Some(first), Some(last)) = (positions.first(), positions.last())
    {
        painter.line_segment([*last, *first], Stroke::new(2.0, ACCENT));
    }
}

fn draw_coordinate_label(painter: &egui::Painter, canvas_rect: Rect, x: f64, y: f64) {
    painter.text(
        canvas_rect.left_bottom() + Vec2::new(14.0, -14.0),
        egui::Align2::LEFT_BOTTOM,
        format!("x={x:.4}  y={y:.4}"),
        FontId::monospace(12.0),
        MUTED_TEXT,
    );
}

fn node_lookup(model: &Model2D) -> HashMap<usize, (f64, f64)> {
    model.nodes().iter().map(|node| (node.id(), (node.x(), node.y()))).collect()
}

fn constrained_node_set(model: &Model2D) -> HashSet<usize> {
    model.constraints().iter().map(|constraint| constraint.node_id()).collect()
}

fn loaded_edge_set(model: &Model2D) -> HashSet<[usize; 2]> {
    let mut edges = HashSet::new();

    for load in model.element_loads() {
        if let ElementLoad2D::EdgeTraction(load) = load {
            let Some(element) = model.elements().iter().find(|element| element.id() == load.element_id()) else {
                continue;
            };

            for edge in element_edge_segments(element, load.edge_node_ids()) {
                edges.insert(edge);
            }
        }
    }

    edges
}

fn boundary_edge_set(model: &Model2D) -> HashSet<[usize; 2]> {
    let mut edge_counts: HashMap<[usize; 2], usize> = HashMap::new();

    for element in model.elements() {
        for edge in element_outline_node_ids(element).windows(2) {
            *edge_counts.entry(normalized_edge(edge[0], edge[1])).or_default() += 1;
        }
    }

    edge_counts.into_iter().filter_map(|(edge, count)| (count == 1).then_some(edge)).collect()
}

fn draw_model_edges(
    painter: &egui::Painter, model: &Model2D, node_lookup: &HashMap<usize, (f64, f64)>, projector: &CanvasProjector,
    edges_to_draw: &HashSet<[usize; 2]>, stroke: Stroke,
) {
    for element in model.elements() {
        for edge in element_outline_node_ids(element).windows(2) {
            if !edges_to_draw.is_empty() && !edges_to_draw.contains(&normalized_edge(edge[0], edge[1])) {
                continue;
            }

            draw_edge_segment(painter, edge[0], edge[1], node_lookup, projector, stroke);
        }
    }
}

fn draw_scalar_fill(
    painter: &egui::Painter, projector: &CanvasProjector, results: Option<&AnalysisResults>, field: ResultField,
    range: Option<ScalarRange>,
) {
    let Some(results) = results else {
        return;
    };
    let Some(range) = range else {
        return;
    };

    let mut mesh = egui::Mesh::default();
    mesh.reserve_triangles(results.contour_triangles().len());

    for triangle in results.contour_triangles() {
        let first_index = mesh.vertices.len() as u32;
        let mut complete = true;

        for vertex in triangle.vertices {
            let Some(value) = vertex.state.scalar(field) else {
                complete = false;
                break;
            };
            mesh.colored_vertex(
                projector.project(vertex.position[0], vertex.position[1]),
                contour_color(range.normalized(value)),
            );
        }

        if complete {
            mesh.add_triangle(first_index, first_index + 1, first_index + 2);
        }
    }

    painter.add(egui::Shape::mesh(mesh));
}

fn draw_deformed_mesh(
    painter: &egui::Painter, model: &Model2D, projector: &CanvasProjector, results: Option<&AnalysisResults>,
    deformation_scale: f64,
) {
    let Some(results) = results else {
        return;
    };

    for element in model.elements() {
        let positions = element_outline_node_ids(element)
            .into_iter()
            .filter_map(|node_id| {
                let node = model.nodes().iter().find(|node| node.id() == node_id)?;
                let [ux, uy] = results.displacement(node_id);
                Some(projector.project(node.x() + deformation_scale * ux, node.y() + deformation_scale * uy))
            })
            .collect::<Vec<_>>();

        for points in positions.windows(2) {
            painter.line_segment([points[0], points[1]], Stroke::new(2.0, ACCENT));
        }
    }
}

fn draw_result_legend(
    painter: &egui::Painter, canvas_rect: Rect, result_field: ResultField, results: &AnalysisResults,
    scalar_range: Option<ScalarRange>, deformation_scale: f64,
) {
    if result_field == ResultField::Model {
        return;
    }

    if result_field == ResultField::Displacement {
        let maximum = format_legend_value(results.max_displacement, results.max_displacement.abs());
        let text = format!("Displacement  max {maximum}  scale {deformation_scale:.3}");
        let legend_rect =
            Rect::from_min_size(canvas_rect.right_top() + Vec2::new(-344.0, 14.0), Vec2::new(330.0, 38.0));

        draw_legend_background(painter, legend_rect);
        painter.text(legend_rect.center(), egui::Align2::CENTER_CENTER, text, FontId::proportional(13.0), TEXT);
        return;
    }

    let Some(range) = scalar_range else {
        return;
    };
    let legend_width = (canvas_rect.width() - 28.0).clamp(180.0, 420.0);
    let legend_rect = Rect::from_min_size(
        Pos2::new(canvas_rect.right() - legend_width - 14.0, canvas_rect.top() + 14.0),
        Vec2::new(legend_width, 82.0),
    );

    draw_legend_background(painter, legend_rect);
    painter.text(
        legend_rect.left_top() + Vec2::new(14.0, 8.0),
        egui::Align2::LEFT_TOP,
        result_field.label(),
        FontId::proportional(13.0),
        TEXT,
    );

    let bar = Rect::from_min_max(
        legend_rect.left_top() + Vec2::new(14.0, 29.0),
        legend_rect.right_top() + Vec2::new(-14.0, 43.0),
    );
    let segments = 96;

    for index in 0..segments {
        let first = index as f32 / segments as f32;
        let second = (index + 1) as f32 / segments as f32;
        let segment = Rect::from_min_max(
            Pos2::new(egui::lerp(bar.x_range(), first), bar.top()),
            Pos2::new(egui::lerp(bar.x_range(), second), bar.bottom()),
        );
        painter.rect_filled(segment, 0.0, contour_color(first));
    }

    painter.rect_stroke(bar, 1.0, Stroke::new(1.0, BORDER), egui::StrokeKind::Inside);

    let tick_fractions: &[f32] = if legend_width >= 340.0 { &[0.0, 0.25, 0.5, 0.75, 1.0] } else { &[0.0, 0.5, 1.0] };
    let magnitude = range.minimum.abs().max(range.maximum.abs());

    for (index, fraction) in tick_fractions.iter().copied().enumerate() {
        let x = egui::lerp(bar.x_range(), fraction);
        let value = egui::lerp(range.minimum..=range.maximum, fraction as f64);
        let formatted = format_legend_value(value, magnitude);
        let (alignment, label) = if index == 0 {
            (egui::Align2::LEFT_TOP, format!("min {formatted}"))
        } else if index + 1 == tick_fractions.len() {
            (egui::Align2::RIGHT_TOP, format!("max {formatted}"))
        } else {
            (egui::Align2::CENTER_TOP, formatted)
        };

        painter.line_segment([Pos2::new(x, bar.bottom()), Pos2::new(x, bar.bottom() + 5.0)], Stroke::new(1.0, BORDER));
        painter.text(Pos2::new(x, bar.bottom() + 8.0), alignment, label, FontId::monospace(11.0), TEXT);
    }
}

fn draw_legend_background(painter: &egui::Painter, legend_rect: Rect) {
    painter.rect_filled(legend_rect, 4.0, Color32::from_rgba_premultiplied(255, 255, 255, 235));
    painter.rect_stroke(legend_rect, 4.0, Stroke::new(1.0, BORDER), egui::StrokeKind::Inside);
}

fn format_legend_value(value: f64, magnitude: f64) -> String {
    if value.abs() <= f64::EPSILON {
        return "0".to_owned();
    }

    if !(0.01..10_000.0).contains(&magnitude) {
        return format!("{value:.3e}");
    }

    let integer_digits = magnitude.max(1.0).log10().floor() as i32 + 1;
    let decimal_places = (4 - integer_digits).clamp(0, 5) as usize;

    format!("{value:.decimal_places$}")
}

fn contour_color(value: f32) -> Color32 {
    let [red, green, blue] = scalar_color_rgb(value as f64);
    Color32::from_rgb(red, green, blue)
}

fn draw_selected_entity(
    painter: &egui::Painter, selection: Option<&SelectedEntity>, model: &Model2D,
    node_lookup: &HashMap<usize, (f64, f64)>, projector: &CanvasProjector,
) {
    let Some(selection) = selection else {
        return;
    };

    match selection {
        SelectedEntity::Node(node_id) => {
            if let Some(&(x, y)) = node_lookup.get(node_id) {
                painter.circle_filled(projector.project(x, y), 7.0, ACCENT);
            }
        }
        SelectedEntity::Element(element_id) => {
            if let Some(element) = model.elements().iter().find(|element| element.id() == *element_id) {
                for edge in element_outline_node_ids(element).windows(2) {
                    draw_edge_segment(painter, edge[0], edge[1], node_lookup, projector, selected_stroke());
                }
            }
        }
        SelectedEntity::Edge { node_ids, .. } => {
            draw_edge_segment(painter, node_ids[0], node_ids[1], node_lookup, projector, selected_stroke());
        }
    }
}

fn mesh_stroke() -> Stroke {
    Stroke::new(0.8, Color32::from_rgb(148, 163, 184))
}

fn boundary_stroke() -> Stroke {
    Stroke::new(1.8, TEXT)
}

fn load_stroke() -> Stroke {
    Stroke::new(2.4, ERROR)
}

fn selected_stroke() -> Stroke {
    Stroke::new(3.2, ACCENT)
}

fn draw_edge_segment(
    painter: &egui::Painter, first: usize, second: usize, node_lookup: &HashMap<usize, (f64, f64)>,
    projector: &CanvasProjector, stroke: Stroke,
) {
    let Some(&(first_x, first_y)) = node_lookup.get(&first) else {
        return;
    };
    let Some(&(second_x, second_y)) = node_lookup.get(&second) else {
        return;
    };

    painter.line_segment([projector.project(first_x, first_y), projector.project(second_x, second_y)], stroke);
}

fn draw_nodes(
    painter: &egui::Painter, model: &Model2D, node_lookup: &HashMap<usize, (f64, f64)>, projector: &CanvasProjector,
    constrained_nodes: &HashSet<usize>, view_options: &ViewOptions, selection: Option<&SelectedEntity>,
) {
    for node in model.nodes() {
        let is_constrained = constrained_nodes.contains(&node.id());
        let is_selected = matches!(selection, Some(SelectedEntity::Node(selected_node)) if *selected_node == node.id());

        if !view_options.show_nodes && !(view_options.show_constraints && is_constrained) && !is_selected {
            continue;
        }

        let Some(&(x, y)) = node_lookup.get(&node.id()) else {
            continue;
        };
        let position = projector.project(x, y);

        if is_selected {
            painter.circle_filled(position, 6.8, ACCENT);
            painter.circle_stroke(position, 8.2, Stroke::new(1.0, Color32::WHITE));
        } else if view_options.show_constraints && is_constrained {
            painter.circle_filled(position, 4.2, Color32::from_rgb(51, 65, 85));
            painter.circle_stroke(position, 5.4, Stroke::new(1.0, Color32::WHITE));
        } else {
            painter.circle_filled(position, 2.1, Color32::from_rgb(71, 85, 105));
        }
    }
}

fn draw_nodal_loads(
    painter: &egui::Painter, model: &Model2D, node_lookup: &HashMap<usize, (f64, f64)>, projector: &CanvasProjector,
) {
    for load in model.loads() {
        if let Some(&(x, y)) = node_lookup.get(&load.node_id()) {
            let start = projector.project(x, y);
            let direction = load_direction(load.dof(), load.value());
            draw_arrow(painter, start, direction, ERROR);
        }
    }
}

fn draw_node_labels(painter: &egui::Painter, model: &Model2D, projector: &CanvasProjector) {
    for node in model.nodes() {
        let position = projector.project(node.x(), node.y()) + Vec2::new(4.0, -4.0);
        painter.text(position, egui::Align2::LEFT_BOTTOM, node.id().to_string(), FontId::monospace(10.0), TEXT);
    }
}

fn draw_element_labels(
    painter: &egui::Painter, model: &Model2D, node_lookup: &HashMap<usize, (f64, f64)>, projector: &CanvasProjector,
) {
    for element in model.elements() {
        if let Some(position) = element_centroid(element, node_lookup, projector) {
            painter.text(
                position,
                egui::Align2::CENTER_CENTER,
                element.id().to_string(),
                FontId::monospace(10.0),
                TEXT,
            );
        }
    }
}

fn draw_hover_label(painter: &egui::Painter, canvas_rect: Rect, text: &str) {
    painter.rect_filled(
        Rect::from_min_size(canvas_rect.left_top() + Vec2::new(16.0, 16.0), Vec2::new(300.0, 32.0)),
        6.0,
        Color32::from_rgba_premultiplied(15, 23, 42, 215),
    );
    painter.text(
        canvas_rect.left_top() + Vec2::new(28.0, 32.0),
        egui::Align2::LEFT_CENTER,
        text,
        FontId::proportional(14.0),
        Color32::WHITE,
    );
}

fn picked_entity(
    model: &Model2D, hover_pos: Pos2, node_lookup: &HashMap<usize, (f64, f64)>, projector: &CanvasProjector,
) -> Option<PickedEntity> {
    hovered_node(model, hover_pos, projector)
        .or_else(|| hovered_edge(model, hover_pos, node_lookup, projector))
        .or_else(|| hovered_element(model, hover_pos, node_lookup, projector))
}

fn hovered_node(model: &Model2D, hover_pos: Pos2, projector: &CanvasProjector) -> Option<PickedEntity> {
    let nearest = model
        .nodes()
        .iter()
        .map(|node| (node, projector.project(node.x(), node.y()).distance(hover_pos)))
        .min_by(|(_, left_distance), (_, right_distance)| left_distance.total_cmp(right_distance))?;

    (nearest.1 <= 8.0).then(|| {
        let node = nearest.0;

        PickedEntity {
            entity: SelectedEntity::Node(node.id()),
            label: format!("node {}  x={:.4}, y={:.4}", node.id(), node.x(), node.y()),
        }
    })
}

fn hovered_edge(
    model: &Model2D, hover_pos: Pos2, node_lookup: &HashMap<usize, (f64, f64)>, projector: &CanvasProjector,
) -> Option<PickedEntity> {
    let nearest = model
        .elements()
        .iter()
        .flat_map(|element| {
            element_outline_node_ids(element)
                .windows(2)
                .map(|edge| (element.id(), [edge[0], edge[1]]))
                .collect::<Vec<_>>()
        })
        .filter_map(|(element_id, edge)| {
            let first = node_lookup.get(&edge[0])?;
            let second = node_lookup.get(&edge[1])?;
            let start = projector.project(first.0, first.1);
            let end = projector.project(second.0, second.1);

            Some((element_id, edge, distance_to_segment(hover_pos, start, end)))
        })
        .min_by(|(_, _, left_distance), (_, _, right_distance)| left_distance.total_cmp(right_distance))?;

    (nearest.2 <= 6.0).then(|| PickedEntity {
        entity: SelectedEntity::Edge { element_id: nearest.0, node_ids: nearest.1 },
        label: format!("edge {}-{}  element {}", nearest.1[0], nearest.1[1], nearest.0),
    })
}

fn hovered_element(
    model: &Model2D, hover_pos: Pos2, node_lookup: &HashMap<usize, (f64, f64)>, projector: &CanvasProjector,
) -> Option<PickedEntity> {
    let nearest = model
        .elements()
        .iter()
        .filter_map(|element| {
            let centroid = element_centroid(element, node_lookup, projector)?;
            Some((element, centroid.distance(hover_pos)))
        })
        .min_by(|(_, left_distance), (_, right_distance)| left_distance.total_cmp(right_distance))?;

    (nearest.1 <= 24.0).then(|| PickedEntity {
        entity: SelectedEntity::Element(nearest.0.id()),
        label: format!("element {}  {}", nearest.0.id(), nearest.0.element_type()),
    })
}

fn element_centroid(
    element: &Element2D, node_lookup: &HashMap<usize, (f64, f64)>, projector: &CanvasProjector,
) -> Option<Pos2> {
    let mut sum = Vec2::ZERO;
    let mut count = 0.0_f32;

    for &node_id in element.node_ids() {
        let &(x, y) = node_lookup.get(&node_id)?;
        sum += projector.project(x, y).to_vec2();
        count += 1.0;
    }

    Some(Pos2::new(sum.x / count, sum.y / count))
}

fn distance_to_segment(point: Pos2, start: Pos2, end: Pos2) -> f32 {
    let segment = end - start;
    let segment_length_squared = segment.length_sq();

    if segment_length_squared <= f32::EPSILON {
        return point.distance(start);
    }

    let t = ((point - start).dot(segment) / segment_length_squared).clamp(0.0, 1.0);
    let projection = start + segment * t;

    point.distance(projection)
}

fn load_direction(dof: Dof2D, value: f64) -> Vec2 {
    let sign = if value >= 0.0 { 1.0 } else { -1.0 };

    match dof {
        Dof2D::Ux => Vec2::new(sign * 22.0, 0.0),
        Dof2D::Uy => Vec2::new(0.0, -sign * 22.0),
        Dof2D::Rz => Vec2::new(sign * 16.0, -16.0),
    }
}

fn draw_arrow(painter: &egui::Painter, start: Pos2, direction: Vec2, color: Color32) {
    let end = start + direction;
    painter.line_segment([start, end], Stroke::new(2.0, color));

    if direction.length_sq() <= f32::EPSILON {
        return;
    }

    let unit = direction.normalized();
    let normal = Vec2::new(-unit.y, unit.x);
    let first = end - unit * 7.0 + normal * 4.0;
    let second = end - unit * 7.0 - normal * 4.0;

    painter.line_segment([end, first], Stroke::new(2.0, color));
    painter.line_segment([end, second], Stroke::new(2.0, color));
}

#[cfg(test)]
mod tests {
    use super::{MeshViewState, format_legend_value};
    use eframe::egui::Pos2;

    #[test]
    fn formats_legend_values_for_engineering_scales() {
        assert_eq!(format_legend_value(0.0, 44.511), "0");
        assert_eq!(format_legend_value(44.511, 44.511), "44.51");
        assert_eq!(format_legend_value(-3.25, 44.511), "-3.25");
        assert!(format_legend_value(0.000_125, 0.000_5).contains('e'));
    }

    #[test]
    fn zoom_keeps_the_model_point_under_the_pointer_stationary() {
        let mut state = MeshViewState::default();
        let canvas_center = Pos2::new(100.0, 100.0);
        let anchor = Pos2::new(150.0, 80.0);
        let before = (anchor - (canvas_center + state.pan)) / state.zoom;

        state.zoom_at(2.0, anchor, canvas_center);

        let after = (anchor - (canvas_center + state.pan)) / state.zoom;
        assert!((state.zoom - 2.0).abs() <= f32::EPSILON);
        assert!((before - after).length() <= f32::EPSILON);
    }
}
