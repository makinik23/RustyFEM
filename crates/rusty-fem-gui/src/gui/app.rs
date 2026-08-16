//! Top-level native GUI application state and layout.

use super::DEFAULT_MODEL_PATH;
use super::canvas::{CanvasAction, CanvasEditor, CanvasView, MeshViewState, draw_model_canvas};
use super::document::{EditHistory, write_model_json};
use super::loaded_model::LoadedModel;
use super::model_browser::discover_json_models;
use super::results::{AnalysisPhase, AnalysisProgress, AnalysisResults, ResultField, ScalarRange};
use super::selection::{SelectedEntity, ViewOptions};
use super::theme::{ACCENT, APP_BACKGROUND, BORDER, ERROR, MUTED_TEXT, PANEL_BACKGROUND, TEXT, WARNING};
use super::topology::{element_edge_segments, normalized_edge};
use super::workflow::{DrawTool, ElementDraftKind, WorkMode};
use eframe::egui::{self, Color32, Frame, Grid, Margin, RichText, Stroke};
use rusty_fem::elements::{Beam2D, Element2D, QuadQ4, QuadQ8, TriangleT3, TriangleT6, Truss2D};
use rusty_fem::io::{
    DisplacementConstraint2DInput, Dof2DInput, ElementLoad2DInput, ElementLoad2DInputKind, ElementType2DInput,
    LoadCoordinateSystem2DInput, Model2DInput, NodalLoad2DInput, Node2DInput,
};
use rusty_fem::model::{
    BeamSection2D, Dof2D, ElementLoad2D, Material2D, Model2D, PlaneStressSection2D, Section2D, SolverKind2D,
    TrussSection2D,
};
use std::path::{Path, PathBuf};
use std::sync::mpsc::{self, Receiver, TryRecvError};
use std::time::{Duration, Instant};

const SIDEBAR_WIDTH: f32 = 320.0;
const PROPERTIES_WIDTH: f32 = 330.0;
const EMPTY_CANVAS_MODEL_PATH: &str = "untitled-fem-model.json";
const DEFAULT_DRAWING_MATERIAL_ID: usize = 1;
const DEFAULT_DRAWING_SECTION_ID: usize = 1;
const DEFAULT_TRUSS_SECTION_ID: usize = 2;
const DEFAULT_BEAM_SECTION_ID: usize = 3;
const NODE_MERGE_TOLERANCE: f64 = 1e-9;
const SUCCESS: Color32 = Color32::from_rgb(21, 128, 61);

/// Main `eframe` app object.
///
/// It owns the loaded model and the view state, then delegates actual mesh
/// drawing to the canvas module. Keeping this type focused on application
/// layout makes it easier to add later panels such as selection properties,
/// solver controls, and result views.
pub(super) struct RustyFemGuiApp {
    loaded_model: Option<LoadedModel>,
    load_error: Option<String>,
    model_path_input: String,
    available_json_models: Vec<PathBuf>,
    selected_json_model: Option<usize>,
    view_state: MeshViewState,
    view_options: ViewOptions,
    selection: Option<SelectedEntity>,
    work_mode: WorkMode,
    draw_tool: DrawTool,
    element_kind: ElementDraftKind,
    element_material_id_input: String,
    element_section_id_input: String,
    element_draft_points: Vec<[f64; 2]>,
    show_grid: bool,
    snap_to_grid: bool,
    grid_spacing: f64,
    constraint_ux: bool,
    constraint_uy: bool,
    constraint_rz: bool,
    constraint_value_input: String,
    load_dof: Dof2D,
    load_value_input: String,
    traction_x_input: String,
    traction_y_input: String,
    history: EditHistory,
    saved_snapshot: Option<Model2DInput>,
    results: Option<AnalysisResults>,
    analysis_task: Option<AnalysisTask>,
    result_field: ResultField,
    contour_auto_range: bool,
    contour_manual_minimum: f64,
    contour_manual_maximum: f64,
    deformation_scale: f64,
    cursor_coordinates: Option<[f64; 2]>,
    model_dirty: bool,
    editor_message: Option<EditorMessage>,
}

struct EditorMessage {
    text: String,
    kind: EditorMessageKind,
}

struct AnalysisTask {
    receiver: Receiver<AnalysisTaskMessage>,
    model_snapshot: Model2DInput,
    phase: AnalysisPhase,
    iterations: usize,
    residual_norm: Option<f64>,
    relative_residual: Option<f64>,
    recovered_elements: usize,
    total_elements: usize,
    tolerance: f64,
    progress: f32,
    started: Instant,
}

enum AnalysisTaskMessage {
    Progress(AnalysisProgress),
    Finished(Box<Result<AnalysisResults, String>>),
}

impl AnalysisTask {
    fn new(receiver: Receiver<AnalysisTaskMessage>, model_snapshot: Model2DInput, tolerance: f64) -> Self {
        Self {
            receiver,
            model_snapshot,
            phase: AnalysisPhase::Assembling,
            iterations: 0,
            residual_norm: None,
            relative_residual: None,
            recovered_elements: 0,
            total_elements: 0,
            tolerance,
            progress: 0.05,
            started: Instant::now(),
        }
    }

    fn apply_progress(&mut self, progress: AnalysisProgress) {
        match progress {
            AnalysisProgress::Phase(phase) => {
                self.phase = phase;
                let phase_progress = match phase {
                    AnalysisPhase::Assembling => 0.05,
                    AnalysisPhase::Solving => 0.10,
                    AnalysisPhase::Recovering => 0.92,
                };
                self.progress = self.progress.max(phase_progress);
            }
            AnalysisProgress::Iteration(sample) => {
                self.phase = AnalysisPhase::Solving;
                self.iterations = sample.iterations;
                self.residual_norm = Some(sample.residual_norm);
                self.relative_residual = Some(sample.relative_residual_norm);
                let convergence = residual_convergence(sample.relative_residual_norm, self.tolerance);
                self.progress = self.progress.max(0.10 + 0.80 * convergence);
            }
            AnalysisProgress::Recovery { completed, total } => {
                self.phase = AnalysisPhase::Recovering;
                self.recovered_elements = completed;
                self.total_elements = total;
                let recovery_progress = if total == 0 { 1.0 } else { completed as f32 / total as f32 };
                self.progress = self.progress.max(0.92 + 0.08 * recovery_progress);
            }
        }
    }
}

fn residual_convergence(relative_residual: f64, tolerance: f64) -> f32 {
    if !relative_residual.is_finite() || !tolerance.is_finite() || tolerance <= 0.0 {
        return 0.0;
    }
    if relative_residual <= tolerance || tolerance >= 1.0 {
        return 1.0;
    }

    let target_orders = -tolerance.log10();
    let completed_orders = -relative_residual.min(1.0).log10();

    (completed_orders / target_orders).clamp(0.0, 1.0) as f32
}

enum EditorMessageKind {
    Info,
    Error,
}

impl RustyFemGuiApp {
    pub(super) fn new(input_path: Option<PathBuf>) -> Self {
        let model_path_input =
            input_path.as_ref().map(|path| path.display().to_string()).unwrap_or_else(|| DEFAULT_MODEL_PATH.to_owned());
        let mut app = Self {
            loaded_model: None,
            load_error: None,
            model_path_input,
            available_json_models: discover_json_models(),
            selected_json_model: None,
            view_state: MeshViewState::default(),
            view_options: ViewOptions::default(),
            selection: None,
            work_mode: WorkMode::Inspect,
            draw_tool: DrawTool::Select,
            element_kind: ElementDraftKind::T3,
            element_material_id_input: "1".to_owned(),
            element_section_id_input: "1".to_owned(),
            element_draft_points: Vec::new(),
            show_grid: true,
            snap_to_grid: true,
            grid_spacing: 1.0,
            constraint_ux: true,
            constraint_uy: true,
            constraint_rz: false,
            constraint_value_input: "0.0".to_owned(),
            load_dof: Dof2D::Uy,
            load_value_input: "-1.0".to_owned(),
            traction_x_input: "0.0".to_owned(),
            traction_y_input: "-1.0".to_owned(),
            history: EditHistory::default(),
            saved_snapshot: None,
            results: None,
            analysis_task: None,
            result_field: ResultField::Model,
            contour_auto_range: true,
            contour_manual_minimum: 0.0,
            contour_manual_maximum: 1.0,
            deformation_scale: 1.0,
            cursor_coordinates: None,
            model_dirty: false,
            editor_message: None,
        };

        let current_path = PathBuf::from(app.model_path_input.clone());
        app.select_available_json_model(&current_path);

        if let Some(path) = input_path {
            app.load_model(path);
        }

        app
    }

    fn load_model(&mut self, path: PathBuf) {
        let input = path.display().to_string();

        match LoadedModel::from_json_file(&path) {
            Ok(model) => {
                self.saved_snapshot = Some(Model2DInput::from_model(&model.model));
                self.loaded_model = Some(model);
                self.load_error = None;
                self.model_path_input = input;
                self.view_state.fit();
                self.selection = None;
                self.model_dirty = false;
                self.history.clear();
                self.results = None;
                self.analysis_task = None;
                self.result_field = ResultField::Model;
                self.element_draft_points.clear();
                self.sync_editor_defaults_from_loaded_model();
                self.select_available_json_model(&path);
            }
            Err(error) => {
                self.load_error = Some(error);
            }
        }
    }

    fn load_model_from_input(&mut self) {
        let trimmed = self.model_path_input.trim();

        if trimmed.is_empty() {
            self.load_error = Some("Model path cannot be empty.".to_owned());
            return;
        }

        self.load_model(PathBuf::from(trimmed));
    }

    fn open_model_dialog(&mut self) {
        if let Some(path) = rfd::FileDialog::new().add_filter("RustyFEM JSON", &["json"]).pick_file() {
            self.load_model(path);
        }
    }

    fn save_model(&mut self) {
        let Some(loaded_model) = &self.loaded_model else {
            self.set_editor_error("No model to save.");
            return;
        };

        if loaded_model.path == Path::new(EMPTY_CANVAS_MODEL_PATH) {
            self.save_model_as();
            return;
        }

        self.save_model_to(loaded_model.path.clone());
    }

    fn save_model_as(&mut self) {
        let suggested_name = self
            .loaded_model
            .as_ref()
            .and_then(|model| model.path.file_name())
            .and_then(|name| name.to_str())
            .unwrap_or(EMPTY_CANVAS_MODEL_PATH);
        let Some(path) =
            rfd::FileDialog::new().add_filter("RustyFEM JSON", &["json"]).set_file_name(suggested_name).save_file()
        else {
            return;
        };

        self.save_model_to(path);
    }

    fn save_model_to(&mut self, path: PathBuf) {
        let result = self
            .loaded_model
            .as_ref()
            .ok_or_else(|| "No model to save.".to_owned())
            .and_then(|loaded_model| write_model_json(&path, &loaded_model.model));

        match result {
            Ok(()) => {
                if let Some(loaded_model) = self.loaded_model.as_mut() {
                    loaded_model.path = path.clone();
                }
                self.model_path_input = path.display().to_string();
                self.model_dirty = false;
                self.saved_snapshot = self.loaded_model.as_ref().map(|model| Model2DInput::from_model(&model.model));
                self.set_editor_info(format!("Saved {}.", path.display()));
            }
            Err(error) => self.set_editor_error(error),
        }
    }

    fn clear_canvas(&mut self) {
        let mut model = Model2D::new();
        let material = Material2D::new(210_000.0, 0.3, 0.0).expect("default drawing material should be valid");
        let section = Section2D::PlaneStress(
            PlaneStressSection2D::new(1.0).expect("default drawing plane-stress section should be valid"),
        );
        let truss_section =
            Section2D::Truss(TrussSection2D::new(1.0).expect("default drawing truss section should be valid"));
        let beam_section = Section2D::Beam(
            BeamSection2D::new_with_section_height(1.0, 1.0, 1.0)
                .expect("default drawing beam section should be valid"),
        );

        model
            .add_material(DEFAULT_DRAWING_MATERIAL_ID, material)
            .expect("empty drawing model should accept default material");
        model
            .add_section(DEFAULT_DRAWING_SECTION_ID, section)
            .expect("empty drawing model should accept default section");
        model
            .add_section(DEFAULT_TRUSS_SECTION_ID, truss_section)
            .expect("empty drawing model should accept default truss section");
        model
            .add_section(DEFAULT_BEAM_SECTION_ID, beam_section)
            .expect("empty drawing model should accept default beam section");

        self.loaded_model = Some(LoadedModel::from_model(PathBuf::from(EMPTY_CANVAS_MODEL_PATH), model));
        self.load_error = None;
        self.model_path_input = EMPTY_CANVAS_MODEL_PATH.to_owned();
        self.selected_json_model = None;
        self.view_state.fit();
        self.selection = None;
        self.model_dirty = true;
        self.saved_snapshot = None;
        self.history.clear();
        self.results = None;
        self.analysis_task = None;
        self.result_field = ResultField::Model;
        self.element_draft_points.clear();
        self.sync_editor_defaults_from_loaded_model();
        self.set_editor_info(
            "Cleared canvas. Started a new model with default plane-stress, truss, and beam sections.",
        );
    }

    fn refresh_available_json_models(&mut self) {
        self.available_json_models = discover_json_models();
        let current_path = PathBuf::from(self.model_path_input.clone());
        self.select_available_json_model(&current_path);
    }

    fn select_available_json_model(&mut self, path: &Path) {
        let path_text = path.display().to_string();

        self.selected_json_model = self
            .available_json_models
            .iter()
            .position(|candidate| candidate == path || candidate.display().to_string() == path_text);
    }

    fn sync_editor_defaults_from_loaded_model(&mut self) {
        let Some(loaded_model) = &self.loaded_model else {
            return;
        };

        let selected_material_exists = self.element_material_id_input.parse::<usize>().ok().is_some_and(|selected| {
            loaded_model.model.materials().materials().iter().any(|(material_id, _)| *material_id == selected)
        });

        if !selected_material_exists && let Some((material_id, _)) = loaded_model.model.materials().materials().first()
        {
            self.element_material_id_input = material_id.to_string();
        }

        self.select_compatible_section();
    }

    fn select_compatible_section(&mut self) {
        let Some(loaded_model) = &self.loaded_model else {
            return;
        };
        let selected_section_is_compatible =
            self.element_section_id_input.parse::<usize>().ok().is_some_and(|selected| {
                loaded_model
                    .model
                    .sections()
                    .sections()
                    .iter()
                    .any(|(id, section)| *id == selected && section_is_compatible(self.element_kind, section))
            });

        if selected_section_is_compatible {
            return;
        }

        let section_id = loaded_model
            .model
            .sections()
            .sections()
            .iter()
            .find_map(|(id, section)| section_is_compatible(self.element_kind, section).then_some(*id));

        if let Some(section_id) = section_id {
            self.element_section_id_input = section_id.to_string();
        }
    }

    fn sidebar(&mut self, ui: &mut egui::Ui) {
        ui.add_space(4.0);
        ui.heading(RichText::new("RustyFEM").color(TEXT));
        ui.horizontal_wrapped(|ui| {
            status_chip(ui, self.work_mode.label(), ChipTone::Accent);

            if self.model_dirty {
                status_chip(ui, "unsaved", ChipTone::Warning);
            }
        });

        sidebar_section(ui, "Workflow", |ui| self.work_mode_panel(ui));
        sidebar_section(ui, "Model", |ui| self.model_loading_panel(ui));
        sidebar_section(ui, "Edit", |ui| self.edit_panel(ui));
        sidebar_section(ui, "View", |ui| self.view_panel(ui));

        if self.work_mode == WorkMode::DrawFem {
            sidebar_section(ui, "Draw", |ui| self.draw_fem_panel(ui));
        }

        sidebar_section(ui, "Analysis", |ui| self.analysis_panel(ui));

        sidebar_section(ui, "Summary", |ui| self.model_summary_panel(ui));
    }

    fn work_mode_panel(&mut self, ui: &mut egui::Ui) {
        ui.horizontal(|ui| {
            ui.selectable_value(&mut self.work_mode, WorkMode::Inspect, WorkMode::Inspect.label());
            ui.selectable_value(&mut self.work_mode, WorkMode::DrawFem, WorkMode::DrawFem.label());
        });
    }

    fn view_panel(&mut self, ui: &mut egui::Ui) {
        ui.horizontal_wrapped(|ui| {
            if ui.button("Fit view").clicked() {
                self.view_state.fit();
            }

            if ui.button("Clear selection").clicked() {
                self.selection = None;
            }
        });

        ui.add_space(6.0);
        ui.checkbox(&mut self.show_grid, "Grid");
        ui.checkbox(&mut self.snap_to_grid, "Snap to grid");
        ui.horizontal(|ui| {
            ui.label("Spacing");
            ui.add(egui::DragValue::new(&mut self.grid_spacing).range(1e-9..=f64::MAX).speed(0.1));
        });
        ui.add_space(4.0);
        ui.checkbox(&mut self.view_options.show_mesh, "Mesh");
        ui.checkbox(&mut self.view_options.show_nodes, "Nodes");
        ui.checkbox(&mut self.view_options.show_boundary_edges, "Boundary");
        ui.checkbox(&mut self.view_options.show_constraints, "Constraints");
        ui.checkbox(&mut self.view_options.show_loads, "Loads");
        ui.checkbox(&mut self.view_options.show_node_ids, "Node IDs");
        ui.checkbox(&mut self.view_options.show_element_ids, "Element IDs");
    }

    fn edit_panel(&mut self, ui: &mut egui::Ui) {
        ui.horizontal_wrapped(|ui| {
            if ui.add_enabled(self.history.can_undo(), egui::Button::new("Undo")).clicked() {
                self.undo();
            }
            if ui.add_enabled(self.history.can_redo(), egui::Button::new("Redo")).clicked() {
                self.redo();
            }
            if ui.add_enabled(self.selection.is_some(), egui::Button::new("Delete")).clicked() {
                self.delete_selection();
            }
        });
    }

    fn model_loading_panel(&mut self, ui: &mut egui::Ui) {
        muted_label(ui, "JSON path");
        let path_response =
            ui.add_sized([ui.available_width(), 24.0], egui::TextEdit::singleline(&mut self.model_path_input));
        let enter_pressed = path_response.lost_focus() && ui.input(|input| input.key_pressed(egui::Key::Enter));
        let mut load_requested = enter_pressed;
        let mut selected_load_path = None;

        ui.horizontal_wrapped(|ui| {
            if ui.button("Open...").clicked() {
                self.open_model_dialog();
            }

            if ui.button("Load").clicked() {
                load_requested = true;
            }

            if ui.add_enabled(self.loaded_model.is_some(), egui::Button::new("Save")).clicked() {
                self.save_model();
            }

            if ui.add_enabled(self.loaded_model.is_some(), egui::Button::new("Save as...")).clicked() {
                self.save_model_as();
            }

            if ui.button("Clear canvas").clicked() {
                self.clear_canvas();
            }

            if ui.button("Scan examples").clicked() {
                self.refresh_available_json_models();
            }
        });

        if self.available_json_models.is_empty() {
            ui.add_space(4.0);
            ui.label(RichText::new("No JSON files found in examples/.").color(MUTED_TEXT));
        } else {
            let selected_text = self
                .selected_json_model
                .and_then(|index| self.available_json_models.get(index))
                .map(|path| short_path_label(path))
                .unwrap_or_else(|| "Choose example".to_owned());

            ui.add_space(6.0);
            muted_label(ui, "Examples");
            egui::ComboBox::from_id_salt("available_json_models")
                .selected_text(selected_text)
                .width(ui.available_width())
                .show_ui(ui, |ui| {
                    for (index, path) in self.available_json_models.iter().enumerate() {
                        ui.selectable_value(&mut self.selected_json_model, Some(index), short_path_label(path))
                            .on_hover_text(path.display().to_string());
                    }
                });

            if ui.button("Load selected").clicked() {
                selected_load_path =
                    self.selected_json_model.and_then(|index| self.available_json_models.get(index).cloned());
            }
        }

        if load_requested {
            self.load_model_from_input();
        }

        if let Some(path) = selected_load_path {
            self.model_path_input = path.display().to_string();
            self.load_model(path);
        }

        if let Some(error) = &self.load_error {
            ui.add_space(8.0);
            ui.colored_label(ERROR, error);
        }
    }

    fn model_summary_panel(&self, ui: &mut egui::Ui) {
        let Some(model) = &self.loaded_model else {
            ui.label(RichText::new("No model loaded.").color(MUTED_TEXT));

            if self.load_error.is_none() {
                ui.add_space(4.0);
                ui.label(RichText::new(format!("Default: {DEFAULT_MODEL_PATH}")).color(MUTED_TEXT));
            }

            return;
        };

        ui.label(RichText::new(short_path_label(&model.path)).strong().color(TEXT))
            .on_hover_text(model.path.display().to_string());

        ui.add_space(6.0);
        Grid::new("model_summary").num_columns(2).spacing([18.0, 4.0]).striped(true).show(ui, |ui| {
            stat_row(ui, "Nodes", model.model.nodes().len());
            stat_row(ui, "Elements", model.model.elements().len());
            stat_row(ui, "Constraints", model.model.constraints().len());
            stat_row(ui, "Nodal loads", model.model.loads().len());
            stat_row(ui, "Element loads", model.model.element_loads().len());
            stat_row(ui, "Materials", model.model.materials().materials().len());
            stat_row(ui, "Sections", model.model.sections().sections().len());
            property_row(ui, "Solver", format!("{:?}", model.model.analysis_settings().solver()));
        });

        if let Some(bounds) = model.bounds {
            ui.add_space(8.0);
            Grid::new("model_bounds").num_columns(2).spacing([18.0, 4.0]).show(ui, |ui| {
                property_row(ui, "x", format!("{:.4} .. {:.4}", bounds.min_x, bounds.max_x));
                property_row(ui, "y", format!("{:.4} .. {:.4}", bounds.min_y, bounds.max_y));
            });
        }
    }

    fn analysis_panel(&mut self, ui: &mut egui::Ui) {
        let Some(current_solver) = self.loaded_model.as_ref().map(|model| model.model.analysis_settings().solver())
        else {
            ui.add_enabled(false, egui::Button::new("Solve"));
            return;
        };
        let mut selected_solver = current_solver;
        let analysis_running = self.analysis_task.is_some();

        ui.add_enabled_ui(!analysis_running, |ui| {
            egui::ComboBox::from_id_salt("solver_kind")
                .selected_text(match selected_solver {
                    SolverKind2D::Dense => "Dense",
                    SolverKind2D::Sparse => "Sparse",
                })
                .width(ui.available_width())
                .show_ui(ui, |ui| {
                    ui.selectable_value(&mut selected_solver, SolverKind2D::Dense, "Dense");
                    ui.selectable_value(&mut selected_solver, SolverKind2D::Sparse, "Sparse");
                });
        });

        if !analysis_running && selected_solver != current_solver {
            self.apply_input_edit("Changed solver settings.", |input| {
                input.analysis_settings.solver = Some(match selected_solver {
                    SolverKind2D::Dense => rusty_fem::io::SolverKind2DInput::Dense,
                    SolverKind2D::Sparse => rusty_fem::io::SolverKind2DInput::Sparse,
                });
                Ok(())
            });
        }

        if ui.add_enabled(!analysis_running, egui::Button::new("Solve")).clicked() {
            self.start_analysis(ui.ctx().clone());
        }

        if let Some(task) = &self.analysis_task {
            ui.add_space(8.0);
            let progress_text = format!("{}  {:.0}%", task.phase.label(), task.progress * 100.0);
            ui.add(
                egui::ProgressBar::new(task.progress)
                    .desired_width(ui.available_width())
                    .animate(true)
                    .text(progress_text),
            );
            ui.add_space(4.0);
            Grid::new("analysis_progress").num_columns(2).spacing([18.0, 4.0]).show(ui, |ui| {
                property_row(ui, "elapsed", format!("{:.1} s", task.started.elapsed().as_secs_f64()));
                property_row(ui, "iteration", task.iterations.to_string());
                if let Some(residual) = task.residual_norm {
                    property_row(ui, "residual", format_value(residual));
                }
                if let Some(relative_residual) = task.relative_residual {
                    property_row(ui, "relative residual", format_value(relative_residual));
                }
                property_row(ui, "target", format_value(task.tolerance));
                if task.phase == AnalysisPhase::Recovering {
                    property_row(ui, "elements", format!("{} / {}", task.recovered_elements, task.total_elements));
                }
            });
        }

        if self.results.is_some() {
            ui.add_space(8.0);
            let previous_field = self.result_field;
            egui::ComboBox::from_id_salt("result_field")
                .selected_text(self.result_field.label())
                .width(ui.available_width())
                .show_ui(ui, |ui| {
                    for field in ResultField::ALL {
                        ui.selectable_value(&mut self.result_field, field, field.label());
                    }
                });

            if self.result_field != previous_field {
                self.use_automatic_range_as_manual();
            }

            if self.result_field == ResultField::Displacement {
                ui.horizontal(|ui| {
                    ui.label("Scale");
                    ui.add(egui::DragValue::new(&mut self.deformation_scale).range(0.0..=1e12).speed(0.1));
                });
            } else if self.result_field.is_scalar_contour() {
                let auto_response = ui.checkbox(&mut self.contour_auto_range, "Automatic contour range");

                if auto_response.changed() && !self.contour_auto_range {
                    self.use_automatic_range_as_manual();
                }
                if !self.contour_auto_range {
                    Grid::new("manual_contour_range").num_columns(2).spacing([18.0, 4.0]).show(ui, |ui| {
                        ui.label("minimum");
                        ui.add(egui::DragValue::new(&mut self.contour_manual_minimum).speed(0.1));
                        ui.end_row();
                        ui.label("maximum");
                        ui.add(egui::DragValue::new(&mut self.contour_manual_maximum).speed(0.1));
                        ui.end_row();
                    });

                    if self.contour_manual_maximum <= self.contour_manual_minimum {
                        ui.label(RichText::new("Maximum must be greater than minimum.").color(ERROR));
                    }
                }
            }

            let results = self.results.as_ref().expect("results were checked above");
            ui.add_space(6.0);
            Grid::new("analysis_summary").num_columns(2).spacing([18.0, 4.0]).show(ui, |ui| {
                property_row(ui, "max |u|", format_value(results.max_displacement));
                property_row(ui, "max VM stress", format_value(results.max_von_mises_stress));
                property_row(ui, "max VM strain", format_value(results.max_equivalent_strain));

                if let Some(iterations) = results.iterations {
                    property_row(ui, "iterations", iterations.to_string());
                }
                if let Some(residual) = results.relative_residual {
                    property_row(ui, "relative residual", format_value(residual));
                }
                if let Some(termination) = results.termination {
                    property_row(ui, "status", termination);
                }
            });
        }
    }

    fn draw_fem_panel(&mut self, ui: &mut egui::Ui) {
        let previous_tool = self.draw_tool;
        muted_label(ui, "Tool");
        egui::ComboBox::from_id_salt("draw_tool")
            .selected_text(self.draw_tool.label())
            .width(ui.available_width())
            .show_ui(ui, |ui| {
                for tool in DrawTool::ALL {
                    ui.selectable_value(&mut self.draw_tool, tool, tool.label());
                }
            });

        if self.draw_tool != previous_tool {
            self.element_draft_points.clear();
        }

        ui.add_space(8.0);

        match self.draw_tool {
            DrawTool::Select | DrawTool::InsertNode | DrawTool::MoveNode => {}
            DrawTool::InsertElement => self.insert_element_panel(ui),
            DrawTool::Constraint => self.constraint_panel(ui),
            DrawTool::NodalLoad => self.nodal_load_panel(ui),
            DrawTool::EdgeTraction => self.edge_traction_panel(ui),
        }

        if let Some(message) = &self.editor_message {
            ui.add_space(8.0);
            match message.kind {
                EditorMessageKind::Info => {
                    ui.colored_label(SUCCESS, &message.text);
                }
                EditorMessageKind::Error => {
                    ui.colored_label(ERROR, &message.text);
                }
            }
        }
    }

    fn insert_element_panel(&mut self, ui: &mut egui::Ui) {
        muted_label(ui, "Element");
        let previous_kind = self.element_kind;
        egui::ComboBox::from_id_salt("element_kind")
            .selected_text(self.element_kind.label())
            .width(ui.available_width())
            .show_ui(ui, |ui| {
                for kind in ElementDraftKind::ALL {
                    ui.selectable_value(&mut self.element_kind, kind, kind.label());
                }
            });

        if self.element_kind != previous_kind {
            self.element_draft_points.clear();
            self.select_compatible_section();
        }

        let material_ids = self
            .loaded_model
            .as_ref()
            .map(|model| model.model.materials().materials().iter().map(|(id, _)| *id).collect::<Vec<_>>())
            .unwrap_or_default();
        let section_ids = self
            .loaded_model
            .as_ref()
            .map(|model| {
                model
                    .model
                    .sections()
                    .sections()
                    .iter()
                    .filter_map(|(id, section)| section_is_compatible(self.element_kind, section).then_some(*id))
                    .collect::<Vec<_>>()
            })
            .unwrap_or_default();

        ui.horizontal(|ui| {
            egui::ComboBox::from_id_salt("element_material")
                .selected_text(format!("Material {}", self.element_material_id_input))
                .show_ui(ui, |ui| {
                    for material_id in material_ids {
                        ui.selectable_value(
                            &mut self.element_material_id_input,
                            material_id.to_string(),
                            format!("Material {material_id}"),
                        );
                    }
                });
            egui::ComboBox::from_id_salt("element_section")
                .selected_text(format!("Section {}", self.element_section_id_input))
                .show_ui(ui, |ui| {
                    for section_id in section_ids {
                        ui.selectable_value(
                            &mut self.element_section_id_input,
                            section_id.to_string(),
                            format!("Section {section_id}"),
                        );
                    }
                });
        });

        let point_count = self.element_kind.placement_point_count();
        ui.add_space(4.0);
        ui.horizontal_wrapped(|ui| {
            ui.label(format!("Corners {} / {point_count}", self.element_draft_points.len()));
            for index in 0..self.element_draft_points.len() {
                status_chip(ui, &(index + 1).to_string(), ChipTone::Accent);
            }
        });

        if !self.element_draft_points.is_empty() {
            ui.horizontal_wrapped(|ui| {
                if ui.button("Cancel").clicked() {
                    self.element_draft_points.clear();
                }
            });
        }
    }

    fn constraint_panel(&mut self, ui: &mut egui::Ui) {
        let ux_available = self.dof_available_for_current_target(Dof2D::Ux);
        let uy_available = self.dof_available_for_current_target(Dof2D::Uy);
        let rz_available = self.dof_available_for_current_target(Dof2D::Rz);

        if !ux_available {
            self.constraint_ux = false;
        }
        if !uy_available {
            self.constraint_uy = false;
        }
        if !rz_available {
            self.constraint_rz = false;
        }

        ui.horizontal(|ui| {
            ui.add_enabled(ux_available, egui::Checkbox::new(&mut self.constraint_ux, "Ux"));
            ui.add_enabled(uy_available, egui::Checkbox::new(&mut self.constraint_uy, "Uy"));
            ui.add_enabled(rz_available, egui::Checkbox::new(&mut self.constraint_rz, "Rz"));
        });
        ui.horizontal(|ui| {
            ui.label("Value");
            ui.text_edit_singleline(&mut self.constraint_value_input);
        });
    }

    fn nodal_load_panel(&mut self, ui: &mut egui::Ui) {
        let ux_available = self.dof_available_for_current_target(Dof2D::Ux);
        let uy_available = self.dof_available_for_current_target(Dof2D::Uy);
        let rz_available = self.dof_available_for_current_target(Dof2D::Rz);

        if !self.dof_available_for_current_target(self.load_dof) {
            self.load_dof = if uy_available { Dof2D::Uy } else { Dof2D::Ux };
        }

        ui.horizontal(|ui| {
            egui::ComboBox::from_id_salt("load_dof").selected_text(self.load_dof.name()).show_ui(ui, |ui| {
                if ux_available {
                    ui.selectable_value(&mut self.load_dof, Dof2D::Ux, "Ux");
                }
                if uy_available {
                    ui.selectable_value(&mut self.load_dof, Dof2D::Uy, "Uy");
                }
                if rz_available {
                    ui.selectable_value(&mut self.load_dof, Dof2D::Rz, "Rz");
                }
            });
            ui.label("Value");
            ui.text_edit_singleline(&mut self.load_value_input);
        });
    }

    fn edge_traction_panel(&mut self, ui: &mut egui::Ui) {
        ui.horizontal(|ui| {
            ui.label("tx");
            ui.text_edit_singleline(&mut self.traction_x_input);
            ui.label("ty");
            ui.text_edit_singleline(&mut self.traction_y_input);
        });
    }

    fn dof_available_for_current_target(&self, dof: Dof2D) -> bool {
        let Some(model) = self.loaded_model.as_ref().map(|loaded| &loaded.model) else {
            return false;
        };

        match self.selection {
            Some(SelectedEntity::Node(node_id)) => node_has_dof(model, node_id, dof),
            Some(SelectedEntity::Element(_)) | Some(SelectedEntity::Edge { .. }) | None => {
                model.elements().iter().any(|element| element.dofs_per_node().contains(&dof))
            }
        }
    }

    fn apply_input_edit(
        &mut self, description: impl Into<String>, edit: impl FnOnce(&mut Model2DInput) -> Result<(), String>,
    ) -> bool {
        let Some(loaded_model) = &self.loaded_model else {
            self.set_editor_error("Load or create a model first.");
            return false;
        };
        let previous = Model2DInput::from_model(&loaded_model.model);
        let mut next = previous.clone();

        if let Err(error) = edit(&mut next) {
            self.set_editor_error(error);
            return false;
        }

        match next.into_model() {
            Ok(model) => {
                self.history.record_snapshot(previous);
                if let Some(loaded_model) = self.loaded_model.as_mut() {
                    loaded_model.model = model;
                    loaded_model.refresh_bounds();
                }
                self.after_model_change();
                self.set_editor_info(description);
                true
            }
            Err(error) => {
                self.set_editor_error(format!("Model edit rejected: {error}"));
                false
            }
        }
    }

    fn after_model_change(&mut self) {
        self.results = None;
        self.result_field = ResultField::Model;
        self.model_dirty = self.loaded_model.as_ref().is_some_and(|loaded_model| {
            self.saved_snapshot.as_ref().is_none_or(|saved| Model2DInput::from_model(&loaded_model.model) != *saved)
        });
        self.sync_editor_defaults_from_loaded_model();
    }

    fn handle_canvas_action(&mut self, action: CanvasAction) {
        match action {
            CanvasAction::InsertNode { x, y } => self.insert_node_at(x, y),
            CanvasAction::PlaceElementPoint { x, y } => self.place_element_point(x, y),
            CanvasAction::MoveNode { node_id, x, y } => self.move_node(node_id, x, y),
            CanvasAction::AddConstraint(node_id) => self.add_constraints(node_id),
            CanvasAction::AddNodalLoad(node_id) => self.add_nodal_load(node_id),
            CanvasAction::AddEdgeTraction { element_id, node_ids } => self.add_edge_traction(element_id, node_ids),
        }
    }

    fn insert_node_at(&mut self, x: f64, y: f64) {
        let id = self.next_node_id();

        if self.apply_input_edit(format!("Inserted node {id}."), |input| {
            input.nodes.push(Node2DInput { id, x, y });
            Ok(())
        }) {
            self.selection = Some(SelectedEntity::Node(id));
        }
    }

    fn place_element_point(&mut self, x: f64, y: f64) {
        let point = [x, y];

        if self.element_draft_points.iter().any(|candidate| coordinates_match(*candidate, point)) {
            self.set_editor_error("This corner is already part of the element draft.");
            return;
        }

        self.element_draft_points.push(point);

        if self.element_draft_points.len() == self.element_kind.placement_point_count() {
            let points = std::mem::take(&mut self.element_draft_points);
            self.insert_element_from_points(points);
        }
    }

    fn insert_element_from_points(&mut self, mut points: Vec<[f64; 2]>) {
        let id = self.next_element_id();
        let material_id = match parse_usize_field(&self.element_material_id_input, "material ID") {
            Ok(id) => id,
            Err(error) => {
                self.set_editor_error(error);
                return;
            }
        };
        let section_id = match parse_usize_field(&self.element_section_id_input, "section ID") {
            Ok(id) => id,
            Err(error) => {
                self.set_editor_error(error);
                return;
            }
        };
        let kind = self.element_kind;
        normalize_corner_order(kind, &mut points);

        if self.apply_input_edit(format!("Inserted {} element {id}.", kind.label()), |input| {
            let node_ids = connectivity_from_points(input, kind, &points)?;
            let element = element_from_draft(id, kind, node_ids, material_id, section_id)?;
            let element_input = ElementType2DInput::from(&element);
            input.elements.push(element_input);
            Ok(())
        }) {
            self.selection = Some(SelectedEntity::Element(id));
        }
    }

    fn move_node(&mut self, node_id: usize, x: f64, y: f64) {
        if self.apply_input_edit(format!("Moved node {node_id}."), |input| {
            let node = input
                .nodes
                .iter_mut()
                .find(|node| node.id == node_id)
                .ok_or_else(|| format!("Unknown node {node_id}."))?;
            node.x = x;
            node.y = y;
            Ok(())
        }) {
            self.selection = Some(SelectedEntity::Node(node_id));
        }
    }

    fn add_constraints(&mut self, node_id: usize) {
        let value = match parse_f64_field(&self.constraint_value_input, "constraint value") {
            Ok(value) => value,
            Err(error) => {
                self.set_editor_error(error);
                return;
            }
        };
        let selected = [
            (self.constraint_ux, Dof2D::Ux, Dof2DInput::Ux),
            (self.constraint_uy, Dof2D::Uy, Dof2DInput::Uy),
            (self.constraint_rz, Dof2D::Rz, Dof2DInput::Rz),
        ];

        if let Some((_, unsupported, _)) = selected.iter().find(|(enabled, dof, _)| {
            *enabled && !self.loaded_model.as_ref().is_some_and(|loaded| node_has_dof(&loaded.model, node_id, *dof))
        }) {
            self.set_editor_error(format!("Node {node_id} does not have degree of freedom {}.", unsupported.name()));
            return;
        }

        self.apply_input_edit(format!("Updated constraints on node {node_id}."), |input| {
            let previous_count = input.constraints.len();
            input.constraints.retain(|constraint| constraint.node != node_id);
            let removed = previous_count - input.constraints.len();

            for (_, _, dof) in selected.iter().filter(|(enabled, _, _)| *enabled) {
                input.constraints.push(DisplacementConstraint2DInput { node: node_id, dof: *dof, value });
            }

            (removed > 0 || selected.iter().any(|(enabled, _, _)| *enabled))
                .then_some(())
                .ok_or_else(|| "Node has no constraints to clear.".to_owned())
        });
    }

    fn add_nodal_load(&mut self, node_id: usize) {
        if !self.loaded_model.as_ref().is_some_and(|loaded| node_has_dof(&loaded.model, node_id, self.load_dof)) {
            self.set_editor_error(format!("Node {node_id} does not have degree of freedom {}.", self.load_dof.name()));
            return;
        }

        let value = match parse_f64_field(&self.load_value_input, "load value") {
            Ok(value) => value,
            Err(error) => {
                self.set_editor_error(error);
                return;
            }
        };
        let dof = match self.load_dof {
            Dof2D::Ux => Dof2DInput::Ux,
            Dof2D::Uy => Dof2DInput::Uy,
            Dof2D::Rz => Dof2DInput::Rz,
        };

        self.apply_input_edit(format!("Applied nodal load to node {node_id}."), |input| {
            input.nodal_loads.push(NodalLoad2DInput { node: node_id, dof, value });
            Ok(())
        });
    }

    fn add_edge_traction(&mut self, element_id: usize, node_ids: [usize; 2]) {
        let tx = match parse_f64_field(&self.traction_x_input, "traction tx") {
            Ok(value) => value,
            Err(error) => {
                self.set_editor_error(error);
                return;
            }
        };
        let ty = match parse_f64_field(&self.traction_y_input, "traction ty") {
            Ok(value) => value,
            Err(error) => {
                self.set_editor_error(error);
                return;
            }
        };

        self.apply_input_edit(format!("Applied edge traction to element {element_id}."), |input| {
            input.loads.push(ElementLoad2DInput {
                kind: ElementLoad2DInputKind::EdgeTraction {
                    element: element_id,
                    edge: node_ids,
                    coordinate_system: LoadCoordinateSystem2DInput::Global,
                    tx,
                    ty,
                },
            });
            Ok(())
        });
    }

    fn delete_selection(&mut self) {
        let Some(selection) = self.selection else {
            return;
        };

        let changed = match selection {
            SelectedEntity::Node(node_id) => self.apply_input_edit(format!("Deleted node {node_id}."), |input| {
                if input.elements.iter().any(|element| element_input_nodes(element).contains(&node_id)) {
                    return Err(format!("Node {node_id} is still used by an element."));
                }
                input.nodes.retain(|node| node.id != node_id);
                input.constraints.retain(|constraint| constraint.node != node_id);
                input.nodal_loads.retain(|load| load.node != node_id);
                input
                    .loads
                    .retain(|load| !matches!(load.kind, ElementLoad2DInputKind::Nodal { node, .. } if node == node_id));
                Ok(())
            }),
            SelectedEntity::Element(element_id) => {
                self.apply_input_edit(format!("Deleted element {element_id}."), |input| {
                    input.elements.retain(|element| element_input_id(element) != element_id);
                    input.loads.retain(|load| element_load_input_id(load) != Some(element_id));
                    Ok(())
                })
            }
            SelectedEntity::Edge { .. } => {
                self.set_editor_error("Select a node or element to delete it.");
                false
            }
        };

        if changed {
            self.selection = None;
        }
    }

    fn undo(&mut self) {
        let snapshot = self.loaded_model.as_ref().and_then(|model| self.history.undo(&model.model));
        self.restore_snapshot(snapshot, "Undo");
    }

    fn redo(&mut self) {
        let snapshot = self.loaded_model.as_ref().and_then(|model| self.history.redo(&model.model));
        self.restore_snapshot(snapshot, "Redo");
    }

    fn restore_snapshot(&mut self, snapshot: Option<Model2DInput>, label: &str) {
        let Some(snapshot) = snapshot else {
            return;
        };

        match snapshot.into_model() {
            Ok(model) => {
                if let Some(loaded_model) = self.loaded_model.as_mut() {
                    loaded_model.model = model;
                    loaded_model.refresh_bounds();
                }
                self.selection = None;
                self.element_draft_points.clear();
                self.after_model_change();
                self.set_editor_info(label);
            }
            Err(error) => self.set_editor_error(format!("Could not restore history: {error}")),
        }
    }

    fn start_analysis(&mut self, repaint_context: egui::Context) {
        if self.analysis_task.is_some() {
            return;
        }
        let Some(model) = self.loaded_model.as_ref() else {
            self.set_editor_error("No model to solve.");
            return;
        };
        let model_snapshot = Model2DInput::from_model(&model.model);
        let worker_input = model_snapshot.clone();
        let tolerance = model.model.analysis_settings().cg_tolerance();
        let (sender, receiver) = mpsc::channel();
        let spawn_result = std::thread::Builder::new().name("rusty-fem-analysis".to_owned()).spawn(move || {
            let result = worker_input.into_model().map_err(|error| error.to_string()).and_then(|worker_model| {
                let mut last_progress_update = Instant::now() - Duration::from_millis(100);

                AnalysisResults::solve_with_progress(&worker_model, |progress| {
                    let should_send = match progress {
                        AnalysisProgress::Phase(_) => true,
                        AnalysisProgress::Iteration(sample) => {
                            sample.iterations == 0
                                || sample.relative_residual_norm <= tolerance
                                || last_progress_update.elapsed() >= Duration::from_millis(50)
                        }
                        AnalysisProgress::Recovery { completed, total } => {
                            completed == 0
                                || completed == total
                                || last_progress_update.elapsed() >= Duration::from_millis(50)
                        }
                    };

                    if should_send {
                        if !matches!(progress, AnalysisProgress::Phase(_)) {
                            last_progress_update = Instant::now();
                        }
                        if sender.send(AnalysisTaskMessage::Progress(progress)).is_ok() {
                            repaint_context.request_repaint();
                        }
                    }
                })
            });

            let _ = sender.send(AnalysisTaskMessage::Finished(Box::new(result)));
            repaint_context.request_repaint();
        });

        match spawn_result {
            Ok(_) => {
                self.results = None;
                self.result_field = ResultField::Model;
                self.analysis_task = Some(AnalysisTask::new(receiver, model_snapshot, tolerance));
                self.set_editor_info("Analysis started.");
            }
            Err(error) => self.set_editor_error(format!("Could not start analysis: {error}")),
        }
    }

    fn poll_analysis(&mut self) {
        let mut completion = None;
        let mut disconnected = false;

        if let Some(task) = &mut self.analysis_task {
            loop {
                match task.receiver.try_recv() {
                    Ok(AnalysisTaskMessage::Progress(progress)) => task.apply_progress(progress),
                    Ok(AnalysisTaskMessage::Finished(result)) => {
                        completion = Some(result);
                        break;
                    }
                    Err(TryRecvError::Empty) => break,
                    Err(TryRecvError::Disconnected) => {
                        disconnected = true;
                        break;
                    }
                }
            }
        }

        if disconnected {
            self.analysis_task = None;
            self.set_editor_error("Analysis worker stopped unexpectedly.");
            return;
        }

        let Some(result) = completion else {
            return;
        };
        let Some(task) = self.analysis_task.take() else {
            return;
        };
        let model_is_unchanged = self
            .loaded_model
            .as_ref()
            .is_some_and(|model| Model2DInput::from_model(&model.model) == task.model_snapshot);

        match *result {
            Ok(results) if model_is_unchanged => {
                self.results = Some(results);
                self.result_field = ResultField::VonMisesStress;
                self.use_automatic_range_as_manual();
                self.set_editor_info(format!("Analysis completed in {:.1} s.", task.started.elapsed().as_secs_f64()));
            }
            Ok(_) => self.set_editor_error("Analysis result was discarded because the model changed."),
            Err(error) => self.set_editor_error(format!("Analysis failed: {error}")),
        }
    }

    fn next_node_id(&self) -> usize {
        self.loaded_model.as_ref().and_then(|model| model.model.nodes().iter().map(|node| node.id()).max()).unwrap_or(0)
            + 1
    }

    fn next_element_id(&self) -> usize {
        self.loaded_model
            .as_ref()
            .and_then(|model| model.model.elements().iter().map(|element| element.id()).max())
            .unwrap_or(0)
            + 1
    }

    fn handle_shortcuts(&mut self, ui: &egui::Ui) {
        let wants_keyboard = ui.ctx().egui_wants_keyboard_input();
        let (undo, redo, save, save_as, delete, cancel) = ui.input(|input| {
            let command = input.modifiers.command;
            (
                command && !input.modifiers.shift && input.key_pressed(egui::Key::Z),
                command
                    && (input.key_pressed(egui::Key::Y) || (input.modifiers.shift && input.key_pressed(egui::Key::Z))),
                command && !input.modifiers.shift && input.key_pressed(egui::Key::S),
                command && input.modifiers.shift && input.key_pressed(egui::Key::S),
                !wants_keyboard && (input.key_pressed(egui::Key::Delete) || input.key_pressed(egui::Key::Backspace)),
                input.key_pressed(egui::Key::Escape),
            )
        });

        if undo {
            self.undo();
        } else if redo {
            self.redo();
        }

        if save_as {
            self.save_model_as();
        } else if save {
            self.save_model();
        }

        if delete {
            self.delete_selection();
        }

        if cancel {
            self.element_draft_points.clear();
        }
    }

    fn set_editor_info(&mut self, text: impl Into<String>) {
        self.editor_message = Some(EditorMessage { text: text.into(), kind: EditorMessageKind::Info });
    }

    fn set_editor_error(&mut self, text: impl Into<String>) {
        self.editor_message = Some(EditorMessage { text: text.into(), kind: EditorMessageKind::Error });
    }

    fn active_scalar_range(&self) -> Option<ScalarRange> {
        let automatic = self.results.as_ref()?.scalar_range(self.result_field)?;

        if self.contour_auto_range
            || !self.contour_manual_minimum.is_finite()
            || !self.contour_manual_maximum.is_finite()
            || self.contour_manual_maximum <= self.contour_manual_minimum
        {
            return Some(automatic);
        }

        Some(ScalarRange { minimum: self.contour_manual_minimum, maximum: self.contour_manual_maximum })
    }

    fn use_automatic_range_as_manual(&mut self) {
        let Some(range) = self.results.as_ref().and_then(|results| results.scalar_range(self.result_field)) else {
            return;
        };

        self.contour_manual_minimum = range.minimum;
        self.contour_manual_maximum = range.maximum;
    }

    fn properties_panel(&self, ui: &mut egui::Ui) {
        ui.heading(RichText::new("Selection").color(TEXT));
        ui.add_space(8.0);

        let Some(loaded_model) = &self.loaded_model else {
            ui.label(RichText::new("No model loaded.").color(MUTED_TEXT));
            return;
        };

        let Some(selection) = self.selection else {
            ui.label(RichText::new("Nothing selected.").color(MUTED_TEXT));
            return;
        };

        status_chip(ui, &selection.label(), ChipTone::Accent);
        ui.add_space(8.0);
        ui.separator();
        ui.add_space(8.0);

        match selection {
            SelectedEntity::Node(node_id) => draw_node_properties(ui, &loaded_model.model, node_id),
            SelectedEntity::Element(element_id) => draw_element_properties(ui, &loaded_model.model, element_id),
            SelectedEntity::Edge { element_id, node_ids } => {
                draw_edge_properties(ui, &loaded_model.model, element_id, node_ids);
            }
        }

        if let Some(results) = &self.results {
            match selection {
                SelectedEntity::Node(node_id) => {
                    let [ux, uy] = results.displacement(node_id);
                    detail_section(ui, "Analysis");
                    Grid::new("selected_node_results").num_columns(2).spacing([18.0, 4.0]).show(ui, |ui| {
                        property_row(ui, "ux", format_value(ux));
                        property_row(ui, "uy", format_value(uy));
                        property_row(ui, "|u|", format_value(ux.hypot(uy)));
                        for field in ResultField::ALL.into_iter().filter(|field| field.is_scalar_contour()) {
                            if let Some(value) = results.nodal_scalar_value(field, node_id) {
                                property_row(ui, field.label(), format_value(value));
                            }
                        }
                    });
                }
                SelectedEntity::Element(element_id) => {
                    if results.von_mises_stress(element_id).is_some() {
                        detail_section(ui, "Analysis");
                        Grid::new("selected_element_results").num_columns(2).spacing([18.0, 4.0]).show(ui, |ui| {
                            for field in ResultField::ALL.into_iter().filter(|field| field.is_scalar_contour()) {
                                if let Some(value) = results.scalar_value(field, element_id) {
                                    property_row(ui, field.label(), format_value(value));
                                }
                            }
                        });
                    }
                }
                SelectedEntity::Edge { .. } => {}
            }
        }
    }
}

impl eframe::App for RustyFemGuiApp {
    fn ui(&mut self, ui: &mut egui::Ui, _frame: &mut eframe::Frame) {
        self.poll_analysis();
        if self.analysis_task.is_some() {
            ui.ctx().request_repaint_after(Duration::from_millis(100));
        }
        self.handle_shortcuts(ui);

        egui::Panel::left("sidebar_panel")
            .exact_size(SIDEBAR_WIDTH)
            .resizable(false)
            .show_separator_line(false)
            .frame(side_panel_frame())
            .show(ui, |ui| {
                egui::ScrollArea::vertical()
                    .id_salt("sidebar_scroll")
                    .auto_shrink([false, false])
                    .show(ui, |ui| self.sidebar(ui));
            });

        egui::Panel::right("properties_panel")
            .exact_size(PROPERTIES_WIDTH)
            .resizable(false)
            .show_separator_line(false)
            .frame(side_panel_frame())
            .show(ui, |ui| {
                egui::ScrollArea::vertical()
                    .id_salt("properties_scroll")
                    .auto_shrink([false, false])
                    .show(ui, |ui| self.properties_panel(ui));
            });

        egui::CentralPanel::default().frame(Frame::new().fill(APP_BACKGROUND).inner_margin(Margin::same(10))).show(
            ui,
            |ui| {
                let scalar_range = self.active_scalar_range();
                let view = CanvasView {
                    options: &self.view_options,
                    editor: CanvasEditor {
                        work_mode: self.work_mode,
                        tool: self.draw_tool,
                        show_grid: self.show_grid,
                        snap_to_grid: self.snap_to_grid,
                        grid_spacing: self.grid_spacing,
                        element_kind: self.element_kind,
                        element_draft_points: &self.element_draft_points,
                    },
                    results: self.results.as_ref(),
                    result_field: self.result_field,
                    scalar_range,
                    deformation_scale: self.deformation_scale,
                };
                let output =
                    draw_model_canvas(ui, self.loaded_model.as_ref(), &mut self.view_state, &mut self.selection, &view);
                self.cursor_coordinates = output.cursor_coordinates;

                if let Some(action) = output.action {
                    self.handle_canvas_action(action);
                }
            },
        );
    }
}

enum ChipTone {
    Accent,
    Warning,
}

fn sidebar_section(ui: &mut egui::Ui, title: &str, add_contents: impl FnOnce(&mut egui::Ui)) {
    ui.add_space(14.0);
    ui.separator();
    ui.add_space(8.0);
    ui.label(RichText::new(title.to_ascii_uppercase()).small().strong().color(MUTED_TEXT));
    ui.add_space(4.0);
    add_contents(ui);
}

fn side_panel_frame() -> Frame {
    Frame::new()
        .fill(PANEL_BACKGROUND)
        .stroke(Stroke::new(1.0, BORDER))
        .corner_radius(6)
        .inner_margin(Margin::same(12))
        .outer_margin(Margin::same(10))
}

fn detail_section(ui: &mut egui::Ui, title: &str) {
    ui.add_space(14.0);
    ui.label(RichText::new(title).strong().color(TEXT));
    ui.add_space(4.0);
}

fn muted_label(ui: &mut egui::Ui, text: &str) {
    ui.label(RichText::new(text).small().color(MUTED_TEXT));
}

fn status_chip(ui: &mut egui::Ui, text: &str, tone: ChipTone) {
    let (fill, stroke, color) = match tone {
        ChipTone::Accent => (Color32::from_rgb(236, 245, 248), Color32::from_rgb(125, 165, 180), ACCENT),
        ChipTone::Warning => (Color32::from_rgb(255, 247, 237), Color32::from_rgb(251, 191, 36), WARNING),
    };

    Frame::new()
        .fill(fill)
        .stroke(Stroke::new(1.0, stroke))
        .corner_radius(4)
        .inner_margin(Margin::symmetric(7, 2))
        .show(ui, |ui| {
            ui.label(RichText::new(text).small().strong().color(color));
        });
}

fn stat_row(ui: &mut egui::Ui, label: &str, value: usize) {
    property_row(ui, label, value.to_string());
}

fn property_row(ui: &mut egui::Ui, label: &str, value: impl Into<String>) {
    ui.label(RichText::new(label).color(MUTED_TEXT));
    ui.label(RichText::new(value.into()).monospace().color(TEXT));
    ui.end_row();
}

fn short_path_label(path: &Path) -> String {
    path.file_name()
        .and_then(|file_name| file_name.to_str())
        .map(str::to_owned)
        .unwrap_or_else(|| path.display().to_string())
}

fn draw_node_properties(ui: &mut egui::Ui, model: &Model2D, node_id: usize) {
    let Some(node) = model.nodes().iter().find(|node| node.id() == node_id) else {
        ui.colored_label(ERROR, "Selected node is not present in the model.");
        return;
    };

    Grid::new("node_properties").num_columns(2).spacing([18.0, 4.0]).show(ui, |ui| {
        property_row(ui, "x", format_value(node.x()));
        property_row(ui, "y", format_value(node.y()));
    });

    detail_section(ui, "Constraints");
    let constraints: Vec<_> = model.constraints().iter().filter(|constraint| constraint.node_id() == node_id).collect();

    if constraints.is_empty() {
        ui.label(RichText::new("none").color(MUTED_TEXT));
    } else {
        Grid::new("node_constraints").num_columns(2).spacing([18.0, 4.0]).show(ui, |ui| {
            for constraint in constraints {
                property_row(ui, constraint.dof().name(), format_value(constraint.displacement()));
            }
        });
    }

    detail_section(ui, "Nodal loads");
    let loads: Vec<_> = model.loads().iter().filter(|load| load.node_id() == node_id).collect();

    if loads.is_empty() {
        ui.label(RichText::new("none").color(MUTED_TEXT));
    } else {
        Grid::new("node_loads").num_columns(2).spacing([18.0, 4.0]).show(ui, |ui| {
            for load in loads {
                property_row(ui, load.dof().name(), format_value(load.value()));
            }
        });
    }
}

fn draw_element_properties(ui: &mut egui::Ui, model: &Model2D, element_id: usize) {
    let Some(element) = model.elements().iter().find(|element| element.id() == element_id) else {
        ui.colored_label(ERROR, "Selected element is not present in the model.");
        return;
    };

    Grid::new("element_properties").num_columns(2).spacing([18.0, 4.0]).show(ui, |ui| {
        property_row(ui, "Type", element.element_type());
        property_row(ui, "Interpolation", format!("{:?}", element.interpolation()));
        property_row(ui, "Nodes", format_node_ids(element.node_ids()));
        property_row(ui, "Material ID", element.material_id().to_string());
        property_row(ui, "Section ID", element.section_id().to_string());
    });

    if let Ok(material) = model.material(element.material_id()) {
        detail_section(ui, "Material");
        Grid::new("element_material").num_columns(2).spacing([18.0, 4.0]).show(ui, |ui| {
            property_row(ui, "E", format_value(material.young_modulus()));
            property_row(ui, "nu", format_value(material.poisson_ratio()));
            property_row(ui, "density", format_value(material.density()));
        });
    }

    if let Ok(section) = model.section(element.section_id()) {
        detail_section(ui, "Section");
        draw_section_properties(ui, section);
    }

    detail_section(ui, "Element loads");
    let loads: Vec<_> = model.element_loads().iter().filter(|load| load.element_id() == element_id).collect();

    if loads.is_empty() {
        ui.label(RichText::new("none").color(MUTED_TEXT));
    } else {
        for load in loads {
            ui.label(format_element_load(load));
        }
    }
}

fn draw_edge_properties(ui: &mut egui::Ui, model: &Model2D, element_id: usize, node_ids: [usize; 2]) {
    Grid::new("edge_properties").num_columns(2).spacing([18.0, 4.0]).show(ui, |ui| {
        property_row(ui, "Element ID", element_id.to_string());
        property_row(ui, "Node IDs", format!("{} - {}", node_ids[0], node_ids[1]));
    });

    let normalized_selected = normalized_edge(node_ids[0], node_ids[1]);
    let selected_element = model.elements().iter().find(|element| element.id() == element_id);
    let matching_loads: Vec<_> = model
        .element_loads()
        .iter()
        .filter_map(|load| match load {
            ElementLoad2D::EdgeTraction(edge_load)
                if edge_load.element_id() == element_id
                    && selected_element
                        .map(|element| {
                            element_edge_segments(element, edge_load.edge_node_ids()).contains(&normalized_selected)
                        })
                        .unwrap_or(false) =>
            {
                Some(edge_load)
            }
            _ => None,
        })
        .collect();

    detail_section(ui, "Edge traction");

    if matching_loads.is_empty() {
        ui.label(RichText::new("none").color(MUTED_TEXT));
    } else {
        for load in matching_loads {
            ui.label(format!(
                "{:?}: tx = {}, ty = {}",
                load.coordinate_system(),
                format_value(load.x_component()),
                format_value(load.y_component())
            ));
        }
    }
}

fn draw_section_properties(ui: &mut egui::Ui, section: &Section2D) {
    Grid::new("section_properties").num_columns(2).spacing([18.0, 4.0]).show(ui, |ui| match section {
        Section2D::Truss(section) => {
            property_row(ui, "Type", "truss");
            property_row(ui, "A", format_value(section.cross_section_area()));
        }
        Section2D::Beam(section) => {
            property_row(ui, "Type", "beam");
            property_row(ui, "A", format_value(section.cross_section_area()));
            property_row(ui, "I", format_value(section.second_moment_of_area()));

            if let Some(height) = section.section_height() {
                property_row(ui, "height", format_value(height));
            }
        }
        Section2D::PlaneStress(section) => {
            property_row(ui, "Type", "plane_stress");
            property_row(ui, "thickness", format_value(section.thickness()));
        }
    });
}

fn format_element_load(load: &ElementLoad2D) -> String {
    match load {
        ElementLoad2D::BeamUniformLine(load) => format!(
            "beam line {:?}: qx = {}, qy = {}",
            load.coordinate_system(),
            format_value(load.x_component()),
            format_value(load.y_component())
        ),
        ElementLoad2D::EdgeTraction(load) => format!(
            "edge {}-{} {:?}: tx = {}, ty = {}",
            load.edge_node_ids()[0],
            load.edge_node_ids()[1],
            load.coordinate_system(),
            format_value(load.x_component()),
            format_value(load.y_component())
        ),
        ElementLoad2D::BodyForce(load) => {
            format!("body force: bx = {}, by = {}", format_value(load.x_component()), format_value(load.y_component()))
        }
        ElementLoad2D::SelfWeight(load) => format!(
            "self weight: ax = {}, ay = {}",
            format_value(load.x_acceleration()),
            format_value(load.y_acceleration())
        ),
    }
}

fn format_node_ids(node_ids: &[usize]) -> String {
    node_ids.iter().map(|node_id| node_id.to_string()).collect::<Vec<_>>().join(", ")
}

fn format_value(value: f64) -> String {
    if value == 0.0 {
        "0".to_owned()
    } else if (1e-4..1e6).contains(&value.abs()) {
        format!("{value:.6}")
    } else {
        format!("{value:.6e}")
    }
}

fn element_input_id(element: &ElementType2DInput) -> usize {
    match element {
        ElementType2DInput::Truss { id, .. }
        | ElementType2DInput::Beam { id, .. }
        | ElementType2DInput::TriangleT3 { id, .. }
        | ElementType2DInput::T6 { id, .. }
        | ElementType2DInput::Q4 { id, .. }
        | ElementType2DInput::Q8 { id, .. } => *id,
    }
}

fn element_input_nodes(element: &ElementType2DInput) -> &[usize] {
    match element {
        ElementType2DInput::Truss { nodes, .. } | ElementType2DInput::Beam { nodes, .. } => nodes,
        ElementType2DInput::TriangleT3 { nodes, .. } => nodes,
        ElementType2DInput::T6 { nodes, .. } => nodes,
        ElementType2DInput::Q4 { nodes, .. } => nodes,
        ElementType2DInput::Q8 { nodes, .. } => nodes,
    }
}

fn element_load_input_id(load: &ElementLoad2DInput) -> Option<usize> {
    match load.kind {
        ElementLoad2DInputKind::BeamUniform { element, .. }
        | ElementLoad2DInputKind::EdgeTraction { element, .. }
        | ElementLoad2DInputKind::BodyForce { element, .. }
        | ElementLoad2DInputKind::SelfWeight { element, .. } => Some(element),
        ElementLoad2DInputKind::Nodal { .. } => None,
    }
}

fn section_is_compatible(kind: ElementDraftKind, section: &Section2D) -> bool {
    match kind {
        ElementDraftKind::Truss => matches!(section, Section2D::Truss(_)),
        ElementDraftKind::Beam => matches!(section, Section2D::Beam(_)),
        ElementDraftKind::T3 | ElementDraftKind::T6 | ElementDraftKind::Q4 | ElementDraftKind::Q8 => {
            matches!(section, Section2D::PlaneStress(_))
        }
    }
}

fn node_has_dof(model: &Model2D, node_id: usize, dof: Dof2D) -> bool {
    model
        .elements()
        .iter()
        .any(|element| element.node_ids().contains(&node_id) && element.dofs_per_node().contains(&dof))
}

fn connectivity_from_points(
    input: &mut Model2DInput, kind: ElementDraftKind, points: &[[f64; 2]],
) -> Result<Vec<usize>, String> {
    if points.len() != kind.placement_point_count() {
        return Err(format!(
            "{} requires {} corner points, got {}.",
            kind.label(),
            kind.placement_point_count(),
            points.len()
        ));
    }

    let mut node_ids = points.iter().map(|&point| find_or_insert_node(input, point)).collect::<Vec<_>>();
    let midside_edges: &[[usize; 2]] = match kind {
        ElementDraftKind::T6 => &[[0, 1], [1, 2], [2, 0]],
        ElementDraftKind::Q8 => &[[0, 1], [1, 2], [2, 3], [3, 0]],
        ElementDraftKind::Truss | ElementDraftKind::Beam | ElementDraftKind::T3 | ElementDraftKind::Q4 => &[],
    };

    for &[first, second] in midside_edges {
        let midpoint = [0.5 * (points[first][0] + points[second][0]), 0.5 * (points[first][1] + points[second][1])];
        node_ids.push(find_or_insert_node(input, midpoint));
    }

    Ok(node_ids)
}

fn find_or_insert_node(input: &mut Model2DInput, point: [f64; 2]) -> usize {
    if let Some(node) = input.nodes.iter().find(|node| coordinates_match([node.x, node.y], point)) {
        return node.id;
    }

    let id = input.nodes.iter().map(|node| node.id).max().unwrap_or(0) + 1;
    input.nodes.push(Node2DInput { id, x: point[0], y: point[1] });
    id
}

fn coordinates_match(first: [f64; 2], second: [f64; 2]) -> bool {
    (first[0] - second[0]).abs() <= NODE_MERGE_TOLERANCE && (first[1] - second[1]).abs() <= NODE_MERGE_TOLERANCE
}

fn normalize_corner_order(kind: ElementDraftKind, points: &mut [[f64; 2]]) {
    if !kind.is_surface() || signed_polygon_area(points) >= 0.0 {
        return;
    }

    points[1..].reverse();
}

fn signed_polygon_area(points: &[[f64; 2]]) -> f64 {
    points
        .iter()
        .zip(points.iter().cycle().skip(1))
        .take(points.len())
        .map(|(first, second)| first[0] * second[1] - second[0] * first[1])
        .sum::<f64>()
        * 0.5
}

fn parse_usize_field(input: &str, field_name: &str) -> Result<usize, String> {
    input.trim().parse::<usize>().map_err(|_| format!("{field_name} must be an integer ID."))
}

fn parse_f64_field(input: &str, field_name: &str) -> Result<f64, String> {
    let value = input.trim().parse::<f64>().map_err(|_| format!("{field_name} must be a number."))?;

    value.is_finite().then_some(value).ok_or_else(|| format!("{field_name} must be finite."))
}

fn element_from_draft(
    id: usize, kind: ElementDraftKind, node_ids: Vec<usize>, material_id: usize, section_id: usize,
) -> Result<Element2D, String> {
    match kind {
        ElementDraftKind::Truss => Ok(Element2D::Truss(
            Truss2D::new(id, node_array(node_ids, kind)?, material_id, section_id)
                .map_err(|error| error.to_string())?,
        )),
        ElementDraftKind::Beam => Ok(Element2D::Beam(
            Beam2D::new(id, node_array(node_ids, kind)?, material_id, section_id).map_err(|error| error.to_string())?,
        )),
        ElementDraftKind::T3 => Ok(Element2D::TriangleT3(
            TriangleT3::new(id, node_array(node_ids, kind)?, material_id, section_id)
                .map_err(|error| error.to_string())?,
        )),
        ElementDraftKind::T6 => Ok(Element2D::TriangleT6(
            TriangleT6::new(id, node_array(node_ids, kind)?, material_id, section_id)
                .map_err(|error| error.to_string())?,
        )),
        ElementDraftKind::Q4 => Ok(Element2D::QuadQ4(
            QuadQ4::new(id, node_array(node_ids, kind)?, material_id, section_id).map_err(|error| error.to_string())?,
        )),
        ElementDraftKind::Q8 => Ok(Element2D::QuadQ8(
            QuadQ8::new(id, node_array(node_ids, kind)?, material_id, section_id).map_err(|error| error.to_string())?,
        )),
    }
}

fn node_array<const COUNT: usize>(node_ids: Vec<usize>, kind: ElementDraftKind) -> Result<[usize; COUNT], String> {
    node_ids
        .try_into()
        .map_err(|node_ids: Vec<usize>| format!("{} requires {COUNT} node IDs, got {}.", kind.label(), node_ids.len()))
}

#[cfg(test)]
mod tests {
    use super::{
        DEFAULT_BEAM_SECTION_ID, DEFAULT_DRAWING_MATERIAL_ID, DEFAULT_DRAWING_SECTION_ID, DEFAULT_TRUSS_SECTION_ID,
        Element2D, ElementDraftKind, RustyFemGuiApp, element_from_draft, parse_usize_field, residual_convergence,
    };
    use std::collections::HashSet;

    #[test]
    fn accepts_zero_as_an_integer_id() {
        assert_eq!(parse_usize_field("0", "material ID").expect("zero material ID should parse"), 0);
    }

    #[test]
    fn residual_progress_tracks_completed_orders_of_magnitude() {
        assert_eq!(residual_convergence(1.0, 1e-8), 0.0);
        assert_eq!(residual_convergence(1e-4, 1e-8), 0.5);
        assert_eq!(residual_convergence(1e-8, 1e-8), 1.0);
    }

    #[test]
    fn creates_t3_element_from_editor_draft() {
        let element = element_from_draft(10, ElementDraftKind::T3, vec![1, 2, 3], 0, 1)
            .expect("valid T3 draft should create an element");

        assert!(matches!(element, Element2D::TriangleT3(_)));
        assert_eq!(element.id(), 10);
        assert_eq!(element.material_id(), 0);
        assert_eq!(element.section_id(), 1);
    }

    #[test]
    fn clear_canvas_starts_empty_dirty_model_with_drawing_defaults() {
        let mut app = RustyFemGuiApp::new(None);

        app.clear_canvas();

        let model = &app.loaded_model.as_ref().expect("clear canvas should create a model").model;
        assert!(model.nodes().is_empty());
        assert!(model.elements().is_empty());
        assert!(model.constraints().is_empty());
        assert!(model.loads().is_empty());
        assert!(model.element_loads().is_empty());
        assert!(model.material(DEFAULT_DRAWING_MATERIAL_ID).is_ok());
        assert!(model.section(DEFAULT_DRAWING_SECTION_ID).is_ok());
        assert!(model.truss_section(DEFAULT_TRUSS_SECTION_ID).is_ok());
        assert!(model.beam_section(DEFAULT_BEAM_SECTION_ID).is_ok());
        assert!(app.model_dirty);
        assert!(app.selection.is_none());
        assert_eq!(app.next_node_id(), 1);
        assert_eq!(app.next_element_id(), 1);
    }

    #[test]
    fn canvas_workflow_builds_edits_and_solves_a_t3_model() {
        let mut app = RustyFemGuiApp::new(None);
        app.clear_canvas();
        app.place_element_point(0.0, 0.0);
        app.place_element_point(1.0, 0.0);
        app.place_element_point(0.0, 1.0);

        app.add_constraints(1);
        app.constraint_ux = false;
        app.constraint_uy = true;
        app.add_constraints(2);
        app.load_dof = rusty_fem::model::Dof2D::Ux;
        app.load_value_input = "1.0".to_owned();
        app.add_nodal_load(3);
        app.start_analysis(eframe::egui::Context::default());
        let deadline = std::time::Instant::now() + std::time::Duration::from_secs(5);

        while app.analysis_task.is_some() && std::time::Instant::now() < deadline {
            app.poll_analysis();
            std::thread::sleep(std::time::Duration::from_millis(1));
        }

        let model = &app.loaded_model.as_ref().expect("canvas should hold a model").model;
        assert_eq!(model.nodes().len(), 3);
        assert_eq!(model.elements().len(), 1);
        assert_eq!(model.constraints().len(), 3);
        assert_eq!(model.loads().len(), 1);
        assert!(app.analysis_task.is_none());
        assert!(app.results.as_ref().expect("analysis should complete").max_displacement > 0.0);

        app.undo();
        assert!(app.results.is_none());
        assert!(app.loaded_model.as_ref().expect("model should remain loaded").model.loads().is_empty());
        app.redo();
        assert_eq!(app.loaded_model.as_ref().expect("model should remain loaded").model.loads().len(), 1);
    }

    #[test]
    fn t6_placement_creates_midside_nodes_and_undoes_atomically() {
        let mut app = RustyFemGuiApp::new(None);
        app.clear_canvas();
        app.element_kind = ElementDraftKind::T6;

        app.place_element_point(0.0, 0.0);
        app.place_element_point(2.0, 0.0);
        app.place_element_point(0.0, 2.0);

        let model = &app.loaded_model.as_ref().expect("canvas should hold a model").model;
        assert_eq!(model.nodes().len(), 6);
        assert_eq!(model.elements().len(), 1);
        assert_eq!(model.elements()[0].node_ids(), &[1, 2, 3, 4, 5, 6]);
        assert!(model.nodes().iter().any(|node| node.x() == 1.0 && node.y() == 0.0));
        assert!(model.nodes().iter().any(|node| node.x() == 1.0 && node.y() == 1.0));
        assert!(model.nodes().iter().any(|node| node.x() == 0.0 && node.y() == 1.0));

        app.undo();
        let model = &app.loaded_model.as_ref().expect("canvas should hold a model").model;
        assert!(model.nodes().is_empty());
        assert!(model.elements().is_empty());
    }

    #[test]
    fn adjacent_t6_elements_share_corner_and_midside_nodes() {
        let mut app = RustyFemGuiApp::new(None);
        app.clear_canvas();
        app.element_kind = ElementDraftKind::T6;

        for [x, y] in [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]] {
            app.place_element_point(x, y);
        }
        for [x, y] in [[1.0, 0.0], [1.0, 1.0], [0.0, 1.0]] {
            app.place_element_point(x, y);
        }

        let model = &app.loaded_model.as_ref().expect("canvas should hold a model").model;
        assert_eq!(model.nodes().len(), 9);
        assert_eq!(model.elements().len(), 2);

        let first = model.elements()[0].node_ids().iter().copied().collect::<HashSet<_>>();
        let second = model.elements()[1].node_ids().iter().copied().collect::<HashSet<_>>();
        assert_eq!(first.intersection(&second).count(), 3);
    }

    #[test]
    fn q8_placement_uses_four_corners_and_creates_four_midside_nodes() {
        let mut app = RustyFemGuiApp::new(None);
        app.clear_canvas();
        app.element_kind = ElementDraftKind::Q8;

        for [x, y] in [[0.0, 0.0], [0.0, 1.0], [1.0, 1.0], [1.0, 0.0]] {
            app.place_element_point(x, y);
        }

        let model = &app.loaded_model.as_ref().expect("canvas should hold a model").model;
        assert_eq!(model.nodes().len(), 8);
        assert_eq!(model.elements().len(), 1);
        assert_eq!(model.elements()[0].node_ids().len(), 8);
    }

    #[test]
    fn constraint_and_load_tools_reject_unavailable_rotational_dof() {
        let mut app = RustyFemGuiApp::new(None);
        app.clear_canvas();
        for [x, y] in [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]] {
            app.place_element_point(x, y);
        }

        app.constraint_ux = false;
        app.constraint_uy = false;
        app.constraint_rz = true;
        app.add_constraints(1);
        app.load_dof = rusty_fem::model::Dof2D::Rz;
        app.add_nodal_load(1);

        let model = &app.loaded_model.as_ref().expect("canvas should hold a model").model;
        assert!(model.constraints().is_empty());
        assert!(model.loads().is_empty());
    }

    #[test]
    fn constraint_tool_replaces_existing_node_constraints() {
        let mut app = RustyFemGuiApp::new(None);
        app.clear_canvas();
        for [x, y] in [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]] {
            app.place_element_point(x, y);
        }
        app.apply_input_edit("Inserted legacy invalid constraint.", |input| {
            input.constraints.push(rusty_fem::io::DisplacementConstraint2DInput {
                node: 1,
                dof: rusty_fem::io::Dof2DInput::Rz,
                value: 0.0,
            });
            Ok(())
        });

        app.constraint_ux = true;
        app.constraint_uy = true;
        app.constraint_rz = false;
        app.add_constraints(1);

        let constraints = app
            .loaded_model
            .as_ref()
            .expect("canvas should hold a model")
            .model
            .constraints()
            .iter()
            .filter(|constraint| constraint.node_id() == 1)
            .map(|constraint| constraint.dof())
            .collect::<HashSet<_>>();
        assert_eq!(constraints, HashSet::from([rusty_fem::model::Dof2D::Ux, rusty_fem::model::Dof2D::Uy]));
    }
}
