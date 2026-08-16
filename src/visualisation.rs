//! Utilities for exporting simple visualisations of FEM models.

use crate::elements::Element2D;
use crate::model::{Model2D, Node2D};
use std::collections::HashMap;
use std::fmt::Write as _;
use std::fs;
use std::io;
use std::path::Path;

const SCALAR_LEGEND_GAP: f64 = 24.0;
const SCALAR_LEGEND_WIDTH: f64 = 220.0;
const SCALAR_LEGEND_TICK_COUNT: usize = 6;

/// Controls the appearance of an exported 2D mesh SVG.
#[derive(Debug, Clone)]
pub struct SvgMeshOptions {
    /// Total SVG canvas width in pixels.
    pub canvas_width: f64,
    /// Empty space around the model in pixels.
    pub margin: f64,
    /// Draw small markers at model nodes.
    pub show_nodes: bool,
    /// Draw the external and internal free boundaries with a stronger stroke.
    pub show_boundary: bool,
}

/// Controls the appearance of an exported element scalar field SVG.
#[derive(Debug, Clone)]
pub struct SvgElementScalarFieldOptions {
    /// Total SVG canvas width in pixels.
    pub canvas_width: f64,
    /// Empty space around the model in pixels.
    pub margin: f64,
    /// Label shown in the scalar legend.
    pub legend_label: String,
    /// Value mapped to blue.
    pub min_value: f64,
    /// Optional value mapped to red. When unset, the maximum finite field value is used.
    pub max_value: Option<f64>,
    /// Draw element edges over the scalar field.
    pub show_mesh: bool,
    /// Draw the external and internal free boundaries with a stronger stroke.
    pub show_boundary: bool,
}

/// One filled triangular patch used to render a sampled scalar field.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TriangleScalarPatch2D {
    /// Source element ID used only for SVG metadata.
    pub element_id: usize,
    /// Physical patch vertices, in model coordinates.
    pub points: [(f64, f64); 3],
    /// Scalar value represented by this patch.
    pub value: f64,
}

impl Default for SvgMeshOptions {
    fn default() -> Self {
        Self { canvas_width: 1000.0, margin: 32.0, show_nodes: false, show_boundary: true }
    }
}

impl Default for SvgElementScalarFieldOptions {
    fn default() -> Self {
        Self {
            canvas_width: 1000.0,
            margin: 32.0,
            legend_label: "value".to_string(),
            min_value: 0.0,
            max_value: None,
            show_mesh: true,
            show_boundary: true,
        }
    }
}

/// Writes a 2D model mesh to an SVG file using default visualisation options.
pub fn write_model_2d_mesh_svg<P: AsRef<Path>>(model: &Model2D, path: P) -> io::Result<()> {
    write_model_2d_mesh_svg_with_options(model, path, &SvgMeshOptions::default())
}

/// Writes an element scalar field to an SVG file using default visualisation options.
pub fn write_model_2d_element_scalar_svg<P: AsRef<Path>>(
    model: &Model2D, element_values: &[(usize, f64)], path: P,
) -> io::Result<()> {
    write_model_2d_element_scalar_svg_with_options(
        model,
        element_values,
        path,
        &SvgElementScalarFieldOptions::default(),
    )
}

/// Writes a 2D model mesh to an SVG file using custom visualisation options.
pub fn write_model_2d_mesh_svg_with_options<P: AsRef<Path>>(
    model: &Model2D, path: P, options: &SvgMeshOptions,
) -> io::Result<()> {
    let svg = model_2d_mesh_svg(model, options)?;
    let path = path.as_ref();

    if let Some(parent) = path.parent()
        && !parent.as_os_str().is_empty()
    {
        fs::create_dir_all(parent)?;
    }

    fs::write(path, svg)
}

/// Writes an element scalar field to an SVG file using custom visualisation options.
pub fn write_model_2d_element_scalar_svg_with_options<P: AsRef<Path>>(
    model: &Model2D, element_values: &[(usize, f64)], path: P, options: &SvgElementScalarFieldOptions,
) -> io::Result<()> {
    let svg = model_2d_element_scalar_svg(model, element_values, options)?;
    let path = path.as_ref();

    if let Some(parent) = path.parent()
        && !parent.as_os_str().is_empty()
    {
        fs::create_dir_all(parent)?;
    }

    fs::write(path, svg)
}

/// Writes sampled triangular scalar patches to an SVG file using custom visualisation options.
pub fn write_model_2d_triangle_scalar_patches_svg_with_options<P: AsRef<Path>>(
    model: &Model2D, patches: &[TriangleScalarPatch2D], path: P, options: &SvgElementScalarFieldOptions,
) -> io::Result<()> {
    let svg = model_2d_triangle_scalar_patches_svg(model, patches, options)?;
    let path = path.as_ref();

    if let Some(parent) = path.parent()
        && !parent.as_os_str().is_empty()
    {
        fs::create_dir_all(parent)?;
    }

    fs::write(path, svg)
}

/// Renders a 2D model mesh as an SVG string.
pub fn model_2d_mesh_svg(model: &Model2D, options: &SvgMeshOptions) -> io::Result<String> {
    let bounds = Bounds::from_nodes(model.nodes())?;
    let projector = Projector::new(bounds, options)?;
    let node_lookup = node_lookup(model.nodes());
    let boundary_edges = if options.show_boundary { boundary_edges(model.elements()) } else { Vec::new() };
    let mut svg = String::new();

    writeln!(
        svg,
        r#"<svg xmlns="http://www.w3.org/2000/svg" width="{:.0}" height="{:.0}" viewBox="0 0 {:.3} {:.3}" role="img">"#,
        projector.canvas_width, projector.canvas_height, projector.canvas_width, projector.canvas_height
    )
    .expect("writing to String should not fail");
    writeln!(svg, "<title>2D FEM mesh</title>").expect("writing to String should not fail");
    writeln!(
        svg,
        "<desc>Mesh with {} nodes and {} elements. Model bounds: x=[{:.6}, {:.6}], y=[{:.6}, {:.6}].</desc>",
        model.nodes().len(),
        model.elements().len(),
        bounds.min_x,
        bounds.max_x,
        bounds.min_y,
        bounds.max_y
    )
    .expect("writing to String should not fail");
    write_styles(&mut svg);
    writeln!(svg, r#"<rect class="background" x="0" y="0" width="100%" height="100%"/>"#)
        .expect("writing to String should not fail");

    writeln!(svg, r#"<g id="mesh">"#).expect("writing to String should not fail");
    for element in model.elements() {
        let node_ids = element_outline_node_ids(element);
        write_polyline(&mut svg, "mesh-edge", Some(element.id()), &node_ids, &node_lookup, &projector)?;
    }
    writeln!(svg, "</g>").expect("writing to String should not fail");

    if options.show_boundary {
        writeln!(svg, r#"<g id="boundary">"#).expect("writing to String should not fail");
        for edge in boundary_edges {
            write_polyline(&mut svg, "boundary-edge", None, &edge.node_ids, &node_lookup, &projector)?;
        }
        writeln!(svg, "</g>").expect("writing to String should not fail");
    }

    if options.show_nodes {
        writeln!(svg, r#"<g id="nodes">"#).expect("writing to String should not fail");
        for node in model.nodes() {
            let (x, y) = projector.project_node(node);
            writeln!(svg, r#"<circle class="node" cx="{x:.3}" cy="{y:.3}" r="1.8" data-node-id="{}"/>"#, node.id())
                .expect("writing to String should not fail");
        }
        writeln!(svg, "</g>").expect("writing to String should not fail");
    }

    writeln!(svg, "</svg>").expect("writing to String should not fail");

    Ok(svg)
}

/// Renders an element scalar field as an SVG string.
pub fn model_2d_element_scalar_svg(
    model: &Model2D, element_values: &[(usize, f64)], options: &SvgElementScalarFieldOptions,
) -> io::Result<String> {
    let bounds = Bounds::from_nodes(model.nodes())?;
    let mesh_options = SvgMeshOptions {
        canvas_width: options.canvas_width,
        margin: options.margin,
        show_nodes: false,
        show_boundary: options.show_boundary,
    };
    let projector = Projector::new(bounds, &mesh_options)?;
    let node_lookup = node_lookup(model.nodes());
    let value_lookup = scalar_value_lookup(element_values)?;
    let max_value = scalar_max_value(element_values, options)?;
    let boundary_edges = if options.show_boundary { boundary_edges(model.elements()) } else { Vec::new() };
    let canvas_width = projector.canvas_width + SCALAR_LEGEND_GAP + SCALAR_LEGEND_WIDTH + options.margin;
    let mut svg = String::new();

    writeln!(
        svg,
        r#"<svg xmlns="http://www.w3.org/2000/svg" width="{:.0}" height="{:.0}" viewBox="0 0 {:.3} {:.3}" role="img">"#,
        canvas_width, projector.canvas_height, canvas_width, projector.canvas_height
    )
    .expect("writing to String should not fail");
    writeln!(svg, "<title>{}</title>", escape_xml(&options.legend_label)).expect("writing to String should not fail");
    writeln!(
        svg,
        "<desc>Element scalar field for {} nodes and {} elements. Values are clamped to [{:.6}, {:.6}].</desc>",
        model.nodes().len(),
        model.elements().len(),
        options.min_value,
        max_value
    )
    .expect("writing to String should not fail");
    write_styles(&mut svg);
    writeln!(svg, r#"<rect class="background" x="0" y="0" width="100%" height="100%"/>"#)
        .expect("writing to String should not fail");

    writeln!(svg, r#"<g id="scalar-field">"#).expect("writing to String should not fail");
    for element in model.elements() {
        let value = value_lookup.get(&element.id()).copied().unwrap_or(options.min_value);
        let color = scalar_color(value, options.min_value, max_value);
        let node_ids = element_outline_node_ids(element);

        write_polygon(&mut svg, "scalar-element", element.id(), &node_ids, &node_lookup, &projector, &color)?;
    }
    writeln!(svg, "</g>").expect("writing to String should not fail");

    if options.show_mesh {
        writeln!(svg, r#"<g id="mesh">"#).expect("writing to String should not fail");
        for element in model.elements() {
            let node_ids = element_outline_node_ids(element);
            write_polyline(
                &mut svg,
                "mesh-edge scalar-mesh-edge",
                Some(element.id()),
                &node_ids,
                &node_lookup,
                &projector,
            )?;
        }
        writeln!(svg, "</g>").expect("writing to String should not fail");
    }

    if options.show_boundary {
        writeln!(svg, r#"<g id="boundary">"#).expect("writing to String should not fail");
        for edge in boundary_edges {
            write_polyline(&mut svg, "boundary-edge", None, &edge.node_ids, &node_lookup, &projector)?;
        }
        writeln!(svg, "</g>").expect("writing to String should not fail");
    }

    write_scalar_legend(&mut svg, options, max_value, &projector);
    writeln!(svg, "</svg>").expect("writing to String should not fail");

    Ok(svg)
}

/// Renders sampled triangular scalar patches as an SVG string.
pub fn model_2d_triangle_scalar_patches_svg(
    model: &Model2D, patches: &[TriangleScalarPatch2D], options: &SvgElementScalarFieldOptions,
) -> io::Result<String> {
    let bounds = Bounds::from_nodes(model.nodes())?;
    let mesh_options = SvgMeshOptions {
        canvas_width: options.canvas_width,
        margin: options.margin,
        show_nodes: false,
        show_boundary: options.show_boundary,
    };
    let projector = Projector::new(bounds, &mesh_options)?;
    let node_lookup = node_lookup(model.nodes());
    let max_value = scalar_patch_max_value(patches, options)?;
    let boundary_edges = if options.show_boundary { boundary_edges(model.elements()) } else { Vec::new() };
    let canvas_width = projector.canvas_width + SCALAR_LEGEND_GAP + SCALAR_LEGEND_WIDTH + options.margin;
    let mut svg = String::new();

    writeln!(
        svg,
        r#"<svg xmlns="http://www.w3.org/2000/svg" width="{:.0}" height="{:.0}" viewBox="0 0 {:.3} {:.3}" role="img">"#,
        canvas_width, projector.canvas_height, canvas_width, projector.canvas_height
    )
    .expect("writing to String should not fail");
    writeln!(svg, "<title>{}</title>", escape_xml(&options.legend_label)).expect("writing to String should not fail");
    writeln!(
        svg,
        "<desc>Sampled triangular scalar field for {} nodes and {} elements. Values are clamped to [{:.6}, {:.6}].</desc>",
        model.nodes().len(),
        model.elements().len(),
        options.min_value,
        max_value
    )
    .expect("writing to String should not fail");
    write_styles(&mut svg);
    writeln!(svg, r#"<rect class="background" x="0" y="0" width="100%" height="100%"/>"#)
        .expect("writing to String should not fail");

    writeln!(svg, r#"<g id="scalar-field">"#).expect("writing to String should not fail");
    for (patch_index, patch) in patches.iter().enumerate() {
        let color = scalar_color(patch.value, options.min_value, max_value);

        write_triangle_patch_polygon(&mut svg, patch_index, patch, &projector, &color)?;
    }
    writeln!(svg, "</g>").expect("writing to String should not fail");

    if options.show_mesh {
        writeln!(svg, r#"<g id="mesh">"#).expect("writing to String should not fail");
        for element in model.elements() {
            let node_ids = element_outline_node_ids(element);
            write_polyline(
                &mut svg,
                "mesh-edge scalar-mesh-edge",
                Some(element.id()),
                &node_ids,
                &node_lookup,
                &projector,
            )?;
        }
        writeln!(svg, "</g>").expect("writing to String should not fail");
    }

    if options.show_boundary {
        writeln!(svg, r#"<g id="boundary">"#).expect("writing to String should not fail");
        for edge in boundary_edges {
            write_polyline(&mut svg, "boundary-edge", None, &edge.node_ids, &node_lookup, &projector)?;
        }
        writeln!(svg, "</g>").expect("writing to String should not fail");
    }

    write_scalar_legend(&mut svg, options, max_value, &projector);
    writeln!(svg, "</svg>").expect("writing to String should not fail");

    Ok(svg)
}

#[derive(Clone, Copy)]
struct Bounds {
    min_x: f64,
    max_x: f64,
    min_y: f64,
    max_y: f64,
}

impl Bounds {
    fn from_nodes(nodes: &[Node2D]) -> io::Result<Self> {
        let first = nodes.first().ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "model has no nodes"))?;
        let mut bounds = Self { min_x: first.x(), max_x: first.x(), min_y: first.y(), max_y: first.y() };

        for node in nodes.iter().skip(1) {
            bounds.min_x = bounds.min_x.min(node.x());
            bounds.max_x = bounds.max_x.max(node.x());
            bounds.min_y = bounds.min_y.min(node.y());
            bounds.max_y = bounds.max_y.max(node.y());
        }

        Ok(bounds)
    }

    fn width(self) -> f64 {
        (self.max_x - self.min_x).max(1.0)
    }

    fn height(self) -> f64 {
        (self.max_y - self.min_y).max(1.0)
    }
}

struct Projector {
    bounds: Bounds,
    scale: f64,
    margin: f64,
    canvas_width: f64,
    canvas_height: f64,
}

impl Projector {
    fn new(bounds: Bounds, options: &SvgMeshOptions) -> io::Result<Self> {
        if !options.canvas_width.is_finite() || options.canvas_width <= 0.0 {
            return Err(io::Error::new(io::ErrorKind::InvalidInput, "SVG canvas width must be positive and finite"));
        }

        if !options.margin.is_finite() || options.margin < 0.0 || 2.0 * options.margin >= options.canvas_width {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "SVG margin must be finite, non-negative, and smaller than half of the canvas width",
            ));
        }

        let drawing_width = options.canvas_width - 2.0 * options.margin;
        let scale = drawing_width / bounds.width();
        let canvas_height = bounds.height() * scale + 2.0 * options.margin;

        Ok(Self { bounds, scale, margin: options.margin, canvas_width: options.canvas_width, canvas_height })
    }

    fn project_node(&self, node: &Node2D) -> (f64, f64) {
        self.project_point(node.x(), node.y())
    }

    fn project_point(&self, x: f64, y: f64) -> (f64, f64) {
        (self.margin + (x - self.bounds.min_x) * self.scale, self.margin + (self.bounds.max_y - y) * self.scale)
    }
}

#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
struct EdgeKey {
    first: usize,
    second: usize,
}

impl EdgeKey {
    fn new(first: usize, second: usize) -> Self {
        if first <= second { Self { first, second } } else { Self { first: second, second: first } }
    }
}

#[derive(Clone, Debug)]
struct BoundaryEdge {
    node_ids: Vec<usize>,
}

fn node_lookup(nodes: &[Node2D]) -> HashMap<usize, &Node2D> {
    nodes.iter().map(|node| (node.id(), node)).collect()
}

fn element_outline_node_ids(element: &Element2D) -> Vec<usize> {
    let node_ids = element.node_ids();

    match element {
        Element2D::Truss(_) | Element2D::Beam(_) => vec![node_ids[0], node_ids[1]],
        Element2D::TriangleT3(_) => vec![node_ids[0], node_ids[1], node_ids[2], node_ids[0]],
        Element2D::TriangleT6(_) => {
            vec![node_ids[0], node_ids[3], node_ids[1], node_ids[4], node_ids[2], node_ids[5], node_ids[0]]
        }
        Element2D::QuadQ4(_) => vec![node_ids[0], node_ids[1], node_ids[2], node_ids[3], node_ids[0]],
        Element2D::QuadQ8(_) => {
            vec![
                node_ids[0],
                node_ids[4],
                node_ids[1],
                node_ids[5],
                node_ids[2],
                node_ids[6],
                node_ids[3],
                node_ids[7],
                node_ids[0],
            ]
        }
    }
}

fn boundary_edges(elements: &[Element2D]) -> Vec<BoundaryEdge> {
    let mut edges = HashMap::<EdgeKey, (usize, Vec<usize>)>::new();

    for element in elements {
        for node_ids in element_boundary_edge_node_ids(element) {
            let key = EdgeKey::new(node_ids[0], *node_ids.last().expect("boundary edge should not be empty"));
            let entry = edges.entry(key).or_insert_with(|| (0, node_ids.clone()));

            entry.0 += 1;
        }
    }

    edges.into_values().filter_map(|(count, node_ids)| (count == 1).then_some(BoundaryEdge { node_ids })).collect()
}

fn element_boundary_edge_node_ids(element: &Element2D) -> Vec<Vec<usize>> {
    let node_ids = element.node_ids();

    match element {
        Element2D::Truss(_) | Element2D::Beam(_) => Vec::new(),
        Element2D::TriangleT3(_) => {
            vec![vec![node_ids[0], node_ids[1]], vec![node_ids[1], node_ids[2]], vec![node_ids[2], node_ids[0]]]
        }
        Element2D::TriangleT6(_) => {
            vec![
                vec![node_ids[0], node_ids[3], node_ids[1]],
                vec![node_ids[1], node_ids[4], node_ids[2]],
                vec![node_ids[2], node_ids[5], node_ids[0]],
            ]
        }
        Element2D::QuadQ4(_) => {
            vec![
                vec![node_ids[0], node_ids[1]],
                vec![node_ids[1], node_ids[2]],
                vec![node_ids[2], node_ids[3]],
                vec![node_ids[3], node_ids[0]],
            ]
        }
        Element2D::QuadQ8(_) => {
            vec![
                vec![node_ids[0], node_ids[4], node_ids[1]],
                vec![node_ids[1], node_ids[5], node_ids[2]],
                vec![node_ids[2], node_ids[6], node_ids[3]],
                vec![node_ids[3], node_ids[7], node_ids[0]],
            ]
        }
    }
}

fn write_styles(svg: &mut String) {
    writeln!(
        svg,
        r#"<style>
  .background {{ fill: #ffffff; }}
  .mesh-edge {{ fill: none; stroke: #9ca3af; stroke-width: 0.7; vector-effect: non-scaling-stroke; }}
  .scalar-element {{ stroke: none; }}
  .scalar-patch {{ stroke: none; shape-rendering: geometricPrecision; }}
  .scalar-mesh-edge {{ stroke: rgba(17, 24, 39, 0.28); stroke-width: 0.45; }}
  .boundary-edge {{ fill: none; stroke: #111827; stroke-width: 1.8; vector-effect: non-scaling-stroke; stroke-linecap: round; stroke-linejoin: round; }}
  .node {{ fill: #2563eb; opacity: 0.75; }}
  .legend-background {{ fill: #f8fafc; stroke: #cbd5e1; stroke-width: 1; }}
  .legend-title {{ fill: #0f172a; font-family: system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif; font-size: 14px; font-weight: 650; }}
  .legend-text {{ fill: #334155; font-family: system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif; font-size: 12px; }}
  .legend-tick {{ stroke: #64748b; stroke-width: 0.8; vector-effect: non-scaling-stroke; }}
  .legend-axis {{ stroke: #94a3b8; stroke-width: 0.8; vector-effect: non-scaling-stroke; }}
</style>"#
    )
    .expect("writing to String should not fail");
}

fn scalar_value_lookup(element_values: &[(usize, f64)]) -> io::Result<HashMap<usize, f64>> {
    let mut values = HashMap::with_capacity(element_values.len());

    for &(element_id, value) in element_values {
        if !value.is_finite() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("scalar value for element {element_id} is not finite"),
            ));
        }

        values.insert(element_id, value);
    }

    Ok(values)
}

fn scalar_max_value(element_values: &[(usize, f64)], options: &SvgElementScalarFieldOptions) -> io::Result<f64> {
    if !options.min_value.is_finite() {
        return Err(io::Error::new(io::ErrorKind::InvalidInput, "scalar minimum must be finite"));
    }

    let max_value = if let Some(max_value) = options.max_value {
        max_value
    } else {
        element_values.iter().map(|(_, value)| *value).fold(options.min_value, f64::max)
    };

    if !max_value.is_finite() || max_value <= options.min_value {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "scalar maximum must be finite and greater than the scalar minimum",
        ));
    }

    Ok(max_value)
}

fn scalar_patch_max_value(
    patches: &[TriangleScalarPatch2D], options: &SvgElementScalarFieldOptions,
) -> io::Result<f64> {
    if !options.min_value.is_finite() {
        return Err(io::Error::new(io::ErrorKind::InvalidInput, "scalar minimum must be finite"));
    }

    for patch in patches {
        if !patch.value.is_finite() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("scalar patch value for element {} is not finite", patch.element_id),
            ));
        }

        for &(x, y) in &patch.points {
            if !x.is_finite() || !y.is_finite() {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!("scalar patch point for element {} is not finite", patch.element_id),
                ));
            }
        }
    }

    let max_value = if let Some(max_value) = options.max_value {
        max_value
    } else {
        patches.iter().map(|patch| patch.value).fold(options.min_value, f64::max)
    };

    if !max_value.is_finite() || max_value <= options.min_value {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "scalar maximum must be finite and greater than the scalar minimum",
        ));
    }

    Ok(max_value)
}

fn scalar_color(value: f64, min_value: f64, max_value: f64) -> String {
    let fraction = ((value - min_value) / (max_value - min_value)).clamp(0.0, 1.0);
    let [red, green, blue] = scalar_color_rgb(fraction);

    format!("#{red:02x}{green:02x}{blue:02x}")
}

/// Returns the RGB color used by RustyFEM scalar contours.
///
/// `fraction` is clamped to the range from the blue minimum to the red
/// maximum. The GUI and SVG renderer share this function so their result
/// fields use an identical palette.
pub fn scalar_color_rgb(fraction: f64) -> [u8; 3] {
    let fraction = fraction.clamp(0.0, 1.0);
    let color_stops = [
        (0.0, [0, 0, 255]),
        (0.2, [0, 180, 255]),
        (0.4, [0, 190, 80]),
        (0.6, [255, 240, 0]),
        (0.8, [255, 145, 0]),
        (1.0, [255, 0, 0]),
    ];
    let (start_fraction, start_color, end_fraction, end_color) = color_stops
        .windows(2)
        .find_map(|window| {
            let (start_fraction, start_color) = window[0];
            let (end_fraction, end_color) = window[1];

            (fraction <= end_fraction).then_some((start_fraction, start_color, end_fraction, end_color))
        })
        .expect("color stops should cover the full scalar range");
    let local_fraction = (fraction - start_fraction) / (end_fraction - start_fraction);
    let red = interpolate_color_channel(start_color[0], end_color[0], local_fraction);
    let green = interpolate_color_channel(start_color[1], end_color[1], local_fraction);
    let blue = interpolate_color_channel(start_color[2], end_color[2], local_fraction);

    [red, green, blue]
}

fn interpolate_color_channel(start: u8, end: u8, fraction: f64) -> u8 {
    (start as f64 + fraction * (end as f64 - start as f64)).round() as u8
}

fn write_polygon(
    svg: &mut String, class_name: &str, element_id: usize, node_ids: &[usize], node_lookup: &HashMap<usize, &Node2D>,
    projector: &Projector, fill: &str,
) -> io::Result<()> {
    write!(svg, r#"<polygon class="{class_name}" data-element-id="{element_id}" fill="{fill}" points=""#)
        .expect("writing to String should not fail");

    for (index, node_id) in node_ids.iter().copied().enumerate() {
        let node = node_lookup.get(&node_id).ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidData, format!("element references unknown node ID {node_id}"))
        })?;
        let (x, y) = projector.project_node(node);

        if index > 0 {
            svg.push(' ');
        }

        write!(svg, "{x:.3},{y:.3}").expect("writing to String should not fail");
    }

    writeln!(svg, r#""/>"#).expect("writing to String should not fail");

    Ok(())
}

fn write_triangle_patch_polygon(
    svg: &mut String, patch_index: usize, patch: &TriangleScalarPatch2D, projector: &Projector, fill: &str,
) -> io::Result<()> {
    write!(
        svg,
        r#"<polygon class="scalar-patch" data-element-id="{}" data-patch-index="{patch_index}" fill="{fill}" points=""#,
        patch.element_id
    )
    .expect("writing to String should not fail");

    for (index, &(physical_x, physical_y)) in patch.points.iter().enumerate() {
        let (x, y) = projector.project_point(physical_x, physical_y);

        if index > 0 {
            svg.push(' ');
        }

        write!(svg, "{x:.3},{y:.3}").expect("writing to String should not fail");
    }

    writeln!(svg, r#""/>"#).expect("writing to String should not fail");

    Ok(())
}

fn write_polyline(
    svg: &mut String, class_name: &str, element_id: Option<usize>, node_ids: &[usize],
    node_lookup: &HashMap<usize, &Node2D>, projector: &Projector,
) -> io::Result<()> {
    write!(svg, r#"<polyline class="{class_name}""#).expect("writing to String should not fail");

    if let Some(element_id) = element_id {
        write!(svg, r#" data-element-id="{element_id}""#).expect("writing to String should not fail");
    }

    write!(svg, r#" points=""#).expect("writing to String should not fail");

    for (index, node_id) in node_ids.iter().copied().enumerate() {
        let node = node_lookup.get(&node_id).ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidData, format!("element references unknown node ID {node_id}"))
        })?;
        let (x, y) = projector.project_node(node);

        if index > 0 {
            svg.push(' ');
        }

        write!(svg, "{x:.3},{y:.3}").expect("writing to String should not fail");
    }

    writeln!(svg, r#""/>"#).expect("writing to String should not fail");

    Ok(())
}

fn write_scalar_legend(
    svg: &mut String, options: &SvgElementScalarFieldOptions, max_value: f64, projector: &Projector,
) {
    let x = projector.canvas_width + SCALAR_LEGEND_GAP;
    let y = options.margin;
    let width = SCALAR_LEGEND_WIDTH;
    let height = (projector.canvas_height - 2.0 * options.margin).max(220.0);
    let ramp_x = x + 24.0;
    let ramp_y = y + 54.0;
    let ramp_width = 28.0;
    let ramp_height = height - 96.0;
    let tick_x = ramp_x + ramp_width;
    let label_x = tick_x + 12.0;
    let steps = 80;

    writeln!(svg, r#"<g id="legend">"#).expect("writing to String should not fail");
    writeln!(
        svg,
        r#"<rect class="legend-background" x="{x:.3}" y="{y:.3}" width="{width:.3}" height="{height:.3}" rx="6"/>"#
    )
    .expect("writing to String should not fail");
    writeln!(
        svg,
        r#"<text class="legend-title" x="{:.3}" y="{:.3}">{}</text>"#,
        x + 16.0,
        y + 25.0,
        escape_xml(&options.legend_label)
    )
    .expect("writing to String should not fail");

    for step in 0..steps {
        let fraction = 1.0 - step as f64 / (steps - 1) as f64;
        let color =
            scalar_color(options.min_value + fraction * (max_value - options.min_value), options.min_value, max_value);
        let rect_y = ramp_y + step as f64 * ramp_height / steps as f64;
        let rect_height = ramp_height / steps as f64 + 0.5;

        writeln!(
            svg,
            r#"<rect x="{ramp_x:.3}" y="{rect_y:.3}" width="{ramp_width:.3}" height="{rect_height:.3}" fill="{color}"/>"#
        )
        .expect("writing to String should not fail");
    }

    writeln!(
        svg,
        r#"<line class="legend-axis" x1="{tick_x:.3}" y1="{ramp_y:.3}" x2="{tick_x:.3}" y2="{:.3}"/>"#,
        ramp_y + ramp_height
    )
    .expect("writing to String should not fail");

    for tick in 0..SCALAR_LEGEND_TICK_COUNT {
        let fraction = tick as f64 / (SCALAR_LEGEND_TICK_COUNT - 1) as f64;
        let y_position = ramp_y + fraction * ramp_height;
        let value = options.min_value + (1.0 - fraction) * (max_value - options.min_value);

        writeln!(
            svg,
            r#"<line class="legend-tick" x1="{:.3}" y1="{y_position:.3}" x2="{:.3}" y2="{y_position:.3}"/>"#,
            tick_x,
            tick_x + 7.0
        )
        .expect("writing to String should not fail");
        writeln!(
            svg,
            r#"<text class="legend-text" x="{label_x:.3}" y="{:.3}">{}</text>"#,
            y_position + 4.0,
            escape_xml(&format_scalar_tick(value))
        )
        .expect("writing to String should not fail");
    }

    writeln!(svg, "</g>").expect("writing to String should not fail");
}

fn format_scalar_tick(value: f64) -> String {
    let absolute = value.abs();

    if absolute == 0.0 {
        "0".to_string()
    } else if !(0.01..10_000.0).contains(&absolute) {
        format!("{value:.3e}")
    } else if absolute < 1.0 {
        format!("{value:.4}")
    } else if absolute < 100.0 {
        format!("{value:.2}")
    } else {
        format!("{value:.1}")
    }
}

fn escape_xml(text: &str) -> String {
    text.replace('&', "&amp;").replace('<', "&lt;").replace('>', "&gt;")
}

#[cfg(test)]
mod tests {
    use super::{
        SvgElementScalarFieldOptions, SvgMeshOptions, TriangleScalarPatch2D, element_outline_node_ids,
        model_2d_element_scalar_svg, model_2d_mesh_svg, model_2d_triangle_scalar_patches_svg, scalar_color,
    };
    use crate::elements::{Element2D, QuadQ8};
    use crate::model::{DEFAULT_MATERIAL_ID, Material2D, Model2D, Node2D, PlaneStressSection2D, Section2D};

    #[test]
    fn q8_outline_passes_through_midside_nodes() {
        let element = Element2D::QuadQ8(QuadQ8::new(1, [1, 2, 3, 4, 5, 6, 7, 8], DEFAULT_MATERIAL_ID, 1).unwrap());

        assert_eq!(element_outline_node_ids(&element), vec![1, 5, 2, 6, 3, 7, 4, 8, 1]);
    }

    #[test]
    fn renders_q8_mesh_svg() {
        let model = sample_q8_model();

        let svg = model_2d_mesh_svg(&model, &SvgMeshOptions::default()).unwrap();

        assert!(svg.contains(r#"<svg "#));
        assert!(svg.contains(r#"data-element-id="1""#));
        assert!(svg.contains(r#"<g id="boundary">"#));
    }

    #[test]
    fn renders_element_scalar_field_svg() {
        let model = sample_q8_model();
        let options = SvgElementScalarFieldOptions {
            legend_label: "sigma_y [MPa]".to_string(),
            max_value: Some(2.0),
            ..Default::default()
        };

        let svg = model_2d_element_scalar_svg(&model, &[(1, 2.0)], &options).unwrap();

        assert!(svg.contains(r#"<g id="scalar-field">"#));
        assert!(svg.contains(r##"fill="#ff0000""##));
        assert!(svg.contains("sigma_y [MPa]"));
    }

    #[test]
    fn renders_triangle_scalar_patches_svg() {
        let model = sample_q8_model();
        let options = SvgElementScalarFieldOptions {
            legend_label: "von Mises [MPa]".to_string(),
            max_value: Some(2.0),
            ..Default::default()
        };
        let patches =
            [TriangleScalarPatch2D { element_id: 1, points: [(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)], value: 2.0 }];

        let svg = model_2d_triangle_scalar_patches_svg(&model, &patches, &options).unwrap();

        assert!(svg.contains(r#"<g id="scalar-field">"#));
        assert!(svg.contains(r#"class="scalar-patch""#));
        assert!(svg.contains(r#"data-patch-index="0""#));
        assert!(svg.contains(r##"fill="#ff0000""##));
        assert!(svg.contains("von Mises [MPa]"));
    }

    #[test]
    fn scalar_color_maps_min_to_blue_and_max_to_red() {
        assert_eq!(scalar_color(0.0, 0.0, 10.0), "#0000ff");
        assert_eq!(scalar_color(5.0, 0.0, 10.0), "#80d728");
        assert_eq!(scalar_color(10.0, 0.0, 10.0), "#ff0000");
    }

    fn sample_q8_model() -> Model2D {
        let mut model = Model2D::new();

        model.set_material(Material2D::new(1.0, 0.3, 1.0).unwrap());
        model.add_section(1, Section2D::PlaneStress(PlaneStressSection2D::new(1.0).unwrap())).unwrap();

        for (id, x, y) in [
            (1, 0.0, 0.0),
            (2, 1.0, 0.0),
            (3, 1.0, 1.0),
            (4, 0.0, 1.0),
            (5, 0.5, 0.0),
            (6, 1.0, 0.5),
            (7, 0.5, 1.0),
            (8, 0.0, 0.5),
        ] {
            model.add_node(Node2D::new(id, x, y).unwrap()).unwrap();
        }

        model
            .add_element(Element2D::QuadQ8(QuadQ8::new(1, [1, 2, 3, 4, 5, 6, 7, 8], DEFAULT_MATERIAL_ID, 1).unwrap()))
            .unwrap();

        model
    }
}
