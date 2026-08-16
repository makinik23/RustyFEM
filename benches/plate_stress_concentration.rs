//! One-shot plate stress-concentration benchmark.
//! Run the T3 benchmark with: `cargo bench --bench plate_stress_concentration -- --t3`.
//! Generate the T3 mesh SVG with:
//! `cargo bench --bench plate_stress_concentration -- --t3-mesh-only`.

use rusty_fem::analysis::iterative_solver::CgOptions;
use rusty_fem::analysis::solver::{solve_sparse_with_options, solve_with_settings};
use rusty_fem::analysis::stress_recovery::{
    recover_nodal_plane_stress_responses, recover_quad_q8_response, recover_triangle_response,
    recover_triangle_t6_response,
};
use rusty_fem::elements::{
    Element2D, QuadQ8, TriangleT3, TriangleT6, quad_q8_shape_function_derivatives,
    triangle_t6_shape_function_derivatives, triangle_t6_shape_functions,
};
use rusty_fem::io::{
    AnalysisSettings2DInput, BeamSection2DInput, DisplacementConstraint2DInput, Dof2DInput, ElementLoad2DInput,
    ElementLoad2DInputKind, ElementType2DInput, LoadCoordinateSystem2DInput, Material2DInput, Model2DInput,
    NodalLoad2DInput, Node2DInput, PlaneStressSection2DInput, Section2DInput, Section2DInputKind, SolverKind2DInput,
    TrussSection2DInput,
};
use rusty_fem::model::{
    DEFAULT_MATERIAL_ID, DisplacementConstraint2D, Dof2D, DofNumbering2D, EdgeTraction2D, ElementLoad2D,
    LoadCoordinateSystem2D, Material2D, Model2D, Node2D, PlaneStressSection2D, Section2D,
};
use rusty_fem::visualisation::{
    SvgElementScalarFieldOptions, TriangleScalarPatch2D, write_model_2d_element_scalar_svg_with_options,
    write_model_2d_mesh_svg, write_model_2d_triangle_scalar_patches_svg_with_options,
};
use std::collections::HashMap;
use std::fmt::Write as _;
use std::fs;
use std::hint::black_box;
use std::time::Instant;

const WIDTH: f64 = 484.0;
const HEIGHT: f64 = 860.0;
const THICKNESS: f64 = 3.0;
const YOUNG_MODULUS: f64 = 70_000.0;
const POISSON_RATIO: f64 = 0.34;
const NOTCH_DEPTH: f64 = 60.0;
const LEFT_NOTCH_RADIUS: f64 = 25.0;
const RIGHT_NOTCH_RADIUS: f64 = 50.0;
const REFERENCE_MASS: f64 = 1_500.0;
const GRAVITY: f64 = 9.81;
const EXPERIMENTAL_LEFT_SIGMA_Y: f64 = 42.95;
const EXPERIMENTAL_RIGHT_SIGMA_Y: f64 = 35.12;
const EXPERIMENTAL_LEFT_ALPHA: f64 = 3.19;
const EXPERIMENTAL_RIGHT_ALPHA: f64 = 2.61;
const MESH_SVG_PATH: &str = "target/plate_stress_concentration_mesh.svg";
const SIGMA_X_SVG_PATH: &str = "target/plate_stress_concentration_sigma_x.svg";
const SIGMA_Y_SVG_PATH: &str = "target/plate_stress_concentration_sigma_y.svg";
const T3_JSON_EXAMPLE_PATH: &str = "examples/plate_stress_concentration_t3.json";
const T6_JSON_EXAMPLE_PATH: &str = "examples/plate_stress_concentration_t6.json";
const T6_IMPROVED_JSON_EXAMPLE_PATH: &str = "examples/plate_stress_concentration_t6_improved_mesh.json";
const T6_DENSE_JSON_EXAMPLE_PATH: &str = "examples/plate_stress_concentration_t6_dense.json";
const T3_MESH_SVG_PATH: &str = "target/plate_stress_concentration/plate_stress_concentration_t3_mesh.svg";
const T3_SIGMA_X_SVG_PATH: &str = "target/plate_stress_concentration/plate_stress_concentration_t3_sigma_x.svg";
const T3_SIGMA_Y_SVG_PATH: &str = "target/plate_stress_concentration/plate_stress_concentration_t3_sigma_y.svg";
const T3_VON_MISES_SVG_PATH: &str = "target/plate_stress_concentration/plate_stress_concentration_t3_von_mises.svg";
const T6_VON_MISES_SVG_PATH: &str = "target/plate_stress_concentration/plate_stress_concentration_t6_von_mises.svg";
const T6_VON_MISES_VISUAL_SUBDIVISIONS: usize = 4;
const T3_DISPLACEMENT_SVG_PATH: &str =
    "target/plate_stress_concentration/plate_stress_concentration_t3_displacement.svg";
const T3_STRAIN_SVG_PATH: &str = "target/plate_stress_concentration/plate_stress_concentration_t3_strain.svg";
const T3_STRAIN_GAUGE_SAMPLING_CSV_PATH: &str =
    "target/plate_stress_concentration/plate_stress_concentration_t3_strain_gauge_sampling.csv";
const T3_VISUAL_ANGLE_COUNT: usize = 128;
const T3_VISUAL_RADIAL_FRACTIONS: [f64; 13] =
    [0.0, 0.012, 0.027, 0.047, 0.075, 0.115, 0.175, 0.26, 0.38, 0.54, 0.73, 0.90, 1.0];
const STRAIN_GAUGE_SAMPLES: [StrainGaugeSampleDefinition; 14] = [
    StrainGaugeSampleDefinition {
        notch: "A_R25",
        distance: 4.0,
        epsilon_y_per_mille: 0.487,
        epsilon_x_per_mille: -0.104,
    },
    StrainGaugeSampleDefinition {
        notch: "A_R25",
        distance: 14.0,
        epsilon_y_per_mille: 0.326,
        epsilon_x_per_mille: -0.017,
    },
    StrainGaugeSampleDefinition {
        notch: "A_R25",
        distance: 24.0,
        epsilon_y_per_mille: 0.248,
        epsilon_x_per_mille: -0.006,
    },
    StrainGaugeSampleDefinition {
        notch: "A_R25",
        distance: 44.0,
        epsilon_y_per_mille: 0.205,
        epsilon_x_per_mille: -0.017,
    },
    StrainGaugeSampleDefinition {
        notch: "A_R25",
        distance: 80.0,
        epsilon_y_per_mille: 0.188,
        epsilon_x_per_mille: -0.035,
    },
    StrainGaugeSampleDefinition {
        notch: "A_R25",
        distance: 134.0,
        epsilon_y_per_mille: 0.173,
        epsilon_x_per_mille: -0.052,
    },
    StrainGaugeSampleDefinition {
        notch: "A_R25",
        distance: 174.0,
        epsilon_y_per_mille: 0.162,
        epsilon_x_per_mille: -0.030,
    },
    StrainGaugeSampleDefinition {
        notch: "B_R50",
        distance: 4.0,
        epsilon_y_per_mille: 0.437,
        epsilon_x_per_mille: -0.115,
    },
    StrainGaugeSampleDefinition {
        notch: "B_R50",
        distance: 14.0,
        epsilon_y_per_mille: 0.328,
        epsilon_x_per_mille: -0.047,
    },
    StrainGaugeSampleDefinition {
        notch: "B_R50",
        distance: 24.0,
        epsilon_y_per_mille: 0.263,
        epsilon_x_per_mille: -0.028,
    },
    StrainGaugeSampleDefinition {
        notch: "B_R50",
        distance: 44.0,
        epsilon_y_per_mille: 0.221,
        epsilon_x_per_mille: -0.012,
    },
    StrainGaugeSampleDefinition {
        notch: "B_R50",
        distance: 80.0,
        epsilon_y_per_mille: 0.194,
        epsilon_x_per_mille: -0.035,
    },
    StrainGaugeSampleDefinition {
        notch: "B_R50",
        distance: 134.0,
        epsilon_y_per_mille: 0.181,
        epsilon_x_per_mille: -0.048,
    },
    StrainGaugeSampleDefinition {
        notch: "B_R50",
        distance: 174.0,
        epsilon_y_per_mille: 0.140,
        epsilon_x_per_mille: -0.049,
    },
];

#[derive(Clone, Copy)]
struct PlateMeshSpec {
    right_upper_elements: usize,
    top_elements: usize,
    left_upper_elements: usize,
    left_lower_elements: usize,
    bottom_elements: usize,
    right_lower_elements: usize,
    radial_elements: usize,
    radial_bias: f64,
}

impl PlateMeshSpec {
    fn smoke_q8() -> Self {
        Self {
            right_upper_elements: 2,
            top_elements: 6,
            left_upper_elements: 2,
            left_lower_elements: 2,
            bottom_elements: 6,
            right_lower_elements: 2,
            radial_elements: 2,
            radial_bias: 1.5,
        }
    }

    fn reference_q8() -> Self {
        Self {
            right_upper_elements: 18,
            top_elements: 40,
            left_upper_elements: 22,
            left_lower_elements: 22,
            bottom_elements: 40,
            right_lower_elements: 18,
            radial_elements: 12,
            radial_bias: 1.7,
        }
    }

    fn angular_elements(self) -> usize {
        self.right_upper_elements
            + self.top_elements
            + self.left_upper_elements
            + self.left_lower_elements
            + self.bottom_elements
            + self.right_lower_elements
    }
}

#[derive(Clone, Copy)]
struct NotchProbe {
    element_before: usize,
    element_after: usize,
}

#[derive(Clone, Copy)]
struct PlateBenchmarkProbes {
    left: NotchProbe,
    right: NotchProbe,
}

#[derive(Debug, Clone, Copy)]
struct NotchBenchmarkResult {
    sigma_y: f64,
    alpha: f64,
    first_sample_sigma_y: f64,
    second_sample_sigma_y: f64,
}

struct PlateBenchmarkModel {
    model: Model2D,
    probes: PlateBenchmarkProbes,
    sigma_nominal: f64,
}

struct T3BenchmarkModel {
    model: Model2D,
    probes: PlateBenchmarkProbes,
    sigma_nominal: f64,
}

struct T6BenchmarkModel {
    model: Model2D,
    probes: PlateBenchmarkProbes,
    sigma_nominal: f64,
}

struct T6ElementResponseData<'a> {
    element_id: usize,
    node_ids: [usize; 6],
    nodes: [&'a Node2D; 6],
    element_displacements: [f64; 12],
    constitutive_matrix: [[f64; 3]; 3],
}

struct PlateBenchmarkOutput {
    element_count: usize,
    node_count: usize,
    sigma_nominal: f64,
    left: NotchBenchmarkResult,
    right: NotchBenchmarkResult,
    max_von_mises: f64,
    max_sigma_x: f64,
    max_sigma_y: f64,
    iterations: usize,
    relative_residual: f64,
}

struct StressComponentFields {
    sigma_x: Vec<(usize, f64)>,
    sigma_y: Vec<(usize, f64)>,
    max_sigma_x: f64,
    max_sigma_y: f64,
    max_von_mises: f64,
}

struct VonMisesFieldOutput {
    max_von_mises: f64,
}

struct T3FieldOutput {
    sigma_x: Vec<(usize, f64)>,
    sigma_y: Vec<(usize, f64)>,
    von_mises: Vec<(usize, f64)>,
    displacement: Vec<(usize, f64)>,
    strain: Vec<(usize, f64)>,
    max_sigma_x: f64,
    max_sigma_y: f64,
    max_von_mises: f64,
    max_displacement: f64,
    max_strain: f64,
}

#[derive(Clone, Copy)]
struct StrainGaugeSampleDefinition {
    notch: &'static str,
    distance: f64,
    epsilon_y_per_mille: f64,
    epsilon_x_per_mille: f64,
}

#[derive(Debug, Clone, Copy)]
struct StrainGaugeSampleResult {
    notch: &'static str,
    distance: f64,
    x: f64,
    y: f64,
    sampled_element_count: usize,
    fem_sigma_x: f64,
    fem_sigma_y: f64,
    experimental_sigma_x: f64,
    experimental_sigma_y: f64,
    sigma_x_error_percent: f64,
    sigma_y_error_percent: f64,
}

struct StrainGaugeSamplingOutput {
    samples: Vec<StrainGaugeSampleResult>,
    mean_abs_sigma_y_error_percent: f64,
    max_abs_sigma_y_error_percent: f64,
}

struct T3BenchmarkOutput {
    element_count: usize,
    node_count: usize,
    sigma_nominal: f64,
    left: NotchBenchmarkResult,
    right: NotchBenchmarkResult,
    max_sigma_x: f64,
    max_sigma_y: f64,
    max_von_mises: f64,
    max_displacement: f64,
    max_strain: f64,
    mean_abs_strain_gauge_sigma_y_error_percent: f64,
    max_abs_strain_gauge_sigma_y_error_percent: f64,
    iterations: usize,
    relative_residual: f64,
}

struct T3T6VonMisesComparisonOutput {
    t3_element_count: usize,
    t3_node_count: usize,
    t3_left: NotchBenchmarkResult,
    t3_right: NotchBenchmarkResult,
    t3_max_von_mises: f64,
    t3_iterations: usize,
    t3_relative_residual: f64,
    t6_element_count: usize,
    t6_node_count: usize,
    t6_left: NotchBenchmarkResult,
    t6_right: NotchBenchmarkResult,
    t6_max_von_mises: f64,
    t6_iterations: usize,
    t6_relative_residual: f64,
    sigma_nominal: f64,
    shared_von_mises_scale: f64,
}

struct T6JsonConvergenceOutput {
    mesh: &'static str,
    element_count: usize,
    node_count: usize,
    left_sigma_y: f64,
    left_alpha: f64,
    right_sigma_y: f64,
    right_alpha: f64,
    max_von_mises: f64,
    iterations: usize,
    relative_residual: f64,
    elapsed_seconds: f64,
}

fn main() {
    run_mesh_smoke_check();

    if std::env::args().any(|argument| argument == "--mesh-only") {
        let benchmark = plate_q8_benchmark_model(PlateMeshSpec::reference_q8());

        validate_q8_corner_jacobians(&benchmark.model);
        write_model_2d_mesh_svg(&benchmark.model, MESH_SVG_PATH).expect("plate mesh SVG should be written");

        println!(
            "plate Q8 mesh SVG written to {MESH_SVG_PATH}: elements = {}, nodes = {}",
            benchmark.model.elements().len(),
            benchmark.model.nodes().len()
        );

        return;
    }

    if std::env::args().any(|argument| argument == "--t3-mesh-only") {
        let benchmark = plate_t3_benchmark_model();

        write_model_2d_mesh_svg(&benchmark.model, T3_MESH_SVG_PATH).expect("T3 plate mesh SVG should be written");

        println!(
            "plate T3 mesh SVG written to {T3_MESH_SVG_PATH}: elements = {}, nodes = {}",
            benchmark.model.elements().len(),
            benchmark.model.nodes().len()
        );

        return;
    }

    if std::env::args().any(|argument| argument == "--t3-json") {
        let benchmark = plate_t3_benchmark_model();

        write_model_2d_json_example(&benchmark.model, T3_JSON_EXAMPLE_PATH)
            .expect("T3 plate JSON example should be written");

        println!(
            "plate T3 JSON example written to {T3_JSON_EXAMPLE_PATH}: elements = {}, nodes = {}",
            benchmark.model.elements().len(),
            benchmark.model.nodes().len()
        );

        return;
    }

    if std::env::args().any(|argument| argument == "--t6-json") {
        let benchmark = plate_t6_benchmark_model();

        write_model_2d_json_example(&benchmark.model, T6_JSON_EXAMPLE_PATH)
            .expect("T6 plate JSON example should be written");

        println!(
            "plate T6 JSON example written to {T6_JSON_EXAMPLE_PATH}: elements = {}, nodes = {}",
            benchmark.model.elements().len(),
            benchmark.model.nodes().len()
        );

        return;
    }

    if std::env::args().any(|argument| argument == "--t3") {
        let started = Instant::now();
        let output = black_box(run_t3_benchmark());
        let elapsed = started.elapsed();

        println!(
            "plate T3 stress concentration benchmark: elapsed = {:.3} s, elements = {}, nodes = {}, sigma_nom = {:.6} MPa, left A(r=25) sigma_y = {:.6} MPa, alpha = {:.6} (samples {:.6}, {:.6}; strain-gauge sigma_y {:.6}, alpha {:.6}), right B(r=50) sigma_y = {:.6} MPa, alpha = {:.6} (samples {:.6}, {:.6}; strain-gauge sigma_y {:.6}, alpha {:.6}), max sigma_x = {:.6} MPa, max sigma_y = {:.6} MPa, max von Mises = {:.6} MPa, max displacement = {:.6} mm, max strain = {:.6e}, strain-gauge sampling mean |sigma_y err| = {:.3}%, max |sigma_y err| = {:.3}%, cg iters = {}, rel residual = {:.3e}, SVGs = {}, {}, {}, {}, {}, {}, sampling CSV = {}",
            elapsed.as_secs_f64(),
            output.element_count,
            output.node_count,
            output.sigma_nominal,
            output.left.sigma_y,
            output.left.alpha,
            output.left.first_sample_sigma_y,
            output.left.second_sample_sigma_y,
            EXPERIMENTAL_LEFT_SIGMA_Y,
            EXPERIMENTAL_LEFT_ALPHA,
            output.right.sigma_y,
            output.right.alpha,
            output.right.first_sample_sigma_y,
            output.right.second_sample_sigma_y,
            EXPERIMENTAL_RIGHT_SIGMA_Y,
            EXPERIMENTAL_RIGHT_ALPHA,
            output.max_sigma_x,
            output.max_sigma_y,
            output.max_von_mises,
            output.max_displacement,
            output.max_strain,
            output.mean_abs_strain_gauge_sigma_y_error_percent,
            output.max_abs_strain_gauge_sigma_y_error_percent,
            output.iterations,
            output.relative_residual,
            T3_MESH_SVG_PATH,
            T3_SIGMA_X_SVG_PATH,
            T3_SIGMA_Y_SVG_PATH,
            T3_VON_MISES_SVG_PATH,
            T3_DISPLACEMENT_SVG_PATH,
            T3_STRAIN_SVG_PATH,
            T3_STRAIN_GAUGE_SAMPLING_CSV_PATH
        );

        return;
    }

    if std::env::args().any(|argument| argument == "--t3-t6") {
        let started = Instant::now();
        let output = black_box(run_t3_t6_von_mises_comparison());
        let elapsed = started.elapsed();

        println!(
            "plate T3/T6 von Mises comparison: elapsed = {:.3} s, sigma_nom = {:.6} MPa, shared VM scale = {:.6} MPa, SVGs = {}, {}",
            elapsed.as_secs_f64(),
            output.sigma_nominal,
            output.shared_von_mises_scale,
            T3_VON_MISES_SVG_PATH,
            T6_VON_MISES_SVG_PATH
        );
        println!(
            "| element | elements | nodes | A/R25 sigma_y [MPa] | A/R25 alpha | B/R50 sigma_y [MPa] | B/R50 alpha | max von Mises [MPa] | cg iters | rel residual |"
        );
        println!("| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |");
        println!(
            "| T3 | {} | {} | {:.6} | {:.6} | {:.6} | {:.6} | {:.6} | {} | {:.3e} |",
            output.t3_element_count,
            output.t3_node_count,
            output.t3_left.sigma_y,
            output.t3_left.alpha,
            output.t3_right.sigma_y,
            output.t3_right.alpha,
            output.t3_max_von_mises,
            output.t3_iterations,
            output.t3_relative_residual
        );
        println!(
            "| T6 | {} | {} | {:.6} | {:.6} | {:.6} | {:.6} | {:.6} | {} | {:.3e} |",
            output.t6_element_count,
            output.t6_node_count,
            output.t6_left.sigma_y,
            output.t6_left.alpha,
            output.t6_right.sigma_y,
            output.t6_right.alpha,
            output.t6_max_von_mises,
            output.t6_iterations,
            output.t6_relative_residual
        );

        return;
    }

    if std::env::args().any(|argument| argument == "--t6-json-convergence") {
        let outputs = black_box(run_t6_json_mesh_convergence());

        println!(
            "| mesh | elements | nodes | A/R25 sigma_y [MPa] | A/R25 error | A/R25 alpha | B/R50 sigma_y [MPa] | B/R50 error | B/R50 alpha | max von Mises [MPa] | CG iterations | relative residual | elapsed [s] |"
        );
        println!("| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |");

        for output in outputs {
            println!(
                "| {} | {} | {} | {:.6} | {:+.2}% | {:.6} | {:.6} | {:+.2}% | {:.6} | {:.6} | {} | {:.3e} | {:.3} |",
                output.mesh,
                output.element_count,
                output.node_count,
                output.left_sigma_y,
                relative_error_percent(output.left_sigma_y, EXPERIMENTAL_LEFT_SIGMA_Y),
                output.left_alpha,
                output.right_sigma_y,
                relative_error_percent(output.right_sigma_y, EXPERIMENTAL_RIGHT_SIGMA_Y),
                output.right_alpha,
                output.max_von_mises,
                output.iterations,
                output.relative_residual,
                output.elapsed_seconds,
            );
        }

        return;
    }

    let started = Instant::now();
    let output = black_box(run_reference_benchmark());
    let elapsed = started.elapsed();

    println!(
        "plate Q8 stress concentration benchmark: elapsed = {:.3} s, elements = {}, nodes = {}, sigma_nom = {:.6} MPa, left A(r=25) sigma_y = {:.6} MPa, alpha = {:.6} (samples {:.6}, {:.6}; strain-gauge sigma_y {:.6}, alpha {:.6}), right B(r=50) sigma_y = {:.6} MPa, alpha = {:.6} (samples {:.6}, {:.6}; strain-gauge sigma_y {:.6}, alpha {:.6}), max sigma_x map = {:.6} MPa, max sigma_y map = {:.6} MPa, max von Mises = {:.6} MPa, cg iters = {}, rel residual = {:.3e}, stress SVGs = {}, {}",
        elapsed.as_secs_f64(),
        output.element_count,
        output.node_count,
        output.sigma_nominal,
        output.left.sigma_y,
        output.left.alpha,
        output.left.first_sample_sigma_y,
        output.left.second_sample_sigma_y,
        EXPERIMENTAL_LEFT_SIGMA_Y,
        EXPERIMENTAL_LEFT_ALPHA,
        output.right.sigma_y,
        output.right.alpha,
        output.right.first_sample_sigma_y,
        output.right.second_sample_sigma_y,
        EXPERIMENTAL_RIGHT_SIGMA_Y,
        EXPERIMENTAL_RIGHT_ALPHA,
        output.max_sigma_x,
        output.max_sigma_y,
        output.max_von_mises,
        output.iterations,
        output.relative_residual,
        SIGMA_X_SVG_PATH,
        SIGMA_Y_SVG_PATH
    );
}

fn run_t6_json_mesh_convergence() -> Vec<T6JsonConvergenceOutput> {
    [
        ("reference", T6_JSON_EXAMPLE_PATH),
        ("improved", T6_IMPROVED_JSON_EXAMPLE_PATH),
        ("dense", T6_DENSE_JSON_EXAMPLE_PATH),
    ]
    .into_iter()
    .map(|(mesh, path)| run_t6_json_mesh(mesh, path))
    .collect()
}

fn run_t6_json_mesh(mesh: &'static str, path: &str) -> T6JsonConvergenceOutput {
    let input: Model2DInput = serde_json::from_str(&fs::read_to_string(path).expect("T6 JSON should be readable"))
        .expect("T6 JSON should parse");
    let model = input.into_model().expect("T6 JSON model should build");

    assert!(model.elements().iter().all(|element| matches!(element, Element2D::TriangleT6(_))));

    let started = Instant::now();
    let result = solve_with_settings(&model).expect("T6 JSON model should solve");
    let responses = recover_nodal_plane_stress_responses(&model, result.displacements())
        .expect("T6 nodal responses should be recovered");
    let response_by_node = responses.iter().map(|response| (response.node_id(), response)).collect::<HashMap<_, _>>();
    let left_node = nearest_node(&model, -NOTCH_DEPTH, 0.0);
    let right_node = nearest_node(&model, NOTCH_DEPTH, 0.0);
    let left_sigma_y = response_by_node[&left_node.id()].stress()[1];
    let right_sigma_y = response_by_node[&right_node.id()].stress()[1];
    let max_von_mises = responses.iter().map(|response| response.von_mises_stress()).fold(0.0_f64, f64::max);
    let report = result.solver_report().expect("JSON benchmark should use the sparse solver");
    let sigma_nominal = REFERENCE_MASS * GRAVITY / ((WIDTH - 2.0 * NOTCH_DEPTH) * THICKNESS);

    assert!((left_node.x() + NOTCH_DEPTH).abs() < 1e-9 && left_node.y().abs() < 1e-9);
    assert!((right_node.x() - NOTCH_DEPTH).abs() < 1e-9 && right_node.y().abs() < 1e-9);
    assert!(left_sigma_y.is_finite() && right_sigma_y.is_finite() && max_von_mises.is_finite());

    T6JsonConvergenceOutput {
        mesh,
        element_count: model.elements().len(),
        node_count: model.nodes().len(),
        left_sigma_y,
        left_alpha: left_sigma_y / sigma_nominal,
        right_sigma_y,
        right_alpha: right_sigma_y / sigma_nominal,
        max_von_mises,
        iterations: report.iterations,
        relative_residual: report.relative_residual_norm,
        elapsed_seconds: started.elapsed().as_secs_f64(),
    }
}

fn nearest_node(model: &Model2D, x: f64, y: f64) -> &Node2D {
    model
        .nodes()
        .iter()
        .min_by(|first, second| {
            let first_distance = (first.x() - x).powi(2) + (first.y() - y).powi(2);
            let second_distance = (second.x() - x).powi(2) + (second.y() - y).powi(2);
            first_distance.total_cmp(&second_distance)
        })
        .expect("benchmark model should contain nodes")
}

fn run_mesh_smoke_check() {
    let spec = PlateMeshSpec::smoke_q8();
    let benchmark = plate_q8_benchmark_model(spec);

    assert_eq!(benchmark.model.elements().len(), spec.angular_elements() * spec.radial_elements);
    assert_eq!(benchmark.model.constraints().len(), 3);
    assert_eq!(benchmark.model.element_loads().len(), spec.top_elements + spec.bottom_elements);
    assert_relative_error_below(benchmark.sigma_nominal, 13.475274725274724, 1e-12, "nominal stress");
    validate_q8_corner_jacobians(&benchmark.model);
}

fn write_model_2d_json_example(model: &Model2D, path: &str) -> Result<(), Box<dyn std::error::Error>> {
    let input = model_2d_input_from_model(model);
    let contents = serde_json::to_string_pretty(&input)?;

    fs::write(path, format!("{contents}\n"))?;

    Ok(())
}

fn model_2d_input_from_model(model: &Model2D) -> Model2DInput {
    Model2DInput {
        analysis_settings: AnalysisSettings2DInput {
            solver: Some(SolverKind2DInput::Sparse),
            cg_tolerance: Some(1e-8),
            cg_max_iterations: Some(100_000),
            cg_stagnation_window: Some(0),
            cg_stagnation_tolerance: Some(1e-12),
        },
        materials: model
            .materials()
            .materials()
            .iter()
            .map(|(id, material)| Material2DInput {
                id: *id,
                young_modulus: material.young_modulus(),
                poisson_ratio: material.poisson_ratio(),
                density: material.density(),
            })
            .collect(),
        sections: model
            .sections()
            .sections()
            .iter()
            .map(|(id, section)| Section2DInput { id: *id, kind: section_input_kind(section) })
            .collect(),
        nodes: model.nodes().iter().map(|node| Node2DInput { id: node.id(), x: node.x(), y: node.y() }).collect(),
        elements: model.elements().iter().map(element_input).collect(),
        constraints: model
            .constraints()
            .iter()
            .map(|constraint| DisplacementConstraint2DInput {
                node: constraint.node_id(),
                dof: dof_input(constraint.dof()),
                value: constraint.displacement(),
            })
            .collect(),
        loads: model.element_loads().iter().map(element_load_input).collect(),
        nodal_loads: model
            .loads()
            .iter()
            .map(|load| NodalLoad2DInput { node: load.node_id(), dof: dof_input(load.dof()), value: load.value() })
            .collect(),
    }
}

fn section_input_kind(section: &Section2D) -> Section2DInputKind {
    match section {
        Section2D::Truss(section) => {
            Section2DInputKind::Truss(TrussSection2DInput { area: section.cross_section_area() })
        }
        Section2D::Beam(section) => Section2DInputKind::Beam(BeamSection2DInput {
            area: section.cross_section_area(),
            second_moment_of_area: section.second_moment_of_area(),
            height: section.section_height(),
        }),
        Section2D::PlaneStress(section) => {
            Section2DInputKind::PlaneStress(PlaneStressSection2DInput { thickness: section.thickness() })
        }
    }
}

fn element_input(element: &Element2D) -> ElementType2DInput {
    match element {
        Element2D::Truss(_) => ElementType2DInput::Truss {
            id: element.id(),
            nodes: element.node_ids().try_into().expect("truss should have two nodes"),
            material: element.material_id(),
            section: element.section_id(),
        },
        Element2D::Beam(_) => ElementType2DInput::Beam {
            id: element.id(),
            nodes: element.node_ids().try_into().expect("beam should have two nodes"),
            material: element.material_id(),
            section: element.section_id(),
        },
        Element2D::TriangleT3(_) => ElementType2DInput::TriangleT3 {
            id: element.id(),
            nodes: element.node_ids().try_into().expect("T3 should have three nodes"),
            material: element.material_id(),
            section: element.section_id(),
        },
        Element2D::TriangleT6(_) => ElementType2DInput::T6 {
            id: element.id(),
            nodes: element.node_ids().try_into().expect("T6 should have six nodes"),
            material: element.material_id(),
            section: element.section_id(),
        },
        Element2D::QuadQ4(_) => ElementType2DInput::Q4 {
            id: element.id(),
            nodes: element.node_ids().try_into().expect("Q4 should have four nodes"),
            material: element.material_id(),
            section: element.section_id(),
        },
        Element2D::QuadQ8(_) => ElementType2DInput::Q8 {
            id: element.id(),
            nodes: element.node_ids().try_into().expect("Q8 should have eight nodes"),
            material: element.material_id(),
            section: element.section_id(),
        },
    }
}

fn element_load_input(load: &ElementLoad2D) -> ElementLoad2DInput {
    let kind = match load {
        ElementLoad2D::BeamUniformLine(load) => ElementLoad2DInputKind::BeamUniform {
            element: load.element_id(),
            coordinate_system: load_coordinate_system_input(load.coordinate_system()),
            qx: load.x_component(),
            qy: load.y_component(),
        },
        ElementLoad2D::EdgeTraction(load) => ElementLoad2DInputKind::EdgeTraction {
            element: load.element_id(),
            edge: load.edge_node_ids(),
            coordinate_system: load_coordinate_system_input(load.coordinate_system()),
            tx: load.x_component(),
            ty: load.y_component(),
        },
        ElementLoad2D::BodyForce(load) => ElementLoad2DInputKind::BodyForce {
            element: load.element_id(),
            bx: load.x_component(),
            by: load.y_component(),
        },
        ElementLoad2D::SelfWeight(load) => ElementLoad2DInputKind::SelfWeight {
            element: load.element_id(),
            ax: load.x_acceleration(),
            ay: load.y_acceleration(),
        },
    };

    ElementLoad2DInput { kind }
}

fn dof_input(dof: Dof2D) -> Dof2DInput {
    match dof {
        Dof2D::Ux => Dof2DInput::Ux,
        Dof2D::Uy => Dof2DInput::Uy,
        Dof2D::Rz => Dof2DInput::Rz,
    }
}

fn load_coordinate_system_input(coordinate_system: LoadCoordinateSystem2D) -> LoadCoordinateSystem2DInput {
    match coordinate_system {
        LoadCoordinateSystem2D::Global => LoadCoordinateSystem2DInput::Global,
        LoadCoordinateSystem2D::Local => LoadCoordinateSystem2DInput::Local,
    }
}

fn run_reference_benchmark() -> PlateBenchmarkOutput {
    let benchmark = plate_q8_benchmark_model(PlateMeshSpec::reference_q8());
    validate_q8_corner_jacobians(&benchmark.model);
    let options =
        CgOptions { max_iterations: 100_000, tolerance: 1e-8, stagnation_window: 0, stagnation_tolerance: 1e-12 };
    let element_count = benchmark.model.elements().len();
    let node_count = benchmark.model.nodes().len();

    let result = solve_sparse_with_options(&benchmark.model, options).expect("plate Q8 benchmark should solve");
    let left = notch_result(&benchmark.model, result.displacements(), benchmark.probes.left, benchmark.sigma_nominal);
    let right = notch_result(&benchmark.model, result.displacements(), benchmark.probes.right, benchmark.sigma_nominal);
    let stress_fields = q8_stress_component_fields(&benchmark.model, result.displacements());
    write_stress_field_svg(
        &benchmark.model,
        &stress_fields.sigma_x,
        SIGMA_X_SVG_PATH,
        "positive sigma_x [MPa]",
        stress_fields.max_sigma_x,
    );
    write_stress_field_svg(
        &benchmark.model,
        &stress_fields.sigma_y,
        SIGMA_Y_SVG_PATH,
        "positive sigma_y [MPa]",
        stress_fields.max_sigma_y,
    );
    let report = result.solver_report().expect("sparse solver should report diagnostics");
    let iterations = report.iterations;
    let relative_residual = report.relative_residual_norm;

    assert!(left.sigma_y.is_finite());
    assert!(right.sigma_y.is_finite());
    assert!(stress_fields.max_von_mises.is_finite());
    assert!(
        left.alpha > right.alpha,
        "sharper left notch should have larger concentration: left={left:?}, right={right:?}"
    );
    assert_relative_error_below(left.alpha, EXPERIMENTAL_LEFT_ALPHA, 0.12, "left notch alpha");
    assert_relative_error_below(right.alpha, EXPERIMENTAL_RIGHT_ALPHA, 0.12, "right notch alpha");
    assert_relative_error_below(left.sigma_y, EXPERIMENTAL_LEFT_SIGMA_Y, 0.12, "left notch sigma_y");
    assert_relative_error_below(right.sigma_y, EXPERIMENTAL_RIGHT_SIGMA_Y, 0.12, "right notch sigma_y");

    PlateBenchmarkOutput {
        element_count,
        node_count,
        sigma_nominal: benchmark.sigma_nominal,
        left,
        right,
        max_von_mises: stress_fields.max_von_mises,
        max_sigma_x: stress_fields.max_sigma_x,
        max_sigma_y: stress_fields.max_sigma_y,
        iterations,
        relative_residual,
    }
}

fn run_t3_benchmark() -> T3BenchmarkOutput {
    let benchmark = plate_t3_benchmark_model();
    let options =
        CgOptions { max_iterations: 50_000, tolerance: 1e-6, stagnation_window: 0, stagnation_tolerance: 1e-12 };
    let element_count = benchmark.model.elements().len();
    let node_count = benchmark.model.nodes().len();

    write_model_2d_mesh_svg(&benchmark.model, T3_MESH_SVG_PATH).expect("T3 plate mesh SVG should be written");

    let result = solve_sparse_with_options(&benchmark.model, options).expect("plate T3 benchmark should solve");
    let left =
        t3_notch_result(&benchmark.model, result.displacements(), benchmark.probes.left, benchmark.sigma_nominal);
    let right =
        t3_notch_result(&benchmark.model, result.displacements(), benchmark.probes.right, benchmark.sigma_nominal);
    let fields = t3_field_output(&benchmark.model, result.displacements());
    let strain_gauge_sampling = t3_strain_gauge_sampling(&benchmark.model, result.displacements());

    write_stress_field_svg(
        &benchmark.model,
        &fields.sigma_x,
        T3_SIGMA_X_SVG_PATH,
        "T3 positive sigma_x [MPa]",
        fields.max_sigma_x,
    );
    write_stress_field_svg(
        &benchmark.model,
        &fields.sigma_y,
        T3_SIGMA_Y_SVG_PATH,
        "T3 positive sigma_y [MPa]",
        fields.max_sigma_y,
    );
    write_stress_field_svg(
        &benchmark.model,
        &fields.von_mises,
        T3_VON_MISES_SVG_PATH,
        "T3 von Mises stress [MPa]",
        fields.max_von_mises,
    );
    write_stress_field_svg(
        &benchmark.model,
        &fields.displacement,
        T3_DISPLACEMENT_SVG_PATH,
        "T3 displacement magnitude [mm]",
        fields.max_displacement,
    );
    write_stress_field_svg(
        &benchmark.model,
        &fields.strain,
        T3_STRAIN_SVG_PATH,
        "T3 strain magnitude [-]",
        fields.max_strain,
    );
    write_strain_gauge_sampling_csv(&strain_gauge_sampling.samples, T3_STRAIN_GAUGE_SAMPLING_CSV_PATH)
        .expect("T3 strain-gauge sampling CSV should be written");

    let report = result.solver_report().expect("sparse solver should report diagnostics");

    T3BenchmarkOutput {
        element_count,
        node_count,
        sigma_nominal: benchmark.sigma_nominal,
        left,
        right,
        max_sigma_x: fields.max_sigma_x,
        max_sigma_y: fields.max_sigma_y,
        max_von_mises: fields.max_von_mises,
        max_displacement: fields.max_displacement,
        max_strain: fields.max_strain,
        mean_abs_strain_gauge_sigma_y_error_percent: strain_gauge_sampling.mean_abs_sigma_y_error_percent,
        max_abs_strain_gauge_sigma_y_error_percent: strain_gauge_sampling.max_abs_sigma_y_error_percent,
        iterations: report.iterations,
        relative_residual: report.relative_residual_norm,
    }
}

fn run_t3_t6_von_mises_comparison() -> T3T6VonMisesComparisonOutput {
    let options =
        CgOptions { max_iterations: 100_000, tolerance: 1e-6, stagnation_window: 0, stagnation_tolerance: 1e-12 };
    let t3_benchmark = plate_t3_benchmark_model();
    let t6_benchmark = plate_t6_benchmark_model();
    let t3_element_count = t3_benchmark.model.elements().len();
    let t3_node_count = t3_benchmark.model.nodes().len();
    let t6_element_count = t6_benchmark.model.elements().len();
    let t6_node_count = t6_benchmark.model.nodes().len();

    let t3_result = solve_sparse_with_options(&t3_benchmark.model, options).expect("plate T3 benchmark should solve");
    let t6_result = solve_sparse_with_options(&t6_benchmark.model, options).expect("plate T6 benchmark should solve");
    let t3_fields = t3_field_output(&t3_benchmark.model, t3_result.displacements());
    let t6_response_data = t6_element_response_data(&t6_benchmark.model, t6_result.displacements());
    let t6_fields = t6_von_mises_field_output(&t6_response_data);
    let t6_von_mises_patches = t6_von_mises_gradient_patches(&t6_response_data, T6_VON_MISES_VISUAL_SUBDIVISIONS);
    let shared_von_mises_scale = t3_fields.max_von_mises.max(t6_fields.max_von_mises);
    let t3_left = t3_notch_result(
        &t3_benchmark.model,
        t3_result.displacements(),
        t3_benchmark.probes.left,
        t3_benchmark.sigma_nominal,
    );
    let t3_right = t3_notch_result(
        &t3_benchmark.model,
        t3_result.displacements(),
        t3_benchmark.probes.right,
        t3_benchmark.sigma_nominal,
    );
    let t6_left = t6_notch_result(
        &t6_benchmark.model,
        t6_result.displacements(),
        t6_benchmark.probes.left,
        t6_benchmark.sigma_nominal,
    );
    let t6_right = t6_notch_result(
        &t6_benchmark.model,
        t6_result.displacements(),
        t6_benchmark.probes.right,
        t6_benchmark.sigma_nominal,
    );

    write_stress_field_svg(
        &t3_benchmark.model,
        &t3_fields.von_mises,
        T3_VON_MISES_SVG_PATH,
        "T3 von Mises stress [MPa]",
        shared_von_mises_scale,
    );
    write_stress_patch_svg(
        &t6_benchmark.model,
        &t6_von_mises_patches,
        T6_VON_MISES_SVG_PATH,
        "T6 von Mises stress [MPa]",
        shared_von_mises_scale,
    );

    let t3_report = t3_result.solver_report().expect("sparse solver should report diagnostics");
    let t6_report = t6_result.solver_report().expect("sparse solver should report diagnostics");

    T3T6VonMisesComparisonOutput {
        t3_element_count,
        t3_node_count,
        t3_left,
        t3_right,
        t3_max_von_mises: t3_fields.max_von_mises,
        t3_iterations: t3_report.iterations,
        t3_relative_residual: t3_report.relative_residual_norm,
        t6_element_count,
        t6_node_count,
        t6_left,
        t6_right,
        t6_max_von_mises: t6_fields.max_von_mises,
        t6_iterations: t6_report.iterations,
        t6_relative_residual: t6_report.relative_residual_norm,
        sigma_nominal: t3_benchmark.sigma_nominal,
        shared_von_mises_scale,
    }
}

fn plate_q8_benchmark_model(spec: PlateMeshSpec) -> PlateBenchmarkModel {
    assert!(spec.angular_elements() > 0);
    assert!(spec.radial_elements > 0);
    assert!(spec.radial_bias > 1.0);

    let load = REFERENCE_MASS * GRAVITY;
    let sigma_nominal = load / ((WIDTH - 2.0 * NOTCH_DEPTH) * THICKNESS);
    let edge_traction = load / (WIDTH * THICKNESS);
    let angles = angular_node_angles(spec);
    let angular_node_count = angles.len();
    let angular_element_count = spec.angular_elements();
    let radial_node_count = 2 * spec.radial_elements + 1;
    let mut model = Model2D::new();

    model.set_material(Material2D::new(YOUNG_MODULUS, POISSON_RATIO, 1.0).expect("valid material"));
    model
        .add_section(1, Section2D::PlaneStress(PlaneStressSection2D::new(THICKNESS).expect("valid section")))
        .expect("section should be added");

    for radial_index in 0..radial_node_count {
        let t = radial_fraction(radial_index, spec.radial_elements, spec.radial_bias);

        for (angular_index, angle) in angles.iter().copied().enumerate() {
            let (x, y) = mapped_point(angle, t);
            let node_id = node_id(radial_index, angular_index, angular_node_count);

            model.add_node(Node2D::new(node_id, x, y).expect("valid node")).expect("node should be added");
        }
    }

    let bottom_left = closest_node_id(model.nodes(), -WIDTH / 2.0, -HEIGHT / 2.0);
    let bottom_right = closest_node_id(model.nodes(), WIDTH / 2.0, -HEIGHT / 2.0);

    model
        .add_constraint(DisplacementConstraint2D::new(bottom_left, Dof2D::Ux, 0.0).expect("valid constraint"))
        .expect("constraint should be added");
    model
        .add_constraint(DisplacementConstraint2D::new(bottom_left, Dof2D::Uy, 0.0).expect("valid constraint"))
        .expect("constraint should be added");
    model
        .add_constraint(DisplacementConstraint2D::new(bottom_right, Dof2D::Uy, 0.0).expect("valid constraint"))
        .expect("constraint should be added");

    for radial_element in 0..spec.radial_elements {
        for angular_element in 0..angular_element_count {
            let element_id = element_id(radial_element, angular_element, angular_element_count);
            let nodes = q8_element_node_ids(radial_element, angular_element, angular_node_count);
            let quad = QuadQ8::new(element_id, nodes, DEFAULT_MATERIAL_ID, 1).expect("valid Q8");

            model.add_element(Element2D::QuadQ8(quad)).expect("Q8 element should be added");
        }
    }

    let outer_radial_element = spec.radial_elements - 1;

    for angular_element in 0..angular_element_count {
        let mid_angle = angles[(2 * angular_element + 1) % angular_node_count];
        let (_, outer_y) = outer_boundary(mid_angle);
        let traction_y = if (outer_y - HEIGHT / 2.0).abs() < 1e-8 {
            Some(edge_traction)
        } else if (outer_y + HEIGHT / 2.0).abs() < 1e-8 {
            Some(-edge_traction)
        } else {
            None
        };

        if let Some(ty) = traction_y {
            let element_id = element_id(outer_radial_element, angular_element, angular_element_count);
            let nodes = q8_element_node_ids(outer_radial_element, angular_element, angular_node_count);
            let load = ElementLoad2D::EdgeTraction(
                EdgeTraction2D::new(element_id, [nodes[1], nodes[2]], LoadCoordinateSystem2D::Global, 0.0, ty)
                    .expect("valid edge traction"),
            );

            model.add_element_load(load).expect("edge traction should be added");
        }
    }

    let right_tip_angular_element = 0;
    let left_tip_angular_node = angles
        .iter()
        .position(|angle| (*angle - std::f64::consts::PI).abs() < 1e-12)
        .expect("angular grid should contain the left notch tip");
    assert_eq!(left_tip_angular_node % 2, 0);
    let left_tip_angular_element = left_tip_angular_node / 2;
    let probes = PlateBenchmarkProbes {
        left: notch_probe(left_tip_angular_element, angular_element_count),
        right: notch_probe(right_tip_angular_element, angular_element_count),
    };

    PlateBenchmarkModel { model, probes, sigma_nominal }
}

fn plate_t3_benchmark_model() -> T3BenchmarkModel {
    let load = REFERENCE_MASS * GRAVITY;
    let sigma_nominal = load / ((WIDTH - 2.0 * NOTCH_DEPTH) * THICKNESS);
    let edge_traction = load / (WIDTH * THICKNESS);
    let angles = t3_visual_angular_samples();
    let angular_node_count = angles.len();
    let radial_node_count = T3_VISUAL_RADIAL_FRACTIONS.len();
    let mut model = Model2D::new();

    model.set_material(Material2D::new(YOUNG_MODULUS, POISSON_RATIO, 1.0).expect("valid material"));
    model
        .add_section(1, Section2D::PlaneStress(PlaneStressSection2D::new(THICKNESS).expect("valid section")))
        .expect("section should be added");

    for (radial_index, t) in T3_VISUAL_RADIAL_FRACTIONS.iter().copied().enumerate() {
        for (angular_index, angle) in angles.iter().copied().enumerate() {
            let (x, y) = mapped_point(angle, t);
            let node_id = node_id(radial_index, angular_index, angular_node_count);

            model.add_node(Node2D::new(node_id, x, y).expect("valid node")).expect("node should be added");
        }
    }

    let bottom_left = closest_node_id(model.nodes(), -WIDTH / 2.0, -HEIGHT / 2.0);
    let bottom_right = closest_node_id(model.nodes(), WIDTH / 2.0, -HEIGHT / 2.0);

    model
        .add_constraint(DisplacementConstraint2D::new(bottom_left, Dof2D::Ux, 0.0).expect("valid constraint"))
        .expect("constraint should be added");
    model
        .add_constraint(DisplacementConstraint2D::new(bottom_left, Dof2D::Uy, 0.0).expect("valid constraint"))
        .expect("constraint should be added");
    model
        .add_constraint(DisplacementConstraint2D::new(bottom_right, Dof2D::Uy, 0.0).expect("valid constraint"))
        .expect("constraint should be added");

    for radial_element in 0..radial_node_count - 1 {
        for angular_element in 0..angular_node_count {
            let next_angular = (angular_element + 1) % angular_node_count;
            let inner_first = node_id(radial_element, angular_element, angular_node_count);
            let outer_first = node_id(radial_element + 1, angular_element, angular_node_count);
            let outer_second = node_id(radial_element + 1, next_angular, angular_node_count);
            let inner_second = node_id(radial_element, next_angular, angular_node_count);
            let (first_element_id, second_element_id) =
                t3_element_ids(radial_element, angular_element, angular_node_count);
            let (first_nodes, second_nodes) = if (radial_element + angular_element).is_multiple_of(2) {
                ([inner_first, outer_first, outer_second], [inner_first, outer_second, inner_second])
            } else {
                ([inner_first, outer_first, inner_second], [inner_second, outer_first, outer_second])
            };

            let first_triangle = TriangleT3::new(first_element_id, first_nodes, DEFAULT_MATERIAL_ID, 1)
                .expect("valid T3 element should be created");
            let second_triangle = TriangleT3::new(second_element_id, second_nodes, DEFAULT_MATERIAL_ID, 1)
                .expect("valid T3 element should be created");

            model.add_element(Element2D::TriangleT3(first_triangle)).expect("first T3 element should be added");
            model.add_element(Element2D::TriangleT3(second_triangle)).expect("second T3 element should be added");
        }
    }

    let outer_radial_element = radial_node_count - 2;

    for angular_element in 0..angular_node_count {
        let next_angular = (angular_element + 1) % angular_node_count;
        let first_node_id = node_id(radial_node_count - 1, angular_element, angular_node_count);
        let second_node_id = node_id(radial_node_count - 1, next_angular, angular_node_count);
        let first_node = find_node(&model, first_node_id);
        let second_node = find_node(&model, second_node_id);
        let traction_y =
            if (first_node.y() - HEIGHT / 2.0).abs() < 1e-8 && (second_node.y() - HEIGHT / 2.0).abs() < 1e-8 {
                Some(edge_traction)
            } else if (first_node.y() + HEIGHT / 2.0).abs() < 1e-8 && (second_node.y() + HEIGHT / 2.0).abs() < 1e-8 {
                Some(-edge_traction)
            } else {
                None
            };

        if let Some(ty) = traction_y {
            let (first_element_id, second_element_id) =
                t3_element_ids(outer_radial_element, angular_element, angular_node_count);
            let element_id = if (outer_radial_element + angular_element).is_multiple_of(2) {
                first_element_id
            } else {
                second_element_id
            };
            let load = ElementLoad2D::EdgeTraction(
                EdgeTraction2D::new(
                    element_id,
                    [first_node_id, second_node_id],
                    LoadCoordinateSystem2D::Global,
                    0.0,
                    ty,
                )
                .expect("valid T3 edge traction"),
            );

            model.add_element_load(load).expect("T3 edge traction should be added");
        }
    }

    let right_tip_angular_node = angles
        .iter()
        .position(|angle| (*angle).abs() < 1e-12)
        .expect("angular grid should contain the right notch tip");
    let left_tip_angular_node = angles
        .iter()
        .position(|angle| (*angle - std::f64::consts::PI).abs() < 1e-12)
        .expect("angular grid should contain the left notch tip");
    let probes = PlateBenchmarkProbes {
        left: t3_notch_probe(left_tip_angular_node, angular_node_count),
        right: t3_notch_probe(right_tip_angular_node, angular_node_count),
    };

    T3BenchmarkModel { model, probes, sigma_nominal }
}

fn plate_t6_benchmark_model() -> T6BenchmarkModel {
    let load = REFERENCE_MASS * GRAVITY;
    let sigma_nominal = load / ((WIDTH - 2.0 * NOTCH_DEPTH) * THICKNESS);
    let edge_traction = load / (WIDTH * THICKNESS);
    let angles = t3_visual_angular_samples();
    let angular_node_count = angles.len();
    let radial_node_count = T3_VISUAL_RADIAL_FRACTIONS.len();
    let mut model = Model2D::new();

    model.set_material(Material2D::new(YOUNG_MODULUS, POISSON_RATIO, 1.0).expect("valid material"));
    model
        .add_section(1, Section2D::PlaneStress(PlaneStressSection2D::new(THICKNESS).expect("valid section")))
        .expect("section should be added");

    for (radial_index, t) in T3_VISUAL_RADIAL_FRACTIONS.iter().copied().enumerate() {
        for (angular_index, angle) in angles.iter().copied().enumerate() {
            let (x, y) = mapped_point(angle, t);
            let node_id = node_id(radial_index, angular_index, angular_node_count);

            model.add_node(Node2D::new(node_id, x, y).expect("valid node")).expect("node should be added");
        }
    }

    let bottom_left = closest_node_id(model.nodes(), -WIDTH / 2.0, -HEIGHT / 2.0);
    let bottom_right = closest_node_id(model.nodes(), WIDTH / 2.0, -HEIGHT / 2.0);

    model
        .add_constraint(DisplacementConstraint2D::new(bottom_left, Dof2D::Ux, 0.0).expect("valid constraint"))
        .expect("constraint should be added");
    model
        .add_constraint(DisplacementConstraint2D::new(bottom_left, Dof2D::Uy, 0.0).expect("valid constraint"))
        .expect("constraint should be added");
    model
        .add_constraint(DisplacementConstraint2D::new(bottom_right, Dof2D::Uy, 0.0).expect("valid constraint"))
        .expect("constraint should be added");

    let mut midside_nodes = HashMap::new();
    let mut next_midside_node_id = radial_node_count * angular_node_count + 1;

    for radial_element in 0..radial_node_count - 1 {
        for angular_element in 0..angular_node_count {
            let next_angular = (angular_element + 1) % angular_node_count;
            let inner_first = node_id(radial_element, angular_element, angular_node_count);
            let outer_first = node_id(radial_element + 1, angular_element, angular_node_count);
            let outer_second = node_id(radial_element + 1, next_angular, angular_node_count);
            let inner_second = node_id(radial_element, next_angular, angular_node_count);
            let (first_element_id, second_element_id) =
                t3_element_ids(radial_element, angular_element, angular_node_count);
            let (first_nodes, second_nodes) = if (radial_element + angular_element).is_multiple_of(2) {
                ([inner_first, outer_first, outer_second], [inner_first, outer_second, inner_second])
            } else {
                ([inner_first, outer_first, inner_second], [inner_second, outer_first, outer_second])
            };

            add_t6_triangle_from_corner_nodes(
                &mut model,
                &mut midside_nodes,
                &mut next_midside_node_id,
                &angles,
                angular_node_count,
                first_element_id,
                first_nodes,
            );
            add_t6_triangle_from_corner_nodes(
                &mut model,
                &mut midside_nodes,
                &mut next_midside_node_id,
                &angles,
                angular_node_count,
                second_element_id,
                second_nodes,
            );
        }
    }

    let outer_radial_element = radial_node_count - 2;

    for angular_element in 0..angular_node_count {
        let next_angular = (angular_element + 1) % angular_node_count;
        let first_node_id = node_id(radial_node_count - 1, angular_element, angular_node_count);
        let second_node_id = node_id(radial_node_count - 1, next_angular, angular_node_count);
        let first_node = find_node(&model, first_node_id);
        let second_node = find_node(&model, second_node_id);
        let traction_y =
            if (first_node.y() - HEIGHT / 2.0).abs() < 1e-8 && (second_node.y() - HEIGHT / 2.0).abs() < 1e-8 {
                Some(edge_traction)
            } else if (first_node.y() + HEIGHT / 2.0).abs() < 1e-8 && (second_node.y() + HEIGHT / 2.0).abs() < 1e-8 {
                Some(-edge_traction)
            } else {
                None
            };

        if let Some(ty) = traction_y {
            let (first_element_id, second_element_id) =
                t3_element_ids(outer_radial_element, angular_element, angular_node_count);
            let element_id = if (outer_radial_element + angular_element).is_multiple_of(2) {
                first_element_id
            } else {
                second_element_id
            };
            let load = ElementLoad2D::EdgeTraction(
                EdgeTraction2D::new(
                    element_id,
                    [first_node_id, second_node_id],
                    LoadCoordinateSystem2D::Global,
                    0.0,
                    ty,
                )
                .expect("valid T6 edge traction"),
            );

            model.add_element_load(load).expect("T6 edge traction should be added");
        }
    }

    let right_tip_angular_node = angles
        .iter()
        .position(|angle| (*angle).abs() < 1e-12)
        .expect("angular grid should contain the right notch tip");
    let left_tip_angular_node = angles
        .iter()
        .position(|angle| (*angle - std::f64::consts::PI).abs() < 1e-12)
        .expect("angular grid should contain the left notch tip");
    let probes = PlateBenchmarkProbes {
        left: t3_notch_probe(left_tip_angular_node, angular_node_count),
        right: t3_notch_probe(right_tip_angular_node, angular_node_count),
    };

    T6BenchmarkModel { model, probes, sigma_nominal }
}

fn angular_node_angles(spec: PlateMeshSpec) -> Vec<f64> {
    let corner_angle = (HEIGHT / WIDTH).atan();
    let segments = [
        AngularSegment {
            end: corner_angle,
            element_count: spec.right_upper_elements,
            spacing: AngularSpacing::Start(1.25),
        },
        AngularSegment {
            end: std::f64::consts::PI - corner_angle,
            element_count: spec.top_elements,
            spacing: AngularSpacing::Middle(1.2),
        },
        AngularSegment {
            end: std::f64::consts::PI,
            element_count: spec.left_upper_elements,
            spacing: AngularSpacing::End(1.25),
        },
        AngularSegment {
            end: std::f64::consts::PI + corner_angle,
            element_count: spec.left_lower_elements,
            spacing: AngularSpacing::Start(1.25),
        },
        AngularSegment {
            end: 2.0 * std::f64::consts::PI - corner_angle,
            element_count: spec.bottom_elements,
            spacing: AngularSpacing::Middle(1.2),
        },
        AngularSegment {
            end: 2.0 * std::f64::consts::PI,
            element_count: spec.right_lower_elements,
            spacing: AngularSpacing::End(1.25),
        },
    ];
    let mut angles = Vec::with_capacity(2 * spec.angular_elements());
    let mut start = 0.0;

    for segment in segments {
        let interval_count = 2 * segment.element_count;

        for interval in 0..interval_count {
            let fraction = interval as f64 / interval_count as f64;
            let spaced_fraction = segment.spacing.map_fraction(fraction);

            angles.push(start + (segment.end - start) * spaced_fraction);
        }

        start = segment.end;
    }

    angles
}

fn t3_visual_angular_samples() -> Vec<f64> {
    const INTEGRATION_STEPS: usize = 2400;

    let mut cumulative = Vec::with_capacity(INTEGRATION_STEPS + 1);
    let mut total = 0.0;

    cumulative.push((0.0, 0.0));

    for index in 1..=INTEGRATION_STEPS {
        let previous_angle = (index - 1) as f64 / INTEGRATION_STEPS as f64 * 2.0 * std::f64::consts::PI;
        let angle = index as f64 / INTEGRATION_STEPS as f64 * 2.0 * std::f64::consts::PI;

        total += 0.5 * (t3_visual_weight(previous_angle) + t3_visual_weight(angle)) * (angle - previous_angle);
        cumulative.push((angle, total));
    }

    let mut angles = Vec::with_capacity(T3_VISUAL_ANGLE_COUNT);

    for index in 0..T3_VISUAL_ANGLE_COUNT {
        let target = index as f64 / T3_VISUAL_ANGLE_COUNT as f64 * total;
        let mut low = 0;
        let mut high = cumulative.len() - 1;

        while high - low > 1 {
            let mid = (low + high) / 2;

            if cumulative[mid].1 < target {
                low = mid;
            } else {
                high = mid;
            }
        }

        let (left_angle, left_value) = cumulative[low];
        let (right_angle, right_value) = cumulative[high];
        let ratio = (target - left_value) / (right_value - left_value).max(1e-12);

        angles.push(left_angle + ratio * (right_angle - left_angle));
    }

    enforce_t3_visual_landmark_angles(&mut angles);

    angles
}

fn enforce_t3_visual_landmark_angles(angles: &mut [f64]) {
    let mut locked = vec![false; angles.len()];

    for landmark in t3_visual_landmark_angles() {
        let nearest_index = angles
            .iter()
            .enumerate()
            .filter(|(index, _)| !locked[*index])
            .min_by(|(_, first), (_, second)| {
                wrapped_angular_distance(**first, landmark).total_cmp(&wrapped_angular_distance(**second, landmark))
            })
            .map(|(index, _)| index)
            .expect("visual T3 mesh should contain angular samples");

        angles[nearest_index] = landmark;
        locked[nearest_index] = true;
    }

    angles.sort_by(f64::total_cmp);
}

fn t3_visual_landmark_angles() -> [f64; 8] {
    let corner_angle = (HEIGHT / WIDTH).atan();

    [
        0.0,
        corner_angle,
        std::f64::consts::FRAC_PI_2,
        std::f64::consts::PI - corner_angle,
        std::f64::consts::PI,
        std::f64::consts::PI + corner_angle,
        3.0 * std::f64::consts::FRAC_PI_2,
        2.0 * std::f64::consts::PI - corner_angle,
    ]
}

fn t3_visual_weight(angle: f64) -> f64 {
    1.0 + gaussian_angular_weight(angle, 0.0, 0.26, 2.7)
        + gaussian_angular_weight(angle, std::f64::consts::PI, 0.24, 4.0)
        + gaussian_angular_weight(angle, std::f64::consts::FRAC_PI_2, 0.30, 1.75)
        + gaussian_angular_weight(angle, 3.0 * std::f64::consts::FRAC_PI_2, 0.30, 1.75)
}

fn gaussian_angular_weight(angle: f64, target: f64, width: f64, amplitude: f64) -> f64 {
    let distance = wrapped_angular_distance(angle, target);

    amplitude * (-(distance * distance) / (2.0 * width * width)).exp()
}

fn wrapped_angular_distance(angle: f64, target: f64) -> f64 {
    let distance = (angle - target).abs() % (2.0 * std::f64::consts::PI);

    distance.min(2.0 * std::f64::consts::PI - distance)
}

#[derive(Clone, Copy)]
struct AngularSegment {
    end: f64,
    element_count: usize,
    spacing: AngularSpacing,
}

#[derive(Clone, Copy)]
enum AngularSpacing {
    Start(f64),
    End(f64),
    Middle(f64),
}

impl AngularSpacing {
    fn map_fraction(self, fraction: f64) -> f64 {
        match self {
            Self::Start(strength) => fraction.powf(strength),
            Self::End(strength) => 1.0 - (1.0 - fraction).powf(strength),
            Self::Middle(strength) => {
                if fraction <= 0.5 {
                    0.5 * (1.0 - (1.0 - 2.0 * fraction).powf(strength))
                } else {
                    0.5 + 0.5 * (2.0 * fraction - 1.0).powf(strength)
                }
            }
        }
    }
}

fn radial_fraction(radial_index: usize, radial_elements: usize, bias: f64) -> f64 {
    (radial_index as f64 / (2 * radial_elements) as f64).powf(bias)
}

fn mapped_point(angle: f64, radial_fraction: f64) -> (f64, f64) {
    let (inner_x, inner_y) = inner_notch_boundary(angle);
    let (outer_x, outer_y) = outer_boundary(angle);

    (inner_x + radial_fraction * (outer_x - inner_x), inner_y + radial_fraction * (outer_y - inner_y))
}

fn inner_notch_boundary(angle: f64) -> (f64, f64) {
    let left_semi_height = (NOTCH_DEPTH * LEFT_NOTCH_RADIUS).sqrt();
    let right_semi_height = (NOTCH_DEPTH * RIGHT_NOTCH_RADIUS).sqrt();
    let mean_semi_height = 0.5 * (left_semi_height + right_semi_height);
    let half_difference = 0.5 * (right_semi_height - left_semi_height);
    let local_semi_height = mean_semi_height + half_difference * angle.cos();

    (NOTCH_DEPTH * angle.cos(), local_semi_height * angle.sin())
}

fn outer_boundary(angle: f64) -> (f64, f64) {
    let cosine = angle.cos();
    let sine = angle.sin();
    let x_scale = if cosine.abs() < 1e-12 { f64::INFINITY } else { WIDTH / 2.0 / cosine.abs() };
    let y_scale = if sine.abs() < 1e-12 { f64::INFINITY } else { HEIGHT / 2.0 / sine.abs() };
    let scale = x_scale.min(y_scale);

    (scale * cosine, scale * sine)
}

fn node_id(radial_index: usize, angular_index: usize, angular_node_count: usize) -> usize {
    radial_index * angular_node_count + angular_index + 1
}

#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
struct MeshEdgeKey {
    first: usize,
    second: usize,
}

impl MeshEdgeKey {
    fn new(first: usize, second: usize) -> Self {
        if first <= second { Self { first, second } } else { Self { first: second, second: first } }
    }
}

fn add_t6_triangle_from_corner_nodes(
    model: &mut Model2D, midside_nodes: &mut HashMap<MeshEdgeKey, usize>, next_midside_node_id: &mut usize,
    angles: &[f64], angular_node_count: usize, element_id: usize, corner_nodes: [usize; 3],
) {
    let midside_12 = t6_midside_node_id(
        model,
        midside_nodes,
        next_midside_node_id,
        angles,
        angular_node_count,
        corner_nodes[0],
        corner_nodes[1],
    );
    let midside_23 = t6_midside_node_id(
        model,
        midside_nodes,
        next_midside_node_id,
        angles,
        angular_node_count,
        corner_nodes[1],
        corner_nodes[2],
    );
    let midside_31 = t6_midside_node_id(
        model,
        midside_nodes,
        next_midside_node_id,
        angles,
        angular_node_count,
        corner_nodes[2],
        corner_nodes[0],
    );
    let triangle = TriangleT6::new(
        element_id,
        [corner_nodes[0], corner_nodes[1], corner_nodes[2], midside_12, midside_23, midside_31],
        DEFAULT_MATERIAL_ID,
        1,
    )
    .expect("valid T6 element should be created");

    model.add_element(Element2D::TriangleT6(triangle)).expect("T6 element should be added");
}

fn t6_midside_node_id(
    model: &mut Model2D, midside_nodes: &mut HashMap<MeshEdgeKey, usize>, next_midside_node_id: &mut usize,
    angles: &[f64], angular_node_count: usize, first_node_id: usize, second_node_id: usize,
) -> usize {
    let key = MeshEdgeKey::new(first_node_id, second_node_id);

    if let Some(node_id) = midside_nodes.get(&key) {
        return *node_id;
    }

    let (x, y) = t6_midside_coordinates(model, angles, angular_node_count, first_node_id, second_node_id);
    let node_id = *next_midside_node_id;

    *next_midside_node_id += 1;
    model.add_node(Node2D::new(node_id, x, y).expect("valid T6 midside node")).expect("midside node should be added");
    midside_nodes.insert(key, node_id);

    node_id
}

fn t6_midside_coordinates(
    model: &Model2D, angles: &[f64], angular_node_count: usize, first_node_id: usize, second_node_id: usize,
) -> (f64, f64) {
    let (first_radial, first_angular) = corner_grid_indices(first_node_id, angular_node_count);
    let (second_radial, second_angular) = corner_grid_indices(second_node_id, angular_node_count);

    if first_radial == 0 && second_radial == 0 {
        let midpoint_angle = angular_midpoint(angles[first_angular], angles[second_angular]);

        return mapped_point(midpoint_angle, 0.0);
    }

    let first_node = find_node(model, first_node_id);
    let second_node = find_node(model, second_node_id);

    (0.5 * (first_node.x() + second_node.x()), 0.5 * (first_node.y() + second_node.y()))
}

fn corner_grid_indices(node_id: usize, angular_node_count: usize) -> (usize, usize) {
    let zero_based = node_id - 1;

    (zero_based / angular_node_count, zero_based % angular_node_count)
}

fn angular_midpoint(first: f64, second: f64) -> f64 {
    let full_turn = 2.0 * std::f64::consts::PI;
    let mut delta = second - first;

    if delta > std::f64::consts::PI {
        delta -= full_turn;
    } else if delta < -std::f64::consts::PI {
        delta += full_turn;
    }

    (first + 0.5 * delta).rem_euclid(full_turn)
}

fn element_id(radial_element: usize, angular_element: usize, angular_element_count: usize) -> usize {
    radial_element * angular_element_count + angular_element + 1
}

fn t3_element_ids(radial_element: usize, angular_element: usize, angular_node_count: usize) -> (usize, usize) {
    let first = 2 * (radial_element * angular_node_count + angular_element) + 1;

    (first, first + 1)
}

fn t3_inner_edge_element_id(radial_element: usize, angular_element: usize, angular_node_count: usize) -> usize {
    let (first_element_id, second_element_id) = t3_element_ids(radial_element, angular_element, angular_node_count);

    if (radial_element + angular_element).is_multiple_of(2) { second_element_id } else { first_element_id }
}

fn q8_element_node_ids(radial_element: usize, angular_element: usize, angular_node_count: usize) -> [usize; 8] {
    let first_radial = 2 * radial_element;
    let second_radial = first_radial + 2;
    let middle_radial = first_radial + 1;
    let first_angle = 2 * angular_element;
    let second_angle = (first_angle + 2) % angular_node_count;
    let middle_angle = (first_angle + 1) % angular_node_count;

    [
        node_id(first_radial, first_angle, angular_node_count),
        node_id(second_radial, first_angle, angular_node_count),
        node_id(second_radial, second_angle, angular_node_count),
        node_id(first_radial, second_angle, angular_node_count),
        node_id(middle_radial, first_angle, angular_node_count),
        node_id(second_radial, middle_angle, angular_node_count),
        node_id(middle_radial, second_angle, angular_node_count),
        node_id(first_radial, middle_angle, angular_node_count),
    ]
}

fn closest_node_id(nodes: &[Node2D], target_x: f64, target_y: f64) -> usize {
    nodes
        .iter()
        .min_by(|first, second| {
            let first_distance = squared_distance(first, target_x, target_y);
            let second_distance = squared_distance(second, target_x, target_y);

            first_distance.total_cmp(&second_distance)
        })
        .expect("model should contain nodes")
        .id()
}

fn squared_distance(node: &Node2D, target_x: f64, target_y: f64) -> f64 {
    (node.x() - target_x).powi(2) + (node.y() - target_y).powi(2)
}

fn notch_probe(tip_angular_element: usize, angular_element_count: usize) -> NotchProbe {
    let element_before = if tip_angular_element == 0 { angular_element_count } else { tip_angular_element };
    let element_after = tip_angular_element + 1;

    NotchProbe { element_before, element_after }
}

fn t3_notch_probe(tip_angular_node: usize, angular_node_count: usize) -> NotchProbe {
    let element_before = t3_inner_edge_element_id(
        0,
        if tip_angular_node == 0 { angular_node_count - 1 } else { tip_angular_node - 1 },
        angular_node_count,
    );
    let element_after = t3_inner_edge_element_id(0, tip_angular_node, angular_node_count);

    NotchProbe { element_before, element_after }
}

fn notch_result(
    model: &Model2D, displacements: &nalgebra::DVector<f64>, probe: NotchProbe, sigma_nominal: f64,
) -> NotchBenchmarkResult {
    let first_sample_sigma_y = q8_sigma_y(model, displacements, probe.element_before, -1.0, 1.0);
    let second_sample_sigma_y = q8_sigma_y(model, displacements, probe.element_after, -1.0, -1.0);
    let sigma_y = first_sample_sigma_y.max(second_sample_sigma_y);
    let alpha = sigma_y / sigma_nominal;

    NotchBenchmarkResult { sigma_y, alpha, first_sample_sigma_y, second_sample_sigma_y }
}

fn t3_notch_result(
    model: &Model2D, displacements: &nalgebra::DVector<f64>, probe: NotchProbe, sigma_nominal: f64,
) -> NotchBenchmarkResult {
    let first_sample_sigma_y = t3_sigma_y(model, displacements, probe.element_before);
    let second_sample_sigma_y = t3_sigma_y(model, displacements, probe.element_after);
    let sigma_y = first_sample_sigma_y.max(second_sample_sigma_y);
    let alpha = sigma_y / sigma_nominal;

    NotchBenchmarkResult { sigma_y, alpha, first_sample_sigma_y, second_sample_sigma_y }
}

fn t6_notch_result(
    model: &Model2D, displacements: &nalgebra::DVector<f64>, probe: NotchProbe, sigma_nominal: f64,
) -> NotchBenchmarkResult {
    let first_sample_sigma_y = t6_element_max_corner_sigma_y(model, displacements, probe.element_before);
    let second_sample_sigma_y = t6_element_max_corner_sigma_y(model, displacements, probe.element_after);
    let sigma_y = first_sample_sigma_y.max(second_sample_sigma_y);
    let alpha = sigma_y / sigma_nominal;

    NotchBenchmarkResult { sigma_y, alpha, first_sample_sigma_y, second_sample_sigma_y }
}

fn t3_sigma_y(model: &Model2D, displacements: &nalgebra::DVector<f64>, element_id: usize) -> f64 {
    let element = model.elements().iter().find(|element| element.id() == element_id).expect("element should exist");
    let Element2D::TriangleT3(triangle) = element else {
        panic!("T3 plate benchmark should contain only T3 elements");
    };
    let response = recover_triangle_response(model, triangle, displacements).expect("T3 stress should be recovered");

    response.stress()[1]
}

fn t6_element_max_corner_sigma_y(model: &Model2D, displacements: &nalgebra::DVector<f64>, element_id: usize) -> f64 {
    let element = model.elements().iter().find(|element| element.id() == element_id).expect("element should exist");
    let Element2D::TriangleT6(triangle) = element else {
        panic!("T6 plate benchmark should contain only T6 elements");
    };
    let mut max_sigma_y = f64::NEG_INFINITY;

    for (xi, eta) in t6_corner_points() {
        let response = recover_triangle_t6_response(model, triangle, displacements, xi, eta)
            .expect("T6 stress should be recovered");

        max_sigma_y = max_sigma_y.max(response.stress()[1]);
    }

    max_sigma_y
}

fn q8_sigma_y(model: &Model2D, displacements: &nalgebra::DVector<f64>, element_id: usize, xi: f64, eta: f64) -> f64 {
    let element = model.elements().iter().find(|element| element.id() == element_id).expect("element should exist");
    let Element2D::QuadQ8(quad) = element else {
        panic!("plate benchmark should contain only Q8 elements");
    };
    let response =
        recover_quad_q8_response(model, quad, displacements, xi, eta).expect("Q8 stress should be recovered");

    response.stress()[1]
}

fn validate_q8_corner_jacobians(model: &Model2D) {
    for element in model.elements() {
        let Element2D::QuadQ8(_) = element else {
            continue;
        };
        let nodes = q8_nodes(model, element);

        for (xi, eta) in [(-1.0, -1.0), (1.0, -1.0), (1.0, 1.0), (-1.0, 1.0)] {
            let jacobian_determinant = q8_jacobian_determinant(nodes, xi, eta);

            assert!(
                jacobian_determinant.is_finite() && jacobian_determinant > 0.0,
                "invalid Q8 corner Jacobian for element {} at ({xi}, {eta}): {jacobian_determinant}",
                element.id()
            );
        }
    }
}

fn q8_jacobian_determinant(nodes: [&Node2D; 8], xi: f64, eta: f64) -> f64 {
    let derivatives = quad_q8_shape_function_derivatives(xi, eta);
    let dndxi = derivatives[0];
    let dndeta = derivatives[1];
    let mut dx_dxi = 0.0;
    let mut dx_deta = 0.0;
    let mut dy_dxi = 0.0;
    let mut dy_deta = 0.0;

    for node_index in 0..8 {
        dx_dxi += dndxi[node_index] * nodes[node_index].x();
        dx_deta += dndeta[node_index] * nodes[node_index].x();
        dy_dxi += dndxi[node_index] * nodes[node_index].y();
        dy_deta += dndeta[node_index] * nodes[node_index].y();
    }

    dx_dxi * dy_deta - dx_deta * dy_dxi
}

fn q8_nodes<'a>(model: &'a Model2D, element: &Element2D) -> [&'a Node2D; 8] {
    let node_ids = element.node_ids();

    [
        find_node(model, node_ids[0]),
        find_node(model, node_ids[1]),
        find_node(model, node_ids[2]),
        find_node(model, node_ids[3]),
        find_node(model, node_ids[4]),
        find_node(model, node_ids[5]),
        find_node(model, node_ids[6]),
        find_node(model, node_ids[7]),
    ]
}

fn find_node(model: &Model2D, node_id: usize) -> &Node2D {
    model.nodes().iter().find(|node| node.id() == node_id).expect("element should reference existing nodes")
}

fn write_stress_field_svg(model: &Model2D, values: &[(usize, f64)], path: &str, label: &str, max_value: f64) {
    let options = SvgElementScalarFieldOptions {
        legend_label: label.to_string(),
        min_value: 0.0,
        max_value: Some(max_value.max(1e-12)),
        show_mesh: true,
        show_boundary: true,
        ..Default::default()
    };

    write_model_2d_element_scalar_svg_with_options(model, values, path, &options)
        .unwrap_or_else(|error| panic!("stress field SVG should be written to {path}: {error}"));
}

fn write_stress_patch_svg(model: &Model2D, patches: &[TriangleScalarPatch2D], path: &str, label: &str, max_value: f64) {
    let options = SvgElementScalarFieldOptions {
        legend_label: label.to_string(),
        min_value: 0.0,
        max_value: Some(max_value.max(1e-12)),
        show_mesh: true,
        show_boundary: true,
        ..Default::default()
    };

    write_model_2d_triangle_scalar_patches_svg_with_options(model, patches, path, &options)
        .unwrap_or_else(|error| panic!("stress patch SVG should be written to {path}: {error}"));
}

fn q8_stress_component_fields(model: &Model2D, displacements: &nalgebra::DVector<f64>) -> StressComponentFields {
    let mut sigma_x = Vec::with_capacity(model.elements().len());
    let mut sigma_y = Vec::with_capacity(model.elements().len());
    let mut max_sigma_x: f64 = 0.0;
    let mut max_sigma_y: f64 = 0.0;
    let mut max_von_mises: f64 = 0.0;

    for element in model.elements() {
        let Element2D::QuadQ8(quad) = element else {
            panic!("plate benchmark should contain only Q8 elements");
        };
        let mut element_sigma_x: f64 = 0.0;
        let mut element_sigma_y: f64 = 0.0;

        for (xi, eta) in q8_gauss_points() {
            let response =
                recover_quad_q8_response(model, quad, displacements, xi, eta).expect("Q8 stress should be recovered");
            let stress = response.stress();

            element_sigma_x = element_sigma_x.max(stress[0].max(0.0));
            element_sigma_y = element_sigma_y.max(stress[1].max(0.0));

            max_von_mises = max_von_mises.max(response.von_mises_stress());
        }

        for (xi, eta) in q8_corner_points() {
            let response =
                recover_quad_q8_response(model, quad, displacements, xi, eta).expect("Q8 stress should be recovered");
            let stress = response.stress();

            element_sigma_x = element_sigma_x.max(stress[0].max(0.0));
            element_sigma_y = element_sigma_y.max(stress[1].max(0.0));
        }

        max_sigma_x = max_sigma_x.max(element_sigma_x);
        max_sigma_y = max_sigma_y.max(element_sigma_y);
        sigma_x.push((element.id(), element_sigma_x));
        sigma_y.push((element.id(), element_sigma_y));
    }

    StressComponentFields { sigma_x, sigma_y, max_sigma_x, max_sigma_y, max_von_mises }
}

fn t3_field_output(model: &Model2D, displacements: &nalgebra::DVector<f64>) -> T3FieldOutput {
    let numbering = DofNumbering2D::from_model(model).expect("T3 DOF numbering should be created");
    let mut sigma_x = Vec::with_capacity(model.elements().len());
    let mut sigma_y = Vec::with_capacity(model.elements().len());
    let mut von_mises = Vec::with_capacity(model.elements().len());
    let mut displacement = Vec::with_capacity(model.elements().len());
    let mut strain = Vec::with_capacity(model.elements().len());
    let mut max_sigma_x: f64 = 0.0;
    let mut max_sigma_y: f64 = 0.0;
    let mut max_von_mises: f64 = 0.0;
    let mut max_displacement: f64 = 0.0;
    let mut max_strain: f64 = 0.0;

    for element in model.elements() {
        let Element2D::TriangleT3(triangle) = element else {
            panic!("T3 plate benchmark should contain only T3 elements");
        };
        let response =
            recover_triangle_response(model, triangle, displacements).expect("T3 stress should be recovered");
        let stress = response.stress();
        let strain_values = response.strain();
        let element_sigma_x = stress[0].max(0.0);
        let element_sigma_y = stress[1].max(0.0);
        let element_von_mises = response.von_mises_stress();
        let element_displacement = element_max_displacement_magnitude(&numbering, displacements, element);
        let element_strain = vector_magnitude(strain_values);

        max_sigma_x = max_sigma_x.max(element_sigma_x);
        max_sigma_y = max_sigma_y.max(element_sigma_y);
        max_von_mises = max_von_mises.max(element_von_mises);
        max_displacement = max_displacement.max(element_displacement);
        max_strain = max_strain.max(element_strain);

        sigma_x.push((element.id(), element_sigma_x));
        sigma_y.push((element.id(), element_sigma_y));
        von_mises.push((element.id(), element_von_mises));
        displacement.push((element.id(), element_displacement));
        strain.push((element.id(), element_strain));
    }

    T3FieldOutput {
        sigma_x,
        sigma_y,
        von_mises,
        displacement,
        strain,
        max_sigma_x,
        max_sigma_y,
        max_von_mises,
        max_displacement,
        max_strain,
    }
}

fn t6_element_response_data<'a>(
    model: &'a Model2D, displacements: &nalgebra::DVector<f64>,
) -> Vec<T6ElementResponseData<'a>> {
    let numbering = DofNumbering2D::from_model(model).expect("T6 DOF numbering should be created");

    assert_eq!(displacements.len(), numbering.count(), "T6 displacement vector should match the model DOF count",);

    model
        .elements()
        .iter()
        .map(|element| {
            let Element2D::TriangleT6(triangle) = element else {
                panic!("T6 plate benchmark should contain only T6 elements");
            };
            let material = model.material(triangle.material_id()).expect("T6 material should exist");
            model.plane_stress_section(triangle.section_id()).expect("T6 plane-stress section should exist");
            let node_ids = element.node_ids();
            let node_ids = [node_ids[0], node_ids[1], node_ids[2], node_ids[3], node_ids[4], node_ids[5]];
            let indices = numbering.element_dof_indices(element).expect("T6 element DOFs should be numbered");
            let element_displacements = [
                displacements[indices[0]],
                displacements[indices[1]],
                displacements[indices[2]],
                displacements[indices[3]],
                displacements[indices[4]],
                displacements[indices[5]],
                displacements[indices[6]],
                displacements[indices[7]],
                displacements[indices[8]],
                displacements[indices[9]],
                displacements[indices[10]],
                displacements[indices[11]],
            ];

            T6ElementResponseData {
                element_id: element.id(),
                node_ids,
                nodes: [
                    find_node(model, node_ids[0]),
                    find_node(model, node_ids[1]),
                    find_node(model, node_ids[2]),
                    find_node(model, node_ids[3]),
                    find_node(model, node_ids[4]),
                    find_node(model, node_ids[5]),
                ],
                element_displacements,
                constitutive_matrix: plane_stress_constitutive_matrix(material),
            }
        })
        .collect()
}

fn t6_von_mises_field_output(elements: &[T6ElementResponseData<'_>]) -> VonMisesFieldOutput {
    let stress_plot_points = t6_stress_plot_points();
    let mut max_von_mises: f64 = 0.0;

    for element in elements {
        let mut element_von_mises: f64 = 0.0;

        for &(xi, eta) in &stress_plot_points {
            element_von_mises = element_von_mises.max(von_mises_stress(t6_stress(element, xi, eta)));
        }

        max_von_mises = max_von_mises.max(element_von_mises);
    }

    VonMisesFieldOutput { max_von_mises }
}

fn t6_von_mises_gradient_patches(
    elements: &[T6ElementResponseData<'_>], subdivisions: usize,
) -> Vec<TriangleScalarPatch2D> {
    assert!(subdivisions > 0);

    let natural_patches = t6_natural_subtriangles(subdivisions);
    let mut patches = Vec::with_capacity(elements.len() * natural_patches.len());

    for element in elements {
        for natural_patch in &natural_patches {
            let centroid_xi = (natural_patch[0].0 + natural_patch[1].0 + natural_patch[2].0) / 3.0;
            let centroid_eta = (natural_patch[0].1 + natural_patch[1].1 + natural_patch[2].1) / 3.0;
            let points = [
                t6_physical_point(element, natural_patch[0].0, natural_patch[0].1),
                t6_physical_point(element, natural_patch[1].0, natural_patch[1].1),
                t6_physical_point(element, natural_patch[2].0, natural_patch[2].1),
            ];
            let value = von_mises_stress(t6_stress(element, centroid_xi, centroid_eta));

            patches.push(TriangleScalarPatch2D { element_id: element.element_id, points, value });
        }
    }

    patches
}

fn t6_natural_subtriangles(subdivisions: usize) -> Vec<[(f64, f64); 3]> {
    assert!(subdivisions > 0);

    let step = 1.0 / subdivisions as f64;
    let mut subtriangles = Vec::with_capacity(subdivisions * subdivisions);

    for i in 0..subdivisions {
        for j in 0..(subdivisions - i) {
            let p00 = (i as f64 * step, j as f64 * step);
            let p10 = ((i + 1) as f64 * step, j as f64 * step);
            let p01 = (i as f64 * step, (j + 1) as f64 * step);

            subtriangles.push([p00, p10, p01]);

            if i + j + 1 < subdivisions {
                let p11 = ((i + 1) as f64 * step, (j + 1) as f64 * step);

                subtriangles.push([p10, p11, p01]);
            }
        }
    }

    subtriangles
}

fn t6_physical_point(element: &T6ElementResponseData<'_>, xi: f64, eta: f64) -> (f64, f64) {
    let shape_functions = triangle_t6_shape_functions(xi, eta);
    let mut x = 0.0;
    let mut y = 0.0;

    for (node, shape_function) in element.nodes.iter().zip(shape_functions) {
        x += shape_function * node.x();
        y += shape_function * node.y();
    }

    (x, y)
}

fn t6_stress(element: &T6ElementResponseData<'_>, xi: f64, eta: f64) -> [f64; 3] {
    let strain_displacement_matrix = t6_strain_displacement_matrix(element, xi, eta);
    let mut strain = [0.0; 3];
    let mut stress = [0.0; 3];

    for (row, strain_value) in strain.iter_mut().enumerate() {
        *strain_value = strain_displacement_matrix[row]
            .iter()
            .zip(element.element_displacements)
            .map(|(coefficient, displacement)| coefficient * displacement)
            .sum();
    }

    for (row, stress_value) in stress.iter_mut().enumerate() {
        *stress_value = element.constitutive_matrix[row]
            .iter()
            .zip(strain)
            .map(|(coefficient, strain_value)| coefficient * strain_value)
            .sum();
    }

    stress
}

fn t6_strain_displacement_matrix(element: &T6ElementResponseData<'_>, xi: f64, eta: f64) -> [[f64; 12]; 3] {
    let derivatives = triangle_t6_shape_function_derivatives(xi, eta);
    let dndxi = derivatives[0];
    let dndeta = derivatives[1];
    let mut dx_dxi = 0.0;
    let mut dx_deta = 0.0;
    let mut dy_dxi = 0.0;
    let mut dy_deta = 0.0;

    for node_index in 0..6 {
        dx_dxi += dndxi[node_index] * element.nodes[node_index].x();
        dx_deta += dndeta[node_index] * element.nodes[node_index].x();
        dy_dxi += dndxi[node_index] * element.nodes[node_index].y();
        dy_deta += dndeta[node_index] * element.nodes[node_index].y();
    }

    let jacobian_determinant = dx_dxi * dy_deta - dx_deta * dy_dxi;

    assert!(
        jacobian_determinant.is_finite() && jacobian_determinant > 0.0,
        "invalid T6 Jacobian for element {} with nodes {:?}: {jacobian_determinant}",
        element.element_id,
        element.node_ids,
    );

    let mut dndx = [0.0; 6];
    let mut dndy = [0.0; 6];

    for node_index in 0..6 {
        dndx[node_index] = (dy_deta * dndxi[node_index] - dy_dxi * dndeta[node_index]) / jacobian_determinant;
        dndy[node_index] = (-dx_deta * dndxi[node_index] + dx_dxi * dndeta[node_index]) / jacobian_determinant;
    }

    let mut matrix = [[0.0; 12]; 3];

    for node_index in 0..6 {
        let x_dof = 2 * node_index;
        let y_dof = x_dof + 1;

        matrix[0][x_dof] = dndx[node_index];
        matrix[1][y_dof] = dndy[node_index];
        matrix[2][x_dof] = dndy[node_index];
        matrix[2][y_dof] = dndx[node_index];
    }

    matrix
}

fn plane_stress_constitutive_matrix(material: &Material2D) -> [[f64; 3]; 3] {
    let constitutive_factor = material.young_modulus() / (1.0 - material.poisson_ratio() * material.poisson_ratio());

    [
        [constitutive_factor, constitutive_factor * material.poisson_ratio(), 0.0],
        [constitutive_factor * material.poisson_ratio(), constitutive_factor, 0.0],
        [0.0, 0.0, constitutive_factor * (1.0 - material.poisson_ratio()) / 2.0],
    ]
}

fn von_mises_stress(stress: [f64; 3]) -> f64 {
    let [sigma_x, sigma_y, tau_xy] = stress;

    (sigma_x.powi(2) - sigma_x * sigma_y + sigma_y.powi(2) + 3.0 * tau_xy.powi(2)).sqrt()
}

fn t6_stress_plot_points() -> Vec<(f64, f64)> {
    let mut points = t6_gauss_points().iter().map(|(xi, eta, _)| (*xi, *eta)).collect::<Vec<_>>();

    points.extend(t6_corner_points());
    points
}

fn t6_gauss_points() -> [(f64, f64, f64); 6] {
    let vertex = 0.816_847_572_980_459;
    let near_edge = 0.091_576_213_509_771;
    let near_vertex = 0.108_103_018_168_070;
    let mid = 0.445_948_490_915_965;
    let vertex_weight = 0.054_975_871_827_661;
    let mid_weight = 0.111_690_794_839_005;

    [
        (near_edge, near_edge, vertex_weight),
        (vertex, near_edge, vertex_weight),
        (near_edge, vertex, vertex_weight),
        (mid, mid, mid_weight),
        (near_vertex, mid, mid_weight),
        (mid, near_vertex, mid_weight),
    ]
}

fn t6_corner_points() -> [(f64, f64); 3] {
    [(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)]
}

fn element_max_displacement_magnitude(
    numbering: &DofNumbering2D, displacements: &nalgebra::DVector<f64>, element: &Element2D,
) -> f64 {
    element
        .node_ids()
        .iter()
        .map(|node_id| {
            let ux_index = numbering.index(*node_id, Dof2D::Ux).expect("element node should have Ux DOF");
            let uy_index = numbering.index(*node_id, Dof2D::Uy).expect("element node should have Uy DOF");

            displacements[ux_index].hypot(displacements[uy_index])
        })
        .fold(0.0, f64::max)
}

fn vector_magnitude(values: &[f64; 3]) -> f64 {
    values.iter().map(|value| value.powi(2)).sum::<f64>().sqrt()
}

fn t3_strain_gauge_sampling(model: &Model2D, displacements: &nalgebra::DVector<f64>) -> StrainGaugeSamplingOutput {
    let samples = STRAIN_GAUGE_SAMPLES
        .iter()
        .copied()
        .map(|definition| {
            let (x, y) = strain_gauge_sample_point(definition);
            let (sampled_element_count, stress) = t3_stress_at_point(model, displacements, x, y);
            let (experimental_sigma_x, experimental_sigma_y) = experimental_stress_from_strain_gauge(definition);
            let sigma_x_error_percent = relative_error_percent(stress[0], experimental_sigma_x);
            let sigma_y_error_percent = relative_error_percent(stress[1], experimental_sigma_y);

            StrainGaugeSampleResult {
                notch: definition.notch,
                distance: definition.distance,
                x,
                y,
                sampled_element_count,
                fem_sigma_x: stress[0],
                fem_sigma_y: stress[1],
                experimental_sigma_x,
                experimental_sigma_y,
                sigma_x_error_percent,
                sigma_y_error_percent,
            }
        })
        .collect::<Vec<_>>();
    let mean_abs_sigma_y_error_percent =
        samples.iter().map(|sample| sample.sigma_y_error_percent.abs()).sum::<f64>() / samples.len() as f64;
    let max_abs_sigma_y_error_percent =
        samples.iter().map(|sample| sample.sigma_y_error_percent.abs()).fold(0.0, f64::max);

    StrainGaugeSamplingOutput { samples, mean_abs_sigma_y_error_percent, max_abs_sigma_y_error_percent }
}

fn strain_gauge_sample_point(definition: StrainGaugeSampleDefinition) -> (f64, f64) {
    match definition.notch {
        "A_R25" => (-NOTCH_DEPTH - definition.distance, 0.0),
        "B_R50" => (NOTCH_DEPTH + definition.distance, 0.0),
        _ => panic!("unknown strain-gauge notch label: {}", definition.notch),
    }
}

fn experimental_stress_from_strain_gauge(definition: StrainGaugeSampleDefinition) -> (f64, f64) {
    let epsilon_y = definition.epsilon_y_per_mille * 1e-3;
    let epsilon_x = definition.epsilon_x_per_mille * 1e-3;
    let coefficient = YOUNG_MODULUS / (1.0 - POISSON_RATIO.powi(2));
    let sigma_x = coefficient * (epsilon_x + POISSON_RATIO * epsilon_y);
    let sigma_y = coefficient * (epsilon_y + POISSON_RATIO * epsilon_x);

    (sigma_x, sigma_y)
}

fn t3_stress_at_point(model: &Model2D, displacements: &nalgebra::DVector<f64>, x: f64, y: f64) -> (usize, [f64; 3]) {
    let mut stress = [0.0; 3];
    let mut count = 0;

    for element in model.elements() {
        let Element2D::TriangleT3(triangle) = element else {
            panic!("T3 plate benchmark should contain only T3 elements");
        };

        if t3_contains_point(model, element, x, y) {
            let response =
                recover_triangle_response(model, triangle, displacements).expect("T3 stress should be recovered");
            let element_stress = response.stress();

            for index in 0..3 {
                stress[index] += element_stress[index];
            }

            count += 1;
        }
    }

    assert!(count > 0, "strain-gauge sample point ({x}, {y}) should lie inside the T3 mesh");

    for value in &mut stress {
        *value /= count as f64;
    }

    (count, stress)
}

fn t3_contains_point(model: &Model2D, element: &Element2D, x: f64, y: f64) -> bool {
    let node_ids = element.node_ids();
    let first = find_node(model, node_ids[0]);
    let second = find_node(model, node_ids[1]);
    let third = find_node(model, node_ids[2]);
    let denominator =
        (second.y() - third.y()) * (first.x() - third.x()) + (third.x() - second.x()) * (first.y() - third.y());

    if denominator.abs() < 1e-12 {
        return false;
    }

    let alpha = ((second.y() - third.y()) * (x - third.x()) + (third.x() - second.x()) * (y - third.y())) / denominator;
    let beta = ((third.y() - first.y()) * (x - third.x()) + (first.x() - third.x()) * (y - third.y())) / denominator;
    let gamma = 1.0 - alpha - beta;
    let tolerance = 1e-9;

    alpha >= -tolerance
        && beta >= -tolerance
        && gamma >= -tolerance
        && alpha <= 1.0 + tolerance
        && beta <= 1.0 + tolerance
        && gamma <= 1.0 + tolerance
}

fn relative_error_percent(actual: f64, expected: f64) -> f64 {
    100.0 * (actual - expected) / expected.abs().max(1e-12)
}

fn write_strain_gauge_sampling_csv(samples: &[StrainGaugeSampleResult], path: &str) -> std::io::Result<()> {
    let mut csv = String::from(
        "notch,distance_from_notch_mm,x_mm,y_mm,sampled_element_count,fem_sigma_x_mpa,experimental_sigma_x_mpa,sigma_x_error_percent,fem_sigma_y_mpa,experimental_sigma_y_mpa,sigma_y_error_percent\n",
    );

    for sample in samples {
        writeln!(
            csv,
            "{},{:.6},{:.6},{:.6},{},{:.12},{:.12},{:.6},{:.12},{:.12},{:.6}",
            sample.notch,
            sample.distance,
            sample.x,
            sample.y,
            sample.sampled_element_count,
            sample.fem_sigma_x,
            sample.experimental_sigma_x,
            sample.sigma_x_error_percent,
            sample.fem_sigma_y,
            sample.experimental_sigma_y,
            sample.sigma_y_error_percent
        )
        .expect("writing to String should not fail");
    }

    if let Some(parent) = std::path::Path::new(path).parent()
        && !parent.as_os_str().is_empty()
    {
        fs::create_dir_all(parent)?;
    }

    fs::write(path, csv)
}

fn q8_gauss_points() -> [(f64, f64); 9] {
    let point = (3.0_f64 / 5.0).sqrt();

    [
        (-point, -point),
        (0.0, -point),
        (point, -point),
        (-point, 0.0),
        (0.0, 0.0),
        (point, 0.0),
        (-point, point),
        (0.0, point),
        (point, point),
    ]
}

fn q8_corner_points() -> [(f64, f64); 4] {
    [(-1.0, -1.0), (1.0, -1.0), (1.0, 1.0), (-1.0, 1.0)]
}

fn assert_relative_error_below(actual: f64, expected: f64, tolerance: f64, label: &str) {
    let relative_error = (actual - expected).abs() / expected.abs();

    assert!(
        relative_error < tolerance,
        "{label} relative error is too large: actual = {actual}, expected = {expected}, relative error = {relative_error}"
    );
}
