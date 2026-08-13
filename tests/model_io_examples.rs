use rusty_fem::analysis::solver::solve_with_settings;
use rusty_fem::io::{AnalysisResult2DOutput, Model2DInput};
use rusty_fem::model::{Dof2D, DofNumbering2D, SolverKind2D};

#[test]
fn t3_cantilever_json_example_solves_to_expected_tip_displacement() {
    let input: Model2DInput =
        serde_json::from_str(include_str!("../examples/t3_cantilever.json")).expect("example JSON should parse");
    let model = input.into_model().expect("example model should be built");

    assert_eq!(model.analysis_settings().solver(), SolverKind2D::Sparse);
    assert_eq!(model.nodes().len(), 15);
    assert_eq!(model.elements().len(), 16);
    assert_eq!(model.constraints().len(), 6);
    assert_eq!(model.element_loads().len(), 2);

    let result = solve_with_settings(&model).expect("example model should solve");
    let numbering = DofNumbering2D::from_model(&model).expect("DOF numbering should be created");
    let middle_right_uy =
        result.displacements()[numbering.index(10, Dof2D::Uy).expect("middle right node should have Uy")];

    assert!(
        (middle_right_uy - -0.473711211345).abs() < 1e-10,
        "unexpected middle-right displacement: {middle_right_uy}"
    );

    let output =
        AnalysisResult2DOutput::from_model_and_result(&model, &result).expect("example output should be created");
    let serialized = serde_json::to_string(&output).expect("example output should serialize");

    assert!(serialized.contains(r#""solver_report""#));
    assert_eq!(output.displacements.len(), 30);
    assert_eq!(output.reactions.len(), 30);
    assert_eq!(output.element_responses.len(), 16);
}

#[test]
#[ignore]
fn t3_and_q4_64x32_json_files_solve_with_sparse() {
    let nx = 64;
    let ny = 32;
    let analytical_tip_deflection = -4.0;
    let t3 = solve_cantilever_json_file(include_str!("../examples/t3_cantilever_64x32.json"), nx, ny);
    let q4 = solve_cantilever_json_file(include_str!("../examples/q4_cantilever_64x32.json"), nx, ny);
    let t3_relative_error = (t3.middle_right_uy - analytical_tip_deflection).abs() / analytical_tip_deflection.abs();
    let q4_relative_error = (q4.middle_right_uy - analytical_tip_deflection).abs() / analytical_tip_deflection.abs();

    eprintln!(
        "64x32 JSON sparse cantilever: T3 uy = {:.12} (rel err {:.6}), Q4 uy = {:.12} (rel err {:.6})",
        t3.middle_right_uy, t3_relative_error, q4.middle_right_uy, q4_relative_error
    );

    assert_eq!(t3.node_count, 2145);
    assert_eq!(q4.node_count, 2145);
    assert_eq!(t3.element_count, 4096);
    assert_eq!(q4.element_count, 2048);
    assert_eq!(t3.displacement_count, 4290);
    assert_eq!(q4.displacement_count, 4290);
    assert_eq!(t3.reaction_count, 4290);
    assert_eq!(q4.reaction_count, 4290);
    assert_eq!(t3.element_response_count, t3.element_count);
    assert_eq!(q4.element_response_count, q4.element_count);
    assert!(t3.middle_right_uy < 0.0);
    assert!(q4.middle_right_uy < 0.0);
    assert!(
        t3_relative_error < 0.03,
        "64x32 T3 JSON file should stay within 3% of the beam reference: uy = {}, error = {}",
        t3.middle_right_uy,
        t3_relative_error
    );
    assert!(
        q4_relative_error < 0.01,
        "64x32 Q4 JSON file should stay within 1% of the beam reference: uy = {}, error = {}",
        q4.middle_right_uy,
        q4_relative_error
    );
    assert!(
        q4_relative_error < t3_relative_error,
        "Q4 should be closer than T3 on the same 64x32 JSON grid: T3 error = {t3_relative_error}, Q4 error = {q4_relative_error}"
    );
}

struct CantileverJsonResult {
    node_count: usize,
    element_count: usize,
    displacement_count: usize,
    reaction_count: usize,
    element_response_count: usize,
    middle_right_uy: f64,
}

fn solve_cantilever_json_file(contents: &str, nx: usize, ny: usize) -> CantileverJsonResult {
    let input: Model2DInput = serde_json::from_str(contents).expect("cantilever JSON should parse");
    let model = input.into_model().expect("cantilever model should be built");

    assert_eq!(model.analysis_settings().solver(), SolverKind2D::Sparse);
    assert_eq!(model.element_loads().len(), ny);

    let result = solve_with_settings(&model).expect("sparse cantilever model should solve");
    let report = result.solver_report().expect("sparse solve should report diagnostics");

    assert!(report.relative_residual_norm < 1e-8);

    let numbering = DofNumbering2D::from_model(&model).expect("DOF numbering should be created");
    let middle_right_node = (ny / 2) * (nx + 1) + nx + 1;
    let middle_right_uy = result.displacements()
        [numbering.index(middle_right_node, Dof2D::Uy).expect("middle right node should have Uy")];
    let output =
        AnalysisResult2DOutput::from_model_and_result(&model, &result).expect("cantilever output should be created");
    let output_json = serde_json::to_string(&output).expect("cantilever output should serialize");

    assert!(output_json.contains(r#""element_responses""#));

    CantileverJsonResult {
        node_count: model.nodes().len(),
        element_count: model.elements().len(),
        displacement_count: output.displacements.len(),
        reaction_count: output.reactions.len(),
        element_response_count: output.element_responses.len(),
        middle_right_uy,
    }
}
