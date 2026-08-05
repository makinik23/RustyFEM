//! Integration tests for the public dense and sparse solver APIs.

use approx::assert_relative_eq;
use rusty_fem::analysis::iterative_solver::CgTerminationReason;
use rusty_fem::analysis::solver::{AnalysisResult2D, solve, solve_sparse};
use rusty_fem::elements::{Beam2D, Element2D, TriangleT3};
use rusty_fem::model::{
    BeamSection2D, DEFAULT_MATERIAL_ID, DisplacementConstraint2D, Dof2D, Material2D, Model2D, NodalLoad2D, Node2D,
    PlaneStressSection2D, Section2D,
};

#[test]
fn dense_and_sparse_solvers_agree_for_a_cantilever_beam() {
    let model = cantilever_beam_model();

    let dense_result = solve(&model).expect("dense beam system should be solvable");
    let sparse_result = solve_sparse(&model).expect("sparse beam system should be solvable");

    assert!(dense_result.solver_report().is_none());
    assert_sparse_solver_converged(&sparse_result);
    assert_results_are_equal(&dense_result, &sparse_result);
}

#[test]
fn dense_and_sparse_solvers_agree_for_a_t3_mesh() {
    let model = two_triangle_cantilever_model();

    let dense_result = solve(&model).expect("dense T3 system should be solvable");
    let sparse_result = solve_sparse(&model).expect("sparse T3 system should be solvable");

    assert!(dense_result.solver_report().is_none());
    assert_sparse_solver_converged(&sparse_result);
    assert_results_are_equal(&dense_result, &sparse_result);
}

fn assert_sparse_solver_converged(result: &AnalysisResult2D) {
    let report = result.solver_report().expect("sparse result should contain a solver report");

    assert_eq!(report.termination_reason, CgTerminationReason::Converged);
    assert!(report.iterations > 0);
    assert!(report.relative_residual_norm < 1e-8);
}

/// Builds a two-node Euler-Bernoulli cantilever beam.
fn cantilever_beam_model() -> Model2D {
    let mut model = Model2D::new();

    model.set_material(Material2D::new(200.0, 0.3, 1.0).expect("test material should be valid"));

    model
        .add_node(Node2D::new(1, 0.0, 0.0).expect("first test node should be valid"))
        .expect("first test node should be added");

    model
        .add_node(Node2D::new(2, 1.0, 0.0).expect("second test node should be valid"))
        .expect("second test node should be added");

    let beam = Beam2D::new(1, [1, 2], DEFAULT_MATERIAL_ID, 100).expect("test beam should be valid");
    let section = Section2D::Beam(BeamSection2D::new(1.0, 2.0).expect("test section should be valid"));

    model.add_element_with_section(Element2D::Beam(beam), section).expect("test beam should be added");

    for dof in [Dof2D::Ux, Dof2D::Uy, Dof2D::Rz] {
        let constraint = DisplacementConstraint2D::new(1, dof, 0.0).expect("beam constraint should be valid");

        model.add_constraint(constraint).expect("beam constraint should be added");
    }

    let load = NodalLoad2D::new(2, Dof2D::Uy, -12.0).expect("beam load should be valid");

    model.add_load(load).expect("beam load should be added");

    model
}

/// Builds a rectangular cantilever consisting of two T3 triangles.
fn two_triangle_cantilever_model() -> Model2D {
    let mut model = Model2D::new();

    model.set_material(Material2D::new(1_000.0, 0.3, 1.0).expect("test material should be valid"));

    for (id, x, y) in [(1, 0.0, 0.0), (2, 1.0, 0.0), (3, 1.0, 1.0), (4, 0.0, 1.0)] {
        model.add_node(Node2D::new(id, x, y).expect("test node should be valid")).expect("test node should be added");
    }

    let lower_right_triangle =
        TriangleT3::new(1, [1, 2, 3], DEFAULT_MATERIAL_ID, 100).expect("first test triangle should be valid");
    let lower_right_section =
        Section2D::PlaneStress(PlaneStressSection2D::new(1.0).expect("first test section should be valid"));

    model
        .add_element_with_section(Element2D::TriangleT3(lower_right_triangle), lower_right_section)
        .expect("first test triangle should be added");

    let upper_left_triangle =
        TriangleT3::new(2, [1, 3, 4], DEFAULT_MATERIAL_ID, 200).expect("second test triangle should be valid");
    let upper_left_section =
        Section2D::PlaneStress(PlaneStressSection2D::new(1.0).expect("second test section should be valid"));

    model
        .add_element_with_section(Element2D::TriangleT3(upper_left_triangle), upper_left_section)
        .expect("second test triangle should be added");

    for node_id in [1, 4] {
        for dof in [Dof2D::Ux, Dof2D::Uy] {
            let constraint = DisplacementConstraint2D::new(node_id, dof, 0.0).expect("T3 constraint should be valid");

            model.add_constraint(constraint).expect("T3 constraint should be added");
        }
    }

    let load = NodalLoad2D::new(3, Dof2D::Uy, -1.0).expect("T3 load should be valid");

    model.add_load(load).expect("T3 load should be added");

    model
}

/// Compares displacements and reactions produced by both public solvers.
fn assert_results_are_equal(dense: &AnalysisResult2D, sparse: &AnalysisResult2D) {
    assert_eq!(dense.displacements().len(), sparse.displacements().len());
    assert_eq!(dense.reactions().len(), sparse.reactions().len());

    for (dense, sparse) in dense.displacements().iter().zip(sparse.displacements().iter()) {
        assert_relative_eq!(dense, sparse, epsilon = 1e-8);
    }

    for (dense, sparse) in dense.reactions().iter().zip(sparse.reactions().iter()) {
        assert_relative_eq!(dense, sparse, epsilon = 1e-8);
    }
}
