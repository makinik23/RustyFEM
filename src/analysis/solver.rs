//! Solves constrained linear FEM systems.

use nalgebra::{DMatrix, DVector};

use crate::analysis::assembly::{assemble_sparse_stiffness_matrix, assemble_stiffness_matrix};
use crate::analysis::boundary_conditions::{apply_displacement_constraints, apply_displacement_constraints_sparse};
use crate::analysis::iterative_solver::{
    CgOptions, CgTerminationReason, JacobiPreconditioner, preconditioned_conjugate_gradient,
};
use crate::analysis::load_vector::assemble_load_vector;
use crate::analysis::sparse::CsrMatrix;
use crate::error::FemError;
use crate::model::{AnalysisSettings2D, DofNumbering2D, Model2D, SolverKind2D};

/// Summarizes the iterative solver execution.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SolverReport {
    /// Number of iterations performed by the iterative solver.
    pub iterations: usize,

    /// Euclidean norm of the final residual.
    pub residual_norm: f64,

    /// Final residual divided by the norm of the right-hand side.
    pub relative_residual_norm: f64,

    /// Reason why the iterative solver stopped.
    pub termination_reason: CgTerminationReason,
}

/// Contains the displacements, reactions, and optional solver diagnostics.
#[derive(Debug, Clone, PartialEq)]
pub struct AnalysisResult2D {
    displacements: DVector<f64>,
    reactions: DVector<f64>,
    solver_report: Option<SolverReport>,
}

impl AnalysisResult2D {
    /// Returns the displacement vector in global DOF order.
    #[must_use]
    pub fn displacements(&self) -> &DVector<f64> {
        &self.displacements
    }

    /// Returns the residual force vector in global DOF order.
    ///
    /// For a solved system, non-zero entries correspond to constrained DOFs
    /// and represent support reactions.
    #[must_use]
    pub fn reactions(&self) -> &DVector<f64> {
        &self.reactions
    }

    /// Returns iterative-solver diagnostics when a sparse solver was used.
    #[must_use]
    pub fn solver_report(&self) -> Option<&SolverReport> {
        self.solver_report.as_ref()
    }
}

/// Solves a 2D model using the solver selected in its analysis settings.
pub fn solve_with_settings(model: &Model2D) -> Result<AnalysisResult2D, FemError> {
    match model.analysis_settings().solver() {
        SolverKind2D::Dense => solve(model),
        SolverKind2D::Sparse => solve_sparse_with_options(model, cg_options_from_settings(model.analysis_settings())),
    }
}

fn cg_options_from_settings(settings: &AnalysisSettings2D) -> CgOptions {
    CgOptions {
        max_iterations: settings.cg_max_iterations(),
        tolerance: settings.cg_tolerance(),
        stagnation_window: settings.cg_stagnation_window(),
        stagnation_tolerance: settings.cg_stagnation_tolerance(),
    }
}

/// Assembles and solves the constrained linear system for a 2D model.
pub fn solve(model: &Model2D) -> Result<AnalysisResult2D, FemError> {
    let numbering = DofNumbering2D::from_model(model)?;
    let original_stiffness_matrix = assemble_stiffness_matrix(model)?;
    let original_load_vector = assemble_load_vector(model)?;
    let mut constrained_stiffness_matrix = original_stiffness_matrix.clone();
    let mut constrained_load_vector = original_load_vector.clone();
    let constraints = numbering.constraint_dof_indices(model)?;

    apply_displacement_constraints(&mut constrained_stiffness_matrix, &mut constrained_load_vector, &constraints)?;

    let displacements = solve_linear_system(constrained_stiffness_matrix, constrained_load_vector)?;
    let reactions = calculate_reactions(&original_stiffness_matrix, &original_load_vector, &displacements)?;

    Ok(AnalysisResult2D { displacements, reactions, solver_report: None })
}

/// Calculates residual forces from the original, unconstrained system.
pub fn calculate_reactions(
    stiffness_matrix: &DMatrix<f64>, load_vector: &DVector<f64>, displacements: &DVector<f64>,
) -> Result<DVector<f64>, FemError> {
    validate_linear_system_dimensions(stiffness_matrix, load_vector.len())?;

    if displacements.len() != stiffness_matrix.nrows() {
        return Err(FemError::IncompatibleLinearSystem {
            stiffness_rows: stiffness_matrix.nrows(),
            stiffness_columns: stiffness_matrix.ncols(),
            load_vector_length: displacements.len(),
        });
    }

    Ok(stiffness_matrix * displacements - load_vector)
}

/// Solves an already constrained linear system.
pub fn solve_linear_system(
    stiffness_matrix: DMatrix<f64>, load_vector: DVector<f64>,
) -> Result<DVector<f64>, FemError> {
    validate_linear_system_dimensions(&stiffness_matrix, load_vector.len())?;

    stiffness_matrix.lu().solve(&load_vector).ok_or(FemError::SingularSystem)
}

fn validate_linear_system_dimensions(
    stiffness_matrix: &DMatrix<f64>, load_vector_length: usize,
) -> Result<(), FemError> {
    if stiffness_matrix.nrows() == 0 && stiffness_matrix.ncols() == 0 && load_vector_length == 0 {
        return Err(FemError::SingularSystem);
    }

    if stiffness_matrix.nrows() != stiffness_matrix.ncols() || stiffness_matrix.nrows() != load_vector_length {
        return Err(FemError::IncompatibleLinearSystem {
            stiffness_rows: stiffness_matrix.nrows(),
            stiffness_columns: stiffness_matrix.ncols(),
            load_vector_length,
        });
    }

    Ok(())
}

/// Assembles and solves a constrained FEM system using sparse matrices.
///
/// The stiffness matrix is assembled as COO, converted to CSR, modified
/// using sparse boundary conditions, and solved with Conjugate Gradient.
pub fn solve_sparse(model: &Model2D) -> Result<AnalysisResult2D, FemError> {
    solve_sparse_with_options(model, CgOptions::default())
}

/// Sparse solver variant with explicit iterative-solver settings.
pub fn solve_sparse_with_options(model: &Model2D, options: CgOptions) -> Result<AnalysisResult2D, FemError> {
    let numbering = DofNumbering2D::from_model(model)?;

    let original_stiffness_matrix = assemble_sparse_stiffness_matrix(model)?;

    let original_load_vector = assemble_load_vector(model)?;

    let mut constrained_stiffness_matrix = original_stiffness_matrix.clone();

    let mut constrained_load_vector = original_load_vector.as_slice().to_vec();

    let constraints = numbering.constraint_dof_indices(model)?;

    apply_displacement_constraints_sparse(
        &mut constrained_stiffness_matrix,
        &mut constrained_load_vector,
        &constraints,
    )?;

    let preconditioner = JacobiPreconditioner::from_matrix(&constrained_stiffness_matrix)?;

    let cg_result = preconditioned_conjugate_gradient(
        &constrained_stiffness_matrix,
        &constrained_load_vector,
        options,
        &preconditioner,
    )?;

    if !cg_result.converged {
        let error = match cg_result.termination_reason {
            CgTerminationReason::MaxIterations => FemError::IterativeSolverDidNotConverge {
                iterations: cg_result.iterations,
                residual_norm: cg_result.residual_norm,
            },
            CgTerminationReason::Stagnated => FemError::IterativeSolverStagnated {
                iterations: cg_result.iterations,
                residual_norm: cg_result.residual_norm,
            },
            CgTerminationReason::Converged => unreachable!("a converged result was marked as non-converged"),
        };

        return Err(error);
    }

    let displacements = DVector::from_vec(cg_result.solution);

    let solver_report = SolverReport {
        iterations: cg_result.iterations,
        residual_norm: cg_result.residual_norm,
        relative_residual_norm: cg_result.relative_residual_norm,
        termination_reason: cg_result.termination_reason,
    };

    // Reactions must use the original, unconstrained system.
    let reactions = calculate_sparse_reactions(
        &original_stiffness_matrix,
        original_load_vector.as_slice(),
        displacements.as_slice(),
    )?;

    Ok(AnalysisResult2D { displacements, reactions, solver_report: Some(solver_report) })
}

/// Calculates residual forces using an original sparse system.
///
/// Reactions are computed before boundary-condition modification:
///
/// ```text
/// reactions = K_original * displacements - f_original
/// ```
pub fn calculate_sparse_reactions(
    stiffness_matrix: &CsrMatrix, load_vector: &[f64], displacements: &[f64],
) -> Result<DVector<f64>, FemError> {
    let size = stiffness_matrix.nrows();

    if stiffness_matrix.ncols() != size {
        return Err(FemError::IncompatibleLinearSystem {
            stiffness_rows: stiffness_matrix.nrows(),
            stiffness_columns: stiffness_matrix.ncols(),
            load_vector_length: load_vector.len(),
        });
    }

    if load_vector.len() != size {
        return Err(FemError::InvalidVectorLength { vector: "load", expected: size, actual: load_vector.len() });
    }

    if displacements.len() != size {
        return Err(FemError::InvalidVectorLength {
            vector: "displacements",
            expected: size,
            actual: displacements.len(),
        });
    }

    let mut internal_forces = vec![0.0; size];
    stiffness_matrix.mul_vector(displacements, &mut internal_forces)?;

    for index in 0..size {
        internal_forces[index] -= load_vector[index];
    }

    Ok(DVector::from_vec(internal_forces))
}

#[cfg(test)]
mod tests {
    use super::{solve, solve_linear_system, solve_with_settings};
    use crate::elements::{Beam2D, Element2D, QuadQ4, QuadQ8, TriangleT3, Truss2D};
    use crate::error::FemError;
    use crate::model::{
        BeamSection2D, BeamUniformLineLoad2D, BodyForce2D, DEFAULT_MATERIAL_ID, DisplacementConstraint2D, Dof2D,
        EdgeTraction2D, ElementLoad2D, LoadCoordinateSystem2D, Material2D, Model2D, NodalLoad2D, Node2D,
        PlaneStressSection2D, Section2D, SelfWeight2D, SolverKind2D, TrussSection2D,
    };
    use approx::assert_relative_eq;
    use nalgebra::{DMatrix, DVector};

    #[test]
    fn solves_an_algebraic_constrained_system() {
        let mut stiffness_matrix = DMatrix::from_row_slice(2, 2, &[2.0, 1.0, 1.0, 2.0]);
        let mut load_vector = DVector::from_row_slice(&[5.0, 4.0]);

        crate::analysis::boundary_conditions::apply_displacement_constraints(
            &mut stiffness_matrix,
            &mut load_vector,
            &[(0, 1.0)],
        )
        .expect("constraint should be applied");

        let solution = solve_linear_system(stiffness_matrix, load_vector).expect("system should be solved");

        assert_eq!(solution, DVector::from_row_slice(&[1.0, 1.5]));
    }

    #[test]
    fn rejects_singular_system() {
        let stiffness_matrix = DMatrix::from_row_slice(2, 2, &[1.0, 0.0, 0.0, 0.0]);
        let load_vector = DVector::from_row_slice(&[1.0, 0.0]);

        let result = solve_linear_system(stiffness_matrix, load_vector);

        assert!(matches!(result, Err(FemError::SingularSystem)));
    }

    #[test]
    fn rejects_incompatible_system_dimensions() {
        let stiffness_matrix = DMatrix::from_element(2, 3, 0.0);
        let load_vector = DVector::zeros(2);

        let result = solve_linear_system(stiffness_matrix, load_vector);

        assert!(matches!(
            result,
            Err(FemError::IncompatibleLinearSystem { stiffness_rows: 2, stiffness_columns: 3, load_vector_length: 2 })
        ));
    }

    #[test]
    fn model_solver_rejects_empty_model() {
        let model = Model2D::new();

        let result = solve(&model);

        assert!(matches!(result, Err(FemError::SingularSystem)));
    }

    #[test]
    fn solves_horizontal_truss_with_known_displacement() {
        let mut model = Model2D::new();

        model.set_material(Material2D::new(200.0, 0.3, 1.0).expect("valid material should be created"));

        model.add_node(Node2D::new(1, 0.0, 0.0).expect("valid node should be created")).expect("node should be added");
        model.add_node(Node2D::new(2, 1.0, 0.0).expect("valid node should be created")).expect("node should be added");

        let truss = Truss2D::new(10, [1, 2], DEFAULT_MATERIAL_ID, 100).expect("valid truss should be created");
        let section = Section2D::Truss(TrussSection2D::new(2.0).expect("valid section should be created"));

        model.add_element_with_section(Element2D::Truss(truss), section).expect("element should be added");

        for (node_id, dof) in [(1, Dof2D::Ux), (1, Dof2D::Uy), (2, Dof2D::Uy)] {
            let constraint =
                DisplacementConstraint2D::new(node_id, dof, 0.0).expect("valid constraint should be created");

            model.add_constraint(constraint).expect("constraint should be added");
        }

        let load = NodalLoad2D::new(2, Dof2D::Ux, 10.0).expect("valid load should be created");
        model.add_load(load).expect("load should be added");

        let result = solve(&model).expect("system should be solved");
        let solution = result.displacements();
        let reactions = result.reactions();

        assert_relative_eq!(solution[0], 0.0, epsilon = 1e-12);
        assert_relative_eq!(solution[1], 0.0, epsilon = 1e-12);
        assert_relative_eq!(solution[2], 0.025, epsilon = 1e-12);
        assert_relative_eq!(solution[3], 0.0, epsilon = 1e-12);

        assert_relative_eq!(reactions[0], -10.0, epsilon = 1e-12);
        assert_relative_eq!(reactions[1], 0.0, epsilon = 1e-12);
        assert_relative_eq!(reactions[2], 0.0, epsilon = 1e-12);
        assert_relative_eq!(reactions[3], 0.0, epsilon = 1e-12);
    }

    #[test]
    fn sparse_solver_matches_dense_solver_on_constrained_truss() {
        let mut model = Model2D::new();

        model.set_material(Material2D::new(200.0, 0.3, 1.0).expect("valid material should be created"));

        model.add_node(Node2D::new(1, 0.0, 0.0).expect("valid node should be created")).expect("node should be added");

        model.add_node(Node2D::new(2, 1.0, 0.0).expect("valid node should be created")).expect("node should be added");

        model
            .add_element_with_section(
                Element2D::Truss(
                    Truss2D::new(10, [1, 2], DEFAULT_MATERIAL_ID, 100).expect("valid truss should be created"),
                ),
                Section2D::Truss(TrussSection2D::new(2.0).expect("valid section should be created")),
            )
            .expect("element should be added");

        for (node_id, dof) in [(1, Dof2D::Ux), (1, Dof2D::Uy), (2, Dof2D::Uy)] {
            model
                .add_constraint(
                    DisplacementConstraint2D::new(node_id, dof, 0.0).expect("valid constraint should be created"),
                )
                .expect("constraint should be added");
        }

        model
            .add_load(NodalLoad2D::new(2, Dof2D::Ux, 10.0).expect("valid load should be created"))
            .expect("load should be added");

        let dense_result = solve(&model).expect("dense system should be solved");

        model.analysis_settings_mut().set_solver(SolverKind2D::Sparse);
        let sparse_result = solve_with_settings(&model).expect("sparse system should be solved");

        assert!(sparse_result.solver_report().is_some());

        for (sparse, dense) in sparse_result.displacements().iter().zip(dense_result.displacements().iter()) {
            assert_relative_eq!(sparse, dense, epsilon = 1e-10);
        }

        for (sparse, dense) in sparse_result.reactions().iter().zip(dense_result.reactions().iter()) {
            assert_relative_eq!(sparse, dense, epsilon = 1e-10);
        }
    }

    #[test]
    fn solves_cantilever_beam_with_known_tip_deflection() {
        let mut model = Model2D::new();

        model.set_material(Material2D::new(200.0, 0.3, 1.0).expect("valid material should be created"));

        model.add_node(Node2D::new(1, 0.0, 0.0).expect("valid node should be created")).expect("node should be added");
        model.add_node(Node2D::new(2, 1.0, 0.0).expect("valid node should be created")).expect("node should be added");

        let beam = Beam2D::new(10, [1, 2], DEFAULT_MATERIAL_ID, 100).expect("valid beam should be created");
        let section = Section2D::Beam(BeamSection2D::new(1.0, 2.0).expect("valid section should be created"));

        model.add_element_with_section(Element2D::Beam(beam), section).expect("element should be added");

        for dof in [Dof2D::Ux, Dof2D::Uy, Dof2D::Rz] {
            let constraint = DisplacementConstraint2D::new(1, dof, 0.0).expect("valid constraint should be created");

            model.add_constraint(constraint).expect("constraint should be added");
        }

        let load = NodalLoad2D::new(2, Dof2D::Uy, -12.0).expect("valid load should be created");
        model.add_load(load).expect("load should be added");

        let result = solve(&model).expect("system should be solved");
        let solution = result.displacements();
        let reactions = result.reactions();

        assert_relative_eq!(solution[0], 0.0, epsilon = 1e-12);
        assert_relative_eq!(solution[1], 0.0, epsilon = 1e-12);
        assert_relative_eq!(solution[2], 0.0, epsilon = 1e-12);
        assert_relative_eq!(solution[3], 0.0, epsilon = 1e-12);
        assert_relative_eq!(solution[4], -0.01, epsilon = 1e-12);
        assert_relative_eq!(solution[5], -0.015, epsilon = 1e-12);

        assert_relative_eq!(reactions[0], 0.0, epsilon = 1e-12);
        assert_relative_eq!(reactions[1], 12.0, epsilon = 1e-12);
        assert_relative_eq!(reactions[2], 12.0, epsilon = 1e-12);
        assert_relative_eq!(reactions[3], 0.0, epsilon = 1e-12);
        assert_relative_eq!(reactions[4], 0.0, epsilon = 1e-12);
        assert_relative_eq!(reactions[5], 0.0, epsilon = 1e-12);
    }

    #[test]
    fn solves_cantilever_beam_with_uniform_line_load() {
        let mut model = Model2D::new();

        model.set_material(Material2D::new(200.0, 0.3, 1.0).expect("valid material should be created"));

        model.add_node(Node2D::new(1, 0.0, 0.0).expect("valid node should be created")).expect("node should be added");
        model.add_node(Node2D::new(2, 1.0, 0.0).expect("valid node should be created")).expect("node should be added");

        let beam = Beam2D::new(10, [1, 2], DEFAULT_MATERIAL_ID, 100).expect("valid beam should be created");
        let section = Section2D::Beam(BeamSection2D::new(1.0, 2.0).expect("valid section should be created"));

        model.add_element_with_section(Element2D::Beam(beam), section).expect("element should be added");

        for dof in [Dof2D::Ux, Dof2D::Uy, Dof2D::Rz] {
            let constraint = DisplacementConstraint2D::new(1, dof, 0.0).expect("valid constraint should be created");

            model.add_constraint(constraint).expect("constraint should be added");
        }

        let load = ElementLoad2D::BeamUniformLine(
            BeamUniformLineLoad2D::new(10, LoadCoordinateSystem2D::Local, 0.0, -12.0)
                .expect("valid load should be created"),
        );
        model.add_element_load(load).expect("load should be added");

        let result = solve(&model).expect("system should be solved");
        let solution = result.displacements();
        let reactions = result.reactions();

        assert_relative_eq!(solution[0], 0.0, epsilon = 1e-12);
        assert_relative_eq!(solution[1], 0.0, epsilon = 1e-12);
        assert_relative_eq!(solution[2], 0.0, epsilon = 1e-12);
        assert_relative_eq!(solution[3], 0.0, epsilon = 1e-12);
        assert_relative_eq!(solution[4], -0.00375, epsilon = 1e-12);
        assert_relative_eq!(solution[5], -0.005, epsilon = 1e-12);

        assert_relative_eq!(reactions[0], 0.0, epsilon = 1e-12);
        assert_relative_eq!(reactions[1], 12.0, epsilon = 1e-12);
        assert_relative_eq!(reactions[2], 6.0, epsilon = 1e-12);
        assert_relative_eq!(reactions[3], 0.0, epsilon = 1e-12);
        assert_relative_eq!(reactions[4], 0.0, epsilon = 1e-12);
        assert_relative_eq!(reactions[5], 0.0, epsilon = 1e-12);
    }

    #[test]
    fn solves_triangle_with_edge_traction_and_balances_reactions() {
        let mut model = Model2D::new();

        model.set_material(Material2D::new(200.0, 0.3, 1.0).expect("valid material should be created"));

        model.add_node(Node2D::new(1, 0.0, 0.0).expect("valid node should be created")).expect("node should be added");
        model.add_node(Node2D::new(2, 1.0, 0.0).expect("valid node should be created")).expect("node should be added");
        model.add_node(Node2D::new(3, 0.0, 1.0).expect("valid node should be created")).expect("node should be added");

        let triangle =
            TriangleT3::new(10, [1, 2, 3], DEFAULT_MATERIAL_ID, 100).expect("valid triangle should be created");
        let section = Section2D::PlaneStress(PlaneStressSection2D::new(1.0).expect("valid section should be created"));

        model.add_element_with_section(Element2D::TriangleT3(triangle), section).expect("element should be added");

        for (node_id, dof) in [(1, Dof2D::Ux), (1, Dof2D::Uy), (2, Dof2D::Uy), (3, Dof2D::Ux)] {
            let constraint =
                DisplacementConstraint2D::new(node_id, dof, 0.0).expect("valid constraint should be created");

            model.add_constraint(constraint).expect("constraint should be added");
        }

        let load = ElementLoad2D::EdgeTraction(
            EdgeTraction2D::new(10, [2, 3], LoadCoordinateSystem2D::Global, 10.0, 0.0)
                .expect("valid load should be created"),
        );
        model.add_element_load(load).expect("load should be added");

        let result = solve(&model).expect("system should be solved");
        let reactions = result.reactions();
        let total_edge_force = 10.0 * 2.0_f64.sqrt();
        let x_reaction_sum = reactions[0] + reactions[4];
        let y_reaction_sum = reactions[1] + reactions[3];

        assert_relative_eq!(x_reaction_sum, -total_edge_force, epsilon = 1e-12);
        assert_relative_eq!(y_reaction_sum, 0.0, epsilon = 1e-12);
    }

    #[test]
    fn solves_triangle_with_body_force_and_balances_reactions() {
        let mut model = Model2D::new();

        model.set_material(Material2D::new(200.0, 0.3, 1.0).expect("valid material should be created"));

        model.add_node(Node2D::new(1, 0.0, 0.0).expect("valid node should be created")).expect("node should be added");
        model.add_node(Node2D::new(2, 1.0, 0.0).expect("valid node should be created")).expect("node should be added");
        model.add_node(Node2D::new(3, 0.0, 1.0).expect("valid node should be created")).expect("node should be added");

        let triangle =
            TriangleT3::new(10, [1, 2, 3], DEFAULT_MATERIAL_ID, 100).expect("valid triangle should be created");
        let section = Section2D::PlaneStress(PlaneStressSection2D::new(1.0).expect("valid section should be created"));

        model.add_element_with_section(Element2D::TriangleT3(triangle), section).expect("element should be added");

        for (node_id, dof) in [(1, Dof2D::Ux), (1, Dof2D::Uy), (2, Dof2D::Uy), (3, Dof2D::Ux)] {
            let constraint =
                DisplacementConstraint2D::new(node_id, dof, 0.0).expect("valid constraint should be created");

            model.add_constraint(constraint).expect("constraint should be added");
        }

        let load = ElementLoad2D::BodyForce(BodyForce2D::new(10, 10.0, 0.0).expect("valid load should be created"));
        model.add_element_load(load).expect("load should be added");

        let result = solve(&model).expect("system should be solved");
        let reactions = result.reactions();
        let total_body_force = 10.0 * 0.5;
        let x_reaction_sum = reactions[0] + reactions[4];
        let y_reaction_sum = reactions[1] + reactions[3];

        assert_relative_eq!(x_reaction_sum, -total_body_force, epsilon = 1e-12);
        assert_relative_eq!(y_reaction_sum, 0.0, epsilon = 1e-12);
    }

    #[test]
    fn solves_triangle_with_self_weight_and_balances_reactions() {
        let mut model = Model2D::new();

        model.set_material(Material2D::new(200.0, 0.3, 2.0).expect("valid material should be created"));

        model.add_node(Node2D::new(1, 0.0, 0.0).expect("valid node should be created")).expect("node should be added");
        model.add_node(Node2D::new(2, 1.0, 0.0).expect("valid node should be created")).expect("node should be added");
        model.add_node(Node2D::new(3, 0.0, 1.0).expect("valid node should be created")).expect("node should be added");

        let triangle =
            TriangleT3::new(10, [1, 2, 3], DEFAULT_MATERIAL_ID, 100).expect("valid triangle should be created");
        let section = Section2D::PlaneStress(PlaneStressSection2D::new(1.0).expect("valid section should be created"));

        model.add_element_with_section(Element2D::TriangleT3(triangle), section).expect("element should be added");

        for (node_id, dof) in [(1, Dof2D::Ux), (1, Dof2D::Uy), (2, Dof2D::Uy), (3, Dof2D::Ux)] {
            let constraint =
                DisplacementConstraint2D::new(node_id, dof, 0.0).expect("valid constraint should be created");

            model.add_constraint(constraint).expect("constraint should be added");
        }

        let load = ElementLoad2D::SelfWeight(SelfWeight2D::new(10, 10.0, 0.0).expect("valid load should be created"));
        model.add_element_load(load).expect("load should be added");

        let result = solve(&model).expect("system should be solved");
        let reactions = result.reactions();
        let total_self_weight_force = 2.0 * 10.0 * 0.5;
        let x_reaction_sum = reactions[0] + reactions[4];
        let y_reaction_sum = reactions[1] + reactions[3];

        assert_relative_eq!(x_reaction_sum, -total_self_weight_force, epsilon = 1e-12);
        assert_relative_eq!(y_reaction_sum, 0.0, epsilon = 1e-12);
    }

    #[test]
    fn solves_quad_with_edge_traction_and_balances_reactions() {
        let mut model = fully_constrained_quad_model(1.0);
        let load = ElementLoad2D::EdgeTraction(
            EdgeTraction2D::new(10, [2, 3], LoadCoordinateSystem2D::Global, 8.0, -4.0)
                .expect("valid load should be created"),
        );
        model.add_element_load(load).expect("load should be added");

        let result = solve(&model).expect("system should be solved");
        let reactions = result.reactions();
        let x_reaction_sum = reactions.iter().step_by(2).sum::<f64>();
        let y_reaction_sum = reactions.iter().skip(1).step_by(2).sum::<f64>();

        assert_relative_eq!(x_reaction_sum, -4.0, epsilon = 1e-12);
        assert_relative_eq!(y_reaction_sum, 2.0, epsilon = 1e-12);
    }

    #[test]
    fn solves_quad_with_body_force_and_balances_reactions() {
        let mut model = fully_constrained_quad_model(1.0);
        let load = ElementLoad2D::BodyForce(BodyForce2D::new(10, 8.0, -4.0).expect("valid load should be created"));
        model.add_element_load(load).expect("load should be added");

        let result = solve(&model).expect("system should be solved");
        let reactions = result.reactions();
        let x_reaction_sum = reactions.iter().step_by(2).sum::<f64>();
        let y_reaction_sum = reactions.iter().skip(1).step_by(2).sum::<f64>();

        assert_relative_eq!(x_reaction_sum, -8.0, epsilon = 1e-12);
        assert_relative_eq!(y_reaction_sum, 4.0, epsilon = 1e-12);
    }

    #[test]
    fn solves_quad_q8_with_body_force_and_balances_reactions() {
        let mut model = fully_constrained_quad_q8_model(1.0);
        let load = ElementLoad2D::BodyForce(BodyForce2D::new(10, 8.0, -4.0).expect("valid load should be created"));
        model.add_element_load(load).expect("load should be added");

        let result = solve(&model).expect("system should be solved");
        let reactions = result.reactions();
        let x_reaction_sum = reactions.iter().step_by(2).sum::<f64>();
        let y_reaction_sum = reactions.iter().skip(1).step_by(2).sum::<f64>();

        assert_relative_eq!(x_reaction_sum, -8.0, epsilon = 1e-12);
        assert_relative_eq!(y_reaction_sum, 4.0, epsilon = 1e-12);
    }

    #[test]
    fn solves_quad_q8_with_edge_traction_and_balances_reactions() {
        let mut model = fully_constrained_quad_q8_model(1.0);
        let load = ElementLoad2D::EdgeTraction(
            EdgeTraction2D::new(10, [2, 3], LoadCoordinateSystem2D::Global, 8.0, -4.0)
                .expect("valid load should be created"),
        );
        model.add_element_load(load).expect("load should be added");

        let result = solve(&model).expect("system should be solved");
        let reactions = result.reactions();
        let x_reaction_sum = reactions.iter().step_by(2).sum::<f64>();
        let y_reaction_sum = reactions.iter().skip(1).step_by(2).sum::<f64>();

        assert_relative_eq!(x_reaction_sum, -4.0, epsilon = 1e-12);
        assert_relative_eq!(y_reaction_sum, 2.0, epsilon = 1e-12);
    }

    #[test]
    fn solves_quad_with_self_weight_and_balances_reactions() {
        let mut model = fully_constrained_quad_model(2.0);
        let load = ElementLoad2D::SelfWeight(SelfWeight2D::new(10, 5.0, -10.0).expect("valid load should be created"));
        model.add_element_load(load).expect("load should be added");

        let result = solve(&model).expect("system should be solved");
        let reactions = result.reactions();
        let x_reaction_sum = reactions.iter().step_by(2).sum::<f64>();
        let y_reaction_sum = reactions.iter().skip(1).step_by(2).sum::<f64>();

        assert_relative_eq!(x_reaction_sum, -10.0, epsilon = 1e-12);
        assert_relative_eq!(y_reaction_sum, 20.0, epsilon = 1e-12);
    }

    fn fully_constrained_quad_model(density: f64) -> Model2D {
        let mut model = Model2D::new();

        model.set_material(Material2D::new(200.0, 0.3, density).expect("valid material should be created"));

        for (id, x, y) in [(1, 0.0, 0.0), (2, 2.0, 0.0), (3, 2.0, 1.0), (4, 0.0, 1.0)] {
            model.add_node(Node2D::new(id, x, y).expect("valid node should be created")).expect("node should be added");
        }

        let quad = QuadQ4::new(10, [1, 2, 3, 4], DEFAULT_MATERIAL_ID, 100).expect("valid quad should be created");
        let section = Section2D::PlaneStress(PlaneStressSection2D::new(0.5).expect("valid section should be created"));

        model.add_element_with_section(Element2D::QuadQ4(quad), section).expect("element should be added");

        for node_id in 1..=4 {
            for dof in [Dof2D::Ux, Dof2D::Uy] {
                model
                    .add_constraint(DisplacementConstraint2D::new(node_id, dof, 0.0).expect("valid constraint"))
                    .expect("constraint should be added");
            }
        }

        model
    }

    fn fully_constrained_quad_q8_model(density: f64) -> Model2D {
        let mut model = Model2D::new();

        model.set_material(Material2D::new(200.0, 0.3, density).expect("valid material should be created"));

        for (id, x, y) in [
            (1, 0.0, 0.0),
            (2, 2.0, 0.0),
            (3, 2.0, 1.0),
            (4, 0.0, 1.0),
            (5, 1.0, 0.0),
            (6, 2.0, 0.5),
            (7, 1.0, 1.0),
            (8, 0.0, 0.5),
        ] {
            model.add_node(Node2D::new(id, x, y).expect("valid node should be created")).expect("node should be added");
        }

        let quad =
            QuadQ8::new(10, [1, 2, 3, 4, 5, 6, 7, 8], DEFAULT_MATERIAL_ID, 100).expect("valid quad should be created");
        let section = Section2D::PlaneStress(PlaneStressSection2D::new(0.5).expect("valid section should be created"));

        model.add_element_with_section(Element2D::QuadQ8(quad), section).expect("element should be added");

        for node_id in 1..=8 {
            for dof in [Dof2D::Ux, Dof2D::Uy] {
                model
                    .add_constraint(DisplacementConstraint2D::new(node_id, dof, 0.0).expect("valid constraint"))
                    .expect("constraint should be added");
            }
        }

        model
    }
}
