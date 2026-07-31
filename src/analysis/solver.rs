//! Solves constrained linear FEM systems.

use nalgebra::{DMatrix, DVector};

use crate::analysis::assembly::assemble_stiffness_matrix;
use crate::analysis::boundary_conditions::apply_displacement_constraints;
use crate::analysis::load_vector::assemble_load_vector;
use crate::error::FemError;
use crate::model::{DofNumbering2D, Model2D};

/// Contains the displacements and support reactions obtained from an analysis.
#[derive(Debug, Clone, PartialEq)]
pub struct AnalysisResult2D {
    displacements: DVector<f64>,
    reactions: DVector<f64>,
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

    Ok(AnalysisResult2D { displacements, reactions })
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
    if stiffness_matrix.nrows() != stiffness_matrix.ncols() || stiffness_matrix.nrows() != load_vector_length {
        return Err(FemError::IncompatibleLinearSystem {
            stiffness_rows: stiffness_matrix.nrows(),
            stiffness_columns: stiffness_matrix.ncols(),
            load_vector_length,
        });
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::{solve, solve_linear_system};
    use crate::elements::{Beam2D, Element2D, Truss2D};
    use crate::error::FemError;
    use crate::model::{DisplacementConstraint2D, Dof2D, Material2D, Model2D, NodalLoad2D, Node2D};
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
    fn model_solver_requires_material() {
        let model = Model2D::new();

        let result = solve(&model);

        assert!(matches!(result, Err(FemError::MissingMaterial)));
    }

    #[test]
    fn solves_horizontal_truss_with_known_displacement() {
        let mut model = Model2D::new();

        model.set_material(Material2D::new(200.0, 0.3, 1.0).expect("valid material should be created"));

        model.add_node(Node2D::new(1, 0.0, 0.0).expect("valid node should be created")).expect("node should be added");
        model.add_node(Node2D::new(2, 1.0, 0.0).expect("valid node should be created")).expect("node should be added");

        let truss = Truss2D::new(10, [1, 2], 2.0).expect("valid truss should be created");

        model.add_element(Element2D::Truss(truss)).expect("element should be added");

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
    fn solves_cantilever_beam_with_known_tip_deflection() {
        let mut model = Model2D::new();

        model.set_material(Material2D::new(200.0, 0.3, 1.0).expect("valid material should be created"));

        model.add_node(Node2D::new(1, 0.0, 0.0).expect("valid node should be created")).expect("node should be added");
        model.add_node(Node2D::new(2, 1.0, 0.0).expect("valid node should be created")).expect("node should be added");

        let beam = Beam2D::new(10, [1, 2], 1.0, 2.0).expect("valid beam should be created");

        model.add_element(Element2D::Beam(beam)).expect("element should be added");

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
}
