//! Applies prescribed displacement constraints to a linear FEM system.

use nalgebra::{DMatrix, DVector};

use crate::error::FemError;

/// Applies prescribed displacements by modifying the stiffness matrix and load vector.
pub fn apply_displacement_constraints(
    stiffness_matrix: &mut DMatrix<f64>, load_vector: &mut DVector<f64>, constraints: &[(usize, f64)],
) -> Result<(), FemError> {
    let size = stiffness_matrix.nrows();

    if stiffness_matrix.ncols() != size || load_vector.len() != size {
        return Err(FemError::IncompatibleLinearSystem {
            stiffness_rows: stiffness_matrix.nrows(),
            stiffness_columns: stiffness_matrix.ncols(),
            load_vector_length: load_vector.len(),
        });
    }

    for &(constrained_index, prescribed_displacement) in constraints {
        if constrained_index >= size {
            return Err(FemError::InvalidDofIndex { index: constrained_index });
        }

        for row in 0..size {
            if row != constrained_index {
                load_vector[row] -= stiffness_matrix[(row, constrained_index)] * prescribed_displacement;
            }
        }

        for index in 0..size {
            stiffness_matrix[(constrained_index, index)] = 0.0;
            stiffness_matrix[(index, constrained_index)] = 0.0;
        }

        stiffness_matrix[(constrained_index, constrained_index)] = 1.0;
        load_vector[constrained_index] = prescribed_displacement;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::apply_displacement_constraints;
    use crate::error::FemError;
    use nalgebra::{DMatrix, DVector};

    #[test]
    fn applies_zero_displacement_constraint() {
        let mut stiffness_matrix = DMatrix::from_row_slice(2, 2, &[2.0, 1.0, 1.0, 2.0]);
        let mut load_vector = DVector::from_row_slice(&[5.0, 4.0]);

        apply_displacement_constraints(&mut stiffness_matrix, &mut load_vector, &[(0, 0.0)])
            .expect("constraint should be applied");

        assert_eq!(stiffness_matrix, DMatrix::from_row_slice(2, 2, &[1.0, 0.0, 0.0, 2.0]));
        assert_eq!(load_vector, DVector::from_row_slice(&[0.0, 4.0]));
    }

    #[test]
    fn adjusts_load_vector_for_nonzero_displacement_constraint() {
        let mut stiffness_matrix = DMatrix::from_row_slice(2, 2, &[2.0, 1.0, 1.0, 2.0]);
        let mut load_vector = DVector::from_row_slice(&[5.0, 4.0]);

        apply_displacement_constraints(&mut stiffness_matrix, &mut load_vector, &[(0, 1.0)])
            .expect("constraint should be applied");

        assert_eq!(stiffness_matrix, DMatrix::from_row_slice(2, 2, &[1.0, 0.0, 0.0, 2.0]));
        assert_eq!(load_vector, DVector::from_row_slice(&[1.0, 3.0]));
    }

    #[test]
    fn rejects_incompatible_linear_system_dimensions() {
        let mut stiffness_matrix = DMatrix::from_element(2, 3, 0.0);
        let mut load_vector = DVector::from_element(2, 0.0);

        let result = apply_displacement_constraints(&mut stiffness_matrix, &mut load_vector, &[]);

        assert!(matches!(
            result,
            Err(FemError::IncompatibleLinearSystem { stiffness_rows: 2, stiffness_columns: 3, load_vector_length: 2 })
        ));
    }

    #[test]
    fn rejects_invalid_constraint_index() {
        let mut stiffness_matrix = DMatrix::identity(2, 2);
        let mut load_vector = DVector::zeros(2);

        let result = apply_displacement_constraints(&mut stiffness_matrix, &mut load_vector, &[(2, 0.0)]);

        assert!(matches!(result, Err(FemError::InvalidDofIndex { index: 2 })));
    }
}
