//! Iterative solvers for sparse linear systems.

use crate::analysis::sparse::CsrMatrix;
use crate::error::FemError;

/// Describes a linear operator that can compute `y = A * x`.
pub trait LinearOperator {
    /// Returns the dimension of the square linear system.
    fn dimension(&self) -> usize;

    /// Computes `y = A * x`.
    fn apply(&self, x: &[f64], y: &mut [f64]) -> Result<(), FemError>;
}

/// Transforms a residual into a preconditioned residual.
///
/// A preconditioner approximates the inverse of a linear operator without
/// explicitly constructing a full inverse matrix.
pub trait Preconditioner {
    /// Computes `corrected = M⁻¹ * residual`.
    fn apply(&self, residual: &[f64], corrected: &mut [f64]) -> Result<(), FemError>;
}

/// A Jacobi preconditioner based on the inverse diagonal of a CSR matrix.
#[derive(Debug, Clone)]
pub struct JacobiPreconditioner {
    inverse_diagonal: Vec<f64>,
}

impl JacobiPreconditioner {
    /// Builds a Jacobi preconditioner from a square CSR matrix.
    pub fn from_matrix(matrix: &CsrMatrix) -> Result<Self, FemError> {
        if matrix.nrows() != matrix.ncols() {
            return Err(FemError::IncompatibleLinearSystem {
                stiffness_rows: matrix.nrows(),
                stiffness_columns: matrix.ncols(),
                load_vector_length: matrix.nrows(),
            });
        }

        let mut inverse_diagonal = Vec::with_capacity(matrix.nrows());

        for index in 0..matrix.nrows() {
            let diagonal = matrix.value_at(index, index)?;

            if !diagonal.is_finite() || diagonal <= f64::EPSILON {
                return Err(FemError::InvalidPreconditionerDiagonal { index, value: diagonal });
            }

            inverse_diagonal.push(1.0 / diagonal);
        }

        Ok(Self { inverse_diagonal })
    }
}

impl Preconditioner for JacobiPreconditioner {
    fn apply(&self, residual: &[f64], corrected: &mut [f64]) -> Result<(), FemError> {
        if residual.len() != self.inverse_diagonal.len() {
            return Err(FemError::InvalidVectorLength {
                vector: "residual",
                expected: self.inverse_diagonal.len(),
                actual: residual.len(),
            });
        }

        if corrected.len() != self.inverse_diagonal.len() {
            return Err(FemError::InvalidVectorLength {
                vector: "corrected residual",
                expected: self.inverse_diagonal.len(),
                actual: corrected.len(),
            });
        }

        for index in 0..self.inverse_diagonal.len() {
            corrected[index] = self.inverse_diagonal[index] * residual[index];
        }

        Ok(())
    }
}

const NORM_EPSILON: f64 = f64::EPSILON;

/// Configuration of the Conjugate Gradient solver.
#[derive(Debug, Clone, Copy)]
pub struct CgOptions {
    /// Maximum number of iterations.
    pub max_iterations: usize,

    /// Relative convergence tolerance.
    pub tolerance: f64,

    /// Number of consecutive iterations with insufficient improvement that
    /// cause the solver to report stagnation. Set to zero to disable this
    /// check.
    pub stagnation_window: usize,

    /// Minimum relative residual improvement required to reset the
    /// stagnation counter.
    pub stagnation_tolerance: f64,
}

impl Default for CgOptions {
    fn default() -> Self {
        Self { max_iterations: 1_000, tolerance: 1e-10, stagnation_window: 8, stagnation_tolerance: 1e-12 }
    }
}

/// Explains why an iterative solve stopped.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CgTerminationReason {
    /// The requested relative residual tolerance was reached.
    Converged,

    /// The configured iteration limit was reached first.
    MaxIterations,

    /// The residual stopped improving sufficiently.
    Stagnated,
}

/// Result returned by the Conjugate Gradient solver.
#[derive(Debug, Clone)]
pub struct CgResult {
    /// Computed approximate solution.
    pub solution: Vec<f64>,

    /// Number of iterations performed.
    pub iterations: usize,

    /// Euclidean norm of the final residual.
    pub residual_norm: f64,

    /// Final residual norm divided by the norm of the right-hand side.
    pub relative_residual_norm: f64,

    /// Indicates whether the requested tolerance was reached.
    pub converged: bool,

    /// Describes the reason why the iteration stopped.
    pub termination_reason: CgTerminationReason,
}

/// Solves `A * x = b` using the Conjugate Gradient method.
///
/// The matrix must be symmetric and positive definite. This is normally
/// satisfied by a properly constrained FEM stiffness matrix.
pub fn conjugate_gradient<A: LinearOperator>(
    operator: &A, rhs: &[f64], options: CgOptions,
) -> Result<CgResult, FemError> {
    preconditioned_conjugate_gradient(operator, rhs, options, &IdentityPreconditioner)
}

/// Solves `A * x = b` using the preconditioned Conjugate Gradient method.
///
/// The matrix must be symmetric and positive definite. The preconditioner
/// should also preserve the positive-definite structure of the system.
pub fn preconditioned_conjugate_gradient<A, P>(
    operator: &A, rhs: &[f64], options: CgOptions, preconditioner: &P,
) -> Result<CgResult, FemError>
where
    A: LinearOperator,
    P: Preconditioner,
{
    if !options.tolerance.is_finite() || options.tolerance <= 0.0 {
        return Err(FemError::InvalidSolverTolerance { value: options.tolerance });
    }

    if !options.stagnation_tolerance.is_finite() || options.stagnation_tolerance < 0.0 {
        return Err(FemError::InvalidStagnationTolerance { value: options.stagnation_tolerance });
    }

    let dimension = operator.dimension();

    if rhs.len() != dimension {
        return Err(FemError::InvalidVectorLength {
            vector: "right-hand side",
            expected: dimension,
            actual: rhs.len(),
        });
    }

    let mut solution = vec![0.0; dimension];

    let mut residual = rhs.to_vec();
    let mut corrected_residual = vec![0.0; dimension];
    let mut direction = vec![0.0; dimension];
    let mut operator_direction = vec![0.0; dimension];

    let rhs_norm = euclidean_norm(rhs);
    let mut residual_norm = euclidean_norm(&residual);
    let mut relative_residual = relative_residual_norm(residual_norm, rhs_norm);

    if relative_residual <= options.tolerance {
        return Ok(make_result(solution, 0, residual_norm, relative_residual, CgTerminationReason::Converged));
    }

    if options.max_iterations == 0 {
        return Ok(make_result(solution, 0, residual_norm, relative_residual, CgTerminationReason::MaxIterations));
    }

    preconditioner.apply(&residual, &mut corrected_residual)?;
    direction.copy_from_slice(&corrected_residual);

    let mut residual_dot_corrected = dot_product(&residual, &corrected_residual);

    if !residual_dot_corrected.is_finite() || residual_dot_corrected <= 0.0 {
        return Err(FemError::ConjugateGradientBreakdown);
    }

    let mut previous_residual_norm = residual_norm;
    let mut stagnation_iterations = 0;

    for iteration in 0..options.max_iterations {
        operator.apply(&direction, &mut operator_direction)?;

        let denominator = dot_product(&direction, &operator_direction);

        // For a symmetric positive-definite matrix, this value must be
        // strictly positive.
        if !denominator.is_finite() || denominator <= 0.0 {
            return Err(FemError::ConjugateGradientBreakdown);
        }

        let alpha = residual_dot_corrected / denominator;

        for index in 0..dimension {
            solution[index] += alpha * direction[index];
            residual[index] -= alpha * operator_direction[index];
        }

        residual_norm = euclidean_norm(&residual);
        relative_residual = relative_residual_norm(residual_norm, rhs_norm);

        if relative_residual <= options.tolerance {
            return Ok(make_result(
                solution,
                iteration + 1,
                residual_norm,
                relative_residual,
                CgTerminationReason::Converged,
            ));
        }

        if residual_improvement(previous_residual_norm, residual_norm) <= options.stagnation_tolerance {
            stagnation_iterations += 1;
        } else {
            stagnation_iterations = 0;
        }

        if options.stagnation_window > 0 && stagnation_iterations >= options.stagnation_window {
            return Ok(make_result(
                solution,
                iteration + 1,
                residual_norm,
                relative_residual,
                CgTerminationReason::Stagnated,
            ));
        }

        previous_residual_norm = residual_norm;

        preconditioner.apply(&residual, &mut corrected_residual)?;

        let new_residual_dot_corrected = dot_product(&residual, &corrected_residual);

        if !new_residual_dot_corrected.is_finite() || new_residual_dot_corrected <= 0.0 {
            return Err(FemError::ConjugateGradientBreakdown);
        }

        let beta = new_residual_dot_corrected / residual_dot_corrected;

        for index in 0..dimension {
            direction[index] = corrected_residual[index] + beta * direction[index];
        }

        residual_dot_corrected = new_residual_dot_corrected;
    }

    Ok(make_result(
        solution,
        options.max_iterations,
        residual_norm,
        relative_residual,
        CgTerminationReason::MaxIterations,
    ))
}

/// Creates a solver result while keeping the convergence flag consistent
/// with the termination reason.
fn make_result(
    solution: Vec<f64>, iterations: usize, residual_norm: f64, relative_residual_norm: f64,
    termination_reason: CgTerminationReason,
) -> CgResult {
    CgResult {
        solution,
        iterations,
        residual_norm,
        relative_residual_norm,
        converged: termination_reason == CgTerminationReason::Converged,
        termination_reason,
    }
}

/// Computes `||r|| / max(||b||, epsilon)`.
fn relative_residual_norm(residual_norm: f64, rhs_norm: f64) -> f64 {
    residual_norm / rhs_norm.max(NORM_EPSILON)
}

/// Returns the relative improvement between two residual norms.
fn residual_improvement(previous: f64, current: f64) -> f64 {
    (previous - current) / previous.max(NORM_EPSILON)
}

/// The identity preconditioner leaves the residual unchanged.
struct IdentityPreconditioner;

impl Preconditioner for IdentityPreconditioner {
    fn apply(&self, residual: &[f64], corrected: &mut [f64]) -> Result<(), FemError> {
        if residual.len() != corrected.len() {
            return Err(FemError::InvalidVectorLength {
                vector: "corrected residual",
                expected: residual.len(),
                actual: corrected.len(),
            });
        }

        corrected.copy_from_slice(residual);
        Ok(())
    }
}

impl LinearOperator for CsrMatrix {
    fn dimension(&self) -> usize {
        self.nrows()
    }

    fn apply(&self, x: &[f64], y: &mut [f64]) -> Result<(), FemError> {
        self.mul_vector(x, y)
    }
}

/// Computes the scalar product of two equally-sized vectors.
fn dot_product(left: &[f64], right: &[f64]) -> f64 {
    left.iter().zip(right).map(|(left, right)| left * right).sum()
}

/// Computes the Euclidean norm of a vector.
fn euclidean_norm(vector: &[f64]) -> f64 {
    dot_product(vector, vector).sqrt()
}

#[cfg(test)]
mod tests {
    use super::{
        CgOptions, CgTerminationReason, JacobiPreconditioner, conjugate_gradient, preconditioned_conjugate_gradient,
    };
    use crate::analysis::sparse::CooMatrix;
    use crate::error::FemError;

    #[test]
    fn solves_small_symmetric_positive_definite_system() {
        let mut matrix = CooMatrix::new(2, 2);

        matrix.push(0, 0, 4.0).unwrap();
        matrix.push(0, 1, 1.0).unwrap();
        matrix.push(1, 0, 1.0).unwrap();
        matrix.push(1, 1, 3.0).unwrap();

        let matrix = matrix.into_csr();

        let result = conjugate_gradient(
            &matrix,
            &[1.0, 2.0],
            CgOptions { max_iterations: 100, tolerance: 1e-12, ..CgOptions::default() },
        )
        .unwrap();

        assert!(result.converged);
        assert_eq!(result.termination_reason, CgTerminationReason::Converged);
        assert!(result.relative_residual_norm < 1e-12);

        // Exact solution:
        //
        // [4 1] [x] = [1]
        // [1 3] [y]   [2]
        //
        // x = 1/11, y = 7/11
        assert!((result.solution[0] - 1.0 / 11.0).abs() < 1e-10);
        assert!((result.solution[1] - 7.0 / 11.0).abs() < 1e-10);
    }

    #[test]
    fn jacobi_preconditioner_scales_residual_by_the_diagonal() {
        let mut matrix = CooMatrix::new(2, 2);
        matrix.push(0, 0, 4.0).unwrap();
        matrix.push(1, 1, 2.0).unwrap();

        let matrix = matrix.into_csr();
        let preconditioner = JacobiPreconditioner::from_matrix(&matrix).unwrap();
        let mut corrected = vec![0.0; 2];

        super::Preconditioner::apply(&preconditioner, &[8.0, 6.0], &mut corrected).unwrap();

        assert_eq!(corrected, vec![2.0, 3.0]);
    }

    #[test]
    fn preconditioned_cg_converges_faster_for_ill_scaled_diagonal_system() {
        let mut matrix = CooMatrix::new(2, 2);
        matrix.push(0, 0, 1_000.0).unwrap();
        matrix.push(1, 1, 1.0).unwrap();

        let matrix = matrix.into_csr();
        let options = CgOptions { max_iterations: 10, tolerance: 1e-12, ..CgOptions::default() };
        let rhs = [1_000.0, 1.0];

        let plain_result = conjugate_gradient(&matrix, &rhs, options).unwrap();
        let preconditioner = JacobiPreconditioner::from_matrix(&matrix).unwrap();
        let preconditioned_result = preconditioned_conjugate_gradient(&matrix, &rhs, options, &preconditioner).unwrap();

        assert!(plain_result.converged);
        assert!(preconditioned_result.converged);
        assert!(preconditioned_result.iterations < plain_result.iterations);
    }

    #[test]
    fn reports_immediate_convergence_for_zero_right_hand_side() {
        let mut matrix = CooMatrix::new(2, 2);
        matrix.push(0, 0, 2.0).unwrap();
        matrix.push(1, 1, 3.0).unwrap();

        let matrix = matrix.into_csr();
        let result = conjugate_gradient(&matrix, &[0.0, 0.0], CgOptions::default()).unwrap();

        assert!(result.converged);
        assert_eq!(result.iterations, 0);
        assert_eq!(result.residual_norm, 0.0);
        assert_eq!(result.relative_residual_norm, 0.0);
        assert_eq!(result.termination_reason, CgTerminationReason::Converged);
    }

    #[test]
    fn reports_iteration_limit_without_claiming_convergence() {
        let mut matrix = CooMatrix::new(2, 2);
        matrix.push(0, 0, 4.0).unwrap();
        matrix.push(0, 1, 1.0).unwrap();
        matrix.push(1, 0, 1.0).unwrap();
        matrix.push(1, 1, 3.0).unwrap();

        let matrix = matrix.into_csr();
        let result =
            conjugate_gradient(&matrix, &[1.0, 2.0], CgOptions { max_iterations: 0, ..CgOptions::default() }).unwrap();

        assert!(!result.converged);
        assert_eq!(result.iterations, 0);
        assert_eq!(result.termination_reason, CgTerminationReason::MaxIterations);
        assert_eq!(result.relative_residual_norm, 1.0);
    }

    #[test]
    fn reports_stagnation_when_residual_improvement_is_too_small() {
        let mut matrix = CooMatrix::new(2, 2);
        matrix.push(0, 0, 4.0).unwrap();
        matrix.push(0, 1, 1.0).unwrap();
        matrix.push(1, 0, 1.0).unwrap();
        matrix.push(1, 1, 3.0).unwrap();

        let matrix = matrix.into_csr();
        let result = conjugate_gradient(
            &matrix,
            &[1.0, 2.0],
            CgOptions { tolerance: 1e-14, stagnation_window: 1, stagnation_tolerance: 1.0, ..CgOptions::default() },
        )
        .unwrap();

        assert!(!result.converged);
        assert_eq!(result.termination_reason, CgTerminationReason::Stagnated);
        assert_eq!(result.iterations, 1);
    }

    #[test]
    fn rejects_missing_jacobi_diagonal() {
        let matrix = CooMatrix::new(2, 2).into_csr();

        let result = JacobiPreconditioner::from_matrix(&matrix);

        assert!(matches!(result, Err(FemError::InvalidPreconditionerDiagonal { index: 0, value: 0.0 })));
    }
}
