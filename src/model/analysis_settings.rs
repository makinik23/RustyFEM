//! Analysis settings for a two-dimensional finite element model.

use crate::FemError;

/// Selects the linear solver used for a 2D analysis.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SolverKind2D {
    /// Solves the system with dense LU decomposition.
    Dense,

    /// Solves the system with sparse CSR storage and Conjugate Gradient.
    Sparse,
}

/// Stores solver settings used by a 2D model.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AnalysisSettings2D {
    solver: SolverKind2D,
    cg_max_iterations: usize,
    cg_tolerance: f64,
    cg_stagnation_window: usize,
    cg_stagnation_tolerance: f64,
}

impl Default for AnalysisSettings2D {
    fn default() -> Self {
        Self {
            solver: SolverKind2D::Dense,
            cg_max_iterations: 1_000,
            cg_tolerance: 1e-10,
            cg_stagnation_window: 8,
            cg_stagnation_tolerance: 1e-12,
        }
    }
}

impl AnalysisSettings2D {
    /// Creates default analysis settings.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Returns the selected solver.
    #[must_use]
    pub fn solver(&self) -> SolverKind2D {
        self.solver
    }

    /// Sets the selected solver.
    pub fn set_solver(&mut self, solver: SolverKind2D) {
        self.solver = solver;
    }

    /// Returns the maximum number of Conjugate Gradient iterations.
    #[must_use]
    pub fn cg_max_iterations(&self) -> usize {
        self.cg_max_iterations
    }

    /// Sets the maximum number of Conjugate Gradient iterations.
    pub fn set_cg_max_iterations(&mut self, value: usize) -> Result<(), FemError> {
        if value == 0 {
            return Err(FemError::InvalidSolverIterationLimit { value });
        }

        self.cg_max_iterations = value;

        Ok(())
    }

    /// Returns the relative Conjugate Gradient residual tolerance.
    #[must_use]
    pub fn cg_tolerance(&self) -> f64 {
        self.cg_tolerance
    }

    /// Sets the relative Conjugate Gradient residual tolerance.
    pub fn set_cg_tolerance(&mut self, value: f64) -> Result<(), FemError> {
        if !value.is_finite() || value <= 0.0 {
            return Err(FemError::InvalidSolverTolerance { value });
        }

        self.cg_tolerance = value;

        Ok(())
    }

    /// Returns the number of consecutive low-improvement CG iterations that trigger stagnation.
    #[must_use]
    pub fn cg_stagnation_window(&self) -> usize {
        self.cg_stagnation_window
    }

    /// Sets the number of consecutive low-improvement CG iterations that trigger stagnation.
    pub fn set_cg_stagnation_window(&mut self, value: usize) {
        self.cg_stagnation_window = value;
    }

    /// Returns the minimum relative residual improvement required to reset stagnation counting.
    #[must_use]
    pub fn cg_stagnation_tolerance(&self) -> f64 {
        self.cg_stagnation_tolerance
    }

    /// Sets the minimum relative residual improvement required to reset stagnation counting.
    pub fn set_cg_stagnation_tolerance(&mut self, value: f64) -> Result<(), FemError> {
        if !value.is_finite() || value < 0.0 {
            return Err(FemError::InvalidStagnationTolerance { value });
        }

        self.cg_stagnation_tolerance = value;

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::{AnalysisSettings2D, SolverKind2D};
    use crate::FemError;

    #[test]
    fn creates_default_analysis_settings() {
        let settings = AnalysisSettings2D::new();

        assert_eq!(settings.solver(), SolverKind2D::Dense);
        assert_eq!(settings.cg_max_iterations(), 1_000);
        assert_eq!(settings.cg_tolerance(), 1e-10);
        assert_eq!(settings.cg_stagnation_window(), 8);
        assert_eq!(settings.cg_stagnation_tolerance(), 1e-12);
    }

    #[test]
    fn updates_analysis_settings() {
        let mut settings = AnalysisSettings2D::new();

        settings.set_solver(SolverKind2D::Sparse);
        settings.set_cg_max_iterations(250).expect("iteration limit should be valid");
        settings.set_cg_tolerance(1e-8).expect("tolerance should be valid");
        settings.set_cg_stagnation_window(4);
        settings.set_cg_stagnation_tolerance(1e-9).expect("stagnation tolerance should be valid");

        assert_eq!(settings.solver(), SolverKind2D::Sparse);
        assert_eq!(settings.cg_max_iterations(), 250);
        assert_eq!(settings.cg_tolerance(), 1e-8);
        assert_eq!(settings.cg_stagnation_window(), 4);
        assert_eq!(settings.cg_stagnation_tolerance(), 1e-9);
    }

    #[test]
    fn rejects_invalid_iteration_limit() {
        let mut settings = AnalysisSettings2D::new();
        let result = settings.set_cg_max_iterations(0);

        assert!(matches!(result, Err(FemError::InvalidSolverIterationLimit { value: 0 })));
    }

    #[test]
    fn rejects_invalid_tolerance() {
        let mut settings = AnalysisSettings2D::new();

        for value in [0.0, -1e-10, f64::INFINITY, f64::NAN] {
            let result = settings.set_cg_tolerance(value);

            assert!(
                matches!(result, Err(FemError::InvalidSolverTolerance { value: actual }) if actual == value || actual.is_nan() && value.is_nan())
            );
        }
    }

    #[test]
    fn rejects_invalid_stagnation_tolerance() {
        let mut settings = AnalysisSettings2D::new();

        for value in [-1e-12, f64::INFINITY, f64::NAN] {
            let result = settings.set_cg_stagnation_tolerance(value);

            assert!(
                matches!(result, Err(FemError::InvalidStagnationTolerance { value: actual }) if actual == value || actual.is_nan() && value.is_nan())
            );
        }
    }
}
