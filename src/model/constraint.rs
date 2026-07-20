use serde::Serialize;

use crate::FemError;

/// Prescribed displacement for a one-dimensional bar model.
#[derive(Debug, Clone, Copy, PartialEq, Serialize)]
pub struct DisplacementConstraint1D {
    node_id: usize,
    displacement: f64,
}

impl DisplacementConstraint1D {
    /// Creates a prescribed displacement in meters.
    pub fn new(node_id: usize, displacement: f64) -> Result<Self, FemError> {
        if !displacement.is_finite() {
            return Err(FemError::InvalidDisplacementConstraint { node_id, displacement });
        }

        Ok(Self { node_id, displacement })
    }

    /// External node identifier.
    #[must_use]
    pub fn node_id(&self) -> usize {
        self.node_id
    }

    /// Prescribed displacement in meters.
    #[must_use]
    pub fn displacement(&self) -> f64 {
        self.displacement
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn creates_displacement_constraint() {
        let constraint = DisplacementConstraint1D::new(1, 1.0e-3).unwrap();

        assert_eq!(constraint.node_id(), 1);
        assert_eq!(constraint.displacement(), 1.0e-3);
    }

    #[test]
    fn rejects_non_finite_displacement() {
        let error = DisplacementConstraint1D::new(1, f64::NEG_INFINITY).unwrap_err();

        assert_eq!(error, FemError::InvalidDisplacementConstraint { node_id: 1, displacement: f64::NEG_INFINITY });
    }
}
