use serde::Serialize;

use crate::FemError;

/// Concentrated nodal force for a one-dimensional bar model.
#[derive(Debug, Clone, Copy, PartialEq, Serialize)]
pub struct NodalLoad1D {
    node_id: usize,
    force: f64,
}

impl NodalLoad1D {
    /// Creates a nodal force in newtons.
    pub fn new(node_id: usize, force: f64) -> Result<Self, FemError> {
        if !force.is_finite() {
            return Err(FemError::InvalidNodalLoad { node_id, force });
        }

        Ok(Self { node_id, force })
    }

    /// External node identifier.
    #[must_use]
    pub fn node_id(&self) -> usize {
        self.node_id
    }

    /// Force in newtons.
    #[must_use]
    pub fn force(&self) -> f64 {
        self.force
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn creates_nodal_load() {
        let load = NodalLoad1D::new(1, -12.0).unwrap();

        assert_eq!(load.node_id(), 1);
        assert_eq!(load.force(), -12.0);
    }

    #[test]
    fn rejects_non_finite_force() {
        let error = NodalLoad1D::new(1, f64::INFINITY).unwrap_err();

        assert_eq!(error, FemError::InvalidNodalLoad { node_id: 1, force: f64::INFINITY });
    }
}
