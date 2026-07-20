use serde::Serialize;

use crate::FemError;

/// One-dimensional node expressed in SI units.
#[derive(Debug, Clone, Copy, PartialEq, Serialize)]
pub struct Node1D {
    id: usize,
    x: f64,
}

impl Node1D {
    /// Creates a node with a finite coordinate in meters.
    pub fn new(id: usize, x: f64) -> Result<Self, FemError> {
        if !x.is_finite() {
            return Err(FemError::InvalidNodeCoordinate { node_id: id, x });
        }

        Ok(Self { id, x })
    }

    /// External node identifier.
    #[must_use]
    pub fn id(&self) -> usize {
        self.id
    }

    /// Coordinate in meters.
    #[must_use]
    pub fn x(&self) -> f64 {
        self.x
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn creates_node_with_finite_coordinate() {
        let node = Node1D::new(7, 1.25).unwrap();

        assert_eq!(node.id(), 7);
        assert_eq!(node.x(), 1.25);
    }

    #[test]
    fn rejects_non_finite_coordinate() {
        let error = Node1D::new(7, f64::NAN).unwrap_err();

        assert!(matches!(
            error,
            FemError::InvalidNodeCoordinate {
                node_id: 7,
                x
            } if x.is_nan()
        ));
    }
}
