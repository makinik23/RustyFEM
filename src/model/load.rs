//! Describes load conditions applied to nodes.

use crate::FemError;
use crate::model::Dof2D;

/// Struct representing a nodal load in 2D space.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct NodalLoad2D {
    node_id: usize,
    dof: Dof2D,
    value: f64,
}

impl NodalLoad2D {
    pub fn new(node_id: usize, dof: Dof2D, value: f64) -> Result<Self, FemError> {
        if !value.is_finite() {
            return Err(FemError::InvalidNodalLoad { node_id, value });
        }

        Ok(Self { node_id, dof, value })
    }

    pub fn node_id(&self) -> usize {
        self.node_id
    }

    pub fn dof(&self) -> Dof2D {
        self.dof
    }

    pub fn value(&self) -> f64 {
        self.value
    }
}

#[cfg(test)]
mod tests {
    use super::NodalLoad2D;
    use crate::error::FemError;
    use crate::model::Dof2D;

    #[test]
    fn creates_nodal_load_with_valid_data() {
        let load = NodalLoad2D::new(7, Dof2D::Ux, -12.5).expect("valid load should be created");

        assert_eq!(load.node_id(), 7);
        assert_eq!(load.dof(), Dof2D::Ux);
        assert_eq!(load.value(), -12.5);
    }

    #[test]
    fn rejects_non_finite_nodal_load_values() {
        let cases = [("NaN", f64::NAN), ("positive infinity", f64::INFINITY), ("negative infinity", f64::NEG_INFINITY)];

        for (name, value) in cases {
            let result = NodalLoad2D::new(7, Dof2D::Uy, value);

            assert!(
                matches!(
                    result,
                    Err(FemError::InvalidNodalLoad { node_id: 7, value: actual_value })
                        if actual_value == value || (actual_value.is_nan() && value.is_nan())
                ),
                "failed case: {name}"
            );
        }
    }
}
