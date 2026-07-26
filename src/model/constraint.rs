//! Defines displacement constraints for nodes.

use crate::FemError;
use crate::model::Dof2D;

/// 2D displacement constraint for a node, specifying the degree of freedom (DOF) and the prescribed displacement value.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DisplacementConstraint2D {
    node_id: usize,
    dof: Dof2D,
    displacement: f64,
}

impl DisplacementConstraint2D {
    /// Creates a new 2D displacement constraint for a node with the given ID, degree of freedom (DOF), and prescribed displacement value.
    pub fn new(node_id: usize, dof: Dof2D, displacement: f64) -> Result<Self, FemError> {
        if !displacement.is_finite() {
            return Err(FemError::InvalidDisplacementConstraint { node_id, displacement });
        }

        Ok(Self { node_id, dof, displacement })
    }
    /// Returns the ID of the node associated with this displacement constraint.
    #[must_use]
    pub fn node_id(&self) -> usize {
        self.node_id
    }

    /// Returns the degree of freedom (DOF) associated with this displacement constraint.
    #[must_use]
    pub fn dof(&self) -> Dof2D {
        self.dof
    }

    /// Returns the prescribed displacement value for this constraint.
    #[must_use]
    pub fn displacement(&self) -> f64 {
        self.displacement
    }
}

#[cfg(test)]
mod test {
    use super::DisplacementConstraint2D;
    use crate::model::Dof2D;

    #[test]
    fn creates_correct_constraints() {
        let cases = [
            (1, Dof2D::Ux, 0.0),
            (2, Dof2D::Ux, 1.0),
            (45, Dof2D::Ux, -1123.2),
            (12, Dof2D::Uy, 0.0),
            (245, Dof2D::Uy, 2.32),
            (1, Dof2D::Uy, -42.22),
            (556, Dof2D::Rz, 0.01),
            (5321, Dof2D::Rz, -0.2),
            (32, Dof2D::Rz, 0.0),
        ];

        for (node_id, dof, displacement) in cases {
            let constraint =
                DisplacementConstraint2D::new(node_id, dof, displacement).expect("valid constraint should be created");

            assert_eq!(constraint.node_id(), node_id);
            assert_eq!(constraint.dof(), dof);
            assert_eq!(constraint.displacement(), displacement);
        }
    }

    #[test]
    fn rejects_non_finite_displacement() {
        let result = DisplacementConstraint2D::new(1, Dof2D::Ux, f64::NAN);

        assert!(result.is_err());
    }
}
