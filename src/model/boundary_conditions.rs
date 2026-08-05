//! Boundary condition data for a two-dimensional finite element model.

use super::constraint::DisplacementConstraint2D;

///  Stores the boundary conditions of a 2D model, including displacement constraints.
#[derive(Default)]
pub struct BoundaryConditions2D {
    displacement_constraints: Vec<DisplacementConstraint2D>,
}

impl BoundaryConditions2D {
    /// Creates an empty set of boundary conditions.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Adds a displacement constraint to the boundary conditions.
    pub(crate) fn push_displacement_constraint(&mut self, constraint: DisplacementConstraint2D) {
        self.displacement_constraints.push(constraint);
    }

    /// Returns all displacement constraints in insertion order.
    #[must_use]
    pub fn displacement_constraints(&self) -> &[DisplacementConstraint2D] {
        &self.displacement_constraints
    }
}

#[cfg(test)]
mod tests {
    use super::BoundaryConditions2D;
    use crate::model::{DisplacementConstraint2D, Dof2D};

    #[test]
    fn creates_empty_boundary_conditions() {
        let boundary_conditions = BoundaryConditions2D::new();

        assert!(boundary_conditions.displacement_constraints().is_empty());
    }

    #[test]
    fn stores_displacement_constraints_in_insertion_order() {
        let mut boundary_conditions = BoundaryConditions2D::new();

        let first = DisplacementConstraint2D::new(10, Dof2D::Ux, 0.0).expect("valid constraint");
        let second = DisplacementConstraint2D::new(20, Dof2D::Uy, 1.5).expect("valid constraint");

        boundary_conditions.push_displacement_constraint(first);
        boundary_conditions.push_displacement_constraint(second);

        assert_eq!(boundary_conditions.displacement_constraints(), &[first, second]);
    }
}
