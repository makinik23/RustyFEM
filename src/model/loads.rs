//! Load data for a two-dimensional finite element model.

use super::nodal_load::NodalLoad2D;

/// Stores loads applied to a 2D model.
#[derive(Default)]
pub struct Loads2D {
    nodal_loads: Vec<NodalLoad2D>,
}

impl Loads2D {
    /// Creates an empty set of loads.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Adds a nodal load to the load collection.
    pub(crate) fn push_nodal_load(&mut self, load: NodalLoad2D) {
        self.nodal_loads.push(load);
    }

    /// Returns all nodal loads in insertion order.
    #[must_use]
    pub fn nodal_loads(&self) -> &[NodalLoad2D] {
        &self.nodal_loads
    }
}

#[cfg(test)]
mod tests {
    use super::Loads2D;
    use crate::model::{Dof2D, NodalLoad2D};

    #[test]
    fn creates_empty_loads() {
        let loads = Loads2D::new();

        assert!(loads.nodal_loads().is_empty());
    }

    #[test]
    fn stores_nodal_loads_in_insertion_order() {
        let mut loads = Loads2D::new();

        let first = NodalLoad2D::new(10, Dof2D::Ux, 100.0).expect("valid load");
        let second = NodalLoad2D::new(20, Dof2D::Uy, -50.0).expect("valid load");

        loads.push_nodal_load(first);
        loads.push_nodal_load(second);

        assert_eq!(loads.nodal_loads(), &[first, second]);
    }
}
