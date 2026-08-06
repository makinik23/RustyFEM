//! Load data for a two-dimensional finite element model.

use super::element_load::ElementLoad2D;
use super::nodal_load::NodalLoad2D;

/// Stores loads applied to a 2D model.
#[derive(Default)]
pub struct Loads2D {
    nodal_loads: Vec<NodalLoad2D>,
    element_loads: Vec<ElementLoad2D>,
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

    /// Adds an element load to the load collection.
    pub(crate) fn push_element_load(&mut self, load: ElementLoad2D) {
        self.element_loads.push(load);
    }

    /// Returns all nodal loads in insertion order.
    #[must_use]
    pub fn nodal_loads(&self) -> &[NodalLoad2D] {
        &self.nodal_loads
    }

    /// Returns all element loads in insertion order.
    #[must_use]
    pub fn element_loads(&self) -> &[ElementLoad2D] {
        &self.element_loads
    }
}

#[cfg(test)]
mod tests {
    use super::Loads2D;
    use crate::model::{BeamUniformLineLoad2D, Dof2D, ElementLoad2D, LoadCoordinateSystem2D, NodalLoad2D};

    #[test]
    fn creates_empty_loads() {
        let loads = Loads2D::new();

        assert!(loads.nodal_loads().is_empty());
        assert!(loads.element_loads().is_empty());
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

    #[test]
    fn stores_element_loads_in_insertion_order() {
        let mut loads = Loads2D::new();

        let first = ElementLoad2D::BeamUniformLine(
            BeamUniformLineLoad2D::new(10, LoadCoordinateSystem2D::Local, 0.0, -1.0).expect("valid load"),
        );
        let second = ElementLoad2D::BeamUniformLine(
            BeamUniformLineLoad2D::new(20, LoadCoordinateSystem2D::Global, 2.0, 0.0).expect("valid load"),
        );

        loads.push_element_load(first);
        loads.push_element_load(second);

        assert_eq!(loads.element_loads(), &[first, second]);
    }
}
