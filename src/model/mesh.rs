//! Mesh data for a two-dimensional finite element model.

use super::node::Node2D;
use crate::elements::Element2D;

/// Stores the geometric and topological part of a 2D model.
#[derive(Default)]
pub struct Mesh2D {
    nodes: Vec<Node2D>,
    elements: Vec<Element2D>,
}

impl Mesh2D {
    /// Creates an empty 2D mesh.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Adds a node without validating model-level uniqueness.
    pub(crate) fn push_node(&mut self, node: Node2D) {
        self.nodes.push(node);
    }

    /// Adds an element without validating model-level connectivity.
    pub(crate) fn push_element(&mut self, element: Element2D) {
        self.elements.push(element);
    }

    /// Returns all nodes in insertion order.
    #[must_use]
    pub fn nodes(&self) -> &[Node2D] {
        &self.nodes
    }

    /// Returns all elements in insertion order.
    #[must_use]
    pub fn elements(&self) -> &[Element2D] {
        &self.elements
    }

    /// Finds a node by ID.
    #[must_use]
    pub(crate) fn node(&self, node_id: usize) -> Option<&Node2D> {
        self.nodes.iter().find(|node| node.id() == node_id)
    }

    /// Finds an element by ID.
    #[must_use]
    pub(crate) fn element(&self, element_id: usize) -> Option<&Element2D> {
        self.elements.iter().find(|element| element.id() == element_id)
    }

    /// Checks whether an element ID already exists.
    #[must_use]
    pub(crate) fn contains_element_id(&self, element_id: usize) -> bool {
        self.elements.iter().any(|element| element.id() == element_id)
    }

    /// Checks whether a node ID already exists.
    #[must_use]
    pub(crate) fn contains_node_id(&self, node_id: usize) -> bool {
        self.nodes.iter().any(|node| node.id() == node_id)
    }
}

#[cfg(test)]
mod tests {
    use super::Mesh2D;
    use crate::elements::{Element2D, Truss2D};
    use crate::model::Node2D;

    #[test]
    fn creates_empty_mesh() {
        let mesh = Mesh2D::new();

        assert!(mesh.nodes().is_empty());
        assert!(mesh.elements().is_empty());
    }

    #[test]
    fn stores_nodes_in_insertion_order() {
        let mut mesh = Mesh2D::new();

        mesh.push_node(Node2D::new(10, 1.0, 2.0).expect("valid node"));
        mesh.push_node(Node2D::new(20, 3.0, 4.0).expect("valid node"));

        assert_eq!(mesh.nodes().len(), 2);
        assert_eq!(mesh.nodes()[0].id(), 10);
        assert_eq!(mesh.nodes()[0].coordinates(), (1.0, 2.0));
        assert_eq!(mesh.nodes()[1].id(), 20);
        assert_eq!(mesh.nodes()[1].coordinates(), (3.0, 4.0));
    }

    #[test]
    fn finds_node_by_id() {
        let mut mesh = Mesh2D::new();

        mesh.push_node(Node2D::new(10, 1.0, 2.0).expect("valid node"));

        let node = mesh.node(10).expect("node should exist");

        assert_eq!(node.coordinates(), (1.0, 2.0));
        assert!(mesh.node(99).is_none());
    }

    #[test]
    fn reports_existing_node_ids() {
        let mut mesh = Mesh2D::new();

        mesh.push_node(Node2D::new(10, 1.0, 2.0).expect("valid node"));

        assert!(mesh.contains_node_id(10));
        assert!(!mesh.contains_node_id(99));
    }

    #[test]
    fn stores_elements_in_insertion_order() {
        let mut mesh = Mesh2D::new();

        let first = Element2D::Truss(Truss2D::new(100, [1, 2], 1, 10).expect("valid truss"));
        let second = Element2D::Truss(Truss2D::new(200, [2, 3], 1, 20).expect("valid truss"));

        mesh.push_element(first);
        mesh.push_element(second);

        assert_eq!(mesh.elements(), &[first, second]);
    }

    #[test]
    fn finds_element_by_id() {
        let mut mesh = Mesh2D::new();

        let element = Element2D::Truss(Truss2D::new(100, [1, 2], 1, 10).expect("valid truss"));

        mesh.push_element(element);

        assert_eq!(mesh.element(100), Some(&element));
        assert!(mesh.element(999).is_none());
    }

    #[test]
    fn reports_existing_element_ids() {
        let mut mesh = Mesh2D::new();

        let element = Element2D::Truss(Truss2D::new(100, [1, 2], 1, 10).expect("valid truss"));

        mesh.push_element(element);

        assert!(mesh.contains_element_id(100));
        assert!(!mesh.contains_element_id(999));
    }
}
