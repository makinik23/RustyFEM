//! Defines a 2D model consisting of nodes and displacement constraints.

use super::constraint::DisplacementConstraint2D;
use super::material::Material2D;
use super::node::Node2D;
use crate::elements::Element2D;
use crate::error::FemError;

/// Represents a 2D model consisting of nodes and displacement constraints.
#[derive(Default)]
pub struct Model2D {
    nodes: Vec<Node2D>,
    constraints: Vec<DisplacementConstraint2D>,
    elements: Vec<Element2D>,
    material: Option<Material2D>,
}

impl Model2D {
    /// Creates a new empty 2D model.
    pub fn new() -> Self {
        Self::default()
    }

    /// Adds a node to the model. Returns an error if a node with the same ID already exists.
    pub fn add_node(&mut self, node: Node2D) -> Result<(), FemError> {
        if self.nodes.iter().any(|existing| existing.id() == node.id()) {
            return Err(FemError::DuplicateId { entity: "node", id: node.id() });
        }

        self.nodes.push(node);

        Ok(())
    }

    /// Adds a displacement constraint to the model. Returns an error if the node ID associated with the constraint does not exist in the model.
    pub fn add_constraint(&mut self, constraint: DisplacementConstraint2D) -> Result<(), FemError> {
        if !self.nodes.iter().any(|node| node.id() == constraint.node_id()) {
            return Err(FemError::UnknownId { entity: "node", id: constraint.node_id() });
        }

        self.constraints.push(constraint);

        Ok(())
    }

    /// Adds an element to the model. Returns an error if an element with the same ID already exists
    /// or if any of the node IDs associated with the element do not exist in the model.
    pub fn add_element(&mut self, element: Element2D) -> Result<(), FemError> {
        if self.elements.iter().any(|existing| existing.id() == element.id()) {
            return Err(FemError::DuplicateId { entity: "element", id: element.id() });
        }

        for &node_id in element.node_ids() {
            self.find_node(node_id)?;
        }

        self.validate_element_geometry(&element)?;

        self.elements.push(element);

        Ok(())
    }

    /// Sets the material properties for the model.
    pub fn set_material(&mut self, material: Material2D) {
        self.material = Some(material);
    }

    /// Finds a node in the model by its ID. Returns a reference to the node if found, or an error if the node ID does not exist in the model.
    fn find_node(&self, node_id: usize) -> Result<&Node2D, FemError> {
        self.nodes.iter().find(|node| node.id() == node_id).ok_or(FemError::UnknownId { entity: "node", id: node_id })
    }

    /// Validates the geometry of an element. Returns an error if the element is degenerate (e.g., zero length for trusses and beams, or zero area for triangles).
    fn validate_element_geometry(&self, element: &Element2D) -> Result<(), FemError> {
        let node_ids = element.node_ids();
        let first_node = self.find_node(node_ids[0])?;
        let second_node = self.find_node(node_ids[1])?;
        let element_type = match element {
            Element2D::Truss(_) => "truss",
            Element2D::Beam(_) => "beam",
            Element2D::TriangleT3(_) => "triangle_t3",
        };

        match element {
            Element2D::Truss(_) | Element2D::Beam(_) => {
                let dx = second_node.x() - first_node.x();
                let dy = second_node.y() - first_node.y();
                let length = (dx * dx + dy * dy).sqrt();

                if length == 0.0 {
                    return Err(FemError::DegenerateElement {
                        element_id: element.id(),
                        element_type,
                        node_ids: node_ids.to_vec(),
                        measure_name: "length",
                        measure: length,
                    });
                }
            }
            Element2D::TriangleT3(_) => {
                let third_node = self.find_node(node_ids[2])?;
                let first_to_second_x = second_node.x() - first_node.x();
                let first_to_second_y = second_node.y() - first_node.y();
                let first_to_third_x = third_node.x() - first_node.x();
                let first_to_third_y = third_node.y() - first_node.y();
                let area = 0.5 * (first_to_second_x * first_to_third_y - first_to_second_y * first_to_third_x).abs();

                if area == 0.0 {
                    return Err(FemError::DegenerateElement {
                        element_id: element.id(),
                        element_type,
                        node_ids: node_ids.to_vec(),
                        measure_name: "area",
                        measure: area,
                    });
                }
            }
        }

        Ok(())
    }

    /// Returns a slice of all nodes in the model.
    #[must_use]
    pub fn nodes(&self) -> &[Node2D] {
        &self.nodes
    }

    /// Returns a slice of all displacement constraints in the model.
    #[must_use]
    pub fn constraints(&self) -> &[DisplacementConstraint2D] {
        &self.constraints
    }

    /// Returns a slice of all elements in the model.
    #[must_use]
    pub fn elements(&self) -> &[Element2D] {
        &self.elements
    }

    /// Returns the material properties of the model.
    #[must_use]
    pub fn material(&self) -> Option<&Material2D> {
        self.material.as_ref()
    }
}

#[cfg(test)]
mod tests {
    use super::Model2D;
    use crate::FemError;
    use crate::elements::{Beam2D, Element2D, TriangleT3, Truss2D};
    use crate::model::{DisplacementConstraint2D, Dof2D, Material2D, Node2D};

    #[test]
    fn creates_empty_model() {
        let model = Model2D::new();

        assert!(model.nodes().is_empty());
        assert!(model.constraints().is_empty());
        assert!(model.elements().is_empty());
        assert!(model.material().is_none());
    }

    #[test]
    fn adds_node_to_model() {
        let mut model = Model2D::new();
        let node = Node2D::new(10, 1.5, -2.0).expect("valid node should be created");

        model.add_node(node).expect("node should be added");

        assert_eq!(model.nodes().len(), 1);
        assert_eq!(model.nodes()[0].id(), 10);
        assert_eq!(model.nodes()[0].coordinates(), (1.5, -2.0));
    }

    #[test]
    fn rejects_duplicate_node_id() {
        let mut model = Model2D::new();

        let first_node = Node2D::new(10, 1.0, 2.0).expect("valid node should be created");
        let second_node = Node2D::new(10, 3.0, 4.0).expect("valid node should be created");

        model.add_node(first_node).expect("first node should be added");

        let result = model.add_node(second_node);

        assert!(matches!(result, Err(FemError::DuplicateId { entity: "node", id: 10 })));

        assert_eq!(model.nodes().len(), 1);
    }

    #[test]
    fn adds_constraint_for_existing_node() {
        let mut model = Model2D::new();

        let node = Node2D::new(10, 1.0, 2.0).expect("valid node should be created");
        model.add_node(node).expect("node should be added");

        let constraint = DisplacementConstraint2D::new(10, Dof2D::Uy, 0.0).expect("valid constraint should be created");

        model.add_constraint(constraint).expect("constraint should be added");

        assert_eq!(model.constraints().len(), 1);
        assert_eq!(model.constraints()[0].node_id(), 10);
        assert_eq!(model.constraints()[0].dof(), Dof2D::Uy);
        assert_eq!(model.constraints()[0].displacement(), 0.0);
    }

    #[test]
    fn rejects_constraint_for_unknown_node() {
        let mut model = Model2D::new();

        let constraint = DisplacementConstraint2D::new(99, Dof2D::Ux, 0.0).expect("constraint itself should be valid");

        let result = model.add_constraint(constraint);

        assert!(matches!(result, Err(FemError::UnknownId { entity: "node", id: 99 })));

        assert!(model.constraints().is_empty());
    }

    #[test]
    fn adds_element_for_existing_nodes() {
        let mut model = Model2D::new();

        let first_node = Node2D::new(1, 0.0, 0.0).expect("valid node should be created");

        let second_node = Node2D::new(2, 1.0, 0.0).expect("valid node should be created");

        model.add_node(first_node).expect("node should be added");
        model.add_node(second_node).expect("node should be added");

        let truss = Truss2D::new(10, [1, 2]).expect("valid truss should be created");
        let element = Element2D::Truss(truss);

        model.add_element(element).expect("element should be added");

        assert_eq!(model.elements().len(), 1);
        assert_eq!(model.elements()[0].id(), 10);
        assert_eq!(model.elements()[0].node_ids(), &[1, 2]);
    }

    #[test]
    fn rejects_element_for_unknown_node() {
        let mut model = Model2D::new();

        let node = Node2D::new(1, 0.0, 0.0).expect("valid node should be created");

        model.add_node(node).expect("node should be added");

        let element = Element2D::Truss(Truss2D::new(10, [1, 99]).expect("valid truss should be created"));

        let result = model.add_element(element);

        assert!(matches!(result, Err(FemError::UnknownId { entity: "node", id: 99 })));

        assert!(model.elements().is_empty());
    }

    #[test]
    fn rejects_degenerate_elements() {
        let mut model = Model2D::new();

        for node_id in 1..=3 {
            let node = Node2D::new(node_id, 1.0, 1.0).expect("valid node should be created");

            model.add_node(node).expect("node should be added");
        }

        let cases = [
            (
                "truss",
                Element2D::Truss(Truss2D::new(10, [1, 2]).expect("valid truss connectivity should be created")),
                "length",
            ),
            (
                "beam",
                Element2D::Beam(Beam2D::new(20, [1, 2]).expect("valid beam connectivity should be created")),
                "length",
            ),
            (
                "triangle_t3",
                Element2D::TriangleT3(
                    TriangleT3::new(30, [1, 2, 3]).expect("valid triangle connectivity should be created"),
                ),
                "area",
            ),
        ];

        for (element_type, element, measure_name) in cases {
            let element_id = element.id();
            let result = model.add_element(element);

            assert!(matches!(
                result,
                Err(FemError::DegenerateElement {
                    element_id: actual_id,
                    element_type: actual_type,
                    measure_name: actual_measure_name,
                    measure,
                    ..
                }) if actual_id == element_id
                    && actual_type == element_type
                    && actual_measure_name == measure_name
                    && measure == 0.0
            ));
        }

        assert!(model.elements().is_empty());
    }

    #[test]
    fn replaces_existing_material() {
        let mut model = Model2D::new();

        let first = Material2D::new(210e9, 0.3, 7800.0).expect("valid material should be created");

        let second = Material2D::new(70e9, 0.33, 2700.0).expect("valid material should be created");

        model.set_material(first);
        model.set_material(second);

        assert_eq!(model.material(), Some(&second));
    }
}
