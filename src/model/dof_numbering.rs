//! Defines the numbering of degrees of freedom (DOF) for a 2D finite element model.

use std::collections::HashMap;

use super::{Dof2D, Model2D};
use crate::elements::Element2D;
use crate::error::FemError;

/// Struct representing the numbering of degrees of freedom (DOF) for a 2D finite element model.
pub struct DofNumbering2D {
    indices: HashMap<(usize, Dof2D), usize>,
}

impl DofNumbering2D {
    /// Creates a new `DofNumbering2D` instance from a given `Model2D`.
    pub fn from_model(model: &Model2D) -> Result<Self, FemError> {
        let mut indices = HashMap::new();
        let mut next_index = 0;

        let ordered_dofs = [Dof2D::Ux, Dof2D::Uy, Dof2D::Rz];

        for node in model.nodes() {
            for dof in ordered_dofs {
                let is_required = model
                    .elements()
                    .iter()
                    .any(|element| element.node_ids().contains(&node.id()) && element.dofs_per_node().contains(&dof));

                if is_required {
                    indices.insert((node.id(), dof), next_index);

                    next_index += 1;
                }
            }
        }

        for constraint in model.constraints() {
            let key = (constraint.node_id(), constraint.dof());

            if !indices.contains_key(&key) {
                return Err(FemError::UnknownDof { node_id: constraint.node_id(), dof: constraint.dof().name() });
            }
        }

        Ok(Self { indices })
    }

    /// Returns the index of the specified degree of freedom (DOF) for a given node ID.
    pub fn index(&self, node_id: usize, dof: Dof2D) -> Result<usize, FemError> {
        self.indices.get(&(node_id, dof)).copied().ok_or(FemError::UnknownDof { node_id, dof: dof.name() })
    }

    pub fn element_dof_indices(&self, element: &Element2D) -> Result<Vec<usize>, FemError> {
        let mut indices = Vec::with_capacity(element.dof_count());

        for &node_id in element.node_ids() {
            for &dof in element.dofs_per_node() {
                indices.push(self.index(node_id, dof)?);
            }
        }

        Ok(indices)
    }

    pub fn constraint_dof_indices(&self, model: &Model2D) -> Result<Vec<(usize, f64)>, FemError> {
        model
            .constraints()
            .iter()
            .map(|constraint| {
                let global_index = self.index(constraint.node_id(), constraint.dof())?;

                Ok((global_index, constraint.displacement()))
            })
            .collect()
    }

    #[must_use]
    /// Returns the total count of degrees of freedom (DOF) in the model.
    pub fn count(&self) -> usize {
        self.indices.len()
    }
}

#[cfg(test)]
mod tests {
    use super::DofNumbering2D;
    use crate::elements::{Beam2D, Element2D, TriangleT3, Truss2D};
    use crate::model::DisplacementConstraint2D;
    use crate::model::{Dof2D, Model2D, Node2D};

    #[test]
    fn numbers_beam_degrees_of_freedom() {
        let mut model = Model2D::new();

        model.add_node(Node2D::new(10, 0.0, 0.0).expect("valid node")).expect("node should be added");

        model.add_node(Node2D::new(20, 1.0, 0.0).expect("valid node")).expect("node should be added");

        let beam = Beam2D::new(1, [10, 20]).expect("valid beam");

        model.add_element(Element2D::Beam(beam)).expect("element should be added");

        let numbering = DofNumbering2D::from_model(&model).expect("numbering should be created");

        assert_eq!(numbering.count(), 6);

        assert_eq!(numbering.index(10, Dof2D::Ux).unwrap(), 0);
        assert_eq!(numbering.index(10, Dof2D::Uy).unwrap(), 1);
        assert_eq!(numbering.index(10, Dof2D::Rz).unwrap(), 2);

        assert_eq!(numbering.index(20, Dof2D::Ux).unwrap(), 3);
        assert_eq!(numbering.index(20, Dof2D::Uy).unwrap(), 4);
        assert_eq!(numbering.index(20, Dof2D::Rz).unwrap(), 5);
    }

    #[test]
    fn maps_local_dofs_for_every_element_type() {
        let cases = [
            (
                "truss",
                Element2D::Truss(Truss2D::new(1, [10, 20]).expect("valid truss should be created")),
                vec![0, 1, 2, 3],
            ),
            (
                "beam",
                Element2D::Beam(Beam2D::new(2, [10, 20]).expect("valid beam should be created")),
                vec![0, 1, 2, 3, 4, 5],
            ),
            (
                "triangle_t3",
                Element2D::TriangleT3(TriangleT3::new(3, [10, 20, 30]).expect("valid triangle should be created")),
                vec![0, 1, 2, 3, 4, 5],
            ),
        ];

        for (name, element, expected_indices) in cases {
            let mut model = Model2D::new();

            let nodes = [(10, 0.0, 0.0), (20, 1.0, 0.0), (30, 0.0, 1.0)];

            for (id, x, y) in nodes {
                let node = Node2D::new(id, x, y).expect("valid node should be created");

                model.add_node(node).expect("node should be added");
            }

            model.add_element(element).expect("element should be added");

            let numbering = DofNumbering2D::from_model(&model).expect("numbering should be created");

            let actual_indices = numbering.element_dof_indices(&element).expect("element DOFs should be mapped");

            assert_eq!(actual_indices, expected_indices, "failed case: {name}");
        }
    }

    #[test]
    fn maps_constraints_to_global_dofs() {
        let mut model = Model2D::new();

        model.add_node(Node2D::new(10, 0.0, 0.0).expect("valid node")).expect("node should be added");

        model.add_node(Node2D::new(20, 1.0, 0.0).expect("valid node")).expect("node should be added");

        model
            .add_element(Element2D::Beam(Beam2D::new(1, [10, 20]).expect("valid beam")))
            .expect("element should be added");

        let constraints = [(10, Dof2D::Ux, 0.0), (10, Dof2D::Rz, 0.01), (20, Dof2D::Uy, 0.0)];

        for (node_id, dof, displacement) in constraints {
            model
                .add_constraint(DisplacementConstraint2D::new(node_id, dof, displacement).expect("valid constraint"))
                .expect("constraint should be added");
        }

        let numbering = DofNumbering2D::from_model(&model).expect("numbering should be created");

        let actual = numbering.constraint_dof_indices(&model).expect("constraints should be mapped");

        assert_eq!(actual, vec![(0, 0.0), (2, 0.01), (4, 0.0),]);
    }
}
