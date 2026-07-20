use std::collections::HashMap;

use crate::FemError;
use crate::elements::BarElement1D;

use super::{
    constraint::DisplacementConstraint1D, load::NodalLoad1D, material::LinearElasticMaterial, node::Node1D,
    section::BarSection,
};

/// Complete, validated input model for a one-dimensional bar analysis.
///
/// External identifiers are kept intact, but every node is also mapped to a
/// contiguous zero-based internal index. Because a 1D bar node has one axial
/// displacement degree of freedom, that internal node index is also the global
/// displacement DOF index.
#[derive(Debug, Clone)]
pub struct BarModel1D {
    nodes: Vec<Node1D>,
    materials: Vec<LinearElasticMaterial>,
    sections: Vec<BarSection>,
    elements: Vec<BarElement1D>,
    loads: Vec<NodalLoad1D>,
    constraints: Vec<DisplacementConstraint1D>,
    node_id_to_index: HashMap<usize, usize>,
    material_id_to_index: HashMap<usize, usize>,
    section_id_to_index: HashMap<usize, usize>,
    element_id_to_index: HashMap<usize, usize>,
}

impl BarModel1D {
    /// Creates a validated 1D bar model and builds its external-id maps.
    pub fn new(
        nodes: Vec<Node1D>, materials: Vec<LinearElasticMaterial>, sections: Vec<BarSection>,
        elements: Vec<BarElement1D>, loads: Vec<NodalLoad1D>, constraints: Vec<DisplacementConstraint1D>,
    ) -> Result<Self, FemError> {
        let node_id_to_index = build_index_map("node", nodes.iter().map(Node1D::id))?;
        let material_id_to_index = build_index_map("material", materials.iter().map(LinearElasticMaterial::id))?;
        let section_id_to_index = build_index_map("section", sections.iter().map(BarSection::id))?;
        let element_id_to_index = build_index_map("element", elements.iter().map(BarElement1D::id))?;

        validate_element_references(&elements, &nodes, &node_id_to_index, &material_id_to_index, &section_id_to_index)?;
        validate_load_references(&loads, &node_id_to_index)?;
        validate_constraint_references(&constraints, &node_id_to_index)?;

        Ok(Self {
            nodes,
            materials,
            sections,
            elements,
            loads,
            constraints,
            node_id_to_index,
            material_id_to_index,
            section_id_to_index,
            element_id_to_index,
        })
    }

    /// Nodes in the same order in which they were provided to the model.
    #[must_use]
    pub fn nodes(&self) -> &[Node1D] {
        &self.nodes
    }

    /// Materials in the same order in which they were provided to the model.
    #[must_use]
    pub fn materials(&self) -> &[LinearElasticMaterial] {
        &self.materials
    }

    /// Sections in the same order in which they were provided to the model.
    #[must_use]
    pub fn sections(&self) -> &[BarSection] {
        &self.sections
    }

    /// Elements in the same order in which they were provided to the model.
    #[must_use]
    pub fn elements(&self) -> &[BarElement1D] {
        &self.elements
    }

    /// Nodal loads in the same order in which they were provided to the model.
    #[must_use]
    pub fn loads(&self) -> &[NodalLoad1D] {
        &self.loads
    }

    /// Displacement constraints in the same order in which they were provided to the model.
    #[must_use]
    pub fn constraints(&self) -> &[DisplacementConstraint1D] {
        &self.constraints
    }

    /// Number of displacement degrees of freedom.
    #[must_use]
    pub fn dof_count(&self) -> usize {
        self.nodes.len()
    }

    /// Returns the contiguous internal node index for an external node id.
    pub fn node_index(&self, node_id: usize) -> Result<usize, FemError> {
        lookup_index(&self.node_id_to_index, "node", node_id)
    }

    /// Returns the internal material index for an external material id.
    pub fn material_index(&self, material_id: usize) -> Result<usize, FemError> {
        lookup_index(&self.material_id_to_index, "material", material_id)
    }

    /// Returns the internal section index for an external section id.
    pub fn section_index(&self, section_id: usize) -> Result<usize, FemError> {
        lookup_index(&self.section_id_to_index, "section", section_id)
    }

    /// Returns the internal element index for an external element id.
    pub fn element_index(&self, element_id: usize) -> Result<usize, FemError> {
        lookup_index(&self.element_id_to_index, "element", element_id)
    }

    /// Returns the global displacement DOF index for an external node id.
    pub fn displacement_dof_index(&self, node_id: usize) -> Result<usize, FemError> {
        self.node_index(node_id)
    }

    /// Returns local-to-global DOF indices for an element's local vector `[u1, u2]^T`.
    pub fn element_dof_indices(&self, element_id: usize) -> Result<[usize; 2], FemError> {
        let element_index = self.element_index(element_id)?;
        let element = match self.elements.get(element_index) {
            Some(element) => element,
            None => {
                return Err(FemError::UnknownId { entity: "element", id: element_id });
            }
        };
        let node_ids = element.local_dof_node_ids();

        Ok([self.displacement_dof_index(node_ids[0])?, self.displacement_dof_index(node_ids[1])?])
    }
}

fn build_index_map<I>(entity: &'static str, ids: I) -> Result<HashMap<usize, usize>, FemError>
where
    I: IntoIterator<Item = usize>,
{
    let mut map = HashMap::new();
    for (index, id) in ids.into_iter().enumerate() {
        if map.insert(id, index).is_some() {
            return Err(FemError::DuplicateId { entity, id });
        }
    }

    Ok(map)
}

fn lookup_index(map: &HashMap<usize, usize>, entity: &'static str, id: usize) -> Result<usize, FemError> {
    match map.get(&id) {
        Some(index) => Ok(*index),
        None => Err(FemError::UnknownId { entity, id }),
    }
}

fn validate_element_references(
    elements: &[BarElement1D], nodes: &[Node1D], node_id_to_index: &HashMap<usize, usize>,
    material_id_to_index: &HashMap<usize, usize>, section_id_to_index: &HashMap<usize, usize>,
) -> Result<(), FemError> {
    for element in elements {
        let node_ids = element.node_ids();
        for node_id in node_ids {
            if !node_id_to_index.contains_key(&node_id) {
                return Err(FemError::InvalidReference {
                    owner_entity: "element",
                    owner_id: element.id(),
                    referenced_entity: "node",
                    referenced_id: node_id,
                });
            }
        }

        if !material_id_to_index.contains_key(&element.material_id()) {
            return Err(FemError::InvalidReference {
                owner_entity: "element",
                owner_id: element.id(),
                referenced_entity: "material",
                referenced_id: element.material_id(),
            });
        }

        if !section_id_to_index.contains_key(&element.section_id()) {
            return Err(FemError::InvalidReference {
                owner_entity: "element",
                owner_id: element.id(),
                referenced_entity: "section",
                referenced_id: element.section_id(),
            });
        }

        let first = match node_id_to_index.get(&node_ids[0]).and_then(|index| nodes.get(*index)) {
            Some(node) => node,
            None => {
                return Err(FemError::InvalidReference {
                    owner_entity: "element",
                    owner_id: element.id(),
                    referenced_entity: "node",
                    referenced_id: node_ids[0],
                });
            }
        };
        let second = match node_id_to_index.get(&node_ids[1]).and_then(|index| nodes.get(*index)) {
            Some(node) => node,
            None => {
                return Err(FemError::InvalidReference {
                    owner_entity: "element",
                    owner_id: element.id(),
                    referenced_entity: "node",
                    referenced_id: node_ids[1],
                });
            }
        };

        element.length(first, second)?;
    }

    Ok(())
}

fn validate_load_references(loads: &[NodalLoad1D], node_id_to_index: &HashMap<usize, usize>) -> Result<(), FemError> {
    for load in loads {
        if !node_id_to_index.contains_key(&load.node_id()) {
            return Err(FemError::InvalidReference {
                owner_entity: "nodal load",
                owner_id: load.node_id(),
                referenced_entity: "node",
                referenced_id: load.node_id(),
            });
        }
    }

    Ok(())
}

fn validate_constraint_references(
    constraints: &[DisplacementConstraint1D], node_id_to_index: &HashMap<usize, usize>,
) -> Result<(), FemError> {
    for constraint in constraints {
        if !node_id_to_index.contains_key(&constraint.node_id()) {
            return Err(FemError::InvalidReference {
                owner_entity: "displacement constraint",
                owner_id: constraint.node_id(),
                referenced_entity: "node",
                referenced_id: constraint.node_id(),
            });
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn node(id: usize, x: f64) -> Node1D {
        Node1D::new(id, x).unwrap()
    }

    fn material(id: usize) -> LinearElasticMaterial {
        LinearElasticMaterial::new(id, 210.0e9, 0.3, 7850.0).unwrap()
    }

    fn section(id: usize) -> BarSection {
        BarSection::new(id, 1.0e-4).unwrap()
    }

    fn element(id: usize, node_ids: [usize; 2], material_id: usize, section_id: usize) -> BarElement1D {
        BarElement1D::new(id, node_ids, material_id, section_id).unwrap()
    }

    #[test]
    fn maps_external_node_ids_to_contiguous_dofs_in_input_order() {
        let model = BarModel1D::new(
            vec![node(10, 0.0), node(40, 1.0), node(99, 2.0)],
            vec![material(7)],
            vec![section(3)],
            vec![element(20, [40, 99], 7, 3)],
            vec![NodalLoad1D::new(99, 100.0).unwrap()],
            vec![DisplacementConstraint1D::new(10, 0.0).unwrap()],
        )
        .unwrap();

        assert_eq!(model.dof_count(), 3);
        assert_eq!(model.node_index(10).unwrap(), 0);
        assert_eq!(model.node_index(40).unwrap(), 1);
        assert_eq!(model.displacement_dof_index(99).unwrap(), 2);
        assert_eq!(model.material_index(7).unwrap(), 0);
        assert_eq!(model.section_index(3).unwrap(), 0);
        assert_eq!(model.element_index(20).unwrap(), 0);
        assert_eq!(model.element_dof_indices(20).unwrap(), [1, 2]);
        assert_eq!(model.nodes().len(), 3);
        assert_eq!(model.materials().len(), 1);
        assert_eq!(model.sections().len(), 1);
        assert_eq!(model.elements().len(), 1);
        assert_eq!(model.loads().len(), 1);
        assert_eq!(model.constraints().len(), 1);
    }

    #[test]
    fn rejects_duplicate_node_ids() {
        let error = BarModel1D::new(
            vec![node(10, 0.0), node(10, 1.0)],
            vec![material(7)],
            vec![section(3)],
            Vec::new(),
            Vec::new(),
            Vec::new(),
        )
        .unwrap_err();

        assert_eq!(error, FemError::DuplicateId { entity: "node", id: 10 });
    }

    #[test]
    fn rejects_duplicate_element_ids() {
        let error = BarModel1D::new(
            vec![node(10, 0.0), node(40, 1.0)],
            vec![material(7)],
            vec![section(3)],
            vec![element(20, [10, 40], 7, 3), element(20, [10, 40], 7, 3)],
            Vec::new(),
            Vec::new(),
        )
        .unwrap_err();

        assert_eq!(error, FemError::DuplicateId { entity: "element", id: 20 });
    }

    #[test]
    fn rejects_element_with_missing_node() {
        let error = BarModel1D::new(
            vec![node(10, 0.0), node(40, 1.0)],
            vec![material(7)],
            vec![section(3)],
            vec![element(20, [40, 99], 7, 3)],
            Vec::new(),
            Vec::new(),
        )
        .unwrap_err();

        assert_eq!(
            error,
            FemError::InvalidReference {
                owner_entity: "element",
                owner_id: 20,
                referenced_entity: "node",
                referenced_id: 99,
            }
        );
    }

    #[test]
    fn rejects_element_with_missing_material() {
        let error = BarModel1D::new(
            vec![node(10, 0.0), node(40, 1.0)],
            vec![material(7)],
            vec![section(3)],
            vec![element(20, [10, 40], 99, 3)],
            Vec::new(),
            Vec::new(),
        )
        .unwrap_err();

        assert_eq!(
            error,
            FemError::InvalidReference {
                owner_entity: "element",
                owner_id: 20,
                referenced_entity: "material",
                referenced_id: 99,
            }
        );
    }

    #[test]
    fn rejects_element_with_missing_section() {
        let error = BarModel1D::new(
            vec![node(10, 0.0), node(40, 1.0)],
            vec![material(7)],
            vec![section(3)],
            vec![element(20, [10, 40], 7, 99)],
            Vec::new(),
            Vec::new(),
        )
        .unwrap_err();

        assert_eq!(
            error,
            FemError::InvalidReference {
                owner_entity: "element",
                owner_id: 20,
                referenced_entity: "section",
                referenced_id: 99,
            }
        );
    }

    #[test]
    fn rejects_element_with_zero_length_geometry() {
        let error = BarModel1D::new(
            vec![node(10, 1.0), node(40, 1.0)],
            vec![material(7)],
            vec![section(3)],
            vec![element(20, [10, 40], 7, 3)],
            Vec::new(),
            Vec::new(),
        )
        .unwrap_err();

        assert_eq!(error, FemError::ZeroLengthElement { element_id: 20, node_ids: [10, 40], length: 0.0 });
    }

    #[test]
    fn rejects_load_with_missing_node() {
        let error = BarModel1D::new(
            vec![node(10, 0.0), node(40, 1.0)],
            vec![material(7)],
            vec![section(3)],
            vec![element(20, [10, 40], 7, 3)],
            vec![NodalLoad1D::new(99, 1.0).unwrap()],
            Vec::new(),
        )
        .unwrap_err();

        assert_eq!(
            error,
            FemError::InvalidReference {
                owner_entity: "nodal load",
                owner_id: 99,
                referenced_entity: "node",
                referenced_id: 99,
            }
        );
    }

    #[test]
    fn rejects_constraint_with_missing_node() {
        let error = BarModel1D::new(
            vec![node(10, 0.0), node(40, 1.0)],
            vec![material(7)],
            vec![section(3)],
            vec![element(20, [10, 40], 7, 3)],
            Vec::new(),
            vec![DisplacementConstraint1D::new(99, 0.0).unwrap()],
        )
        .unwrap_err();

        assert_eq!(
            error,
            FemError::InvalidReference {
                owner_entity: "displacement constraint",
                owner_id: 99,
                referenced_entity: "node",
                referenced_id: 99,
            }
        );
    }

    #[test]
    fn returns_unknown_id_when_dof_lookup_misses() {
        let model = BarModel1D::new(
            vec![node(10, 0.0), node(40, 1.0)],
            vec![material(7)],
            vec![section(3)],
            vec![element(20, [10, 40], 7, 3)],
            Vec::new(),
            Vec::new(),
        )
        .unwrap();

        let error = model.displacement_dof_index(99).unwrap_err();

        assert_eq!(error, FemError::UnknownId { entity: "node", id: 99 });
    }
}
