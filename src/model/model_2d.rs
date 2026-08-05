//! Defines a 2D model consisting of nodes and displacement constraints.

use super::AnalysisSettings2D;
use super::BoundaryConditions2D;
use super::Loads2D;
use super::Materials2D;
use super::Mesh2D;
use super::Sections2D;
use super::constraint::DisplacementConstraint2D;
use super::material_model::Material2D;
use super::nodal_load::NodalLoad2D;
use super::node::Node2D;
use super::sections::{BeamSection2D, PlaneStressSection2D, Section2D, TrussSection2D};
use crate::elements::Element2D;
use crate::error::FemError;

/// Represents a 2D model consisting of nodes and displacement constraints.
#[derive(Default)]
pub struct Model2D {
    mesh: Mesh2D,
    boundary_conditions: BoundaryConditions2D,
    materials: Materials2D,
    sections: Sections2D,
    loads: Loads2D,
    analysis_settings: AnalysisSettings2D,
}

impl Model2D {
    /// Creates a new empty 2D model.
    pub fn new() -> Self {
        Self::default()
    }

    /// Adds a node to the model. Returns an error if a node with the same ID already exists.
    pub fn add_node(&mut self, node: Node2D) -> Result<(), FemError> {
        if self.mesh.contains_node_id(node.id()) {
            return Err(FemError::DuplicateId { entity: "node", id: node.id() });
        }

        self.mesh.push_node(node);

        Ok(())
    }

    /// Adds a displacement constraint to the model. Returns an error if the node ID associated with the constraint does not exist in the model.
    pub fn add_constraint(&mut self, constraint: DisplacementConstraint2D) -> Result<(), FemError> {
        if !self.mesh.contains_node_id(constraint.node_id()) {
            return Err(FemError::UnknownId { entity: "node", id: constraint.node_id() });
        }

        self.boundary_conditions.push_displacement_constraint(constraint);

        Ok(())
    }

    /// Adds an element to the model. Returns an error if an element with the same ID already exists
    /// or if any of the node IDs associated with the element do not exist in the model.
    pub fn add_element(&mut self, element: Element2D) -> Result<(), FemError> {
        self.validate_element_for_insert(&element)?;

        self.materials.material(element.material_id())?;
        let section = self.sections.section(element.section_id())?;
        element.validate_section(section)?;

        self.mesh.push_element(element);

        Ok(())
    }

    /// Adds an element together with its referenced section.
    pub fn add_element_with_section(&mut self, element: Element2D, section: Section2D) -> Result<(), FemError> {
        self.validate_element_for_insert(&element)?;
        self.materials.material(element.material_id())?;
        element.validate_section(&section)?;

        self.sections.add_section(element.section_id(), section)?;
        self.mesh.push_element(element);

        Ok(())
    }

    fn validate_element_for_insert(&self, element: &Element2D) -> Result<(), FemError> {
        if self.mesh.contains_element_id(element.id()) {
            return Err(FemError::DuplicateId { entity: "element", id: element.id() });
        }

        for &node_id in element.node_ids() {
            self.find_node(node_id)?;
        }

        self.validate_element_geometry(element)?;

        Ok(())
    }

    /// Adds a section to the model.
    pub fn add_section(&mut self, section_id: usize, section: Section2D) -> Result<(), FemError> {
        self.sections.add_section(section_id, section)
    }

    /// Adds a material to the model.
    pub fn add_material(&mut self, material_id: usize, material: Material2D) -> Result<(), FemError> {
        self.materials.add_material(material_id, material)
    }

    /// Sets the material properties for the model.
    pub fn set_material(&mut self, material: Material2D) {
        self.materials.set_default_material(material);
    }

    /// Finds a node in the model by its ID. Returns a reference to the node if found, or an error if the node ID does not exist in the model.
    fn find_node(&self, node_id: usize) -> Result<&Node2D, FemError> {
        self.mesh.node(node_id).ok_or(FemError::UnknownId { entity: "node", id: node_id })
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

    pub fn add_load(&mut self, load: NodalLoad2D) -> Result<(), FemError> {
        if !self.mesh.contains_node_id(load.node_id()) {
            return Err(FemError::UnknownId { entity: "node", id: load.node_id() });
        }

        self.loads.push_nodal_load(load);

        Ok(())
    }

    /// Returns a slice of all nodes in the model.
    #[must_use]
    pub fn nodes(&self) -> &[Node2D] {
        self.mesh.nodes()
    }

    /// Returns a slice of all displacement constraints in the model.
    #[must_use]
    pub fn constraints(&self) -> &[DisplacementConstraint2D] {
        self.boundary_conditions.displacement_constraints()
    }

    /// Returns a slice of all elements in the model.
    #[must_use]
    pub fn elements(&self) -> &[Element2D] {
        self.mesh.elements()
    }

    /// Returns all section definitions.
    #[must_use]
    pub fn sections(&self) -> &Sections2D {
        &self.sections
    }

    /// Returns all material definitions.
    #[must_use]
    pub fn materials(&self) -> &Materials2D {
        &self.materials
    }

    /// Returns a section by ID.
    pub fn section(&self, section_id: usize) -> Result<&Section2D, FemError> {
        self.sections.section(section_id)
    }

    /// Returns a truss section by ID.
    pub fn truss_section(&self, section_id: usize) -> Result<&TrussSection2D, FemError> {
        self.sections.truss_section(section_id)
    }

    /// Returns a beam section by ID.
    pub fn beam_section(&self, section_id: usize) -> Result<&BeamSection2D, FemError> {
        self.sections.beam_section(section_id)
    }

    /// Returns a plane-stress section by ID.
    pub fn plane_stress_section(&self, section_id: usize) -> Result<&PlaneStressSection2D, FemError> {
        self.sections.plane_stress_section(section_id)
    }

    /// Returns a material by ID.
    pub fn material(&self, material_id: usize) -> Result<&Material2D, FemError> {
        self.materials.material(material_id)
    }

    /// Returns the default material if one was set through `set_material`.
    #[must_use]
    pub fn default_material(&self) -> Option<&Material2D> {
        self.materials.default_material()
    }

    /// Returns a slice of all nodal loads in the model.
    #[must_use]
    pub fn loads(&self) -> &[NodalLoad2D] {
        self.loads.nodal_loads()
    }

    /// Returns the analysis settings for the model.
    #[must_use]
    pub fn analysis_settings(&self) -> &AnalysisSettings2D {
        &self.analysis_settings
    }

    /// Returns mutable analysis settings for the model.
    pub fn analysis_settings_mut(&mut self) -> &mut AnalysisSettings2D {
        &mut self.analysis_settings
    }
}

#[cfg(test)]
mod tests {
    use super::Model2D;
    use crate::FemError;
    use crate::elements::{Beam2D, Element2D, TriangleT3, Truss2D};
    use crate::model::{
        DEFAULT_MATERIAL_ID, DisplacementConstraint2D, Dof2D, Material2D, NodalLoad2D, Node2D, Section2D, SolverKind2D,
        TrussSection2D,
    };

    #[test]
    fn creates_empty_model() {
        let model = Model2D::new();

        assert!(model.nodes().is_empty());
        assert!(model.constraints().is_empty());
        assert!(model.elements().is_empty());
        assert!(model.sections().sections().is_empty());
        assert!(model.materials().materials().is_empty());
        assert!(model.default_material().is_none());
        assert_eq!(model.analysis_settings().solver(), SolverKind2D::Dense);
    }

    #[test]
    fn updates_analysis_settings() {
        let mut model = Model2D::new();

        model.analysis_settings_mut().set_solver(SolverKind2D::Sparse);
        model.analysis_settings_mut().set_cg_max_iterations(250).expect("iteration limit should be valid");
        model.analysis_settings_mut().set_cg_tolerance(1e-8).expect("tolerance should be valid");

        assert_eq!(model.analysis_settings().solver(), SolverKind2D::Sparse);
        assert_eq!(model.analysis_settings().cg_max_iterations(), 250);
        assert_eq!(model.analysis_settings().cg_tolerance(), 1e-8);
    }

    #[test]
    fn adds_node_to_model() {
        let mut model = Model2D::new();
        model.set_material(Material2D::new(200.0, 0.3, 1.0).expect("valid material should be created"));
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
    fn rejects_load_for_unknown_node() {
        let mut model = Model2D::new();
        let load = NodalLoad2D::new(99, Dof2D::Ux, 10.0).expect("valid load should be created");

        let result = model.add_load(load);

        assert!(matches!(result, Err(FemError::UnknownId { entity: "node", id: 99 })));
        assert!(model.loads().is_empty());
    }

    #[test]
    fn allows_multiple_loads_on_the_same_dof() {
        let mut model = Model2D::new();
        let node = Node2D::new(7, 0.0, 0.0).expect("valid node should be created");

        model.add_node(node).expect("node should be added");

        let first_load = NodalLoad2D::new(7, Dof2D::Ux, 10.0).expect("valid load should be created");
        let second_load = NodalLoad2D::new(7, Dof2D::Ux, -3.0).expect("valid load should be created");

        model.add_load(first_load).expect("first load should be added");
        model.add_load(second_load).expect("second load should be added");

        assert_eq!(model.loads().len(), 2);
        assert_eq!(model.loads()[0], first_load);
        assert_eq!(model.loads()[1], second_load);
    }

    #[test]
    fn adds_element_for_existing_nodes() {
        let mut model = Model2D::new();
        model.set_material(Material2D::new(200.0, 0.3, 1.0).expect("valid material should be created"));

        let first_node = Node2D::new(1, 0.0, 0.0).expect("valid node should be created");

        let second_node = Node2D::new(2, 1.0, 0.0).expect("valid node should be created");

        model.add_node(first_node).expect("node should be added");
        model.add_node(second_node).expect("node should be added");

        let truss = Truss2D::new(10, [1, 2], DEFAULT_MATERIAL_ID, 100).expect("valid truss should be created");
        let element = Element2D::Truss(truss);
        let section = Section2D::Truss(TrussSection2D::new(1.0).expect("valid section should be created"));

        model.add_element_with_section(element, section).expect("element should be added");

        assert_eq!(model.elements().len(), 1);
        assert_eq!(model.elements()[0].id(), 10);
        assert_eq!(model.elements()[0].node_ids(), &[1, 2]);
        assert_eq!(model.elements()[0].material_id(), DEFAULT_MATERIAL_ID);
        assert_eq!(model.elements()[0].section_id(), 100);
        assert!(model.truss_section(100).is_ok());
    }

    #[test]
    fn rejects_element_for_unknown_material() {
        let mut model = Model2D::new();

        model.add_node(Node2D::new(1, 0.0, 0.0).expect("valid node should be created")).expect("node should be added");
        model.add_node(Node2D::new(2, 1.0, 0.0).expect("valid node should be created")).expect("node should be added");

        let truss = Truss2D::new(10, [1, 2], 99, 100).expect("valid truss should be created");
        let section = Section2D::Truss(TrussSection2D::new(1.0).expect("valid section should be created"));
        let result = model.add_element_with_section(Element2D::Truss(truss), section);

        assert!(matches!(result, Err(FemError::UnknownId { entity: "material", id: 99 })));
        assert!(model.elements().is_empty());
        assert!(model.sections().sections().is_empty());
    }

    #[test]
    fn rejects_element_for_unknown_node() {
        let mut model = Model2D::new();

        let node = Node2D::new(1, 0.0, 0.0).expect("valid node should be created");

        model.add_node(node).expect("node should be added");

        let element = Element2D::Truss(
            Truss2D::new(10, [1, 99], DEFAULT_MATERIAL_ID, 100).expect("valid truss should be created"),
        );

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
                Element2D::Truss(
                    Truss2D::new(10, [1, 2], DEFAULT_MATERIAL_ID, 100)
                        .expect("valid truss connectivity should be created"),
                ),
                "length",
            ),
            (
                "beam",
                Element2D::Beam(
                    Beam2D::new(20, [1, 2], DEFAULT_MATERIAL_ID, 200)
                        .expect("valid beam connectivity should be created"),
                ),
                "length",
            ),
            (
                "triangle_t3",
                Element2D::TriangleT3(
                    TriangleT3::new(30, [1, 2, 3], DEFAULT_MATERIAL_ID, 300)
                        .expect("valid triangle connectivity should be created"),
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

        assert_eq!(model.default_material(), Some(&second));
        assert_eq!(model.material(DEFAULT_MATERIAL_ID).expect("default material should exist"), &second);
    }
}
