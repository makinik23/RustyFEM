//! Serializable input format for building a `Model2D`.

use serde::{Deserialize, Serialize};

use crate::FemError;
use crate::elements::{Beam2D, Element2D, QuadQ4, QuadQ8, TriangleT3, TriangleT6, Truss2D};
use crate::model::{
    BeamSection2D, BeamUniformLineLoad2D, BodyForce2D, DEFAULT_MATERIAL_ID, DisplacementConstraint2D, Dof2D,
    EdgeTraction2D, ElementLoad2D, LoadCoordinateSystem2D, Material2D, Model2D, NodalLoad2D, Node2D,
    PlaneStressSection2D, Section2D, SelfWeight2D, SolverKind2D, TrussSection2D,
};

/// Complete JSON-friendly definition of a two-dimensional FEM model.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Model2DInput {
    #[serde(default)]
    pub analysis_settings: AnalysisSettings2DInput,
    #[serde(default)]
    pub materials: Vec<Material2DInput>,
    #[serde(default)]
    pub sections: Vec<Section2DInput>,
    #[serde(default)]
    pub nodes: Vec<Node2DInput>,
    #[serde(default)]
    pub elements: Vec<ElementType2DInput>,
    #[serde(default)]
    pub constraints: Vec<DisplacementConstraint2DInput>,
    #[serde(default)]
    pub loads: Vec<ElementLoad2DInput>,
    #[serde(default)]
    pub nodal_loads: Vec<NodalLoad2DInput>,
}

/// Optional analysis settings. Missing values use `AnalysisSettings2D` defaults.
#[derive(Debug, Clone, Copy, Default, PartialEq, Serialize, Deserialize)]
pub struct AnalysisSettings2DInput {
    pub solver: Option<SolverKind2DInput>,
    pub cg_tolerance: Option<f64>,
    pub cg_max_iterations: Option<usize>,
    pub cg_stagnation_window: Option<usize>,
    pub cg_stagnation_tolerance: Option<f64>,
}

/// Solver selector used by the serialized model format.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum SolverKind2DInput {
    Dense,
    Sparse,
}

/// Material definition.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Material2DInput {
    pub id: usize,
    pub young_modulus: f64,
    pub poisson_ratio: f64,
    pub density: f64,
}

/// Section definition with a tagged concrete section kind.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Section2DInput {
    pub id: usize,
    #[serde(flatten)]
    pub kind: Section2DInputKind,
}

/// Concrete section data.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(tag = "type", rename_all = "snake_case")]
pub enum Section2DInputKind {
    Truss(TrussSection2DInput),
    Beam(BeamSection2DInput),
    PlaneStress(PlaneStressSection2DInput),
}

/// Truss section input.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct TrussSection2DInput {
    pub area: f64,
}

/// Beam section input.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct BeamSection2DInput {
    pub area: f64,
    pub second_moment_of_area: f64,
    pub height: Option<f64>,
}

/// Plane-stress section input.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct PlaneStressSection2DInput {
    pub thickness: f64,
}

/// Node definition.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Node2DInput {
    pub id: usize,
    pub x: f64,
    pub y: f64,
}

/// Displacement constraint definition.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct DisplacementConstraint2DInput {
    pub node: usize,
    pub dof: Dof2DInput,
    pub value: f64,
}

/// Nodal load definition.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct NodalLoad2DInput {
    pub node: usize,
    pub dof: Dof2DInput,
    pub value: f64,
}

/// Serialized 2D degree of freedom.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum Dof2DInput {
    #[serde(rename = "ux", alias = "Ux", alias = "UX")]
    Ux,
    #[serde(rename = "uy", alias = "Uy", alias = "UY")]
    Uy,
    #[serde(rename = "rz", alias = "Rz", alias = "RZ")]
    Rz,
}

/// Element definition with a tagged concrete element kind.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(tag = "type", rename_all = "snake_case")]
pub enum ElementType2DInput {
    Truss {
        id: usize,
        nodes: [usize; 2],
        material: usize,
        section: usize,
    },
    Beam {
        id: usize,
        nodes: [usize; 2],
        material: usize,
        section: usize,
    },
    #[serde(rename = "triangle", alias = "t3")]
    TriangleT3 {
        id: usize,
        nodes: [usize; 3],
        material: usize,
        section: usize,
    },
    T6 {
        id: usize,
        nodes: [usize; 6],
        material: usize,
        section: usize,
    },
    Q4 {
        id: usize,
        nodes: [usize; 4],
        material: usize,
        section: usize,
    },
    Q8 {
        id: usize,
        nodes: [usize; 8],
        material: usize,
        section: usize,
    },
}

/// Element load definition with a tagged concrete load kind.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct ElementLoad2DInput {
    #[serde(flatten)]
    pub kind: ElementLoad2DInputKind,
}

/// Concrete element load data.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(tag = "type", rename_all = "snake_case")]
pub enum ElementLoad2DInputKind {
    #[serde(alias = "beam_uniform_line")]
    BeamUniform {
        element: usize,
        coordinate_system: LoadCoordinateSystem2DInput,
        qx: f64,
        qy: f64,
    },
    EdgeTraction {
        element: usize,
        edge: [usize; 2],
        coordinate_system: LoadCoordinateSystem2DInput,
        tx: f64,
        ty: f64,
    },
    BodyForce {
        element: usize,
        bx: f64,
        by: f64,
    },
    SelfWeight {
        element: usize,
        ax: f64,
        ay: f64,
    },
    Nodal {
        node: usize,
        dof: Dof2DInput,
        value: f64,
    },
}

/// Coordinate system used by serialized element loads.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum LoadCoordinateSystem2DInput {
    Global,
    Local,
}

impl Model2DInput {
    /// Creates a serializable snapshot of an existing validated model.
    #[must_use]
    pub fn from_model(model: &Model2D) -> Self {
        Self {
            analysis_settings: AnalysisSettings2DInput::from_model(model),
            materials: model
                .materials()
                .materials()
                .iter()
                .map(|&(id, material)| Material2DInput {
                    id,
                    young_modulus: material.young_modulus(),
                    poisson_ratio: material.poisson_ratio(),
                    density: material.density(),
                })
                .collect(),
            sections: model
                .sections()
                .sections()
                .iter()
                .map(|&(id, section)| Section2DInput { id, kind: section.into() })
                .collect(),
            nodes: model.nodes().iter().map(|node| Node2DInput { id: node.id(), x: node.x(), y: node.y() }).collect(),
            elements: model.elements().iter().map(ElementType2DInput::from).collect(),
            constraints: model
                .constraints()
                .iter()
                .map(|constraint| DisplacementConstraint2DInput {
                    node: constraint.node_id(),
                    dof: constraint.dof().into(),
                    value: constraint.displacement(),
                })
                .collect(),
            loads: model.element_loads().iter().map(ElementLoad2DInput::from).collect(),
            nodal_loads: model
                .loads()
                .iter()
                .map(|load| NodalLoad2DInput { node: load.node_id(), dof: load.dof().into(), value: load.value() })
                .collect(),
        }
    }

    /// Builds a validated `Model2D` using the same insertion API as hand-written models.
    pub fn into_model(self) -> Result<Model2D, FemError> {
        let mut model = Model2D::new();

        self.analysis_settings.apply_to(&mut model)?;

        for material in self.materials {
            model.add_material(material.id, material.try_into()?)?;
        }

        if model.materials().materials().is_empty() {
            model.add_material(DEFAULT_MATERIAL_ID, Material2D::new(1.0, 0.3, 0.0)?)?;
        }

        for section in self.sections {
            model.add_section(section.id, section.kind.try_into()?)?;
        }

        for node in self.nodes {
            model.add_node(node.try_into()?)?;
        }

        for element in self.elements {
            model.add_element(element.try_into()?)?;
        }

        for constraint in self.constraints {
            model.add_constraint(constraint.try_into()?)?;
        }

        for load in self.nodal_loads {
            model.add_load(load.try_into()?)?;
        }

        for load in self.loads {
            match load.kind.try_into()? {
                Load2DInputConversion::Nodal(load) => model.add_load(load)?,
                Load2DInputConversion::Element(load) => model.add_element_load(load)?,
            }
        }

        Ok(model)
    }
}

impl AnalysisSettings2DInput {
    fn from_model(model: &Model2D) -> Self {
        let settings = model.analysis_settings();

        Self {
            solver: Some(settings.solver().into()),
            cg_tolerance: Some(settings.cg_tolerance()),
            cg_max_iterations: Some(settings.cg_max_iterations()),
            cg_stagnation_window: Some(settings.cg_stagnation_window()),
            cg_stagnation_tolerance: Some(settings.cg_stagnation_tolerance()),
        }
    }
}

impl From<SolverKind2D> for SolverKind2DInput {
    fn from(value: SolverKind2D) -> Self {
        match value {
            SolverKind2D::Dense => Self::Dense,
            SolverKind2D::Sparse => Self::Sparse,
        }
    }
}

impl From<Section2D> for Section2DInputKind {
    fn from(value: Section2D) -> Self {
        match value {
            Section2D::Truss(section) => Self::Truss(TrussSection2DInput { area: section.cross_section_area() }),
            Section2D::Beam(section) => Self::Beam(BeamSection2DInput {
                area: section.cross_section_area(),
                second_moment_of_area: section.second_moment_of_area(),
                height: section.section_height(),
            }),
            Section2D::PlaneStress(section) => {
                Self::PlaneStress(PlaneStressSection2DInput { thickness: section.thickness() })
            }
        }
    }
}

impl From<Dof2D> for Dof2DInput {
    fn from(value: Dof2D) -> Self {
        match value {
            Dof2D::Ux => Self::Ux,
            Dof2D::Uy => Self::Uy,
            Dof2D::Rz => Self::Rz,
        }
    }
}

impl From<&Element2D> for ElementType2DInput {
    fn from(value: &Element2D) -> Self {
        let id = value.id();
        let nodes = value.node_ids();
        let material = value.material_id();
        let section = value.section_id();

        match value {
            Element2D::Truss(_) => Self::Truss { id, nodes: [nodes[0], nodes[1]], material, section },
            Element2D::Beam(_) => Self::Beam { id, nodes: [nodes[0], nodes[1]], material, section },
            Element2D::TriangleT3(_) => {
                Self::TriangleT3 { id, nodes: [nodes[0], nodes[1], nodes[2]], material, section }
            }
            Element2D::TriangleT6(_) => {
                Self::T6 { id, nodes: nodes.try_into().expect("T6 connectivity has six nodes"), material, section }
            }
            Element2D::QuadQ4(_) => {
                Self::Q4 { id, nodes: nodes.try_into().expect("Q4 connectivity has four nodes"), material, section }
            }
            Element2D::QuadQ8(_) => {
                Self::Q8 { id, nodes: nodes.try_into().expect("Q8 connectivity has eight nodes"), material, section }
            }
        }
    }
}

impl From<&ElementLoad2D> for ElementLoad2DInput {
    fn from(value: &ElementLoad2D) -> Self {
        let kind = match value {
            ElementLoad2D::BeamUniformLine(load) => ElementLoad2DInputKind::BeamUniform {
                element: load.element_id(),
                coordinate_system: load.coordinate_system().into(),
                qx: load.x_component(),
                qy: load.y_component(),
            },
            ElementLoad2D::EdgeTraction(load) => ElementLoad2DInputKind::EdgeTraction {
                element: load.element_id(),
                edge: load.edge_node_ids(),
                coordinate_system: load.coordinate_system().into(),
                tx: load.x_component(),
                ty: load.y_component(),
            },
            ElementLoad2D::BodyForce(load) => ElementLoad2DInputKind::BodyForce {
                element: load.element_id(),
                bx: load.x_component(),
                by: load.y_component(),
            },
            ElementLoad2D::SelfWeight(load) => ElementLoad2DInputKind::SelfWeight {
                element: load.element_id(),
                ax: load.x_acceleration(),
                ay: load.y_acceleration(),
            },
        };

        Self { kind }
    }
}

impl From<LoadCoordinateSystem2D> for LoadCoordinateSystem2DInput {
    fn from(value: LoadCoordinateSystem2D) -> Self {
        match value {
            LoadCoordinateSystem2D::Global => Self::Global,
            LoadCoordinateSystem2D::Local => Self::Local,
        }
    }
}

impl AnalysisSettings2DInput {
    fn apply_to(self, model: &mut Model2D) -> Result<(), FemError> {
        let settings = model.analysis_settings_mut();

        if let Some(solver) = self.solver {
            settings.set_solver(solver.into());
        }

        if let Some(value) = self.cg_tolerance {
            settings.set_cg_tolerance(value)?;
        }

        if let Some(value) = self.cg_max_iterations {
            settings.set_cg_max_iterations(value)?;
        }

        if let Some(value) = self.cg_stagnation_window {
            settings.set_cg_stagnation_window(value);
        }

        if let Some(value) = self.cg_stagnation_tolerance {
            settings.set_cg_stagnation_tolerance(value)?;
        }

        Ok(())
    }
}

impl From<SolverKind2DInput> for SolverKind2D {
    fn from(value: SolverKind2DInput) -> Self {
        match value {
            SolverKind2DInput::Dense => Self::Dense,
            SolverKind2DInput::Sparse => Self::Sparse,
        }
    }
}

impl TryFrom<Material2DInput> for Material2D {
    type Error = FemError;

    fn try_from(value: Material2DInput) -> Result<Self, Self::Error> {
        Self::new(value.young_modulus, value.poisson_ratio, value.density)
    }
}

impl TryFrom<Section2DInputKind> for Section2D {
    type Error = FemError;

    fn try_from(value: Section2DInputKind) -> Result<Self, Self::Error> {
        match value {
            Section2DInputKind::Truss(section) => Ok(Self::Truss(TrussSection2D::new(section.area)?)),
            Section2DInputKind::Beam(section) => match section.height {
                Some(height) => Ok(Self::Beam(BeamSection2D::new_with_section_height(
                    section.area,
                    section.second_moment_of_area,
                    height,
                )?)),
                None => Ok(Self::Beam(BeamSection2D::new(section.area, section.second_moment_of_area)?)),
            },
            Section2DInputKind::PlaneStress(section) => {
                Ok(Self::PlaneStress(PlaneStressSection2D::new(section.thickness)?))
            }
        }
    }
}

impl TryFrom<Node2DInput> for Node2D {
    type Error = FemError;

    fn try_from(value: Node2DInput) -> Result<Self, Self::Error> {
        Self::new(value.id, value.x, value.y)
    }
}

impl TryFrom<DisplacementConstraint2DInput> for DisplacementConstraint2D {
    type Error = FemError;

    fn try_from(value: DisplacementConstraint2DInput) -> Result<Self, Self::Error> {
        Self::new(value.node, value.dof.into(), value.value)
    }
}

impl TryFrom<NodalLoad2DInput> for NodalLoad2D {
    type Error = FemError;

    fn try_from(value: NodalLoad2DInput) -> Result<Self, Self::Error> {
        Self::new(value.node, value.dof.into(), value.value)
    }
}

impl From<Dof2DInput> for Dof2D {
    fn from(value: Dof2DInput) -> Self {
        match value {
            Dof2DInput::Ux => Self::Ux,
            Dof2DInput::Uy => Self::Uy,
            Dof2DInput::Rz => Self::Rz,
        }
    }
}

impl TryFrom<ElementType2DInput> for Element2D {
    type Error = FemError;

    fn try_from(value: ElementType2DInput) -> Result<Self, Self::Error> {
        match value {
            ElementType2DInput::Truss { id, nodes, material, section } => {
                Ok(Self::Truss(Truss2D::new(id, nodes, material, section)?))
            }
            ElementType2DInput::Beam { id, nodes, material, section } => {
                Ok(Self::Beam(Beam2D::new(id, nodes, material, section)?))
            }
            ElementType2DInput::TriangleT3 { id, nodes, material, section } => {
                Ok(Self::TriangleT3(TriangleT3::new(id, nodes, material, section)?))
            }
            ElementType2DInput::T6 { id, nodes, material, section } => {
                Ok(Self::TriangleT6(TriangleT6::new(id, nodes, material, section)?))
            }
            ElementType2DInput::Q4 { id, nodes, material, section } => {
                Ok(Self::QuadQ4(QuadQ4::new(id, nodes, material, section)?))
            }
            ElementType2DInput::Q8 { id, nodes, material, section } => {
                Ok(Self::QuadQ8(QuadQ8::new(id, nodes, material, section)?))
            }
        }
    }
}

impl TryFrom<ElementLoad2DInputKind> for Load2DInputConversion {
    type Error = FemError;

    fn try_from(value: ElementLoad2DInputKind) -> Result<Self, Self::Error> {
        match value {
            ElementLoad2DInputKind::BeamUniform { element, coordinate_system, qx, qy } => Ok(Self::Element(
                ElementLoad2D::BeamUniformLine(BeamUniformLineLoad2D::new(element, coordinate_system.into(), qx, qy)?),
            )),
            ElementLoad2DInputKind::EdgeTraction { element, edge, coordinate_system, tx, ty } => Ok(Self::Element(
                ElementLoad2D::EdgeTraction(EdgeTraction2D::new(element, edge, coordinate_system.into(), tx, ty)?),
            )),
            ElementLoad2DInputKind::BodyForce { element, bx, by } => {
                Ok(Self::Element(ElementLoad2D::BodyForce(BodyForce2D::new(element, bx, by)?)))
            }
            ElementLoad2DInputKind::SelfWeight { element, ax, ay } => {
                Ok(Self::Element(ElementLoad2D::SelfWeight(SelfWeight2D::new(element, ax, ay)?)))
            }
            ElementLoad2DInputKind::Nodal { node, dof, value } => {
                Ok(Self::Nodal(NodalLoad2D::new(node, dof.into(), value)?))
            }
        }
    }
}

impl From<LoadCoordinateSystem2DInput> for LoadCoordinateSystem2D {
    fn from(value: LoadCoordinateSystem2DInput) -> Self {
        match value {
            LoadCoordinateSystem2DInput::Global => Self::Global,
            LoadCoordinateSystem2DInput::Local => Self::Local,
        }
    }
}

enum Load2DInputConversion {
    Nodal(NodalLoad2D),
    Element(ElementLoad2D),
}

#[cfg(test)]
mod tests {
    use super::{ElementType2DInput, Model2DInput, SolverKind2DInput};
    use crate::elements::Element2D;
    use crate::model::{Dof2D, SolverKind2D};

    #[test]
    fn builds_model_from_minimal_json_input() {
        let input: Model2DInput = serde_json::from_str(
            r#"
            {
              "analysis_settings": { "solver": "sparse", "cg_tolerance": 1e-8 },
              "materials": [
                { "id": 1, "young_modulus": 200.0, "poisson_ratio": 0.3, "density": 1.0 }
              ],
              "sections": [
                { "id": 10, "type": "plane_stress", "thickness": 0.2 }
              ],
              "nodes": [
                { "id": 1, "x": 0.0, "y": 0.0 },
                { "id": 2, "x": 1.0, "y": 0.0 },
                { "id": 3, "x": 0.0, "y": 1.0 }
              ],
              "elements": [
                { "type": "triangle", "id": 100, "nodes": [1, 2, 3], "material": 1, "section": 10 }
              ],
              "constraints": [
                { "node": 1, "dof": "ux", "value": 0.0 },
                { "node": 1, "dof": "Uy", "value": 0.0 }
              ],
              "loads": [
                { "type": "edge_traction", "element": 100, "edge": [2, 3], "coordinate_system": "global", "tx": 0.0, "ty": -5.0 },
                { "type": "nodal", "node": 2, "dof": "uy", "value": -1.0 }
              ]
            }
            "#,
        )
        .expect("JSON input should deserialize");

        assert_eq!(input.analysis_settings.solver, Some(SolverKind2DInput::Sparse));

        let model = input.into_model().expect("model should be built");

        assert_eq!(model.analysis_settings().solver(), SolverKind2D::Sparse);
        assert_eq!(model.nodes().len(), 3);
        assert_eq!(model.elements().len(), 1);
        assert_eq!(model.constraints()[1].dof(), Dof2D::Uy);
        assert_eq!(model.loads().len(), 1);
        assert_eq!(model.element_loads().len(), 1);
    }

    #[test]
    fn supports_t3_alias_when_deserializing_elements() {
        let element: ElementType2DInput =
            serde_json::from_str(r#"{ "type": "t3", "id": 7, "nodes": [1, 2, 3], "material": 1, "section": 2 }"#)
                .expect("T3 alias should deserialize");

        let element: Element2D = element.try_into().expect("element should convert");

        assert!(matches!(element, Element2D::TriangleT3(_)));
    }

    #[test]
    fn snapshots_and_rebuilds_a_validated_model() {
        let input: Model2DInput = serde_json::from_str(include_str!("../../examples/t3_cantilever.json"))
            .expect("example input should deserialize");
        let model = input.into_model().expect("example model should build");

        let snapshot = Model2DInput::from_model(&model);
        let rebuilt = snapshot.clone().into_model().expect("snapshot should rebuild");
        let rebuilt_snapshot = Model2DInput::from_model(&rebuilt);

        assert_eq!(rebuilt.analysis_settings(), model.analysis_settings());
        assert_eq!(rebuilt_snapshot, snapshot);
    }
}
