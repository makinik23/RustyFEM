//! Types of FEM errors.

#[derive(Debug, thiserror::Error)]
pub enum FemError {
    #[error("node {node_id} has invalid {axis} coordinate: {value}")]
    InvalidNodeCoordinate { node_id: usize, axis: &'static str, value: f64 },

    #[error("node {node_id} has invalid prescribed displacement: {displacement}")]
    InvalidDisplacementConstraint { node_id: usize, displacement: f64 },

    #[error("duplicate {entity} ID: {id}")]
    DuplicateId { entity: &'static str, id: usize },

    #[error("unknown {entity} ID: {id}")]
    UnknownId { entity: &'static str, id: usize },

    #[error("element {element_id} has invalid connectivity: {node_ids:?}")]
    InvalidElementConnectivity { element_id: usize, node_ids: Vec<usize> },

    #[error("element {element_id} ({element_type}) has invalid {property}: {value}; {reason}")]
    InvalidElementProperty {
        element_id: usize,
        element_type: &'static str,
        property: &'static str,
        value: f64,
        reason: &'static str,
    },

    #[error("element {element_id} ({element_type}) is degenerate: {measure_name} = {measure}; node IDs: {node_ids:?}")]
    DegenerateElement {
        element_id: usize,
        element_type: &'static str,
        node_ids: Vec<usize>,
        measure_name: &'static str,
        measure: f64,
    },

    #[error("invalid analysis space: {value}")]
    InvalidAnalysisSpace { value: String },

    #[error("node {node_id} does not have degree of freedom {dof}")]
    UnknownDof { node_id: usize, dof: &'static str },

    #[error("material has invalid {property}: {value}; {reason}")]
    InvalidMaterialProperty { property: &'static str, value: f64, reason: &'static str },

    #[error("model has no material")]
    MissingMaterial,

    #[error("nodal load on node {node_id} has invalid value: {value}")]
    InvalidNodalLoad { node_id: usize, value: f64 },

    #[error(
        "linear system has incompatible dimensions: stiffness matrix is {stiffness_rows}x{stiffness_columns}, load vector has length {load_vector_length}"
    )]
    IncompatibleLinearSystem { stiffness_rows: usize, stiffness_columns: usize, load_vector_length: usize },

    #[error("linear system references invalid DOF index: {index}")]
    InvalidDofIndex { index: usize },

    #[error("the stiffness matrix is singular")]
    SingularSystem,

    #[error("displacement vector has invalid length: expected {expected}, got {actual}")]
    InvalidDisplacementVector { expected: usize, actual: usize },

    #[error("interpolation length is invalid: {value}; must be finite and strictly positive")]
    InvalidInterpolationLength { value: f64 },

    #[error("interpolation coordinate is invalid: {coordinate}; expected a value in [0, {length}]")]
    InvalidInterpolationCoordinate { coordinate: f64, length: f64 },

    #[error(
        "triangle natural coordinates are invalid: xi = {xi}, eta = {eta}; expected xi >= 0, eta >= 0, xi + eta <= 1"
    )]
    InvalidTriangleNaturalCoordinates { xi: f64, eta: f64 },

    #[error("invalid interpolation position for {element_type}: expected {expected}")]
    InvalidElementInterpolationPosition { element_type: &'static str, expected: &'static str },
}
