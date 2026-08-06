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

    #[error("section ({section_type}) has invalid {property}: {value}; {reason}")]
    InvalidSectionProperty { section_type: &'static str, property: &'static str, value: f64, reason: &'static str },

    #[error("section {section_id} has type {actual}, expected {expected}")]
    InvalidSectionType { section_id: usize, expected: &'static str, actual: &'static str },

    #[error("nodal load on node {node_id} has invalid value: {value}")]
    InvalidNodalLoad { node_id: usize, value: f64 },

    #[error("element load {load_type} on element {element_id} has invalid {component}: {value}")]
    InvalidElementLoadValue { element_id: usize, load_type: &'static str, component: &'static str, value: f64 },

    #[error("element load {load_type} cannot be applied to {actual} element {element_id}; expected {expected}")]
    InvalidElementLoadType { element_id: usize, load_type: &'static str, expected: &'static str, actual: &'static str },

    #[error("element load {load_type} on element {element_id} has invalid edge {node_ids:?}; expected {expected}")]
    InvalidElementLoadEdge { element_id: usize, load_type: &'static str, node_ids: Vec<usize>, expected: &'static str },

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

    #[error("vector {vector} has invalid length: expected {expected}, got {actual}")]
    InvalidVectorLength { vector: &'static str, expected: usize, actual: usize },

    #[error("conjugate gradient tolerance must be finite and positive: {value}")]
    InvalidSolverTolerance { value: f64 },

    #[error("conjugate gradient iteration limit must be positive: {value}")]
    InvalidSolverIterationLimit { value: usize },

    #[error("stagnation tolerance must be finite and non-negative: {value}")]
    InvalidStagnationTolerance { value: f64 },

    #[error("conjugate gradient method broke down")]
    ConjugateGradientBreakdown,

    #[error("invalid Jacobi preconditioner diagonal at index {index}: {value}")]
    InvalidPreconditionerDiagonal { index: usize, value: f64 },

    #[error("iterative solver did not converge after {iterations} iterations; residual norm = {residual_norm}")]
    IterativeSolverDidNotConverge { iterations: usize, residual_norm: f64 },

    #[error("iterative solver stagnated after {iterations} iterations; residual norm = {residual_norm}")]
    IterativeSolverStagnated { iterations: usize, residual_norm: f64 },
}
