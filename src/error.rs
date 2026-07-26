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
}
