use thiserror::Error;

/// Error type shared by public RustyFEM APIs.
#[derive(Debug, Error, PartialEq)]
pub enum FemError {
    /// A model collection contains the same external identifier more than once.
    #[error("duplicate {entity} id {id}")]
    DuplicateId { entity: &'static str, id: usize },

    /// A requested external identifier cannot be found in a model collection.
    #[error("{entity} with id {id} does not exist")]
    UnknownId { entity: &'static str, id: usize },

    /// A model object references another object that does not exist.
    #[error("{owner_entity} {owner_id} references missing {referenced_entity} {referenced_id}")]
    InvalidReference {
        owner_entity: &'static str,
        owner_id: usize,
        referenced_entity: &'static str,
        referenced_id: usize,
    },

    /// A node coordinate is not finite.
    #[error("node {node_id} has invalid x coordinate {x}; coordinates must be finite")]
    InvalidNodeCoordinate { node_id: usize, x: f64 },

    /// A material property violates the supported linear-elastic assumptions.
    #[error("material {material_id} has invalid {property}={value}: {reason}")]
    InvalidMaterialProperty { material_id: usize, property: &'static str, value: f64, reason: &'static str },

    /// A section property is outside the physically meaningful range.
    #[error("section {section_id} has invalid {property}={value}: {reason}")]
    InvalidSectionProperty { section_id: usize, property: &'static str, value: f64, reason: &'static str },

    /// A nodal load value cannot be represented in the model.
    #[error("nodal load at node {node_id} has invalid force {force}; force must be finite")]
    InvalidNodalLoad { node_id: usize, force: f64 },

    /// A prescribed displacement value cannot be represented in the model.
    #[error(
        "displacement constraint at node {node_id} has invalid displacement {displacement}; \
         displacement must be finite"
    )]
    InvalidDisplacementConstraint { node_id: usize, displacement: f64 },

    /// An element references the same node more than once.
    #[error("element {element_id} has invalid connectivity {node_ids:?}: node ids must differ")]
    InvalidElementConnectivity { element_id: usize, node_ids: [usize; 2] },

    /// The supplied nodes do not match an element's connectivity.
    #[error("element {element_id} expected node ids {expected:?}, but received {received:?}")]
    ElementNodeMismatch { element_id: usize, expected: [usize; 2], received: [usize; 2] },

    /// The supplied material does not match an element's material reference.
    #[error("element {element_id} expected material id {expected}, but received {received}")]
    ElementMaterialMismatch { element_id: usize, expected: usize, received: usize },

    /// The supplied section does not match an element's section reference.
    #[error("element {element_id} expected section id {expected}, but received {received}")]
    ElementSectionMismatch { element_id: usize, expected: usize, received: usize },

    /// An element has zero or numerically negligible length.
    #[error("element {element_id} has zero or near-zero length {length} between nodes {node_ids:?}")]
    ZeroLengthElement { element_id: usize, node_ids: [usize; 2], length: f64 },
}
