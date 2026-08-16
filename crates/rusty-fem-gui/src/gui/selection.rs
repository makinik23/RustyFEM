//! Shared GUI selection and view-layer state.

/// Entity selected in the model canvas.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum SelectedEntity {
    Node(usize),
    Element(usize),
    Edge { element_id: usize, node_ids: [usize; 2] },
}

impl SelectedEntity {
    pub(super) fn label(self) -> String {
        match self {
            Self::Node(node_id) => format!("Node {node_id}"),
            Self::Element(element_id) => format!("Element {element_id}"),
            Self::Edge { element_id, node_ids } => {
                format!("Edge {}-{} on element {element_id}", node_ids[0], node_ids[1])
            }
        }
    }
}

/// Canvas visibility toggles.
#[derive(Debug, Clone, Copy)]
pub(super) struct ViewOptions {
    pub(super) show_mesh: bool,
    pub(super) show_nodes: bool,
    pub(super) show_node_ids: bool,
    pub(super) show_element_ids: bool,
    pub(super) show_constraints: bool,
    pub(super) show_loads: bool,
    pub(super) show_boundary_edges: bool,
}

impl Default for ViewOptions {
    fn default() -> Self {
        Self {
            show_mesh: true,
            show_nodes: true,
            show_node_ids: false,
            show_element_ids: false,
            show_constraints: true,
            show_loads: true,
            show_boundary_edges: true,
        }
    }
}
