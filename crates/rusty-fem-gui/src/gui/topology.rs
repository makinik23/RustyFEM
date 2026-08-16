//! Mesh-topology helpers shared by GUI panels.

use rusty_fem::elements::Element2D;

/// Returns a normalized, direction-independent two-node edge key.
pub(super) fn normalized_edge(first: usize, second: usize) -> [usize; 2] {
    if first < second { [first, second] } else { [second, first] }
}

/// Returns the visible outline node sequence for an element.
///
/// Higher-order element outlines include midside nodes so curved or quadratic
/// edges are visible as segmented paths in the viewer.
pub(super) fn element_outline_node_ids(element: &Element2D) -> Vec<usize> {
    let nodes = element.node_ids();

    match element {
        Element2D::Truss(_) | Element2D::Beam(_) => nodes.to_vec(),
        Element2D::TriangleT3(_) => vec![nodes[0], nodes[1], nodes[2], nodes[0]],
        Element2D::TriangleT6(_) => vec![nodes[0], nodes[3], nodes[1], nodes[4], nodes[2], nodes[5], nodes[0]],
        Element2D::QuadQ4(_) => vec![nodes[0], nodes[1], nodes[2], nodes[3], nodes[0]],
        Element2D::QuadQ8(_) => {
            vec![nodes[0], nodes[4], nodes[1], nodes[5], nodes[2], nodes[6], nodes[3], nodes[7], nodes[0]]
        }
    }
}

/// Expands an element edge into the visible outline segments used by the GUI.
///
/// For T3 and Q4 this usually returns one segment. For T6 and Q8 it returns the
/// two corner-to-midside segments that together represent the loaded edge.
pub(super) fn element_edge_segments(element: &Element2D, edge_node_ids: [usize; 2]) -> Vec<[usize; 2]> {
    let outline = element_outline_node_ids(element);
    let Some(closed_node_id) = outline.last() else {
        return Vec::new();
    };
    let unique_outline_length = outline.len() - usize::from(outline.first() == Some(closed_node_id));
    let outline_nodes = &outline[..unique_outline_length];

    let Some(first_index) = outline_nodes.iter().position(|node_id| *node_id == edge_node_ids[0]) else {
        return Vec::new();
    };
    let Some(second_index) = outline_nodes.iter().position(|node_id| *node_id == edge_node_ids[1]) else {
        return Vec::new();
    };

    let forward = outline_path(outline_nodes, first_index, second_index);
    let backward = outline_path(outline_nodes, second_index, first_index);
    let path = if forward.len() <= backward.len() { forward } else { backward };

    path.windows(2).map(|edge| normalized_edge(edge[0], edge[1])).collect()
}

fn outline_path(outline_nodes: &[usize], start_index: usize, end_index: usize) -> Vec<usize> {
    let mut path = vec![outline_nodes[start_index]];
    let mut index = start_index;

    while index != end_index {
        index = (index + 1) % outline_nodes.len();
        path.push(outline_nodes[index]);

        if path.len() > outline_nodes.len() + 1 {
            return Vec::new();
        }
    }

    path
}

#[cfg(test)]
mod tests {
    use super::{element_edge_segments, element_outline_node_ids};
    use rusty_fem::elements::{Element2D, QuadQ4, QuadQ8, TriangleT6};

    #[test]
    fn q4_edge_has_one_visible_segment() {
        let element = Element2D::QuadQ4(QuadQ4::new(1, [1, 2, 3, 4], 1, 1).expect("valid element"));

        assert_eq!(element_edge_segments(&element, [1, 2]), vec![[1, 2]]);
        assert_eq!(element_edge_segments(&element, [2, 1]), vec![[1, 2]]);
    }

    #[test]
    fn t6_edge_expands_through_midside_node() {
        let element = Element2D::TriangleT6(TriangleT6::new(1, [1, 2, 3, 4, 5, 6], 1, 1).expect("valid element"));

        assert_eq!(element_edge_segments(&element, [1, 2]), vec![[1, 4], [2, 4]]);
        assert_eq!(element_edge_segments(&element, [2, 1]), vec![[1, 4], [2, 4]]);
    }

    #[test]
    fn q8_edge_expands_through_midside_node() {
        let element = Element2D::QuadQ8(QuadQ8::new(1, [1, 2, 3, 4, 5, 6, 7, 8], 1, 1).expect("valid element"));

        assert_eq!(element_edge_segments(&element, [4, 1]), vec![[4, 8], [1, 8]]);
    }

    #[test]
    fn q8_outline_includes_midside_nodes() {
        let element = Element2D::QuadQ8(QuadQ8::new(1, [1, 2, 3, 4, 5, 6, 7, 8], 1, 1).expect("valid element"));

        assert_eq!(element_outline_node_ids(&element), vec![1, 5, 2, 6, 3, 7, 4, 8, 1]);
    }
}
