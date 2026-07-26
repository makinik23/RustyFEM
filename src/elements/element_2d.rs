//! Defines 2D elements used in finite element analysis.

use super::Interpolation;
use crate::error::FemError;
use crate::model::Dof2D;

const TRANSLATIONAL_DOFS: &[Dof2D] = &[Dof2D::Ux, Dof2D::Uy];

const FRAME_DOFS: &[Dof2D] = &[Dof2D::Ux, Dof2D::Uy, Dof2D::Rz];

/// Truss element in 2D space, defined by two nodes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Truss2D {
    id: usize,
    node_ids: [usize; 2],
}

/// Beam element in 2D space, defined by two nodes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Beam2D {
    id: usize,
    node_ids: [usize; 2],
}

/// Triangle element in 2D space, defined by three nodes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct TriangleT3 {
    id: usize,
    node_ids: [usize; 3],
}

/// Enum representing different types of 2D elements used in finite element analysis.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Element2D {
    Truss(Truss2D),
    Beam(Beam2D),
    TriangleT3(TriangleT3),
}

impl Truss2D {
    pub fn new(id: usize, node_ids: [usize; 2]) -> Result<Self, FemError> {
        if node_ids[0] == node_ids[1] {
            return Err(FemError::InvalidElementConnectivity { element_id: id, node_ids: node_ids.to_vec() });
        }

        Ok(Self { id, node_ids })
    }
}

impl Beam2D {
    pub fn new(id: usize, node_ids: [usize; 2]) -> Result<Self, FemError> {
        if node_ids[0] == node_ids[1] {
            return Err(FemError::InvalidElementConnectivity { element_id: id, node_ids: node_ids.to_vec() });
        }

        Ok(Self { id, node_ids })
    }
}

impl TriangleT3 {
    pub fn new(id: usize, node_ids: [usize; 3]) -> Result<Self, FemError> {
        if node_ids[0] == node_ids[1] || node_ids[0] == node_ids[2] || node_ids[1] == node_ids[2] {
            return Err(FemError::InvalidElementConnectivity { element_id: id, node_ids: node_ids.to_vec() });
        }

        Ok(Self { id, node_ids })
    }
}

impl Element2D {
    pub fn id(&self) -> usize {
        match self {
            Self::Truss(element) => element.id,
            Self::Beam(element) => element.id,
            Self::TriangleT3(element) => element.id,
        }
    }

    pub fn node_ids(&self) -> &[usize] {
        match self {
            Self::Truss(element) => &element.node_ids,
            Self::Beam(element) => &element.node_ids,
            Self::TriangleT3(element) => &element.node_ids,
        }
    }

    pub fn interpolation(&self) -> Interpolation {
        match self {
            Self::Truss(_) => Interpolation::LinearLagrange,
            Self::Beam(_) => Interpolation::CubicHermite,
            Self::TriangleT3(_) => Interpolation::LinearTriangleT3,
        }
    }

    pub fn dofs_per_node(&self) -> &'static [Dof2D] {
        match self {
            Self::Truss(_) => TRANSLATIONAL_DOFS,
            Self::Beam(_) => FRAME_DOFS,
            Self::TriangleT3(_) => TRANSLATIONAL_DOFS,
        }
    }

    pub fn dof_count(&self) -> usize {
        self.node_ids().len() * self.dofs_per_node().len()
    }
}

// TODO cases
#[cfg(test)]
mod tests {
    use super::{Beam2D, Element2D, TriangleT3, Truss2D};
    use crate::elements::interpolation::Interpolation;
    use crate::error::FemError;
    use crate::model::Dof2D;

    #[test]
    fn creates_valid_elements() {
        let cases = [
            ("truss", Truss2D::new(10, [1, 2]).map(Element2D::Truss), 10, vec![1, 2]),
            ("beam", Beam2D::new(20, [2, 3]).map(Element2D::Beam), 20, vec![2, 3]),
            ("triangle", TriangleT3::new(30, [3, 4, 5]).map(Element2D::TriangleT3), 30, vec![3, 4, 5]),
        ];

        for (name, result, expected_element_id, expected_node_ids) in cases {
            let element = result.expect("valid element should be created");

            assert_eq!(element.id(), expected_element_id, "failed case: {name}");

            assert_eq!(element.node_ids(), expected_node_ids.as_slice(), "failed case: {name}");
        }
    }

    #[test]
    fn rejects_elements_with_repeated_nodes() {
        let cases = [
            ("truss", Truss2D::new(10, [1, 1]).map(Element2D::Truss), 10, vec![1, 1]),
            ("beam", Beam2D::new(20, [2, 2]).map(Element2D::Beam), 20, vec![2, 2]),
            ("triangle", TriangleT3::new(30, [3, 4, 3]).map(Element2D::TriangleT3), 30, vec![3, 4, 3]),
        ];

        for (name, result, expected_element_id, expected_node_ids) in cases {
            assert!(
                matches!(
                    result,
                    Err(FemError::InvalidElementConnectivity
                    {
                        element_id,
                        node_ids,
                    })
                    if element_id == expected_element_id
                        && node_ids == expected_node_ids
                ),
                "failed case: {name}"
            );
        }
    }

    #[test]
    fn elements_use_correct_interpolation() {
        let cases = [
            (
                "truss",
                Element2D::Truss(Truss2D::new(10, [1, 2]).expect("valid truss should be created")),
                Interpolation::LinearLagrange,
            ),
            (
                "beam",
                Element2D::Beam(Beam2D::new(20, [2, 3]).expect("valid beam should be created")),
                Interpolation::CubicHermite,
            ),
            (
                "triangle",
                Element2D::TriangleT3(TriangleT3::new(30, [3, 4, 5]).expect("valid triangle should be created")),
                Interpolation::LinearTriangleT3,
            ),
        ];

        for (name, element, expected_interpolation) in cases {
            assert_eq!(element.interpolation(), expected_interpolation, "failed case: {name}");
        }
    }

    #[test]
    fn elements_use_correct_degrees_of_freedom() {
        let cases = [
            (
                "truss",
                Element2D::Truss(Truss2D::new(10, [1, 2]).expect("valid truss should be created")),
                vec![Dof2D::Ux, Dof2D::Uy],
                4,
            ),
            (
                "beam",
                Element2D::Beam(Beam2D::new(20, [2, 3]).expect("valid beam should be created")),
                vec![Dof2D::Ux, Dof2D::Uy, Dof2D::Rz],
                6,
            ),
            (
                "triangle",
                Element2D::TriangleT3(TriangleT3::new(30, [3, 4, 5]).expect("valid triangle should be created")),
                vec![Dof2D::Ux, Dof2D::Uy],
                6,
            ),
        ];

        for (name, element, expected_dofs, expected_count) in cases {
            assert_eq!(element.dofs_per_node(), expected_dofs.as_slice(), "failed case: {name}");

            assert_eq!(element.dof_count(), expected_count, "failed case: {name}");
        }
    }
}
