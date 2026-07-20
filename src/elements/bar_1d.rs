use nalgebra::Matrix2;
use serde::Serialize;

use crate::FemError;
use crate::math::tolerance::{DEFAULT_GEOMETRY_TOLERANCE, is_near_zero};
use crate::model::{BarSection, LinearElasticMaterial, Node1D};

/// Two-node axial bar element in one spatial dimension.
#[derive(Debug, Clone, Copy, PartialEq, Serialize)]
pub struct BarElement1D {
    id: usize,
    node_ids: [usize; 2],
    material_id: usize,
    section_id: usize,
}

impl BarElement1D {
    /// Creates a two-node bar element.
    pub fn new(id: usize, node_ids: [usize; 2], material_id: usize, section_id: usize) -> Result<Self, FemError> {
        if node_ids[0] == node_ids[1] {
            return Err(FemError::InvalidElementConnectivity { element_id: id, node_ids });
        }

        Ok(Self { id, node_ids, material_id, section_id })
    }

    /// External element identifier.
    #[must_use]
    pub fn id(&self) -> usize {
        self.id
    }

    /// External identifiers of the element end nodes.
    #[must_use]
    pub fn node_ids(&self) -> [usize; 2] {
        self.node_ids
    }

    /// External material identifier.
    #[must_use]
    pub fn material_id(&self) -> usize {
        self.material_id
    }

    /// External section identifier.
    #[must_use]
    pub fn section_id(&self) -> usize {
        self.section_id
    }

    /// Local-to-global node id mapping for `[u1, u2]^T`.
    #[must_use]
    pub fn local_dof_node_ids(&self) -> [usize; 2] {
        self.node_ids
    }

    /// Computes the element length using the default geometry tolerance.
    pub fn length(&self, first: &Node1D, second: &Node1D) -> Result<f64, FemError> {
        self.length_with_tolerance(first, second, DEFAULT_GEOMETRY_TOLERANCE)
    }

    /// Computes the element length using an explicit absolute tolerance.
    pub fn length_with_tolerance(&self, first: &Node1D, second: &Node1D, tolerance: f64) -> Result<f64, FemError> {
        let received = [first.id(), second.id()];
        let reversed = [self.node_ids[1], self.node_ids[0]];
        if received != self.node_ids && received != reversed {
            return Err(FemError::ElementNodeMismatch { element_id: self.id, expected: self.node_ids, received });
        }

        let length = (second.x() - first.x()).abs();
        if is_near_zero(length, tolerance) {
            return Err(FemError::ZeroLengthElement { element_id: self.id, node_ids: self.node_ids, length });
        }

        Ok(length)
    }

    /// Computes the local 2 x 2 axial stiffness matrix.
    ///
    /// The local displacement order is `[u1, u2]^T`, where `u1` belongs to
    /// `self.node_ids()[0]` and `u2` belongs to `self.node_ids()[1]`.
    pub fn stiffness_matrix(
        &self, first: &Node1D, second: &Node1D, material: &LinearElasticMaterial, section: &BarSection,
    ) -> Result<Matrix2<f64>, FemError> {
        if material.id() != self.material_id {
            return Err(FemError::ElementMaterialMismatch {
                element_id: self.id,
                expected: self.material_id,
                received: material.id(),
            });
        }

        if section.id() != self.section_id {
            return Err(FemError::ElementSectionMismatch {
                element_id: self.id,
                expected: self.section_id,
                received: section.id(),
            });
        }

        let length = self.length(first, second)?;
        let axial_stiffness = material.young_modulus() * section.area() / length;

        Ok(Matrix2::new(axial_stiffness, -axial_stiffness, -axial_stiffness, axial_stiffness))
    }
}

#[cfg(test)]
mod tests {
    use approx::assert_relative_eq;

    use super::*;

    fn material(id: usize, young_modulus: f64) -> LinearElasticMaterial {
        LinearElasticMaterial::new(id, young_modulus, 0.3, 7850.0).unwrap()
    }

    fn section(id: usize, area: f64) -> BarSection {
        BarSection::new(id, area).unwrap()
    }

    fn assert_matrix_relative_eq(actual: &Matrix2<f64>, expected: &Matrix2<f64>) {
        for row in 0..2 {
            for col in 0..2 {
                assert_relative_eq!(actual[(row, col)], expected[(row, col)]);
            }
        }
    }

    #[test]
    fn creates_bar_element() {
        let element = BarElement1D::new(10, [1, 2], 3, 4).unwrap();

        assert_eq!(element.id(), 10);
        assert_eq!(element.node_ids(), [1, 2]);
        assert_eq!(element.material_id(), 3);
        assert_eq!(element.section_id(), 4);
        assert_eq!(element.local_dof_node_ids(), [1, 2]);
    }

    #[test]
    fn rejects_repeated_node_ids() {
        let error = BarElement1D::new(10, [1, 1], 3, 4).unwrap_err();

        assert_eq!(error, FemError::InvalidElementConnectivity { element_id: 10, node_ids: [1, 1] });
    }

    #[test]
    fn computes_length_for_matching_nodes() {
        let element = BarElement1D::new(10, [1, 2], 3, 4).unwrap();
        let first = Node1D::new(1, -0.25).unwrap();
        let second = Node1D::new(2, 0.75).unwrap();

        assert_relative_eq!(element.length(&first, &second).unwrap(), 1.0);
    }

    #[test]
    fn computes_length_for_reversed_nodes() {
        let element = BarElement1D::new(10, [1, 2], 3, 4).unwrap();
        let first = Node1D::new(2, 0.75).unwrap();
        let second = Node1D::new(1, -0.25).unwrap();

        assert_relative_eq!(element.length(&first, &second).unwrap(), 1.0);
    }

    #[test]
    fn rejects_mismatched_nodes() {
        let element = BarElement1D::new(10, [1, 2], 3, 4).unwrap();
        let first = Node1D::new(2, 0.0).unwrap();
        let second = Node1D::new(3, 1.0).unwrap();

        let error = element.length(&first, &second).unwrap_err();

        assert_eq!(error, FemError::ElementNodeMismatch { element_id: 10, expected: [1, 2], received: [2, 3] });
    }

    #[test]
    fn rejects_zero_length_geometry() {
        let element = BarElement1D::new(10, [1, 2], 3, 4).unwrap();
        let first = Node1D::new(1, 1.0).unwrap();
        let second = Node1D::new(2, 1.0).unwrap();

        let error = element.length(&first, &second).unwrap_err();

        assert_eq!(error, FemError::ZeroLengthElement { element_id: 10, node_ids: [1, 2], length: 0.0 });
    }

    #[test]
    fn computes_local_stiffness_matrix() {
        let element = BarElement1D::new(10, [1, 2], 3, 4).unwrap();
        let first = Node1D::new(1, 0.0).unwrap();
        let second = Node1D::new(2, 2.0).unwrap();
        let material = material(3, 200.0);
        let section = section(4, 0.5);

        let stiffness = element.stiffness_matrix(&first, &second, &material, &section).unwrap();

        let expected = Matrix2::new(50.0, -50.0, -50.0, 50.0);
        assert_matrix_relative_eq(&stiffness, &expected);
    }

    #[test]
    fn stiffness_matrix_is_independent_of_node_argument_order() {
        let element = BarElement1D::new(10, [1, 2], 3, 4).unwrap();
        let first = Node1D::new(2, 2.0).unwrap();
        let second = Node1D::new(1, 0.0).unwrap();
        let material = material(3, 200.0);
        let section = section(4, 0.5);

        let stiffness = element.stiffness_matrix(&first, &second, &material, &section).unwrap();

        let expected = Matrix2::new(50.0, -50.0, -50.0, 50.0);
        assert_matrix_relative_eq(&stiffness, &expected);
    }

    #[test]
    fn rejects_material_that_does_not_match_element_reference() {
        let element = BarElement1D::new(10, [1, 2], 3, 4).unwrap();
        let first = Node1D::new(1, 0.0).unwrap();
        let second = Node1D::new(2, 2.0).unwrap();
        let material = material(99, 200.0);
        let section = section(4, 0.5);

        let error = element.stiffness_matrix(&first, &second, &material, &section).unwrap_err();

        assert_eq!(error, FemError::ElementMaterialMismatch { element_id: 10, expected: 3, received: 99 });
    }

    #[test]
    fn rejects_section_that_does_not_match_element_reference() {
        let element = BarElement1D::new(10, [1, 2], 3, 4).unwrap();
        let first = Node1D::new(1, 0.0).unwrap();
        let second = Node1D::new(2, 2.0).unwrap();
        let material = material(3, 200.0);
        let section = section(99, 0.5);

        let error = element.stiffness_matrix(&first, &second, &material, &section).unwrap_err();

        assert_eq!(error, FemError::ElementSectionMismatch { element_id: 10, expected: 4, received: 99 });
    }

    #[test]
    fn stiffness_matrix_reuses_length_validation() {
        let element = BarElement1D::new(10, [1, 2], 3, 4).unwrap();
        let first = Node1D::new(1, 1.0).unwrap();
        let second = Node1D::new(2, 1.0).unwrap();
        let material = material(3, 200.0);
        let section = section(4, 0.5);

        let error = element.stiffness_matrix(&first, &second, &material, &section).unwrap_err();

        assert_eq!(error, FemError::ZeroLengthElement { element_id: 10, node_ids: [1, 2], length: 0.0 });
    }
}
