//! Defines 2D elements used in finite element analysis.

use super::Interpolation;
use crate::error::FemError;
use crate::model::{Dof2D, Material2D, Node2D};
use nalgebra::DMatrix;

const TRANSLATIONAL_DOFS: &[Dof2D] = &[Dof2D::Ux, Dof2D::Uy];

const FRAME_DOFS: &[Dof2D] = &[Dof2D::Ux, Dof2D::Uy, Dof2D::Rz];

/// Truss element in 2D space, defined by two nodes.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Truss2D {
    id: usize,
    node_ids: [usize; 2],
    cross_section_area: f64,
}

/// Beam element in 2D space, defined by two nodes.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Beam2D {
    id: usize,
    node_ids: [usize; 2],
    cross_section_area: f64,
    second_moment_of_area: f64,
}

/// Triangle element in 2D space, defined by three nodes.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TriangleT3 {
    id: usize,
    node_ids: [usize; 3],
    thickness: f64,
}

/// Enum representing different types of 2D elements used in finite element analysis.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Element2D {
    Truss(Truss2D),
    Beam(Beam2D),
    TriangleT3(TriangleT3),
}

impl Truss2D {
    /// Creates a new Truss2D element with the specified ID, node IDs, and cross-sectional area.
    pub fn new(id: usize, node_ids: [usize; 2], cross_section_area: f64) -> Result<Self, FemError> {
        if node_ids[0] == node_ids[1] {
            return Err(FemError::InvalidElementConnectivity { element_id: id, node_ids: node_ids.to_vec() });
        }

        if !cross_section_area.is_finite() || cross_section_area <= 0.0 {
            return Err(FemError::InvalidElementProperty {
                element_id: id,
                element_type: "truss",
                property: "cross-sectional area",
                value: cross_section_area,
                reason: "must be finite and strictly positive",
            });
        }

        Ok(Self { id, node_ids, cross_section_area })
    }

    /// Calculates the stiffness matrix for the truss element based on the provided material properties and node coordinates.
    pub fn stiffness_matrix(
        &self, material: &Material2D, first_node: &Node2D, second_node: &Node2D,
    ) -> Result<[[f64; 4]; 4], FemError> {
        let dx = second_node.x() - first_node.x();
        let dy = second_node.y() - first_node.y();
        let length = (dx * dx + dy * dy).sqrt();

        if length == 0.0 {
            return Err(FemError::DegenerateElement {
                element_id: self.id,
                element_type: "truss",
                node_ids: self.node_ids.to_vec(),
                measure_name: "length",
                measure: length,
            });
        }

        let c = dx / length;
        let s = dy / length;

        let factor = material.young_modulus() * self.cross_section_area() / length;

        let matrix = [
            [factor * c * c, factor * c * s, -factor * c * c, -factor * c * s],
            [factor * c * s, factor * s * s, -factor * c * s, -factor * s * s],
            [-factor * c * c, -factor * c * s, factor * c * c, factor * c * s],
            [-factor * c * s, -factor * s * s, factor * c * s, factor * s * s],
        ];

        Ok(matrix)
    }

    /// Returns the cross-sectional area of the truss.
    #[must_use]
    pub fn cross_section_area(&self) -> f64 {
        self.cross_section_area
    }
}

impl Beam2D {
    /// Creates a new Beam2D element with the specified ID, node IDs, cross-sectional area, and second moment of area.
    pub fn new(
        id: usize, node_ids: [usize; 2], cross_section_area: f64, second_moment_of_area: f64,
    ) -> Result<Self, FemError> {
        if node_ids[0] == node_ids[1] {
            return Err(FemError::InvalidElementConnectivity { element_id: id, node_ids: node_ids.to_vec() });
        }

        if !cross_section_area.is_finite() || cross_section_area <= 0.0 {
            return Err(FemError::InvalidElementProperty {
                element_id: id,
                element_type: "beam",
                property: "cross-sectional area",
                value: cross_section_area,
                reason: "must be finite and strictly positive",
            });
        }

        if !second_moment_of_area.is_finite() || second_moment_of_area <= 0.0 {
            return Err(FemError::InvalidElementProperty {
                element_id: id,
                element_type: "beam",
                property: "second moment of area",
                value: second_moment_of_area,
                reason: "must be finite and strictly positive",
            });
        }

        Ok(Self { id, node_ids, cross_section_area, second_moment_of_area })
    }

    /// Returns the cross-sectional area of the beam.
    #[must_use]
    pub fn cross_section_area(&self) -> f64 {
        self.cross_section_area
    }

    /// Returns the second moment of area of the beam.
    #[must_use]
    pub fn second_moment_of_area(&self) -> f64 {
        self.second_moment_of_area
    }

    /// Calculates the stiffness matrix of the beam element in global coordinates.
    pub fn stiffness_matrix(
        &self, material: &Material2D, first_node: &Node2D, second_node: &Node2D,
    ) -> Result<[[f64; 6]; 6], FemError> {
        let dx = second_node.x() - first_node.x();
        let dy = second_node.y() - first_node.y();
        let length = (dx * dx + dy * dy).sqrt();

        if length == 0.0 {
            return Err(FemError::DegenerateElement {
                element_id: self.id,
                element_type: "beam",
                node_ids: self.node_ids.to_vec(),
                measure_name: "length",
                measure: length,
            });
        }

        let c = dx / length;
        let s = dy / length;

        let e = material.young_modulus();
        let a = self.cross_section_area();
        let i = self.second_moment_of_area();

        let ea_over_l = e * a / length;
        let twelve_ei_over_l3 = 12.0 * e * i / length.powi(3);
        let six_ei_over_l2 = 6.0 * e * i / length.powi(2);
        let four_ei_over_l = 4.0 * e * i / length;
        let two_ei_over_l = 2.0 * e * i / length;

        let local_matrix = [
            [ea_over_l, 0.0, 0.0, -ea_over_l, 0.0, 0.0],
            [0.0, twelve_ei_over_l3, six_ei_over_l2, 0.0, -twelve_ei_over_l3, six_ei_over_l2],
            [0.0, six_ei_over_l2, four_ei_over_l, 0.0, -six_ei_over_l2, two_ei_over_l],
            [-ea_over_l, 0.0, 0.0, ea_over_l, 0.0, 0.0],
            [0.0, -twelve_ei_over_l3, -six_ei_over_l2, 0.0, twelve_ei_over_l3, -six_ei_over_l2],
            [0.0, six_ei_over_l2, two_ei_over_l, 0.0, -six_ei_over_l2, four_ei_over_l],
        ];

        let transformation = [
            [c, s, 0.0, 0.0, 0.0, 0.0],
            [-s, c, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, c, s, 0.0],
            [0.0, 0.0, 0.0, -s, c, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 1.0],
        ];

        let transformed_local_matrix = multiply_6x6(&transpose_6x6(&transformation), &local_matrix);

        Ok(multiply_6x6(&transformed_local_matrix, &transformation))
    }
}

/// Transposes a 6x6 matrix.
fn transpose_6x6(matrix: &[[f64; 6]; 6]) -> [[f64; 6]; 6] {
    let mut transposed = [[0.0; 6]; 6];

    for (row, matrix_row) in matrix.iter().enumerate() {
        for (column, value) in matrix_row.iter().enumerate() {
            transposed[column][row] = *value;
        }
    }

    transposed
}

/// Multiplies two 6x6 matrices and returns the resulting matrix.
fn multiply_6x6(left: &[[f64; 6]; 6], right: &[[f64; 6]; 6]) -> [[f64; 6]; 6] {
    let mut product = [[0.0; 6]; 6];

    for (row, product_row) in product.iter_mut().enumerate() {
        for (column, product_value) in product_row.iter_mut().enumerate() {
            *product_value =
                left[row].iter().enumerate().map(|(index, left_value)| left_value * right[index][column]).sum();
        }
    }

    product
}

impl TriangleT3 {
    /// Creates a three-node triangular element for a plane-stress analysis.
    pub fn new(id: usize, node_ids: [usize; 3], thickness: f64) -> Result<Self, FemError> {
        if node_ids[0] == node_ids[1] || node_ids[0] == node_ids[2] || node_ids[1] == node_ids[2] {
            return Err(FemError::InvalidElementConnectivity { element_id: id, node_ids: node_ids.to_vec() });
        }

        if !thickness.is_finite() || thickness <= 0.0 {
            return Err(FemError::InvalidElementProperty {
                element_id: id,
                element_type: "triangle_t3",
                property: "thickness",
                value: thickness,
                reason: "must be finite and strictly positive",
            });
        }

        Ok(Self { id, node_ids, thickness })
    }

    /// Returns the thickness of the triangular element.
    #[must_use]
    pub fn thickness(&self) -> f64 {
        self.thickness
    }

    /// Calculates the stiffness matrix using the linear plane-stress T3 formulation.
    pub fn stiffness_matrix(
        &self, material: &Material2D, first_node: &Node2D, second_node: &Node2D, third_node: &Node2D,
    ) -> Result<[[f64; 6]; 6], FemError> {
        let x1 = first_node.x();
        let y1 = first_node.y();
        let x2 = second_node.x();
        let y2 = second_node.y();
        let x3 = third_node.x();
        let y3 = third_node.y();

        let twice_signed_area = (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1);
        let area = 0.5 * twice_signed_area.abs();

        if area == 0.0 {
            return Err(FemError::DegenerateElement {
                element_id: self.id,
                element_type: "triangle_t3",
                node_ids: self.node_ids.to_vec(),
                measure_name: "area",
                measure: area,
            });
        }

        let b1 = y2 - y3;
        let b2 = y3 - y1;
        let b3 = y1 - y2;
        let c1 = x3 - x2;
        let c2 = x1 - x3;
        let c3 = x2 - x1;
        let inverse_twice_area = 1.0 / twice_signed_area;

        let strain_displacement_matrix = [
            [b1 * inverse_twice_area, 0.0, b2 * inverse_twice_area, 0.0, b3 * inverse_twice_area, 0.0],
            [0.0, c1 * inverse_twice_area, 0.0, c2 * inverse_twice_area, 0.0, c3 * inverse_twice_area],
            [
                c1 * inverse_twice_area,
                b1 * inverse_twice_area,
                c2 * inverse_twice_area,
                b2 * inverse_twice_area,
                c3 * inverse_twice_area,
                b3 * inverse_twice_area,
            ],
        ];

        let young_modulus = material.young_modulus();
        let poisson_ratio = material.poisson_ratio();
        let constitutive_factor = young_modulus / (1.0 - poisson_ratio * poisson_ratio);
        let constitutive_matrix = [
            [constitutive_factor, constitutive_factor * poisson_ratio, 0.0],
            [constitutive_factor * poisson_ratio, constitutive_factor, 0.0],
            [0.0, 0.0, constitutive_factor * (1.0 - poisson_ratio) / 2.0],
        ];

        let constitutive_times_strain = multiply_3x3_by_3x6(&constitutive_matrix, &strain_displacement_matrix);
        let mut stiffness_matrix =
            multiply_transpose_3x6_by_3x6(&strain_displacement_matrix, &constitutive_times_strain);
        let scale = self.thickness * area;

        for row in &mut stiffness_matrix {
            for value in row {
                *value *= scale;
            }
        }

        Ok(stiffness_matrix)
    }
}

/// Multiplies a 3x3 matrix by a 3x6 matrix.
fn multiply_3x3_by_3x6(left: &[[f64; 3]; 3], right: &[[f64; 6]; 3]) -> [[f64; 6]; 3] {
    let mut product = [[0.0; 6]; 3];

    for (row, product_row) in product.iter_mut().enumerate() {
        for (column, product_value) in product_row.iter_mut().enumerate() {
            *product_value =
                left[row].iter().enumerate().map(|(index, left_value)| left_value * right[index][column]).sum();
        }
    }

    product
}

/// Multiplies the transpose of a 3x6 matrix by a 3x6 matrix.
fn multiply_transpose_3x6_by_3x6(left: &[[f64; 6]; 3], right: &[[f64; 6]; 3]) -> [[f64; 6]; 6] {
    let mut product = [[0.0; 6]; 6];

    for (row, product_row) in product.iter_mut().enumerate() {
        for (column, product_value) in product_row.iter_mut().enumerate() {
            *product_value =
                left.iter().enumerate().map(|(index, left_row)| left_row[row] * right[index][column]).sum();
        }
    }

    product
}

impl Element2D {
    /// Returns the ID of the element.
    pub fn id(&self) -> usize {
        match self {
            Self::Truss(element) => element.id,
            Self::Beam(element) => element.id,
            Self::TriangleT3(element) => element.id,
        }
    }

    /// Returns the node IDs associated with the element.
    pub fn node_ids(&self) -> &[usize] {
        match self {
            Self::Truss(element) => &element.node_ids,
            Self::Beam(element) => &element.node_ids,
            Self::TriangleT3(element) => &element.node_ids,
        }
    }

    /// Returns the interpolation type used for the element.
    pub fn interpolation(&self) -> Interpolation {
        match self {
            Self::Truss(_) => Interpolation::LinearLagrange,
            Self::Beam(_) => Interpolation::CubicHermite,
            Self::TriangleT3(_) => Interpolation::LinearTriangleT3,
        }
    }

    /// Returns the degrees of freedom (DOFs) associated with each node of the element.
    pub fn dofs_per_node(&self) -> &'static [Dof2D] {
        match self {
            Self::Truss(_) => TRANSLATIONAL_DOFS,
            Self::Beam(_) => FRAME_DOFS,
            Self::TriangleT3(_) => TRANSLATIONAL_DOFS,
        }
    }

    /// Returns the total number of degrees of freedom (DOFs) for the element.
    pub fn dof_count(&self) -> usize {
        self.node_ids().len() * self.dofs_per_node().len()
    }

    /// Calculates the stiffness matrix of this element in global coordinates.
    ///
    /// The nodes slice is used to resolve the element's node IDs. The returned
    /// dynamic matrix has the same row and column order as the element DOFs.
    pub fn stiffness_matrix(&self, material: &Material2D, nodes: &[Node2D]) -> Result<DMatrix<f64>, FemError> {
        let element_nodes = self
            .node_ids()
            .iter()
            .map(|&node_id| {
                nodes
                    .iter()
                    .find(|node| node.id() == node_id)
                    .ok_or(FemError::UnknownId { entity: "node", id: node_id })
            })
            .collect::<Result<Vec<_>, _>>()?;

        match self {
            Self::Truss(element) => {
                let matrix = element.stiffness_matrix(material, element_nodes[0], element_nodes[1])?;

                Ok(dynamic_matrix_from_array(matrix))
            }
            Self::Beam(element) => {
                let matrix = element.stiffness_matrix(material, element_nodes[0], element_nodes[1])?;

                Ok(dynamic_matrix_from_array(matrix))
            }
            Self::TriangleT3(element) => {
                let matrix =
                    element.stiffness_matrix(material, element_nodes[0], element_nodes[1], element_nodes[2])?;

                Ok(dynamic_matrix_from_array(matrix))
            }
        }
    }

    /// Returns the cross-sectional area if the element is a truss. TODO
    #[must_use]
    pub fn cross_section_area(&self) -> Option<f64> {
        match self {
            Self::Truss(element) => Some(element.cross_section_area()),
            Self::Beam(element) => Some(element.cross_section_area()),
            _ => None,
        }
    }

    /// Returns the second moment of area if the element is a beam. TODO
    #[must_use]
    pub fn second_moment_of_area(&self) -> Option<f64> {
        match self {
            Self::Beam(element) => Some(element.second_moment_of_area()),
            _ => None,
        }
    }
}

/// Converts a square stack-allocated matrix into nalgebra's dynamic matrix.
fn dynamic_matrix_from_array<const N: usize>(matrix: [[f64; N]; N]) -> DMatrix<f64> {
    let values = matrix.into_iter().flatten().collect::<Vec<_>>();

    DMatrix::from_row_slice(N, N, &values)
}

// TODO cases
#[cfg(test)]
mod tests {
    use super::{Beam2D, Element2D, TriangleT3, Truss2D, dynamic_matrix_from_array};
    use crate::elements::interpolation::Interpolation;
    use crate::error::FemError;
    use crate::model::{Dof2D, Material2D, Node2D};

    #[test]
    fn creates_valid_elements() {
        let cases = [
            ("truss", Truss2D::new(10, [1, 2], 1.0).map(Element2D::Truss), 10, vec![1, 2]),
            ("beam", Beam2D::new(20, [2, 3], 1.0, 1.0).map(Element2D::Beam), 20, vec![2, 3]),
            ("triangle", TriangleT3::new(30, [3, 4, 5], 1.0).map(Element2D::TriangleT3), 30, vec![3, 4, 5]),
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
            ("truss", Truss2D::new(10, [1, 1], 1.0).map(Element2D::Truss), 10, vec![1, 1]),
            ("beam", Beam2D::new(20, [2, 2], 1.0, 1.0).map(Element2D::Beam), 20, vec![2, 2]),
            ("triangle", TriangleT3::new(30, [3, 4, 3], 1.0).map(Element2D::TriangleT3), 30, vec![3, 4, 3]),
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
                Element2D::Truss(Truss2D::new(10, [1, 2], 1.0).expect("valid truss should be created")),
                Interpolation::LinearLagrange,
            ),
            (
                "beam",
                Element2D::Beam(Beam2D::new(20, [2, 3], 1.0, 1.0).expect("valid beam should be created")),
                Interpolation::CubicHermite,
            ),
            (
                "triangle",
                Element2D::TriangleT3(TriangleT3::new(30, [3, 4, 5], 1.0).expect("valid triangle should be created")),
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
                Element2D::Truss(Truss2D::new(10, [1, 2], 1.0).expect("valid truss should be created")),
                vec![Dof2D::Ux, Dof2D::Uy],
                4,
            ),
            (
                "beam",
                Element2D::Beam(Beam2D::new(20, [2, 3], 1.0, 1.0).expect("valid beam should be created")),
                vec![Dof2D::Ux, Dof2D::Uy, Dof2D::Rz],
                6,
            ),
            (
                "triangle",
                Element2D::TriangleT3(TriangleT3::new(30, [3, 4, 5], 1.0).expect("valid triangle should be created")),
                vec![Dof2D::Ux, Dof2D::Uy],
                6,
            ),
        ];

        for (name, element, expected_dofs, expected_count) in cases {
            assert_eq!(element.dofs_per_node(), expected_dofs.as_slice(), "failed case: {name}");

            assert_eq!(element.dof_count(), expected_count, "failed case: {name}");
        }
    }

    #[test]
    fn common_stiffness_matrix_dispatches_for_every_element_type() {
        let material = Material2D::new(200.0, 0.3, 7800.0).expect("valid material should be created");
        let nodes = vec![
            Node2D::new(10, 0.0, 0.0).expect("valid node should be created"),
            Node2D::new(20, 1.0, 0.0).expect("valid node should be created"),
            Node2D::new(30, 0.0, 1.0).expect("valid node should be created"),
        ];
        let cases = [
            ("truss", Element2D::Truss(Truss2D::new(1, [10, 20], 1.0).expect("valid truss should be created"))),
            ("beam", Element2D::Beam(Beam2D::new(2, [10, 20], 1.0, 1.0).expect("valid beam should be created"))),
            (
                "triangle",
                Element2D::TriangleT3(TriangleT3::new(3, [10, 20, 30], 1.0).expect("valid triangle should be created")),
            ),
        ];

        for (name, element) in cases {
            let actual = element.stiffness_matrix(&material, &nodes).expect("stiffness matrix should be calculated");
            let expected = match &element {
                Element2D::Truss(truss) => dynamic_matrix_from_array(
                    truss
                        .stiffness_matrix(&material, &nodes[0], &nodes[1])
                        .expect("stiffness matrix should be calculated"),
                ),
                Element2D::Beam(beam) => dynamic_matrix_from_array(
                    beam.stiffness_matrix(&material, &nodes[0], &nodes[1])
                        .expect("stiffness matrix should be calculated"),
                ),
                Element2D::TriangleT3(triangle) => dynamic_matrix_from_array(
                    triangle
                        .stiffness_matrix(&material, &nodes[0], &nodes[1], &nodes[2])
                        .expect("stiffness matrix should be calculated"),
                ),
            };

            assert_eq!(actual, expected, "failed case: {name}");
        }
    }

    #[test]
    fn common_stiffness_matrix_rejects_unknown_node() {
        let material = Material2D::new(200.0, 0.3, 7800.0).expect("valid material should be created");
        let known_node = Node2D::new(10, 0.0, 0.0).expect("valid node should be created");
        let element =
            Element2D::Truss(Truss2D::new(1, [10, 20], 1.0).expect("valid truss connectivity should be created"));

        let result = element.stiffness_matrix(&material, &[known_node]);

        assert!(matches!(result, Err(FemError::UnknownId { entity: "node", id: 20 })));
    }

    #[test]
    fn rejects_invalid_truss_cross_section_area() {
        let cases = [("zero", 0.0), ("negative", -1.0), ("infinite", f64::INFINITY), ("not a number", f64::NAN)];

        for (name, area) in cases {
            let result = Truss2D::new(10, [1, 2], area);

            assert!(
                matches!(
                    result,
                    Err(FemError::InvalidElementProperty {
                        element_id: 10,
                        element_type: "truss",
                        property: "cross-sectional area",
                        value,
                        reason: "must be finite and strictly positive",
                    }) if value == area || (value.is_nan() && area.is_nan())
                ),
                "failed case: {name}"
            );
        }
    }

    #[test]
    fn rejects_invalid_beam_cross_section_area() {
        let cases = [("zero", 0.0), ("negative", -1.0), ("infinite", f64::INFINITY), ("not a number", f64::NAN)];

        for (name, area) in cases {
            let result = Beam2D::new(20, [2, 3], area, 1.0);

            assert!(
                matches!(
                    result,
                    Err(FemError::InvalidElementProperty {
                        element_id: 20,
                        element_type: "beam",
                        property: "cross-sectional area",
                        value,
                        reason: "must be finite and strictly positive",
                    }) if value == area || (value.is_nan() && area.is_nan())
                ),
                "failed case: {name}"
            );
        }
    }

    #[test]
    fn rejects_invalid_beam_second_moment_of_area() {
        let cases = [("zero", 0.0), ("negative", -1.0), ("infinite", f64::INFINITY), ("not a number", f64::NAN)];

        for (name, second_moment_of_area) in cases {
            let result = Beam2D::new(20, [2, 3], 1.0, second_moment_of_area);

            assert!(
                matches!(
                    result,
                    Err(FemError::InvalidElementProperty {
                        element_id: 20,
                        element_type: "beam",
                        property: "second moment of area",
                        value,
                        reason: "must be finite and strictly positive",
                    }) if value == second_moment_of_area || (value.is_nan() && second_moment_of_area.is_nan())
                ),
                "failed case: {name}"
            );
        }
    }

    #[test]
    fn rejects_invalid_triangle_thickness() {
        let cases = [("zero", 0.0), ("negative", -1.0), ("infinite", f64::INFINITY), ("not a number", f64::NAN)];

        for (name, thickness) in cases {
            let result = TriangleT3::new(30, [3, 4, 5], thickness);

            assert!(
                matches!(
                    result,
                    Err(FemError::InvalidElementProperty {
                        element_id: 30,
                        element_type: "triangle_t3",
                        property: "thickness",
                        value,
                        reason: "must be finite and strictly positive",
                    }) if value == thickness || (value.is_nan() && thickness.is_nan())
                ),
                "failed case: {name}"
            );
        }
    }

    #[test]
    fn calculates_plane_stress_stiffness_matrix_for_right_triangle() {
        let material = Material2D::new(1.0, 0.3, 1.0).expect("valid material should be created");
        let first_node = Node2D::new(1, 0.0, 0.0).expect("valid node should be created");
        let second_node = Node2D::new(2, 1.0, 0.0).expect("valid node should be created");
        let third_node = Node2D::new(3, 0.0, 1.0).expect("valid node should be created");
        let triangle = TriangleT3::new(30, [1, 2, 3], 1.0).expect("valid triangle should be created");

        let matrix = triangle
            .stiffness_matrix(&material, &first_node, &second_node, &third_node)
            .expect("stiffness matrix should be calculated");

        let expected = [
            [135.0 / 182.0, 65.0 / 182.0, -100.0 / 182.0, -35.0 / 182.0, -35.0 / 182.0, -30.0 / 182.0],
            [65.0 / 182.0, 135.0 / 182.0, -30.0 / 182.0, -35.0 / 182.0, -35.0 / 182.0, -100.0 / 182.0],
            [-100.0 / 182.0, -30.0 / 182.0, 100.0 / 182.0, 0.0, 0.0, 30.0 / 182.0],
            [-35.0 / 182.0, -35.0 / 182.0, 0.0, 35.0 / 182.0, 35.0 / 182.0, 0.0],
            [-35.0 / 182.0, -35.0 / 182.0, 0.0, 35.0 / 182.0, 35.0 / 182.0, 0.0],
            [-30.0 / 182.0, -100.0 / 182.0, 30.0 / 182.0, 0.0, 0.0, 100.0 / 182.0],
        ];

        assert_matrix_approximately_equal_6(&matrix, &expected);
    }

    #[test]
    fn triangle_stiffness_matrix_is_symmetric() {
        let material = Material2D::new(200.0, 0.3, 7800.0).expect("valid material should be created");
        let first_node = Node2D::new(1, 0.0, 0.0).expect("valid node should be created");
        let second_node = Node2D::new(2, 2.0, 0.0).expect("valid node should be created");
        let third_node = Node2D::new(3, 0.0, 1.0).expect("valid node should be created");
        let triangle = TriangleT3::new(30, [1, 2, 3], 0.2).expect("valid triangle should be created");

        let matrix = triangle
            .stiffness_matrix(&material, &first_node, &second_node, &third_node)
            .expect("stiffness matrix should be calculated");

        let transposed = [
            [matrix[0][0], matrix[1][0], matrix[2][0], matrix[3][0], matrix[4][0], matrix[5][0]],
            [matrix[0][1], matrix[1][1], matrix[2][1], matrix[3][1], matrix[4][1], matrix[5][1]],
            [matrix[0][2], matrix[1][2], matrix[2][2], matrix[3][2], matrix[4][2], matrix[5][2]],
            [matrix[0][3], matrix[1][3], matrix[2][3], matrix[3][3], matrix[4][3], matrix[5][3]],
            [matrix[0][4], matrix[1][4], matrix[2][4], matrix[3][4], matrix[4][4], matrix[5][4]],
            [matrix[0][5], matrix[1][5], matrix[2][5], matrix[3][5], matrix[4][5], matrix[5][5]],
        ];

        assert_matrix_approximately_equal_6(&matrix, &transposed);
    }

    #[test]
    fn rejects_stiffness_matrix_for_zero_area_triangle() {
        let material = Material2D::new(200.0, 0.3, 7800.0).expect("valid material should be created");
        let first_node = Node2D::new(1, 0.0, 0.0).expect("valid node should be created");
        let second_node = Node2D::new(2, 1.0, 1.0).expect("valid node should be created");
        let third_node = Node2D::new(3, 2.0, 2.0).expect("valid node should be created");
        let triangle = TriangleT3::new(30, [1, 2, 3], 0.2).expect("valid connectivity should be created");

        let result = triangle.stiffness_matrix(&material, &first_node, &second_node, &third_node);

        assert!(matches!(
            result,
            Err(FemError::DegenerateElement {
                element_id: 30,
                element_type: "triangle_t3",
                measure_name: "area",
                measure: 0.0,
                ..
            })
        ));
    }

    #[test]
    fn exposes_element_cross_section_area() {
        let truss = Truss2D::new(10, [1, 2], 0.02).expect("valid truss should be created");
        let beam = Beam2D::new(20, [2, 3], 0.03, 0.001).expect("valid beam should be created");

        let cases = [(Element2D::Truss(truss), Some(0.02)), (Element2D::Beam(beam), Some(0.03))];

        assert_eq!(truss.cross_section_area(), 0.02);
        assert_eq!(beam.cross_section_area(), 0.03);
        assert_eq!(beam.second_moment_of_area(), 0.001);

        for (element, expected_area) in cases {
            assert_eq!(element.cross_section_area(), expected_area);
        }
    }

    #[test]
    fn calculates_stiffness_matrix_for_horizontal_truss() {
        let material = Material2D::new(200.0, 0.3, 7800.0).expect("valid material should be created");
        let first_node = Node2D::new(1, 0.0, 0.0).expect("valid node should be created");
        let second_node = Node2D::new(2, 2.0, 0.0).expect("valid node should be created");
        let truss = Truss2D::new(10, [1, 2], 0.5).expect("valid truss should be created");

        let matrix = truss
            .stiffness_matrix(&material, &first_node, &second_node)
            .expect("stiffness matrix should be calculated");

        assert_eq!(
            matrix,
            [[50.0, 0.0, -50.0, 0.0], [0.0, 0.0, 0.0, 0.0], [-50.0, 0.0, 50.0, 0.0], [0.0, 0.0, 0.0, 0.0],]
        );
    }

    #[test]
    fn calculates_stiffness_matrix_for_horizontal_beam() {
        let material = Material2D::new(200.0, 0.3, 7800.0).expect("valid material should be created");
        let first_node = Node2D::new(1, 0.0, 0.0).expect("valid node should be created");
        let second_node = Node2D::new(2, 1.0, 0.0).expect("valid node should be created");
        let beam = Beam2D::new(20, [1, 2], 2.0, 3.0).expect("valid beam should be created");

        let matrix =
            beam.stiffness_matrix(&material, &first_node, &second_node).expect("stiffness matrix should be calculated");

        assert_eq!(
            matrix,
            [
                [400.0, 0.0, 0.0, -400.0, 0.0, 0.0],
                [0.0, 7200.0, 3600.0, 0.0, -7200.0, 3600.0],
                [0.0, 3600.0, 2400.0, 0.0, -3600.0, 1200.0],
                [-400.0, 0.0, 0.0, 400.0, 0.0, 0.0],
                [0.0, -7200.0, -3600.0, 0.0, 7200.0, -3600.0],
                [0.0, 3600.0, 1200.0, 0.0, -3600.0, 2400.0],
            ]
        );
    }

    #[test]
    fn calculates_stiffness_matrix_for_longer_beam() {
        let material = Material2D::new(200.0, 0.3, 7800.0).expect("valid material should be created");
        let first_node = Node2D::new(1, 0.0, 0.0).expect("valid node should be created");
        let second_node = Node2D::new(2, 2.0, 0.0).expect("valid node should be created");
        let beam = Beam2D::new(20, [1, 2], 2.0, 3.0).expect("valid beam should be created");

        let matrix =
            beam.stiffness_matrix(&material, &first_node, &second_node).expect("stiffness matrix should be calculated");

        let expected = [
            [200.0, 0.0, 0.0, -200.0, 0.0, 0.0],
            [0.0, 900.0, 900.0, 0.0, -900.0, 900.0],
            [0.0, 900.0, 1200.0, 0.0, -900.0, 600.0],
            [-200.0, 0.0, 0.0, 200.0, 0.0, 0.0],
            [0.0, -900.0, -900.0, 0.0, 900.0, -900.0],
            [0.0, 900.0, 600.0, 0.0, -900.0, 1200.0],
        ];

        assert_matrix_approximately_equal_6(&matrix, &expected);
    }

    #[test]
    fn calculates_global_stiffness_matrix_for_diagonal_beam() {
        let material = Material2D::new(1.0, 0.3, 1.0).expect("valid material should be created");
        let first_node = Node2D::new(1, 0.0, 0.0).expect("valid node should be created");
        let second_node = Node2D::new(2, 1.0, 1.0).expect("valid node should be created");
        let beam = Beam2D::new(20, [1, 2], 1.0, 1.0).expect("valid beam should be created");

        let matrix =
            beam.stiffness_matrix(&material, &first_node, &second_node).expect("stiffness matrix should be calculated");

        let length = 2.0_f64.sqrt();
        let c = 1.0 / length;
        let s = c;
        let axial = 1.0 / length;
        let bending = 12.0 / length.powi(3);
        let coupling = 6.0 / length.powi(2);
        let rotation = 4.0 / length;
        let rotation_coupling = 2.0 / length;
        let k00 = axial * c * c + bending * s * s;
        let k01 = (axial - bending) * c * s;
        let k02 = -coupling * s;
        let k11 = axial * s * s + bending * c * c;
        let k12 = coupling * c;

        let expected = [
            [k00, k01, k02, -k00, -k01, k02],
            [k01, k11, k12, -k01, -k11, k12],
            [k02, k12, rotation, -k02, -k12, rotation_coupling],
            [-k00, -k01, -k02, k00, k01, -k02],
            [-k01, -k11, -k12, k01, k11, -k12],
            [k02, k12, rotation_coupling, -k02, -k12, rotation],
        ];

        assert_matrix_approximately_equal_6(&matrix, &expected);
    }

    #[test]
    fn beam_stiffness_matrix_is_symmetric() {
        let material = Material2D::new(200.0, 0.3, 7800.0).expect("valid material should be created");
        let first_node = Node2D::new(1, 0.0, 0.0).expect("valid node should be created");
        let second_node = Node2D::new(2, 2.0, 0.0).expect("valid node should be created");
        let beam = Beam2D::new(20, [1, 2], 2.0, 3.0).expect("valid beam should be created");

        let matrix =
            beam.stiffness_matrix(&material, &first_node, &second_node).expect("stiffness matrix should be calculated");

        let transposed = [
            [matrix[0][0], matrix[1][0], matrix[2][0], matrix[3][0], matrix[4][0], matrix[5][0]],
            [matrix[0][1], matrix[1][1], matrix[2][1], matrix[3][1], matrix[4][1], matrix[5][1]],
            [matrix[0][2], matrix[1][2], matrix[2][2], matrix[3][2], matrix[4][2], matrix[5][2]],
            [matrix[0][3], matrix[1][3], matrix[2][3], matrix[3][3], matrix[4][3], matrix[5][3]],
            [matrix[0][4], matrix[1][4], matrix[2][4], matrix[3][4], matrix[4][4], matrix[5][4]],
            [matrix[0][5], matrix[1][5], matrix[2][5], matrix[3][5], matrix[4][5], matrix[5][5]],
        ];

        assert_matrix_approximately_equal_6(&matrix, &transposed);
    }

    #[test]
    fn rejects_stiffness_matrix_for_zero_length_beam() {
        let material = Material2D::new(200.0, 0.3, 7800.0).expect("valid material should be created");
        let first_node = Node2D::new(1, 0.0, 0.0).expect("valid node should be created");
        let second_node = Node2D::new(2, 0.0, 0.0).expect("valid node should be created");
        let beam = Beam2D::new(20, [1, 2], 2.0, 3.0).expect("valid connectivity should be created");

        let result = beam.stiffness_matrix(&material, &first_node, &second_node);

        assert!(matches!(
            result,
            Err(FemError::DegenerateElement {
                element_id: 20,
                element_type: "beam",
                measure_name: "length",
                measure: 0.0,
                ..
            })
        ));
    }

    #[test]
    fn calculates_stiffness_matrix_for_diagonal_truss() {
        let material = Material2D::new(200.0, 0.3, 7800.0).expect("valid material should be created");
        let first_node = Node2D::new(1, 0.0, 0.0).expect("valid node should be created");
        let second_node = Node2D::new(2, 1.0, 1.0).expect("valid node should be created");
        let truss = Truss2D::new(10, [1, 2], 0.5).expect("valid truss should be created");

        let matrix = truss
            .stiffness_matrix(&material, &first_node, &second_node)
            .expect("stiffness matrix should be calculated");

        let value = 25.0 * 2.0_f64.sqrt();
        let expected = [
            [value, value, -value, -value],
            [value, value, -value, -value],
            [-value, -value, value, value],
            [-value, -value, value, value],
        ];

        assert_matrix_approximately_equal(&matrix, &expected);
    }

    #[test]
    fn stiffness_matrix_is_symmetric() {
        let material = Material2D::new(200.0, 0.3, 7800.0).expect("valid material should be created");
        let first_node = Node2D::new(1, 0.0, 0.0).expect("valid node should be created");
        let second_node = Node2D::new(2, 1.0, 1.0).expect("valid node should be created");
        let truss = Truss2D::new(10, [1, 2], 0.5).expect("valid truss should be created");

        let matrix = truss
            .stiffness_matrix(&material, &first_node, &second_node)
            .expect("stiffness matrix should be calculated");

        let transposed = [
            [matrix[0][0], matrix[1][0], matrix[2][0], matrix[3][0]],
            [matrix[0][1], matrix[1][1], matrix[2][1], matrix[3][1]],
            [matrix[0][2], matrix[1][2], matrix[2][2], matrix[3][2]],
            [matrix[0][3], matrix[1][3], matrix[2][3], matrix[3][3]],
        ];

        assert_matrix_approximately_equal(&matrix, &transposed);
    }

    #[test]
    fn rejects_stiffness_matrix_for_zero_length_truss() {
        let material = Material2D::new(200.0, 0.3, 7800.0).expect("valid material should be created");
        let first_node = Node2D::new(1, 0.0, 0.0).expect("valid node should be created");
        let second_node = Node2D::new(2, 0.0, 0.0).expect("valid node should be created");
        let truss = Truss2D::new(10, [1, 2], 0.5).expect("valid connectivity should be created");

        let result = truss.stiffness_matrix(&material, &first_node, &second_node);

        assert!(matches!(
            result,
            Err(FemError::DegenerateElement {
                element_id: 10,
                element_type: "truss",
                measure_name: "length",
                measure: 0.0,
                ..
            })
        ));
    }

    fn assert_matrix_approximately_equal(actual: &[[f64; 4]; 4], expected: &[[f64; 4]; 4]) {
        for (row, (actual_row, expected_row)) in actual.iter().zip(expected).enumerate() {
            for (column, (actual_value, expected_value)) in actual_row.iter().zip(expected_row).enumerate() {
                assert!(
                    (actual_value - expected_value).abs() < 1e-12,
                    "different matrix entry at row {row}, column {column}: actual = {}, expected = {}",
                    actual_value,
                    expected_value
                );
            }
        }
    }

    fn assert_matrix_approximately_equal_6(actual: &[[f64; 6]; 6], expected: &[[f64; 6]; 6]) {
        for (row, (actual_row, expected_row)) in actual.iter().zip(expected).enumerate() {
            for (column, (actual_value, expected_value)) in actual_row.iter().zip(expected_row).enumerate() {
                assert!(
                    (actual_value - expected_value).abs() < 1e-12,
                    "different matrix entry at row {row}, column {column}: actual = {}, expected = {}",
                    actual_value,
                    expected_value
                );
            }
        }
    }
}
