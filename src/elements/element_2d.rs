//! Defines 2D elements used in finite element analysis.

use super::Interpolation;
use super::interpolation::{
    quad_q4_shape_function_derivatives, quad_q8_shape_function_derivatives, triangle_t6_shape_function_derivatives,
};
use crate::error::FemError;
use crate::model::{BeamSection2D, Dof2D, Material2D, Node2D, PlaneStressSection2D, Section2D, TrussSection2D};
use nalgebra::DMatrix;

const TRANSLATIONAL_DOFS: &[Dof2D] = &[Dof2D::Ux, Dof2D::Uy];

const FRAME_DOFS: &[Dof2D] = &[Dof2D::Ux, Dof2D::Uy, Dof2D::Rz];

/// Truss element in 2D space, defined by two nodes.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Truss2D {
    id: usize,
    node_ids: [usize; 2],
    material_id: usize,
    section_id: usize,
}

/// Beam element in 2D space, defined by two nodes.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Beam2D {
    id: usize,
    node_ids: [usize; 2],
    material_id: usize,
    section_id: usize,
}

/// Triangle element in 2D space, defined by three nodes.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TriangleT3 {
    id: usize,
    node_ids: [usize; 3],
    material_id: usize,
    section_id: usize,
}

/// Quadratic triangle element in 2D space, defined by six nodes.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TriangleT6 {
    id: usize,
    node_ids: [usize; 6],
    material_id: usize,
    section_id: usize,
}

/// Quadrilateral element in 2D space, defined by four nodes.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct QuadQ4 {
    id: usize,
    node_ids: [usize; 4],
    material_id: usize,
    section_id: usize,
}

/// Quadratic serendipity quadrilateral element in 2D space, defined by eight nodes.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct QuadQ8 {
    id: usize,
    node_ids: [usize; 8],
    material_id: usize,
    section_id: usize,
}

/// Enum representing different types of 2D elements used in finite element analysis.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Element2D {
    Truss(Truss2D),
    Beam(Beam2D),
    TriangleT3(TriangleT3),
    TriangleT6(TriangleT6),
    QuadQ4(QuadQ4),
    QuadQ8(QuadQ8),
}

impl Truss2D {
    /// Creates a new Truss2D element with the specified ID, node IDs, material ID, and section ID.
    pub fn new(id: usize, node_ids: [usize; 2], material_id: usize, section_id: usize) -> Result<Self, FemError> {
        if node_ids[0] == node_ids[1] {
            return Err(FemError::InvalidElementConnectivity { element_id: id, node_ids: node_ids.to_vec() });
        }

        Ok(Self { id, node_ids, material_id, section_id })
    }

    /// Calculates the stiffness matrix for the truss element based on the provided material properties and node coordinates.
    pub fn stiffness_matrix(
        &self, material: &Material2D, section: &TrussSection2D, first_node: &Node2D, second_node: &Node2D,
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

        let factor = material.young_modulus() * section.cross_section_area() / length;

        let matrix = [
            [factor * c * c, factor * c * s, -factor * c * c, -factor * c * s],
            [factor * c * s, factor * s * s, -factor * c * s, -factor * s * s],
            [-factor * c * c, -factor * c * s, factor * c * c, factor * c * s],
            [-factor * c * s, -factor * s * s, factor * c * s, factor * s * s],
        ];

        Ok(matrix)
    }

    /// Returns the material ID used by this element.
    #[must_use]
    pub fn material_id(&self) -> usize {
        self.material_id
    }

    /// Returns the section ID used by this element.
    #[must_use]
    pub fn section_id(&self) -> usize {
        self.section_id
    }

    /// Returns the distance between the truss's two nodes.
    pub fn length(&self, first_node: &Node2D, second_node: &Node2D) -> Result<f64, FemError> {
        let dx = second_node.x() - first_node.x();
        let dy = second_node.y() - first_node.y();
        let length = (dx * dx + dy * dy).sqrt();

        if !length.is_finite() || length == 0.0 {
            return Err(FemError::DegenerateElement {
                element_id: self.id,
                element_type: "truss",
                node_ids: self.node_ids.to_vec(),
                measure_name: "length",
                measure: length,
            });
        }

        Ok(length)
    }
}

impl Beam2D {
    /// Creates a new Beam2D element with the specified ID, node IDs, material ID, and section ID.
    pub fn new(id: usize, node_ids: [usize; 2], material_id: usize, section_id: usize) -> Result<Self, FemError> {
        if node_ids[0] == node_ids[1] {
            return Err(FemError::InvalidElementConnectivity { element_id: id, node_ids: node_ids.to_vec() });
        }

        Ok(Self { id, node_ids, material_id, section_id })
    }

    /// Returns the material ID used by this element.
    #[must_use]
    pub fn material_id(&self) -> usize {
        self.material_id
    }

    /// Returns the section ID used by this element.
    #[must_use]
    pub fn section_id(&self) -> usize {
        self.section_id
    }

    pub(crate) fn geometry(&self, first_node: &Node2D, second_node: &Node2D) -> Result<(f64, f64, f64), FemError> {
        let dx = second_node.x() - first_node.x();
        let dy = second_node.y() - first_node.y();
        let length = (dx * dx + dy * dy).sqrt();

        if !length.is_finite() || length == 0.0 {
            return Err(FemError::DegenerateElement {
                element_id: self.id,
                element_type: "beam",
                node_ids: self.node_ids.to_vec(),
                measure_name: "length",
                measure: length,
            });
        }

        Ok((length, dx / length, dy / length))
    }

    /// Returns the distance between the beam's two nodes.
    pub fn length(&self, first_node: &Node2D, second_node: &Node2D) -> Result<f64, FemError> {
        self.geometry(first_node, second_node).map(|(length, _, _)| length)
    }

    pub(crate) fn local_stiffness_matrix(
        &self, material: &Material2D, section: &BeamSection2D, length: f64,
    ) -> [[f64; 6]; 6] {
        let ea_over_l = material.young_modulus() * section.cross_section_area() / length;
        let twelve_ei_over_l3 = 12.0 * material.young_modulus() * section.second_moment_of_area() / length.powi(3);
        let six_ei_over_l2 = 6.0 * material.young_modulus() * section.second_moment_of_area() / length.powi(2);
        let four_ei_over_l = 4.0 * material.young_modulus() * section.second_moment_of_area() / length;
        let two_ei_over_l = 2.0 * material.young_modulus() * section.second_moment_of_area() / length;

        [
            [ea_over_l, 0.0, 0.0, -ea_over_l, 0.0, 0.0],
            [0.0, twelve_ei_over_l3, six_ei_over_l2, 0.0, -twelve_ei_over_l3, six_ei_over_l2],
            [0.0, six_ei_over_l2, four_ei_over_l, 0.0, -six_ei_over_l2, two_ei_over_l],
            [-ea_over_l, 0.0, 0.0, ea_over_l, 0.0, 0.0],
            [0.0, -twelve_ei_over_l3, -six_ei_over_l2, 0.0, twelve_ei_over_l3, -six_ei_over_l2],
            [0.0, six_ei_over_l2, two_ei_over_l, 0.0, -six_ei_over_l2, four_ei_over_l],
        ]
    }

    /// Calculates the stiffness matrix of the beam element in global coordinates.
    pub fn stiffness_matrix(
        &self, material: &Material2D, section: &BeamSection2D, first_node: &Node2D, second_node: &Node2D,
    ) -> Result<[[f64; 6]; 6], FemError> {
        let (length, c, s) = self.geometry(first_node, second_node)?;
        let local_matrix = self.local_stiffness_matrix(material, section, length);

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
    pub fn new(id: usize, node_ids: [usize; 3], material_id: usize, section_id: usize) -> Result<Self, FemError> {
        if node_ids[0] == node_ids[1] || node_ids[0] == node_ids[2] || node_ids[1] == node_ids[2] {
            return Err(FemError::InvalidElementConnectivity { element_id: id, node_ids: node_ids.to_vec() });
        }

        Ok(Self { id, node_ids, material_id, section_id })
    }

    /// Returns the material ID used by this element.
    #[must_use]
    pub fn material_id(&self) -> usize {
        self.material_id
    }

    /// Returns the section ID used by this element.
    #[must_use]
    pub fn section_id(&self) -> usize {
        self.section_id
    }

    pub(crate) fn strain_displacement_matrix(
        &self, first_node: &Node2D, second_node: &Node2D, third_node: &Node2D,
    ) -> Result<([[f64; 6]; 3], f64), FemError> {
        let x1 = first_node.x();
        let y1 = first_node.y();
        let x2 = second_node.x();
        let y2 = second_node.y();
        let x3 = third_node.x();
        let y3 = third_node.y();

        let twice_signed_area = (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1);
        let area = 0.5 * twice_signed_area.abs();

        if !area.is_finite() || area == 0.0 {
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

        let matrix = [
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

        Ok((matrix, area))
    }

    pub(crate) fn constitutive_matrix(material: &Material2D) -> [[f64; 3]; 3] {
        let constitutive_factor =
            material.young_modulus() / (1.0 - material.poisson_ratio() * material.poisson_ratio());

        [
            [constitutive_factor, constitutive_factor * material.poisson_ratio(), 0.0],
            [constitutive_factor * material.poisson_ratio(), constitutive_factor, 0.0],
            [0.0, 0.0, constitutive_factor * (1.0 - material.poisson_ratio()) / 2.0],
        ]
    }

    /// Calculates the stiffness matrix using the linear plane-stress T3 formulation.
    pub fn stiffness_matrix(
        &self, material: &Material2D, section: &PlaneStressSection2D, first_node: &Node2D, second_node: &Node2D,
        third_node: &Node2D,
    ) -> Result<[[f64; 6]; 6], FemError> {
        let (strain_displacement_matrix, area) =
            self.strain_displacement_matrix(first_node, second_node, third_node)?;
        let constitutive_matrix = Self::constitutive_matrix(material);

        let constitutive_times_strain = multiply_3x3_by_3x6(&constitutive_matrix, &strain_displacement_matrix);
        let mut stiffness_matrix =
            multiply_transpose_3x6_by_3x6(&strain_displacement_matrix, &constitutive_times_strain);
        let scale = section.thickness() * area;

        for row in &mut stiffness_matrix {
            for value in row {
                *value *= scale;
            }
        }

        Ok(stiffness_matrix)
    }
}

impl TriangleT6 {
    /// Creates a six-node quadratic triangular element for a plane-stress analysis.
    pub fn new(id: usize, node_ids: [usize; 6], material_id: usize, section_id: usize) -> Result<Self, FemError> {
        for first_index in 0..node_ids.len() {
            for second_index in (first_index + 1)..node_ids.len() {
                if node_ids[first_index] == node_ids[second_index] {
                    return Err(FemError::InvalidElementConnectivity { element_id: id, node_ids: node_ids.to_vec() });
                }
            }
        }

        Ok(Self { id, node_ids, material_id, section_id })
    }

    /// Returns the material ID used by this element.
    #[must_use]
    pub fn material_id(&self) -> usize {
        self.material_id
    }

    /// Returns the section ID used by this element.
    #[must_use]
    pub fn section_id(&self) -> usize {
        self.section_id
    }

    pub(crate) fn gauss_points() -> [(f64, f64, f64); 6] {
        let vertex = 0.816_847_572_980_459;
        let near_edge = 0.091_576_213_509_771;
        let near_vertex = 0.108_103_018_168_070;
        let mid = 0.445_948_490_915_965;
        let vertex_weight = 0.054_975_871_827_661;
        let mid_weight = 0.111_690_794_839_005;

        [
            (near_edge, near_edge, vertex_weight),
            (vertex, near_edge, vertex_weight),
            (near_edge, vertex, vertex_weight),
            (mid, mid, mid_weight),
            (near_vertex, mid, mid_weight),
            (mid, near_vertex, mid_weight),
        ]
    }

    pub(crate) fn strain_displacement_matrix(
        &self, nodes: [&Node2D; 6], xi: f64, eta: f64,
    ) -> Result<([[f64; 12]; 3], f64), FemError> {
        let derivatives = triangle_t6_shape_function_derivatives(xi, eta);
        let dndxi = derivatives[0];
        let dndeta = derivatives[1];
        let mut dx_dxi = 0.0;
        let mut dx_deta = 0.0;
        let mut dy_dxi = 0.0;
        let mut dy_deta = 0.0;

        for node_index in 0..6 {
            dx_dxi += dndxi[node_index] * nodes[node_index].x();
            dx_deta += dndeta[node_index] * nodes[node_index].x();
            dy_dxi += dndxi[node_index] * nodes[node_index].y();
            dy_deta += dndeta[node_index] * nodes[node_index].y();
        }

        let jacobian_determinant = dx_dxi * dy_deta - dx_deta * dy_dxi;

        if !jacobian_determinant.is_finite() || jacobian_determinant <= 0.0 {
            return Err(FemError::DegenerateElement {
                element_id: self.id,
                element_type: "triangle_t6",
                node_ids: self.node_ids.to_vec(),
                measure_name: "jacobian determinant",
                measure: jacobian_determinant,
            });
        }

        let mut dndx = [0.0; 6];
        let mut dndy = [0.0; 6];

        for node_index in 0..6 {
            dndx[node_index] = (dy_deta * dndxi[node_index] - dy_dxi * dndeta[node_index]) / jacobian_determinant;
            dndy[node_index] = (-dx_deta * dndxi[node_index] + dx_dxi * dndeta[node_index]) / jacobian_determinant;
        }

        let mut matrix = [[0.0; 12]; 3];

        for node_index in 0..6 {
            let x_dof = 2 * node_index;
            let y_dof = x_dof + 1;

            matrix[0][x_dof] = dndx[node_index];
            matrix[1][y_dof] = dndy[node_index];
            matrix[2][x_dof] = dndy[node_index];
            matrix[2][y_dof] = dndx[node_index];
        }

        Ok((matrix, jacobian_determinant))
    }

    /// Calculates the stiffness matrix using a quadratic T6 plane-stress formulation.
    pub fn stiffness_matrix(
        &self, material: &Material2D, section: &PlaneStressSection2D, nodes: [&Node2D; 6],
    ) -> Result<[[f64; 12]; 12], FemError> {
        let constitutive_matrix = TriangleT3::constitutive_matrix(material);
        let mut stiffness_matrix = [[0.0; 12]; 12];

        for (xi, eta, weight) in Self::gauss_points() {
            let (strain_displacement_matrix, jacobian_determinant) = self.strain_displacement_matrix(nodes, xi, eta)?;
            let constitutive_times_strain = multiply_3x3_by_3x12(&constitutive_matrix, &strain_displacement_matrix);
            let contribution = multiply_transpose_3x12_by_3x12(&strain_displacement_matrix, &constitutive_times_strain);
            let scale = section.thickness() * jacobian_determinant * weight;

            for (row_index, stiffness_row) in stiffness_matrix.iter_mut().enumerate() {
                for (column_index, stiffness_value) in stiffness_row.iter_mut().enumerate() {
                    *stiffness_value += contribution[row_index][column_index] * scale;
                }
            }
        }

        Ok(stiffness_matrix)
    }
}

impl QuadQ4 {
    /// Creates a four-node quadrilateral element for a plane-stress analysis.
    pub fn new(id: usize, node_ids: [usize; 4], material_id: usize, section_id: usize) -> Result<Self, FemError> {
        for first_index in 0..node_ids.len() {
            for second_index in (first_index + 1)..node_ids.len() {
                if node_ids[first_index] == node_ids[second_index] {
                    return Err(FemError::InvalidElementConnectivity { element_id: id, node_ids: node_ids.to_vec() });
                }
            }
        }

        Ok(Self { id, node_ids, material_id, section_id })
    }

    /// Returns the material ID used by this element.
    #[must_use]
    pub fn material_id(&self) -> usize {
        self.material_id
    }

    /// Returns the section ID used by this element.
    #[must_use]
    pub fn section_id(&self) -> usize {
        self.section_id
    }

    pub(crate) fn gauss_points() -> [(f64, f64); 4] {
        let point = 1.0 / 3.0_f64.sqrt();

        [(-point, -point), (point, -point), (point, point), (-point, point)]
    }

    pub(crate) fn strain_displacement_matrix(
        &self, nodes: [&Node2D; 4], xi: f64, eta: f64,
    ) -> Result<([[f64; 8]; 3], f64), FemError> {
        let derivatives = quad_q4_shape_function_derivatives(xi, eta);
        let dndxi = derivatives[0];
        let dndeta = derivatives[1];
        let mut dx_dxi = 0.0;
        let mut dx_deta = 0.0;
        let mut dy_dxi = 0.0;
        let mut dy_deta = 0.0;

        for node_index in 0..4 {
            dx_dxi += dndxi[node_index] * nodes[node_index].x();
            dx_deta += dndeta[node_index] * nodes[node_index].x();
            dy_dxi += dndxi[node_index] * nodes[node_index].y();
            dy_deta += dndeta[node_index] * nodes[node_index].y();
        }

        let jacobian_determinant = dx_dxi * dy_deta - dx_deta * dy_dxi;

        if !jacobian_determinant.is_finite() || jacobian_determinant <= 0.0 {
            return Err(FemError::DegenerateElement {
                element_id: self.id,
                element_type: "quad_q4",
                node_ids: self.node_ids.to_vec(),
                measure_name: "jacobian determinant",
                measure: jacobian_determinant,
            });
        }

        let mut dndx = [0.0; 4];
        let mut dndy = [0.0; 4];

        for node_index in 0..4 {
            dndx[node_index] = (dy_deta * dndxi[node_index] - dy_dxi * dndeta[node_index]) / jacobian_determinant;
            dndy[node_index] = (-dx_deta * dndxi[node_index] + dx_dxi * dndeta[node_index]) / jacobian_determinant;
        }

        let mut matrix = [[0.0; 8]; 3];

        for node_index in 0..4 {
            let x_dof = 2 * node_index;
            let y_dof = x_dof + 1;

            matrix[0][x_dof] = dndx[node_index];
            matrix[1][y_dof] = dndy[node_index];
            matrix[2][x_dof] = dndy[node_index];
            matrix[2][y_dof] = dndx[node_index];
        }

        Ok((matrix, jacobian_determinant))
    }

    /// Calculates the stiffness matrix using a bilinear Q4 plane-stress formulation.
    pub fn stiffness_matrix(
        &self, material: &Material2D, section: &PlaneStressSection2D, first_node: &Node2D, second_node: &Node2D,
        third_node: &Node2D, fourth_node: &Node2D,
    ) -> Result<[[f64; 8]; 8], FemError> {
        let nodes = [first_node, second_node, third_node, fourth_node];
        let constitutive_matrix = TriangleT3::constitutive_matrix(material);
        let mut stiffness_matrix = [[0.0; 8]; 8];

        for (xi, eta) in Self::gauss_points() {
            let (strain_displacement_matrix, jacobian_determinant) = self.strain_displacement_matrix(nodes, xi, eta)?;
            let constitutive_times_strain = multiply_3x3_by_3x8(&constitutive_matrix, &strain_displacement_matrix);
            let contribution = multiply_transpose_3x8_by_3x8(&strain_displacement_matrix, &constitutive_times_strain);
            let scale = section.thickness() * jacobian_determinant;

            for (row_index, stiffness_row) in stiffness_matrix.iter_mut().enumerate() {
                for (column_index, stiffness_value) in stiffness_row.iter_mut().enumerate() {
                    *stiffness_value += contribution[row_index][column_index] * scale;
                }
            }
        }

        Ok(stiffness_matrix)
    }
}

impl QuadQ8 {
    /// Creates an eight-node serendipity quadrilateral element for a plane-stress analysis.
    pub fn new(id: usize, node_ids: [usize; 8], material_id: usize, section_id: usize) -> Result<Self, FemError> {
        for first_index in 0..node_ids.len() {
            for second_index in (first_index + 1)..node_ids.len() {
                if node_ids[first_index] == node_ids[second_index] {
                    return Err(FemError::InvalidElementConnectivity { element_id: id, node_ids: node_ids.to_vec() });
                }
            }
        }

        Ok(Self { id, node_ids, material_id, section_id })
    }

    /// Returns the material ID used by this element.
    #[must_use]
    pub fn material_id(&self) -> usize {
        self.material_id
    }

    /// Returns the section ID used by this element.
    #[must_use]
    pub fn section_id(&self) -> usize {
        self.section_id
    }

    pub(crate) fn gauss_points() -> [(f64, f64, f64); 9] {
        let point = (3.0_f64 / 5.0).sqrt();
        let edge_weight = 5.0 / 9.0;
        let center_weight = 8.0 / 9.0;

        [
            (-point, -point, edge_weight * edge_weight),
            (0.0, -point, center_weight * edge_weight),
            (point, -point, edge_weight * edge_weight),
            (-point, 0.0, edge_weight * center_weight),
            (0.0, 0.0, center_weight * center_weight),
            (point, 0.0, edge_weight * center_weight),
            (-point, point, edge_weight * edge_weight),
            (0.0, point, center_weight * edge_weight),
            (point, point, edge_weight * edge_weight),
        ]
    }

    pub(crate) fn strain_displacement_matrix(
        &self, nodes: [&Node2D; 8], xi: f64, eta: f64,
    ) -> Result<([[f64; 16]; 3], f64), FemError> {
        let derivatives = quad_q8_shape_function_derivatives(xi, eta);
        let dndxi = derivatives[0];
        let dndeta = derivatives[1];
        let mut dx_dxi = 0.0;
        let mut dx_deta = 0.0;
        let mut dy_dxi = 0.0;
        let mut dy_deta = 0.0;

        for node_index in 0..8 {
            dx_dxi += dndxi[node_index] * nodes[node_index].x();
            dx_deta += dndeta[node_index] * nodes[node_index].x();
            dy_dxi += dndxi[node_index] * nodes[node_index].y();
            dy_deta += dndeta[node_index] * nodes[node_index].y();
        }

        let jacobian_determinant = dx_dxi * dy_deta - dx_deta * dy_dxi;

        if !jacobian_determinant.is_finite() || jacobian_determinant <= 0.0 {
            return Err(FemError::DegenerateElement {
                element_id: self.id,
                element_type: "quad_q8",
                node_ids: self.node_ids.to_vec(),
                measure_name: "jacobian determinant",
                measure: jacobian_determinant,
            });
        }

        let mut dndx = [0.0; 8];
        let mut dndy = [0.0; 8];

        for node_index in 0..8 {
            dndx[node_index] = (dy_deta * dndxi[node_index] - dy_dxi * dndeta[node_index]) / jacobian_determinant;
            dndy[node_index] = (-dx_deta * dndxi[node_index] + dx_dxi * dndeta[node_index]) / jacobian_determinant;
        }

        let mut matrix = [[0.0; 16]; 3];

        for node_index in 0..8 {
            let x_dof = 2 * node_index;
            let y_dof = x_dof + 1;

            matrix[0][x_dof] = dndx[node_index];
            matrix[1][y_dof] = dndy[node_index];
            matrix[2][x_dof] = dndy[node_index];
            matrix[2][y_dof] = dndx[node_index];
        }

        Ok((matrix, jacobian_determinant))
    }

    /// Calculates the stiffness matrix using a quadratic serendipity Q8 plane-stress formulation.
    pub fn stiffness_matrix(
        &self, material: &Material2D, section: &PlaneStressSection2D, nodes: [&Node2D; 8],
    ) -> Result<[[f64; 16]; 16], FemError> {
        let constitutive_matrix = TriangleT3::constitutive_matrix(material);
        let mut stiffness_matrix = [[0.0; 16]; 16];

        for (xi, eta, weight) in Self::gauss_points() {
            let (strain_displacement_matrix, jacobian_determinant) = self.strain_displacement_matrix(nodes, xi, eta)?;
            let constitutive_times_strain = multiply_3x3_by_3x16(&constitutive_matrix, &strain_displacement_matrix);
            let contribution = multiply_transpose_3x16_by_3x16(&strain_displacement_matrix, &constitutive_times_strain);
            let scale = section.thickness() * jacobian_determinant * weight;

            for (row_index, stiffness_row) in stiffness_matrix.iter_mut().enumerate() {
                for (column_index, stiffness_value) in stiffness_row.iter_mut().enumerate() {
                    *stiffness_value += contribution[row_index][column_index] * scale;
                }
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

/// Multiplies a 3x3 matrix by a 3x8 matrix.
fn multiply_3x3_by_3x8(left: &[[f64; 3]; 3], right: &[[f64; 8]; 3]) -> [[f64; 8]; 3] {
    let mut product = [[0.0; 8]; 3];

    for (row, product_row) in product.iter_mut().enumerate() {
        for (column, product_value) in product_row.iter_mut().enumerate() {
            *product_value =
                left[row].iter().enumerate().map(|(index, left_value)| left_value * right[index][column]).sum();
        }
    }

    product
}

/// Multiplies a 3x3 matrix by a 3x12 matrix.
fn multiply_3x3_by_3x12(left: &[[f64; 3]; 3], right: &[[f64; 12]; 3]) -> [[f64; 12]; 3] {
    let mut product = [[0.0; 12]; 3];

    for (row, product_row) in product.iter_mut().enumerate() {
        for (column, product_value) in product_row.iter_mut().enumerate() {
            *product_value =
                left[row].iter().enumerate().map(|(index, left_value)| left_value * right[index][column]).sum();
        }
    }

    product
}

/// Multiplies a 3x3 matrix by a 3x16 matrix.
fn multiply_3x3_by_3x16(left: &[[f64; 3]; 3], right: &[[f64; 16]; 3]) -> [[f64; 16]; 3] {
    let mut product = [[0.0; 16]; 3];

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

/// Multiplies the transpose of a 3x8 matrix by a 3x8 matrix.
fn multiply_transpose_3x8_by_3x8(left: &[[f64; 8]; 3], right: &[[f64; 8]; 3]) -> [[f64; 8]; 8] {
    let mut product = [[0.0; 8]; 8];

    for (row, product_row) in product.iter_mut().enumerate() {
        for (column, product_value) in product_row.iter_mut().enumerate() {
            *product_value =
                left.iter().enumerate().map(|(index, left_row)| left_row[row] * right[index][column]).sum();
        }
    }

    product
}

/// Multiplies the transpose of a 3x12 matrix by a 3x12 matrix.
fn multiply_transpose_3x12_by_3x12(left: &[[f64; 12]; 3], right: &[[f64; 12]; 3]) -> [[f64; 12]; 12] {
    let mut product = [[0.0; 12]; 12];

    for (row, product_row) in product.iter_mut().enumerate() {
        for (column, product_value) in product_row.iter_mut().enumerate() {
            *product_value =
                left.iter().enumerate().map(|(index, left_row)| left_row[row] * right[index][column]).sum();
        }
    }

    product
}

/// Multiplies the transpose of a 3x16 matrix by a 3x16 matrix.
fn multiply_transpose_3x16_by_3x16(left: &[[f64; 16]; 3], right: &[[f64; 16]; 3]) -> [[f64; 16]; 16] {
    let mut product = [[0.0; 16]; 16];

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
            Self::TriangleT6(element) => element.id,
            Self::QuadQ4(element) => element.id,
            Self::QuadQ8(element) => element.id,
        }
    }

    /// Returns the element type name for diagnostics.
    #[must_use]
    pub fn element_type(&self) -> &'static str {
        match self {
            Self::Truss(_) => "truss",
            Self::Beam(_) => "beam",
            Self::TriangleT3(_) => "triangle_t3",
            Self::TriangleT6(_) => "triangle_t6",
            Self::QuadQ4(_) => "quad_q4",
            Self::QuadQ8(_) => "quad_q8",
        }
    }

    /// Returns the node IDs associated with the element.
    pub fn node_ids(&self) -> &[usize] {
        match self {
            Self::Truss(element) => &element.node_ids,
            Self::Beam(element) => &element.node_ids,
            Self::TriangleT3(element) => &element.node_ids,
            Self::TriangleT6(element) => &element.node_ids,
            Self::QuadQ4(element) => &element.node_ids,
            Self::QuadQ8(element) => &element.node_ids,
        }
    }

    /// Returns the section ID referenced by this element.
    #[must_use]
    pub fn section_id(&self) -> usize {
        match self {
            Self::Truss(element) => element.section_id(),
            Self::Beam(element) => element.section_id(),
            Self::TriangleT3(element) => element.section_id(),
            Self::TriangleT6(element) => element.section_id(),
            Self::QuadQ4(element) => element.section_id(),
            Self::QuadQ8(element) => element.section_id(),
        }
    }

    /// Returns the material ID referenced by this element.
    #[must_use]
    pub fn material_id(&self) -> usize {
        match self {
            Self::Truss(element) => element.material_id(),
            Self::Beam(element) => element.material_id(),
            Self::TriangleT3(element) => element.material_id(),
            Self::TriangleT6(element) => element.material_id(),
            Self::QuadQ4(element) => element.material_id(),
            Self::QuadQ8(element) => element.material_id(),
        }
    }

    /// Returns the interpolation type used for the element.
    pub fn interpolation(&self) -> Interpolation {
        match self {
            Self::Truss(_) => Interpolation::LinearLagrange,
            Self::Beam(_) => Interpolation::CubicHermite,
            Self::TriangleT3(_) => Interpolation::LinearTriangleT3,
            Self::TriangleT6(_) => Interpolation::QuadraticTriangleT6,
            Self::QuadQ4(_) => Interpolation::BilinearQuadQ4,
            Self::QuadQ8(_) => Interpolation::SerendipityQuadQ8,
        }
    }

    /// Returns the degrees of freedom (DOFs) associated with each node of the element.
    pub fn dofs_per_node(&self) -> &'static [Dof2D] {
        match self {
            Self::Truss(_) => TRANSLATIONAL_DOFS,
            Self::Beam(_) => FRAME_DOFS,
            Self::TriangleT3(_) => TRANSLATIONAL_DOFS,
            Self::TriangleT6(_) => TRANSLATIONAL_DOFS,
            Self::QuadQ4(_) => TRANSLATIONAL_DOFS,
            Self::QuadQ8(_) => TRANSLATIONAL_DOFS,
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
    pub fn stiffness_matrix(
        &self, material: &Material2D, section: &Section2D, nodes: &[Node2D],
    ) -> Result<DMatrix<f64>, FemError> {
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

        match (self, section) {
            (Self::Truss(element), Section2D::Truss(section)) => {
                let matrix = element.stiffness_matrix(material, section, element_nodes[0], element_nodes[1])?;

                Ok(dynamic_matrix_from_array(matrix))
            }
            (Self::Beam(element), Section2D::Beam(section)) => {
                let matrix = element.stiffness_matrix(material, section, element_nodes[0], element_nodes[1])?;

                Ok(dynamic_matrix_from_array(matrix))
            }
            (Self::TriangleT3(element), Section2D::PlaneStress(section)) => {
                let matrix = element.stiffness_matrix(
                    material,
                    section,
                    element_nodes[0],
                    element_nodes[1],
                    element_nodes[2],
                )?;

                Ok(dynamic_matrix_from_array(matrix))
            }
            (Self::TriangleT6(element), Section2D::PlaneStress(section)) => {
                let matrix = element.stiffness_matrix(
                    material,
                    section,
                    [
                        element_nodes[0],
                        element_nodes[1],
                        element_nodes[2],
                        element_nodes[3],
                        element_nodes[4],
                        element_nodes[5],
                    ],
                )?;

                Ok(dynamic_matrix_from_array(matrix))
            }
            (Self::QuadQ4(element), Section2D::PlaneStress(section)) => {
                let matrix = element.stiffness_matrix(
                    material,
                    section,
                    element_nodes[0],
                    element_nodes[1],
                    element_nodes[2],
                    element_nodes[3],
                )?;

                Ok(dynamic_matrix_from_array(matrix))
            }
            (Self::QuadQ8(element), Section2D::PlaneStress(section)) => {
                let matrix = element.stiffness_matrix(
                    material,
                    section,
                    [
                        element_nodes[0],
                        element_nodes[1],
                        element_nodes[2],
                        element_nodes[3],
                        element_nodes[4],
                        element_nodes[5],
                        element_nodes[6],
                        element_nodes[7],
                    ],
                )?;

                Ok(dynamic_matrix_from_array(matrix))
            }
            _ => Err(FemError::InvalidSectionType {
                section_id: self.section_id(),
                expected: self.expected_section_type(),
                actual: section.section_type(),
            }),
        }
    }

    /// Checks that a section has the expected type for this element.
    pub fn validate_section(&self, section: &Section2D) -> Result<(), FemError> {
        match (self, section) {
            (Self::Truss(_), Section2D::Truss(_))
            | (Self::Beam(_), Section2D::Beam(_))
            | (Self::TriangleT3(_), Section2D::PlaneStress(_))
            | (Self::TriangleT6(_), Section2D::PlaneStress(_))
            | (Self::QuadQ4(_), Section2D::PlaneStress(_))
            | (Self::QuadQ8(_), Section2D::PlaneStress(_)) => Ok(()),
            _ => Err(FemError::InvalidSectionType {
                section_id: self.section_id(),
                expected: self.expected_section_type(),
                actual: section.section_type(),
            }),
        }
    }

    /// Returns the expected section type name for diagnostics.
    #[must_use]
    pub fn expected_section_type(&self) -> &'static str {
        match self {
            Self::Truss(_) => "truss",
            Self::Beam(_) => "beam",
            Self::TriangleT3(_) => "plane_stress",
            Self::TriangleT6(_) => "plane_stress",
            Self::QuadQ4(_) => "plane_stress",
            Self::QuadQ8(_) => "plane_stress",
        }
    }
}

/// Converts a square stack-allocated matrix into nalgebra's dynamic matrix.
fn dynamic_matrix_from_array<const N: usize>(matrix: [[f64; N]; N]) -> DMatrix<f64> {
    let values = matrix.into_iter().flatten().collect::<Vec<_>>();

    DMatrix::from_row_slice(N, N, &values)
}

#[cfg(test)]
mod tests {
    use super::{Beam2D, Element2D, QuadQ4, QuadQ8, TriangleT3, TriangleT6, Truss2D, dynamic_matrix_from_array};
    use crate::elements::interpolation::Interpolation;
    use crate::error::FemError;
    use crate::model::{BeamSection2D, Dof2D, Material2D, Node2D, PlaneStressSection2D, Section2D, TrussSection2D};

    #[test]
    fn creates_valid_elements() {
        let cases = [
            ("truss", Truss2D::new(10, [1, 2], 90, 100).map(Element2D::Truss), "truss", 10, 90, 100, vec![1, 2]),
            ("beam", Beam2D::new(20, [2, 3], 91, 200).map(Element2D::Beam), "beam", 20, 91, 200, vec![2, 3]),
            (
                "triangle",
                TriangleT3::new(30, [3, 4, 5], 92, 300).map(Element2D::TriangleT3),
                "triangle_t3",
                30,
                92,
                300,
                vec![3, 4, 5],
            ),
            (
                "triangle_t6",
                TriangleT6::new(35, [3, 4, 5, 6, 7, 8], 95, 350).map(Element2D::TriangleT6),
                "triangle_t6",
                35,
                95,
                350,
                vec![3, 4, 5, 6, 7, 8],
            ),
            (
                "quad",
                QuadQ4::new(40, [4, 5, 6, 7], 93, 400).map(Element2D::QuadQ4),
                "quad_q4",
                40,
                93,
                400,
                vec![4, 5, 6, 7],
            ),
            (
                "quad_q8",
                QuadQ8::new(50, [4, 5, 6, 7, 8, 9, 10, 11], 94, 500).map(Element2D::QuadQ8),
                "quad_q8",
                50,
                94,
                500,
                vec![4, 5, 6, 7, 8, 9, 10, 11],
            ),
        ];

        for (
            name,
            result,
            expected_element_type,
            expected_element_id,
            expected_material_id,
            expected_section_id,
            expected_node_ids,
        ) in cases
        {
            let element = result.expect("valid element should be created");

            assert_eq!(element.id(), expected_element_id, "failed case: {name}");
            assert_eq!(element.element_type(), expected_element_type, "failed case: {name}");
            assert_eq!(element.material_id(), expected_material_id, "failed case: {name}");
            assert_eq!(element.section_id(), expected_section_id, "failed case: {name}");

            assert_eq!(element.node_ids(), expected_node_ids.as_slice(), "failed case: {name}");
        }
    }

    #[test]
    fn rejects_elements_with_repeated_nodes() {
        let cases = [
            ("truss", Truss2D::new(10, [1, 1], 0, 100).map(Element2D::Truss), 10, vec![1, 1]),
            ("beam", Beam2D::new(20, [2, 2], 0, 200).map(Element2D::Beam), 20, vec![2, 2]),
            ("triangle", TriangleT3::new(30, [3, 4, 3], 0, 300).map(Element2D::TriangleT3), 30, vec![3, 4, 3]),
            (
                "triangle_t6",
                TriangleT6::new(35, [3, 4, 5, 6, 7, 3], 0, 350).map(Element2D::TriangleT6),
                35,
                vec![3, 4, 5, 6, 7, 3],
            ),
            ("quad", QuadQ4::new(40, [4, 5, 6, 4], 0, 400).map(Element2D::QuadQ4), 40, vec![4, 5, 6, 4]),
            (
                "quad_q8",
                QuadQ8::new(50, [4, 5, 6, 7, 8, 9, 10, 4], 0, 500).map(Element2D::QuadQ8),
                50,
                vec![4, 5, 6, 7, 8, 9, 10, 4],
            ),
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
                Element2D::Truss(Truss2D::new(10, [1, 2], 0, 100).expect("valid truss should be created")),
                Interpolation::LinearLagrange,
            ),
            (
                "beam",
                Element2D::Beam(Beam2D::new(20, [2, 3], 0, 200).expect("valid beam should be created")),
                Interpolation::CubicHermite,
            ),
            (
                "triangle",
                Element2D::TriangleT3(
                    TriangleT3::new(30, [3, 4, 5], 0, 300).expect("valid triangle should be created"),
                ),
                Interpolation::LinearTriangleT3,
            ),
            (
                "triangle_t6",
                Element2D::TriangleT6(
                    TriangleT6::new(35, [3, 4, 5, 6, 7, 8], 0, 350).expect("valid triangle should be created"),
                ),
                Interpolation::QuadraticTriangleT6,
            ),
            (
                "quad",
                Element2D::QuadQ4(QuadQ4::new(40, [4, 5, 6, 7], 0, 400).expect("valid quad should be created")),
                Interpolation::BilinearQuadQ4,
            ),
            (
                "quad_q8",
                Element2D::QuadQ8(
                    QuadQ8::new(50, [4, 5, 6, 7, 8, 9, 10, 11], 0, 500).expect("valid quad should be created"),
                ),
                Interpolation::SerendipityQuadQ8,
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
                Element2D::Truss(Truss2D::new(10, [1, 2], 0, 100).expect("valid truss should be created")),
                vec![Dof2D::Ux, Dof2D::Uy],
                4,
            ),
            (
                "beam",
                Element2D::Beam(Beam2D::new(20, [2, 3], 0, 200).expect("valid beam should be created")),
                vec![Dof2D::Ux, Dof2D::Uy, Dof2D::Rz],
                6,
            ),
            (
                "triangle",
                Element2D::TriangleT3(
                    TriangleT3::new(30, [3, 4, 5], 0, 300).expect("valid triangle should be created"),
                ),
                vec![Dof2D::Ux, Dof2D::Uy],
                6,
            ),
            (
                "triangle_t6",
                Element2D::TriangleT6(
                    TriangleT6::new(35, [3, 4, 5, 6, 7, 8], 0, 350).expect("valid triangle should be created"),
                ),
                vec![Dof2D::Ux, Dof2D::Uy],
                12,
            ),
            (
                "quad",
                Element2D::QuadQ4(QuadQ4::new(40, [4, 5, 6, 7], 0, 400).expect("valid quad should be created")),
                vec![Dof2D::Ux, Dof2D::Uy],
                8,
            ),
            (
                "quad_q8",
                Element2D::QuadQ8(
                    QuadQ8::new(50, [4, 5, 6, 7, 8, 9, 10, 11], 0, 500).expect("valid quad should be created"),
                ),
                vec![Dof2D::Ux, Dof2D::Uy],
                16,
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
            Node2D::new(40, 1.0, 1.0).expect("valid node should be created"),
            Node2D::new(50, 0.5, 0.0).expect("valid node should be created"),
            Node2D::new(60, 1.0, 0.5).expect("valid node should be created"),
            Node2D::new(70, 0.5, 1.0).expect("valid node should be created"),
            Node2D::new(80, 0.0, 0.5).expect("valid node should be created"),
        ];
        let cases = [
            (
                "truss",
                Element2D::Truss(Truss2D::new(1, [10, 20], 0, 100).expect("valid truss should be created")),
                Section2D::Truss(TrussSection2D::new(1.0).expect("valid section")),
            ),
            (
                "beam",
                Element2D::Beam(Beam2D::new(2, [10, 20], 0, 200).expect("valid beam should be created")),
                Section2D::Beam(BeamSection2D::new(1.0, 1.0).expect("valid section")),
            ),
            (
                "triangle",
                Element2D::TriangleT3(
                    TriangleT3::new(3, [10, 20, 30], 0, 300).expect("valid triangle should be created"),
                ),
                Section2D::PlaneStress(PlaneStressSection2D::new(1.0).expect("valid section")),
            ),
            (
                "triangle_t6",
                Element2D::TriangleT6(
                    TriangleT6::new(6, [10, 20, 30, 50, 60, 80], 0, 600).expect("valid triangle should be created"),
                ),
                Section2D::PlaneStress(PlaneStressSection2D::new(1.0).expect("valid section")),
            ),
            (
                "quad",
                Element2D::QuadQ4(QuadQ4::new(4, [10, 20, 40, 30], 0, 400).expect("valid quad should be created")),
                Section2D::PlaneStress(PlaneStressSection2D::new(1.0).expect("valid section")),
            ),
            (
                "quad_q8",
                Element2D::QuadQ8(
                    QuadQ8::new(5, [10, 20, 40, 30, 50, 60, 70, 80], 0, 500).expect("valid quad should be created"),
                ),
                Section2D::PlaneStress(PlaneStressSection2D::new(1.0).expect("valid section")),
            ),
        ];

        for (name, element, section) in cases {
            let actual =
                element.stiffness_matrix(&material, &section, &nodes).expect("stiffness matrix should be calculated");
            let expected = match (&element, &section) {
                (Element2D::Truss(truss), Section2D::Truss(section)) => dynamic_matrix_from_array(
                    truss
                        .stiffness_matrix(&material, section, &nodes[0], &nodes[1])
                        .expect("stiffness matrix should be calculated"),
                ),
                (Element2D::Beam(beam), Section2D::Beam(section)) => dynamic_matrix_from_array(
                    beam.stiffness_matrix(&material, section, &nodes[0], &nodes[1])
                        .expect("stiffness matrix should be calculated"),
                ),
                (Element2D::TriangleT3(triangle), Section2D::PlaneStress(section)) => dynamic_matrix_from_array(
                    triangle
                        .stiffness_matrix(&material, section, &nodes[0], &nodes[1], &nodes[2])
                        .expect("stiffness matrix should be calculated"),
                ),
                (Element2D::TriangleT6(triangle), Section2D::PlaneStress(section)) => dynamic_matrix_from_array(
                    triangle
                        .stiffness_matrix(
                            &material,
                            section,
                            [&nodes[0], &nodes[1], &nodes[2], &nodes[4], &nodes[5], &nodes[7]],
                        )
                        .expect("stiffness matrix should be calculated"),
                ),
                (Element2D::QuadQ4(quad), Section2D::PlaneStress(section)) => dynamic_matrix_from_array(
                    quad.stiffness_matrix(&material, section, &nodes[0], &nodes[1], &nodes[3], &nodes[2])
                        .expect("stiffness matrix should be calculated"),
                ),
                (Element2D::QuadQ8(quad), Section2D::PlaneStress(section)) => dynamic_matrix_from_array(
                    quad.stiffness_matrix(
                        &material,
                        section,
                        [&nodes[0], &nodes[1], &nodes[3], &nodes[2], &nodes[4], &nodes[5], &nodes[6], &nodes[7]],
                    )
                    .expect("stiffness matrix should be calculated"),
                ),
                _ => unreachable!("test cases should pair compatible elements and sections"),
            };

            assert_eq!(actual, expected, "failed case: {name}");
        }
    }

    #[test]
    fn common_stiffness_matrix_rejects_unknown_node() {
        let material = Material2D::new(200.0, 0.3, 7800.0).expect("valid material should be created");
        let known_node = Node2D::new(10, 0.0, 0.0).expect("valid node should be created");
        let element =
            Element2D::Truss(Truss2D::new(1, [10, 20], 0, 100).expect("valid truss connectivity should be created"));
        let section = Section2D::Truss(TrussSection2D::new(1.0).expect("valid section"));

        let result = element.stiffness_matrix(&material, &section, &[known_node]);

        assert!(matches!(result, Err(FemError::UnknownId { entity: "node", id: 20 })));
    }

    #[test]
    fn calculates_plane_stress_stiffness_matrix_for_right_triangle() {
        let material = Material2D::new(1.0, 0.3, 1.0).expect("valid material should be created");
        let first_node = Node2D::new(1, 0.0, 0.0).expect("valid node should be created");
        let second_node = Node2D::new(2, 1.0, 0.0).expect("valid node should be created");
        let third_node = Node2D::new(3, 0.0, 1.0).expect("valid node should be created");
        let triangle = TriangleT3::new(30, [1, 2, 3], 0, 300).expect("valid triangle should be created");
        let section = PlaneStressSection2D::new(1.0).expect("valid section should be created");

        let matrix = triangle
            .stiffness_matrix(&material, &section, &first_node, &second_node, &third_node)
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
        let triangle = TriangleT3::new(30, [1, 2, 3], 0, 300).expect("valid triangle should be created");
        let section = PlaneStressSection2D::new(0.2).expect("valid section should be created");

        let matrix = triangle
            .stiffness_matrix(&material, &section, &first_node, &second_node, &third_node)
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
    fn triangle_t6_stiffness_matrix_is_symmetric_and_balanced() {
        let material = Material2D::new(200.0, 0.3, 7800.0).expect("valid material should be created");
        let first_node = Node2D::new(1, 0.0, 0.0).expect("valid node should be created");
        let second_node = Node2D::new(2, 2.0, 0.0).expect("valid node should be created");
        let third_node = Node2D::new(3, 0.0, 1.0).expect("valid node should be created");
        let fourth_node = Node2D::new(4, 1.0, 0.0).expect("valid node should be created");
        let fifth_node = Node2D::new(5, 1.0, 0.5).expect("valid node should be created");
        let sixth_node = Node2D::new(6, 0.0, 0.5).expect("valid node should be created");
        let triangle = TriangleT6::new(35, [1, 2, 3, 4, 5, 6], 0, 350).expect("valid triangle should be created");
        let section = PlaneStressSection2D::new(0.2).expect("valid section should be created");

        let matrix = triangle
            .stiffness_matrix(
                &material,
                &section,
                [&first_node, &second_node, &third_node, &fourth_node, &fifth_node, &sixth_node],
            )
            .expect("stiffness matrix should be calculated");

        for (row_index, row) in matrix.iter().enumerate() {
            for (column_index, value) in row.iter().enumerate() {
                assert!(
                    (*value - matrix[column_index][row_index]).abs() < 1e-10,
                    "matrix is not symmetric at row {row_index}, column {column_index}"
                );
            }
        }

        for (row_index, row) in matrix.iter().enumerate() {
            let sum = row.iter().sum::<f64>();

            assert!(sum.abs() < 1e-10, "row {row_index} is not balanced: {sum}");
        }
    }

    #[test]
    fn quad_q4_stiffness_matrix_is_symmetric_and_balanced() {
        let material = Material2D::new(200.0, 0.3, 7800.0).expect("valid material should be created");
        let first_node = Node2D::new(1, 0.0, 0.0).expect("valid node should be created");
        let second_node = Node2D::new(2, 2.0, 0.0).expect("valid node should be created");
        let third_node = Node2D::new(3, 2.0, 1.0).expect("valid node should be created");
        let fourth_node = Node2D::new(4, 0.0, 1.0).expect("valid node should be created");
        let quad = QuadQ4::new(40, [1, 2, 3, 4], 0, 400).expect("valid quad should be created");
        let section = PlaneStressSection2D::new(0.2).expect("valid section should be created");

        let matrix = quad
            .stiffness_matrix(&material, &section, &first_node, &second_node, &third_node, &fourth_node)
            .expect("stiffness matrix should be calculated");

        for (row_index, row) in matrix.iter().enumerate() {
            for (column_index, value) in row.iter().enumerate() {
                assert!(
                    (*value - matrix[column_index][row_index]).abs() < 1e-10,
                    "matrix is not symmetric at row {row_index}, column {column_index}"
                );
            }
        }

        for (row_index, row) in matrix.iter().enumerate() {
            let sum = row.iter().sum::<f64>();

            assert!(sum.abs() < 1e-10, "row {row_index} is not balanced: {sum}");
        }
    }

    #[test]
    fn quad_q8_stiffness_matrix_is_symmetric_and_balanced() {
        let material = Material2D::new(200.0, 0.3, 7800.0).expect("valid material should be created");
        let first_node = Node2D::new(1, 0.0, 0.0).expect("valid node should be created");
        let second_node = Node2D::new(2, 2.0, 0.0).expect("valid node should be created");
        let third_node = Node2D::new(3, 2.0, 1.0).expect("valid node should be created");
        let fourth_node = Node2D::new(4, 0.0, 1.0).expect("valid node should be created");
        let fifth_node = Node2D::new(5, 1.0, 0.0).expect("valid node should be created");
        let sixth_node = Node2D::new(6, 2.0, 0.5).expect("valid node should be created");
        let seventh_node = Node2D::new(7, 1.0, 1.0).expect("valid node should be created");
        let eighth_node = Node2D::new(8, 0.0, 0.5).expect("valid node should be created");
        let quad = QuadQ8::new(50, [1, 2, 3, 4, 5, 6, 7, 8], 0, 500).expect("valid quad should be created");
        let section = PlaneStressSection2D::new(0.2).expect("valid section should be created");

        let matrix = quad
            .stiffness_matrix(
                &material,
                &section,
                [
                    &first_node,
                    &second_node,
                    &third_node,
                    &fourth_node,
                    &fifth_node,
                    &sixth_node,
                    &seventh_node,
                    &eighth_node,
                ],
            )
            .expect("stiffness matrix should be calculated");

        for (row_index, row) in matrix.iter().enumerate() {
            for (column_index, value) in row.iter().enumerate() {
                assert!(
                    (*value - matrix[column_index][row_index]).abs() < 1e-10,
                    "matrix is not symmetric at row {row_index}, column {column_index}"
                );
            }
        }

        for (row_index, row) in matrix.iter().enumerate() {
            let sum = row.iter().sum::<f64>();

            assert!(sum.abs() < 1e-9, "row {row_index} is not balanced: {sum}");
        }
    }

    #[test]
    fn rejects_stiffness_matrix_for_inverted_quad_q4() {
        let material = Material2D::new(200.0, 0.3, 7800.0).expect("valid material should be created");
        let first_node = Node2D::new(1, 0.0, 0.0).expect("valid node should be created");
        let second_node = Node2D::new(2, 0.0, 1.0).expect("valid node should be created");
        let third_node = Node2D::new(3, 1.0, 1.0).expect("valid node should be created");
        let fourth_node = Node2D::new(4, 1.0, 0.0).expect("valid node should be created");
        let quad = QuadQ4::new(40, [1, 2, 3, 4], 0, 400).expect("valid connectivity should be created");
        let section = PlaneStressSection2D::new(0.2).expect("valid section should be created");

        let result = quad.stiffness_matrix(&material, &section, &first_node, &second_node, &third_node, &fourth_node);

        assert!(matches!(
            result,
            Err(FemError::DegenerateElement {
                element_id: 40,
                element_type: "quad_q4",
                measure_name: "jacobian determinant",
                ..
            })
        ));
    }

    #[test]
    fn rejects_stiffness_matrix_for_inverted_quad_q8() {
        let material = Material2D::new(200.0, 0.3, 7800.0).expect("valid material should be created");
        let first_node = Node2D::new(1, 0.0, 0.0).expect("valid node should be created");
        let second_node = Node2D::new(2, 0.0, 1.0).expect("valid node should be created");
        let third_node = Node2D::new(3, 1.0, 1.0).expect("valid node should be created");
        let fourth_node = Node2D::new(4, 1.0, 0.0).expect("valid node should be created");
        let fifth_node = Node2D::new(5, 0.0, 0.5).expect("valid node should be created");
        let sixth_node = Node2D::new(6, 0.5, 1.0).expect("valid node should be created");
        let seventh_node = Node2D::new(7, 1.0, 0.5).expect("valid node should be created");
        let eighth_node = Node2D::new(8, 0.5, 0.0).expect("valid node should be created");
        let quad = QuadQ8::new(50, [1, 2, 3, 4, 5, 6, 7, 8], 0, 500).expect("valid connectivity should be created");
        let section = PlaneStressSection2D::new(0.2).expect("valid section should be created");

        let result = quad.stiffness_matrix(
            &material,
            &section,
            [
                &first_node,
                &second_node,
                &third_node,
                &fourth_node,
                &fifth_node,
                &sixth_node,
                &seventh_node,
                &eighth_node,
            ],
        );

        assert!(matches!(
            result,
            Err(FemError::DegenerateElement {
                element_id: 50,
                element_type: "quad_q8",
                measure_name: "jacobian determinant",
                ..
            })
        ));
    }

    #[test]
    fn rejects_stiffness_matrix_for_zero_area_triangle() {
        let material = Material2D::new(200.0, 0.3, 7800.0).expect("valid material should be created");
        let first_node = Node2D::new(1, 0.0, 0.0).expect("valid node should be created");
        let second_node = Node2D::new(2, 1.0, 1.0).expect("valid node should be created");
        let third_node = Node2D::new(3, 2.0, 2.0).expect("valid node should be created");
        let triangle = TriangleT3::new(30, [1, 2, 3], 0, 300).expect("valid connectivity should be created");
        let section = PlaneStressSection2D::new(0.2).expect("valid section should be created");

        let result = triangle.stiffness_matrix(&material, &section, &first_node, &second_node, &third_node);

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
    fn calculates_stiffness_matrix_for_horizontal_truss() {
        let material = Material2D::new(200.0, 0.3, 7800.0).expect("valid material should be created");
        let first_node = Node2D::new(1, 0.0, 0.0).expect("valid node should be created");
        let second_node = Node2D::new(2, 2.0, 0.0).expect("valid node should be created");
        let truss = Truss2D::new(10, [1, 2], 0, 100).expect("valid truss should be created");
        let section = TrussSection2D::new(0.5).expect("valid section should be created");

        let matrix = truss
            .stiffness_matrix(&material, &section, &first_node, &second_node)
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
        let beam = Beam2D::new(20, [1, 2], 0, 200).expect("valid beam should be created");
        let section = BeamSection2D::new(2.0, 3.0).expect("valid section should be created");

        let matrix = beam
            .stiffness_matrix(&material, &section, &first_node, &second_node)
            .expect("stiffness matrix should be calculated");

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
        let beam = Beam2D::new(20, [1, 2], 0, 200).expect("valid beam should be created");
        let section = BeamSection2D::new(2.0, 3.0).expect("valid section should be created");

        let matrix = beam
            .stiffness_matrix(&material, &section, &first_node, &second_node)
            .expect("stiffness matrix should be calculated");

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
        let beam = Beam2D::new(20, [1, 2], 0, 200).expect("valid beam should be created");
        let section = BeamSection2D::new(1.0, 1.0).expect("valid section should be created");

        let matrix = beam
            .stiffness_matrix(&material, &section, &first_node, &second_node)
            .expect("stiffness matrix should be calculated");

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
        let beam = Beam2D::new(20, [1, 2], 0, 200).expect("valid beam should be created");
        let section = BeamSection2D::new(2.0, 3.0).expect("valid section should be created");

        let matrix = beam
            .stiffness_matrix(&material, &section, &first_node, &second_node)
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
    fn rejects_stiffness_matrix_for_zero_length_beam() {
        let material = Material2D::new(200.0, 0.3, 7800.0).expect("valid material should be created");
        let first_node = Node2D::new(1, 0.0, 0.0).expect("valid node should be created");
        let second_node = Node2D::new(2, 0.0, 0.0).expect("valid node should be created");
        let beam = Beam2D::new(20, [1, 2], 0, 200).expect("valid connectivity should be created");
        let section = BeamSection2D::new(2.0, 3.0).expect("valid section should be created");

        let result = beam.stiffness_matrix(&material, &section, &first_node, &second_node);

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
        let truss = Truss2D::new(10, [1, 2], 0, 100).expect("valid truss should be created");
        let section = TrussSection2D::new(0.5).expect("valid section should be created");

        let matrix = truss
            .stiffness_matrix(&material, &section, &first_node, &second_node)
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
        let truss = Truss2D::new(10, [1, 2], 0, 100).expect("valid truss should be created");
        let section = TrussSection2D::new(0.5).expect("valid section should be created");

        let matrix = truss
            .stiffness_matrix(&material, &section, &first_node, &second_node)
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
        let truss = Truss2D::new(10, [1, 2], 0, 100).expect("valid connectivity should be created");
        let section = TrussSection2D::new(0.5).expect("valid section should be created");

        let result = truss.stiffness_matrix(&material, &section, &first_node, &second_node);

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
