//! Postprocessing of solved FEM results.

use nalgebra::DVector;

use crate::elements::interpolation::{
    cubic_hermite_first_derivatives, cubic_hermite_second_derivatives, cubic_hermite_shape_functions,
    linear_lagrange_shape_functions, triangle_t3_shape_functions,
};
use crate::elements::{Beam2D, Element2D, TriangleT3, Truss2D};
use crate::error::FemError;
use crate::model::{BeamSection2D, DofNumbering2D, Material2D, Model2D, Node2D, TrussSection2D};

/// Position used when interpolating a displacement field inside a 2D element.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ElementPosition2D {
    /// Physical distance from the first node of a two-node element.
    Line { position: f64 },

    /// Natural coordinates inside a T3 triangle.
    Triangle { xi: f64, eta: f64 },
}

/// Stress and strain results recovered for a 2D truss element.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TrussResponse2D {
    strain: f64,
    stress: f64,
    axial_force: f64,
}

/// Interpolated translational displacement inside a 2D truss element.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TrussDisplacement2D {
    ux: f64,
    uy: f64,
}

impl TrussDisplacement2D {
    /// Returns the interpolated x-displacement.
    #[must_use]
    pub fn ux(&self) -> f64 {
        self.ux
    }

    /// Returns the interpolated y-displacement.
    #[must_use]
    pub fn uy(&self) -> f64 {
        self.uy
    }
}

/// Interpolates a truss displacement at a physical position measured from its first node.
pub fn interpolate_truss_displacement(
    model: &Model2D, truss: &Truss2D, global_displacements: &DVector<f64>, position: f64,
) -> Result<TrussDisplacement2D, FemError> {
    let numbering = DofNumbering2D::from_model(model)?;

    if global_displacements.len() != numbering.count() {
        return Err(FemError::InvalidDisplacementVector {
            expected: numbering.count(),
            actual: global_displacements.len(),
        });
    }

    let element = Element2D::Truss(*truss);
    let node_ids = element.node_ids();
    let first_node = find_node(model.nodes(), node_ids[0])?;
    let second_node = find_node(model.nodes(), node_ids[1])?;
    let length = truss.length(first_node, second_node)?;
    let xi = normalized_two_node_position(position, length)?;
    let shape_functions = linear_lagrange_shape_functions(xi);
    let indices = numbering.element_dof_indices(&element)?;
    let first_ux = global_displacements[indices[0]];
    let first_uy = global_displacements[indices[1]];
    let second_ux = global_displacements[indices[2]];
    let second_uy = global_displacements[indices[3]];

    Ok(TrussDisplacement2D {
        ux: shape_functions[0] * first_ux + shape_functions[1] * second_ux,
        uy: shape_functions[0] * first_uy + shape_functions[1] * second_uy,
    })
}

impl TrussResponse2D {
    /// Returns the axial strain of the truss element.
    #[must_use]
    pub fn strain(&self) -> f64 {
        self.strain
    }

    /// Returns the normal stress of the truss element.
    #[must_use]
    pub fn stress(&self) -> f64 {
        self.stress
    }

    /// Returns the signed axial force of the truss element.
    #[must_use]
    pub fn axial_force(&self) -> f64 {
        self.axial_force
    }
}

/// Recovers the axial response of a truss from the global displacement vector.
pub fn recover_truss_response(
    model: &Model2D, truss: &Truss2D, global_displacements: &DVector<f64>,
) -> Result<TrussResponse2D, FemError> {
    let material = model.material(truss.material_id())?;
    let section = model.truss_section(truss.section_id())?;
    let numbering = DofNumbering2D::from_model(model)?;

    if global_displacements.len() != numbering.count() {
        return Err(FemError::InvalidDisplacementVector {
            expected: numbering.count(),
            actual: global_displacements.len(),
        });
    }

    let element = Element2D::Truss(*truss);
    let node_ids = element.node_ids();
    let first_node = find_node(model.nodes(), node_ids[0])?;
    let second_node = find_node(model.nodes(), node_ids[1])?;
    let indices = numbering.element_dof_indices(&element)?;

    let element_displacements = [
        global_displacements[indices[0]],
        global_displacements[indices[1]],
        global_displacements[indices[2]],
        global_displacements[indices[3]],
    ];

    calculate_truss_response(truss, material, section, first_node, second_node, element_displacements)
}

fn find_node(nodes: &[Node2D], node_id: usize) -> Result<&Node2D, FemError> {
    nodes.iter().find(|node| node.id() == node_id).ok_or(FemError::UnknownId { entity: "node", id: node_id })
}

fn calculate_truss_response(
    truss: &Truss2D, material: &Material2D, section: &TrussSection2D, first_node: &Node2D, second_node: &Node2D,
    [first_x, first_y, second_x, second_y]: [f64; 4],
) -> Result<TrussResponse2D, FemError> {
    let dx = second_node.x() - first_node.x();
    let dy = second_node.y() - first_node.y();
    let length = (dx * dx + dy * dy).sqrt();

    if !length.is_finite() || length == 0.0 {
        return Err(FemError::DegenerateElement {
            element_id: truss_element_id(truss),
            element_type: "truss",
            node_ids: vec![first_node.id(), second_node.id()],
            measure_name: "length",
            measure: length,
        });
    }

    let cosine = dx / length;
    let sine = dy / length;
    let strain = (-cosine * first_x - sine * first_y + cosine * second_x + sine * second_y) / length;
    let stress = material.young_modulus() * strain;
    let axial_force = stress * section.cross_section_area();

    Ok(TrussResponse2D { strain, stress, axial_force })
}

fn truss_element_id(truss: &Truss2D) -> usize {
    Element2D::Truss(*truss).id()
}

/// Signed local end forces recovered for a 2D beam element.
///
/// The entries are ordered as `[N1, V1, M1, N2, V2, M2]`, where `N` is the
/// axial force, `V` is the shear force, and `M` is the bending moment.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BeamResponse2D {
    end_forces: [f64; 6],
}

impl BeamResponse2D {
    /// Returns the local end-force vector `[N1, V1, M1, N2, V2, M2]`.
    #[must_use]
    pub fn end_forces(&self) -> &[f64; 6] {
        &self.end_forces
    }

    /// Returns the axial force at the first node.
    #[must_use]
    pub fn first_axial_force(&self) -> f64 {
        self.end_forces[0]
    }

    /// Returns the shear force at the first node.
    #[must_use]
    pub fn first_shear_force(&self) -> f64 {
        self.end_forces[1]
    }

    /// Returns the bending moment at the first node.
    #[must_use]
    pub fn first_bending_moment(&self) -> f64 {
        self.end_forces[2]
    }

    /// Returns the axial force at the second node.
    #[must_use]
    pub fn second_axial_force(&self) -> f64 {
        self.end_forces[3]
    }

    /// Returns the shear force at the second node.
    #[must_use]
    pub fn second_shear_force(&self) -> f64 {
        self.end_forces[4]
    }

    /// Returns the bending moment at the second node.
    #[must_use]
    pub fn second_bending_moment(&self) -> f64 {
        self.end_forces[5]
    }
}

/// Recovers the local end forces of a beam from the global displacement vector.
pub fn recover_beam_response(
    model: &Model2D, beam: &Beam2D, global_displacements: &DVector<f64>,
) -> Result<BeamResponse2D, FemError> {
    let material = model.material(beam.material_id())?;
    let section = model.beam_section(beam.section_id())?;
    let (local_displacements, length) = extract_beam_local_displacements(model, beam, global_displacements)?;

    Ok(calculate_beam_end_forces(beam, material, section, length, local_displacements))
}

fn extract_beam_local_displacements(
    model: &Model2D, beam: &Beam2D, global_displacements: &DVector<f64>,
) -> Result<([f64; 6], f64), FemError> {
    let numbering = DofNumbering2D::from_model(model)?;

    if global_displacements.len() != numbering.count() {
        return Err(FemError::InvalidDisplacementVector {
            expected: numbering.count(),
            actual: global_displacements.len(),
        });
    }

    let element = Element2D::Beam(*beam);
    let node_ids = element.node_ids();
    let first_node = find_node(model.nodes(), node_ids[0])?;
    let second_node = find_node(model.nodes(), node_ids[1])?;
    let indices = numbering.element_dof_indices(&element)?;
    let global_displacements = [
        global_displacements[indices[0]],
        global_displacements[indices[1]],
        global_displacements[indices[2]],
        global_displacements[indices[3]],
        global_displacements[indices[4]],
        global_displacements[indices[5]],
    ];
    let (length, cosine, sine) = beam.geometry(first_node, second_node)?;
    let local_displacements = [
        cosine * global_displacements[0] + sine * global_displacements[1],
        -sine * global_displacements[0] + cosine * global_displacements[1],
        global_displacements[2],
        cosine * global_displacements[3] + sine * global_displacements[4],
        -sine * global_displacements[3] + cosine * global_displacements[4],
        global_displacements[5],
    ];

    Ok((local_displacements, length))
}

fn calculate_beam_end_forces(
    beam: &Beam2D, material: &Material2D, section: &BeamSection2D, length: f64, local_displacements: [f64; 6],
) -> BeamResponse2D {
    let local_stiffness_matrix = beam.local_stiffness_matrix(material, section, length);
    let mut end_forces = [0.0; 6];

    for (row, end_force) in end_forces.iter_mut().enumerate() {
        *end_force = local_stiffness_matrix[row]
            .iter()
            .zip(local_displacements)
            .map(|(stiffness, displacement)| stiffness * displacement)
            .sum();
    }

    BeamResponse2D { end_forces }
}

/// Interpolated local beam displacement at a physical position measured from
/// the first node.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BeamDisplacement2D {
    axial: f64,
    transverse: f64,
    rotation: f64,
}

impl BeamDisplacement2D {
    /// Returns the axial displacement in the local beam coordinate system.
    #[must_use]
    pub fn axial(&self) -> f64 {
        self.axial
    }

    /// Returns the transverse displacement in the local beam coordinate system.
    #[must_use]
    pub fn transverse(&self) -> f64 {
        self.transverse
    }

    /// Returns the cross-section rotation in the local beam coordinate system.
    #[must_use]
    pub fn rotation(&self) -> f64 {
        self.rotation
    }
}

/// Section-level bending results at a point inside a beam element.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BeamSectionResponse2D {
    axial_strain: f64,
    axial_force: f64,
    curvature: f64,
    bending_moment: f64,
    cross_section_area: f64,
    second_moment_of_area: f64,
}

impl BeamSectionResponse2D {
    /// Returns the axial strain in the beam section.
    #[must_use]
    pub fn axial_strain(&self) -> f64 {
        self.axial_strain
    }

    /// Returns the signed axial force in the beam section.
    #[must_use]
    pub fn axial_force(&self) -> f64 {
        self.axial_force
    }

    /// Returns the curvature `d2v/dx2` in the local beam coordinate system.
    #[must_use]
    pub fn curvature(&self) -> f64 {
        self.curvature
    }

    /// Returns the bending moment calculated as `E * I * curvature`.
    #[must_use]
    pub fn bending_moment(&self) -> f64 {
        self.bending_moment
    }

    /// Returns normal stress at distance `y` from the section's neutral axis.
    ///
    /// The sign convention is `sigma(y) = N/A - M*y/I`.
    #[must_use]
    pub fn normal_stress(&self, y: f64) -> f64 {
        self.axial_force / self.cross_section_area - self.bending_moment * y / self.second_moment_of_area
    }
}

/// Recovers curvature and bending moment at a physical position inside a beam.
pub fn recover_beam_section_response(
    model: &Model2D, beam: &Beam2D, global_displacements: &DVector<f64>, position: f64,
) -> Result<BeamSectionResponse2D, FemError> {
    let material = model.material(beam.material_id())?;
    let section = model.beam_section(beam.section_id())?;
    let (local_displacements, length) = extract_beam_local_displacements(model, beam, global_displacements)?;
    let axial_strain = (local_displacements[3] - local_displacements[0]) / length;
    let axial_force = material.young_modulus() * section.cross_section_area() * axial_strain;
    let xi = normalized_two_node_position(position, length)?;
    let second_derivatives = cubic_hermite_second_derivatives(xi, length)?;
    let curvature = second_derivatives[0] * local_displacements[1]
        + second_derivatives[1] * local_displacements[2]
        + second_derivatives[2] * local_displacements[4]
        + second_derivatives[3] * local_displacements[5];
    let bending_moment = material.young_modulus() * section.second_moment_of_area() * curvature;

    Ok(BeamSectionResponse2D {
        axial_strain,
        axial_force,
        curvature,
        bending_moment,
        cross_section_area: section.cross_section_area(),
        second_moment_of_area: section.second_moment_of_area(),
    })
}

/// Interpolates a beam's local displacement field at `position`.
pub fn interpolate_beam_displacement(
    model: &Model2D, beam: &Beam2D, global_displacements: &DVector<f64>, position: f64,
) -> Result<BeamDisplacement2D, FemError> {
    let (local_displacements, length) = extract_beam_local_displacements(model, beam, global_displacements)?;
    let xi = normalized_two_node_position(position, length)?;
    let axial_shape_functions = linear_lagrange_shape_functions(xi);
    let transverse_shape_functions = cubic_hermite_shape_functions(xi, length)?;
    let transverse_derivatives = cubic_hermite_first_derivatives(xi, length)?;
    let axial = axial_shape_functions[0] * local_displacements[0] + axial_shape_functions[1] * local_displacements[3];
    let transverse = transverse_shape_functions[0] * local_displacements[1]
        + transverse_shape_functions[1] * local_displacements[2]
        + transverse_shape_functions[2] * local_displacements[4]
        + transverse_shape_functions[3] * local_displacements[5];
    let rotation = transverse_derivatives[0] * local_displacements[1]
        + transverse_derivatives[1] * local_displacements[2]
        + transverse_derivatives[2] * local_displacements[4]
        + transverse_derivatives[3] * local_displacements[5];

    Ok(BeamDisplacement2D { axial, transverse, rotation })
}

fn normalized_two_node_position(position: f64, length: f64) -> Result<f64, FemError> {
    if !position.is_finite() || position < 0.0 || position > length {
        return Err(FemError::InvalidInterpolationCoordinate { coordinate: position, length });
    }

    Ok(position / length)
}

/// Interpolated translational displacement inside a T3 triangle.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TriangleDisplacement2D {
    ux: f64,
    uy: f64,
}

/// Displacement result returned by the common interpolation interface.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ElementDisplacement2D {
    /// Interpolated displacement inside a truss.
    Truss(TrussDisplacement2D),

    /// Interpolated displacement inside a beam.
    Beam(BeamDisplacement2D),

    /// Interpolated displacement inside a T3 triangle.
    Triangle(TriangleDisplacement2D),
}

impl TriangleDisplacement2D {
    /// Returns the interpolated x-displacement.
    #[must_use]
    pub fn ux(&self) -> f64 {
        self.ux
    }

    /// Returns the interpolated y-displacement.
    #[must_use]
    pub fn uy(&self) -> f64 {
        self.uy
    }
}

/// Interpolates a T3 displacement at natural coordinates `(xi, eta)`.
pub fn interpolate_triangle_displacement(
    model: &Model2D, triangle: &TriangleT3, global_displacements: &DVector<f64>, xi: f64, eta: f64,
) -> Result<TriangleDisplacement2D, FemError> {
    let numbering = DofNumbering2D::from_model(model)?;

    if global_displacements.len() != numbering.count() {
        return Err(FemError::InvalidDisplacementVector {
            expected: numbering.count(),
            actual: global_displacements.len(),
        });
    }

    validate_triangle_natural_coordinates(xi, eta)?;

    let element = Element2D::TriangleT3(*triangle);
    let shape_functions = triangle_t3_shape_functions(xi, eta);
    let indices = numbering.element_dof_indices(&element)?;
    let mut ux = 0.0;
    let mut uy = 0.0;

    for (node_index, shape_function) in shape_functions.iter().enumerate() {
        ux += shape_function * global_displacements[indices[2 * node_index]];
        uy += shape_function * global_displacements[indices[2 * node_index + 1]];
    }

    Ok(TriangleDisplacement2D { ux, uy })
}

fn validate_triangle_natural_coordinates(xi: f64, eta: f64) -> Result<(), FemError> {
    if !xi.is_finite() || !eta.is_finite() || xi < 0.0 || eta < 0.0 || xi + eta > 1.0 {
        return Err(FemError::InvalidTriangleNaturalCoordinates { xi, eta });
    }

    Ok(())
}

/// Interpolates the displacement field of any supported 2D element.
pub fn interpolate_element_displacement(
    model: &Model2D, element: &Element2D, global_displacements: &DVector<f64>, position: ElementPosition2D,
) -> Result<ElementDisplacement2D, FemError> {
    match (element, position) {
        (Element2D::Truss(truss), ElementPosition2D::Line { position }) => {
            let displacement = interpolate_truss_displacement(model, truss, global_displacements, position)?;

            Ok(ElementDisplacement2D::Truss(displacement))
        }
        (Element2D::Beam(beam), ElementPosition2D::Line { position }) => {
            let displacement = interpolate_beam_displacement(model, beam, global_displacements, position)?;

            Ok(ElementDisplacement2D::Beam(displacement))
        }
        (Element2D::TriangleT3(triangle), ElementPosition2D::Triangle { xi, eta }) => {
            let displacement = interpolate_triangle_displacement(model, triangle, global_displacements, xi, eta)?;

            Ok(ElementDisplacement2D::Triangle(displacement))
        }
        (Element2D::Truss(_), ElementPosition2D::Triangle { .. }) => {
            Err(FemError::InvalidElementInterpolationPosition { element_type: "truss", expected: "a line position" })
        }
        (Element2D::Beam(_), ElementPosition2D::Triangle { .. }) => {
            Err(FemError::InvalidElementInterpolationPosition { element_type: "beam", expected: "a line position" })
        }
        (Element2D::TriangleT3(_), ElementPosition2D::Line { .. }) => {
            Err(FemError::InvalidElementInterpolationPosition {
                element_type: "triangle",
                expected: "triangle natural coordinates",
            })
        }
    }
}

/// Constant strain and plane-stress results recovered for a T3 triangle.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TriangleResponse2D {
    strain: [f64; 3],
    stress: [f64; 3],
}

impl TriangleResponse2D {
    /// Returns `[epsilon_x, epsilon_y, gamma_xy]`.
    #[must_use]
    pub fn strain(&self) -> &[f64; 3] {
        &self.strain
    }

    /// Returns `[sigma_x, sigma_y, tau_xy]`.
    #[must_use]
    pub fn stress(&self) -> &[f64; 3] {
        &self.stress
    }

    /// Returns the equivalent von Mises stress for plane stress.
    #[must_use]
    pub fn von_mises_stress(&self) -> f64 {
        let [sigma_x, sigma_y, tau_xy] = self.stress;

        (sigma_x.powi(2) - sigma_x * sigma_y + sigma_y.powi(2) + 3.0 * tau_xy.powi(2)).sqrt()
    }
}

/// Recovers the constant plane-stress response of a T3 triangle.
pub fn recover_triangle_response(
    model: &Model2D, triangle: &TriangleT3, global_displacements: &DVector<f64>,
) -> Result<TriangleResponse2D, FemError> {
    let material = model.material(triangle.material_id())?;
    model.plane_stress_section(triangle.section_id())?;
    let numbering = DofNumbering2D::from_model(model)?;

    if global_displacements.len() != numbering.count() {
        return Err(FemError::InvalidDisplacementVector {
            expected: numbering.count(),
            actual: global_displacements.len(),
        });
    }

    let element = Element2D::TriangleT3(*triangle);
    let node_ids = element.node_ids();
    let first_node = find_node(model.nodes(), node_ids[0])?;
    let second_node = find_node(model.nodes(), node_ids[1])?;
    let third_node = find_node(model.nodes(), node_ids[2])?;
    let indices = numbering.element_dof_indices(&element)?;
    let element_displacements = [
        global_displacements[indices[0]],
        global_displacements[indices[1]],
        global_displacements[indices[2]],
        global_displacements[indices[3]],
        global_displacements[indices[4]],
        global_displacements[indices[5]],
    ];

    calculate_triangle_response(triangle, material, first_node, second_node, third_node, element_displacements)
}

fn calculate_triangle_response(
    triangle: &TriangleT3, material: &Material2D, first_node: &Node2D, second_node: &Node2D, third_node: &Node2D,
    element_displacements: [f64; 6],
) -> Result<TriangleResponse2D, FemError> {
    let (strain_displacement_matrix, _) = triangle.strain_displacement_matrix(first_node, second_node, third_node)?;
    let constitutive_matrix = TriangleT3::constitutive_matrix(material);
    let mut strain = [0.0; 3];
    let mut stress = [0.0; 3];

    for (row, strain_value) in strain.iter_mut().enumerate() {
        *strain_value = strain_displacement_matrix[row]
            .iter()
            .zip(element_displacements)
            .map(|(coefficient, displacement)| coefficient * displacement)
            .sum();
    }

    for (row, stress_value) in stress.iter_mut().enumerate() {
        *stress_value = constitutive_matrix[row]
            .iter()
            .zip(strain)
            .map(|(coefficient, strain_value)| coefficient * strain_value)
            .sum();
    }

    Ok(TriangleResponse2D { strain, stress })
}

/// Postprocessing result for any supported 2D element.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ElementResponse2D {
    /// Recovered response of a truss element.
    Truss(TrussResponse2D),

    /// Recovered response of a beam element.
    Beam(BeamResponse2D),

    /// Recovered response of a T3 triangle.
    Triangle(TriangleResponse2D),
}

/// Recovers the response of one element using its concrete element type.
pub fn recover_element_response(
    model: &Model2D, element: &Element2D, global_displacements: &DVector<f64>,
) -> Result<ElementResponse2D, FemError> {
    match element {
        Element2D::Truss(truss) => {
            let response = recover_truss_response(model, truss, global_displacements)?;

            Ok(ElementResponse2D::Truss(response))
        }
        Element2D::Beam(beam) => {
            let response = recover_beam_response(model, beam, global_displacements)?;

            Ok(ElementResponse2D::Beam(response))
        }
        Element2D::TriangleT3(triangle) => {
            let response = recover_triangle_response(model, triangle, global_displacements)?;

            Ok(ElementResponse2D::Triangle(response))
        }
    }
}

/// Recovers responses for all elements in a model, preserving their IDs.
pub fn recover_model_responses(
    model: &Model2D, global_displacements: &DVector<f64>,
) -> Result<Vec<(usize, ElementResponse2D)>, FemError> {
    model
        .elements()
        .iter()
        .map(|element| {
            let response = recover_element_response(model, element, global_displacements)?;

            Ok((element.id(), response))
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::{
        ElementDisplacement2D, ElementPosition2D, ElementResponse2D, interpolate_beam_displacement,
        interpolate_element_displacement, interpolate_triangle_displacement, interpolate_truss_displacement,
        recover_beam_response, recover_beam_section_response, recover_element_response, recover_model_responses,
        recover_triangle_response, recover_truss_response,
    };
    use crate::analysis::solver::solve;
    use crate::elements::{Beam2D, Element2D, TriangleT3, Truss2D};
    use crate::error::FemError;
    use crate::model::{
        BeamSection2D, DEFAULT_MATERIAL_ID, DisplacementConstraint2D, Dof2D, Material2D, Model2D, NodalLoad2D, Node2D,
        PlaneStressSection2D, Section2D, TrussSection2D,
    };
    use approx::assert_relative_eq;
    use nalgebra::DVector;

    fn add_truss_element(model: &mut Model2D, id: usize, node_ids: [usize; 2], cross_section_area: f64) -> Truss2D {
        let truss = Truss2D::new(id, node_ids, DEFAULT_MATERIAL_ID, id).expect("valid truss should be created");
        let section =
            Section2D::Truss(TrussSection2D::new(cross_section_area).expect("valid truss section should be created"));

        model.add_element_with_section(Element2D::Truss(truss), section).expect("truss should be added");

        truss
    }

    fn add_beam_element(
        model: &mut Model2D, id: usize, node_ids: [usize; 2], cross_section_area: f64, second_moment_of_area: f64,
    ) -> Beam2D {
        let beam = Beam2D::new(id, node_ids, DEFAULT_MATERIAL_ID, id).expect("valid beam should be created");
        let section = Section2D::Beam(
            BeamSection2D::new(cross_section_area, second_moment_of_area)
                .expect("valid beam section should be created"),
        );

        model.add_element_with_section(Element2D::Beam(beam), section).expect("beam should be added");

        beam
    }

    fn add_triangle_element(model: &mut Model2D, id: usize, node_ids: [usize; 3], thickness: f64) -> TriangleT3 {
        let triangle =
            TriangleT3::new(id, node_ids, DEFAULT_MATERIAL_ID, id).expect("valid triangle should be created");
        let section = Section2D::PlaneStress(
            PlaneStressSection2D::new(thickness).expect("valid plane-stress section should be created"),
        );

        model.add_element_with_section(Element2D::TriangleT3(triangle), section).expect("triangle should be added");

        triangle
    }

    fn horizontal_truss_model() -> Model2D {
        let mut model = Model2D::new();

        model.set_material(Material2D::new(200.0, 0.3, 1.0).expect("valid material should be created"));
        model.add_node(Node2D::new(1, 0.0, 0.0).expect("valid node should be created")).expect("node should be added");
        model.add_node(Node2D::new(2, 1.0, 0.0).expect("valid node should be created")).expect("node should be added");
        add_truss_element(&mut model, 10, [1, 2], 2.0);

        model
    }

    #[test]
    fn recovers_tensile_response_for_horizontal_truss() {
        let model = horizontal_truss_model();
        let truss = match model.elements()[0] {
            Element2D::Truss(truss) => truss,
            _ => panic!("expected a truss element"),
        };
        let displacements = DVector::from_row_slice(&[0.0, 0.0, 0.025, 0.0]);

        let response = recover_truss_response(&model, &truss, &displacements).expect("response should be recovered");

        assert_relative_eq!(response.strain(), 0.025, epsilon = 1e-12);
        assert_relative_eq!(response.stress(), 5.0, epsilon = 1e-12);
        assert_relative_eq!(response.axial_force(), 10.0, epsilon = 1e-12);
    }

    #[test]
    fn interpolates_horizontal_truss_displacement() {
        let model = horizontal_truss_model();
        let truss = match model.elements()[0] {
            Element2D::Truss(truss) => truss,
            _ => panic!("expected a truss element"),
        };
        let displacements = DVector::from_row_slice(&[0.0, 1.0, 2.0, 3.0]);

        let actual = interpolate_truss_displacement(&model, &truss, &displacements, 0.25)
            .expect("truss displacement should be interpolated");

        assert_relative_eq!(actual.ux(), 0.5, epsilon = 1e-12);
        assert_relative_eq!(actual.uy(), 1.5, epsilon = 1e-12);
    }

    #[test]
    fn rejects_truss_interpolation_outside_element() {
        let model = horizontal_truss_model();
        let truss = match model.elements()[0] {
            Element2D::Truss(truss) => truss,
            _ => panic!("expected a truss element"),
        };
        let displacements = DVector::zeros(4);

        for position in [-0.1, 1.1, f64::NAN, f64::INFINITY] {
            let result = interpolate_truss_displacement(&model, &truss, &displacements, position);

            assert!(matches!(result, Err(FemError::InvalidInterpolationCoordinate { .. })));
        }
    }

    #[test]
    fn common_interpolation_dispatches_for_every_element_type() {
        let truss_model = horizontal_truss_model();
        let truss = &truss_model.elements()[0];
        let truss_displacement = interpolate_element_displacement(
            &truss_model,
            truss,
            &DVector::from_row_slice(&[0.0, 0.0, 2.0, 4.0]),
            ElementPosition2D::Line { position: 0.25 },
        )
        .expect("truss displacement should be interpolated");

        let ElementDisplacement2D::Truss(truss_displacement) = truss_displacement else {
            panic!("expected a truss displacement");
        };

        assert_relative_eq!(truss_displacement.ux(), 0.5, epsilon = 1e-12);
        assert_relative_eq!(truss_displacement.uy(), 1.0, epsilon = 1e-12);

        let beam_model = horizontal_beam_model();
        let beam = &beam_model.elements()[0];
        let beam_displacement = interpolate_element_displacement(
            &beam_model,
            beam,
            &DVector::from_row_slice(&[0.0, 0.0, 0.0, 0.1, 0.2, 0.3]),
            ElementPosition2D::Line { position: 0.5 },
        )
        .expect("beam displacement should be interpolated");

        let ElementDisplacement2D::Beam(beam_displacement) = beam_displacement else {
            panic!("expected a beam displacement");
        };

        assert_relative_eq!(beam_displacement.axial(), 0.05, epsilon = 1e-12);
        assert_relative_eq!(beam_displacement.transverse(), 0.0625, epsilon = 1e-12);

        let triangle_model = right_triangle_model([1, 2, 3]);
        let triangle = &triangle_model.elements()[0];
        let triangle_displacement = interpolate_element_displacement(
            &triangle_model,
            triangle,
            &DVector::from_row_slice(&[1.0, -1.0, 3.0, 3.0, 4.0, 4.0]),
            ElementPosition2D::Triangle { xi: 0.2, eta: 0.3 },
        )
        .expect("triangle displacement should be interpolated");

        let ElementDisplacement2D::Triangle(triangle_displacement) = triangle_displacement else {
            panic!("expected a triangle displacement");
        };

        assert_relative_eq!(triangle_displacement.ux(), 2.3, epsilon = 1e-12);
        assert_relative_eq!(triangle_displacement.uy(), 1.3, epsilon = 1e-12);
    }

    #[test]
    fn common_interpolation_rejects_position_for_wrong_element_shape() {
        let truss_model = horizontal_truss_model();
        let truss = &truss_model.elements()[0];
        let truss_result = interpolate_element_displacement(
            &truss_model,
            truss,
            &DVector::zeros(4),
            ElementPosition2D::Triangle { xi: 0.2, eta: 0.3 },
        );

        assert!(matches!(
            truss_result,
            Err(FemError::InvalidElementInterpolationPosition { element_type: "truss", expected: "a line position" })
        ));

        let triangle_model = right_triangle_model([1, 2, 3]);
        let triangle = &triangle_model.elements()[0];
        let triangle_result = interpolate_element_displacement(
            &triangle_model,
            triangle,
            &DVector::zeros(6),
            ElementPosition2D::Line { position: 0.2 },
        );

        assert!(matches!(
            triangle_result,
            Err(FemError::InvalidElementInterpolationPosition {
                element_type: "triangle",
                expected: "triangle natural coordinates"
            })
        ));
    }

    #[test]
    fn recovers_compressive_response_for_horizontal_truss() {
        let model = horizontal_truss_model();
        let truss = match model.elements()[0] {
            Element2D::Truss(truss) => truss,
            _ => panic!("expected a truss element"),
        };
        let displacements = DVector::from_row_slice(&[0.0, 0.0, -0.025, 0.0]);

        let response = recover_truss_response(&model, &truss, &displacements).expect("response should be recovered");

        assert_relative_eq!(response.strain(), -0.025, epsilon = 1e-12);
        assert_relative_eq!(response.stress(), -5.0, epsilon = 1e-12);
        assert_relative_eq!(response.axial_force(), -10.0, epsilon = 1e-12);
    }

    #[test]
    fn recovers_response_for_diagonal_truss() {
        let mut model = Model2D::new();

        model.set_material(Material2D::new(200.0, 0.3, 1.0).expect("valid material should be created"));
        model.add_node(Node2D::new(1, 0.0, 0.0).expect("valid node should be created")).expect("node should be added");
        model.add_node(Node2D::new(2, 1.0, 1.0).expect("valid node should be created")).expect("node should be added");
        add_truss_element(&mut model, 10, [1, 2], 1.0);

        let truss = match model.elements()[0] {
            Element2D::Truss(truss) => truss,
            _ => panic!("expected a truss element"),
        };
        let displacements = DVector::from_row_slice(&[0.0, 0.0, 0.1, 0.1]);

        let response = recover_truss_response(&model, &truss, &displacements).expect("response should be recovered");

        assert_relative_eq!(response.strain(), 0.1, epsilon = 1e-12);
        assert_relative_eq!(response.stress(), 20.0, epsilon = 1e-12);
        assert_relative_eq!(response.axial_force(), 20.0, epsilon = 1e-12);
    }

    #[test]
    fn rejects_displacement_vector_with_wrong_length() {
        let model = horizontal_truss_model();
        let truss = match model.elements()[0] {
            Element2D::Truss(truss) => truss,
            _ => panic!("expected a truss element"),
        };
        let displacements = DVector::from_row_slice(&[0.0, 0.0]);

        let result = recover_truss_response(&model, &truss, &displacements);

        assert!(matches!(result, Err(FemError::InvalidDisplacementVector { expected: 4, actual: 2 })));
    }

    #[test]
    fn recovers_response_from_solver_result() {
        let mut model = horizontal_truss_model();

        for (node_id, dof) in [(1, Dof2D::Ux), (1, Dof2D::Uy), (2, Dof2D::Uy)] {
            model
                .add_constraint(DisplacementConstraint2D::new(node_id, dof, 0.0).expect("valid constraint"))
                .expect("constraint should be added");
        }

        model.add_load(NodalLoad2D::new(2, Dof2D::Ux, 10.0).expect("valid load")).expect("load should be added");

        let result = solve(&model).expect("system should be solved");
        let truss = match model.elements()[0] {
            Element2D::Truss(truss) => truss,
            _ => panic!("expected a truss element"),
        };

        let response =
            recover_truss_response(&model, &truss, result.displacements()).expect("response should be recovered");

        assert_relative_eq!(response.strain(), 0.025, epsilon = 1e-12);
        assert_relative_eq!(response.stress(), 5.0, epsilon = 1e-12);
        assert_relative_eq!(response.axial_force(), 10.0, epsilon = 1e-12);
    }

    fn horizontal_beam_model() -> Model2D {
        let mut model = Model2D::new();

        model.set_material(Material2D::new(200.0, 0.3, 1.0).expect("valid material should be created"));
        model.add_node(Node2D::new(1, 0.0, 0.0).expect("valid node should be created")).expect("node should be added");
        model.add_node(Node2D::new(2, 1.0, 0.0).expect("valid node should be created")).expect("node should be added");
        add_beam_element(&mut model, 10, [1, 2], 1.0, 2.0);

        model
    }

    #[test]
    fn recovers_axial_force_for_horizontal_beam() {
        let model = horizontal_beam_model();
        let beam = match model.elements()[0] {
            Element2D::Beam(beam) => beam,
            _ => panic!("expected a beam element"),
        };
        let displacements = DVector::from_row_slice(&[0.0, 0.0, 0.0, 0.01, 0.0, 0.0]);

        let response = recover_beam_response(&model, &beam, &displacements).expect("response should be recovered");

        assert_relative_eq!(response.first_axial_force(), -2.0, epsilon = 1e-12);
        assert_relative_eq!(response.second_axial_force(), 2.0, epsilon = 1e-12);
        assert_relative_eq!(response.first_shear_force(), 0.0, epsilon = 1e-12);
        assert_relative_eq!(response.first_bending_moment(), 0.0, epsilon = 1e-12);

        let section = recover_beam_section_response(&model, &beam, &displacements, 0.5)
            .expect("section response should be recovered");

        assert_relative_eq!(section.axial_strain(), 0.01, epsilon = 1e-12);
        assert_relative_eq!(section.axial_force(), 2.0, epsilon = 1e-12);
        assert_relative_eq!(section.curvature(), 0.0, epsilon = 1e-12);
        assert_relative_eq!(section.bending_moment(), 0.0, epsilon = 1e-12);
        assert_relative_eq!(section.normal_stress(0.25), 2.0, epsilon = 1e-12);
    }

    #[test]
    fn recovers_axial_force_for_diagonal_beam() {
        let mut model = Model2D::new();

        model.set_material(Material2D::new(200.0, 0.3, 1.0).expect("valid material should be created"));
        model.add_node(Node2D::new(1, 0.0, 0.0).expect("valid node should be created")).expect("node should be added");
        model.add_node(Node2D::new(2, 1.0, 1.0).expect("valid node should be created")).expect("node should be added");
        add_beam_element(&mut model, 10, [1, 2], 1.0, 2.0);

        let beam = match model.elements()[0] {
            Element2D::Beam(beam) => beam,
            _ => panic!("expected a beam element"),
        };
        let cosine = 1.0 / 2.0_f64.sqrt();
        let sine = cosine;
        let displacement_along_beam = 0.01;
        let displacements = DVector::from_row_slice(&[
            0.0,
            0.0,
            0.0,
            cosine * displacement_along_beam,
            sine * displacement_along_beam,
            0.0,
        ]);

        let response = recover_beam_response(&model, &beam, &displacements).expect("response should be recovered");
        let expected_force = 200.0 * displacement_along_beam / 2.0_f64.sqrt();

        assert_relative_eq!(response.first_axial_force(), -expected_force, epsilon = 1e-12);
        assert_relative_eq!(response.second_axial_force(), expected_force, epsilon = 1e-12);
        assert_relative_eq!(response.first_shear_force(), 0.0, epsilon = 1e-12);
        assert_relative_eq!(response.first_bending_moment(), 0.0, epsilon = 1e-12);
    }

    #[test]
    fn recovers_cantilever_end_forces_from_solver_result() {
        let mut model = horizontal_beam_model();

        for dof in [Dof2D::Ux, Dof2D::Uy, Dof2D::Rz] {
            model
                .add_constraint(DisplacementConstraint2D::new(1, dof, 0.0).expect("valid constraint"))
                .expect("constraint should be added");
        }

        model.add_load(NodalLoad2D::new(2, Dof2D::Uy, -12.0).expect("valid load")).expect("load should be added");

        let result = solve(&model).expect("system should be solved");
        let beam = match model.elements()[0] {
            Element2D::Beam(beam) => beam,
            _ => panic!("expected a beam element"),
        };
        let response =
            recover_beam_response(&model, &beam, result.displacements()).expect("response should be recovered");

        assert_relative_eq!(response.first_axial_force(), 0.0, epsilon = 1e-12);
        assert_relative_eq!(response.first_shear_force(), 12.0, epsilon = 1e-12);
        assert_relative_eq!(response.first_bending_moment(), 12.0, epsilon = 1e-12);
        assert_relative_eq!(response.second_axial_force(), 0.0, epsilon = 1e-12);
        assert_relative_eq!(response.second_shear_force(), -12.0, epsilon = 1e-12);
        assert_relative_eq!(response.second_bending_moment(), 0.0, epsilon = 1e-12);
    }

    #[test]
    fn interpolates_cantilever_deflection_against_analytical_solution() {
        let mut model = horizontal_beam_model();

        for dof in [Dof2D::Ux, Dof2D::Uy, Dof2D::Rz] {
            model
                .add_constraint(DisplacementConstraint2D::new(1, dof, 0.0).expect("valid constraint"))
                .expect("constraint should be added");
        }

        let load = -12.0;
        let young_modulus = 200.0;
        let second_moment_of_area = 2.0;
        model.add_load(NodalLoad2D::new(2, Dof2D::Uy, load).expect("valid load")).expect("load should be added");

        let result = solve(&model).expect("system should be solved");
        let beam = match model.elements()[0] {
            Element2D::Beam(beam) => beam,
            _ => panic!("expected a beam element"),
        };

        for position in [0.0, 0.25, 0.5, 0.75, 1.0] {
            let actual = interpolate_beam_displacement(&model, &beam, result.displacements(), position)
                .expect("beam displacement should be interpolated");
            let expected_transverse =
                load * position.powi(2) * (3.0 - position) / (6.0 * young_modulus * second_moment_of_area);
            let expected_rotation = load * position * (2.0 - position) / (2.0 * young_modulus * second_moment_of_area);

            assert_relative_eq!(actual.axial(), 0.0, epsilon = 1e-12);
            assert_relative_eq!(actual.transverse(), expected_transverse, epsilon = 1e-12);
            assert_relative_eq!(actual.rotation(), expected_rotation, epsilon = 1e-12);
        }
    }

    #[test]
    fn recovers_cantilever_curvature_and_moment_against_analytical_solution() {
        let mut model = horizontal_beam_model();

        for dof in [Dof2D::Ux, Dof2D::Uy, Dof2D::Rz] {
            model
                .add_constraint(DisplacementConstraint2D::new(1, dof, 0.0).expect("valid constraint"))
                .expect("constraint should be added");
        }

        let load = -12.0;
        let young_modulus = 200.0;
        let second_moment_of_area = 2.0;
        model.add_load(NodalLoad2D::new(2, Dof2D::Uy, load).expect("valid load")).expect("load should be added");

        let result = solve(&model).expect("system should be solved");
        let beam = match model.elements()[0] {
            Element2D::Beam(beam) => beam,
            _ => panic!("expected a beam element"),
        };

        for position in [0.0, 0.25, 0.5, 0.75, 1.0] {
            let actual = recover_beam_section_response(&model, &beam, result.displacements(), position)
                .expect("section response should be recovered");
            let expected_curvature = load * (1.0 - position) / (young_modulus * second_moment_of_area);
            let expected_moment = load * (1.0 - position);

            assert_relative_eq!(actual.curvature(), expected_curvature, epsilon = 1e-12);
            assert_relative_eq!(actual.bending_moment(), expected_moment, epsilon = 1e-12);
            assert_relative_eq!(actual.axial_strain(), 0.0, epsilon = 1e-12);
            assert_relative_eq!(actual.axial_force(), 0.0, epsilon = 1e-12);
            assert_relative_eq!(actual.normal_stress(1.0), -expected_moment / second_moment_of_area, epsilon = 1e-12);
        }
    }

    fn simply_supported_two_element_beam_model() -> Model2D {
        let mut model = Model2D::new();

        model.set_material(Material2D::new(200.0, 0.3, 1.0).expect("valid material should be created"));
        for (id, x) in [(1, 0.0), (2, 1.0), (3, 2.0)] {
            model
                .add_node(Node2D::new(id, x, 0.0).expect("valid node should be created"))
                .expect("node should be added");
        }
        add_beam_element(&mut model, 10, [1, 2], 1.0, 2.0);
        add_beam_element(&mut model, 20, [2, 3], 1.0, 2.0);

        for (node_id, dof) in [(1, Dof2D::Ux), (1, Dof2D::Uy), (3, Dof2D::Uy)] {
            model
                .add_constraint(DisplacementConstraint2D::new(node_id, dof, 0.0).expect("valid constraint"))
                .expect("constraint should be added");
        }

        model.add_load(NodalLoad2D::new(2, Dof2D::Uy, -12.0).expect("valid load")).expect("load should be added");

        model
    }

    #[test]
    fn interpolates_simply_supported_beam_against_analytical_solution() {
        let model = simply_supported_two_element_beam_model();
        let result = solve(&model).expect("system should be solved");
        let first_beam = match model.elements()[0] {
            Element2D::Beam(beam) => beam,
            _ => panic!("expected a beam element"),
        };
        let second_beam = match model.elements()[1] {
            Element2D::Beam(beam) => beam,
            _ => panic!("expected a beam element"),
        };
        let load_magnitude = 12.0;
        let young_modulus = 200.0;
        let second_moment_of_area = 2.0;
        let span = 2.0;

        for position in [0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0] {
            let (beam, local_position) =
                if position <= 1.0 { (&first_beam, position) } else { (&second_beam, position - 1.0) };
            let actual = interpolate_beam_displacement(&model, beam, result.displacements(), local_position)
                .expect("beam displacement should be interpolated");
            let distance_from_support = position.min(span - position);
            let expected_transverse =
                -load_magnitude * distance_from_support * (3.0 * span.powi(2) - 4.0 * distance_from_support.powi(2))
                    / (48.0 * young_modulus * second_moment_of_area);

            assert_relative_eq!(actual.transverse(), expected_transverse, epsilon = 1e-12);
        }
    }

    #[test]
    fn recovers_simply_supported_beam_moment_against_analytical_solution() {
        let model = simply_supported_two_element_beam_model();
        let result = solve(&model).expect("system should be solved");
        let first_beam = match model.elements()[0] {
            Element2D::Beam(beam) => beam,
            _ => panic!("expected a beam element"),
        };
        let second_beam = match model.elements()[1] {
            Element2D::Beam(beam) => beam,
            _ => panic!("expected a beam element"),
        };
        let load_magnitude = 12.0;
        let young_modulus = 200.0;
        let second_moment_of_area = 2.0;
        let span = 2.0;

        for position in [0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0] {
            let (beam, local_position) =
                if position <= 1.0 { (&first_beam, position) } else { (&second_beam, position - 1.0) };
            let actual = recover_beam_section_response(&model, beam, result.displacements(), local_position)
                .expect("section response should be recovered");
            let distance_from_support = position.min(span - position);
            let expected_moment = load_magnitude * distance_from_support / 2.0;
            let expected_curvature = expected_moment / (young_modulus * second_moment_of_area);

            assert_relative_eq!(actual.bending_moment(), expected_moment, epsilon = 1e-12);
            assert_relative_eq!(actual.curvature(), expected_curvature, epsilon = 1e-12);
        }
    }

    fn right_triangle_model(node_ids: [usize; 3]) -> Model2D {
        let mut model = Model2D::new();

        model.set_material(Material2D::new(200.0, 0.3, 1.0).expect("valid material should be created"));
        model.add_node(Node2D::new(1, 0.0, 0.0).expect("valid node should be created")).expect("node should be added");
        model.add_node(Node2D::new(2, 1.0, 0.0).expect("valid node should be created")).expect("node should be added");
        model.add_node(Node2D::new(3, 0.0, 1.0).expect("valid node should be created")).expect("node should be added");
        add_triangle_element(&mut model, 10, node_ids, 1.0);

        model
    }

    fn rectangular_t3_tension_model() -> Model2D {
        let mut model = Model2D::new();
        let material = Material2D::new(200.0, 0.3, 7800.0).expect("valid material should be created");

        model.set_material(material);

        for (id, x, y) in [(1, 0.0, 0.0), (2, 2.0, 0.0), (3, 2.0, 1.0), (4, 0.0, 1.0)] {
            model.add_node(Node2D::new(id, x, y).expect("valid node should be created")).expect("node should be added");
        }

        add_triangle_element(&mut model, 10, [1, 2, 3], 1.0);
        add_triangle_element(&mut model, 20, [1, 3, 4], 1.0);

        for (node_id, dof) in [(1, Dof2D::Ux), (1, Dof2D::Uy), (4, Dof2D::Ux)] {
            model
                .add_constraint(DisplacementConstraint2D::new(node_id, dof, 0.0).expect("valid constraint"))
                .expect("constraint should be added");
        }

        for node_id in [2, 3] {
            model
                .add_load(NodalLoad2D::new(node_id, Dof2D::Ux, 5.0).expect("valid load"))
                .expect("load should be added");
        }

        model
    }

    fn rectangular_cantilever_t3_model(nx: usize, ny: usize) -> Model2D {
        assert!(nx > 0);
        assert!(ny > 0);

        let length = 10.0;
        let height = 1.0;
        let thickness = 1.0;
        let young_modulus = 1_000.0;
        let poisson_ratio = 0.3;
        let mut model = Model2D::new();

        model.set_material(
            Material2D::new(young_modulus, poisson_ratio, 1.0).expect("valid material should be created"),
        );

        for row in 0..=ny {
            for column in 0..=nx {
                let node_id = row * (nx + 1) + column + 1;
                let x = length * column as f64 / nx as f64;
                let y = height * row as f64 / ny as f64;

                model
                    .add_node(Node2D::new(node_id, x, y).expect("valid node should be created"))
                    .expect("node should be added");
            }
        }

        for row in 0..ny {
            for column in 0..nx {
                let lower_left = row * (nx + 1) + column + 1;
                let lower_right = lower_left + 1;
                let upper_left = lower_left + nx + 1;
                let upper_right = upper_left + 1;
                let element_id = 100 + 2 * (row * nx + column);

                add_triangle_element(&mut model, element_id, [lower_left, lower_right, upper_right], thickness);
                add_triangle_element(&mut model, element_id + 1, [lower_left, upper_right, upper_left], thickness);
            }
        }

        for row in 0..=ny {
            let node_id = row * (nx + 1) + 1;

            for dof in [Dof2D::Ux, Dof2D::Uy] {
                model
                    .add_constraint(DisplacementConstraint2D::new(node_id, dof, 0.0).expect("valid constraint"))
                    .expect("constraint should be added");
            }
        }

        let load_per_right_edge_node = -1.0 / (ny + 1) as f64;

        for row in 0..=ny {
            let node_id = row * (nx + 1) + nx + 1;

            model
                .add_load(NodalLoad2D::new(node_id, Dof2D::Uy, load_per_right_edge_node).expect("valid load"))
                .expect("load should be added");
        }

        model
    }

    #[test]
    fn recovers_uniaxial_plane_stress_for_triangle() {
        let model = right_triangle_model([1, 2, 3]);
        let triangle = match model.elements()[0] {
            Element2D::TriangleT3(triangle) => triangle,
            _ => panic!("expected a triangle element"),
        };
        let poisson_ratio = 0.3;
        let axial_strain = 0.01;
        let displacements = DVector::from_row_slice(&[0.0, 0.0, axial_strain, 0.0, 0.0, -poisson_ratio * axial_strain]);

        let response =
            recover_triangle_response(&model, &triangle, &displacements).expect("response should be recovered");

        assert_relative_eq!(response.strain()[0], axial_strain, epsilon = 1e-12);
        assert_relative_eq!(response.strain()[1], -poisson_ratio * axial_strain, epsilon = 1e-12);
        assert_relative_eq!(response.strain()[2], 0.0, epsilon = 1e-12);
        assert_relative_eq!(response.stress()[0], 2.0, epsilon = 1e-12);
        assert_relative_eq!(response.stress()[1], 0.0, epsilon = 1e-12);
        assert_relative_eq!(response.stress()[2], 0.0, epsilon = 1e-12);
        assert_relative_eq!(response.von_mises_stress(), 2.0, epsilon = 1e-12);
    }

    #[test]
    fn recovers_pure_shear_for_triangle() {
        let model = right_triangle_model([1, 2, 3]);
        let triangle = match model.elements()[0] {
            Element2D::TriangleT3(triangle) => triangle,
            _ => panic!("expected a triangle element"),
        };
        let engineering_shear = 0.02;
        let half_shear = engineering_shear / 2.0;
        let displacements = DVector::from_row_slice(&[0.0, 0.0, 0.0, half_shear, half_shear, 0.0]);

        let response =
            recover_triangle_response(&model, &triangle, &displacements).expect("response should be recovered");
        let expected_shear_stress = 200.0 / (2.0 * (1.0 + 0.3)) * engineering_shear;

        assert_relative_eq!(response.strain()[0], 0.0, epsilon = 1e-12);
        assert_relative_eq!(response.strain()[1], 0.0, epsilon = 1e-12);
        assert_relative_eq!(response.strain()[2], engineering_shear, epsilon = 1e-12);
        assert_relative_eq!(response.stress()[0], 0.0, epsilon = 1e-12);
        assert_relative_eq!(response.stress()[1], 0.0, epsilon = 1e-12);
        assert_relative_eq!(response.stress()[2], expected_shear_stress, epsilon = 1e-12);
        assert_relative_eq!(response.von_mises_stress(), 3.0_f64.sqrt() * expected_shear_stress, epsilon = 1e-12);
    }

    #[test]
    fn interpolates_triangle_displacement_and_passes_affine_patch_test() {
        let model = right_triangle_model([1, 2, 3]);
        let triangle = match model.elements()[0] {
            Element2D::TriangleT3(triangle) => triangle,
            _ => panic!("expected a triangle element"),
        };
        let displacements = DVector::from_row_slice(&[1.0, -1.0, 3.0, 3.0, 4.0, 4.0]);

        let actual = interpolate_triangle_displacement(&model, &triangle, &displacements, 0.2, 0.3)
            .expect("triangle displacement should be interpolated");

        assert_relative_eq!(actual.ux(), 2.3, epsilon = 1e-12);
        assert_relative_eq!(actual.uy(), 1.3, epsilon = 1e-12);

        let response = recover_triangle_response(&model, &triangle, &displacements)
            .expect("triangle response should be recovered");

        assert_relative_eq!(response.strain()[0], 2.0, epsilon = 1e-12);
        assert_relative_eq!(response.strain()[1], 5.0, epsilon = 1e-12);
        assert_relative_eq!(response.strain()[2], 7.0, epsilon = 1e-12);
    }

    #[test]
    fn rejects_triangle_interpolation_outside_natural_domain() {
        let model = right_triangle_model([1, 2, 3]);
        let triangle = match model.elements()[0] {
            Element2D::TriangleT3(triangle) => triangle,
            _ => panic!("expected a triangle element"),
        };
        let displacements = DVector::zeros(6);

        for coordinates in [(-0.1, 0.0), (0.0, -0.1), (0.6, 0.5), (f64::NAN, 0.0), (0.0, f64::INFINITY)] {
            let result =
                interpolate_triangle_displacement(&model, &triangle, &displacements, coordinates.0, coordinates.1);

            assert!(matches!(result, Err(FemError::InvalidTriangleNaturalCoordinates { .. })));
        }
    }

    #[test]
    fn recovers_same_response_for_reversed_triangle_orientation() {
        let model = right_triangle_model([1, 3, 2]);
        let triangle = match model.elements()[0] {
            Element2D::TriangleT3(triangle) => triangle,
            _ => panic!("expected a triangle element"),
        };
        let displacements = DVector::from_row_slice(&[0.0, 0.0, 0.01, 0.0, 0.0, -0.003]);

        let response =
            recover_triangle_response(&model, &triangle, &displacements).expect("response should be recovered");

        assert_relative_eq!(response.strain()[0], 0.01, epsilon = 1e-12);
        assert_relative_eq!(response.strain()[1], -0.003, epsilon = 1e-12);
        assert_relative_eq!(response.strain()[2], 0.0, epsilon = 1e-12);
        assert_relative_eq!(response.stress()[0], 2.0, epsilon = 1e-12);
        assert_relative_eq!(response.stress()[1], 0.0, epsilon = 1e-12);
        assert_relative_eq!(response.stress()[2], 0.0, epsilon = 1e-12);
    }

    #[test]
    fn recovers_triangle_response_from_solver_result() {
        let mut model = right_triangle_model([1, 2, 3]);
        let poisson_ratio = 0.3;
        let axial_strain = 0.01;
        let prescribed_displacements = [
            (1, Dof2D::Ux, 0.0),
            (1, Dof2D::Uy, 0.0),
            (2, Dof2D::Ux, axial_strain),
            (2, Dof2D::Uy, 0.0),
            (3, Dof2D::Ux, 0.0),
            (3, Dof2D::Uy, -poisson_ratio * axial_strain),
        ];

        for (node_id, dof, displacement) in prescribed_displacements {
            model
                .add_constraint(DisplacementConstraint2D::new(node_id, dof, displacement).expect("valid constraint"))
                .expect("constraint should be added");
        }

        let result = solve(&model).expect("system should be solved");
        let triangle = match model.elements()[0] {
            Element2D::TriangleT3(triangle) => triangle,
            _ => panic!("expected a triangle element"),
        };
        let response =
            recover_triangle_response(&model, &triangle, result.displacements()).expect("response should be recovered");

        assert_relative_eq!(response.stress()[0], 2.0, epsilon = 1e-12);
        assert_relative_eq!(response.stress()[1], 0.0, epsilon = 1e-12);
        assert_relative_eq!(response.stress()[2], 0.0, epsilon = 1e-12);
    }

    #[test]
    fn solves_rectangular_t3_membrane_against_analytical_tension_solution() {
        let model = rectangular_t3_tension_model();
        let result = solve(&model).expect("membrane model should be solved");
        let displacements = result.displacements();
        let numbering = crate::model::DofNumbering2D::from_model(&model).expect("DOF numbering should be created");

        let young_modulus = 200.0;
        let poisson_ratio = 0.3;
        let applied_traction = 10.0;
        let expected_epsilon_x = applied_traction / young_modulus;
        let expected_epsilon_y = -poisson_ratio * expected_epsilon_x;
        let expected_nodal_displacements = [
            (1, 0.0, 0.0),
            (2, 2.0 * expected_epsilon_x, 0.0),
            (3, 2.0 * expected_epsilon_x, expected_epsilon_y),
            (4, 0.0, expected_epsilon_y),
        ];

        for (node_id, expected_ux, expected_uy) in expected_nodal_displacements {
            assert_relative_eq!(
                displacements[numbering.index(node_id, Dof2D::Ux).expect("node should have Ux")],
                expected_ux,
                epsilon = 1e-12
            );
            assert_relative_eq!(
                displacements[numbering.index(node_id, Dof2D::Uy).expect("node should have Uy")],
                expected_uy,
                epsilon = 1e-12
            );
        }

        for element in model.elements() {
            let Element2D::TriangleT3(triangle) = element else {
                panic!("expected a triangle element");
            };

            let response = recover_triangle_response(&model, triangle, displacements)
                .expect("triangle response should be recovered");

            assert_relative_eq!(response.strain()[0], expected_epsilon_x, epsilon = 1e-12);
            assert_relative_eq!(response.strain()[1], expected_epsilon_y, epsilon = 1e-12);
            assert_relative_eq!(response.strain()[2], 0.0, epsilon = 1e-12);
            assert_relative_eq!(response.stress()[0], applied_traction, epsilon = 1e-12);
            assert_relative_eq!(response.stress()[1], 0.0, epsilon = 1e-12);
            assert_relative_eq!(response.stress()[2], 0.0, epsilon = 1e-12);
            assert_relative_eq!(response.von_mises_stress(), applied_traction, epsilon = 1e-12);
        }

        assert_relative_eq!(
            result.reactions()[numbering.index(1, Dof2D::Ux).expect("node should have Ux")],
            -5.0,
            epsilon = 1e-12
        );
        assert_relative_eq!(
            result.reactions()[numbering.index(4, Dof2D::Ux).expect("node should have Ux")],
            -5.0,
            epsilon = 1e-12
        );
        assert_relative_eq!(
            result.reactions()[numbering.index(1, Dof2D::Uy).expect("node should have Uy")],
            0.0,
            epsilon = 1e-12
        );
    }

    #[test]
    fn t3_cantilever_deflection_converges_towards_beam_solution() {
        let length: f64 = 10.0;
        let height: f64 = 1.0;
        let thickness: f64 = 1.0;
        let young_modulus: f64 = 1_000.0;
        let applied_load: f64 = -1.0;
        let second_moment_of_area = thickness * height.powi(3) / 12.0;
        let analytical_tip_deflection = applied_load * length.powi(3) / (3.0 * young_modulus * second_moment_of_area);
        let mut relative_errors = Vec::new();

        for (nx, ny) in [(4, 2), (8, 4), (16, 8), (32, 16)] {
            let model = rectangular_cantilever_t3_model(nx, ny);
            let result = solve(&model).expect("cantilever membrane should be solved");
            let numbering = crate::model::DofNumbering2D::from_model(&model).expect("DOF numbering should be created");
            let middle_right_node = (ny / 2) * (nx + 1) + nx + 1;
            let numerical_tip_deflection =
                result.displacements()[numbering.index(middle_right_node, Dof2D::Uy).expect("node should have Uy")];
            let relative_error =
                (numerical_tip_deflection - analytical_tip_deflection).abs() / analytical_tip_deflection.abs();

            assert!(numerical_tip_deflection < 0.0);
            relative_errors.push(relative_error);
        }

        assert!(
            relative_errors[1] < relative_errors[0],
            "mesh refinement should reduce the error: {relative_errors:?}"
        );
        assert!(
            relative_errors[2] < relative_errors[1],
            "mesh refinement should reduce the error: {relative_errors:?}"
        );
        assert!(relative_errors[3] < 0.11, "the dense mesh should bring the error below 11%: {relative_errors:?}");
    }

    #[test]
    #[ignore = "large dense benchmark; run explicitly with cargo test -- --ignored"]
    fn t3_cantilever_64x32_dense_mesh_benchmark() {
        let length: f64 = 10.0;
        let height: f64 = 1.0;
        let thickness: f64 = 1.0;
        let young_modulus: f64 = 1_000.0;
        let applied_load: f64 = -1.0;
        let second_moment_of_area = thickness * height.powi(3) / 12.0;
        let analytical_tip_deflection = applied_load * length.powi(3) / (3.0 * young_modulus * second_moment_of_area);
        let nx = 64;
        let ny = 32;
        let model = rectangular_cantilever_t3_model(nx, ny);
        let result = solve(&model).expect("64x32 cantilever membrane should be solved");
        let numbering = crate::model::DofNumbering2D::from_model(&model).expect("DOF numbering should be created");
        let middle_right_node = (ny / 2) * (nx + 1) + nx + 1;
        let numerical_tip_deflection =
            result.displacements()[numbering.index(middle_right_node, Dof2D::Uy).expect("node should have Uy")];
        let relative_error =
            (numerical_tip_deflection - analytical_tip_deflection).abs() / analytical_tip_deflection.abs();

        eprintln!(
            "64x32 mesh: numerical = {numerical_tip_deflection:.12}, analytical = {analytical_tip_deflection:.12}, relative error = {relative_error:.6}"
        );

        assert!(numerical_tip_deflection < 0.0);
    }

    #[test]
    fn dispatches_truss_response_through_common_interface() {
        let model = horizontal_truss_model();
        let element = &model.elements()[0];
        let displacements = DVector::from_row_slice(&[0.0, 0.0, 0.025, 0.0]);

        let response =
            recover_element_response(&model, element, &displacements).expect("element response should be recovered");

        let ElementResponse2D::Truss(response) = response else {
            panic!("expected a truss response");
        };

        assert_relative_eq!(response.axial_force(), 10.0, epsilon = 1e-12);
    }

    #[test]
    fn dispatches_beam_response_through_common_interface() {
        let model = horizontal_beam_model();
        let element = &model.elements()[0];
        let displacements = DVector::from_row_slice(&[0.0, 0.0, 0.0, 0.01, 0.0, 0.0]);

        let response =
            recover_element_response(&model, element, &displacements).expect("element response should be recovered");

        let ElementResponse2D::Beam(response) = response else {
            panic!("expected a beam response");
        };

        assert_relative_eq!(response.second_axial_force(), 2.0, epsilon = 1e-12);
    }

    #[test]
    fn dispatches_triangle_response_through_common_interface() {
        let model = right_triangle_model([1, 2, 3]);
        let element = &model.elements()[0];
        let displacements = DVector::from_row_slice(&[0.0, 0.0, 0.01, 0.0, 0.0, -0.003]);

        let response =
            recover_element_response(&model, element, &displacements).expect("element response should be recovered");

        let ElementResponse2D::Triangle(response) = response else {
            panic!("expected a triangle response");
        };

        assert_relative_eq!(response.stress()[0], 2.0, epsilon = 1e-12);
    }

    #[test]
    fn recovers_responses_for_all_model_elements_and_preserves_ids() {
        let mut model = right_triangle_model([1, 2, 3]);
        add_beam_element(&mut model, 20, [1, 2], 1.0, 2.0);
        add_truss_element(&mut model, 30, [2, 3], 1.0);

        let responses = recover_model_responses(&model, &DVector::zeros(8)).expect("responses should be recovered");

        assert_eq!(responses.len(), 3);
        assert_eq!(responses[0].0, 10);
        assert_eq!(responses[1].0, 20);
        assert_eq!(responses[2].0, 30);
        assert!(matches!(responses[0].1, ElementResponse2D::Triangle(_)));
        assert!(matches!(responses[1].1, ElementResponse2D::Beam(_)));
        assert!(matches!(responses[2].1, ElementResponse2D::Truss(_)));
    }
}
