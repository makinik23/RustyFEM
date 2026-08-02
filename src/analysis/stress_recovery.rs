//! Postprocessing of solved FEM results.

use nalgebra::DVector;

use crate::elements::{Beam2D, Element2D, TriangleT3, Truss2D};
use crate::error::FemError;
use crate::model::{DofNumbering2D, Material2D, Model2D, Node2D};

/// Stress and strain results recovered for a 2D truss element.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TrussResponse2D {
    strain: f64,
    stress: f64,
    axial_force: f64,
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
    let material = model.material().ok_or(FemError::MissingMaterial)?;
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

    calculate_truss_response(truss, material, first_node, second_node, element_displacements)
}

fn find_node(nodes: &[Node2D], node_id: usize) -> Result<&Node2D, FemError> {
    nodes.iter().find(|node| node.id() == node_id).ok_or(FemError::UnknownId { entity: "node", id: node_id })
}

fn calculate_truss_response(
    truss: &Truss2D, material: &Material2D, first_node: &Node2D, second_node: &Node2D,
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
    let axial_force = stress * truss.cross_section_area();

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
    let material = model.material().ok_or(FemError::MissingMaterial)?;
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
    let global_element_displacements = [
        global_displacements[indices[0]],
        global_displacements[indices[1]],
        global_displacements[indices[2]],
        global_displacements[indices[3]],
        global_displacements[indices[4]],
        global_displacements[indices[5]],
    ];

    calculate_beam_response(beam, material, first_node, second_node, global_element_displacements)
}

fn calculate_beam_response(
    beam: &Beam2D, material: &Material2D, first_node: &Node2D, second_node: &Node2D,
    [first_x, first_y, first_rotation, second_x, second_y, second_rotation]: [f64; 6],
) -> Result<BeamResponse2D, FemError> {
    let (length, cosine, sine) = beam.geometry(first_node, second_node)?;
    let local_displacements = [
        cosine * first_x + sine * first_y,
        -sine * first_x + cosine * first_y,
        first_rotation,
        cosine * second_x + sine * second_y,
        -sine * second_x + cosine * second_y,
        second_rotation,
    ];
    let local_stiffness_matrix = beam.local_stiffness_matrix(material, length);
    let mut end_forces = [0.0; 6];

    for (row, end_force) in end_forces.iter_mut().enumerate() {
        *end_force = local_stiffness_matrix[row]
            .iter()
            .zip(local_displacements)
            .map(|(stiffness, displacement)| stiffness * displacement)
            .sum();
    }

    Ok(BeamResponse2D { end_forces })
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
}

/// Recovers the constant plane-stress response of a T3 triangle.
pub fn recover_triangle_response(
    model: &Model2D, triangle: &TriangleT3, global_displacements: &DVector<f64>,
) -> Result<TriangleResponse2D, FemError> {
    let material = model.material().ok_or(FemError::MissingMaterial)?;
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
        ElementResponse2D, recover_beam_response, recover_element_response, recover_model_responses,
        recover_triangle_response, recover_truss_response,
    };
    use crate::analysis::solver::solve;
    use crate::elements::{Beam2D, Element2D, TriangleT3, Truss2D};
    use crate::error::FemError;
    use crate::model::{DisplacementConstraint2D, Dof2D, Material2D, Model2D, NodalLoad2D, Node2D};
    use approx::assert_relative_eq;
    use nalgebra::DVector;

    fn horizontal_truss_model() -> Model2D {
        let mut model = Model2D::new();

        model.set_material(Material2D::new(200.0, 0.3, 1.0).expect("valid material should be created"));
        model.add_node(Node2D::new(1, 0.0, 0.0).expect("valid node should be created")).expect("node should be added");
        model.add_node(Node2D::new(2, 1.0, 0.0).expect("valid node should be created")).expect("node should be added");
        model
            .add_element(Element2D::Truss(Truss2D::new(10, [1, 2], 2.0).expect("valid truss should be created")))
            .expect("element should be added");

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
        model
            .add_element(Element2D::Truss(Truss2D::new(10, [1, 2], 1.0).expect("valid truss should be created")))
            .expect("element should be added");

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
        model
            .add_element(Element2D::Beam(Beam2D::new(10, [1, 2], 1.0, 2.0).expect("valid beam should be created")))
            .expect("element should be added");

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
    }

    #[test]
    fn recovers_axial_force_for_diagonal_beam() {
        let mut model = Model2D::new();

        model.set_material(Material2D::new(200.0, 0.3, 1.0).expect("valid material should be created"));
        model.add_node(Node2D::new(1, 0.0, 0.0).expect("valid node should be created")).expect("node should be added");
        model.add_node(Node2D::new(2, 1.0, 1.0).expect("valid node should be created")).expect("node should be added");
        model
            .add_element(Element2D::Beam(Beam2D::new(10, [1, 2], 1.0, 2.0).expect("valid beam should be created")))
            .expect("element should be added");

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

    fn right_triangle_model(node_ids: [usize; 3]) -> Model2D {
        let mut model = Model2D::new();

        model.set_material(Material2D::new(200.0, 0.3, 1.0).expect("valid material should be created"));
        model.add_node(Node2D::new(1, 0.0, 0.0).expect("valid node should be created")).expect("node should be added");
        model.add_node(Node2D::new(2, 1.0, 0.0).expect("valid node should be created")).expect("node should be added");
        model.add_node(Node2D::new(3, 0.0, 1.0).expect("valid node should be created")).expect("node should be added");
        model
            .add_element(Element2D::TriangleT3(
                TriangleT3::new(10, node_ids, 1.0).expect("valid triangle should be created"),
            ))
            .expect("element should be added");

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
        let beam = Beam2D::new(20, [1, 2], 1.0, 2.0).expect("valid beam should be created");
        let truss = Truss2D::new(30, [2, 3], 1.0).expect("valid truss should be created");

        model.add_element(Element2D::Beam(beam)).expect("beam should be added");
        model.add_element(Element2D::Truss(truss)).expect("truss should be added");

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
