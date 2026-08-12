use nalgebra::DVector;

use crate::elements::interpolation::{quad_q4_shape_functions, quad_q8_shape_functions, triangle_t6_shape_functions};
use crate::elements::{Element2D, QuadQ4, QuadQ8};
use crate::error::FemError;
use crate::model::{BeamUniformLineLoad2D, BodyForce2D, DofNumbering2D, EdgeTraction2D, Model2D, Node2D, SelfWeight2D};

pub fn assemble_load_vector(model: &Model2D) -> Result<DVector<f64>, FemError> {
    let numbering = DofNumbering2D::from_model(model)?;
    let mut load_vector = DVector::zeros(numbering.count());

    for load in model.loads() {
        let global_index = numbering.index(load.node_id(), load.dof())?;

        load_vector[global_index] += load.value();
    }

    for load in model.element_loads() {
        match load {
            crate::model::ElementLoad2D::BeamUniformLine(load) => {
                assemble_beam_uniform_line_load(&mut load_vector, &numbering, model, load)?;
            }
            crate::model::ElementLoad2D::EdgeTraction(load) => {
                assemble_edge_traction_load(&mut load_vector, &numbering, model, load)?;
            }
            crate::model::ElementLoad2D::BodyForce(load) => {
                assemble_body_force_load(&mut load_vector, &numbering, model, load)?;
            }
            crate::model::ElementLoad2D::SelfWeight(load) => {
                assemble_self_weight_load(&mut load_vector, &numbering, model, load)?;
            }
        }
    }

    Ok(load_vector)
}

fn assemble_beam_uniform_line_load(
    load_vector: &mut DVector<f64>, numbering: &DofNumbering2D, model: &Model2D, load: &BeamUniformLineLoad2D,
) -> Result<(), FemError> {
    let element = model
        .elements()
        .iter()
        .find(|element| element.id() == load.element_id())
        .ok_or(FemError::UnknownId { entity: "element", id: load.element_id() })?;
    let Element2D::Beam(beam) = element else {
        return Err(FemError::InvalidElementLoadType {
            element_id: load.element_id(),
            load_type: BeamUniformLineLoad2D::LOAD_TYPE,
            expected: "beam",
            actual: element.element_type(),
        });
    };
    let node_ids = element.node_ids();
    let first_node = find_node(model.nodes(), node_ids[0])?;
    let second_node = find_node(model.nodes(), node_ids[1])?;
    let (length, cosine, sine) = beam.geometry(first_node, second_node)?;
    let local_load_vector = load.local_equivalent_nodal_load(length, cosine, sine);
    let element_load_vector = transform_beam_load_to_global(local_load_vector, cosine, sine);
    let indices = numbering.element_dof_indices(element)?;

    for (global_index, value) in indices.iter().zip(element_load_vector) {
        load_vector[*global_index] += value;
    }

    Ok(())
}

fn assemble_edge_traction_load(
    load_vector: &mut DVector<f64>, numbering: &DofNumbering2D, model: &Model2D, load: &EdgeTraction2D,
) -> Result<(), FemError> {
    let element = model
        .elements()
        .iter()
        .find(|element| element.id() == load.element_id())
        .ok_or(FemError::UnknownId { entity: "element", id: load.element_id() })?;
    if !matches!(
        element,
        Element2D::TriangleT3(_) | Element2D::TriangleT6(_) | Element2D::QuadQ4(_) | Element2D::QuadQ8(_)
    ) {
        return Err(FemError::InvalidElementLoadType {
            element_id: load.element_id(),
            load_type: EdgeTraction2D::LOAD_TYPE,
            expected: "plane_stress",
            actual: element.element_type(),
        });
    }
    let edge_node_ids = load.edge_node_ids();
    let Some(edge_positions) = element_edge_positions(element, edge_node_ids) else {
        return Err(FemError::InvalidElementLoadEdge {
            element_id: load.element_id(),
            load_type: EdgeTraction2D::LOAD_TYPE,
            node_ids: edge_node_ids.to_vec(),
            expected: "one of the element's boundary edges",
        });
    };
    let section = model.plane_stress_section(element.section_id())?;
    let first_node = find_node(model.nodes(), edge_node_ids[0])?;
    let second_node = find_node(model.nodes(), edge_node_ids[1])?;
    let dx = second_node.x() - first_node.x();
    let dy = second_node.y() - first_node.y();
    let edge_length = (dx * dx + dy * dy).sqrt();

    if !edge_length.is_finite() || edge_length == 0.0 {
        return Err(FemError::DegenerateElement {
            element_id: load.element_id(),
            element_type: element.element_type(),
            node_ids: element.node_ids().to_vec(),
            measure_name: "edge length",
            measure: edge_length,
        });
    }

    let (x_component, y_component) = load.global_components(dx / edge_length, dy / edge_length);
    let mut element_load_vector = vec![0.0; element.dof_count()];

    match edge_positions {
        EdgeNodePositions::Linear(edge_positions) => {
            let nodal_factor = section.thickness() * edge_length / 2.0;

            for edge_position in edge_positions {
                element_load_vector[2 * edge_position] += x_component * nodal_factor;
                element_load_vector[2 * edge_position + 1] += y_component * nodal_factor;
            }
        }
        EdgeNodePositions::Quadratic(edge_positions) => {
            let element_node_ids = element.node_ids();
            let edge_nodes = [
                find_node(model.nodes(), element_node_ids[edge_positions[0]])?,
                find_node(model.nodes(), element_node_ids[edge_positions[1]])?,
                find_node(model.nodes(), element_node_ids[edge_positions[2]])?,
            ];

            for (coordinate, weight) in quadratic_line_gauss_points() {
                let shape_functions = quadratic_line_shape_functions(coordinate);
                let derivatives = quadratic_line_shape_function_derivatives(coordinate);
                let mut dx_dcoordinate = 0.0;
                let mut dy_dcoordinate = 0.0;

                for node_index in 0..3 {
                    dx_dcoordinate += derivatives[node_index] * edge_nodes[node_index].x();
                    dy_dcoordinate += derivatives[node_index] * edge_nodes[node_index].y();
                }

                let edge_jacobian = (dx_dcoordinate * dx_dcoordinate + dy_dcoordinate * dy_dcoordinate).sqrt();

                if !edge_jacobian.is_finite() || edge_jacobian == 0.0 {
                    return Err(FemError::DegenerateElement {
                        element_id: load.element_id(),
                        element_type: element.element_type(),
                        node_ids: element.node_ids().to_vec(),
                        measure_name: "edge jacobian",
                        measure: edge_jacobian,
                    });
                }

                let scale = section.thickness() * edge_jacobian * weight;

                for (edge_node_index, edge_position) in edge_positions.iter().enumerate() {
                    let nodal_factor = shape_functions[edge_node_index] * scale;

                    element_load_vector[2 * edge_position] += x_component * nodal_factor;
                    element_load_vector[2 * edge_position + 1] += y_component * nodal_factor;
                }
            }
        }
    }

    let indices = numbering.element_dof_indices(element)?;

    for (global_index, value) in indices.iter().zip(element_load_vector) {
        load_vector[*global_index] += value;
    }

    Ok(())
}

fn assemble_body_force_load(
    load_vector: &mut DVector<f64>, numbering: &DofNumbering2D, model: &Model2D, load: &BodyForce2D,
) -> Result<(), FemError> {
    assemble_plane_stress_volume_load(
        load_vector,
        numbering,
        model,
        load.element_id(),
        BodyForce2D::LOAD_TYPE,
        load.x_component(),
        load.y_component(),
    )
}

fn assemble_self_weight_load(
    load_vector: &mut DVector<f64>, numbering: &DofNumbering2D, model: &Model2D, load: &SelfWeight2D,
) -> Result<(), FemError> {
    let element = model
        .elements()
        .iter()
        .find(|element| element.id() == load.element_id())
        .ok_or(FemError::UnknownId { entity: "element", id: load.element_id() })?;
    let material = model.material(element.material_id())?;
    let density = material.density();

    assemble_plane_stress_volume_load(
        load_vector,
        numbering,
        model,
        load.element_id(),
        SelfWeight2D::LOAD_TYPE,
        density * load.x_acceleration(),
        density * load.y_acceleration(),
    )
}

fn assemble_plane_stress_volume_load(
    load_vector: &mut DVector<f64>, numbering: &DofNumbering2D, model: &Model2D, element_id: usize,
    load_type: &'static str, x_force_per_volume: f64, y_force_per_volume: f64,
) -> Result<(), FemError> {
    let element = model
        .elements()
        .iter()
        .find(|element| element.id() == element_id)
        .ok_or(FemError::UnknownId { entity: "element", id: element_id })?;
    let element_load_vector = match element {
        Element2D::TriangleT3(triangle) => {
            let node_ids = element.node_ids();
            let first_node = find_node(model.nodes(), node_ids[0])?;
            let second_node = find_node(model.nodes(), node_ids[1])?;
            let third_node = find_node(model.nodes(), node_ids[2])?;
            let (_, area) = triangle.strain_displacement_matrix(first_node, second_node, third_node)?;
            let section = model.plane_stress_section(element.section_id())?;
            let nodal_factor = area * section.thickness() / 3.0;

            vec![
                x_force_per_volume * nodal_factor,
                y_force_per_volume * nodal_factor,
                x_force_per_volume * nodal_factor,
                y_force_per_volume * nodal_factor,
                x_force_per_volume * nodal_factor,
                y_force_per_volume * nodal_factor,
            ]
        }
        Element2D::TriangleT6(triangle) => {
            let node_ids = element.node_ids();
            let nodes = [
                find_node(model.nodes(), node_ids[0])?,
                find_node(model.nodes(), node_ids[1])?,
                find_node(model.nodes(), node_ids[2])?,
                find_node(model.nodes(), node_ids[3])?,
                find_node(model.nodes(), node_ids[4])?,
                find_node(model.nodes(), node_ids[5])?,
            ];
            let section = model.plane_stress_section(element.section_id())?;
            let mut element_load_vector = vec![0.0; element.dof_count()];

            for (xi, eta, weight) in crate::elements::TriangleT6::gauss_points() {
                let (_, jacobian_determinant) = triangle.strain_displacement_matrix(nodes, xi, eta)?;
                let shape_functions = triangle_t6_shape_functions(xi, eta);
                let scale = section.thickness() * jacobian_determinant * weight;

                for (node_index, shape_function) in shape_functions.iter().enumerate() {
                    let nodal_factor = shape_function * scale;

                    element_load_vector[2 * node_index] += x_force_per_volume * nodal_factor;
                    element_load_vector[2 * node_index + 1] += y_force_per_volume * nodal_factor;
                }
            }

            element_load_vector
        }
        Element2D::QuadQ4(quad) => {
            let node_ids = element.node_ids();
            let nodes = [
                find_node(model.nodes(), node_ids[0])?,
                find_node(model.nodes(), node_ids[1])?,
                find_node(model.nodes(), node_ids[2])?,
                find_node(model.nodes(), node_ids[3])?,
            ];
            let section = model.plane_stress_section(element.section_id())?;
            let mut element_load_vector = vec![0.0; element.dof_count()];

            for (xi, eta) in QuadQ4::gauss_points() {
                let (_, jacobian_determinant) = quad.strain_displacement_matrix(nodes, xi, eta)?;
                let shape_functions = quad_q4_shape_functions(xi, eta);
                let scale = section.thickness() * jacobian_determinant;

                for (node_index, shape_function) in shape_functions.iter().enumerate() {
                    let nodal_factor = shape_function * scale;

                    element_load_vector[2 * node_index] += x_force_per_volume * nodal_factor;
                    element_load_vector[2 * node_index + 1] += y_force_per_volume * nodal_factor;
                }
            }

            element_load_vector
        }
        Element2D::QuadQ8(quad) => {
            let node_ids = element.node_ids();
            let nodes = [
                find_node(model.nodes(), node_ids[0])?,
                find_node(model.nodes(), node_ids[1])?,
                find_node(model.nodes(), node_ids[2])?,
                find_node(model.nodes(), node_ids[3])?,
                find_node(model.nodes(), node_ids[4])?,
                find_node(model.nodes(), node_ids[5])?,
                find_node(model.nodes(), node_ids[6])?,
                find_node(model.nodes(), node_ids[7])?,
            ];
            let section = model.plane_stress_section(element.section_id())?;
            let mut element_load_vector = vec![0.0; element.dof_count()];

            for (xi, eta, weight) in QuadQ8::gauss_points() {
                let (_, jacobian_determinant) = quad.strain_displacement_matrix(nodes, xi, eta)?;
                let shape_functions = quad_q8_shape_functions(xi, eta);
                let scale = section.thickness() * jacobian_determinant * weight;

                for (node_index, shape_function) in shape_functions.iter().enumerate() {
                    let nodal_factor = shape_function * scale;

                    element_load_vector[2 * node_index] += x_force_per_volume * nodal_factor;
                    element_load_vector[2 * node_index + 1] += y_force_per_volume * nodal_factor;
                }
            }

            element_load_vector
        }
        _ => {
            return Err(FemError::InvalidElementLoadType {
                element_id,
                load_type,
                expected: "plane_stress",
                actual: element.element_type(),
            });
        }
    };
    let indices = numbering.element_dof_indices(element)?;

    for (global_index, value) in indices.iter().zip(element_load_vector) {
        load_vector[*global_index] += value;
    }

    Ok(())
}

fn transform_beam_load_to_global(local_load_vector: [f64; 6], cosine: f64, sine: f64) -> [f64; 6] {
    [
        cosine * local_load_vector[0] - sine * local_load_vector[1],
        sine * local_load_vector[0] + cosine * local_load_vector[1],
        local_load_vector[2],
        cosine * local_load_vector[3] - sine * local_load_vector[4],
        sine * local_load_vector[3] + cosine * local_load_vector[4],
        local_load_vector[5],
    ]
}

enum EdgeNodePositions {
    Linear([usize; 2]),
    Quadratic([usize; 3]),
}

fn element_edge_positions(element: &Element2D, edge_node_ids: [usize; 2]) -> Option<EdgeNodePositions> {
    match element {
        Element2D::TriangleT6(_) => {
            triangle_t6_edge_positions(element.node_ids(), edge_node_ids).map(EdgeNodePositions::Quadratic)
        }
        Element2D::QuadQ8(_) => {
            quad_q8_edge_positions(element.node_ids(), edge_node_ids).map(EdgeNodePositions::Quadratic)
        }
        _ => linear_edge_positions(element.node_ids(), edge_node_ids).map(EdgeNodePositions::Linear),
    }
}

fn linear_edge_positions(element_node_ids: &[usize], edge_node_ids: [usize; 2]) -> Option<[usize; 2]> {
    let first_position = element_node_ids.iter().position(|&node_id| node_id == edge_node_ids[0])?;
    let second_position = element_node_ids.iter().position(|&node_id| node_id == edge_node_ids[1])?;
    let node_count = element_node_ids.len();
    let is_edge =
        (first_position + 1) % node_count == second_position || (second_position + 1) % node_count == first_position;

    is_edge.then_some([first_position, second_position])
}

fn triangle_t6_edge_positions(element_node_ids: &[usize], edge_node_ids: [usize; 2]) -> Option<[usize; 3]> {
    let edges = [([0, 1], 3), ([1, 2], 4), ([2, 0], 5)];

    for ([first_corner, second_corner], midside) in edges {
        let first_node_id = element_node_ids[first_corner];
        let second_node_id = element_node_ids[second_corner];

        if edge_node_ids == [first_node_id, second_node_id] {
            return Some([first_corner, midside, second_corner]);
        }

        if edge_node_ids == [second_node_id, first_node_id] {
            return Some([second_corner, midside, first_corner]);
        }
    }

    None
}

fn quad_q8_edge_positions(element_node_ids: &[usize], edge_node_ids: [usize; 2]) -> Option<[usize; 3]> {
    let edges = [([0, 1], 4), ([1, 2], 5), ([2, 3], 6), ([3, 0], 7)];

    for ([first_corner, second_corner], midside) in edges {
        let first_node_id = element_node_ids[first_corner];
        let second_node_id = element_node_ids[second_corner];

        if edge_node_ids == [first_node_id, second_node_id] {
            return Some([first_corner, midside, second_corner]);
        }

        if edge_node_ids == [second_node_id, first_node_id] {
            return Some([second_corner, midside, first_corner]);
        }
    }

    None
}

fn quadratic_line_gauss_points() -> [(f64, f64); 3] {
    let point = (3.0_f64 / 5.0).sqrt();

    [(-point, 5.0 / 9.0), (0.0, 8.0 / 9.0), (point, 5.0 / 9.0)]
}

fn quadratic_line_shape_functions(coordinate: f64) -> [f64; 3] {
    [0.5 * coordinate * (coordinate - 1.0), 1.0 - coordinate.powi(2), 0.5 * coordinate * (coordinate + 1.0)]
}

fn quadratic_line_shape_function_derivatives(coordinate: f64) -> [f64; 3] {
    [coordinate - 0.5, -2.0 * coordinate, coordinate + 0.5]
}

fn find_node(nodes: &[Node2D], node_id: usize) -> Result<&Node2D, FemError> {
    nodes.iter().find(|node| node.id() == node_id).ok_or(FemError::UnknownId { entity: "node", id: node_id })
}

#[cfg(test)]
mod tests {
    use super::assemble_load_vector;
    use crate::elements::{Beam2D, Element2D, QuadQ4, QuadQ8, TriangleT3, TriangleT6, Truss2D};
    use crate::error::FemError;
    use crate::model::{
        BeamSection2D, BeamUniformLineLoad2D, BodyForce2D, DEFAULT_MATERIAL_ID, Dof2D, EdgeTraction2D, ElementLoad2D,
        LoadCoordinateSystem2D, Material2D, Model2D, NodalLoad2D, Node2D, PlaneStressSection2D, Section2D,
        SelfWeight2D, TrussSection2D,
    };
    use nalgebra::DVector;

    #[test]
    fn creates_zero_load_vector_when_model_has_no_loads() {
        let model = model_with_beam();

        let actual = assemble_load_vector(&model).expect("load vector should be assembled");
        let expected = DVector::from_row_slice(&[0.0, 0.0, 0.0, 0.0, 0.0, 0.0]);

        assert_eq!(actual, expected);
    }

    #[test]
    fn maps_nodal_loads_to_global_dof_indices() {
        let mut model = model_with_beam();

        model
            .add_load(NodalLoad2D::new(1, Dof2D::Ux, 5.0).expect("valid load should be created"))
            .expect("load should be added");
        model
            .add_load(NodalLoad2D::new(2, Dof2D::Rz, -2.0).expect("valid load should be created"))
            .expect("load should be added");

        let actual = assemble_load_vector(&model).expect("load vector should be assembled");
        let expected = DVector::from_row_slice(&[5.0, 0.0, 0.0, 0.0, 0.0, -2.0]);

        assert_eq!(actual, expected);
    }

    #[test]
    fn sums_loads_applied_to_the_same_dof() {
        let mut model = model_with_truss();

        model
            .add_load(NodalLoad2D::new(2, Dof2D::Ux, 10.0).expect("valid load should be created"))
            .expect("first load should be added");
        model
            .add_load(NodalLoad2D::new(2, Dof2D::Ux, -3.0).expect("valid load should be created"))
            .expect("second load should be added");

        let actual = assemble_load_vector(&model).expect("load vector should be assembled");
        let expected = DVector::from_row_slice(&[0.0, 0.0, 7.0, 0.0]);

        assert_eq!(actual, expected);
    }

    #[test]
    fn maps_uniform_beam_line_load_to_equivalent_nodal_loads() {
        let mut model = model_with_beam();
        let load = ElementLoad2D::BeamUniformLine(
            BeamUniformLineLoad2D::new(10, LoadCoordinateSystem2D::Local, 2.0, -6.0)
                .expect("valid load should be created"),
        );

        model.add_element_load(load).expect("element load should be added");

        let actual = assemble_load_vector(&model).expect("load vector should be assembled");
        let expected = DVector::from_row_slice(&[1.0, -3.0, -0.5, 1.0, -3.0, 0.5]);

        assert_eq!(actual, expected);
    }

    #[test]
    fn transforms_global_uniform_beam_line_load_to_global_dofs() {
        let mut model = model_with_diagonal_beam();
        let load = ElementLoad2D::BeamUniformLine(
            BeamUniformLineLoad2D::new(10, LoadCoordinateSystem2D::Global, 0.0, -10.0)
                .expect("valid load should be created"),
        );

        model.add_element_load(load).expect("element load should be added");

        let actual = assemble_load_vector(&model).expect("load vector should be assembled");
        let expected = DVector::from_row_slice(&[0.0, -25.0, -12.5, 0.0, -25.0, 12.5]);

        assert_vector_approximately_equal(&actual, &expected);
    }

    #[test]
    fn maps_global_edge_traction_to_triangle_edge_nodes() {
        let mut model = model_with_triangle();
        let load = ElementLoad2D::EdgeTraction(
            EdgeTraction2D::new(10, [2, 3], LoadCoordinateSystem2D::Global, 0.0, -10.0)
                .expect("valid load should be created"),
        );

        model.add_element_load(load).expect("edge traction should be added");

        let actual = assemble_load_vector(&model).expect("load vector should be assembled");
        let nodal_force = -10.0 * 0.5 * 5.0_f64.sqrt() / 2.0;
        let expected = DVector::from_row_slice(&[0.0, 0.0, 0.0, nodal_force, 0.0, nodal_force]);

        assert_vector_approximately_equal(&actual, &expected);
    }

    #[test]
    fn maps_local_edge_traction_to_triangle_edge_nodes() {
        let mut model = model_with_triangle();
        let load = ElementLoad2D::EdgeTraction(
            EdgeTraction2D::new(10, [1, 2], LoadCoordinateSystem2D::Local, 4.0, -8.0)
                .expect("valid load should be created"),
        );

        model.add_element_load(load).expect("edge traction should be added");

        let actual = assemble_load_vector(&model).expect("load vector should be assembled");
        let expected = DVector::from_row_slice(&[2.0, -4.0, 2.0, -4.0, 0.0, 0.0]);

        assert_vector_approximately_equal(&actual, &expected);
    }

    #[test]
    fn maps_edge_traction_to_quad_edge_nodes() {
        let mut model = model_with_quad();
        let load = ElementLoad2D::EdgeTraction(
            EdgeTraction2D::new(10, [2, 3], LoadCoordinateSystem2D::Global, 8.0, -4.0)
                .expect("valid load should be created"),
        );

        model.add_element_load(load).expect("edge traction should be added");

        let actual = assemble_load_vector(&model).expect("load vector should be assembled");
        let expected = DVector::from_row_slice(&[0.0, 0.0, 2.0, -1.0, 2.0, -1.0, 0.0, 0.0]);

        assert_vector_approximately_equal(&actual, &expected);
    }

    #[test]
    fn maps_edge_traction_to_quad_q8_edge_nodes() {
        let mut model = model_with_quad_q8();
        let load = ElementLoad2D::EdgeTraction(
            EdgeTraction2D::new(10, [2, 3], LoadCoordinateSystem2D::Global, 8.0, -4.0)
                .expect("valid load should be created"),
        );

        model.add_element_load(load).expect("edge traction should be added");

        let actual = assemble_load_vector(&model).expect("load vector should be assembled");
        let expected = DVector::from_row_slice(&[
            0.0,
            0.0,
            2.0 / 3.0,
            -1.0 / 3.0,
            2.0 / 3.0,
            -1.0 / 3.0,
            0.0,
            0.0,
            0.0,
            0.0,
            8.0 / 3.0,
            -4.0 / 3.0,
            0.0,
            0.0,
            0.0,
            0.0,
        ]);

        assert_vector_approximately_equal(&actual, &expected);
    }

    #[test]
    fn maps_edge_traction_to_triangle_t6_edge_nodes() {
        let mut model = model_with_triangle_t6();
        let load = ElementLoad2D::EdgeTraction(
            EdgeTraction2D::new(10, [2, 3], LoadCoordinateSystem2D::Global, 0.0, -10.0)
                .expect("valid load should be created"),
        );

        model.add_element_load(load).expect("edge traction should be added");

        let actual = assemble_load_vector(&model).expect("load vector should be assembled");
        let corner_force = -10.0 * 0.5 * 5.0_f64.sqrt() / 6.0;
        let midside_force = -10.0 * 0.5 * 2.0 * 5.0_f64.sqrt() / 3.0;
        let expected = DVector::from_row_slice(&[
            0.0,
            0.0,
            0.0,
            corner_force,
            0.0,
            corner_force,
            0.0,
            0.0,
            0.0,
            midside_force,
            0.0,
            0.0,
        ]);

        assert_vector_approximately_equal(&actual, &expected);
    }

    #[test]
    fn maps_body_force_to_all_triangle_nodes() {
        let mut model = model_with_triangle();
        let load = ElementLoad2D::BodyForce(BodyForce2D::new(10, 6.0, -12.0).expect("valid load should be created"));

        model.add_element_load(load).expect("body force should be added");

        let actual = assemble_load_vector(&model).expect("load vector should be assembled");
        let expected = DVector::from_row_slice(&[1.0, -2.0, 1.0, -2.0, 1.0, -2.0]);

        assert_vector_approximately_equal(&actual, &expected);
    }

    #[test]
    fn maps_body_force_to_all_quad_nodes_with_gauss_integration() {
        let mut model = model_with_quad();
        let load = ElementLoad2D::BodyForce(BodyForce2D::new(10, 8.0, -4.0).expect("valid load should be created"));

        model.add_element_load(load).expect("body force should be added");

        let actual = assemble_load_vector(&model).expect("load vector should be assembled");
        let expected = DVector::from_row_slice(&[2.0, -1.0, 2.0, -1.0, 2.0, -1.0, 2.0, -1.0]);

        assert_vector_approximately_equal(&actual, &expected);
    }

    #[test]
    fn maps_body_force_to_quad_q8_nodes_with_gauss_integration() {
        let mut model = model_with_quad_q8();
        let load = ElementLoad2D::BodyForce(BodyForce2D::new(10, 8.0, -4.0).expect("valid load should be created"));

        model.add_element_load(load).expect("body force should be added");

        let actual = assemble_load_vector(&model).expect("load vector should be assembled");
        let expected = DVector::from_row_slice(&[
            -2.0 / 3.0,
            1.0 / 3.0,
            -2.0 / 3.0,
            1.0 / 3.0,
            -2.0 / 3.0,
            1.0 / 3.0,
            -2.0 / 3.0,
            1.0 / 3.0,
            8.0 / 3.0,
            -4.0 / 3.0,
            8.0 / 3.0,
            -4.0 / 3.0,
            8.0 / 3.0,
            -4.0 / 3.0,
            8.0 / 3.0,
            -4.0 / 3.0,
        ]);

        assert_vector_approximately_equal(&actual, &expected);
    }

    #[test]
    fn maps_body_force_to_triangle_t6_nodes_with_gauss_integration() {
        let mut model = model_with_triangle_t6();
        let load = ElementLoad2D::BodyForce(BodyForce2D::new(10, 6.0, -12.0).expect("valid load should be created"));

        model.add_element_load(load).expect("body force should be added");

        let actual = assemble_load_vector(&model).expect("load vector should be assembled");
        let expected = DVector::from_row_slice(&[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -2.0, 1.0, -2.0, 1.0, -2.0]);

        assert_vector_approximately_equal(&actual, &expected);
    }

    #[test]
    fn maps_self_weight_to_all_triangle_nodes_using_material_density() {
        let mut model = model_with_triangle_density(2.0);
        let load = ElementLoad2D::SelfWeight(SelfWeight2D::new(10, 6.0, -12.0).expect("valid load should be created"));

        model.add_element_load(load).expect("self-weight load should be added");

        let actual = assemble_load_vector(&model).expect("load vector should be assembled");
        let expected = DVector::from_row_slice(&[2.0, -4.0, 2.0, -4.0, 2.0, -4.0]);

        assert_vector_approximately_equal(&actual, &expected);
    }

    #[test]
    fn maps_self_weight_to_all_quad_nodes_using_material_density() {
        let mut model = model_with_quad_density(2.0);
        let load = ElementLoad2D::SelfWeight(SelfWeight2D::new(10, 5.0, -10.0).expect("valid load should be created"));

        model.add_element_load(load).expect("self-weight load should be added");

        let actual = assemble_load_vector(&model).expect("load vector should be assembled");
        let expected = DVector::from_row_slice(&[2.5, -5.0, 2.5, -5.0, 2.5, -5.0, 2.5, -5.0]);

        assert_vector_approximately_equal(&actual, &expected);
    }

    #[test]
    fn maps_self_weight_to_quad_q8_nodes_using_material_density() {
        let mut model = model_with_quad_q8_density(2.0);
        let load = ElementLoad2D::SelfWeight(SelfWeight2D::new(10, 5.0, -10.0).expect("valid load should be created"));

        model.add_element_load(load).expect("self-weight load should be added");

        let actual = assemble_load_vector(&model).expect("load vector should be assembled");
        let expected = DVector::from_row_slice(&[
            -5.0 / 6.0,
            5.0 / 3.0,
            -5.0 / 6.0,
            5.0 / 3.0,
            -5.0 / 6.0,
            5.0 / 3.0,
            -5.0 / 6.0,
            5.0 / 3.0,
            10.0 / 3.0,
            -20.0 / 3.0,
            10.0 / 3.0,
            -20.0 / 3.0,
            10.0 / 3.0,
            -20.0 / 3.0,
            10.0 / 3.0,
            -20.0 / 3.0,
        ]);

        assert_vector_approximately_equal(&actual, &expected);
    }

    #[test]
    fn rejects_load_on_dof_not_used_by_element() {
        let mut model = model_with_truss();

        model
            .add_load(NodalLoad2D::new(1, Dof2D::Rz, 1.0).expect("valid load should be created"))
            .expect("load should be added");

        let result = assemble_load_vector(&model);

        assert!(matches!(result, Err(FemError::UnknownDof { node_id: 1, dof: "Rz" })));
    }

    fn model_with_beam() -> Model2D {
        let mut model = Model2D::new();
        model.set_material(Material2D::new(200.0, 0.3, 1.0).expect("valid material should be created"));

        model.add_node(Node2D::new(1, 0.0, 0.0).expect("valid node should be created")).expect("node should be added");
        model.add_node(Node2D::new(2, 1.0, 0.0).expect("valid node should be created")).expect("node should be added");

        let beam = Beam2D::new(10, [1, 2], DEFAULT_MATERIAL_ID, 100).expect("valid beam should be created");
        let section = Section2D::Beam(BeamSection2D::new(1.0, 1.0).expect("valid section should be created"));
        model.add_element_with_section(Element2D::Beam(beam), section).expect("beam should be added");

        model
    }

    fn model_with_diagonal_beam() -> Model2D {
        let mut model = Model2D::new();
        model.set_material(Material2D::new(200.0, 0.3, 1.0).expect("valid material should be created"));

        model.add_node(Node2D::new(1, 0.0, 0.0).expect("valid node should be created")).expect("node should be added");
        model.add_node(Node2D::new(2, 3.0, 4.0).expect("valid node should be created")).expect("node should be added");

        let beam = Beam2D::new(10, [1, 2], DEFAULT_MATERIAL_ID, 100).expect("valid beam should be created");
        let section = Section2D::Beam(BeamSection2D::new(1.0, 1.0).expect("valid section should be created"));
        model.add_element_with_section(Element2D::Beam(beam), section).expect("beam should be added");

        model
    }

    fn model_with_triangle() -> Model2D {
        model_with_triangle_density(1.0)
    }

    fn model_with_triangle_density(density: f64) -> Model2D {
        let mut model = Model2D::new();
        model.set_material(Material2D::new(200.0, 0.3, density).expect("valid material should be created"));

        model.add_node(Node2D::new(1, 0.0, 0.0).expect("valid node should be created")).expect("node should be added");
        model.add_node(Node2D::new(2, 2.0, 0.0).expect("valid node should be created")).expect("node should be added");
        model.add_node(Node2D::new(3, 0.0, 1.0).expect("valid node should be created")).expect("node should be added");

        let triangle =
            TriangleT3::new(10, [1, 2, 3], DEFAULT_MATERIAL_ID, 100).expect("valid triangle should be created");
        let section = Section2D::PlaneStress(PlaneStressSection2D::new(0.5).expect("valid section should be created"));
        model.add_element_with_section(Element2D::TriangleT3(triangle), section).expect("triangle should be added");

        model
    }

    fn model_with_triangle_t6() -> Model2D {
        let mut model = Model2D::new();
        model.set_material(Material2D::new(200.0, 0.3, 1.0).expect("valid material should be created"));

        for (id, x, y) in [(1, 0.0, 0.0), (2, 2.0, 0.0), (3, 0.0, 1.0), (4, 1.0, 0.0), (5, 1.0, 0.5), (6, 0.0, 0.5)] {
            model.add_node(Node2D::new(id, x, y).expect("valid node should be created")).expect("node should be added");
        }

        let triangle = TriangleT6::new(10, [1, 2, 3, 4, 5, 6], DEFAULT_MATERIAL_ID, 100)
            .expect("valid triangle should be created");
        let section = Section2D::PlaneStress(PlaneStressSection2D::new(0.5).expect("valid section should be created"));
        model.add_element_with_section(Element2D::TriangleT6(triangle), section).expect("triangle should be added");

        model
    }

    fn model_with_quad() -> Model2D {
        model_with_quad_density(1.0)
    }

    fn model_with_quad_density(density: f64) -> Model2D {
        let mut model = Model2D::new();
        model.set_material(Material2D::new(200.0, 0.3, density).expect("valid material should be created"));

        for (id, x, y) in [(1, 0.0, 0.0), (2, 2.0, 0.0), (3, 2.0, 1.0), (4, 0.0, 1.0)] {
            model.add_node(Node2D::new(id, x, y).expect("valid node should be created")).expect("node should be added");
        }

        let quad = QuadQ4::new(10, [1, 2, 3, 4], DEFAULT_MATERIAL_ID, 100).expect("valid quad should be created");
        let section = Section2D::PlaneStress(PlaneStressSection2D::new(0.5).expect("valid section should be created"));
        model.add_element_with_section(Element2D::QuadQ4(quad), section).expect("quad should be added");

        model
    }

    fn model_with_quad_q8() -> Model2D {
        model_with_quad_q8_density(1.0)
    }

    fn model_with_quad_q8_density(density: f64) -> Model2D {
        let mut model = Model2D::new();
        model.set_material(Material2D::new(200.0, 0.3, density).expect("valid material should be created"));

        for (id, x, y) in [
            (1, 0.0, 0.0),
            (2, 2.0, 0.0),
            (3, 2.0, 1.0),
            (4, 0.0, 1.0),
            (5, 1.0, 0.0),
            (6, 2.0, 0.5),
            (7, 1.0, 1.0),
            (8, 0.0, 0.5),
        ] {
            model.add_node(Node2D::new(id, x, y).expect("valid node should be created")).expect("node should be added");
        }

        let quad =
            QuadQ8::new(10, [1, 2, 3, 4, 5, 6, 7, 8], DEFAULT_MATERIAL_ID, 100).expect("valid quad should be created");
        let section = Section2D::PlaneStress(PlaneStressSection2D::new(0.5).expect("valid section should be created"));
        model.add_element_with_section(Element2D::QuadQ8(quad), section).expect("quad should be added");

        model
    }

    fn model_with_truss() -> Model2D {
        let mut model = Model2D::new();
        model.set_material(Material2D::new(200.0, 0.3, 1.0).expect("valid material should be created"));

        model.add_node(Node2D::new(1, 0.0, 0.0).expect("valid node should be created")).expect("node should be added");
        model.add_node(Node2D::new(2, 1.0, 0.0).expect("valid node should be created")).expect("node should be added");

        let truss = Truss2D::new(10, [1, 2], DEFAULT_MATERIAL_ID, 100).expect("valid truss should be created");
        let section = Section2D::Truss(TrussSection2D::new(1.0).expect("valid section should be created"));
        model.add_element_with_section(Element2D::Truss(truss), section).expect("truss should be added");

        model
    }

    fn assert_vector_approximately_equal(actual: &DVector<f64>, expected: &DVector<f64>) {
        assert_eq!(actual.len(), expected.len());

        for (index, (actual, expected)) in actual.iter().zip(expected.iter()).enumerate() {
            assert!(
                (actual - expected).abs() < 1e-12,
                "different vector entry at index {index}: actual = {actual}, expected = {expected}"
            );
        }
    }
}
