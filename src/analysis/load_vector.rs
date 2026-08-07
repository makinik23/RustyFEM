use nalgebra::DVector;

use crate::elements::Element2D;
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
    let Element2D::TriangleT3(_) = element else {
        return Err(FemError::InvalidElementLoadType {
            element_id: load.element_id(),
            load_type: EdgeTraction2D::LOAD_TYPE,
            expected: "triangle_t3",
            actual: element.element_type(),
        });
    };
    let edge_node_ids = load.edge_node_ids();
    let Some(edge_positions) = element_edge_positions(element.node_ids(), edge_node_ids) else {
        return Err(FemError::InvalidElementLoadEdge {
            element_id: load.element_id(),
            load_type: EdgeTraction2D::LOAD_TYPE,
            node_ids: edge_node_ids.to_vec(),
            expected: "one of the triangle's three edges",
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
            element_type: "triangle_t3",
            node_ids: element.node_ids().to_vec(),
            measure_name: "edge length",
            measure: edge_length,
        });
    }

    let (x_component, y_component) = load.global_components(dx / edge_length, dy / edge_length);
    let nodal_factor = section.thickness() * edge_length / 2.0;
    let mut element_load_vector = vec![0.0; element.dof_count()];

    for edge_position in edge_positions {
        element_load_vector[2 * edge_position] += x_component * nodal_factor;
        element_load_vector[2 * edge_position + 1] += y_component * nodal_factor;
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
    assemble_triangle_volume_load(
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

    assemble_triangle_volume_load(
        load_vector,
        numbering,
        model,
        load.element_id(),
        SelfWeight2D::LOAD_TYPE,
        density * load.x_acceleration(),
        density * load.y_acceleration(),
    )
}

fn assemble_triangle_volume_load(
    load_vector: &mut DVector<f64>, numbering: &DofNumbering2D, model: &Model2D, element_id: usize,
    load_type: &'static str, x_force_per_volume: f64, y_force_per_volume: f64,
) -> Result<(), FemError> {
    let element = model
        .elements()
        .iter()
        .find(|element| element.id() == element_id)
        .ok_or(FemError::UnknownId { entity: "element", id: element_id })?;
    let Element2D::TriangleT3(triangle) = element else {
        return Err(FemError::InvalidElementLoadType {
            element_id,
            load_type,
            expected: "triangle_t3",
            actual: element.element_type(),
        });
    };
    let node_ids = element.node_ids();
    let first_node = find_node(model.nodes(), node_ids[0])?;
    let second_node = find_node(model.nodes(), node_ids[1])?;
    let third_node = find_node(model.nodes(), node_ids[2])?;
    let (_, area) = triangle.strain_displacement_matrix(first_node, second_node, third_node)?;
    let section = model.plane_stress_section(element.section_id())?;
    let nodal_factor = area * section.thickness() / 3.0;
    let element_load_vector = [
        x_force_per_volume * nodal_factor,
        y_force_per_volume * nodal_factor,
        x_force_per_volume * nodal_factor,
        y_force_per_volume * nodal_factor,
        x_force_per_volume * nodal_factor,
        y_force_per_volume * nodal_factor,
    ];
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

fn element_edge_positions(element_node_ids: &[usize], edge_node_ids: [usize; 2]) -> Option<[usize; 2]> {
    let first_position = element_node_ids.iter().position(|&node_id| node_id == edge_node_ids[0])?;
    let second_position = element_node_ids.iter().position(|&node_id| node_id == edge_node_ids[1])?;
    let node_count = element_node_ids.len();
    let is_edge =
        (first_position + 1) % node_count == second_position || (second_position + 1) % node_count == first_position;

    is_edge.then_some([first_position, second_position])
}

fn find_node(nodes: &[Node2D], node_id: usize) -> Result<&Node2D, FemError> {
    nodes.iter().find(|node| node.id() == node_id).ok_or(FemError::UnknownId { entity: "node", id: node_id })
}

#[cfg(test)]
mod tests {
    use super::assemble_load_vector;
    use crate::elements::{Beam2D, Element2D, TriangleT3, Truss2D};
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
    fn maps_body_force_to_all_triangle_nodes() {
        let mut model = model_with_triangle();
        let load = ElementLoad2D::BodyForce(BodyForce2D::new(10, 6.0, -12.0).expect("valid load should be created"));

        model.add_element_load(load).expect("body force should be added");

        let actual = assemble_load_vector(&model).expect("load vector should be assembled");
        let expected = DVector::from_row_slice(&[1.0, -2.0, 1.0, -2.0, 1.0, -2.0]);

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
