use nalgebra::DVector;

use crate::error::FemError;
use crate::model::{DofNumbering2D, Model2D};

pub fn assemble_load_vector(model: &Model2D) -> Result<DVector<f64>, FemError> {
    let numbering = DofNumbering2D::from_model(model)?;
    let mut load_vector = DVector::zeros(numbering.count());

    for load in model.loads() {
        let global_index = numbering.index(load.node_id(), load.dof())?;

        load_vector[global_index] += load.value();
    }

    Ok(load_vector)
}

#[cfg(test)]
mod tests {
    use super::assemble_load_vector;
    use crate::elements::{Beam2D, Element2D, Truss2D};
    use crate::error::FemError;
    use crate::model::{
        BeamSection2D, DEFAULT_MATERIAL_ID, Dof2D, Material2D, Model2D, NodalLoad2D, Node2D, Section2D, TrussSection2D,
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
}
