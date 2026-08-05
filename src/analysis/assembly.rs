//! Contains functions for assembling the global stiffness matrix for a 2D finite element model.

use nalgebra::DMatrix;

use crate::analysis::sparse::{CooMatrix, CsrMatrix};
use crate::error::FemError;
use crate::model::{DofNumbering2D, Model2D};

/// Assembles the global stiffness matrix for a given 2D finite element model.
pub fn assemble_stiffness_matrix(model: &Model2D) -> Result<DMatrix<f64>, FemError> {
    let numbering = DofNumbering2D::from_model(model)?;

    let size = numbering.count();
    let mut global_matrix = DMatrix::<f64>::zeros(size, size);

    for element in model.elements() {
        let material = model.material(element.material_id())?;
        let section = model.section(element.section_id())?;
        let element_matrix = element.stiffness_matrix(material, section, model.nodes())?;

        let indices = numbering.element_dof_indices(element)?;

        for (local_row, &global_row) in indices.iter().enumerate() {
            for (local_column, &global_column) in indices.iter().enumerate() {
                global_matrix[(global_row, global_column)] += element_matrix[(local_row, local_column)];
            }
        }
    }

    Ok(global_matrix)
}

/// Assembles the global stiffness matrix in sparse CSR format.
pub fn assemble_sparse_stiffness_matrix(model: &Model2D) -> Result<CsrMatrix, FemError> {
    let numbering = DofNumbering2D::from_model(model)?;
    let size = numbering.count();

    // COO is convenient during assembly because it accepts duplicate entries.
    let mut global_matrix = CooMatrix::new(size, size);

    for element in model.elements() {
        let material = model.material(element.material_id())?;
        let section = model.section(element.section_id())?;
        let element_matrix = element.stiffness_matrix(material, section, model.nodes())?;
        let indices = numbering.element_dof_indices(element)?;

        for (local_row, &global_row) in indices.iter().enumerate() {
            for (local_column, &global_column) in indices.iter().enumerate() {
                let value = element_matrix[(local_row, local_column)];

                global_matrix.push(global_row, global_column, value)?;
            }
        }
    }

    // This sums duplicate entries and creates the efficient CSR structure.
    Ok(global_matrix.into_csr())
}

#[cfg(test)]
mod tests {
    use super::{assemble_sparse_stiffness_matrix, assemble_stiffness_matrix};
    use crate::elements::{Beam2D, Element2D, TriangleT3, Truss2D};
    use crate::model::{
        BeamSection2D, DEFAULT_MATERIAL_ID, Material2D, Model2D, Node2D, PlaneStressSection2D, Section2D,
        TrussSection2D,
    };
    use nalgebra::{DMatrix, DVector};

    #[test]
    fn assembles_empty_model_without_material() {
        let model = Model2D::new();

        let matrix = assemble_stiffness_matrix(&model).expect("empty matrix should be assembled");

        assert_eq!(matrix.shape(), (0, 0));
    }

    #[test]
    fn assembles_global_matrix_for_one_horizontal_truss() {
        let mut model = Model2D::new();
        let material = Material2D::new(2.0, 0.3, 1.0).expect("valid material should be created");

        model.set_material(material);
        model.add_node(Node2D::new(1, 0.0, 0.0).expect("valid node should be created")).expect("node should be added");
        model.add_node(Node2D::new(2, 1.0, 0.0).expect("valid node should be created")).expect("node should be added");

        let truss = Truss2D::new(10, [1, 2], DEFAULT_MATERIAL_ID, 100).expect("valid truss should be created");
        let section = Section2D::Truss(TrussSection2D::new(3.0).expect("valid section should be created"));
        model.add_element_with_section(Element2D::Truss(truss), section).expect("element should be added");

        let actual = assemble_stiffness_matrix(&model).expect("global matrix should be assembled");
        let expected = DMatrix::from_row_slice(
            4,
            4,
            &[6.0, 0.0, -6.0, 0.0, 0.0, 0.0, 0.0, 0.0, -6.0, 0.0, 6.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        );

        assert_matrix_approximately_equal(&actual, &expected);
    }

    #[test]
    fn sums_contributions_from_two_trusses_sharing_a_node() {
        let mut model = Model2D::new();
        let material = Material2D::new(2.0, 0.3, 1.0).expect("valid material should be created");

        model.set_material(material);

        for (id, x) in [(1, 0.0), (2, 1.0), (3, 2.0)] {
            model
                .add_node(Node2D::new(id, x, 0.0).expect("valid node should be created"))
                .expect("node should be added");
        }

        let section = Section2D::Truss(TrussSection2D::new(3.0).expect("valid section should be created"));
        model.add_section(100, section).expect("section should be added");

        let first_truss = Truss2D::new(10, [1, 2], DEFAULT_MATERIAL_ID, 100).expect("valid truss should be created");
        let second_truss = Truss2D::new(20, [2, 3], DEFAULT_MATERIAL_ID, 100).expect("valid truss should be created");

        model.add_element(Element2D::Truss(first_truss)).expect("first element should be added");
        model.add_element(Element2D::Truss(second_truss)).expect("second element should be added");

        let actual = assemble_stiffness_matrix(&model).expect("global matrix should be assembled");
        let expected = DMatrix::from_row_slice(
            6,
            6,
            &[
                6.0, 0.0, -6.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -6.0, 0.0, 12.0, 0.0, -6.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -6.0, 0.0, 6.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            ],
        );

        assert_matrix_approximately_equal(&actual, &expected);
    }

    #[test]
    fn assembles_mixed_elements_with_correct_global_size() {
        let mut model = Model2D::new();
        let material = Material2D::new(2.0, 0.3, 1.0).expect("valid material should be created");

        model.set_material(material);

        for (id, x, y) in [(1, 0.0, 0.0), (2, 1.0, 0.0), (3, 0.0, 1.0)] {
            model.add_node(Node2D::new(id, x, y).expect("valid node should be created")).expect("node should be added");
        }

        let beam = Beam2D::new(10, [1, 2], DEFAULT_MATERIAL_ID, 100).expect("valid beam should be created");
        let beam_section = Section2D::Beam(BeamSection2D::new(1.0, 1.0).expect("valid section should be created"));
        let triangle =
            TriangleT3::new(20, [1, 2, 3], DEFAULT_MATERIAL_ID, 200).expect("valid triangle should be created");
        let triangle_section =
            Section2D::PlaneStress(PlaneStressSection2D::new(1.0).expect("valid section should be created"));

        model.add_element_with_section(Element2D::Beam(beam), beam_section).expect("beam should be added");
        model
            .add_element_with_section(Element2D::TriangleT3(triangle), triangle_section)
            .expect("triangle should be added");

        let actual = assemble_stiffness_matrix(&model).expect("global matrix should be assembled");

        assert_eq!(actual.shape(), (8, 8));

        for row in 0..actual.nrows() {
            for column in 0..actual.ncols() {
                assert!(
                    (actual[(row, column)] - actual[(column, row)]).abs() < 1e-12,
                    "matrix is not symmetric at row {row}, column {column}"
                );
            }
        }
    }

    fn assert_matrix_approximately_equal(actual: &DMatrix<f64>, expected: &DMatrix<f64>) {
        assert_eq!(actual.shape(), expected.shape());

        for row in 0..actual.nrows() {
            for column in 0..actual.ncols() {
                assert!(
                    (actual[(row, column)] - expected[(row, column)]).abs() < 1e-12,
                    "different matrix entry at row {row}, column {column}: actual = {}, expected = {}",
                    actual[(row, column)],
                    expected[(row, column)]
                );
            }
        }
    }

    #[test]
    fn sparse_and_dense_assembly_produce_the_same_result() {
        let mut model = Model2D::new();

        let material = Material2D::new(2.0, 0.3, 1.0).expect("valid material should be created");

        model.set_material(material);

        for (id, x) in [(1, 0.0), (2, 1.0), (3, 2.0)] {
            model
                .add_node(Node2D::new(id, x, 0.0).expect("valid node should be created"))
                .expect("node should be added");
        }

        let section = Section2D::Truss(TrussSection2D::new(3.0).expect("valid section should be created"));
        model.add_section(100, section).expect("section should be added");

        let first_truss = Truss2D::new(10, [1, 2], DEFAULT_MATERIAL_ID, 100).expect("valid truss should be created");

        let second_truss = Truss2D::new(20, [2, 3], DEFAULT_MATERIAL_ID, 100).expect("valid truss should be created");

        model.add_element(Element2D::Truss(first_truss)).expect("first element should be added");

        model.add_element(Element2D::Truss(second_truss)).expect("second element should be added");

        let dense = assemble_stiffness_matrix(&model).expect("dense matrix should be assembled");

        let sparse = assemble_sparse_stiffness_matrix(&model).expect("sparse matrix should be assembled");

        let x = DVector::from_vec(vec![1.0; dense.ncols()]);
        let dense_result = &dense * &x;

        let mut sparse_result = vec![0.0; sparse.nrows()];

        sparse
            .mul_vector(x.as_slice(), &mut sparse_result)
            .expect("sparse matrix-vector multiplication should succeed");

        for (actual, expected) in sparse_result.iter().zip(dense_result.iter()) {
            assert!((actual - expected).abs() < 1e-12, "different result: sparse = {actual}, dense = {expected}");
        }
    }
}
