//! Benchmarks dense LU and sparse preconditioned CG on regular T3 meshes.
//! Run with: cargo bench --bench solver_benchmark

use std::hint::black_box;

use criterion::{BenchmarkId, Criterion, Throughput, criterion_group, criterion_main};
use rusty_fem::analysis::solver::{solve, solve_sparse};
use rusty_fem::elements::{Element2D, TriangleT3};
use rusty_fem::model::{
    DEFAULT_MATERIAL_ID, DisplacementConstraint2D, Dof2D, Material2D, Model2D, NodalLoad2D, Node2D,
    PlaneStressSection2D, Section2D,
};

/// Benchmarks a complete analysis for several rectangular T3 meshes.
fn benchmark_solver_scaling(criterion: &mut Criterion) {
    let mut group = criterion.benchmark_group("t3_cantilever_solver");

    for &(nx, ny) in &[(4, 2), (8, 4), (16, 8), (32, 16), (64, 32)] {
        let model = rectangular_t3_cantilever(nx, ny);
        let element_count = 2 * nx * ny;
        let mesh_name = format!("{nx}x{ny}");

        group.throughput(Throughput::Elements(element_count as u64));

        group.bench_with_input(BenchmarkId::new("dense_lu", &mesh_name), &model, |benchmark, model| {
            benchmark.iter(|| {
                let result = solve(black_box(model));

                black_box(result.expect("dense benchmark model should be solvable"));
            });
        });

        group.bench_with_input(BenchmarkId::new("sparse_pcg_jacobi", &mesh_name), &model, |benchmark, model| {
            benchmark.iter(|| {
                let result = solve_sparse(black_box(model));

                black_box(result.expect("sparse benchmark model should be solvable"));
            });
        });
    }

    group.finish();
}

/// Creates a rectangular cantilever made from two T3 elements per cell.
///
/// The left edge is fixed in both translational degrees of freedom. A unit
/// vertical load is applied at the upper-right node.
fn rectangular_t3_cantilever(nx: usize, ny: usize) -> Model2D {
    assert!(nx > 0, "the mesh must contain at least one cell in x");
    assert!(ny > 0, "the mesh must contain at least one cell in y");

    let mut model = Model2D::new();
    let material = Material2D::new(200.0, 0.3, 1.0).expect("benchmark material should be valid");

    model.set_material(material);

    let node_id = |i: usize, j: usize| -> usize { j * (nx + 1) + i + 1 };

    for j in 0..=ny {
        for i in 0..=nx {
            let node = Node2D::new(i + j * (nx + 1) + 1, i as f64, j as f64).expect("benchmark node should be valid");

            model.add_node(node).expect("benchmark node should be added");
        }
    }

    let mut element_id = 1;

    for j in 0..ny {
        for i in 0..nx {
            let lower_left = node_id(i, j);
            let lower_right = node_id(i + 1, j);
            let upper_right = node_id(i + 1, j + 1);
            let upper_left = node_id(i, j + 1);

            let first_triangle =
                TriangleT3::new(element_id, [lower_left, lower_right, upper_right], DEFAULT_MATERIAL_ID, element_id)
                    .expect("benchmark triangle should be valid");
            let first_section =
                Section2D::PlaneStress(PlaneStressSection2D::new(1.0).expect("benchmark section should be valid"));

            model
                .add_element_with_section(Element2D::TriangleT3(first_triangle), first_section)
                .expect("benchmark triangle should be added");

            element_id += 1;

            let second_triangle =
                TriangleT3::new(element_id, [lower_left, upper_right, upper_left], DEFAULT_MATERIAL_ID, element_id)
                    .expect("benchmark triangle should be valid");
            let second_section =
                Section2D::PlaneStress(PlaneStressSection2D::new(1.0).expect("benchmark section should be valid"));

            model
                .add_element_with_section(Element2D::TriangleT3(second_triangle), second_section)
                .expect("benchmark triangle should be added");

            element_id += 1;
        }
    }

    for j in 0..=ny {
        let node = node_id(0, j);

        for dof in [Dof2D::Ux, Dof2D::Uy] {
            let constraint =
                DisplacementConstraint2D::new(node, dof, 0.0).expect("benchmark constraint should be valid");

            model.add_constraint(constraint).expect("benchmark constraint should be added");
        }
    }

    let loaded_node = node_id(nx, ny);
    let load = NodalLoad2D::new(loaded_node, Dof2D::Uy, -1.0).expect("benchmark load should be valid");

    model.add_load(load).expect("benchmark load should be added");

    model
}

criterion_group!(benches, benchmark_solver_scaling);
criterion_main!(benches);
