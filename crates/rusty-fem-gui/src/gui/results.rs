//! Solver execution and compact postprocessing state for the GUI.

use rusty_fem::analysis::iterative_solver::{CgProgress, CgTerminationReason};
use rusty_fem::analysis::solver::solve_with_settings_and_progress;
use rusty_fem::analysis::{NodalPlaneStressResponse2D, recover_nodal_plane_stress_responses_with_progress};
use rusty_fem::elements::{Element2D, quad_q4_shape_functions, quad_q8_shape_functions, triangle_t6_shape_functions};
use rusty_fem::model::{Dof2D, DofNumbering2D, Model2D, SolverKind2D};
use std::collections::HashMap;

const T6_CONTOUR_SUBDIVISIONS: usize = 3;
const Q8_CONTOUR_SUBDIVISIONS: usize = 2;

/// Result layer drawn on the model canvas.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum ResultField {
    Model,
    Displacement,
    StressX,
    StressY,
    ShearStress,
    VonMisesStress,
    StrainX,
    StrainY,
    ShearStrain,
    EquivalentStrain,
}

impl ResultField {
    pub(super) const ALL: [Self; 10] = [
        Self::Model,
        Self::Displacement,
        Self::StressX,
        Self::StressY,
        Self::ShearStress,
        Self::VonMisesStress,
        Self::StrainX,
        Self::StrainY,
        Self::ShearStrain,
        Self::EquivalentStrain,
    ];

    pub(super) fn label(self) -> &'static str {
        match self {
            Self::Model => "Model",
            Self::Displacement => "Displacement",
            Self::StressX => "sigma x",
            Self::StressY => "sigma y",
            Self::ShearStress => "tau xy",
            Self::VonMisesStress => "von Mises stress",
            Self::StrainX => "epsilon x",
            Self::StrainY => "epsilon y",
            Self::ShearStrain => "gamma xy",
            Self::EquivalentStrain => "equivalent strain",
        }
    }

    pub(super) fn is_scalar_contour(self) -> bool {
        !matches!(self, Self::Model | Self::Displacement)
    }

    pub(super) fn is_signed(self) -> bool {
        matches!(
            self,
            Self::StressX | Self::StressY | Self::ShearStress | Self::StrainX | Self::StrainY | Self::ShearStrain
        )
    }
}

/// Long-running analysis stage reported to the GUI.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum AnalysisPhase {
    Assembling,
    Solving,
    Recovering,
}

impl AnalysisPhase {
    pub(super) fn label(self) -> &'static str {
        match self {
            Self::Assembling => "Assembling system",
            Self::Solving => "Solving",
            Self::Recovering => "Recovering results",
        }
    }
}

/// Progress message produced by a background analysis.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(super) enum AnalysisProgress {
    Phase(AnalysisPhase),
    Iteration(CgProgress),
    Recovery { completed: usize, total: usize },
}

#[derive(Debug, Clone, Copy, Default)]
pub(super) struct PlaneStressState {
    strain: [f64; 3],
    stress: [f64; 3],
    equivalent_strain: f64,
}

impl PlaneStressState {
    fn from_nodal(response: &NodalPlaneStressResponse2D) -> Self {
        Self { strain: *response.strain(), stress: *response.stress(), equivalent_strain: response.equivalent_strain() }
    }

    fn weighted(states: &[Self], weights: &[f64]) -> Self {
        let mut result = Self::default();

        for (state, weight) in states.iter().zip(weights) {
            for component in 0..3 {
                result.strain[component] += weight * state.strain[component];
                result.stress[component] += weight * state.stress[component];
            }
            result.equivalent_strain += weight * state.equivalent_strain;
        }

        result
    }

    fn mean(states: &[Self]) -> Self {
        let weight = 1.0 / states.len() as f64;
        Self::weighted(states, &vec![weight; states.len()])
    }

    fn von_mises_stress(self) -> f64 {
        let [sigma_x, sigma_y, tau_xy] = self.stress;

        (sigma_x.powi(2) - sigma_x * sigma_y + sigma_y.powi(2) + 3.0 * tau_xy.powi(2)).sqrt()
    }

    pub(super) fn scalar(self, field: ResultField) -> Option<f64> {
        match field {
            ResultField::StressX => Some(self.stress[0]),
            ResultField::StressY => Some(self.stress[1]),
            ResultField::ShearStress => Some(self.stress[2]),
            ResultField::VonMisesStress => Some(self.von_mises_stress()),
            ResultField::StrainX => Some(self.strain[0]),
            ResultField::StrainY => Some(self.strain[1]),
            ResultField::ShearStrain => Some(self.strain[2]),
            ResultField::EquivalentStrain => Some(self.equivalent_strain),
            ResultField::Model | ResultField::Displacement => None,
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub(super) struct ContourVertex2D {
    pub(super) position: [f64; 2],
    pub(super) state: PlaneStressState,
}

#[derive(Debug, Clone, Copy)]
pub(super) struct ContourTriangle2D {
    pub(super) vertices: [ContourVertex2D; 3],
}

#[derive(Debug, Clone, Copy)]
pub(super) struct ScalarRange {
    pub(super) minimum: f64,
    pub(super) maximum: f64,
}

impl ScalarRange {
    pub(super) fn normalized(self, value: f64) -> f32 {
        ((value - self.minimum) / (self.maximum - self.minimum)).clamp(0.0, 1.0) as f32
    }
}

/// Data needed to draw and summarize one completed analysis.
pub(super) struct AnalysisResults {
    nodal_displacements: HashMap<usize, [f64; 2]>,
    nodal_plane_stress: HashMap<usize, PlaneStressState>,
    element_plane_stress: HashMap<usize, PlaneStressState>,
    contour_triangles: Vec<ContourTriangle2D>,
    pub(super) max_displacement: f64,
    pub(super) max_von_mises_stress: f64,
    pub(super) max_equivalent_strain: f64,
    pub(super) iterations: Option<usize>,
    pub(super) relative_residual: Option<f64>,
    pub(super) termination: Option<&'static str>,
}

impl AnalysisResults {
    pub(super) fn solve_with_progress<F>(model: &Model2D, mut on_progress: F) -> Result<Self, String>
    where
        F: FnMut(AnalysisProgress),
    {
        on_progress(AnalysisProgress::Phase(AnalysisPhase::Assembling));
        let mut solving_reported = model.analysis_settings().solver() == SolverKind2D::Dense;

        if solving_reported {
            on_progress(AnalysisProgress::Phase(AnalysisPhase::Solving));
        }

        let analysis = solve_with_settings_and_progress(model, |sample| {
            if !solving_reported {
                on_progress(AnalysisProgress::Phase(AnalysisPhase::Solving));
                solving_reported = true;
            }
            on_progress(AnalysisProgress::Iteration(sample));
        })
        .map_err(|error| error.to_string())?;
        on_progress(AnalysisProgress::Phase(AnalysisPhase::Recovering));
        let numbering = DofNumbering2D::from_model(model).map_err(|error| error.to_string())?;
        let nodal_displacements = model
            .nodes()
            .iter()
            .map(|node| {
                let ux =
                    numbering.index(node.id(), Dof2D::Ux).map(|index| analysis.displacements()[index]).unwrap_or(0.0);
                let uy =
                    numbering.index(node.id(), Dof2D::Uy).map(|index| analysis.displacements()[index]).unwrap_or(0.0);
                (node.id(), [ux, uy])
            })
            .collect::<HashMap<_, _>>();
        let max_displacement = nodal_displacements.values().map(|[ux, uy]| ux.hypot(*uy)).fold(0.0_f64, f64::max);

        let nodal_plane_stress =
            recover_nodal_plane_stress_responses_with_progress(model, analysis.displacements(), |completed, total| {
                on_progress(AnalysisProgress::Recovery { completed, total })
            })
            .map_err(|error| error.to_string())?
            .into_iter()
            .map(|response| (response.node_id(), PlaneStressState::from_nodal(&response)))
            .collect::<HashMap<_, _>>();
        let element_plane_stress = build_element_states(model, &nodal_plane_stress);
        let contour_triangles = build_contour_triangles(model, &nodal_plane_stress);
        let max_von_mises_stress =
            nodal_plane_stress.values().copied().map(PlaneStressState::von_mises_stress).fold(0.0_f64, f64::max);
        let max_equivalent_strain =
            nodal_plane_stress.values().map(|state| state.equivalent_strain).fold(0.0_f64, f64::max);
        let (iterations, relative_residual, termination) =
            analysis.solver_report().map_or((None, None, None), |report| {
                let termination = match report.termination_reason {
                    CgTerminationReason::Converged => "converged",
                    CgTerminationReason::MaxIterations => "maximum iterations",
                    CgTerminationReason::Stagnated => "stagnated",
                };

                (Some(report.iterations), Some(report.relative_residual_norm), Some(termination))
            });

        Ok(Self {
            nodal_displacements,
            nodal_plane_stress,
            element_plane_stress,
            contour_triangles,
            max_displacement,
            max_von_mises_stress,
            max_equivalent_strain,
            iterations,
            relative_residual,
            termination,
        })
    }

    pub(super) fn displacement(&self, node_id: usize) -> [f64; 2] {
        self.nodal_displacements.get(&node_id).copied().unwrap_or([0.0; 2])
    }

    pub(super) fn nodal_scalar_value(&self, field: ResultField, node_id: usize) -> Option<f64> {
        self.nodal_plane_stress.get(&node_id).and_then(|state| state.scalar(field))
    }

    pub(super) fn scalar_value(&self, field: ResultField, element_id: usize) -> Option<f64> {
        self.element_plane_stress.get(&element_id).and_then(|state| state.scalar(field))
    }

    pub(super) fn von_mises_stress(&self, element_id: usize) -> Option<f64> {
        self.scalar_value(ResultField::VonMisesStress, element_id)
    }

    pub(super) fn contour_triangles(&self) -> &[ContourTriangle2D] {
        &self.contour_triangles
    }

    pub(super) fn scalar_range(&self, field: ResultField) -> Option<ScalarRange> {
        let mut minimum = f64::INFINITY;
        let mut maximum = f64::NEG_INFINITY;

        for value in self
            .contour_triangles
            .iter()
            .flat_map(|triangle| triangle.vertices)
            .filter_map(|vertex| vertex.state.scalar(field))
            .filter(|value| value.is_finite())
        {
            minimum = minimum.min(value);
            maximum = maximum.max(value);
        }

        if !minimum.is_finite() || !maximum.is_finite() {
            return None;
        }
        if !field.is_signed() {
            minimum = 0.0;
        }
        if (maximum - minimum).abs() <= f64::EPSILON {
            let padding = maximum.abs().max(1.0) * 1e-6;
            minimum -= padding;
            maximum += padding;
        }

        Some(ScalarRange { minimum, maximum })
    }
}

fn build_element_states(
    model: &Model2D, nodal_states: &HashMap<usize, PlaneStressState>,
) -> HashMap<usize, PlaneStressState> {
    model
        .elements()
        .iter()
        .filter_map(|element| {
            if matches!(element, Element2D::Truss(_) | Element2D::Beam(_)) {
                return None;
            }

            let states =
                element.node_ids().iter().filter_map(|node_id| nodal_states.get(node_id).copied()).collect::<Vec<_>>();
            (!states.is_empty()).then(|| (element.id(), PlaneStressState::mean(&states)))
        })
        .collect()
}

fn build_contour_triangles(model: &Model2D, nodal_states: &HashMap<usize, PlaneStressState>) -> Vec<ContourTriangle2D> {
    let positions = model.nodes().iter().map(|node| (node.id(), [node.x(), node.y()])).collect::<HashMap<_, _>>();
    let mut triangles = Vec::new();

    for element in model.elements() {
        match element {
            Element2D::TriangleT3(_) => {
                let Some((positions, states)) = element_nodal_data::<3>(element, &positions, nodal_states) else {
                    continue;
                };
                triangles.push(ContourTriangle2D {
                    vertices: std::array::from_fn(|index| ContourVertex2D {
                        position: positions[index],
                        state: states[index],
                    }),
                });
            }
            Element2D::TriangleT6(_) => {
                let Some((positions, states)) = element_nodal_data::<6>(element, &positions, nodal_states) else {
                    continue;
                };

                for natural_triangle in triangle_subdivision(T6_CONTOUR_SUBDIVISIONS) {
                    triangles.push(ContourTriangle2D {
                        vertices: natural_triangle.map(|[xi, eta]| {
                            let weights = triangle_t6_shape_functions(xi, eta);
                            ContourVertex2D {
                                position: weighted_position(&positions, &weights),
                                state: PlaneStressState::weighted(&states, &weights),
                            }
                        }),
                    });
                }
            }
            Element2D::QuadQ4(_) => {
                let Some((positions, states)) = element_nodal_data::<4>(element, &positions, nodal_states) else {
                    continue;
                };
                append_quadrilateral_contours(&mut triangles, &positions, &states, 1, quad_q4_shape_functions);
            }
            Element2D::QuadQ8(_) => {
                let Some((positions, states)) = element_nodal_data::<8>(element, &positions, nodal_states) else {
                    continue;
                };
                append_quadrilateral_contours(
                    &mut triangles,
                    &positions,
                    &states,
                    Q8_CONTOUR_SUBDIVISIONS,
                    quad_q8_shape_functions,
                );
            }
            Element2D::Truss(_) | Element2D::Beam(_) => {}
        }
    }

    triangles
}

fn element_nodal_data<const N: usize>(
    element: &Element2D, positions: &HashMap<usize, [f64; 2]>, states: &HashMap<usize, PlaneStressState>,
) -> Option<([[f64; 2]; N], [PlaneStressState; N])> {
    let node_positions = element
        .node_ids()
        .iter()
        .take(N)
        .map(|node_id| positions.get(node_id).copied())
        .collect::<Option<Vec<_>>>()?
        .try_into()
        .ok()?;
    let node_states = element
        .node_ids()
        .iter()
        .take(N)
        .map(|node_id| states.get(node_id).copied())
        .collect::<Option<Vec<_>>>()?
        .try_into()
        .ok()?;

    Some((node_positions, node_states))
}

fn triangle_subdivision(subdivisions: usize) -> Vec<[[f64; 2]; 3]> {
    let mut triangles = Vec::with_capacity(subdivisions * subdivisions);
    let step = 1.0 / subdivisions as f64;

    for row in 0..subdivisions {
        for column in 0..(subdivisions - row) {
            let lower_left = [column as f64 * step, row as f64 * step];
            let lower_right = [(column + 1) as f64 * step, row as f64 * step];
            let upper_left = [column as f64 * step, (row + 1) as f64 * step];
            triangles.push([lower_left, lower_right, upper_left]);

            if column + row + 1 < subdivisions {
                let upper_right = [(column + 1) as f64 * step, (row + 1) as f64 * step];
                triangles.push([lower_right, upper_right, upper_left]);
            }
        }
    }

    triangles
}

fn append_quadrilateral_contours<const N: usize, F>(
    triangles: &mut Vec<ContourTriangle2D>, positions: &[[f64; 2]; N], states: &[PlaneStressState; N],
    subdivisions: usize, shape_functions: F,
) where
    F: Fn(f64, f64) -> [f64; N],
{
    let step = 2.0 / subdivisions as f64;

    for row in 0..subdivisions {
        for column in 0..subdivisions {
            let xi = -1.0 + column as f64 * step;
            let eta = -1.0 + row as f64 * step;
            let natural_triangles = [
                [[xi, eta], [xi + step, eta], [xi + step, eta + step]],
                [[xi, eta], [xi + step, eta + step], [xi, eta + step]],
            ];

            for natural_triangle in natural_triangles {
                triangles.push(ContourTriangle2D {
                    vertices: natural_triangle.map(|[xi, eta]| {
                        let weights = shape_functions(xi, eta);
                        ContourVertex2D {
                            position: weighted_position(positions, &weights),
                            state: PlaneStressState::weighted(states, &weights),
                        }
                    }),
                });
            }
        }
    }
}

fn weighted_position<const N: usize>(positions: &[[f64; 2]; N], weights: &[f64; N]) -> [f64; 2] {
    let mut result = [0.0; 2];

    for (position, weight) in positions.iter().zip(weights) {
        result[0] += weight * position[0];
        result[1] += weight * position[1];
    }

    result
}

#[cfg(test)]
mod tests {
    use super::{AnalysisResults, ResultField, triangle_subdivision};
    use rusty_fem::io::Model2DInput;

    #[test]
    fn solves_example_and_prepares_displacement_stress_and_strain_layers() {
        let input: Model2DInput = serde_json::from_str(include_str!("../../../../examples/t3_cantilever.json"))
            .expect("example input should deserialize");
        let model = input.into_model().expect("example model should build");

        let results = AnalysisResults::solve_with_progress(&model, |_| {}).expect("GUI analysis should complete");

        assert!(results.max_displacement > 0.0);
        assert!(results.max_von_mises_stress > 0.0);
        assert!(results.max_equivalent_strain > 0.0);
        assert!(results.scalar_range(ResultField::StressX).is_some());
        assert!(results.scalar_range(ResultField::StrainY).is_some());
        assert!(results.nodal_scalar_value(ResultField::VonMisesStress, 1).is_some());
        assert!(results.von_mises_stress(100).is_some());
        assert!(results.scalar_value(ResultField::EquivalentStrain, 100).is_some());
        assert_eq!(results.displacement(1), [0.0, 0.0]);
    }

    #[test]
    fn t6_contour_subdivision_has_n_squared_triangles() {
        assert_eq!(triangle_subdivision(1).len(), 1);
        assert_eq!(triangle_subdivision(3).len(), 9);
        assert_eq!(triangle_subdivision(4).len(), 16);
    }
}
