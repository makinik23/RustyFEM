//! Serializable output format for solved `Model2D` analyses.

use serde::{Deserialize, Serialize};

use crate::FemError;
use crate::analysis::iterative_solver::CgTerminationReason;
use crate::analysis::solver::{AnalysisResult2D, SolverReport};
use crate::analysis::{ElementResponse2D, recover_model_responses};
use crate::model::{Dof2D, DofNumbering2D, Model2D};

/// JSON-friendly solved-result payload for a two-dimensional model.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct AnalysisResult2DOutput {
    pub displacements: Vec<Displacement2DOutput>,
    pub reactions: Vec<Reaction2DOutput>,
    pub solver_report: Option<SolverReport2DOutput>,
    pub element_responses: Vec<ElementResponse2DOutput>,
}

/// One displacement component in global node/DOF notation.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Displacement2DOutput {
    pub node: usize,
    pub dof: String,
    pub value: f64,
}

/// One reaction component in global node/DOF notation.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Reaction2DOutput {
    pub node: usize,
    pub dof: String,
    pub value: f64,
}

/// Sparse solver diagnostics in a serializable shape.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct SolverReport2DOutput {
    pub iterations: usize,
    pub residual_norm: f64,
    pub relative_residual_norm: f64,
    pub termination_reason: SolverTerminationReason2DOutput,
}

/// Serialized sparse solver termination reason.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum SolverTerminationReason2DOutput {
    Converged,
    MaxIterations,
    Stagnated,
}

/// One recovered element response.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(tag = "type", rename_all = "snake_case")]
pub enum ElementResponse2DOutput {
    Truss { element: usize, response: TrussResponse2DOutput },
    Beam { element: usize, response: BeamResponse2DOutput },
    PlaneStress { element: usize, response: PlaneStressResponse2DOutput },
}

/// Serialized truss response.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct TrussResponse2DOutput {
    pub strain: f64,
    pub stress: f64,
    pub axial_force: f64,
}

/// Serialized beam end-force response.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct BeamResponse2DOutput {
    pub end_forces: [f64; 6],
}

/// Serialized plane-stress response at one recovery point.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct PlaneStressResponse2DOutput {
    pub natural_coordinates: Option<[f64; 2]>,
    pub strain: [f64; 3],
    pub stress: [f64; 3],
    pub von_mises: f64,
}

impl AnalysisResult2DOutput {
    /// Creates a serializable result payload from a solved model.
    pub fn from_model_and_result(model: &Model2D, result: &AnalysisResult2D) -> Result<Self, FemError> {
        let numbering = DofNumbering2D::from_model(model)?;
        let mut displacements = Vec::new();
        let mut reactions = Vec::new();

        for node in model.nodes() {
            for dof in [Dof2D::Ux, Dof2D::Uy, Dof2D::Rz] {
                if let Ok(index) = numbering.index(node.id(), dof) {
                    displacements.push(Displacement2DOutput {
                        node: node.id(),
                        dof: dof.name().to_owned(),
                        value: result.displacements()[index],
                    });
                    reactions.push(Reaction2DOutput {
                        node: node.id(),
                        dof: dof.name().to_owned(),
                        value: result.reactions()[index],
                    });
                }
            }
        }

        let element_responses = recover_model_responses(model, result.displacements())?
            .into_iter()
            .map(ElementResponse2DOutput::from)
            .collect();

        Ok(Self {
            displacements,
            reactions,
            solver_report: result.solver_report().copied().map(SolverReport2DOutput::from),
            element_responses,
        })
    }
}

impl From<SolverReport> for SolverReport2DOutput {
    fn from(value: SolverReport) -> Self {
        Self {
            iterations: value.iterations,
            residual_norm: value.residual_norm,
            relative_residual_norm: value.relative_residual_norm,
            termination_reason: value.termination_reason.into(),
        }
    }
}

impl From<CgTerminationReason> for SolverTerminationReason2DOutput {
    fn from(value: CgTerminationReason) -> Self {
        match value {
            CgTerminationReason::Converged => Self::Converged,
            CgTerminationReason::MaxIterations => Self::MaxIterations,
            CgTerminationReason::Stagnated => Self::Stagnated,
        }
    }
}

impl From<(usize, ElementResponse2D)> for ElementResponse2DOutput {
    fn from((element, response): (usize, ElementResponse2D)) -> Self {
        match response {
            ElementResponse2D::Truss(response) => Self::Truss {
                element,
                response: TrussResponse2DOutput {
                    strain: response.strain(),
                    stress: response.stress(),
                    axial_force: response.axial_force(),
                },
            },
            ElementResponse2D::Beam(response) => {
                Self::Beam { element, response: BeamResponse2DOutput { end_forces: *response.end_forces() } }
            }
            ElementResponse2D::Triangle(response) => Self::PlaneStress {
                element,
                response: PlaneStressResponse2DOutput {
                    natural_coordinates: None,
                    strain: *response.strain(),
                    stress: *response.stress(),
                    von_mises: response.von_mises_stress(),
                },
            },
            ElementResponse2D::Quadrilateral(response) => {
                let (xi, eta) = response.natural_coordinates();

                Self::PlaneStress {
                    element,
                    response: PlaneStressResponse2DOutput {
                        natural_coordinates: Some([xi, eta]),
                        strain: *response.strain(),
                        stress: *response.stress(),
                        von_mises: response.von_mises_stress(),
                    },
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::{AnalysisResult2DOutput, ElementResponse2DOutput};
    use crate::analysis::solver::solve;
    use crate::elements::{Element2D, Truss2D};
    use crate::model::{
        DEFAULT_MATERIAL_ID, DisplacementConstraint2D, Dof2D, Material2D, Model2D, NodalLoad2D, Section2D,
        TrussSection2D,
    };

    #[test]
    fn creates_serializable_output_from_solved_model() {
        let mut model = Model2D::new();
        model.set_material(Material2D::new(200.0, 0.3, 1.0).expect("valid material"));
        model.add_node(crate::model::Node2D::new(1, 0.0, 0.0).expect("valid node")).expect("node should be added");
        model.add_node(crate::model::Node2D::new(2, 1.0, 0.0).expect("valid node")).expect("node should be added");
        model
            .add_element_with_section(
                Element2D::Truss(Truss2D::new(1, [1, 2], DEFAULT_MATERIAL_ID, 1).expect("valid truss")),
                Section2D::Truss(TrussSection2D::new(1.0).expect("valid section")),
            )
            .expect("element should be added");
        model
            .add_constraint(DisplacementConstraint2D::new(1, Dof2D::Ux, 0.0).expect("valid constraint"))
            .expect("constraint should be added");
        model
            .add_constraint(DisplacementConstraint2D::new(1, Dof2D::Uy, 0.0).expect("valid constraint"))
            .expect("constraint should be added");
        model
            .add_constraint(DisplacementConstraint2D::new(2, Dof2D::Uy, 0.0).expect("valid constraint"))
            .expect("constraint should be added");
        model.add_load(NodalLoad2D::new(2, Dof2D::Ux, 10.0).expect("valid load")).expect("load should be added");

        let result = solve(&model).expect("model should solve");
        let output = AnalysisResult2DOutput::from_model_and_result(&model, &result).expect("output should be created");
        let serialized = serde_json::to_string(&output).expect("output should serialize");

        assert!(serialized.contains(r#""displacements""#));
        assert_eq!(output.displacements.len(), 4);
        assert_eq!(output.reactions.len(), 4);
        assert!(matches!(output.element_responses[0], ElementResponse2DOutput::Truss { .. }));
    }
}
