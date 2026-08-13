//! JSON-friendly input and output data transfer objects.

pub mod model_2d_input;
pub mod result_2d_output;

pub use model_2d_input::{
    AnalysisSettings2DInput, BeamSection2DInput, Dof2DInput, ElementLoad2DInput, ElementLoad2DInputKind,
    ElementType2DInput, LoadCoordinateSystem2DInput, Material2DInput, Model2DInput, NodalLoad2DInput, Node2DInput,
    PlaneStressSection2DInput, Section2DInput, Section2DInputKind, SolverKind2DInput, TrussSection2DInput,
};
pub use result_2d_output::{
    AnalysisResult2DOutput, BeamResponse2DOutput, Displacement2DOutput, ElementResponse2DOutput,
    PlaneStressResponse2DOutput, Reaction2DOutput, SolverReport2DOutput, SolverTerminationReason2DOutput,
    TrussResponse2DOutput,
};
