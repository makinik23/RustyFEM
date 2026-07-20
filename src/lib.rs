//! Educational finite element method solver.
//!
//! The first implemented domain is a one-dimensional axial bar model. The
//! current crate exposes only validated domain types; assembly and solving are
//! intentionally left for later, smaller increments.

pub mod elements;
pub mod error;
pub mod math;
pub mod model;

pub use error::FemError;
