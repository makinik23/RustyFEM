//! Defines interpolation methods used to approximate a displacement field inside an element.

/// Interpolation used to approximate a displacement field inside an element.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Interpolation {
    /// Two-node linear Lagrange interpolation.
    LinearLagrange,

    /// Two-node cubic Hermite interpolation used by beam elements.
    CubicHermite,

    /// Linear interpolation over a three-node triangle.
    LinearTriangleT3,
}
