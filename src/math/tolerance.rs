//! Tolerance values used by validation and numerical tests.

/// Absolute tolerance used for simple geometry checks in SI units.
pub const DEFAULT_GEOMETRY_TOLERANCE: f64 = 1.0e-12;

/// Returns `true` when `value` is finite and close enough to zero.
#[must_use]
pub fn is_near_zero(value: f64, tolerance: f64) -> bool {
    value.is_finite() && value.abs() <= tolerance
}
