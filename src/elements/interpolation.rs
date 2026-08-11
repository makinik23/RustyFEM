//! Defines interpolation methods used to approximate a displacement field inside an element.

use crate::error::FemError;

/// Interpolation used to approximate a displacement field inside an element.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Interpolation {
    /// Two-node linear Lagrange interpolation.
    LinearLagrange,

    /// Two-node cubic Hermite interpolation used by beam elements.
    CubicHermite,

    /// Linear interpolation over a three-node triangle.
    LinearTriangleT3,

    /// Quadratic interpolation over a six-node triangle.
    QuadraticTriangleT6,

    /// Bilinear interpolation over a four-node quadrilateral.
    BilinearQuadQ4,

    /// Quadratic serendipity interpolation over an eight-node quadrilateral.
    SerendipityQuadQ8,
}

/// Returns the two linear Lagrange shape functions evaluated at `xi`.
///
/// The natural coordinate is usually in `[0, 1]`, but the polynomial is also
/// defined outside that interval for extrapolation.
#[must_use]
pub fn linear_lagrange_shape_functions(xi: f64) -> [f64; 2] {
    [1.0 - xi, xi]
}

/// Returns the four cubic Hermite shape functions for a beam's transverse field.
///
/// The returned functions multiply `[v1, theta1, v2, theta2]`, respectively.
/// `xi` is the normalized coordinate `x / length`.
pub fn cubic_hermite_shape_functions(xi: f64, length: f64) -> Result<[f64; 4], FemError> {
    validate_interpolation_length(length)?;

    Ok([
        1.0 - 3.0 * xi.powi(2) + 2.0 * xi.powi(3),
        length * (xi - 2.0 * xi.powi(2) + xi.powi(3)),
        3.0 * xi.powi(2) - 2.0 * xi.powi(3),
        length * (-xi.powi(2) + xi.powi(3)),
    ])
}

/// Returns the four bilinear shape functions of a Q4 quadrilateral.
///
/// The natural node order is `(-1, -1)`, `(1, -1)`, `(1, 1)`, and `(-1, 1)`.
#[must_use]
pub fn quad_q4_shape_functions(xi: f64, eta: f64) -> [f64; 4] {
    [
        0.25 * (1.0 - xi) * (1.0 - eta),
        0.25 * (1.0 + xi) * (1.0 - eta),
        0.25 * (1.0 + xi) * (1.0 + eta),
        0.25 * (1.0 - xi) * (1.0 + eta),
    ]
}

/// Returns derivatives of the Q4 shape functions with respect to `xi` and `eta`.
///
/// The first row contains `dN/dxi`; the second row contains `dN/deta`.
#[must_use]
pub fn quad_q4_shape_function_derivatives(xi: f64, eta: f64) -> [[f64; 4]; 2] {
    [
        [-0.25 * (1.0 - eta), 0.25 * (1.0 - eta), 0.25 * (1.0 + eta), -0.25 * (1.0 + eta)],
        [-0.25 * (1.0 - xi), -0.25 * (1.0 + xi), 0.25 * (1.0 + xi), 0.25 * (1.0 - xi)],
    ]
}

/// Returns the eight serendipity shape functions of a Q8 quadrilateral.
///
/// The natural node order is `(-1, -1)`, `(1, -1)`, `(1, 1)`, `(-1, 1)`,
/// `(0, -1)`, `(1, 0)`, `(0, 1)`, and `(-1, 0)`.
#[must_use]
pub fn quad_q8_shape_functions(xi: f64, eta: f64) -> [f64; 8] {
    [
        -0.25 * (1.0 - xi) * (1.0 - eta) * (1.0 + xi + eta),
        -0.25 * (1.0 + xi) * (1.0 - eta) * (1.0 - xi + eta),
        -0.25 * (1.0 + xi) * (1.0 + eta) * (1.0 - xi - eta),
        -0.25 * (1.0 - xi) * (1.0 + eta) * (1.0 + xi - eta),
        0.5 * (1.0 - xi.powi(2)) * (1.0 - eta),
        0.5 * (1.0 + xi) * (1.0 - eta.powi(2)),
        0.5 * (1.0 - xi.powi(2)) * (1.0 + eta),
        0.5 * (1.0 - xi) * (1.0 - eta.powi(2)),
    ]
}

/// Returns derivatives of the Q8 shape functions with respect to `xi` and `eta`.
///
/// The first row contains `dN/dxi`; the second row contains `dN/deta`.
#[must_use]
pub fn quad_q8_shape_function_derivatives(xi: f64, eta: f64) -> [[f64; 8]; 2] {
    [
        [
            0.25 * (1.0 - eta) * (2.0 * xi + eta),
            0.25 * (1.0 - eta) * (2.0 * xi - eta),
            0.25 * (1.0 + eta) * (2.0 * xi + eta),
            0.25 * (1.0 + eta) * (2.0 * xi - eta),
            -xi * (1.0 - eta),
            0.5 * (1.0 - eta.powi(2)),
            -xi * (1.0 + eta),
            -0.5 * (1.0 - eta.powi(2)),
        ],
        [
            0.25 * (1.0 - xi) * (xi + 2.0 * eta),
            0.25 * (1.0 + xi) * (2.0 * eta - xi),
            0.25 * (1.0 + xi) * (xi + 2.0 * eta),
            0.25 * (1.0 - xi) * (2.0 * eta - xi),
            -0.5 * (1.0 - xi.powi(2)),
            -eta * (1.0 + xi),
            0.5 * (1.0 - xi.powi(2)),
            -eta * (1.0 - xi),
        ],
    ]
}

/// Returns first derivatives with respect to the physical coordinate `x` of
/// the four cubic Hermite shape functions.
pub fn cubic_hermite_first_derivatives(xi: f64, length: f64) -> Result<[f64; 4], FemError> {
    validate_interpolation_length(length)?;

    Ok([
        (-6.0 * xi + 6.0 * xi.powi(2)) / length,
        1.0 - 4.0 * xi + 3.0 * xi.powi(2),
        (6.0 * xi - 6.0 * xi.powi(2)) / length,
        -2.0 * xi + 3.0 * xi.powi(2),
    ])
}

/// Returns second derivatives with respect to the physical coordinate `x` of
/// the four cubic Hermite shape functions.
pub fn cubic_hermite_second_derivatives(xi: f64, length: f64) -> Result<[f64; 4], FemError> {
    validate_interpolation_length(length)?;

    Ok([
        (-6.0 + 12.0 * xi) / length.powi(2),
        (-4.0 + 6.0 * xi) / length,
        (6.0 - 12.0 * xi) / length.powi(2),
        (-2.0 + 6.0 * xi) / length,
    ])
}

/// Returns the three linear shape functions of a T3 triangle.
///
/// `xi` and `eta` are the two natural coordinates. The physical triangle is
/// represented by `xi >= 0`, `eta >= 0`, and `xi + eta <= 1`.
#[must_use]
pub fn triangle_t3_shape_functions(xi: f64, eta: f64) -> [f64; 3] {
    [1.0 - xi - eta, xi, eta]
}

/// Returns the six quadratic shape functions of a T6 triangle.
///
/// The natural node order is corner `L1`, corner `L2`, corner `L3`, midside
/// `L1-L2`, midside `L2-L3`, and midside `L3-L1`, where
/// `L1 = 1 - xi - eta`, `L2 = xi`, and `L3 = eta`.
#[must_use]
pub fn triangle_t6_shape_functions(xi: f64, eta: f64) -> [f64; 6] {
    let l1 = 1.0 - xi - eta;
    let l2 = xi;
    let l3 = eta;

    [l1 * (2.0 * l1 - 1.0), l2 * (2.0 * l2 - 1.0), l3 * (2.0 * l3 - 1.0), 4.0 * l1 * l2, 4.0 * l2 * l3, 4.0 * l3 * l1]
}

/// Returns derivatives of the T6 shape functions with respect to `xi` and `eta`.
///
/// The first row contains `dN/dxi`; the second row contains `dN/deta`.
#[must_use]
pub fn triangle_t6_shape_function_derivatives(xi: f64, eta: f64) -> [[f64; 6]; 2] {
    let l1 = 1.0 - xi - eta;
    let l2 = xi;
    let l3 = eta;

    [
        [1.0 - 4.0 * l1, 4.0 * l2 - 1.0, 0.0, 4.0 * (l1 - l2), 4.0 * l3, -4.0 * l3],
        [1.0 - 4.0 * l1, 0.0, 4.0 * l3 - 1.0, -4.0 * l2, 4.0 * l2, 4.0 * (l1 - l3)],
    ]
}

fn validate_interpolation_length(length: f64) -> Result<(), FemError> {
    if !length.is_finite() || length <= 0.0 {
        return Err(FemError::InvalidInterpolationLength { value: length });
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::{
        cubic_hermite_first_derivatives, cubic_hermite_second_derivatives, cubic_hermite_shape_functions,
        linear_lagrange_shape_functions, quad_q4_shape_function_derivatives, quad_q4_shape_functions,
        quad_q8_shape_function_derivatives, quad_q8_shape_functions, triangle_t3_shape_functions,
        triangle_t6_shape_function_derivatives, triangle_t6_shape_functions,
    };
    use crate::FemError;
    use approx::assert_relative_eq;

    #[test]
    fn linear_lagrange_functions_have_partition_of_unity() {
        let functions = linear_lagrange_shape_functions(0.25);

        assert_relative_eq!(functions[0] + functions[1], 1.0, epsilon = 1e-12);
        assert_relative_eq!(functions[0], 0.75, epsilon = 1e-12);
        assert_relative_eq!(functions[1], 0.25, epsilon = 1e-12);
    }

    #[test]
    fn cubic_hermite_functions_reproduce_nodal_values() {
        let at_first_node = cubic_hermite_shape_functions(0.0, 2.0).expect("valid length");
        let at_second_node = cubic_hermite_shape_functions(1.0, 2.0).expect("valid length");

        assert_eq!(at_first_node, [1.0, 0.0, 0.0, 0.0]);
        assert_eq!(at_second_node, [0.0, 0.0, 1.0, 0.0]);
    }

    #[test]
    fn cubic_hermite_first_derivatives_reproduce_nodal_rotations() {
        let at_first_node = cubic_hermite_first_derivatives(0.0, 2.0).expect("valid length");
        let at_second_node = cubic_hermite_first_derivatives(1.0, 2.0).expect("valid length");

        assert_relative_eq!(at_first_node[0], 0.0, epsilon = 1e-12);
        assert_relative_eq!(at_first_node[1], 1.0, epsilon = 1e-12);
        assert_relative_eq!(at_first_node[2], 0.0, epsilon = 1e-12);
        assert_relative_eq!(at_first_node[3], 0.0, epsilon = 1e-12);
        assert_relative_eq!(at_second_node[0], 0.0, epsilon = 1e-12);
        assert_relative_eq!(at_second_node[1], 0.0, epsilon = 1e-12);
        assert_relative_eq!(at_second_node[2], 0.0, epsilon = 1e-12);
        assert_relative_eq!(at_second_node[3], 1.0, epsilon = 1e-12);
    }

    #[test]
    fn cubic_hermite_second_derivatives_reproduce_a_rigid_linear_field() {
        let derivatives = cubic_hermite_second_derivatives(0.5, 2.0).expect("valid length");
        let curvature = derivatives[0] * 0.0 + derivatives[1] * 1.0 + derivatives[2] * 2.0 + derivatives[3] * 1.0;

        assert_relative_eq!(curvature, 0.0, epsilon = 1e-12);
    }

    #[test]
    fn triangle_functions_have_partition_of_unity() {
        let functions = triangle_t3_shape_functions(0.2, 0.3);

        assert_relative_eq!(functions.iter().sum::<f64>(), 1.0, epsilon = 1e-12);
        assert_relative_eq!(functions[0], 0.5, epsilon = 1e-12);
        assert_relative_eq!(functions[1], 0.2, epsilon = 1e-12);
        assert_relative_eq!(functions[2], 0.3, epsilon = 1e-12);
    }

    #[test]
    fn triangle_t6_functions_have_partition_of_unity() {
        let functions = triangle_t6_shape_functions(0.2, 0.3);

        assert_relative_eq!(functions.iter().sum::<f64>(), 1.0, epsilon = 1e-12);
    }

    #[test]
    fn triangle_t6_functions_reproduce_nodal_values() {
        let natural_nodes = [(0.0, 0.0), (1.0, 0.0), (0.0, 1.0), (0.5, 0.0), (0.5, 0.5), (0.0, 0.5)];

        for (node_index, (xi, eta)) in natural_nodes.into_iter().enumerate() {
            let functions = triangle_t6_shape_functions(xi, eta);

            for (function_index, function) in functions.into_iter().enumerate() {
                let expected = if node_index == function_index { 1.0 } else { 0.0 };

                assert_relative_eq!(function, expected, epsilon = 1e-12);
            }
        }
    }

    #[test]
    fn triangle_t6_derivatives_sum_to_zero() {
        let derivatives = triangle_t6_shape_function_derivatives(0.2, 0.3);

        assert_relative_eq!(derivatives[0].iter().sum::<f64>(), 0.0, epsilon = 1e-12);
        assert_relative_eq!(derivatives[1].iter().sum::<f64>(), 0.0, epsilon = 1e-12);
    }

    #[test]
    fn quad_q4_functions_have_partition_of_unity() {
        let functions = quad_q4_shape_functions(0.2, -0.4);

        assert_relative_eq!(functions.iter().sum::<f64>(), 1.0, epsilon = 1e-12);
    }

    #[test]
    fn quad_q4_functions_reproduce_nodal_values() {
        let natural_nodes = [(-1.0, -1.0), (1.0, -1.0), (1.0, 1.0), (-1.0, 1.0)];

        for (node_index, (xi, eta)) in natural_nodes.into_iter().enumerate() {
            let functions = quad_q4_shape_functions(xi, eta);

            for (function_index, function) in functions.into_iter().enumerate() {
                let expected = if node_index == function_index { 1.0 } else { 0.0 };

                assert_relative_eq!(function, expected, epsilon = 1e-12);
            }
        }
    }

    #[test]
    fn quad_q4_derivatives_match_center_values() {
        let derivatives = quad_q4_shape_function_derivatives(0.0, 0.0);

        assert_eq!(derivatives[0], [-0.25, 0.25, 0.25, -0.25]);
        assert_eq!(derivatives[1], [-0.25, -0.25, 0.25, 0.25]);
    }

    #[test]
    fn quad_q4_derivatives_sum_to_zero() {
        let derivatives = quad_q4_shape_function_derivatives(0.2, -0.4);

        assert_relative_eq!(derivatives[0].iter().sum::<f64>(), 0.0, epsilon = 1e-12);
        assert_relative_eq!(derivatives[1].iter().sum::<f64>(), 0.0, epsilon = 1e-12);
    }

    #[test]
    fn quad_q8_functions_have_partition_of_unity() {
        let functions = quad_q8_shape_functions(0.2, -0.4);

        assert_relative_eq!(functions.iter().sum::<f64>(), 1.0, epsilon = 1e-12);
    }

    #[test]
    fn quad_q8_functions_reproduce_nodal_values() {
        let natural_nodes =
            [(-1.0, -1.0), (1.0, -1.0), (1.0, 1.0), (-1.0, 1.0), (0.0, -1.0), (1.0, 0.0), (0.0, 1.0), (-1.0, 0.0)];

        for (node_index, (xi, eta)) in natural_nodes.into_iter().enumerate() {
            let functions = quad_q8_shape_functions(xi, eta);

            for (function_index, function) in functions.into_iter().enumerate() {
                let expected = if node_index == function_index { 1.0 } else { 0.0 };

                assert_relative_eq!(function, expected, epsilon = 1e-12);
            }
        }
    }

    #[test]
    fn quad_q8_derivatives_match_center_values() {
        let derivatives = quad_q8_shape_function_derivatives(0.0, 0.0);

        assert_eq!(derivatives[0], [0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.0, -0.5]);
        assert_eq!(derivatives[1], [0.0, 0.0, 0.0, 0.0, -0.5, 0.0, 0.5, 0.0]);
    }

    #[test]
    fn quad_q8_derivatives_sum_to_zero() {
        let derivatives = quad_q8_shape_function_derivatives(0.2, -0.4);

        assert_relative_eq!(derivatives[0].iter().sum::<f64>(), 0.0, epsilon = 1e-12);
        assert_relative_eq!(derivatives[1].iter().sum::<f64>(), 0.0, epsilon = 1e-12);
    }

    #[test]
    fn quad_q8_derivatives_match_finite_differences() {
        let xi = 0.17;
        let eta = -0.31;
        let step = 1e-6;
        let derivatives = quad_q8_shape_function_derivatives(xi, eta);
        let xi_plus = quad_q8_shape_functions(xi + step, eta);
        let xi_minus = quad_q8_shape_functions(xi - step, eta);
        let eta_plus = quad_q8_shape_functions(xi, eta + step);
        let eta_minus = quad_q8_shape_functions(xi, eta - step);

        for index in 0..8 {
            let finite_difference_xi = (xi_plus[index] - xi_minus[index]) / (2.0 * step);
            let finite_difference_eta = (eta_plus[index] - eta_minus[index]) / (2.0 * step);

            assert_relative_eq!(derivatives[0][index], finite_difference_xi, epsilon = 1e-9);
            assert_relative_eq!(derivatives[1][index], finite_difference_eta, epsilon = 1e-9);
        }
    }

    #[test]
    fn rejects_invalid_hermite_length() {
        let cases = [0.0, -1.0, f64::INFINITY, f64::NAN];

        for length in cases {
            let result = cubic_hermite_shape_functions(0.5, length);

            assert!(
                matches!(result, Err(FemError::InvalidInterpolationLength { value }) if value == length || (value.is_nan() && length.is_nan()))
            );
        }
    }
}
