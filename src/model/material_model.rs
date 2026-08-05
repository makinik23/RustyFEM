//! Defines material properties used by two-dimensional FEM elements.

use crate::FemError;

/// Describes a homogeneous, isotropic, linear-elastic material.
/// The material model assumes small strains and displacements.
/// Plasticity, material damage, anisotropy, thermal effects, and
/// time-dependent behavior are not supported yet.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Material2D {
    /// Young's modulus in units consistent with the model.
    young_modulus: f64,

    /// Poisson's ratio of the material.
    poisson_ratio: f64,

    /// Mass density per unit volume.
    density: f64,
}

impl Material2D {
    /// Creates a new material after validating its physical properties.
    pub fn new(young_modulus: f64, poisson_ratio: f64, density: f64) -> Result<Self, FemError> {
        if !young_modulus.is_finite() || young_modulus <= 0.0 {
            return Err(FemError::InvalidMaterialProperty {
                property: "Young's modulus",
                value: young_modulus,
                reason: "must be finite and strictly positive",
            });
        }

        if !poisson_ratio.is_finite() || !(-1.0 < poisson_ratio && poisson_ratio < 0.5) {
            return Err(FemError::InvalidMaterialProperty {
                property: "Poisson's ratio",
                value: poisson_ratio,
                reason: "must be finite and in the open interval (-1, 0.5)",
            });
        }

        if !density.is_finite() || density < 0.0 {
            return Err(FemError::InvalidMaterialProperty {
                property: "Density",
                value: density,
                reason: "must be finite and non-negative",
            });
        }

        Ok(Self { young_modulus, poisson_ratio, density })
    }

    /// Returns the Young's modulus of the material.
    #[must_use]
    pub fn young_modulus(&self) -> f64 {
        self.young_modulus
    }

    /// Returns the Poisson's ratio of the material.
    #[must_use]
    pub fn poisson_ratio(&self) -> f64 {
        self.poisson_ratio
    }

    /// Returns the mass density of the material.
    #[must_use]
    pub fn density(&self) -> f64 {
        self.density
    }
}

#[cfg(test)]
mod tests {
    use super::Material2D;
    use crate::FemError;

    fn assert_invalid_property(
        result: Result<Material2D, FemError>, case_name: &str, expected_property: &'static str, expected_value: f64,
        expected_reason: &'static str,
    ) {
        match result {
            Err(FemError::InvalidMaterialProperty { property, value, reason }) => {
                assert_eq!(property, expected_property, "failed case: {case_name}");
                assert!(
                    value == expected_value || (value.is_nan() && expected_value.is_nan()),
                    "failed case: {case_name}"
                );
                assert_eq!(reason, expected_reason, "failed case: {case_name}");
            }
            other => panic!("failed case: {case_name}; expected an invalid material property error, got {other:?}"),
        }
    }

    #[test]
    fn creates_material_with_valid_properties() {
        let material = Material2D::new(210e9, 0.3, 7800.0).expect("valid material should be created");

        assert_eq!(material.young_modulus(), 210e9);
        assert_eq!(material.poisson_ratio(), 0.3);
        assert_eq!(material.density(), 7800.0);
    }

    #[test]
    fn rejects_invalid_young_modulus_values() {
        let cases = [("zero", 0.0), ("negative", -210e9), ("infinite", f64::INFINITY), ("not a number", f64::NAN)];

        for (name, value) in cases {
            let result = Material2D::new(value, 0.3, 7800.0);

            assert_invalid_property(result, name, "Young's modulus", value, "must be finite and strictly positive");
        }
    }

    #[test]
    fn rejects_invalid_poisson_ratio_values() {
        let cases = [
            ("lower boundary", -1.0),
            ("upper boundary", 0.5),
            ("negative infinity", f64::NEG_INFINITY),
            ("positive infinity", f64::INFINITY),
            ("not a number", f64::NAN),
        ];

        for (name, value) in cases {
            let result = Material2D::new(210e9, value, 7800.0);

            assert_invalid_property(
                result,
                name,
                "Poisson's ratio",
                value,
                "must be finite and in the open interval (-1, 0.5)",
            );
        }
    }

    #[test]
    fn rejects_invalid_density_values() {
        let cases = [
            ("negative", -7800.0),
            ("negative infinity", f64::NEG_INFINITY),
            ("positive infinity", f64::INFINITY),
            ("not a number", f64::NAN),
        ];

        for (name, value) in cases {
            let result = Material2D::new(210e9, 0.3, value);

            assert_invalid_property(result, name, "Density", value, "must be finite and non-negative");
        }
    }

    #[test]
    fn accepts_zero_density() {
        let result = Material2D::new(210e9, 0.3, 0.0);

        assert!(result.is_ok());
    }

    #[test]
    fn formats_invalid_material_property_error() {
        let error = Material2D::new(0.0, 0.3, 7800.0).expect_err("zero Young's modulus should be rejected");

        assert_eq!(error.to_string(), "material has invalid Young's modulus: 0; must be finite and strictly positive");
    }
}
