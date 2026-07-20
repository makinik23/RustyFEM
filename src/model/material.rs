use serde::Serialize;

use crate::FemError;

/// Isotropic linear-elastic material in SI units.
#[derive(Debug, Clone, Copy, PartialEq, Serialize)]
pub struct LinearElasticMaterial {
    id: usize,
    young_modulus: f64,
    poisson_ratio: f64,
    density: f64,
}

impl LinearElasticMaterial {
    /// Creates a material with physically meaningful elastic properties.
    ///
    /// `young_modulus` is measured in pascals, `density` in kg/m^3. The
    /// supported Poisson ratio range is the open interval `(-1, 0.5)`, matching
    /// stable isotropic linear elasticity.
    pub fn new(id: usize, young_modulus: f64, poisson_ratio: f64, density: f64) -> Result<Self, FemError> {
        validate_positive_material_property(id, "young_modulus", young_modulus)?;

        if !poisson_ratio.is_finite() || poisson_ratio <= -1.0 || poisson_ratio >= 0.5 {
            return Err(FemError::InvalidMaterialProperty {
                material_id: id,
                property: "poisson_ratio",
                value: poisson_ratio,
                reason: "must be finite and in the open interval (-1, 0.5)",
            });
        }

        if !density.is_finite() || density < 0.0 {
            return Err(FemError::InvalidMaterialProperty {
                material_id: id,
                property: "density",
                value: density,
                reason: "must be finite and non-negative",
            });
        }

        Ok(Self { id, young_modulus, poisson_ratio, density })
    }

    /// External material identifier.
    #[must_use]
    pub fn id(&self) -> usize {
        self.id
    }

    /// Young's modulus in pascals.
    #[must_use]
    pub fn young_modulus(&self) -> f64 {
        self.young_modulus
    }

    /// Poisson ratio.
    #[must_use]
    pub fn poisson_ratio(&self) -> f64 {
        self.poisson_ratio
    }

    /// Density in kg/m^3.
    #[must_use]
    pub fn density(&self) -> f64 {
        self.density
    }
}

fn validate_positive_material_property(material_id: usize, property: &'static str, value: f64) -> Result<(), FemError> {
    if !value.is_finite() || value <= 0.0 {
        return Err(FemError::InvalidMaterialProperty {
            material_id,
            property,
            value,
            reason: "must be finite and positive",
        });
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn creates_linear_elastic_material() {
        let material = LinearElasticMaterial::new(2, 70.0e9, 0.33, 2700.0).unwrap();

        assert_eq!(material.id(), 2);
        assert_eq!(material.young_modulus(), 70.0e9);
        assert_eq!(material.poisson_ratio(), 0.33);
        assert_eq!(material.density(), 2700.0);
    }

    #[test]
    fn rejects_non_positive_young_modulus() {
        let error = LinearElasticMaterial::new(2, 0.0, 0.3, 1.0).unwrap_err();

        assert_eq!(
            error,
            FemError::InvalidMaterialProperty {
                material_id: 2,
                property: "young_modulus",
                value: 0.0,
                reason: "must be finite and positive",
            }
        );
    }

    #[test]
    fn rejects_invalid_poisson_ratio() {
        let error = LinearElasticMaterial::new(2, 1.0, 0.5, 1.0).unwrap_err();

        assert_eq!(
            error,
            FemError::InvalidMaterialProperty {
                material_id: 2,
                property: "poisson_ratio",
                value: 0.5,
                reason: "must be finite and in the open interval (-1, 0.5)",
            }
        );
    }

    #[test]
    fn rejects_negative_density() {
        let error = LinearElasticMaterial::new(2, 1.0, 0.3, -1.0).unwrap_err();

        assert_eq!(
            error,
            FemError::InvalidMaterialProperty {
                material_id: 2,
                property: "density",
                value: -1.0,
                reason: "must be finite and non-negative",
            }
        );
    }
}
