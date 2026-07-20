use serde::Serialize;

use crate::FemError;

/// Constant cross-section for a one-dimensional bar element.
#[derive(Debug, Clone, Copy, PartialEq, Serialize)]
pub struct BarSection {
    id: usize,
    area: f64,
}

impl BarSection {
    /// Creates a bar section with positive area in square meters.
    pub fn new(id: usize, area: f64) -> Result<Self, FemError> {
        if !area.is_finite() || area <= 0.0 {
            return Err(FemError::InvalidSectionProperty {
                section_id: id,
                property: "area",
                value: area,
                reason: "must be finite and positive",
            });
        }

        Ok(Self { id, area })
    }

    /// External section identifier.
    #[must_use]
    pub fn id(&self) -> usize {
        self.id
    }

    /// Cross-sectional area in square meters.
    #[must_use]
    pub fn area(&self) -> f64 {
        self.area
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn creates_bar_section() {
        let section = BarSection::new(4, 2.5e-4).unwrap();

        assert_eq!(section.id(), 4);
        assert_eq!(section.area(), 2.5e-4);
    }

    #[test]
    fn rejects_non_positive_area() {
        let error = BarSection::new(4, -1.0).unwrap_err();

        assert_eq!(
            error,
            FemError::InvalidSectionProperty {
                section_id: 4,
                property: "area",
                value: -1.0,
                reason: "must be finite and positive",
            }
        );
    }
}
