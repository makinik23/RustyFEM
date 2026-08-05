//! Section data for a two-dimensional finite element model.

use crate::FemError;

/// Cross-section properties for a 2D truss element.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TrussSection2D {
    cross_section_area: f64,
}

/// Cross-section properties for a 2D Euler-Bernoulli beam element.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BeamSection2D {
    cross_section_area: f64,
    second_moment_of_area: f64,
    section_height: Option<f64>,
}

/// Thickness properties for a 2D plane-stress element.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PlaneStressSection2D {
    thickness: f64,
}

/// Section definition used by supported 2D elements.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Section2D {
    Truss(TrussSection2D),
    Beam(BeamSection2D),
    PlaneStress(PlaneStressSection2D),
}

/// Stores section definitions used by a 2D model.
#[derive(Default)]
pub struct Sections2D {
    sections: Vec<(usize, Section2D)>,
}

impl TrussSection2D {
    /// Creates a truss section after validating its area.
    pub fn new(cross_section_area: f64) -> Result<Self, FemError> {
        validate_positive_section_property("truss", "cross-sectional area", cross_section_area)?;

        Ok(Self { cross_section_area })
    }

    /// Returns the cross-sectional area.
    #[must_use]
    pub fn cross_section_area(&self) -> f64 {
        self.cross_section_area
    }
}

impl BeamSection2D {
    /// Creates a beam section without a stored section height.
    pub fn new(cross_section_area: f64, second_moment_of_area: f64) -> Result<Self, FemError> {
        Self::new_internal(cross_section_area, second_moment_of_area, None)
    }

    /// Creates a beam section with a height for fiber-stress recovery.
    pub fn new_with_section_height(
        cross_section_area: f64, second_moment_of_area: f64, section_height: f64,
    ) -> Result<Self, FemError> {
        Self::new_internal(cross_section_area, second_moment_of_area, Some(section_height))
    }

    fn new_internal(
        cross_section_area: f64, second_moment_of_area: f64, section_height: Option<f64>,
    ) -> Result<Self, FemError> {
        validate_positive_section_property("beam", "cross-sectional area", cross_section_area)?;
        validate_positive_section_property("beam", "second moment of area", second_moment_of_area)?;

        if let Some(value) = section_height {
            validate_positive_section_property("beam", "section height", value)?;
        }

        Ok(Self { cross_section_area, second_moment_of_area, section_height })
    }

    /// Returns the cross-sectional area.
    #[must_use]
    pub fn cross_section_area(&self) -> f64 {
        self.cross_section_area
    }

    /// Returns the second moment of area.
    #[must_use]
    pub fn second_moment_of_area(&self) -> f64 {
        self.second_moment_of_area
    }

    /// Returns the section height if it was provided.
    #[must_use]
    pub fn section_height(&self) -> Option<f64> {
        self.section_height
    }
}

impl PlaneStressSection2D {
    /// Creates a plane-stress section after validating its thickness.
    pub fn new(thickness: f64) -> Result<Self, FemError> {
        validate_positive_section_property("plane_stress", "thickness", thickness)?;

        Ok(Self { thickness })
    }

    /// Returns the section thickness.
    #[must_use]
    pub fn thickness(&self) -> f64 {
        self.thickness
    }
}

impl Section2D {
    /// Returns a stable section type name for diagnostics.
    #[must_use]
    pub fn section_type(&self) -> &'static str {
        match self {
            Self::Truss(_) => "truss",
            Self::Beam(_) => "beam",
            Self::PlaneStress(_) => "plane_stress",
        }
    }
}

impl Sections2D {
    /// Creates an empty section collection.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Adds a section with a model-unique ID.
    pub fn add_section(&mut self, section_id: usize, section: Section2D) -> Result<(), FemError> {
        if self.contains_section_id(section_id) {
            return Err(FemError::DuplicateId { entity: "section", id: section_id });
        }

        self.sections.push((section_id, section));

        Ok(())
    }

    /// Returns a section by ID.
    pub fn section(&self, section_id: usize) -> Result<&Section2D, FemError> {
        self.sections
            .iter()
            .find(|(id, _)| *id == section_id)
            .map(|(_, section)| section)
            .ok_or(FemError::UnknownId { entity: "section", id: section_id })
    }

    /// Returns a truss section by ID.
    pub fn truss_section(&self, section_id: usize) -> Result<&TrussSection2D, FemError> {
        match self.section(section_id)? {
            Section2D::Truss(section) => Ok(section),
            other => Err(FemError::InvalidSectionType { section_id, expected: "truss", actual: other.section_type() }),
        }
    }

    /// Returns a beam section by ID.
    pub fn beam_section(&self, section_id: usize) -> Result<&BeamSection2D, FemError> {
        match self.section(section_id)? {
            Section2D::Beam(section) => Ok(section),
            other => Err(FemError::InvalidSectionType { section_id, expected: "beam", actual: other.section_type() }),
        }
    }

    /// Returns a plane-stress section by ID.
    pub fn plane_stress_section(&self, section_id: usize) -> Result<&PlaneStressSection2D, FemError> {
        match self.section(section_id)? {
            Section2D::PlaneStress(section) => Ok(section),
            other => {
                Err(FemError::InvalidSectionType { section_id, expected: "plane_stress", actual: other.section_type() })
            }
        }
    }

    /// Returns all sections in insertion order.
    #[must_use]
    pub fn sections(&self) -> &[(usize, Section2D)] {
        &self.sections
    }

    /// Checks whether a section ID already exists.
    #[must_use]
    pub(crate) fn contains_section_id(&self, section_id: usize) -> bool {
        self.sections.iter().any(|(id, _)| *id == section_id)
    }
}

fn validate_positive_section_property(
    section_type: &'static str, property: &'static str, value: f64,
) -> Result<(), FemError> {
    if !value.is_finite() || value <= 0.0 {
        return Err(FemError::InvalidSectionProperty {
            section_type,
            property,
            value,
            reason: "must be finite and strictly positive",
        });
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::{BeamSection2D, PlaneStressSection2D, Section2D, Sections2D, TrussSection2D};
    use crate::FemError;

    #[test]
    fn creates_empty_sections() {
        let sections = Sections2D::new();

        assert!(sections.sections().is_empty());
    }

    #[test]
    fn stores_sections_in_insertion_order() {
        let mut sections = Sections2D::new();
        let truss = Section2D::Truss(TrussSection2D::new(0.02).expect("valid section"));
        let beam = Section2D::Beam(BeamSection2D::new(0.03, 0.001).expect("valid section"));

        sections.add_section(10, truss).expect("section should be added");
        sections.add_section(20, beam).expect("section should be added");

        assert_eq!(sections.sections(), &[(10, truss), (20, beam)]);
        assert_eq!(sections.section(10).expect("section should exist"), &truss);
        assert_eq!(sections.section(20).expect("section should exist"), &beam);
    }

    #[test]
    fn rejects_duplicate_section_ids() {
        let mut sections = Sections2D::new();
        let section = Section2D::Truss(TrussSection2D::new(0.02).expect("valid section"));

        sections.add_section(10, section).expect("section should be added");
        let result = sections.add_section(10, section);

        assert!(matches!(result, Err(FemError::DuplicateId { entity: "section", id: 10 })));
    }

    #[test]
    fn rejects_unknown_section_ids() {
        let sections = Sections2D::new();

        assert!(matches!(sections.section(99), Err(FemError::UnknownId { entity: "section", id: 99 })));
    }

    #[test]
    fn returns_typed_sections() {
        let mut sections = Sections2D::new();
        let truss = TrussSection2D::new(0.02).expect("valid truss section");
        let beam = BeamSection2D::new_with_section_height(0.03, 0.001, 0.2).expect("valid beam section");
        let plane_stress = PlaneStressSection2D::new(0.1).expect("valid plane-stress section");

        sections.add_section(10, Section2D::Truss(truss)).expect("section should be added");
        sections.add_section(20, Section2D::Beam(beam)).expect("section should be added");
        sections.add_section(30, Section2D::PlaneStress(plane_stress)).expect("section should be added");

        assert_eq!(sections.truss_section(10).expect("section should exist"), &truss);
        assert_eq!(sections.beam_section(20).expect("section should exist"), &beam);
        assert_eq!(sections.plane_stress_section(30).expect("section should exist"), &plane_stress);
    }

    #[test]
    fn rejects_invalid_typed_section_request() {
        let mut sections = Sections2D::new();
        let section = Section2D::Truss(TrussSection2D::new(0.02).expect("valid section"));

        sections.add_section(10, section).expect("section should be added");
        let result = sections.beam_section(10);

        assert!(matches!(
            result,
            Err(FemError::InvalidSectionType { section_id: 10, expected: "beam", actual: "truss" })
        ));
    }

    #[test]
    fn rejects_invalid_section_properties() {
        let cases = [("zero", 0.0), ("negative", -1.0), ("infinite", f64::INFINITY), ("not a number", f64::NAN)];

        for (name, value) in cases {
            let result = TrussSection2D::new(value);

            assert!(
                matches!(
                    result,
                    Err(FemError::InvalidSectionProperty {
                        section_type: "truss",
                        property: "cross-sectional area",
                        value: actual_value,
                        reason: "must be finite and strictly positive",
                    }) if actual_value == value || (actual_value.is_nan() && value.is_nan())
                ),
                "failed case: {name}"
            );
        }
    }
}
