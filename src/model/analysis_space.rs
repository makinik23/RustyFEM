//! Types of FEM analysis spaces.

use std::str::FromStr;

use crate::FemError;

/// Describes whether the FEM model uses two or three spatial dimensions.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AnalysisSpace {
    TwoDimensional,
    ThreeDimensional,
}

impl AnalysisSpace {
    /// Returns the number of spatial dimensions for a given analysis space.
    pub fn spatial_dimension(self) -> usize {
        match self {
            Self::TwoDimensional => 2,
            Self::ThreeDimensional => 3,
        }
    }
}

impl FromStr for AnalysisSpace {
    type Err = FemError;

    fn from_str(value: &str) -> Result<Self, Self::Err> {
        let normalized = value.trim().to_ascii_lowercase();

        match normalized.as_str() {
            "2" | "2d" => Ok(Self::TwoDimensional),
            "3" | "3d" => Ok(Self::ThreeDimensional),
            _ => Err(FemError::InvalidAnalysisSpace { value: value.to_owned() }),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::AnalysisSpace;
    use crate::FemError;

    #[test]
    fn two_dimensional_space_has_two_dimensions() {
        assert_eq!(AnalysisSpace::TwoDimensional.spatial_dimension(), 2);
    }

    #[test]
    fn three_dimensional_space_has_three_dimensions() {
        assert_eq!(AnalysisSpace::ThreeDimensional.spatial_dimension(), 3);
    }

    #[test]
    fn parses_two_dimensional_space() {
        assert_eq!("2d".parse::<AnalysisSpace>().expect("valid space"), AnalysisSpace::TwoDimensional);
    }

    #[test]
    fn parses_three_dimensional_space() {
        assert_eq!("3d".parse::<AnalysisSpace>().expect("valid space"), AnalysisSpace::ThreeDimensional);
    }

    #[test]
    fn rejects_unknown_analysis_space() {
        assert!(matches!("4d".parse::<AnalysisSpace>(), Err(FemError::InvalidAnalysisSpace { .. })));
    }
}
