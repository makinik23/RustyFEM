//! Defines the degrees of freedom (DOF) for used types of nodes.

/// 2D degrees of freedom (DOF) for nodes.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Dof2D {
    /// Displacement in the x-direction.
    Ux,

    /// Displacement in the y-direction.
    Uy,

    /// Rotation about the z-axis.
    Rz,
}

impl Dof2D {
    pub fn name(self) -> &'static str {
        match self {
            Self::Ux => "Ux",
            Self::Uy => "Uy",
            Self::Rz => "Rz",
        }
    }
}
