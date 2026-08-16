//! GUI workflow modes and direct FEM drawing tool state.

/// Top-level way of using the GUI.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum WorkMode {
    Inspect,
    DrawFem,
}

impl WorkMode {
    pub(super) fn label(self) -> &'static str {
        match self {
            Self::Inspect => "Inspect",
            Self::DrawFem => "Draw FEM",
        }
    }
}

/// Active tool in the direct FEM drawing workflow.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum DrawTool {
    Select,
    InsertNode,
    InsertElement,
    MoveNode,
    Constraint,
    NodalLoad,
    EdgeTraction,
}

impl DrawTool {
    pub(super) const ALL: [Self; 7] = [
        Self::Select,
        Self::InsertNode,
        Self::InsertElement,
        Self::MoveNode,
        Self::Constraint,
        Self::NodalLoad,
        Self::EdgeTraction,
    ];

    pub(super) fn label(self) -> &'static str {
        match self {
            Self::Select => "Select",
            Self::InsertNode => "Insert node",
            Self::InsertElement => "Insert element",
            Self::MoveNode => "Move node",
            Self::Constraint => "Constraint",
            Self::NodalLoad => "Nodal load",
            Self::EdgeTraction => "Edge traction",
        }
    }
}

/// Element type selected for direct element insertion.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum ElementDraftKind {
    Truss,
    Beam,
    T3,
    T6,
    Q4,
    Q8,
}

impl ElementDraftKind {
    pub(super) const ALL: [Self; 6] = [Self::Truss, Self::Beam, Self::T3, Self::T6, Self::Q4, Self::Q8];

    pub(super) fn label(self) -> &'static str {
        match self {
            Self::Truss => "Truss",
            Self::Beam => "Beam",
            Self::T3 => "T3",
            Self::T6 => "T6",
            Self::Q4 => "Q4",
            Self::Q8 => "Q8",
        }
    }

    /// Number of geometric points the user places on the canvas.
    ///
    /// Quadratic elements create their midside nodes automatically, so T6
    /// and Q8 require only the same corner points as T3 and Q4.
    pub(super) fn placement_point_count(self) -> usize {
        match self {
            Self::Truss | Self::Beam => 2,
            Self::T3 | Self::T6 => 3,
            Self::Q4 | Self::Q8 => 4,
        }
    }

    pub(super) fn is_surface(self) -> bool {
        matches!(self, Self::T3 | Self::T6 | Self::Q4 | Self::Q8)
    }
}
