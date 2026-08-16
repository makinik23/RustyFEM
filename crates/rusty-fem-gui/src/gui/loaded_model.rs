//! Loading and lightweight metadata for GUI-inspected models.
//!
//! The GUI does not define its own model format. JSON files are deserialized
//! into `Model2DInput` and then converted into the canonical `Model2D`, so all
//! existing validation rules remain in force.

use rusty_fem::io::Model2DInput;
use rusty_fem::model::Model2D;
use std::fs;
use std::path::{Path, PathBuf};

/// A model loaded into the GUI, together with file and view metadata.
pub(super) struct LoadedModel {
    pub(super) path: PathBuf,
    pub(super) model: Model2D,
    pub(super) bounds: Option<ModelBounds>,
}

impl LoadedModel {
    pub(super) fn from_json_file(path: &Path) -> Result<Self, String> {
        let contents =
            fs::read_to_string(path).map_err(|error| format!("Could not read {}: {error}", path.display()))?;
        let input: Model2DInput =
            serde_json::from_str(&contents).map_err(|error| format!("Could not parse {}: {error}", path.display()))?;
        let model = input.into_model().map_err(|error| format!("Could not build model: {error}"))?;
        let bounds = ModelBounds::from_model(&model);

        Ok(Self { path: path.to_path_buf(), model, bounds })
    }

    pub(super) fn from_model(path: PathBuf, model: Model2D) -> Self {
        let bounds = ModelBounds::from_model(&model);

        Self { path, model, bounds }
    }

    pub(super) fn refresh_bounds(&mut self) {
        self.bounds = ModelBounds::from_model(&self.model);
    }
}

/// Axis-aligned model extents used to fit model coordinates into the canvas.
#[derive(Debug, Clone, Copy)]
pub(super) struct ModelBounds {
    pub(super) min_x: f64,
    pub(super) max_x: f64,
    pub(super) min_y: f64,
    pub(super) max_y: f64,
}

impl ModelBounds {
    pub(super) fn drawing_default() -> Self {
        Self { min_x: -10.0, max_x: 10.0, min_y: -10.0, max_y: 10.0 }
    }

    fn from_model(model: &Model2D) -> Option<Self> {
        let first = model.nodes().first()?;
        let mut bounds = Self { min_x: first.x(), max_x: first.x(), min_y: first.y(), max_y: first.y() };

        for node in model.nodes().iter().skip(1) {
            bounds.min_x = bounds.min_x.min(node.x());
            bounds.max_x = bounds.max_x.max(node.x());
            bounds.min_y = bounds.min_y.min(node.y());
            bounds.max_y = bounds.max_y.max(node.y());
        }

        Some(bounds)
    }
}
