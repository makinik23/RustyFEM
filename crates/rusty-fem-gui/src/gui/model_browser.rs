//! Discovery of local JSON models for the GUI load panel.

use std::fs;
use std::path::{Path, PathBuf};

/// Finds JSON models that are useful to expose in the GUI without a file picker.
pub(super) fn discover_json_models() -> Vec<PathBuf> {
    let mut models = Vec::new();

    collect_json_files(Path::new("examples"), &mut models);
    models.sort_by_key(|path| path.display().to_string());
    models.dedup();

    models
}

fn collect_json_files(directory: &Path, models: &mut Vec<PathBuf>) {
    let Ok(entries) = fs::read_dir(directory) else {
        return;
    };

    for entry in entries.flatten() {
        let path = entry.path();

        if path.is_file() && path.extension().is_some_and(|extension| extension == "json") {
            models.push(path);
        }
    }
}
