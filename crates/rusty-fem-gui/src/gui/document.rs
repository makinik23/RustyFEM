//! Editable document history and JSON persistence for the GUI.

use rusty_fem::io::Model2DInput;
use rusty_fem::model::Model2D;
use std::fs;
use std::path::Path;

const HISTORY_LIMIT: usize = 100;

/// Full-model snapshots used by Undo and Redo.
#[derive(Default)]
pub(super) struct EditHistory {
    undo: Vec<Model2DInput>,
    redo: Vec<Model2DInput>,
}

impl EditHistory {
    pub(super) fn clear(&mut self) {
        self.undo.clear();
        self.redo.clear();
    }

    pub(super) fn record_snapshot(&mut self, snapshot: Model2DInput) {
        self.undo.push(snapshot);

        if self.undo.len() > HISTORY_LIMIT {
            self.undo.remove(0);
        }

        self.redo.clear();
    }

    pub(super) fn undo(&mut self, current: &Model2D) -> Option<Model2DInput> {
        let previous = self.undo.pop()?;
        self.redo.push(Model2DInput::from_model(current));
        Some(previous)
    }

    pub(super) fn redo(&mut self, current: &Model2D) -> Option<Model2DInput> {
        let next = self.redo.pop()?;
        self.undo.push(Model2DInput::from_model(current));
        Some(next)
    }

    pub(super) fn can_undo(&self) -> bool {
        !self.undo.is_empty()
    }

    pub(super) fn can_redo(&self) -> bool {
        !self.redo.is_empty()
    }
}

pub(super) fn write_model_json(path: &Path, model: &Model2D) -> Result<(), String> {
    if let Some(parent) = path.parent().filter(|parent| !parent.as_os_str().is_empty()) {
        fs::create_dir_all(parent).map_err(|error| format!("Could not create {}: {error}", parent.display()))?;
    }

    let input = Model2DInput::from_model(model);
    let contents =
        serde_json::to_string_pretty(&input).map_err(|error| format!("Could not serialize model: {error}"))?;

    fs::write(path, format!("{contents}\n")).map_err(|error| format!("Could not write {}: {error}", path.display()))
}

#[cfg(test)]
mod tests {
    use super::EditHistory;
    use rusty_fem::io::Model2DInput;
    use rusty_fem::model::{Model2D, Node2D};

    #[test]
    fn undo_and_redo_exchange_complete_model_snapshots() {
        let mut model = Model2D::new();
        let mut history = EditHistory::default();
        history.record_snapshot(Model2DInput::from_model(&model));
        model.add_node(Node2D::new(1, 2.0, 3.0).expect("valid node")).expect("node should be added");

        let previous = history.undo(&model).expect("undo snapshot should exist");
        assert!(previous.nodes.is_empty());

        let restored = previous.into_model().expect("previous snapshot should build");
        let next = history.redo(&restored).expect("redo snapshot should exist");
        assert_eq!(next.nodes.len(), 1);
        assert_eq!(next.nodes[0].id, 1);
    }
}
