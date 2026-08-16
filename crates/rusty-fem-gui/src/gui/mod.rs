//! Native GUI for preprocessing and inspecting RustyFEM 2D models.
//!
//! The GUI is a separate workspace crate that consumes the public `rusty-fem`
//! API. The solver core therefore does not depend on native windowing or
//! rendering crates.
//!
//! The GUI loads and saves JSON models through the same `Model2DInput ->
//! Model2D` path as the CLI. Its direct FEM editor creates numbered nodes and
//! elements from canvas gestures, while the result view renders solved
//! displacement and stress fields.

mod app;
mod canvas;
mod document;
mod loaded_model;
mod model_browser;
mod results;
mod selection;
mod theme;
mod topology;
mod workflow;

use app::RustyFemGuiApp;
use eframe::egui;
use std::path::PathBuf;

pub(super) const DEFAULT_MODEL_PATH: &str = "examples/t3_cantilever.json";

/// Runs the native RustyFEM GUI.
pub fn run() -> eframe::Result {
    let input_path = std::env::args_os().nth(1).map(PathBuf::from).or_else(default_model_path);
    let native_options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default().with_title("RustyFEM").with_inner_size([1280.0, 820.0]),
        ..Default::default()
    };

    eframe::run_native(
        "rusty-fem-gui",
        native_options,
        Box::new(move |creation_context| {
            theme::apply(&creation_context.egui_ctx);

            Ok(Box::new(RustyFemGuiApp::new(input_path)))
        }),
    )
}

fn default_model_path() -> Option<PathBuf> {
    let path = PathBuf::from(DEFAULT_MODEL_PATH);

    path.exists().then_some(path)
}
