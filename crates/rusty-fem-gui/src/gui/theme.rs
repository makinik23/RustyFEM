//! Minimal visual theme for the native GUI.

use eframe::egui::{self, Color32, CornerRadius, Margin, Stroke, Vec2};

pub(super) const APP_BACKGROUND: Color32 = Color32::from_rgb(245, 247, 250);
pub(super) const PANEL_BACKGROUND: Color32 = Color32::from_rgb(250, 252, 255);
pub(super) const CANVAS_BACKGROUND: Color32 = Color32::from_rgb(248, 250, 252);
pub(super) const BORDER: Color32 = Color32::from_rgb(203, 213, 225);
pub(super) const TEXT: Color32 = Color32::from_rgb(15, 23, 42);
pub(super) const MUTED_TEXT: Color32 = Color32::from_rgb(100, 116, 139);
pub(super) const ACCENT: Color32 = Color32::from_rgb(14, 116, 144);
pub(super) const WARNING: Color32 = Color32::from_rgb(180, 83, 9);
pub(super) const ERROR: Color32 = Color32::from_rgb(185, 28, 28);

pub(super) fn apply(context: &egui::Context) {
    context.set_theme(egui::Theme::Light);
    context.style_mut_of(egui::Theme::Light, |style| {
        style.spacing.item_spacing = Vec2::new(8.0, 6.0);
        style.spacing.button_padding = Vec2::new(10.0, 5.0);
        style.spacing.menu_margin = Margin::symmetric(8, 6);
        style.spacing.window_margin = Margin::same(10);

        style.visuals = egui::Visuals::light();
        style.visuals.panel_fill = APP_BACKGROUND;
        style.visuals.window_fill = PANEL_BACKGROUND;
        style.visuals.window_stroke = Stroke::new(1.0, BORDER);
        style.visuals.extreme_bg_color = Color32::WHITE;
        style.visuals.faint_bg_color = Color32::from_rgb(241, 245, 249);
        style.visuals.selection.bg_fill = ACCENT;
        style.visuals.selection.stroke = Stroke::new(1.0, Color32::WHITE);

        for widget in [
            &mut style.visuals.widgets.noninteractive,
            &mut style.visuals.widgets.inactive,
            &mut style.visuals.widgets.hovered,
            &mut style.visuals.widgets.active,
            &mut style.visuals.widgets.open,
        ] {
            widget.corner_radius = CornerRadius::same(4);
        }

        style.visuals.widgets.inactive.bg_stroke = Stroke::new(1.0, BORDER);
        style.visuals.widgets.hovered.bg_fill = Color32::from_rgb(236, 245, 248);
        style.visuals.widgets.hovered.bg_stroke = Stroke::new(1.0, Color32::from_rgb(125, 165, 180));
        style.visuals.widgets.active.bg_fill = Color32::from_rgb(219, 237, 243);
        style.visuals.widgets.active.bg_stroke = Stroke::new(1.0, ACCENT);
    });
}
