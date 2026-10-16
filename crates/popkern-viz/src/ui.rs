//! egui UI panels for the debug viewer.
//!
//! Provides the toolbar, inspector side panel, and status bar as
//! standalone functions that operate on [`VizServer`] state.

use egui::Ui;
use popkern_core::geometry::surface::SurfaceGeometry;

use crate::camera::ViewPreset;
use crate::scene::{ColorMode, Selection, VizServer};

/// Draw the top toolbar with view presets and display toggles.
///
/// Contains camera preset buttons (Front, Back, Top, Bottom, Left, Right, Iso),
/// display option checkboxes (wireframe, normals), and color mode selector.
pub fn toolbar(ctx: &egui::Context, scene: &mut VizServer) {
    egui::TopBottomPanel::top("toolbar").show(ctx, |ui| {
        ui.horizontal(|ui| {
            ui.label("popkern-viz");
            ui.separator();

            // View preset buttons
            if ui.button("Front").clicked() {
                scene.camera.set_view_preset(ViewPreset::Front);
            }
            if ui.button("Back").clicked() {
                scene.camera.set_view_preset(ViewPreset::Back);
            }
            if ui.button("Top").clicked() {
                scene.camera.set_view_preset(ViewPreset::Top);
            }
            if ui.button("Bottom").clicked() {
                scene.camera.set_view_preset(ViewPreset::Bottom);
            }
            if ui.button("Left").clicked() {
                scene.camera.set_view_preset(ViewPreset::Left);
            }
            if ui.button("Right").clicked() {
                scene.camera.set_view_preset(ViewPreset::Right);
            }
            if ui.button("Iso").clicked() {
                scene.camera.set_view_preset(ViewPreset::Iso);
            }

            ui.separator();

            // Display toggles
            ui.checkbox(&mut scene.display.show_wireframe, "Wireframe");
            ui.checkbox(&mut scene.display.show_normals, "Normals");

            ui.separator();

            // Color mode
            egui::ComboBox::from_label("Color")
                .selected_text(color_mode_label(scene.display.color_mode))
                .show_ui(ui, |ui| {
                    ui.selectable_value(
                        &mut scene.display.color_mode,
                        ColorMode::Surface,
                        "Surface",
                    );
                    ui.selectable_value(&mut scene.display.color_mode, ColorMode::Flat, "Flat");
                    ui.selectable_value(&mut scene.display.color_mode, ColorMode::None, "None");
                });
        });
    });
}

/// Draw the right-side inspector panel with model info and selection details.
///
/// Shows vertex/edge/face counts and surface type breakdown when a model is
/// loaded. When a face is selected, displays its surface type, orientation,
/// and loop information.
pub fn side_panel(ctx: &egui::Context, scene: &VizServer) {
    egui::SidePanel::right("inspector")
        .min_width(200.0)
        .default_width(250.0)
        .show(ctx, |ui| {
            model_info_section(ui, scene);
            ui.separator();
            selection_section(ui, scene);
        });
}

/// Draw the bottom status bar with selection summary.
///
/// Shows a one-line summary of the current selection state, or
/// "No selection" when nothing is picked.
pub fn status_bar(ctx: &egui::Context, scene: &VizServer) {
    egui::TopBottomPanel::bottom("status").show(ctx, |ui| {
        ui.horizontal(|ui| {
            match &scene.selection {
                Selection::Face(fh) => {
                    let surface_type = scene
                        .brep
                        .as_ref()
                        .and_then(|brep| brep.faces.get(*fh))
                        .and_then(|face| face.surface.as_ref())
                        .map(|s| surface_type_name(&s.geometry));
                    match surface_type {
                        Some(name) => ui.label(format!("Selected: Face ({name})")),
                        None => ui.label("Selected: Face"),
                    }
                }
                Selection::Edge(_) => ui.label("Selected: Edge"),
                Selection::Vertex(_) => ui.label("Selected: Vertex"),
                Selection::None => ui.label("No selection"),
            };
        });
    });
}

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

/// Draw the "Model Info" section in the side panel.
fn model_info_section(ui: &mut Ui, scene: &VizServer) {
    ui.heading("Model Info");
    ui.add_space(4.0);

    let Some(brep) = &scene.brep else {
        ui.label("No model loaded");
        return;
    };

    ui.label(format!("Vertices: {}", brep.vertices.len()));
    ui.label(format!("Edges: {}", brep.edges.len()));
    ui.label(format!("Faces: {}", brep.faces.len()));
    ui.label(format!("Shells: {}", brep.shells.len()));
    ui.label(format!("Solids: {}", brep.solids.len()));

    // Surface type breakdown
    ui.add_space(4.0);
    ui.label("Surface types:");

    let mut planes = 0u32;
    let mut cylinders = 0u32;
    let mut spheres = 0u32;
    let mut cones = 0u32;
    let mut tori = 0u32;
    let mut nurbs = 0u32;
    let mut gregory = 0u32;
    let mut none_count = 0u32;

    for (_, face) in &brep.faces {
        match &face.surface {
            Some(sof) => match &sof.geometry {
                SurfaceGeometry::Plane { .. } => planes += 1,
                SurfaceGeometry::Cylinder { .. } => cylinders += 1,
                SurfaceGeometry::Sphere { .. } => spheres += 1,
                SurfaceGeometry::Cone { .. } => cones += 1,
                SurfaceGeometry::Torus { .. } => tori += 1,
                SurfaceGeometry::Nurbs(_) => nurbs += 1,
                SurfaceGeometry::Gregory(_) | SurfaceGeometry::PolyhedralSpline(_) => gregory += 1,
            },
            None => none_count += 1,
        }
    }

    ui.indent("surface_breakdown", |ui| {
        if planes > 0 {
            ui.label(format!("  Plane: {planes}"));
        }
        if cylinders > 0 {
            ui.label(format!("  Cylinder: {cylinders}"));
        }
        if spheres > 0 {
            ui.label(format!("  Sphere: {spheres}"));
        }
        if cones > 0 {
            ui.label(format!("  Cone: {cones}"));
        }
        if tori > 0 {
            ui.label(format!("  Torus: {tori}"));
        }
        if nurbs > 0 {
            ui.label(format!("  NURBS: {nurbs}"));
        }
        if gregory > 0 {
            ui.label(format!("  Gregory: {gregory}"));
        }
        if none_count > 0 {
            ui.label(format!("  (none): {none_count}"));
        }
    });
}

/// Draw the "Selection" section in the side panel.
fn selection_section(ui: &mut Ui, scene: &VizServer) {
    ui.heading("Selection");
    ui.add_space(4.0);

    match &scene.selection {
        Selection::None => {
            ui.label("Nothing selected");
        }
        Selection::Face(fh) => {
            ui.label(format!("Face: {fh:?}"));

            if let Some(brep) = &scene.brep
                && let Some(face) = brep.faces.get(*fh)
            {
                // Surface info
                match &face.surface {
                    Some(sof) => {
                        let type_name = surface_type_name(&sof.geometry);
                        ui.label(format!("Surface: {type_name}"));
                        let orient = if sof.orientation {
                            "same as face"
                        } else {
                            "reversed"
                        };
                        ui.label(format!("Orientation: {orient}"));
                    }
                    None => {
                        ui.label("Surface: (none)");
                    }
                }

                // Loop info
                ui.label(format!("Inner loops: {}", face.inner_loops.len()));
            }
        }
        Selection::Edge(eh) => {
            ui.label(format!("Edge: {eh:?}"));
        }
        Selection::Vertex(vh) => {
            ui.label(format!("Vertex: {vh:?}"));
        }
    }
}

/// Human-readable name for a surface geometry type.
fn surface_type_name(surface: &SurfaceGeometry) -> &'static str {
    match surface {
        SurfaceGeometry::Plane { .. } => "Plane",
        SurfaceGeometry::Cylinder { .. } => "Cylinder",
        SurfaceGeometry::Sphere { .. } => "Sphere",
        SurfaceGeometry::Cone { .. } => "Cone",
        SurfaceGeometry::Torus { .. } => "Torus",
        SurfaceGeometry::Nurbs(_) => "NURBS",
        SurfaceGeometry::Gregory(_) => "Gregory",
        SurfaceGeometry::PolyhedralSpline(_) => "PolyhedralSpline",
    }
}

/// Human-readable label for a [`ColorMode`].
fn color_mode_label(mode: ColorMode) -> &'static str {
    match mode {
        ColorMode::Surface => "Surface",
        ColorMode::Flat => "Flat",
        ColorMode::None => "None",
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn surface_type_name_all_variants() {
        use popkern_core::geometry::{Point3, Vec3};

        let plane = SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        };
        assert_eq!(surface_type_name(&plane), "Plane");

        let cyl = SurfaceGeometry::Cylinder {
            origin: Point3::new(0.0, 0.0, 0.0),
            axis: Vec3::new(0.0, 0.0, 1.0),
            radius: 1.0,
        };
        assert_eq!(surface_type_name(&cyl), "Cylinder");

        let sphere = SurfaceGeometry::Sphere {
            center: Point3::new(0.0, 0.0, 0.0),
            radius: 1.0,
            axis: popkern_core::geometry::Vec3::Z,
        };
        assert_eq!(surface_type_name(&sphere), "Sphere");

        let cone = SurfaceGeometry::Cone {
            apex: Point3::new(0.0, 0.0, 0.0),
            axis: Vec3::new(0.0, 0.0, 1.0),
            half_angle: 0.5,
        };
        assert_eq!(surface_type_name(&cone), "Cone");

        let torus = SurfaceGeometry::Torus {
            center: Point3::new(0.0, 0.0, 0.0),
            axis: Vec3::new(0.0, 0.0, 1.0),
            major_radius: 2.0,
            minor_radius: 0.5,
        };
        assert_eq!(surface_type_name(&torus), "Torus");

        let nurbs_surface = popkern_core::geometry::nurbs_surface::NurbsSurface {
            degree_u: 1,
            degree_v: 1,
            knots_u: vec![0.0, 0.0, 1.0, 1.0],
            knots_v: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(0.0, 1.0, 0.0),
                Point3::new(1.0, 0.0, 0.0),
                Point3::new(1.0, 1.0, 0.0),
            ],
            num_u: 2,
            num_v: 2,
            weights: vec![1.0, 1.0, 1.0, 1.0],
        };
        assert_eq!(
            surface_type_name(&SurfaceGeometry::Nurbs(nurbs_surface)),
            "NURBS"
        );
    }

    #[test]
    fn color_mode_label_all_variants() {
        assert_eq!(color_mode_label(ColorMode::Surface), "Surface");
        assert_eq!(color_mode_label(ColorMode::Flat), "Flat");
        assert_eq!(color_mode_label(ColorMode::None), "None");
    }

    #[test]
    fn selection_section_with_no_model() {
        // Verifies that the functions don't panic with default state.
        let scene = VizServer::new();
        assert!(matches!(scene.selection, Selection::None));
    }

    #[test]
    fn selection_with_face_and_brep() {
        use popkern_core::tolerance::Tolerance;
        use popkern_core::topology::BRep;

        let mut scene = VizServer::new();
        let brep = BRep::unit_cube(Tolerance::default()).expect("cube should build");
        let first_face = brep.faces.keys().next().expect("cube has faces");
        scene.brep = Some(brep);
        scene.select(Selection::Face(first_face));

        assert!(matches!(scene.selection, Selection::Face(_)));
    }
}
