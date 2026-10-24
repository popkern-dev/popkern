//! Additional CPU render edge-case tests.

use super::*;
use crate::camera::{Camera, ViewPreset};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_mesh::{TessellationOptions, tessellate};

#[test]
fn render_empty_brep_with_wireframe_is_all_background() {
    let brep = BRep::new(Tolerance::default());
    let mesh = popkern_mesh::TriMesh::default();
    let camera = Camera::default_iso();
    let pixels = render(&brep, &mesh, &camera, 16, 16, true);
    for chunk in pixels.chunks(4) {
        assert_eq!(chunk, &[38, 38, 38, 255]);
    }
}

#[test]
fn render_cube_from_different_views() {
    let brep = BRep::unit_cube(Tolerance::default()).expect("cube builds");
    let mesh = tessellate(&brep, &TessellationOptions::default()).expect("tessellate");

    let presets = [
        ViewPreset::Front,
        ViewPreset::Back,
        ViewPreset::Top,
        ViewPreset::Left,
    ];
    for preset in &presets {
        let mut camera = Camera::default_iso();
        camera.set_view_preset(*preset);
        let pixels = render(&brep, &mesh, &camera, 32, 32, false);
        assert_eq!(pixels.len(), 32 * 32 * 4);
        // Should have some non-background pixels (model visible).
        let has_non_bg = pixels
            .chunks(4)
            .any(|px| px[0] != 38 || px[1] != 38 || px[2] != 38);
        assert!(has_non_bg, "cube should be visible from preset {preset:?}");
    }
}

#[test]
fn render_1x1_does_not_panic() {
    let brep = BRep::unit_cube(Tolerance::default()).expect("cube builds");
    let mesh = tessellate(&brep, &TessellationOptions::default()).expect("tessellate");
    let camera = Camera::default_iso();
    let pixels = render(&brep, &mesh, &camera, 1, 1, false);
    assert_eq!(pixels.len(), 4);
}

#[test]
fn render_large_resolution_is_correct_size() {
    let brep = BRep::new(Tolerance::default());
    let mesh = popkern_mesh::TriMesh::default();
    let camera = Camera::default_iso();
    let pixels = render(&brep, &mesh, &camera, 200, 100, false);
    assert_eq!(pixels.len(), 200 * 100 * 4);
}

#[test]
fn render_sphere_produces_visible_pixels() {
    let brep = BRep::sphere(1.0, Tolerance::default()).expect("sphere builds");
    let mesh = tessellate(&brep, &TessellationOptions::default()).expect("tessellate");
    let mut camera = Camera::default_iso();
    camera.fit_to_bounds(
        glam::Vec3::new(-1.0, -1.0, -1.0),
        glam::Vec3::new(1.0, 1.0, 1.0),
    );
    let pixels = render(&brep, &mesh, &camera, 64, 64, false);
    let non_bg_count = pixels
        .chunks(4)
        .filter(|px| px[0] != 38 || px[1] != 38 || px[2] != 38)
        .count();
    assert!(
        non_bg_count > 0,
        "sphere render should produce visible pixels"
    );
}

/// A closed box, wound outward (counter-clockwise seen from outside, as
/// the tessellator emits faces) or inward, with its +Z face red, its
/// -Z face blue and its sides grey.
/// A box face: its outward normal, two in-plane axes with u x v = n,
/// and its colour.
type BoxFace = ([f32; 3], [f32; 3], [f32; 3], [f32; 4]);

fn coloured_box(outward: bool) -> (Vec<mesh_data::MeshVertex>, Vec<u32>) {
    let faces: [BoxFace; 6] = [
        (
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0, 0.0, 0.0, 1.0],
        ),
        (
            [0.0, 0.0, -1.0],
            [0.0, 1.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 1.0],
        ),
        (
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [0.5, 0.5, 0.5, 1.0],
        ),
        (
            [-1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0],
            [0.0, 1.0, 0.0],
            [0.5, 0.5, 0.5, 1.0],
        ),
        (
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 0.0],
            [0.5, 0.5, 0.5, 1.0],
        ),
        (
            [0.0, -1.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0],
            [0.5, 0.5, 0.5, 1.0],
        ),
    ];
    let mut vertices = Vec::new();
    let mut indices = Vec::new();
    for (n, u, v, color) in faces {
        let base = vertices.len() as u32;
        // u x v = n, so these corners run counter-clockwise about n.
        for (a, b) in [(-1.0, -1.0), (1.0, -1.0), (1.0, 1.0), (-1.0, 1.0)] {
            let p = [
                n[0] + a * u[0] + b * v[0],
                n[1] + a * u[1] + b * v[1],
                n[2] + a * u[2] + b * v[2],
            ];
            vertices.push(mesh_data::MeshVertex {
                position: p,
                normal: n,
                color,
            });
        }
        let quad = if outward {
            [0, 1, 2, 0, 2, 3]
        } else {
            [0, 2, 1, 0, 3, 2]
        };
        indices.extend(quad.iter().map(|i| base + i));
    }
    (vertices, indices)
}

#[test]
fn the_near_face_is_drawn_whatever_the_winding() {
    // Seen from +Z (the Front preset: this camera is Y-up), a box shows
    // its +Z face. Rasterization culled clockwise screen triangles, and
    // an outward mesh lands clockwise on a y-down screen: every
    // thumbnail drew the far side of its parts (the -Z face here, a
    // bolt head's underside through its top).
    let mut camera = Camera::default_iso();
    camera.set_view_preset(ViewPreset::Front);
    camera.fit_to_bounds(glam::Vec3::splat(-1.0), glam::Vec3::splat(1.0));
    for outward in [true, false] {
        let (vertices, indices) = coloured_box(outward);
        let pixels = render_prepared(&vertices, &indices, &[], &camera, 32, 32, false);
        let centre = &pixels[(16 * 32 + 16) * 4..(16 * 32 + 16) * 4 + 4];
        assert!(
            centre[0] > 100 && centre[2] < 60,
            "outward {outward}: the +Z face (red) must face the camera, got {centre:?}"
        );
    }
}
