//! Angular facet bound (`TessellationOptions::max_facet_angle`) and
//! the display preset: acceptance gates on real blend fixtures.
//!
//! - A variable-radius round (0.04 → 0.18, the playground's tapered
//!   fillet) meshed at display quality keeps every triangle's
//!   vertex-normal spread within the bound, watertight.
//! - At DEFAULT options (chord-only) the same strip must still honour
//!   its chord tolerance: before the ruled-strip transpose it declined
//!   into a boundary-only Delaunay of 12 rail-to-rail triangles, 90° of
//!   normal turn each, 7x the chord tolerance off the surface.
//! - Default meshes stay chord-only (the bound is opt-in): triangle
//!   counts of the gate meshes are unchanged by the display preset.
//!
//! Pole-normal repair is gated end-to-end by the gallery test
//! (`popkern-cli/tests/gallery.rs`): the polyhedral dice's pip poles
//! are the only fixture known to invert past v = π.

use std::collections::HashMap;

use popkern_blend::{FilletOptions, FilletRadius, FilletSpec, fillet_edges};
use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_mesh::{DISPLAY_MAX_FACET_ANGLE, TessellationOptions, TriMesh};

/// A unit cube with one top edge along x rounded by `radius`.
fn rounded_bar(radius: FilletRadius) -> BRep {
    let tol = Tolerance::default();
    let bar = BRep::unit_cube(tol).expect("cube");
    let top = bar
        .vertices
        .iter()
        .map(|(_, v)| v.position.z)
        .fold(f64::NEG_INFINITY, f64::max);
    let front = bar
        .vertices
        .iter()
        .map(|(_, v)| v.position.y)
        .fold(f64::INFINITY, f64::min);
    let edge = bar
        .edges
        .iter()
        .find(|(_, d)| {
            let he = &bar.half_edges[d.half_edge];
            let a = bar.vertices[he.origin].position;
            let b = bar.vertices[bar.half_edges[he.twin].origin].position;
            (a.z - top).abs() < 1e-9
                && (b.z - top).abs() < 1e-9
                && (a.y - front).abs() < 1e-9
                && (b.y - front).abs() < 1e-9
        })
        .map(|(e, _)| e)
        .expect("top front x edge");
    fillet_edges(
        &bar,
        &[FilletSpec { edge, radius }],
        &FilletOptions::default(),
        tol,
        None,
    )
    .expect("round serves")
    .brep
}

fn angle(a: Vec3, b: Vec3) -> f64 {
    (a.dot(&b) / (a.length() * b.length()))
        .clamp(-1.0, 1.0)
        .acos()
}

/// Worst per-triangle vertex-normal spread (radians).
fn worst_spread(mesh: &TriMesh) -> f64 {
    mesh.triangles
        .iter()
        .map(|t| {
            (0..3)
                .map(|i| {
                    angle(
                        mesh.vertices[t[i]].normal,
                        mesh.vertices[t[(i + 1) % 3]].normal,
                    )
                })
                .fold(0.0, f64::max)
        })
        .fold(0.0, f64::max)
}

/// Edges used by exactly one triangle after welding at 1e-7.
fn welded_open_edges(mesh: &TriMesh) -> usize {
    let mut key: HashMap<(i64, i64, i64), usize> = HashMap::new();
    let mut id = |p: Point3| {
        let k = (
            (p.x * 1e7).round() as i64,
            (p.y * 1e7).round() as i64,
            (p.z * 1e7).round() as i64,
        );
        let next = key.len();
        *key.entry(k).or_insert(next)
    };
    let mut count: HashMap<(usize, usize), usize> = HashMap::new();
    for t in &mesh.triangles {
        let w = [
            id(mesh.vertices[t[0]].position),
            id(mesh.vertices[t[1]].position),
            id(mesh.vertices[t[2]].position),
        ];
        if w[0] == w[1] || w[1] == w[2] || w[0] == w[2] {
            continue;
        }
        for (a, b) in [(w[0], w[1]), (w[1], w[2]), (w[2], w[0])] {
            *count.entry((a.min(b), a.max(b))).or_default() += 1;
        }
    }
    count.values().filter(|&&c| c == 1).count()
}

/// Largest distance from a triangle centroid of a curved face to that
/// face's surface (closest-point projection).
fn worst_curved_face_deviation(brep: &BRep, mesh: &TriMesh) -> f64 {
    let mut worst = 0.0_f64;
    for fm in &mesh.face_meshes {
        let surface = &brep.faces[fm.face]
            .surface
            .as_ref()
            .expect("surface")
            .geometry;
        if surface.kind_name() == "Plane" {
            continue;
        }
        for t in &fm.triangles {
            let p = [
                fm.vertices[t[0]].position,
                fm.vertices[t[1]].position,
                fm.vertices[t[2]].position,
            ];
            let c = Point3::new(
                (p[0].x + p[1].x + p[2].x) / 3.0,
                (p[0].y + p[1].y + p[2].y) / 3.0,
                (p[0].z + p[1].z + p[2].z) / 3.0,
            );
            let (u, v) = popkern_core::geometry::uv_projection::project_point_to_uv(surface, &c)
                .expect("projects");
            worst = worst.max(surface.point_at(u, v).distance_to(&c));
        }
    }
    worst
}

#[test]
fn tapered_round_meets_the_display_facet_bound_watertight() {
    let brep = rounded_bar(FilletRadius::Variable {
        start: 0.04,
        end: 0.18,
    });
    let mesh = popkern_mesh::tessellate(&brep, &TessellationOptions::display()).expect("mesh");
    assert_eq!(welded_open_edges(&mesh), 0, "watertight");
    let worst = worst_spread(&mesh);
    assert!(
        worst <= DISPLAY_MAX_FACET_ANGLE + 1e-9,
        "worst facet spread {:.2}° over the {:.1}° bound",
        worst.to_degrees(),
        DISPLAY_MAX_FACET_ANGLE.to_degrees()
    );
}

#[test]
fn tapered_round_honours_the_chord_tolerance_at_default_options() {
    let brep = rounded_bar(FilletRadius::Variable {
        start: 0.04,
        end: 0.18,
    });
    let options = TessellationOptions::default();
    let mesh = popkern_mesh::tessellate(&brep, &options).expect("mesh");
    assert_eq!(welded_open_edges(&mesh), 0, "watertight");
    let deviation = worst_curved_face_deviation(&brep, &mesh);
    assert!(
        deviation <= options.chord_height,
        "round mesh deviates {deviation:.4} from its surface (chord tolerance {})",
        options.chord_height
    );
}

#[test]
fn small_constant_round_meets_the_display_facet_bound() {
    let brep = rounded_bar(FilletRadius::Constant(0.03));
    let mesh = popkern_mesh::tessellate(&brep, &TessellationOptions::display()).expect("mesh");
    assert_eq!(welded_open_edges(&mesh), 0, "watertight");
    assert!(worst_spread(&mesh) <= DISPLAY_MAX_FACET_ANGLE + 1e-9);
    // Chord-only default: the same round facets at the chord
    // tolerance alone, far coarser in angle — the bound is opt-in.
    let coarse = popkern_mesh::tessellate(&brep, &TessellationOptions::default()).expect("mesh");
    assert!(worst_spread(&coarse) > 2.0 * DISPLAY_MAX_FACET_ANGLE);
}
