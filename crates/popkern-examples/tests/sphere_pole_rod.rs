//! A rod around a sphere's pole axis, holding both poles (PK-1055),
//! moved around the axis: the sphere is symmetric about it, so union,
//! difference and common keep one volume at every turn of the rod's
//! offset, while the seams the rod's loops cross move.
//!
//! Each loop crosses the sphere's seam meridian and the rod's own
//! seam; the clip cuts it at both, a long arc and a short one between
//! them. Where the two crossings sat inside the consolidation's slop
//! gate (5% of the arcs' median extent: 0.06 on a 0.6 rod), the long
//! arc's gap was welded shut at its midpoint and the short arc traced
//! twice; every boolean declined (non-manifold edge).
//!
//! A rod whose wall passes just beyond the poles meets the sphere in
//! one loop that sweeps half a turn of longitude near each pole. Three
//! samplers read that sweep from too few points: the pcurve refit (23
//! samples, its 4.4e-2 miss taken for a misfit and left unrefined), the
//! winding check (8 samples: the hole read +3.8 for its -13.4) and the
//! walk's chart embedding (16 samples: unfolded the short way through
//! the seam, and the face fell to the 3D walk).

use popkern_bool::{BoolOperation, boolean};
use popkern_core::geometry::Transform3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_measure::{body_volume, face_area};

/// Union, difference and common volumes of the unit sphere and a rod
/// of `radius` along z offset by `distance` at `turn` about the axis.
fn volumes(radius: f64, distance: f64, turn: f64) -> [f64; 3] {
    let tol = Tolerance::default();
    let ops = [
        BoolOperation::Union,
        BoolOperation::Subtraction,
        BoolOperation::Intersection,
    ];
    ops.map(|op| {
        let sphere = BRep::sphere(1.0, tol).expect("sphere");
        let rod = BRep::cylinder(radius, 3.0, tol)
            .and_then(|c| c.transform(&Transform3::translation(distance, 0.0, -1.5)))
            .and_then(|c| c.transform(&Transform3::rotation_z(turn)))
            .expect("rod");
        let at = format!("r {radius} at {distance} turn {turn:.2} {op:?}");
        let body = boolean(op, sphere, rod, tol, None).unwrap_or_else(|e| panic!("{at}: {e}"));
        body.validate()
            .unwrap_or_else(|e| panic!("{at}: invalid: {e:?}"));
        body_volume(&body, &tol).unwrap_or_else(|e| panic!("{at}: {e}"))
    })
}

/// Union and common add to the operands' volumes, difference and
/// common to the sphere's; every turn of the rod's offset keeps each.
fn invariant_over_turns(radius: f64, distance: f64, turns: &[f64]) {
    let sphere = 4.0 / 3.0 * std::f64::consts::PI;
    let rod = std::f64::consts::PI * radius * radius * 3.0;
    let reference = volumes(radius, distance, 0.0);
    let [union, minus, common] = reference;
    assert!(
        (union + common - (sphere + rod)).abs() < 1e-7,
        "r {radius} at {distance}: union {union} + common {common}"
    );
    assert!(
        (minus + common - sphere).abs() < 1e-7,
        "r {radius} at {distance}: difference {minus} + common {common}"
    );
    for &turn in turns {
        let got = volumes(radius, distance, turn);
        for (g, w) in got.iter().zip(&reference) {
            assert!(
                (g - w).abs() < 1e-7,
                "r {radius} at {distance} turn {turn:.2}: {got:?} vs {reference:?}"
            );
        }
    }
}

/// Sphere and rod-wall areas (exact boundary integrals) of each
/// boolean: `[sphere in union, sphere in difference, sphere in common,
/// wall in union, wall in difference, wall in common]`.
fn areas(radius: f64, distance: f64, turn: f64) -> [f64; 6] {
    use popkern_core::geometry::surface::SurfaceGeometry;
    let tol = Tolerance::default();
    let mut out = [0.0; 6];
    let ops = [
        BoolOperation::Union,
        BoolOperation::Subtraction,
        BoolOperation::Intersection,
    ];
    for (k, op) in ops.into_iter().enumerate() {
        let sphere = BRep::sphere(1.0, tol).expect("sphere");
        let rod = BRep::cylinder(radius, 3.0, tol)
            .and_then(|c| c.transform(&Transform3::translation(distance, 0.0, -1.5)))
            .and_then(|c| c.transform(&Transform3::rotation_z(turn)))
            .expect("rod");
        let at = format!("r {radius} at {distance} turn {turn:.2} {op:?}");
        let body = boolean(op, sphere, rod, tol, None).unwrap_or_else(|e| panic!("{at}: {e}"));
        body.validate()
            .unwrap_or_else(|e| panic!("{at}: invalid: {e:?}"));
        for (fh, f) in &body.faces {
            let slot = match f.surface.as_ref().map(|s| &s.geometry) {
                Some(SurfaceGeometry::Sphere { .. }) => k,
                Some(SurfaceGeometry::Cylinder { .. }) => 3 + k,
                _ => continue,
            };
            out[slot] += face_area(&body, fh, &tol).unwrap_or_else(|e| panic!("{at}: {e}"));
        }
    }
    out
}

#[test]
fn a_rod_whose_wall_passes_beside_the_poles_keeps_its_areas_at_every_turn() {
    // The wall 0.05, 0.02 and 0.005 beyond the poles. The sphere's
    // area splits between the difference (outside the rod) and the
    // common (inside); the wall's between the union and the common;
    // the difference keeps the common's wall. Volumes are left to the
    // OCC leg: the fine mesh of the common's sphere face overlaps
    // itself near the poles (PK-1060).
    let turns: Vec<f64> = (1..5).map(|k| 1.3 * f64::from(k)).collect();
    let sphere = 4.0 * std::f64::consts::PI;
    let wall = 2.0 * std::f64::consts::PI * 0.6 * 3.0;
    for distance in [0.55, 0.58, 0.595] {
        let reference = areas(0.6, distance, 0.0);
        let [s_union, s_minus, s_common, w_union, w_minus, w_common] = reference;
        let at = format!("at {distance}: {reference:?}");
        assert!((s_minus + s_common - sphere).abs() < 1e-7, "{at}");
        assert!((s_union - s_minus).abs() < 1e-9, "{at}");
        assert!((w_union + w_common - wall).abs() < 1e-7, "{at}");
        assert!((w_minus - w_common).abs() < 1e-9, "{at}");
        // Trims fit to the 1e-6 linear tolerance along a loop ~4 long
        // move an area by ~4e-7 (measured 3.7e-7 at turn 3.9).
        for &turn in &turns {
            let got = areas(0.6, distance, turn);
            for (g, w) in got.iter().zip(&reference) {
                assert!(
                    (g - w).abs() < 2e-6,
                    "at {distance} turn {turn:.2}: {got:?} vs {reference:?}"
                );
            }
        }
    }
}

/// The common of a rod whose wall passes 0.005 beyond the poles: its
/// sphere face's trim samples near each pole stepped half a turn of
/// longitude, the mesh rings unwrapped them the short way through the
/// seam, and the face meshed its hole: the volume read 3.34 for the
/// exact 1.401867 (PK-1059; OCC 1.401866). The bound is the residual
/// pole overlap (PK-1060), 5.4e-4 here.
#[test]
fn a_rod_beside_the_poles_measures_its_common_volume() {
    let tol = Tolerance::default();
    let sphere = BRep::sphere(1.0, tol).expect("sphere");
    let rod = BRep::cylinder(0.6, 3.0, tol)
        .and_then(|c| c.transform(&Transform3::translation(0.595, 0.0, -1.5)))
        .expect("rod");
    let common = boolean(BoolOperation::Intersection, sphere, rod, tol, None).expect("common");
    let volume = body_volume(&common, &tol).expect("volume");
    assert!((volume - 1.401867).abs() < 1e-3, "{volume}");
    // The display mesh, body-wide and face by face, covers the sphere
    // face itself (it meshed 10.9 for 3.35 on its own: PK-1059).
    let options = popkern_mesh::TessellationOptions::default();
    let body_mesh = popkern_mesh::tessellate(&common, &options).expect("mesh");
    for (fh, _) in &common.faces {
        let exact = face_area(&common, fh, &tol).expect("area");
        let own = popkern_mesh::tessellate_face(&common, fh, &options).expect("face mesh");
        let shared = body_mesh
            .face_meshes
            .iter()
            .find(|m| m.face == fh)
            .expect("face in body mesh");
        for (name, mesh) in [("face", &own), ("body", shared)] {
            let area: f64 = mesh
                .triangles
                .iter()
                .map(|t| {
                    let [a, b, c] = t.map(|i| mesh.vertices[i].position);
                    0.5 * (b - a).cross(&(c - a)).length()
                })
                .sum();
            assert!(
                (area - exact).abs() < 1e-2 * exact,
                "{fh:?} {name} mesh {area} for {exact}"
            );
        }
    }
}

#[test]
fn a_rod_holding_both_poles_keeps_its_volumes_at_every_turn() {
    // Offsets whose seam crossings fall inside the old slop gate: at
    // turn atan(1/2) the 0.112 offset is (0.1, 0.05), the PK-1055
    // fixture. Measured spread over 13 turns of each: 3e-9.
    let turns: Vec<f64> = std::iter::once(0.5_f64.atan())
        .chain((1..7).map(|k| 0.9 * f64::from(k)))
        .collect();
    for (radius, distance) in [(0.6, 0.1_f64.hypot(0.05)), (0.6, 0.2), (0.4, 0.1)] {
        invariant_over_turns(radius, distance, &turns);
    }
}
