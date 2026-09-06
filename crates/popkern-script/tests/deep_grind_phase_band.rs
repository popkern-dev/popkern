//! PK-452: off-station deep coil grinds across the phase band.
//!
//! A cut plane through the 6.5-turn coil used to produce PHASE-LUCK
//! outcomes away from the pinned stations (the station itself is
//! `station_aligned_deep_grind_serves_census_clean`):
//!
//! - z = 0.85/0.875/0.89 served an EMPTY body as success: the chart
//!   embedding's projection run anchored on the ADJACENT COIL TURN
//!   (cold 5-start projection, the PK-446 multi-sheet class), the
//!   seam-crossing cut oval's unfolded span collapsed below one
//!   period, both its endpoints clustered into ONE seam branch, the
//!   walk isolated the oval instead of splitting the wall into
//!   bands, and the implicit-closure shoelace filed it as a HOLE —
//!   one interior sample above the plane then classified the whole
//!   wall Inside and postprocess spurred the remnants to nothing.
//! - z = 0.8 declined with an Euler-Poincaré violation from the same
//!   substrate.
//!
//! Fixed by sheet-correct multi-sheet anchors in the chart embedding
//! and the non-contractible (chart-winding) orbit role rule: a
//! wrapping orbit can never take the hole role. This pin holds the
//! band to full serves: watertight mesh within 2% of the analytic
//! deep-cut volume, primary cellular, validation clean.

use popkern_script::Plane;
use popkern_script::body_sweep::HelixParams;
use popkern_script::prelude::*;

fn mesh_volume_and_open_edges(brep: &popkern_core::topology::BRep) -> (f64, usize) {
    let options = popkern_mesh::TessellationOptions {
        chord_height: 5.0e-4,
        min_subdivisions: 4,
        max_subdivisions: 512,
    };
    let mesh = popkern_mesh::tessellate(brep, &options).expect("mesh");
    let mut vol6 = 0.0;
    for t in &mesh.triangles {
        let a = mesh.vertices[t[0]].position;
        let b = mesh.vertices[t[1]].position;
        let c = mesh.vertices[t[2]].position;
        vol6 += a.x * (b.y * c.z - b.z * c.y) - a.y * (b.x * c.z - b.z * c.x)
            + a.z * (b.x * c.y - b.y * c.x);
    }
    let mut edge_uses: std::collections::HashMap<(u64, u64), i64> =
        std::collections::HashMap::new();
    let key = |p: popkern_core::geometry::Point3| -> u64 {
        let q = |x: f64| ((x * 1.0e7).round() as i64) as u64;
        q(p.x)
            .wrapping_mul(0x9E3779B97F4A7C15)
            .wrapping_add(q(p.y).wrapping_mul(0xC2B2AE3D27D4EB4F))
            .wrapping_add(q(p.z))
    };
    for t in &mesh.triangles {
        for (i, j) in [(0usize, 1usize), (1, 2), (2, 0)] {
            let a = key(mesh.vertices[t[i]].position);
            let b = key(mesh.vertices[t[j]].position);
            *edge_uses.entry((a.min(b), a.max(b))).or_insert(0) += if a < b { 1 } else { -1 };
        }
    }
    let open = edge_uses.values().filter(|&&c| c != 0).count();
    (vol6 / 6.0, open)
}

/// Analytic volume of the 6.5-turn coil with everything above
/// `plane_z` removed AND everything above `plane_z + 0.2` restored
/// when the coil pokes through the plate's top (plate thickness 0.2).
/// Same section-tilt caveat as the loft fidelity suite — well inside
/// the 2% gate.
fn coil_band_cut_volume(plane_z: f64) -> f64 {
    let big_r = 0.3_f64;
    let pitch = 0.15_f64;
    let wire_r = 0.04_f64;
    let turns = 6.5_f64;
    let speed = ((std::f64::consts::TAU * big_r).powi(2) + pitch * pitch).sqrt();
    let full = std::f64::consts::PI * wire_r * wire_r * turns * speed;
    let seg = |d: f64| -> f64 {
        wire_r * wire_r * (d / wire_r).acos() - d * (wire_r * wire_r - d * d).sqrt()
    };
    // Cross-section area removed at station s for a plate occupying
    // z in [plane_z, plane_z + 0.2].
    let removed_at = |s: f64| -> f64 {
        let zc = pitch * s;
        let circle = std::f64::consts::PI * wire_r * wire_r;
        let d_bot = plane_z - zc; // signed distance of plate bottom above centre
        let d_top = plane_z + 0.2 - zc;
        let above = |d: f64| -> f64 {
            if d <= -wire_r {
                circle
            } else if d >= wire_r {
                0.0
            } else {
                seg(d)
            }
        };
        above(d_bot) - above(d_top)
    };
    let n = 8192;
    let mut removed = 0.0;
    for k in 0..n {
        let s = turns * (k as f64 + 0.5) / n as f64;
        removed += removed_at(s) * speed * turns / n as f64;
    }
    full - removed
}

fn grind_serves_watertight(z: f64) {
    let mut m = Model::new();
    let z_axis = AxisSpec::Literal {
        origin: [0.0, 0.0, 0.0],
        direction: [0.0, 0.0, 1.0],
    };
    let spring = m
        .profile(Plane::XY)
        .circle(0.04)
        .sketch("wire_section")
        .expect("wire sketch")
        .helical_sweep(&mut m, z_axis, HelixParams::new(0.3, 0.15, 6.5))
        .expect("coil");
    let top = m
        .profile(Plane::xy_at(z))
        .rect(1.0, 1.0)
        .sketch("grind_top")
        .expect("plate sketch")
        .extrude(&mut m, End::Blind(Scalar::from(0.2)))
        .expect("plate");
    let _ = popkern_bool::merge_fallback::take_merge_fallback_events();
    let ground = spring
        .subtract(&mut m, top)
        .unwrap_or_else(|e| panic!("deep grind at z={z} must serve: {e:?}"));
    let brep = m.brep_of(ground).expect("brep");
    assert!(
        brep.validate().is_ok(),
        "z={z}: invalid: {:?}",
        brep.validate().err()
    );
    assert!(
        brep.validate_curve_surface().is_ok(),
        "z={z}: curve-surface lie: {:?}",
        brep.validate_curve_surface().err()
    );
    let fallbacks = popkern_bool::merge_fallback::take_merge_fallback_events();
    assert!(
        fallbacks.is_empty(),
        "z={z}: fell back to graph-merge: {fallbacks:?}"
    );
    assert!(!brep.faces.is_empty(), "z={z}: EMPTY body served");
    let (vol, open) = mesh_volume_and_open_edges(&brep);
    assert_eq!(open, 0, "z={z}: mesh has {open} open edges");
    let expect = coil_band_cut_volume(z);
    let rel = (vol - expect).abs() / expect;
    assert!(
        rel < 0.02,
        "z={z}: mesh volume {vol:.6} vs analytic {expect:.6} ({:.2}% off)",
        rel * 100.0
    );
}

/// The empty-serve phase: the cut oval's embedding anchored on the
/// wrong coil turn and the whole wall classified away.
#[test]
fn off_station_deep_grind_serves_at_the_former_empty_phase() {
    grind_serves_watertight(0.85);
}

/// The Euler-violation phase: same substrate, declined instead of
/// serving empty. The plate top also pierces the coil tip here, so
/// the serve carries the extra tip piece.
#[test]
fn off_station_deep_grind_serves_at_the_former_euler_phase() {
    grind_serves_watertight(0.8);
}
