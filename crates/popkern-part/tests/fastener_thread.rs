//! PK-460 thread-engine serves: external and internal threads
//! against watertight and analytic grid-CSG volume oracles.
//!
//! The volume model is the ISO basic form evaluated per (θ, z): the
//! material boundary is star-shaped about the axis, so
//! `V = ∫∫ R(θ,z)²/2 dθ dz` integrates it exactly up to grid
//! resolution. Test lengths are exact pitch multiples so the volume
//! is independent of the construction's thread phase (every θ column
//! integrates whole periods). Tolerance 1.5%: covers tessellation
//! sag and the rational helix's ≈0.5%-of-lead axial ripple, while
//! catching any structural defect (a phantom ring, a doubled groove,
//! a glued protrusion) by an order of magnitude.

use popkern_core::tolerance::Tolerance;
use popkern_part::Part;
use popkern_part::fastener::{Handedness, HexNut, ThreadSpec, ThreadedRod};

fn mesh_volume(part: &Part) -> f64 {
    let mesh = part
        .mesh(&popkern_mesh::TessellationOptions::default())
        .expect("tessellate");
    let mut v6 = 0.0;
    for t in &mesh.triangles {
        let a = mesh.vertices[t[0]].position;
        let b = mesh.vertices[t[1]].position;
        let c = mesh.vertices[t[2]].position;
        v6 += a.x * (b.y * c.z - b.z * c.y) - a.y * (b.x * c.z - b.z * c.x)
            + a.z * (b.x * c.y - b.y * c.x);
    }
    (v6 / 6.0).abs()
}

/// External thread surface radius at axial distance `s` from the
/// nearest ridge crest centreline (ISO basic form).
fn external_radius(spec: &ThreadSpec, s: f64) -> f64 {
    let p = spec.pitch();
    let tan_a = (spec.flank_angle() / 2.0).tan();
    let r_major = spec.major_diameter() / 2.0;
    let r_minor = spec.minor_diameter() / 2.0;
    if s <= p / 16.0 {
        r_major
    } else if s < 3.0 * p / 8.0 {
        (r_major - (s - p / 16.0) / tan_a).max(r_minor)
    } else {
        r_minor
    }
}

/// Distance from `z` to the nearest thread crest centreline at
/// azimuth `theta` (crests repeat every pitch across starts).
fn crest_distance(spec: &ThreadSpec, theta: f64, z: f64) -> f64 {
    let p = spec.pitch();
    let dir = match spec.handedness() {
        Handedness::Right => 1.0,
        Handedness::Left => -1.0,
    };
    let u = (z - dir * theta / std::f64::consts::TAU * spec.lead()).rem_euclid(p);
    u.min(p - u)
}

/// Grid-CSG volume of a full-thread rod of `length` (no tip
/// chamfer), phase-independent because `length` is a pitch multiple.
fn expected_rod_volume(spec: &ThreadSpec, length: f64) -> f64 {
    let n_theta = 720;
    let n_z = 1500;
    let d_theta = std::f64::consts::TAU / n_theta as f64;
    let dz = length / n_z as f64;
    let mut v = 0.0;
    for i in 0..n_theta {
        let theta = (i as f64 + 0.5) * d_theta;
        for j in 0..n_z {
            let z = (j as f64 + 0.5) * dz;
            let r = external_radius(spec, crest_distance(spec, theta, z));
            v += 0.5 * r * r * d_theta * dz;
        }
    }
    v
}

/// Internal thread hole radius at axial distance `s` from the
/// nearest ridge centreline, grown by `clearance` (matches
/// `thread_form::internal_ridge`'s in-material profile).
fn internal_hole_radius(spec: &ThreadSpec, clearance: f64, s: f64) -> f64 {
    let p = spec.pitch();
    let half = spec.flank_angle() / 2.0;
    let (sin_a, cos_a) = half.sin_cos();
    let tan_a = half.tan();
    let r_major = spec.major_diameter() / 2.0;
    let r_minor = spec.minor_diameter() / 2.0;
    let r_bore = r_minor + clearance;
    let flank_r0 = r_minor + clearance * sin_a;
    let flank_w0 = 3.0 * p / 8.0 + clearance * cos_a;
    let r_tip = r_major + clearance;
    let w_tip = flank_w0 - tan_a * (r_tip - flank_r0);
    if s <= w_tip {
        r_tip
    } else {
        // Invert the flank: w(r) = flank_w0 − tanα (r − flank_r0).
        (flank_r0 + (flank_w0 - s) / tan_a).max(r_bore)
    }
}

fn expected_nut_volume(spec: &ThreadSpec, clearance: f64, across_flats: f64, m: f64) -> f64 {
    let hex = 3.0_f64.sqrt() / 2.0 * across_flats * across_flats * m;
    let n_theta = 720;
    let n_z = 1500;
    let d_theta = std::f64::consts::TAU / n_theta as f64;
    let dz = m / n_z as f64;
    let mut hole = 0.0;
    for i in 0..n_theta {
        let theta = (i as f64 + 0.5) * d_theta;
        for j in 0..n_z {
            let z = (j as f64 + 0.5) * dz;
            let r = internal_hole_radius(spec, clearance, crest_distance(spec, theta, z));
            hole += 0.5 * r * r * d_theta * dz;
        }
    }
    hex - hole
}

fn serve_and_check_volume(part: &Part, expected: f64, label: &str) {
    let open = part.open_edge_count().expect("mesh");
    assert_eq!(open, 0, "{label}: {open} open mesh edges");
    let volume = mesh_volume(part);
    let rel = (volume - expected).abs() / expected;
    assert!(
        rel < 0.015,
        "{label}: volume {volume:.6} vs expected {expected:.6} ({:.2}% off)",
        rel * 100.0
    );
}

#[test]
fn m6_rod_serves_with_iso_volume() {
    let spec = ThreadSpec::parse("M6").expect("M6");
    let rod = ThreadedRod::new(spec.clone(), 4.0)
        .tip_chamfer(false)
        .build(Tolerance::default())
        .expect("rod");
    serve_and_check_volume(&rod, expected_rod_volume(&spec, 4.0), "M6 rod");
}

#[test]
fn left_hand_rod_mirrors_the_right_hand_volume() {
    let spec = ThreadSpec::parse("M6").expect("M6").left_handed();
    let rod = ThreadedRod::new(spec.clone(), 4.0)
        .tip_chamfer(false)
        .build(Tolerance::default())
        .expect("LH rod");
    serve_and_check_volume(&rod, expected_rod_volume(&spec, 4.0), "M6 LH rod");
}

#[test]
fn two_start_fine_rod_serves() {
    let spec = ThreadSpec::parse("M8x1").expect("M8 fine").with_starts(2);
    let rod = ThreadedRod::new(spec.clone(), 5.0)
        .tip_chamfer(false)
        .build(Tolerance::default())
        .expect("two-start rod");
    serve_and_check_volume(&rod, expected_rod_volume(&spec, 5.0), "M8x1 two-start rod");
}

#[test]
fn rod_with_tip_chamfer_serves_watertight() {
    let spec = ThreadSpec::parse("M6").expect("M6");
    let rod = ThreadedRod::new(spec, 4.0)
        .build(Tolerance::default())
        .expect("chamfered rod");
    assert_eq!(rod.open_edge_count().expect("mesh"), 0);
}

#[test]
fn m6_nut_thread_serves_with_grown_volume() {
    let spec = ThreadSpec::parse("M6").expect("M6");
    let clearance = 0.15;
    // Bare thread geometry for the oracle: chamfers off, thickness a
    // pitch multiple, across-flats from the standard.
    let nut = HexNut::new(spec.clone())
        .clearance(clearance)
        .thickness(5.0)
        .corner_chamfer(false)
        .build(Tolerance::default())
        .expect("nut");
    let s = popkern_part::fastener::standards::hex_head(6.0).across_flats;
    serve_and_check_volume(
        &nut,
        expected_nut_volume(&spec, clearance, s, 5.0),
        "M6 nut",
    );
}

/// Off-pitch rod lengths (PK-474): a cap phase that lands mid-flank
/// puts the blank-seam helix records on the tool's multi-turn NURBS
/// flank where a cold multisheet projection reads a phantom ~pitch
/// deviation — the surface-honesty gate dropped an on-surface seam
/// record, the crest chart's seam ladder lost one link, and the
/// graze filter cascaded 8 half-claimed strands. The piece's own
/// pcurve is the chart-true witness.
#[test]
fn off_pitch_rod_lengths_serve_watertight() {
    for len in [5.25, 5.1, 5.4, 6.75] {
        let spec = ThreadSpec::parse("M6").expect("M6");
        let rod = ThreadedRod::new(spec, len)
            .tip_chamfer(false)
            .build(Tolerance::default())
            .unwrap_or_else(|e| panic!("rod {len}: {e:?}"));
        assert_eq!(rod.open_edge_count().expect("mesh"), 0, "rod {len}");
    }
}
