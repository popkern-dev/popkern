//! PK-327 global self-intersection gates (design §5 Stage 4, rows
//! 478/480).
//!
//! The grooved cylinder's outward grow past half the groove width
//! collides the groove walls: under `TrimPolicy::Trim` the unrealized
//! region (wall/floor offsets + collar undersides) is removed and the
//! two mouth-rim join tori weld along their exact coaxial
//! intersection circle; under `TrimPolicy::Decline` the operation
//! declines with reproducible region reports.

use popkern_core::geometry::Point3;
use popkern_core::tolerance::Tolerance;
use popkern_examples::offset_corpus::{grooved_cylinder, torus_lens};
use popkern_mesh::offset_oracle::{
    OffsetOracleOptions, OffsetSemantics, classify_offset_agreement,
};
use popkern_mesh::point_query::PointMeshIndex;
use popkern_mesh::{TessellationOptions, tessellate};
use popkern_shell::{ShellError, ShellOptions, ShellWarning, TrimPolicy, offset::offset_solid};
use std::f64::consts::PI;

fn tol() -> Tolerance {
    Tolerance::default()
}

fn mesh_volume(brep: &popkern_core::topology::BRep) -> f64 {
    let mesh = tessellate(brep, &TessellationOptions::default()).expect("tessellate");
    let mut vol = 0.0_f64;
    for t in &mesh.triangles {
        let (a, b, c) = (
            mesh.vertices[t[0]].position,
            mesh.vertices[t[1]].position,
            mesh.vertices[t[2]].position,
        );
        vol += a.x * (b.y * c.z - b.z * c.y) - a.y * (b.x * c.z - b.z * c.x)
            + a.z * (b.x * c.y - b.y * c.x);
    }
    (vol / 6.0).abs()
}

/// Rounded-grow closed form for the grooved cylinder at distance `d`
/// past the collision reach (2d > groove width 0.1): the plain
/// cylinder's Steiner grow minus the V-groove notch between the two
/// mouth tori (they cross at z = 0.5, r = 0.5 + sqrt(d² − 0.05²)).
fn grooved_grow_volume(d: f64) -> f64 {
    let (r, h) = (0.5_f64, 1.0_f64);
    let steiner = PI * r * r * h
        + (2.0 * PI * r * h + 2.0 * PI * r * r) * d
        + (PI * h + PI * PI * r) * d * d
        + (4.0 * PI / 3.0) * d * d * d;
    // Notch: 2π ∫₀^{w/2} [(r+d)² − (r + √(d²−s²))²]/… ds — expanded:
    // (r+d)² − (r+g)² with g = √(d²−s²) = 2r(d−g) + d² − g².
    let w2 = 0.05_f64;
    let n = 4000;
    let mut integral = 0.0;
    for i in 0..n {
        let s = w2 * (i as f64 + 0.5) / n as f64;
        let g = (d * d - s * s).sqrt();
        integral += (2.0 * r * (d - g) + d * d - g * g) * (w2 / n as f64);
    }
    // Notch volume of revolution: ΔV = 2π ∫ ((r+d)² − r_boundary²)/2
    // over the groove band; band symmetry folds the two halves into
    // twice the half-band integral, cancelling the /2.
    steiner - 2.0 * PI * integral
}

#[test]
fn grooved_cylinder_grow_past_reach_trims_and_serves() {
    let fixture = grooved_cylinder(tol()).unwrap();
    let d = 0.08;
    let result = offset_solid(&fixture.brep, d, &ShellOptions::default(), tol())
        .expect("outward grow past the groove reach must trim and serve");

    // The trim happened and reported itself.
    let trimmed: Vec<_> = result
        .warnings
        .iter()
        .filter_map(|w| match w {
            ShellWarning::RegionTrimmed {
                removed_faces,
                welded_seams,
                regions,
            } => Some((*removed_faces, *welded_seams, regions.len())),
            _ => None,
        })
        .collect();
    assert_eq!(
        trimmed.len(),
        1,
        "exactly one trim report: {:?}",
        result.warnings
    );
    let (removed, welded, regions) = trimmed[0];
    assert!(
        removed >= 3,
        "wall+floor offsets and collar undersides die: {removed}"
    );
    assert!(
        welded >= 1,
        "the mouth tori weld along their crossing: {welded}"
    );
    assert!(regions >= 1, "reports carried");

    // Validator-clean (the full set, self-intersection included).
    let errors = result.brep.validate_all();
    assert!(errors.is_empty(), "served invalid: {errors:?}");

    // The groove is gone: bands (2) + caps (2) + cap-rim torus halves
    // (4) + mouth collar halves (4).
    assert_eq!(result.brep.faces.len(), 12, "trimmed face census");

    // Exact rounded-grow volume.
    let expected = grooved_grow_volume(d);
    let measured = mesh_volume(&result.brep);
    let rel = (measured - expected).abs() / expected;
    assert!(
        rel < 0.01,
        "volume {measured:.6} vs closed form {expected:.6} (rel {rel:.5})"
    );

    // SDF oracle: the served body IS the Minkowski grow.
    let source_mesh = tessellate(&fixture.brep, &TessellationOptions::default()).unwrap();
    let served_mesh = tessellate(&result.brep, &TessellationOptions::default()).unwrap();
    let source_index = PointMeshIndex::build(&source_mesh);
    let served_index = PointMeshIndex::build(&served_mesh);
    let report = classify_offset_agreement(
        &source_index,
        &served_index,
        OffsetSemantics::Grow { distance: d },
        &OffsetOracleOptions::default(),
        |_| false,
    );
    assert!(report.is_clean(), "SDF oracle disagreements: {report:?}");
    assert!(report.classified_fraction() > 0.5);
}

#[test]
fn grooved_cylinder_below_reach_stays_untouched_by_the_trim_stage() {
    let fixture = grooved_cylinder(tol()).unwrap();
    let result = offset_solid(&fixture.brep, 0.03, &ShellOptions::default(), tol())
        .expect("sub-reach grow serves");
    assert!(
        !result
            .warnings
            .iter()
            .any(|w| matches!(w, ShellWarning::RegionTrimmed { .. })),
        "no trim below the reach: {:?}",
        result.warnings
    );
    assert!(result.brep.validate_all().is_empty());
}

#[test]
fn decline_policy_reports_regions_instead_of_trimming() {
    let fixture = grooved_cylinder(tol()).unwrap();
    let options = ShellOptions::default().with_trim_policy(TrimPolicy::Decline);
    let err = offset_solid(&fixture.brep, 0.08, &options, tol())
        .expect_err("Decline policy must refuse the colliding grow");
    match err {
        ShellError::SelfIntersecting {
            regions,
            policy_note,
        } => {
            assert!(!regions.is_empty(), "reports must name the collisions");
            assert!(policy_note.contains("Decline"));
            // Reports are reproducible: witnesses sit on the collision
            // circle band around z = 0.5.
            for r in &regions {
                let w: Point3 = r.witness;
                assert!(
                    (0.4..=0.6).contains(&w.z),
                    "witness off the groove band: {w:?}"
                );
                assert!(r.curve_length > 0.0);
            }
        }
        other => panic!("expected SelfIntersecting, got {other:?}"),
    }
}

/// Apple volume: revolution of the offset tube disk (radius r'
/// centered at distance R from the axis) clipped to the ρ ≥ 0
/// half-plane: V = 2π ∫∫ ρ dA over the clipped disk.
fn apple_volume(major: f64, minor_grown: f64) -> f64 {
    let (r, rp) = (major, minor_grown);
    let w = (rp * rp - r * r).sqrt();
    2.0 * PI
        * (PI * r * rp * rp + (2.0 / 3.0) * w * w * w + r * r * w
            - r * rp * rp * (PI / 2.0 - (r / rp).asin()))
}

#[test]
fn torus_lens_past_reach_serves_the_apple_with_pinch_poles() {
    let fixture = torus_lens(tol()).unwrap();
    let d = 0.2;
    let result = offset_solid(&fixture.brep, d, &ShellOptions::default(), tol())
        .expect("spindled torus grow must serve the apple");

    let trimmed: Vec<_> = result
        .warnings
        .iter()
        .filter_map(|w| match w {
            ShellWarning::RegionTrimmed {
                removed_faces,
                welded_seams,
                regions,
            } => Some((*removed_faces, *welded_seams, regions.len())),
            _ => None,
        })
        .collect();
    assert_eq!(
        trimmed.len(),
        1,
        "one spindle trim report: {:?}",
        result.warnings
    );
    assert_eq!(trimmed[0].0, 8, "the two dead-window rows die");
    assert_eq!(trimmed[0].2, 1, "one spindle region report");

    let errors = result.brep.validate_all();
    assert!(errors.is_empty(), "served invalid: {errors:?}");

    // 8 surviving band faces + 8 pinch fan sectors, all on the
    // spindle torus; two exact pinch poles on the axis at ±w.
    assert_eq!(result.brep.faces.len(), 16, "apple face census");
    let w = ((0.45_f64 + d).powi(2) - 0.36).sqrt();
    assert!((w - 0.25).abs() < 1e-12);
    let mut pole_hits = 0;
    for (_, v) in result.brep.vertices.iter() {
        let p = v.position;
        if (p.x * p.x + p.y * p.y).sqrt() < 1e-9 && ((p.z.abs() - w).abs() < 1e-9) {
            pole_hits += 1;
        }
    }
    assert_eq!(pole_hits, 2, "two exact pinch poles");

    // Exact apple volume.
    let expected = apple_volume(0.6, 0.45 + d);
    let measured = mesh_volume(&result.brep);
    let rel = (measured - expected).abs() / expected;
    assert!(
        rel < 0.01,
        "volume {measured:.6} vs apple closed form {expected:.6} (rel {rel:.5})"
    );

    // SDF oracle: the apple IS the Minkowski grow of the torus.
    let source_mesh = tessellate(&fixture.brep, &TessellationOptions::default()).unwrap();
    let served_mesh = tessellate(&result.brep, &TessellationOptions::default()).unwrap();
    let source_index = PointMeshIndex::build(&source_mesh);
    let served_index = PointMeshIndex::build(&served_mesh);
    let report = classify_offset_agreement(
        &source_index,
        &served_index,
        OffsetSemantics::Grow { distance: d },
        &OffsetOracleOptions::default(),
        |_| false,
    );
    assert!(report.is_clean(), "SDF oracle disagreements: {report:?}");
    assert!(report.classified_fraction() > 0.5);
}

#[test]
fn torus_lens_decline_policy_reports_the_pinch() {
    let fixture = torus_lens(tol()).unwrap();
    let options = ShellOptions::default().with_trim_policy(TrimPolicy::Decline);
    let err = offset_solid(&fixture.brep, 0.2, &options, tol())
        .expect_err("Decline policy must refuse the spindled grow");
    match err {
        ShellError::SelfIntersecting {
            regions,
            policy_note,
        } => {
            assert_eq!(regions.len(), 1);
            assert!(policy_note.contains("Decline"));
            // Witness = a pinch point on the axis.
            let p = regions[0].witness;
            assert!((p.x * p.x + p.y * p.y).sqrt() < 1e-9 && (p.z.abs() - 0.25).abs() < 1e-9);
        }
        other => panic!("expected SelfIntersecting, got {other:?}"),
    }
}

#[test]
fn distance_map_counter_is_zero_on_analytic_fixtures_and_runs_on_wavy() {
    // Row 480: analytic-only inputs never enter the numeric tier —
    // the FreeFormFoldScan warning (the tier-3 counter) must be absent
    // on all-analytic fixtures and PRESENT on the free-form one.
    use popkern_examples::offset_corpus::{l_bracket, wavy_plate};
    use popkern_shell::thicken::thicken;

    let lb = l_bracket(tol()).unwrap();
    let served = thicken(&lb.brep, 0.02, &ShellOptions::default(), tol()).unwrap();
    assert!(
        !served
            .warnings
            .iter()
            .any(|w| matches!(w, ShellWarning::FreeFormFoldScan { .. })),
        "analytic fixture entered the numeric tier: {:?}",
        served.warnings
    );

    let wavy = wavy_plate(tol()).unwrap();
    let served = thicken(&wavy.brep, 0.02, &ShellOptions::default(), tol()).unwrap();
    let scans: Vec<_> = served
        .warnings
        .iter()
        .filter_map(|w| match w {
            ShellWarning::FreeFormFoldScan { faces } => Some(*faces),
            _ => None,
        })
        .collect();
    assert_eq!(
        scans,
        vec![2],
        "both NURBS skins scanned clean: {:?}",
        served.warnings
    );
    assert!(served.brep.validate_all().is_empty());
}
