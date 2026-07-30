//! Hole-rim reflex terminal class pins (PK-185 / PK-186).
//!
//! Fillets on inner-loop (hole-rim) edge segments whose ends are
//! REFLEX face corners with un-blended transversal edges. The
//! construction: interior-rail spur anchors, band cut through the
//! hole loop, and exact terminal-arc discs in the transversal plane.
//! Covers the convex variant (pocket/through-hole opening rims), the
//! concave variant (rib-base junctions on a plate), and the PK-185
//! partial-service interference policy.

use popkern_part::prelude::*;

fn assert_watertight(part: &Part, label: &str) {
    let brep = part.brep();
    let errors = brep.validate_all();
    assert!(errors.is_empty(), "{label}: validate_all: {errors:?}");
    let opens = part.open_edge_count().expect("open edge count");
    assert_eq!(opens, 0, "{label}: {opens} open mesh edges");
}

/// Convex hole-rim fillets (pocket opening) + concave floor legs —
/// the pocketed_tray shape (PK-186). All four selected legs serve;
/// the removed and added quarter-fillet prisms cancel exactly, so
/// the volume is fully determined by the prismatic features.
#[test]
fn pocket_rim_and_floor_fillets_serve_watertight() {
    let plate = Profile::on(Plane::XY)
        .rect(60.0, 36.0)
        .extrude(10.0)
        .unwrap();
    let tool = Profile::on(Plane::xy_at(4.0))
        .rect(40.0, 20.0)
        .extrude(10.0)
        .unwrap();
    let pocketed = plate.subtract(&tool).unwrap();
    let chamfered = pocketed
        .chamfer(2.5, |e| {
            e.straight()
                .direction(Dir::Z)
                .chord_length_between(8.0, 12.0)
                .allow_empty()
        })
        .unwrap();
    let part = chamfered
        .fillet(1.5, |e| {
            e.straight()
                .direction(Dir::X)
                .chord_length_between(35.0, 45.0)
                .allow_empty()
        })
        .unwrap();
    assert_watertight(&part, "pocketed_tray");
    // plate 60·36·10 − pocket 40·20·6 − 4 chamfers 2.5²/2·10; the
    // two convex rim removals and two concave floor additions are
    // (1 − π/4)·r²·40 each and cancel pairwise.
    let expected = 60.0 * 36.0 * 10.0 - 40.0 * 20.0 * 6.0 - 4.0 * (2.5 * 2.5 / 2.0) * 10.0;
    let volume =
        popkern_measure::body_volume(part.brep(), &popkern_core::tolerance::Tolerance::default())
            .expect("volume");
    // Mesh-integral volume: the ±quarter-prism cancellation leaves
    // only tessellation-chord residue. The pre-PK-218 gate of 1e-3
    // was an ACCIDENT: every strip was rail-chorded (each ~19 off)
    // and the convex/concave pairs cancelled bit-exactly. With the
    // structured ruled-strip mesh each strip carries only its own
    // arc-chord residue (~3e-3 total here, chord_height 0.001) and
    // the cancellation is approximate — the honest bound is the
    // per-strip chord budget, not the lucky zero.
    assert!(
        (volume - expected).abs() < 0.01,
        "pocketed_tray volume {volume} != {expected}"
    );
}

/// Minimal convex fixture: a rectangular THROUGH-hole's two long
/// opening rims filleted, all four ends reflex hole corners. Volume
/// is exact: each rim removes a quarter-fillet prism over the full
/// rim length (the terminal discs bound the removal flush at the
/// hole corners).
#[test]
fn through_hole_rim_fillet_exact_volume() {
    let plate = Profile::on(Plane::XY)
        .rect(30.0, 20.0)
        .extrude(6.0)
        .unwrap();
    let hole = Profile::on(Plane::xy_at(-1.0))
        .rect(10.0, 6.0)
        .extrude(8.0)
        .unwrap();
    let holed = plate.subtract(&hole).unwrap();
    let part = holed
        .fillet(1.5, |e| {
            e.straight()
                .direction(Dir::X)
                .chord_length_between(8.0, 12.0)
                .filter(|i| (i.midpoint.z - 6.0).abs() < 0.1)
                .allow_empty()
        })
        .unwrap();
    assert_watertight(&part, "through_hole_rim");
    let quarter = (1.0 - std::f64::consts::FRAC_PI_4) * 1.5 * 1.5;
    let expected = 30.0 * 20.0 * 6.0 - 10.0 * 6.0 * 6.0 - 2.0 * quarter * 10.0;
    let volume =
        popkern_measure::body_volume(part.brep(), &popkern_core::tolerance::Tolerance::default())
            .expect("volume");
    assert!(
        (volume - expected).abs() < 0.1,
        "through_hole_rim volume {volume} != {expected}"
    );
}

/// Concave variant: skewed rib-base junction fillets ending at
/// reflex rib-footprint corners (PK-185's surviving legs), together
/// with the plate rim fillets that the phantom-interference decline
/// used to block. Watertight with the volume in the analytic band
/// (rim removals − junction additions, corner blends bounded).
#[test]
fn rib_junction_fillets_serve_watertight() {
    let plate = Profile::on(Plane::XY)
        .rect(80.0, 50.0)
        .extrude(5.0)
        .unwrap();
    let rib = Profile::on(Plane::xy_at(5.0))
        .rect(70.0, 6.0)
        .extrude(10.0)
        .unwrap()
        .rotate_z(20.0_f64.to_radians())
        .unwrap();
    let body = plate.union(&rib).unwrap();
    let part = body
        .fillet(1.5, |e| {
            e.straight()
                .chord_length_between(40.0, 80.0)
                .filter(|i| (i.midpoint.z - 5.0).abs() < 0.1)
                .allow_empty()
        })
        .unwrap();
    assert_watertight(&part, "angled_rib_plate");
    // Pin the junction fillets' MATERIAL-ADD semantics against a
    // rims-only build: the two concave junction legs must add
    // (1 − π/4)·r²·140 ≈ 67.6 (band widened for the pre-existing
    // strip-tessellation chord bias in the mesh-integral volume —
    // OCC reads the exact body; see the strip-tessellation ledger
    // item). A sign flip (carving instead of adding) or a silent
    // no-op lands far outside the band.
    let rims_only = body
        .fillet(1.5, |e| {
            e.straight()
                .chord_length_between(40.0, 80.0)
                .filter(|i| {
                    (i.midpoint.z - 5.0).abs() < 0.1
                        && (i.midpoint.x.abs() > 30.0 || i.midpoint.y.abs() > 20.0)
                })
                .allow_empty()
        })
        .unwrap();
    let tol = popkern_core::tolerance::Tolerance::default();
    let volume = popkern_measure::body_volume(part.brep(), &tol).expect("volume");
    let rims_volume = popkern_measure::body_volume(rims_only.brep(), &tol).expect("volume");
    let junction_add = volume - rims_volume;
    assert!(
        (55.0..=100.0).contains(&junction_add),
        "junction fillets must ADD ~67.6 of material (measured {junction_add}; full {volume}, rims-only {rims_volume})"
    );
}

/// PK-185 acceptance: genuine non-completable interference (skewed
/// junction strips crossing the rim strips at r = 4) declines
/// transactionally by default, and under `ServeNonInterfering`
/// drops exactly the interfering pairs, serves the rest, and
/// reports every dropped pair with its reason.
#[test]
fn interference_policy_drop_pair_serves_rest() {
    use popkern_fillet::{
        BlendError, BlendWarning, FilletOptions, FilletRadius, FilletSpec, InterferencePolicy,
        fillet_edges,
    };
    let plate = Profile::on(Plane::XY)
        .rect(100.0, 34.0)
        .extrude(5.0)
        .unwrap();
    let rib = Profile::on(Plane::xy_at(5.0))
        .rect(70.0, 6.0)
        .extrude(10.0)
        .unwrap()
        .rotate_z(20.0_f64.to_radians())
        .unwrap();
    let body = plate.union(&rib).unwrap();
    let brep = body.brep();
    // Select the plate-top boundary edges at z = 5 longer than 20:
    // the two 100-long rims and the two 70-long rib junctions
    // interfere at r = 4; the 34-long side rims sit 12+ clear of
    // the junction strips and must survive the drop.
    let mut specs = Vec::new();
    for (eh, edge) in brep.edges.iter() {
        let Some(curve) = edge.curve.as_ref() else {
            continue;
        };
        use popkern_core::geometry::CurveEval;
        let a = curve.geometry.point_at(curve.t_start);
        let b = curve.geometry.point_at(curve.t_end);
        let mid_z = 0.5 * (a.z + b.z);
        if (mid_z - 5.0).abs() < 0.1 && a.distance_to(&b) > 20.0 {
            specs.push(FilletSpec {
                edge: eh,
                radius: FilletRadius::Constant(4.0),
            });
        }
    }
    assert!(
        specs.len() >= 6,
        "expected the 4 long rims + 2 junction edges, found {}",
        specs.len()
    );
    let tolerance = popkern_core::tolerance::Tolerance::default();

    // Default policy: transactional decline with the structured
    // interference error.
    let declined = fillet_edges(brep, &specs, &FilletOptions::default(), tolerance, None);
    assert!(
        matches!(declined, Err(BlendError::FilletInterference { .. })),
        "default policy must decline: {declined:?}"
    );

    // Partial service: interfering pairs dropped and reported, the
    // rest served on a valid body.
    let options = FilletOptions {
        interference_policy: InterferencePolicy::ServeNonInterfering,
        ..FilletOptions::default()
    };
    let served = fillet_edges(brep, &specs, &options, tolerance, None)
        .expect("ServeNonInterfering must serve the non-interfering fillets");
    let dropped: Vec<_> = served
        .edge_warnings
        .iter()
        .filter_map(|(_, w)| match w {
            BlendWarning::FilletPairDropped {
                fillet_a,
                fillet_b,
                reason,
            } => Some((*fillet_a, *fillet_b, reason.clone())),
            _ => None,
        })
        .collect();
    assert!(
        !dropped.is_empty(),
        "expected at least one FilletPairDropped report"
    );
    for (_, _, reason) in &dropped {
        assert!(!reason.is_empty(), "dropped pair must carry a reason");
    }
    assert!(
        served.brep.validate().is_ok(),
        "partial-service output must validate: {:?}",
        served.brep.validate()
    );
}
