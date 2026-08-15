//! PK-287 — closed-rim pcurve staleness measurement.
//!
//! A boolean hole leaves two closed circular rim edges. Two authoring/
//! measurement defects made their pcurves read permanently stale:
//!
//! 1. The trace-distance fallback was a nearest-of-128-samples scan,
//!    floored at half the sample spacing (~5.4e-3 on a r=0.3 rim) —
//!    too coarse to certify a geometrically coincident pcurve whose
//!    parameterization speed differs from the edge curve's (rational
//!    NURBS arc spans vs angle-uniform). Fixed by refined true
//!    distance-to-trace with a designed direction contract.
//! 2. The plane-face pcurve for a NURBS rim curve came from a
//!    23-sample projection walk + cubic interpolation (~1.4e-4
//!    sagitta, above the stale gate), and `fix_stale_pcurves` could
//!    never measure its own recompute better — 4 stale forever. Fixed
//!    by the exact affine NURBS-in-plane pcurve candidate.

use popkern_bool::BoolOperation;
use popkern_core::geometry::transform::Transform3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

fn holed_cube() -> BRep {
    let tol = Tolerance::default();
    let cube = BRep::unit_cube(tol).expect("cube");
    let cyl = BRep::cylinder(0.3, 2.0, tol)
        .expect("cylinder")
        .transform(&Transform3::translation(0.5, 0.5, -0.5))
        .expect("translate");
    popkern_bool::boolean(BoolOperation::Subtraction, cube, cyl, tol, None).expect("subtraction")
}

/// The boolean output must ship rim pcurves that measure clean — both
/// uses of both closed rim circle edges, with no tolerance widening.
#[test]
fn boolean_hole_rims_ship_no_stale_pcurves() {
    let brep = holed_cube();
    let stale = brep.stale_pcurve_handles();
    assert!(
        stale.is_empty(),
        "boolean output ships {} stale pcurves: {:?}",
        stale.len(),
        stale
            .iter()
            .map(|heh| {
                let report = brep.measure_half_edge_pcurve_deviation(*heh);
                (
                    *heh,
                    report.map(|r| (r.max_deviation, r.trace_direction_agrees)),
                )
            })
            .collect::<Vec<_>>()
    );
}

/// `fix_stale_pcurves` on the reproducer must converge: no repair
/// pass may leave (or reintroduce) staleness it would re-flag on the
/// next pass — the clear+recompute+revert cycle this class caused ran
/// on every whole-model tail pass, forever.
#[test]
fn fix_stale_converges_on_boolean_hole_rims() {
    let mut brep = holed_cube();
    brep.fix_stale_pcurves();
    let stale = brep.stale_pcurve_handles();
    assert!(
        stale.is_empty(),
        "fix_stale_pcurves left {} stale pcurves",
        stale.len()
    );
}

/// The rim pcurves must be genuinely exact — far below the stale
/// gate, not hovering under it. Guards against the fix regressing
/// into a barely-passing fit.
#[test]
fn rim_pcurve_deviation_is_exact_not_marginal() {
    let brep = holed_cube();
    for (eh, edge) in brep.edges.iter() {
        let Some(curve) = edge.curve.as_ref() else {
            continue;
        };
        // Closed rim edges only.
        use popkern_core::geometry::curve::CurveEval;
        let start = curve.geometry.point_at(curve.t_start);
        let end = curve.geometry.point_at(curve.t_end);
        if start.distance_to(&end) > 1e-9 {
            continue;
        }
        let report = brep
            .measure_edge_pcurve_deviation(eh)
            .expect("rim edge report");
        assert!(
            report.max_deviation < 1.0e-9,
            "rim edge {eh:?} deviation {:.3e} should be exact",
            report.max_deviation
        );
    }
}
