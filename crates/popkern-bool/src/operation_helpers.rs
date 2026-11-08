//! Helper functions for the boolean pipeline orchestrator: the BRep
//! copy/combine utilities, containment verification and record
//! tracing. The no-intersection dispatch
//! lives in `operation_no_intersection.rs` (500-line cap).

use std::collections::HashMap;

use popkern_core::topology::BRep;

use crate::error::BoolError;

/// Combine two disjoint BReps into a single BRep by copying all topology
/// from `b` into `a`.
///
/// The two inputs should not share geometry — this is a structural
/// concatenation that copies handles from `b` into fresh slotmap slots
/// in `a` without welding at intersection points. Used by the sheet
/// dispatcher to implement non-overlapping `SheetSheetUnion`.
pub(crate) fn combine_breps(a: BRep, b: BRep) -> Result<BRep, BoolError> {
    copy_brep_into(a, b, false)
}

/// Copy all topology from `b` into `a`, optionally flipping B's face normals.
///
/// When `flip_b_normals` is true, B's material side is mirrored before the
/// copy (loops reversed together with their orientation flags), producing an
/// inner void shell for hollow-solid subtraction. Flipping the flags alone
/// left every void face winding-dishonest — its loop still wound as the
/// tool's outside while its flag said the void's (PK-879).
pub(crate) fn copy_brep_into(a: BRep, b: BRep, flip_b_normals: bool) -> Result<BRep, BoolError> {
    copy_brep_into_mapped(a, b, flip_b_normals).map(|(body, _)| body)
}

/// Where each of `b`'s faces went in a [`copy_brep_into_mapped`].
pub(crate) type FaceMap =
    HashMap<popkern_core::topology::FaceHandle, popkern_core::topology::FaceHandle>;

/// [`copy_brep_into`], also returning where each of `b`'s faces went.
pub(crate) fn copy_brep_into_mapped(
    mut a: BRep,
    mut b: BRep,
    flip_b_normals: bool,
) -> Result<(BRep, FaceMap), BoolError> {
    if flip_b_normals {
        b.reverse_body_orientation();
    }
    let map = a.append_body(&b);
    Ok((a, map.faces))
}

/// the boolean trace topic helper: endpoint/chord dump of every intersection
/// record (a fragmented closed curve shows up as short-chord copies).
pub(crate) fn trace_records(records: &[crate::intersect::IntersectionRecord], n_overlaps: usize) {
    eprintln!("  [exec] {} records, {n_overlaps} overlaps", records.len());
    for (i, r) in records.iter().enumerate() {
        use popkern_core::geometry::curve::CurveEval as _;
        let c = &r.curve.curve_3d;
        let (t0, t1) = c.parameter_range();
        let s = c.point_at(t0);
        let e = c.point_at(t1);
        let uv_span = |c: &popkern_core::geometry::NurbsCurve| -> (f64, f64) {
            let (t0, t1) = c.parameter_range();
            let (mut lo, mut hi) = (f64::MAX, f64::MIN);
            for k in 0..=16 {
                let u = c.point_at(t0 + (t1 - t0) * k as f64 / 16.0).x;
                lo = lo.min(u);
                hi = hi.max(u);
            }
            (lo, hi)
        };
        let (u1lo, u1hi) = uv_span(&r.curve.curve_uv1);
        let (u2lo, u2hi) = uv_span(&r.curve.curve_uv2);
        {
            let (a0, a1) = r.curve.curve_uv1.parameter_range();
            let p0 = r.curve.curve_uv1.point_at(a0);
            let p1 = r.curve.curve_uv1.point_at(a1);
            eprintln!(
                "    [rec {i} uv1-ends] ({:.6},{:.6}) -> ({:.6},{:.6})",
                p0.x, p0.y, p1.x, p1.y
            );
        }
        eprintln!(
            "    [rec {i}] ({:?},{:?}) ({:.6},{:.6},{:.6})->({:.6},{:.6},{:.6}) chord {:.4} u1[{u1lo:.4},{u1hi:.4}] u2[{u2lo:.4},{u2hi:.4}]",
            r.face_a,
            r.face_b,
            s.x,
            s.y,
            s.z,
            e.x,
            e.y,
            e.z,
            s.distance_to(&e)
        );
    }
}

/// Pre-serve retirements and the empty-serve lie gate (PK-465), run
/// after loop-chart normalization and before `post_check`.
pub(super) fn pre_check_retire_and_empty_gate(
    result: &mut BRep,
    certify_empty: &dyn Fn() -> Result<(), crate::error::BoolError>,
    tolerance: &popkern_core::tolerance::Tolerance,
    region: &crate::serve_region::ServeRegion,
) -> Result<usize, crate::error::BoolError> {
    // Weld-residue point-loop slivers retire BEFORE the serve gate:
    // their stale analytic windows measure as loop/shell gaps
    // exactly the weld distance wide (PK-465).
    let point_loops = crate::postprocess::point_loop_retire::retire_point_loop_pairs_in(
        result, tolerance, region,
    );
    if point_loops > 0 && popkern_core::trace::trace_enabled("boolean") {
        eprintln!("[boolean] pre-check: {point_loops} weld point-loop sliver(s) retired");
    }
    let slit_faces = crate::postprocess::slit_retire::retire_pure_slit_faces_in(result, region);
    if slit_faces > 0 && popkern_core::trace::trace_enabled("boolean") {
        eprintln!("[boolean] pre-check: {slit_faces} pure-slit face(s) retired");
    }
    let slit_pairs =
        crate::postprocess::slit_retire::excise_same_face_slit_pairs_in(result, region);
    if slit_pairs > 0 && popkern_core::trace::trace_enabled("boolean") {
        eprintln!("[boolean] pre-check: {slit_pairs} same-face slit pair(s) excised");
    }

    if result.faces.is_empty() {
        certify_empty()?;
    }

    Ok(point_loops + slit_faces + slit_pairs)
}

/// Post-check with ONE junction-reconciliation retry (PK-472): a
/// failing validate re-runs after the bounded geometry repairs in
/// [`crate::postprocess::junction_reconcile`]; a first-time-green
/// result is untouched (bit-identical serves), and a reconciliation
/// that cannot reach the gates declines exactly as before. `carried`
/// edges skip the curve-surface check (see `post_check_gates`).
pub(super) fn validate_with_junction_reconcile_retry(
    result: &mut popkern_core::topology::BRep,
    tolerance: &popkern_core::tolerance::Tolerance,
    carried: &crate::serve_region::CarriedEdges,
    gate: &mut super::postprocess_gate::ValidityGate,
) -> Vec<String> {
    let collect_errors = |report: &crate::validate::post_check::PostCheckReport| -> Vec<String> {
        report
            .findings
            .iter()
            .filter(|f| {
                matches!(
                    f.severity,
                    crate::error::Severity::Error | crate::error::Severity::Critical
                )
            })
            .map(|f| f.message.clone())
            .collect()
    };
    // A body known to validate passes every gate but the curve-surface
    // check without re-measuring them.
    let validated = gate.known_valid();
    let pc_report =
        crate::validate::post_check::post_check_gates(result, tolerance, carried, validated);
    let pc_errors = collect_errors(&pc_report);
    if pc_errors.is_empty() {
        return pc_errors;
    }
    crate::postprocess::junction_reconcile::reconcile_junction_cluster_geometry(result, tolerance);
    gate.edited(1);
    let retry = crate::validate::post_check::post_check_gates(result, tolerance, carried, false);
    let retry_errors = collect_errors(&retry);
    if popkern_core::trace::trace_enabled("boolean") {
        eprintln!(
            "[boolean] junction reconcile retry: {} -> {} validation error(s)",
            pc_errors.len(),
            retry_errors.len()
        );
    }
    retry_errors
}

#[path = "operation_invalid_trace.rs"]
mod operation_invalid_trace;
pub(super) use operation_invalid_trace::trace_invalid_result;
