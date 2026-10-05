//! Crossing-report hygiene for the UV clip (split from
//! `clip_uv_crossings.rs` under the 500-line file cap): drop reports
//! off the query's own parameter window and collapse duplicate
//! reports of one transversal crossing.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::nurbs_curve::NurbsCurve;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};

use super::clip_uv::UvCrossing;

/// Keep the t-sorted `crossings` that lie on `curve_uv`'s parameter
/// window and collapse same-edge duplicate reports of one
/// transversal crossing (see the inline notes).
pub(super) fn collapse_crossing_reports(
    brep: &BRep,
    face: FaceHandle,
    curve_uv: &NurbsCurve,
    mut crossings: Vec<UvCrossing>,
    tolerance: &Tolerance,
) -> Vec<UvCrossing> {
    // Crossings must lie ON the curve: the pcurve-pair intersector
    // works on basis extensions, and hits beyond the curve's own
    // parameter window are not boundary crossings (measured: a fully
    // interior chord got one phantom crossing at t ≈ 4π and was
    // clipped away as Outside, losing socket segments on the second
    // boolean of a part).
    let (ct0, ct1) = curve_uv.parameter_range();
    let span = (ct1 - ct0).abs().max(1e-12);
    let tol_t = span * 1e-6;
    crossings.retain(|c| c.t_curve >= ct0.min(ct1) - tol_t && c.t_curve <= ct0.max(ct1) + tol_t);

    // DUPLICATE-REPORT collapse: the curve×curve intersector can
    // report ONE transversal boundary crossing several times
    // (adjacent spans of a multi-span boundary pcurve each converge
    // onto the shared hit — measured: a thread-groove flank ×
    // end-cap spiral's rim entry came back twice at t = 0.0613 on
    // the 9-CP rational rim circle, the two-crossing hull clipped
    // the record to a 1e-9-wide sliver, and the merge stranded the
    // end-face loop, PK-460). Collapse same-edge same-t reports only
    // when the curve TRANSVERSALLY crosses there — a tangential
    // graze legitimately yields two equal-t crossings and must keep
    // its parity.
    // The grouping gate sits an order above the intersector's own
    // agreement (duplicate reports of one hit land within ~1e-6 of
    // span) and orders below any distinct-crossing spacing (a
    // 1e-3-of-span gate merged two REAL crossings of the radial
    // keyway pin's slot records and stranded the serve). The parity
    // probe keeps the collapse honest: an odd number of true
    // crossings inside the band collapses to one, an even number (a
    // graze) keeps every report, so downstream parity never flips.
    let group_eps = span * 1e-5;
    let mut deduped: Vec<UvCrossing> = Vec::new();
    let mut i = 0;
    while i < crossings.len() {
        let mut j = i + 1;
        while j < crossings.len()
            && crossings[j].boundary_he == crossings[i].boundary_he
            && (crossings[j].t_curve - crossings[i].t_curve).abs() < group_eps
        {
            j += 1;
        }
        if j - i >= 2 {
            let t_lo = crossings[i].t_curve;
            let t_hi = crossings[j - 1].t_curve;
            let lo = ct0.min(ct1);
            let hi = ct0.max(ct1);
            // Interval probes, NOT local offsets: the containment
            // predicate samples the boundary into a coarse polygon
            // whose sagitta band swallows any point near a curved
            // rim (measured 0.23 chart units on an 8-sample circle),
            // so probe at the midpoints between this group and the
            // ADJACENT crossings (or the curve ends) — maximally far
            // from every boundary, the standard interval
            // classification.
            let prev_t = if i > 0 { crossings[i - 1].t_curve } else { lo };
            let next_t = if j < crossings.len() {
                crossings[j].t_curve
            } else {
                hi
            };
            let before = crate::split_inner_loop::point_inside_face_uv(
                brep,
                face,
                curve_uv.point_at(0.5 * (prev_t + t_lo)),
                tolerance,
            );
            let after = crate::split_inner_loop::point_inside_face_uv(
                brep,
                face,
                curve_uv.point_at(0.5 * (t_hi + next_t)),
                tolerance,
            );
            if before != after {
                let t = t_lo;
                if popkern_core::trace::trace_enabled("ring-crossing") {
                    eprintln!(
                        "      [xing-dup] {} duplicate report(s) at t={t:.4} collapsed \
                         (transversal)",
                        j - i - 1
                    );
                }
                deduped.push(crossings[i].clone());
            } else {
                if popkern_core::trace::trace_enabled("ring-crossing") {
                    eprintln!(
                        "      [xing-dup] kept {} report(s) at t {t_lo:.7}..{t_hi:.7} \
                         (no parity flip: before={before} after={after})",
                        j - i
                    );
                }
                deduped.extend_from_slice(&crossings[i..j]);
            }
        } else {
            deduped.push(crossings[i].clone());
        }
        i = j;
    }
    deduped
}
