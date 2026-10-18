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

/// True when the record leaves its crossings' boundary by more than
/// the linear tolerance in 3D between two consecutive crossings: they
/// are two junctions, however close in the chart (PK-788). The
/// degenerate-hull rescue clusters crossings within a chart-relative
/// band (1e-3 of the record's uv extent — record-fit slop), which on a
/// plane chart swallowed real geometry: a cap 1e-5 inside a ring's
/// extreme cuts the ring's end disc in a chord 0.068 long whose ends
/// sat inside a 0.167 band, and the chord dropped as a tangent touch.
/// A boundary the ring cannot place (an inner-loop half-edge) is no
/// evidence either way.
pub(super) fn crossings_part_in_3d(
    brep: &BRep,
    face: FaceHandle,
    curve_uv: &NurbsCurve,
    crossings: &[UvCrossing],
    tolerance: &Tolerance,
) -> bool {
    let Some(surface) = brep
        .faces
        .get(face)
        .and_then(|f| f.surface.as_ref())
        .map(|s| &s.geometry)
    else {
        return false;
    };
    let Some(boundary) = crate::clip::face_uv_boundary(brep, face) else {
        return false;
    };
    crossings.windows(2).any(|w| {
        boundary
            .iter()
            .find(|(h, _)| *h == w[0].boundary_he)
            .is_some_and(|(_, pc)| {
                pair_parts(surface, &tight_nurbs(pc, tolerance), curve_uv, w, tolerance)
            })
    })
}

/// `pcurve` as NURBS within a ten-thousandth of the linear tolerance
/// (circles in a plane chart are elliptical pcurves, fitted adaptively;
/// the default conversion's slop moved shallow-crossing signs by 1e-8).
fn tight_nurbs(
    pcurve: &popkern_core::geometry::pcurve::Pcurve,
    tolerance: &Tolerance,
) -> NurbsCurve {
    pcurve.to_nurbs_within(tolerance.linear() * 1e-4).0
}

/// [`crossings_part_in_3d`] for one consecutive pair on `pcurve`.
fn pair_parts(
    surface: &popkern_core::geometry::surface::SurfaceGeometry,
    pcurve: &NurbsCurve,
    curve_uv: &NurbsCurve,
    pair: &[UvCrossing],
    tolerance: &Tolerance,
) -> bool {
    use popkern_core::geometry::surface::SurfaceEval as _;
    let q = curve_uv.point_at(0.5 * (pair[0].t_curve + pair[1].t_curve));
    let Some(tb) = crate::curve_projection::project_onto_curve_dense(pcurve, &q) else {
        return false;
    };
    let b = pcurve.point_at(tb);
    let (pq, pb) = (surface.point_at(q.x, q.y), surface.point_at(b.x, b.y));
    let scale = [pq, pb]
        .iter()
        .map(|p| p.x.abs().max(p.y.abs()).max(p.z.abs()))
        .fold(0.0, f64::max);
    pq.distance_to(&pb) >= tolerance.tangency_band(scale)
}

/// Bisection steps when locating a shallow crossing's sign change.
const RUN_BISECTIONS: usize = 60;

/// Resolve RUNS of same-edge reports into the contact they describe
/// (PK-887). Where a record crosses a boundary at a shallow angle the
/// intersector reports a hit all along the stretch within tolerance
/// of it (measured: a cap one tolerance inside a ring cuts the ring's
/// end disc in a chord crossing the rim at 5.9e-4 rad, and reports
/// spread 3.6e-3 past each true crossing — the clip kept the chord
/// out to the outermost report and its ends stranded). A run is a
/// sequence of reports between which the record stays within the
/// tangency band of that edge in 3D; when the record's signed side of
/// the edge just beyond the run differs at its two ends, the run is one
/// transversal crossing at the bisected sign change. Same-side runs (a
/// graze, or a record running along the edge) keep their reports for
/// the tuned clip arms that read them. Non-periodic charts
/// only: a periodic boundary's pcurve branch need not share the
/// record's, and its signed side is not a chart-local fact.
pub(super) fn resolve_contact_runs(
    brep: &BRep,
    face: FaceHandle,
    curve_uv: &NurbsCurve,
    crossings: Vec<UvCrossing>,
    tolerance: &Tolerance,
) -> Vec<UvCrossing> {
    let periodic = brep
        .faces
        .get(face)
        .and_then(|f| f.surface.as_ref())
        .map(|s| popkern_core::geometry::surface::SurfaceEval::parameter_periodicity(&s.geometry))
        .is_none_or(|(u, v)| u.is_some() || v.is_some());
    if periodic || crossings.len() < 2 {
        return crossings;
    }
    let Some(surface) = brep
        .faces
        .get(face)
        .and_then(|f| f.surface.as_ref())
        .map(|s| &s.geometry)
    else {
        return crossings;
    };
    let Some(boundary) = crate::clip::face_uv_boundary(brep, face) else {
        return crossings;
    };
    let mut tight: std::collections::HashMap<popkern_core::topology::HalfEdgeHandle, NurbsCurve> =
        std::collections::HashMap::new();
    for c in &crossings {
        if let Some((h, pc)) = boundary.iter().find(|(h, _)| *h == c.boundary_he) {
            tight
                .entry(*h)
                .or_insert_with(|| tight_nurbs(pc, tolerance));
        }
    }
    let (ct0, ct1) = curve_uv.parameter_range();
    let mut out: Vec<UvCrossing> = Vec::with_capacity(crossings.len());
    let mut i = 0;
    while i < crossings.len() {
        let mut j = i + 1;
        while j < crossings.len()
            && crossings[j].boundary_he == crossings[i].boundary_he
            && tight.get(&crossings[j].boundary_he).is_some_and(|pc| {
                !pair_parts(surface, pc, curve_uv, &crossings[j - 1..=j], tolerance)
            })
        {
            j += 1;
        }
        // Reports at one parameter were already judged (collapsed or kept
        // as a graze pair) by the duplicate collapse; only a run SPREAD
        // along the record is a shallow crossing's stretch.
        let spread = crossings[j - 1].t_curve - crossings[i].t_curve >= (ct1 - ct0).abs() * 1e-5;
        let resolved = (j - i >= 2 && spread)
            .then(|| {
                let pcurve = tight.get(&crossings[i].boundary_he)?;
                let side = |t: f64| -> Option<(f64, f64)> {
                    let q = curve_uv.point_at(t);
                    let tb = crate::curve_projection::project_onto_curve_dense(pcurve, &q)?;
                    let (b, tan) = (pcurve.point_at(tb), pcurve.tangent_at(tb));
                    Some((tan.x * (q.y - b.y) - tan.y * (q.x - b.x), tb))
                };
                let (lo, hi) = (crossings[i].t_curve, crossings[j - 1].t_curve);
                let reach = (hi - lo).max((ct1 - ct0).abs() * 1e-9);
                let before = if i > 0 {
                    0.5 * (crossings[i - 1].t_curve + lo)
                } else {
                    ct0
                };
                let after = if j < crossings.len() {
                    0.5 * (hi + crossings[j].t_curve)
                } else {
                    ct1
                };
                let (mut a, mut b) = ((lo - reach).max(before), (hi + reach).min(after));
                let (sa, sb) = (side(a)?.0, side(b)?.0);
                if sa == 0.0 || sb == 0.0 {
                    return None;
                }
                if sa.signum() == sb.signum() {
                    // Same side beyond the run: a graze, or a record
                    // running along the edge — the tuned clip arms own
                    // those reports as they are.
                    return None;
                }
                for _ in 0..RUN_BISECTIONS {
                    let m = 0.5 * (a + b);
                    let sm = side(m)?.0;
                    if sm.signum() == sa.signum() {
                        a = m;
                    } else {
                        b = m;
                    }
                }
                let t = 0.5 * (a + b);
                Some(vec![UvCrossing {
                    t_curve: t,
                    boundary_he: crossings[i].boundary_he,
                    t_boundary: side(t)?.1,
                }])
            })
            .flatten();
        match resolved {
            Some(r) => out.extend(r),
            None => out.extend_from_slice(&crossings[i..j]),
        }
        i = j;
    }
    out
}
