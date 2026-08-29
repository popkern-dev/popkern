//! Terminal ISO-CURVE settle for swept (NURBS) blend strips whose
//! terminal cross-section is neither circular nor conic — the G2
//! degree-5 section class (PK-375).
//!
//! `repair_corner_swept` settles a stale chain-end corner by fitting
//! the terminal circle, then (PK-189) the terminal conic. A
//! curvature-continuous section is neither: measured on a single
//! r=0.1 G2 edge fillet of the unit cube, the circular residual is
//! 9.0e-4 and the conic-fit residual 1.9e-3 against a 1.4e-4 gate, so
//! both arms declined honestly and the corner vertex stayed at the
//! cube corner `(√2−1)·r`-class off the strip — check 9 reported a
//! 3.6e-2 curve-vs-surface lie on all four stub edges, and
//! `Part::fillet_with(G2)` never served on a cube.
//!
//! **Repair, by construction.** The sweep interpolates its sections
//! exactly, so the strip surface's iso curve at the terminal station
//! IS the terminal section, and [`NurbsSurface::iso_curve_at_u`] /
//! [`NurbsSurface::iso_curve_at_v`] extract it as an exact NURBS curve
//! (de Boor on the control grid — no sampling, no fit). The corner is
//! settled onto the iso point at the stubs' shared chart corner (the
//! deepest section station, where the strip-side pcurves already
//! meet), both stub 3-D curves become the exact iso-curve pieces
//! between the trusted anchors and that point, and the cap-side
//! pcurves are the same pieces mapped affinely into the planar cap
//! chart (exact, since the section lies in the cap plane — verified
//! by the caller before any arm runs). Strip-side pcurves are kept:
//! they are the sweep-chart authority that located the section.
//!
//! **Arm ordering.** Circle, conic, then iso: an analytic spelling is
//! preferred whenever it is exact, so this arm only serves sections
//! the analytic arms cannot represent. **Invariants:** topology
//! handles and counts; support trims. **Fallback:** any anchor miss
//! (the iso does not pass through the trusted anchors within the
//! repair band — a stale chart) keeps the honest decline.
//! **Tolerance:** the caller's linear tolerance widened by the shared
//! 10× repair band, local to this settle.

use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::geometry::{
    CurveEval, CurveGeometry, CurveOnEdge, NurbsCurve, NurbsSurface, Point3,
};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle, HalfEdgeHandle};

use super::terminal_arc_author::StagedStubArc;
use super::terminal_arc_trim::CornerCandidate;
use super::terminal_chord_trim::StripPlane;
use crate::types::TerminalArcTrimDeclineReason;

/// The terminal station of a swept strip: which chart direction is
/// held fixed and at what value, plus the free-parameter stations of
/// the two stubs (`free_in` at the in-anchor, `free_deep` at the
/// shared chart corner, `free_out` at the out-anchor).
pub(super) struct TerminalIsoStation {
    /// `true` when `u` is fixed (the iso is parameterized by `v`).
    pub(super) iso_is_u: bool,
    /// The fixed chart coordinate.
    pub(super) fixed: f64,
    /// Free coordinate at the in-stub's anchor end.
    pub(super) free_in: f64,
    /// Free coordinate at the stubs' shared chart corner.
    pub(super) free_deep: f64,
    /// Free coordinate at the out-stub's anchor end.
    pub(super) free_out: f64,
}

/// Settle `c.corner` onto the strip's exact terminal iso curve and
/// re-author both stubs. Transactional: every piece is computed and
/// verified before anything is committed.
#[allow(clippy::too_many_arguments)]
pub(super) fn repair_corner_iso(
    brep: &mut BRep,
    c: &CornerCandidate,
    strip: &NurbsSurface,
    station: &TerminalIsoStation,
    anchor_in: Point3,
    anchor_out: Point3,
    cap: &StripPlane,
    cap_face: FaceHandle,
    tolerance: &Tolerance,
) -> Result<(), TerminalArcTrimDeclineReason> {
    let trace = popkern_core::trace::trace_enabled("terminal-arc");
    let pos_tol = (tolerance.linear() * 10.0).max(1e-9);
    let iso = if station.iso_is_u {
        strip.iso_curve_at_u(station.fixed)
    } else {
        strip.iso_curve_at_v(station.fixed)
    };
    let deep = iso.point_at(station.free_deep);

    // The section must run monotonically in → deep → out along the
    // iso; anything else is a folded or stale chart.
    let (f_in, f_deep, f_out) = (station.free_in, station.free_deep, station.free_out);
    let monotone = (f_in < f_deep && f_deep < f_out) || (f_in > f_deep && f_deep > f_out);
    if !monotone {
        if trace {
            eprintln!("[tat]   iso stations not monotone: {f_in:.4} {f_deep:.4} {f_out:.4}");
        }
        return Err(TerminalArcTrimDeclineReason::ArcBranchUnresolved);
    }
    // Pieces in walk direction (in → deep, deep → out).
    let piece_in = iso_piece(&iso, f_in, f_deep)?;
    let piece_out = iso_piece(&iso, f_deep, f_out)?;

    // The trusted anchors must lie on the iso (a stale chart would
    // put the terminal station elsewhere).
    let miss_in = (piece_in.point_at(piece_range(&piece_in).0) - anchor_in).length();
    let miss_out = (piece_out.point_at(piece_range(&piece_out).1) - anchor_out).length();
    if miss_in > pos_tol || miss_out > pos_tol {
        if trace {
            eprintln!("[tat]   iso anchor miss {miss_in:.3e}/{miss_out:.3e} (band {pos_tol:.3e})");
        }
        return Err(TerminalArcTrimDeclineReason::ArcBranchUnresolved);
    }

    let staged = [
        compute_iso_stub(brep, c.he_in, &piece_in, anchor_in, deep, cap, pos_tol)?,
        compute_iso_stub(brep, c.he_out, &piece_out, deep, anchor_out, cap, pos_tol)?,
    ];
    if trace {
        eprintln!(
            "[tat]   swept terminal section settled on the exact iso curve (degree {}, {} spans); corner -> {deep:?}",
            iso.degree,
            iso.control_points.len()
        );
    }
    if let Some(v) = brep.vertices.get_mut(c.corner) {
        v.position = deep;
    }
    for stub in staged {
        if let Some(e) = brep.edges.get_mut(stub.edge) {
            e.curve = Some(stub.curve);
        }
        if let Some(tw) = brep.half_edges.get_mut(stub.cap_half_edge) {
            tw.pcurve = Some(stub.cap_pcurve.clone());
        }
        brep.set_face_use_pcurve(cap_face, stub.cap_half_edge, stub.cap_pcurve);
    }
    Ok(())
}

/// Clamped parameter range of a curve from its knot vector.
fn piece_range(curve: &NurbsCurve) -> (f64, f64) {
    let d = curve.degree;
    let n = curve.knots.len();
    (curve.knots[d], curve.knots[n - 1 - d])
}

/// The exact sub-curve of `iso` running from `from` to `to` in that
/// walk direction (reversed when `from > to`).
fn iso_piece(
    iso: &NurbsCurve,
    from: f64,
    to: f64,
) -> Result<NurbsCurve, TerminalArcTrimDeclineReason> {
    let (lo, hi) = (from.min(to), from.max(to));
    let (t0, t1) = piece_range(iso);
    let eps = 1e-12 * (t1 - t0).abs().max(1.0);
    let mut piece = iso.clone();
    if lo > t0 + eps {
        piece = piece
            .split_at(lo)
            .map_err(|_| TerminalArcTrimDeclineReason::ArcBranchUnresolved)?
            .1;
    }
    let (_, p1) = piece_range(&piece);
    if hi < p1 - eps {
        piece = piece
            .split_at(hi)
            .map_err(|_| TerminalArcTrimDeclineReason::ArcBranchUnresolved)?
            .0;
    }
    Ok(if from > to { reversed(&piece) } else { piece })
}

/// Exact reversal: control points and weights mirrored, knots
/// reflected about the parameter range (`k ↦ t0 + t1 − k`).
fn reversed(c: &NurbsCurve) -> NurbsCurve {
    let (t0, t1) = piece_range(c);
    NurbsCurve {
        degree: c.degree,
        knots: c.knots.iter().rev().map(|k| t0 + t1 - k).collect(),
        control_points: c.control_points.iter().rev().copied().collect(),
        weights: c.weights.iter().rev().copied().collect(),
    }
}

/// Author one stub as an exact iso-curve piece (walk direction
/// `start_pos → end_pos`) plus the cap-side NURBS pcurve. The
/// strip-side pcurve stays as-is.
fn compute_iso_stub(
    brep: &BRep,
    he_strip: HalfEdgeHandle,
    piece: &NurbsCurve,
    start_pos: Point3,
    end_pos: Point3,
    cap: &StripPlane,
    pos_tol: f64,
) -> Result<StagedStubArc, TerminalArcTrimDeclineReason> {
    let trace = popkern_core::trace::trace_enabled("terminal-arc");
    let (edge, primary, strip_pcurve, cap_half_edge) = {
        let he = brep
            .half_edges
            .get(he_strip)
            .ok_or(TerminalArcTrimDeclineReason::ArcBranchUnresolved)?;
        let edge = he.edge;
        let primary = brep
            .edges
            .get(edge)
            .ok_or(TerminalArcTrimDeclineReason::ArcBranchUnresolved)?
            .half_edge;
        let strip_face = brep
            .loops
            .get(he.loop_)
            .map(|l| l.face)
            .ok_or(TerminalArcTrimDeclineReason::ArcBranchUnresolved)?;
        let pc = brep
            .face_use_pcurve(strip_face, he_strip)
            .or(he.pcurve.as_ref())
            .cloned()
            .ok_or(TerminalArcTrimDeclineReason::ArcBranchUnresolved)?;
        (edge, primary, pc, he.twin)
    };
    // The recorded curve follows the PRIMARY half-edge.
    let curve3 = if primary == he_strip {
        piece.clone()
    } else {
        reversed(piece)
    };
    let (t_start, t_end) = piece_range(&curve3);
    let (expect0, expect1) = if primary == he_strip {
        (start_pos, end_pos)
    } else {
        (end_pos, start_pos)
    };
    let (c0, c1) = (curve3.point_at(t_start), curve3.point_at(t_end));
    if (c0 - expect0).length() > pos_tol || (c1 - expect1).length() > pos_tol {
        if trace {
            eprintln!(
                "[tat]   iso stub endpoint miss {:.3e}/{:.3e}",
                (c0 - expect0).length(),
                (c1 - expect1).length()
            );
        }
        return Err(TerminalArcTrimDeclineReason::ArcBranchUnresolved);
    }
    let coe = CurveOnEdge {
        geometry: CurveGeometry::Nurbs(curve3),
        t_start,
        t_end,
    };

    // Cap-side pcurve: the same piece mapped into the cap chart
    // (u, v, 0) — an affine map, so the rational curve is exact. The
    // cap twin runs opposite the strip side: end → start of the walk.
    let uv_of = |p: Point3| -> Point3 {
        let d = p - cap.origin;
        Point3::new(d.dot(&cap.u_axis), d.dot(&cap.v_axis), 0.0)
    };
    let mut cap_curve = reversed(piece);
    for cp in &mut cap_curve.control_points {
        *cp = uv_of(*cp);
    }
    let cap_pcurve = Pcurve::from_nurbs(cap_curve);
    let cap_p0 = cap_pcurve.point_at(cap_pcurve.t_start);
    let want0 = uv_of(end_pos);
    if ((cap_p0.x - want0.x).powi(2) + (cap_p0.y - want0.y).powi(2)).sqrt() > pos_tol {
        if trace {
            eprintln!("[tat]   iso cap pcurve start miss");
        }
        return Err(TerminalArcTrimDeclineReason::ArcBranchUnresolved);
    }

    Ok(StagedStubArc {
        edge,
        curve: coe,
        strip_half_edge: he_strip,
        strip_pcurve,
        cap_half_edge,
        cap_pcurve,
    })
}
