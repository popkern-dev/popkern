//! Fitted terminal-CONIC trim for swept blend strips (PK-189).
//!
//! The swept-strip terminal settlement (`terminal_arc_fit_trim`)
//! verifies its section is circular and declines honestly when it is
//! not. A conic-profile strip's terminal section IS non-circular by
//! construction, so its chain-end corners were left on the stale
//! sharp corner with off-surface stub curves (measured 1.15e-1
//! CurveSurfaceDeviation on the variable conic cube edge). This
//! module is the conic counterpart: fit the terminal rational
//! quadratic in the planar cap chart — the tangent-triangle apex is
//! the STALE CORNER itself for plane supports, and the shoulder
//! weight comes from the deepest sample along the chord-to-apex
//! direction — verify every sample lies on that conic, settle the
//! corner onto the conic's shoulder point, and author the two stub
//! curves as exact rational-quadratic halves (homogeneous
//! de Casteljau split at t = 1/2). The construction is exact; the
//! only tolerance is the explicit sample residual gate shared with
//! the circular settlement.

use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::geometry::{CurveEval, CurveGeometry, CurveOnEdge, NurbsCurve, Point3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, HalfEdgeHandle, VertexHandle};

use crate::types::TerminalArcTrimDeclineReason;

use super::terminal_arc_author::StagedStubArc;
use super::terminal_chord_trim::StripPlane;

/// The fitted terminal conic: tangent-triangle control points in 3-D
/// (with their cap-chart images) and the shoulder weight.
pub(super) struct TerminalConic {
    p0: Point3,
    apex: Point3,
    p2: Point3,
    w: f64,
    /// Settled corner: the conic's shoulder point C(1/2).
    pub(super) shoulder: Point3,
}

/// Fit the terminal conic through the two anchors with the stale
/// corner as the tangent-triangle apex; the weight comes from the
/// deepest sample. Declines (`None`) when the triangle is degenerate,
/// the fitted rho leaves (0, 1), or any sample strays from the conic
/// beyond the residual gate.
pub(super) fn fit_terminal_conic(
    brep: &BRep,
    corner: VertexHandle,
    anchor_in: Point3,
    anchor_out: Point3,
    samples: &[Point3],
    tolerance: &Tolerance,
) -> Option<TerminalConic> {
    let trace = std::env::var_os("PK_TERMINAL_ARC_TRACE").is_some();
    let apex = brep.vertices.get(corner)?.position;
    let chord = anchor_out - anchor_in;
    let scale = chord.length().max(1e-12);
    let rise = apex - (anchor_in + chord * 0.5);
    let rise_len2 = rise.dot(&rise);
    if rise_len2 <= (1e-6 * scale) * (1e-6 * scale) {
        return None; // apex on the chord — no tangent triangle
    }
    // Per-sample weight via barycentric coordinates in the tangent
    // triangle (P0, A, P2): every point of the conic satisfies
    // τ1² = 4·w²·τ0·τ2, so each interior sample determines
    // w = τ1 / (2·√(τ0·τ2)). The median is robust to the two
    // near-endpoint samples where τ0·τ2 → 0.
    let plane_normal = chord.cross(&rise);
    let e1 = chord * (1.0 / scale);
    let e2 = {
        let n = plane_normal.cross(&e1);
        let l = n.length();
        if l <= 1e-12 {
            return None;
        }
        n * (1.0 / l)
    };
    let to_2d = |p: Point3| -> (f64, f64) {
        let d = p - anchor_in;
        (d.dot(&e1), d.dot(&e2))
    };
    let (a2, b2, c2) = (to_2d(anchor_in), to_2d(apex), to_2d(anchor_out));
    let det = (b2.0 - a2.0) * (c2.1 - a2.1) - (c2.0 - a2.0) * (b2.1 - a2.1);
    if det.abs() <= 1e-12 * scale * scale {
        return None;
    }
    let mut weight_samples: Vec<f64> = Vec::with_capacity(samples.len());
    for p in samples {
        let q = to_2d(*p);
        let t1 = ((q.0 - a2.0) * (c2.1 - a2.1) - (c2.0 - a2.0) * (q.1 - a2.1)) / det;
        let t2 = ((b2.0 - a2.0) * (q.1 - a2.1) - (q.0 - a2.0) * (b2.1 - a2.1)) / det;
        let t0 = 1.0 - t1 - t2;
        if t0 > 1e-3 && t2 > 1e-3 && t1 > 0.0 {
            weight_samples.push(t1 / (2.0 * (t0 * t2).sqrt()));
        }
    }
    if weight_samples.len() < 3 {
        if trace {
            eprintln!("[tat]   conic fit: too few interior samples for the weight");
        }
        return None;
    }
    weight_samples.sort_by(f64::total_cmp);
    let w = weight_samples[weight_samples.len() / 2];
    let rho = w / (1.0 + w);
    if !(0.02..=0.98).contains(&rho) {
        if trace {
            eprintln!("[tat]   conic fit rho {rho:.4} outside (0.02, 0.98)");
        }
        return None;
    }
    let conic = TerminalConic {
        p0: anchor_in,
        apex,
        p2: anchor_out,
        w,
        shoulder: conic_point(anchor_in, apex, anchor_out, w, 0.5),
    };
    // Residual gate (same shape as the circular settlement): every
    // strip sample must lie on the fitted conic.
    let pos_tol = (tolerance.linear() * 10.0).max(1e-9);
    let gate = pos_tol.max(scale * 1.0e-3);
    let dense: Vec<Point3> = (0..=128)
        .map(|k| conic_point(conic.p0, conic.apex, conic.p2, w, k as f64 / 128.0))
        .collect();
    let mut residual = 0.0_f64;
    for p in samples {
        let mut best = f64::INFINITY;
        for pair in dense.windows(2) {
            let (a, b) = (pair[0], pair[1]);
            let ab = b - a;
            let len2 = ab.dot(&ab).max(1e-300);
            let t = ((*p - a).dot(&ab) / len2).clamp(0.0, 1.0);
            let foot = a + ab * t;
            best = best.min((*p - foot).length());
        }
        residual = residual.max(best);
    }
    if residual > gate {
        if trace {
            let planar_worst = samples
                .iter()
                .map(|p| ((*p - anchor_in).dot(&plane_normal) / plane_normal.length()).abs())
                .fold(0.0_f64, f64::max);
            eprintln!(
                "[tat]   conic fit residual {residual:.3e} exceeds gate {gate:.3e} \
                 (w spread {:.6}..{:.6} median {w:.6}, off-plane worst {planar_worst:.3e})",
                weight_samples.first().copied().unwrap_or(f64::NAN),
                weight_samples.last().copied().unwrap_or(f64::NAN),
            );
        }
        return None;
    }
    Some(conic)
}

/// Evaluate the rational quadratic with control points (p0, apex, p2)
/// and middle weight `w` at parameter `t`.
fn conic_point(p0: Point3, apex: Point3, p2: Point3, w: f64, t: f64) -> Point3 {
    let b0 = (1.0 - t) * (1.0 - t);
    let b1 = 2.0 * t * (1.0 - t) * w;
    let b2 = t * t;
    let denom = b0 + b1 + b2;
    Point3::new(
        (p0.x * b0 + apex.x * b1 + p2.x * b2) / denom,
        (p0.y * b0 + apex.y * b1 + p2.y * b2) / denom,
        (p0.z * b0 + apex.z * b1 + p2.z * b2) / denom,
    )
}

/// One half of the terminal conic as an exact rational-quadratic
/// NURBS (homogeneous de Casteljau split at t = 1/2). `first_half`
/// selects [p0 → shoulder]; the other half is [shoulder → p2].
fn conic_half(c: &TerminalConic, first_half: bool) -> NurbsCurve {
    let w = c.w;
    let w_mid = 0.5 * (1.0 + w);
    // Homogeneous midpoint control points.
    let blend = |a: Point3, b: Point3, wb: f64| -> Point3 {
        Point3::new(
            (a.x + wb * b.x) / (1.0 + wb),
            (a.y + wb * b.y) / (1.0 + wb),
            (a.z + wb * b.z) / (1.0 + wb),
        )
    };
    let q_in = blend(c.p0, c.apex, w); // weight w_mid
    let q_out = blend(c.p2, c.apex, w); // weight w_mid
    let (control_points, weights) = if first_half {
        (vec![c.p0, q_in, c.shoulder], vec![1.0, w_mid, w_mid])
    } else {
        (vec![c.shoulder, q_out, c.p2], vec![w_mid, w_mid, 1.0])
    };
    NurbsCurve {
        degree: 2,
        knots: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
        control_points,
        weights,
    }
}

/// Author one swept-strip stub as an exact conic half plus the
/// cap-side NURBS pcurve (the same rational quadratic expressed in
/// the cap plane's chart). The strip-side pcurve stays untouched —
/// it is the sweep-chart authority, exactly as in the circular
/// settlement.
pub(super) fn compute_swept_stub_conic(
    brep: &BRep,
    he_strip: HalfEdgeHandle,
    conic: &TerminalConic,
    first_half: bool,
    cap: &StripPlane,
    tolerance: &Tolerance,
) -> Result<StagedStubArc, TerminalArcTrimDeclineReason> {
    let trace = std::env::var_os("PK_TERMINAL_ARC_TRACE").is_some();
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

    // Walk direction along the loop is (start → end) for the staged
    // pair the caller passes; the recorded curve follows the PRIMARY
    // half-edge, exactly as the circular settlement does.
    let half = conic_half(conic, first_half);
    let (walk_start, walk_end) = if first_half {
        (conic.p0, conic.shoulder)
    } else {
        (conic.shoulder, conic.p2)
    };
    let curve3 = if primary == he_strip {
        half.clone()
    } else {
        reversed_quadratic(&half)
    };
    let pos_tol = (tolerance.linear() * 10.0).max(1e-9);
    let geometry = CurveGeometry::Nurbs(curve3);
    let coe = CurveOnEdge {
        geometry,
        t_start: 0.0,
        t_end: 1.0,
    };
    let (c0, c1) = (coe.geometry.point_at(0.0), coe.geometry.point_at(1.0));
    let (expect0, expect1) = if primary == he_strip {
        (walk_start, walk_end)
    } else {
        (walk_end, walk_start)
    };
    if (c0 - expect0).length() > pos_tol || (c1 - expect1).length() > pos_tol {
        if trace {
            eprintln!(
                "[tat]   conic stub endpoint miss {:.3e}/{:.3e}",
                (c0 - expect0).length(),
                (c1 - expect1).length()
            );
        }
        return Err(TerminalArcTrimDeclineReason::ArcBranchUnresolved);
    }

    // Cap-side pcurve: the same rational quadratic mapped into the
    // cap plane's chart (u, v, 0) — exact, since the section lies in
    // the cap plane. The cap twin runs opposite the strip side.
    let uv_of = |p: Point3| -> Point3 {
        let d = p - cap.origin;
        Point3::new(d.dot(&cap.u_axis), d.dot(&cap.v_axis), 0.0)
    };
    let mut cap_curve = conic_half(conic, first_half);
    for cp in &mut cap_curve.control_points {
        *cp = uv_of(*cp);
    }
    // Twin direction: end → start of the walk.
    let cap_curve = reversed_quadratic(&cap_curve);
    let cap_pcurve = Pcurve::from_nurbs(cap_curve);
    let cap_p0 = cap_pcurve.point_at(cap_pcurve.t_start);
    let want0 = uv_of(walk_end);
    if ((cap_p0.x - want0.x).powi(2) + (cap_p0.y - want0.y).powi(2)).sqrt() > pos_tol {
        if trace {
            eprintln!("[tat]   conic cap pcurve start miss");
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

/// Reverse a single-span rational quadratic (control points and
/// weights mirrored — the exact reversed parameterization).
fn reversed_quadratic(c: &NurbsCurve) -> NurbsCurve {
    NurbsCurve {
        degree: c.degree,
        knots: c.knots.clone(),
        control_points: c.control_points.iter().rev().copied().collect(),
        weights: c.weights.iter().rev().copied().collect(),
    }
}
