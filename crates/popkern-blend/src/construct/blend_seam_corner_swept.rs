//! PK-160 swept-strip pair-corner repair.
//!
//! Variable-radius fillets author swept rational (NURBS) strips, so
//! the pair-corner seam between two adjacent fillets has no analytic
//! form: it is the strip×strip surface intersection. This stage
//! reuses the PK-152 tangle recognition and Euler surgery
//! (`blend_seam_corner` / `blend_seam_corner_repair`) and authors the
//! seam by alternating projection: chord samples between the two
//! trusted anchor vertices are projected onto each strip in turn
//! until they converge onto the intersection curve, then a degree-3
//! interpolating curve and per-chart pcurves (hint-seeded uv
//! projection of the same samples) are committed.
//!
//! **Fallback:** any sample that fails to converge onto both strips
//! within tolerance declines the corner (`MarchDidNotConverge`), the
//! caller restores its snapshot, and the stale-but-valid tangle
//! ships as before. **Tolerance:** the caller's blend tolerance —
//! convergence gate `linear × 10`, same as the analytic stages.

use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{Point3, SurfaceEval};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle, HalfEdgeHandle, VertexHandle};

use crate::types::BlendSeamCornerDeclineReason as Reason;

use super::blend_seam_corner::PairCornerCandidate;
use super::blend_seam_corner_march::authoring::{chart_pcurve, chord_params};
use super::blend_seam_corner_repair::{commit_edge_authoring, perform_flap_surgery};

/// Seam samples between the anchors (endpoints included).
const SAMPLES: usize = 17;
/// Alternating-projection refinement rounds per sample. The
/// iteration is linear and near the support rails the strips meet at
/// a shallow angle (both tangent to the shared support plane), so
/// the crawl needs headroom; each round is two hint-seeded Newton
/// projections and the loop breaks on the residual gate.
const ROUNDS: usize = 400;

pub(super) fn repair_pair_corner_swept(
    brep: &mut BRep,
    c: &PairCornerCandidate,
    tolerance: &Tolerance,
) -> Result<(), Reason> {
    let surf_of = |face: FaceHandle| -> Result<SurfaceGeometry, Reason> {
        brep.faces
            .get(face)
            .and_then(|f| f.surface.as_ref())
            .map(|s| s.geometry.clone())
            .ok_or(Reason::EulerRejected)
    };
    let geom_n = surf_of(c.strip_n)?;
    let geom_p = surf_of(c.strip_p)?;
    let pos_of = |v: VertexHandle| -> Result<Point3, Reason> {
        Ok(brep.vertices.get(v).ok_or(Reason::EulerRejected)?.position)
    };
    let a_pos = pos_of(c.anchor_near)?;
    let b_pos = pos_of(c.anchor_far)?;
    let pos_tol = (tolerance.linear() * 10.0).max(1e-9);
    if popkern_core::trace::trace_enabled("blend-seam") {
        eprintln!(
            "[bsc-swept] a off ({:.3e},{:.3e}) b off ({:.3e},{:.3e})",
            surface_distance(&geom_n, &a_pos, None),
            surface_distance(&geom_p, &a_pos, None),
            surface_distance(&geom_n, &b_pos, None),
            surface_distance(&geom_p, &b_pos, None),
        );
    }

    // Mismatched corner radii: the near anchor sits on ONE strip
    // only and the seam exits through the support plane — the W-exit
    // branch.
    let a_on_n = surface_distance(&geom_n, &a_pos, None) <= pos_tol;
    let a_on_p = surface_distance(&geom_p, &a_pos, None) <= pos_tol;
    match (a_on_n, a_on_p) {
        (true, true) => {}
        (true, false) | (false, true) => {
            return super::blend_seam_corner_swept_exit::repair_swept_mismatch(
                brep, c, &geom_n, &geom_p, a_on_n, tolerance,
            );
        }
        (false, false) => return Err(Reason::StubUnresolved),
    }

    // March the seam: chord samples refined onto BOTH strips by
    // alternating projection. Interior samples only — the anchors
    // are trusted fixed endpoints.
    let mut pts: Vec<Point3> = Vec::with_capacity(SAMPLES);
    let mut hint_n: Option<(f64, f64)> = None;
    let mut hint_p: Option<(f64, f64)> = None;
    for i in 0..SAMPLES {
        let f = i as f64 / (SAMPLES - 1) as f64;
        let mut p = Point3::new(
            a_pos.x + f * (b_pos.x - a_pos.x),
            a_pos.y + f * (b_pos.y - a_pos.y),
            a_pos.z + f * (b_pos.z - a_pos.z),
        );
        if i > 0 && i < SAMPLES - 1 {
            p = refine_onto_both(&geom_n, &geom_p, p, &mut hint_n, &mut hint_p, pos_tol)
                .ok_or(Reason::StubUnresolved)?;
            let off_n = surface_distance(&geom_n, &p, hint_n);
            let off_p = surface_distance(&geom_p, &p, hint_p);
            if off_n > pos_tol || off_p > pos_tol {
                if popkern_core::trace::trace_enabled("blend-seam") {
                    eprintln!(
                        "[bsc-swept] sample {i} off_n={off_n:.3e} off_p={off_p:.3e} p=({:.4},{:.4},{:.4})",
                        p.x, p.y, p.z
                    );
                }
                return Err(Reason::StubUnresolved);
            }
        } else {
            // Seed the hints from the anchors so interior samples
            // stay on the corner-local chart branch.
            hint_n = popkern_core::geometry::uv_projection::project_point_to_uv(&geom_n, &p);
            hint_p = popkern_core::geometry::uv_projection::project_point_to_uv(&geom_p, &p);
        }
        pts.push(p);
    }

    let params = chord_params(&pts);
    // The stored curve follows the primary half-edge's direction.
    let primary = brep
        .edges
        .get(c.seam_edge)
        .ok_or(Reason::EulerRejected)?
        .half_edge;
    let (curve_pts, curve_params) = if primary == c.he_n {
        (pts.clone(), params.clone())
    } else {
        let mut rp: Vec<Point3> = pts.clone();
        rp.reverse();
        let rpar: Vec<f64> = params.iter().rev().map(|t| 1.0 - t).collect();
        (rp, rpar)
    };
    let fitted = popkern_core::geometry::fitting::interpolate_curve(&curve_pts, &curve_params, 3)
        .map_err(|_| Reason::StubUnresolved)?;
    let seam_curve = popkern_core::geometry::CurveOnEdge {
        geometry: popkern_core::geometry::CurveGeometry::Nurbs(fitted),
        t_start: 0.0,
        t_end: 1.0,
    };

    // Chart pcurves: hint-seeded projection of the SAME samples into
    // each strip chart, in that side's half-edge direction.
    let staged = stage_side_pcurves(brep, c, &geom_n, &geom_p, &pts, &params)?;

    perform_flap_surgery(brep, c)?;
    commit_edge_authoring(brep, c.seam_edge, seam_curve, staged);

    let dev = brep.measure_edge_deviation(c.seam_edge);
    if dev > pos_tol {
        if popkern_core::trace::trace_enabled("blend-seam") {
            eprintln!("[bsc-swept] seam deviation {dev:.3e} > gate {pos_tol:.3e}");
        }
        return Err(Reason::ValidationFailed);
    }
    if let Err(e) = brep.validate() {
        if popkern_core::trace::trace_enabled("blend-seam") {
            eprintln!("[bsc-swept] validate failed: {e:?}");
        }
        return Err(Reason::ValidationFailed);
    }
    Ok(())
}

/// Alternating-projection refinement of one point onto BOTH
/// surfaces (residual-gated linear iteration; hint-seeded Newton
/// projections keep the walk on the corner-local chart branch).
/// Returns the refined point; the caller re-checks residuals.
pub(super) fn refine_onto_both(
    geom_a: &SurfaceGeometry,
    geom_b: &SurfaceGeometry,
    start: Point3,
    hint_a: &mut Option<(f64, f64)>,
    hint_b: &mut Option<(f64, f64)>,
    pos_tol: f64,
) -> Option<Point3> {
    let mut p = start;
    for _ in 0..ROUNDS {
        let (ua, va) = project_seeded(geom_a, &p, *hint_a)?;
        p = geom_a.point_at(ua, va);
        let (ub, vb) = project_seeded(geom_b, &p, *hint_b)?;
        let q = geom_b.point_at(ub, vb);
        let moved = p.distance_to(&q);
        p = q;
        *hint_a = Some((ua, va));
        *hint_b = Some((ub, vb));
        if moved < pos_tol * 1e-3 {
            break;
        }
        if moved < pos_tol * 0.05 && surface_distance(geom_a, &p, *hint_a) <= pos_tol * 0.5 {
            break;
        }
    }
    Some(p)
}

pub(super) fn project_seeded(
    geom: &SurfaceGeometry,
    p: &Point3,
    hint: Option<(f64, f64)>,
) -> Option<(f64, f64)> {
    match hint {
        Some((u, v)) => {
            popkern_core::geometry::uv_projection::project_point_to_uv_near(geom, p, u, v)
                .or_else(|| popkern_core::geometry::uv_projection::project_point_to_uv(geom, p))
        }
        None => popkern_core::geometry::uv_projection::project_point_to_uv(geom, p),
    }
}

pub(super) fn surface_distance(
    geom: &SurfaceGeometry,
    p: &Point3,
    hint: Option<(f64, f64)>,
) -> f64 {
    match project_seeded(geom, p, hint) {
        Some((u, v)) => p.distance_to(&geom.point_at(u, v)),
        None => f64::INFINITY,
    }
}

/// Project the marched samples into both strip charts (half-edge
/// direction per side: `he_n` runs anchor_near → anchor_far, `he_p`
/// the reverse) and fit the side pcurves.
fn stage_side_pcurves(
    _brep: &BRep,
    c: &PairCornerCandidate,
    geom_n: &SurfaceGeometry,
    geom_p: &SurfaceGeometry,
    pts: &[Point3],
    params: &[f64],
) -> Result<Vec<(FaceHandle, HalfEdgeHandle, Pcurve)>, Reason> {
    let rev_pts: Vec<Point3> = pts.iter().rev().copied().collect();
    let rev_params: Vec<f64> = params.iter().rev().map(|t| 1.0 - t).collect();
    let mut staged = Vec::with_capacity(2);
    for (face, he, geom, side_pts, side_params) in [
        (c.strip_n, c.he_n, geom_n, pts, params),
        (
            c.strip_p,
            c.he_p,
            geom_p,
            rev_pts.as_slice(),
            rev_params.as_slice(),
        ),
    ] {
        let mut uv: Vec<(f64, f64)> = Vec::with_capacity(side_pts.len());
        let mut hint: Option<(f64, f64)> = None;
        for p in side_pts {
            let (u, v) = project_seeded(geom, p, hint).ok_or(Reason::PcurveClosureFailed)?;
            hint = Some((u, v));
            uv.push((u, v));
        }
        let pcurve = chart_pcurve(&uv, side_params)?;
        staged.push((face, he, pcurve));
    }
    Ok(staged)
}
