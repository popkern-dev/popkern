//! Fan edges of a strip's fan end (PK-530, PK-534): where the strip
//! crosses each fan edge, and the edge re-authored from its crossing.
//!
//! A fan edge is whatever the fan faces meet along:
//! - a plane ∩ plane line at a pyramid apex;
//! - a blend-on-blend seam, where a vertical fillet runs into two top
//!   fillets that already meet along their PK-152 ellipse seam.
//!
//! Lines against cylinder strips cross in closed form. Any other pair
//! is found as the first sign change of the strip's signed distance
//! along the edge's own curve from the corner, then bisected. The root
//! is accepted only where the point genuinely touches the strip. The
//! edge keeps its curve, trimmed to start at the crossing. Its pcurves
//! on both flanks are re-fitted from the trimmed curve (lines get
//! exact linear plane pcurves), and each is gated by the caller's
//! validation.

use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{CurveEval, CurveGeometry, CurveOnEdge, Point3, SurfaceEval};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle, HalfEdgeHandle};

use crate::types::TerminalArcTrimDeclineReason as Reason;

use super::super::blend_seam_corner_swept::{project_seeded, surface_distance};
use super::super::cap_insert_boundary::projected_pcurve_sampled;
use super::super::terminal_arc_trim::StripCylinder;

/// Bracketing samples along a fan edge, and bisection steps.
const BRACKET_SAMPLES: usize = 64;
const BISECTIONS: usize = 60;

fn face_of(brep: &BRep, he: HalfEdgeHandle) -> Option<FaceHandle> {
    Some(brep.loops.get(brep.half_edges.get(he)?.loop_)?.face)
}

/// The fan edge's curve as a function of the fraction `s` from the
/// corner (`s = 0`) to its far end (`s = 1`), with the edge.
fn fan_curve(brep: &BRep, f: HalfEdgeHandle) -> Option<(CurveOnEdge, bool)> {
    let edge = brep.half_edges.get(f)?.edge;
    let data = brep.edges.get(edge)?;
    // The stored curve follows the primary half-edge; `f` leaves the
    // corner, so the curve runs from the corner iff `f` is primary.
    Some((data.curve.clone()?, data.half_edge == f))
}

fn at_fraction(curve: &CurveOnEdge, from_corner: bool, s: f64) -> Point3 {
    let f = if from_corner { s } else { 1.0 - s };
    curve
        .geometry
        .point_at(curve.t_start + (curve.t_end - curve.t_start) * f)
}

/// Signed distance from `p` to the strip.
fn signed_distance(strip: &SurfaceGeometry, p: Point3) -> Option<f64> {
    let (u, v) = project_seeded(strip, &p, None)?;
    let foot = strip.point_at(u, v);
    let n = strip.normal_at(u, v).normalized()?;
    Some((p - foot).dot(&n))
}

/// Where the strip crosses the fan edge `f` (leaving the corner),
/// nearest the corner.
pub(super) fn fan_edge_crossing(
    brep: &BRep,
    f: HalfEdgeHandle,
    corner: Point3,
    cylinder: Option<&StripCylinder>,
    strip: &SurfaceGeometry,
    pos_tol: f64,
) -> Option<Point3> {
    let (curve, from_corner) = fan_curve(brep, f)?;
    if let (Some(cyl), CurveGeometry::Line { .. }) = (cylinder, &curve.geometry) {
        let far = at_fraction(&curve, from_corner, 1.0);
        return line_cylinder_crossing(cyl, corner, far);
    }
    let at = |s: f64| at_fraction(&curve, from_corner, s);
    let mut prev_s = 0.0;
    let mut prev = signed_distance(strip, at(0.0))?;
    for k in 1..=BRACKET_SAMPLES {
        let s = k as f64 / BRACKET_SAMPLES as f64;
        let Some(fs) = signed_distance(strip, at(s)) else {
            continue;
        };
        if prev * fs <= 0.0 {
            let (mut lo, mut hi, mut f_lo) = (prev_s, s, prev);
            for _ in 0..BISECTIONS {
                let mid = 0.5 * (lo + hi);
                let f_mid = signed_distance(strip, at(mid))?;
                if f_lo * f_mid <= 0.0 {
                    hi = mid;
                } else {
                    lo = mid;
                    f_lo = f_mid;
                }
            }
            let q = at(0.5 * (lo + hi));
            return (surface_distance(strip, &q, None) <= pos_tol).then_some(q);
        }
        prev_s = s;
        prev = fs;
    }
    None
}

/// The point nearest `from` along the segment to `to` where the line
/// meets the cylinder.
fn line_cylinder_crossing(cyl: &StripCylinder, from: Point3, to: Point3) -> Option<Point3> {
    let axis = cyl.axis.normalized()?;
    let d = to - from;
    let w = from - cyl.origin;
    let d_perp = d - axis * d.dot(&axis);
    let w_perp = w - axis * w.dot(&axis);
    let a = d_perp.dot(&d_perp);
    let b = 2.0 * d_perp.dot(&w_perp);
    let cc = w_perp.dot(&w_perp) - cyl.radius * cyl.radius;
    let disc = b * b - 4.0 * a * cc;
    if a < 1e-18 || disc < 0.0 {
        return None;
    }
    let sq = disc.sqrt();
    let mut roots = [(-b - sq) / (2.0 * a), (-b + sq) / (2.0 * a)];
    roots.sort_by(f64::total_cmp);
    let t = roots.into_iter().find(|t| *t > 1e-12 && *t < 1.0 - 1e-12)?;
    Some(from + d * t)
}

/// Re-author the fan edge (outgoing half-edge `f`, now starting at its
/// crossing `start`): its curve trimmed to `[start, far end]`, with
/// pcurves on both flanks.
pub(super) fn reauthor_fan_edge(
    brep: &mut BRep,
    f: HalfEdgeHandle,
    start: Point3,
    flanks: (FaceHandle, FaceHandle),
    tolerance: &Tolerance,
) -> Result<(), Reason> {
    let (curve, from_corner) = fan_curve(brep, f).ok_or(Reason::FanEndRejected)?;
    let twin = brep.half_edges.get(f).ok_or(Reason::FanEndRejected)?.twin;
    let edge = brep.half_edges.get(f).ok_or(Reason::FanEndRejected)?.edge;
    let far = at_fraction(&curve, from_corner, 1.0);
    let trimmed = match &curve.geometry {
        CurveGeometry::Line { .. } => {
            let (p0, p1) = if from_corner {
                (start, far)
            } else {
                (far, start)
            };
            CurveOnEdge::line_segment(p0, p1, tolerance.linear()).ok_or(Reason::FanEndRejected)?
        }
        _ => {
            // Same geometry, the corner-side end moved to the crossing's
            // parameter.
            let t_start = curve
                .geometry
                .closest_point(&start)
                .ok_or(Reason::FanEndRejected)?;
            if from_corner {
                CurveOnEdge {
                    t_start,
                    ..curve.clone()
                }
            } else {
                CurveOnEdge {
                    t_end: t_start,
                    ..curve.clone()
                }
            }
        }
    };
    let pos_tol = (tolerance.linear() * 10.0).max(1e-9);
    if (trimmed.geometry.point_at(if from_corner {
        trimmed.t_start
    } else {
        trimmed.t_end
    }) - start)
        .length()
        > pos_tol
    {
        return Err(Reason::FanEndRejected);
    }
    let primary = brep
        .edges
        .get(edge)
        .ok_or(Reason::FanEndRejected)?
        .half_edge;
    // A fitted fan edge (a marched blend seam) carries its own band;
    // the trimmed edge is held to it, never tighter than `pos_tol`.
    let band = brep
        .edges
        .get(edge)
        .and_then(|e| e.tolerance)
        .unwrap_or(0.0)
        .max(pos_tol);
    let mut sides: Vec<(FaceHandle, HalfEdgeHandle, SurfaceGeometry)> = Vec::with_capacity(2);
    for h in [f, twin] {
        let face = face_of(brep, h).ok_or(Reason::FanEndRejected)?;
        if face != flanks.0 && face != flanks.1 {
            return Err(Reason::FanEndRejected);
        }
        let surface = brep
            .faces
            .get(face)
            .and_then(|x| x.surface.as_ref())
            .map(|x| x.geometry.clone())
            .ok_or(Reason::FanEndRejected)?;
        sides.push((face, h, surface));
    }
    if let Some(e) = brep.edges.get_mut(edge) {
        e.curve = Some(trimmed.clone());
    }
    // Projected pcurves densify (17 -> 65 samples) until the edge
    // holds its band; lines on planes are exact at once.
    for samples in [17usize, 33, 65] {
        for (face, h, surface) in &sides {
            let pcurve = match (&trimmed.geometry, surface) {
                (
                    CurveGeometry::Line { .. },
                    SurfaceGeometry::Plane {
                        origin,
                        u_axis,
                        v_axis,
                        ..
                    },
                ) => {
                    let uv = |p: Point3| ((p - *origin).dot(u_axis), (p - *origin).dot(v_axis));
                    let (a, b) = if *h == f { (start, far) } else { (far, start) };
                    Pcurve::linear(uv(a), uv(b), 0.0, 1.0)
                }
                _ => projected_pcurve_sampled(surface, &trimmed, *h != primary, samples)
                    .map_err(|_| Reason::FanEndRejected)?,
            };
            if let Some(x) = brep.half_edges.get_mut(*h) {
                x.pcurve = Some(pcurve.clone());
            }
            brep.set_face_use_pcurve(*face, *h, pcurve);
        }
        if brep.measure_edge_deviation(edge) <= band {
            return Ok(());
        }
    }
    Err(Reason::FanEndRejected)
}
