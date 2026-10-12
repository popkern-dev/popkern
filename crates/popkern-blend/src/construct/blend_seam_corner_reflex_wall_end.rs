//! Asymmetric reflex mitre ends (PK-742).
//!
//! On a planar boss base the two root strips' wall rails reach the
//! boss's sharp corner edge at one point `T`, where the strips' mitre
//! ends (PK-740). On a curved support they do not: a fin's side plane
//! is offset from radial on a cylinder housing, so the side strip's
//! wall rail meets the sharp edge 0.964 r from the corner while the arc
//! strip's cap rail meets it 1.0003 r out. The mitre cannot end on the
//! edge. Measured against OpenCASCADE's fillet of the same loop: the
//! mitre ends at `P`, where the side strip's straight wall rail meets
//! the arc strip, 6.4e-4 below the sharp edge. The arc strip and the
//! side wall then share a short edge `P → Q`, their intersection, to
//! the arc strip's own rail end `Q` on the sharp edge. The wall face
//! dips below its original edge along that curve.
//!
//! The anchor stage builds that topology but lands `P` on the sharp
//! edge (where a planar boss's `T` would be), off the arc strip. This
//! module finds the true `P` by sliding the vertex along its straight
//! rail to the arc strip (a 1-D signed-distance root; the rail is a
//! tangential contact, so a three-surface Newton would be singular
//! there), and re-authors the vertex's other edges once it has moved.

use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::{CurveGeometry, CurveOnEdge, Point3};
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle, VertexHandle};

use crate::types::BlendSeamCornerDeclineReason as Reason;

use super::super::blend_seam_corner_repair::{commit_edge_authoring, neighbour_chart_point};
use super::super::blend_seam_corner_swept::{project_seeded, surface_distance};
use super::super::surface_pair_march as marched;

/// A reflex chain end that lies on one strip's wall rail but off the
/// other strip, and the point on that rail where it meets the other
/// strip.
pub(super) struct WallEnd {
    pub(super) vertex: VertexHandle,
    pub(super) target: Point3,
}

/// Bisection steps for the rail root: enough to reach roundoff from a
/// bracket a few thousandths long.
const BISECTIONS: usize = 80;

/// Samples per direction when bracketing the rail root.
const BRACKET_SAMPLES: usize = 64;

fn face_of(brep: &BRep, he: HalfEdgeHandle) -> Option<FaceHandle> {
    Some(brep.loops.get(brep.half_edges.get(he)?.loop_)?.face)
}

fn geometry(brep: &BRep, face: FaceHandle) -> Option<SurfaceGeometry> {
    Some(brep.faces.get(face)?.surface.as_ref()?.geometry.clone())
}

/// `p`'s signed distance to `g`, positive on the side its normal points.
fn signed_distance(g: &SurfaceGeometry, p: Point3) -> Option<f64> {
    let (u, v) = project_seeded(g, &p, None)?;
    let n = g.normal_at(u, v).normalized()?;
    Some((p - g.point_at(u, v)).dot(&n))
}

/// The wall end at `vertex`, when it is one: an incident straight edge
/// between strip `s` (one of `strips`) and a non-strip wall face, with
/// the vertex off the other strip. The slide is bounded to 32 times the
/// vertex's distance from the other strip, so only a vertex the anchor
/// stage put near its true position moves.
pub(super) fn find(
    brep: &BRep,
    vertex: VertexHandle,
    strips: (FaceHandle, FaceHandle),
    pos_tol: f64,
) -> Option<WallEnd> {
    let pos = brep.vertices.get(vertex)?.position;
    for (_, hd) in brep.half_edges.iter().filter(|(_, h)| h.origin == vertex) {
        let (fa, fb) = (
            face_of(brep, brep.edges.get(hd.edge)?.half_edge)?,
            face_of(
                brep,
                brep.half_edges
                    .get(brep.edges.get(hd.edge)?.half_edge)?
                    .twin,
            )?,
        );
        let is_strip = |f: FaceHandle| f == strips.0 || f == strips.1;
        let (strip, wall) = match (is_strip(fa), is_strip(fb)) {
            (true, false) => (fa, fb),
            (false, true) => (fb, fa),
            _ => continue,
        };
        let other = if strip == strips.0 {
            strips.1
        } else {
            strips.0
        };
        let Some(CurveGeometry::Line { direction, .. }) =
            brep.edges.get(hd.edge)?.curve.as_ref().map(|c| &c.geometry)
        else {
            continue;
        };
        let dir = direction.normalized()?;
        let (g_other, g_strip, g_wall) = (
            geometry(brep, other)?,
            geometry(brep, strip)?,
            geometry(brep, wall)?,
        );
        let d0 = surface_distance(&g_other, &pos, None);
        if !(d0.is_finite() && d0 > pos_tol) {
            continue;
        }
        let along = |t: f64| pos + dir * t;
        let g0 = signed_distance(&g_other, pos)?;
        let reach = 32.0 * d0;
        let mut bracket = None;
        'scan: for k in 1..=BRACKET_SAMPLES {
            for sign in [1.0, -1.0] {
                let t = sign * reach * k as f64 / BRACKET_SAMPLES as f64;
                let g = signed_distance(&g_other, along(t))?;
                if g.signum() != g0.signum() {
                    let t_prev = sign * reach * (k - 1) as f64 / BRACKET_SAMPLES as f64;
                    bracket = Some((t_prev, t, g0.signum()));
                    break 'scan;
                }
            }
        }
        let (mut lo, mut hi, lo_sign) = bracket?;
        for _ in 0..BISECTIONS {
            let mid = 0.5 * (lo + hi);
            if signed_distance(&g_other, along(mid))?.signum() == lo_sign {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        let target = along(0.5 * (lo + hi));
        let on = |g: &SurfaceGeometry| surface_distance(g, &target, None) <= pos_tol;
        if on(&g_other) && on(&g_strip) && on(&g_wall) {
            return Some(WallEnd { vertex, target });
        }
    }
    None
}

/// Re-author every edge at the moved `vertex` except `seam`: a straight
/// edge as the line through its new ends (linear pcurves, the plane and
/// cylinder-ruling case the rail is), any other edge as the marched
/// intersection of its two faces.
pub(super) fn reauthor_incident(
    brep: &mut BRep,
    vertex: VertexHandle,
    seam: EdgeHandle,
    pos_tol: f64,
) -> Result<(), Reason> {
    let mut edges: Vec<EdgeHandle> = brep
        .half_edges
        .iter()
        .filter(|(_, h)| h.origin == vertex && h.edge != seam)
        .map(|(_, h)| h.edge)
        .collect();
    edges.sort();
    edges.dedup();
    for edge in edges {
        let primary = brep.edges.get(edge).ok_or(Reason::EulerRejected)?.half_edge;
        let twin = brep
            .half_edges
            .get(primary)
            .ok_or(Reason::EulerRejected)?
            .twin;
        let origin = |h: HalfEdgeHandle| -> Result<Point3, Reason> {
            let o = brep.half_edges.get(h).ok_or(Reason::EulerRejected)?.origin;
            Ok(brep.vertices.get(o).ok_or(Reason::EulerRejected)?.position)
        };
        let (a, b) = (origin(primary)?, origin(twin)?);
        let sides: Vec<(HalfEdgeHandle, FaceHandle, SurfaceGeometry)> = [primary, twin]
            .into_iter()
            .map(|he| {
                let face = face_of(brep, he).ok_or(Reason::EulerRejected)?;
                Ok((he, face, geometry(brep, face).ok_or(Reason::EulerRejected)?))
            })
            .collect::<Result<_, Reason>>()?;
        // A straight edge only where the segment lies on both faces (the
        // rail, a cylinder ruling on its tangent plane); the stored curve
        // type is not evidence: the anchor stage authored the arc strip's
        // short wall edge as a chord.
        let straight = (1..4).all(|k| {
            let p = a + (b - a) * (k as f64 / 4.0);
            sides
                .iter()
                .all(|(_, _, g)| surface_distance(g, &p, None) <= pos_tol)
        });
        let (curve, staged) = if straight {
            let curve = CurveOnEdge::line_segment(a, b, pos_tol).ok_or(Reason::EulerRejected)?;
            let mut staged = Vec::new();
            for (i, (he, face, g)) in sides.iter().enumerate() {
                let (from, to) = if i == 0 { (a, b) } else { (b, a) };
                // Seed each end from the half-edge's current pcurve, so
                // a periodic chart keeps its branch.
                let ends = brep
                    .half_edges
                    .get(*he)
                    .and_then(|h| h.pcurve.as_ref())
                    .map(|pc| {
                        let (t0, t1) = pc.parameter_range();
                        let (p, q) = (pc.point_at(t0), pc.point_at(t1));
                        ((p.x, p.y), (q.x, q.y))
                    });
                let uv_from =
                    project_seeded(g, &from, ends.map(|e| e.0)).ok_or(Reason::EulerRejected)?;
                let uv_to =
                    project_seeded(g, &to, ends.map(|e| e.1)).ok_or(Reason::EulerRejected)?;
                staged.push((*face, *he, Pcurve::linear(uv_from, uv_to, 0.0, 1.0)));
            }
            (curve, staged)
        } else {
            let seam_m =
                marched::march_intersection(&sides[0].2, &sides[1].2, a, b, pos_tol, pos_tol)?;
            let curve = marched::intersection_curve(&seam_m, false)?;
            let mut staged = Vec::new();
            for (i, (he, face, g)) in sides.iter().enumerate() {
                let h = brep.half_edges.get(*he).ok_or(Reason::EulerRejected)?;
                let seed = if i == 0 {
                    neighbour_chart_point(brep, *face, h.prev, true)
                } else {
                    neighbour_chart_point(brep, *face, h.next, false)
                };
                staged.push((
                    *face,
                    *he,
                    marched::intersection_pcurve(g, &seam_m, i == 1, seed)?,
                ));
            }
            (curve, staged)
        };
        commit_edge_authoring(brep, edge, curve, staged);
        if brep.measure_edge_deviation(edge) > pos_tol {
            return Err(Reason::ValidationFailed);
        }
    }
    Ok(())
}

/// `to` moved by whole periods onto the branch the rail runs along from
/// `from`: the travel the strip half-edge's own pcurve makes (`seeds`),
/// else the shorter way. A strip crossing its own chart seam projects
/// the seam end onto either branch (PK-798: a torus root strip split by
/// the housing seam, where the torus seam also lies), and the iso walk
/// between branches ran the long way round the torus, every sample
/// still on the housing's contact circle.
fn same_branch(
    surface: &SurfaceGeometry,
    from: (f64, f64),
    to: (f64, f64),
    seeds: Option<((f64, f64), (f64, f64))>,
) -> (f64, f64) {
    let (pu, pv) = surface.periodicity();
    let travel = seeds.map(|(a, b)| (b.0 - a.0, b.1 - a.1));
    let unwrap = |f: f64, t: f64, period: Option<f64>, want: Option<f64>| match period {
        Some(p) if p > 0.0 => {
            let target = want.unwrap_or(0.0);
            t - p * ((t - f - target) / p).round()
        }
        _ => t,
    };
    (
        unwrap(from.0, to.0, pu, travel.map(|d| d.0)),
        unwrap(from.1, to.1, pv, travel.map(|d| d.1)),
    )
}

/// Re-author each strip-to-support edge at `vertex` that is off its
/// faces as the strip's contact iso-curve between its ends.
///
/// A reflex pair's support rails end at the mitre point, which lies
/// past the root edge's own trim by design (PK-740). A rail between two
/// mitre points (a fin root arc's housing rail, M⁺ to M⁻) was left a
/// straight chord. On a plate that chord is exact; on a cylinder
/// housing it sagged 5.5e-2 (PK-742). A rolling-ball strip touches its
/// support along an iso-line of its own chart, so the rail is that
/// iso-line. Every sample is checked on the support; an edge whose
/// ends do not share an iso-parameter declines.
pub(super) fn reauthor_contact_rails(
    brep: &mut BRep,
    vertex: VertexHandle,
    strips: (FaceHandle, FaceHandle),
    pos_tol: f64,
) -> Result<(), Reason> {
    const RAIL_SAMPLES: usize = 33;
    let mut edges: Vec<EdgeHandle> = brep
        .half_edges
        .iter()
        .filter(|(_, h)| h.origin == vertex)
        .map(|(_, h)| h.edge)
        .collect();
    edges.sort();
    edges.dedup();
    for edge in edges {
        if brep.measure_edge_deviation(edge) <= pos_tol {
            continue;
        }
        let primary = brep.edges.get(edge).ok_or(Reason::EulerRejected)?.half_edge;
        let twin = brep
            .half_edges
            .get(primary)
            .ok_or(Reason::EulerRejected)?
            .twin;
        let (Some(f0), Some(f1)) = (face_of(brep, primary), face_of(brep, twin)) else {
            continue;
        };
        let is_strip = |f: FaceHandle| f == strips.0 || f == strips.1;
        let strip_side = match (is_strip(f0), is_strip(f1)) {
            (true, false) => 0,
            (false, true) => 1,
            _ => continue,
        };
        let sides = [(primary, f0), (twin, f1)];
        let (strip_he, strip_face) = sides[strip_side];
        let g_strip = geometry(brep, strip_face).ok_or(Reason::EulerRejected)?;
        let g_support = geometry(brep, sides[1 - strip_side].1).ok_or(Reason::EulerRejected)?;
        let origin = |h: HalfEdgeHandle| -> Result<Point3, Reason> {
            let o = brep.half_edges.get(h).ok_or(Reason::EulerRejected)?.origin;
            Ok(brep.vertices.get(o).ok_or(Reason::EulerRejected)?.position)
        };
        // The strip side's own ends, seeded from its pcurve.
        let (from, to) = (origin(strip_he)?, origin(brep.half_edges[strip_he].twin)?);
        let seeds = brep
            .half_edges
            .get(strip_he)
            .and_then(|h| h.pcurve.as_ref())
            .map(|pc| {
                let (t0, t1) = pc.parameter_range();
                let (p, q) = (pc.point_at(t0), pc.point_at(t1));
                ((p.x, p.y), (q.x, q.y))
            });
        let uv_from =
            project_seeded(&g_strip, &from, seeds.map(|s| s.0)).ok_or(Reason::EulerRejected)?;
        let uv_to =
            project_seeded(&g_strip, &to, seeds.map(|s| s.1)).ok_or(Reason::EulerRejected)?;
        let uv_to = same_branch(&g_strip, uv_from, uv_to, seeds);
        let iso_tol = 1e-9 * (1.0 + uv_from.0.abs().max(uv_from.1.abs()));
        let along = |t: f64| -> (f64, f64) {
            if (uv_from.1 - uv_to.1).abs() <= iso_tol {
                (uv_from.0 + (uv_to.0 - uv_from.0) * t, uv_from.1)
            } else {
                (uv_from.0, uv_from.1 + (uv_to.1 - uv_from.1) * t)
            }
        };
        let iso = (uv_from.1 - uv_to.1).abs() <= iso_tol || (uv_from.0 - uv_to.0).abs() <= iso_tol;
        if !iso {
            return Err(Reason::MarchUnresolved);
        }
        let mut samples: Vec<Point3> = (0..RAIL_SAMPLES)
            .map(|k| {
                let (u, v) = along(k as f64 / (RAIL_SAMPLES - 1) as f64);
                g_strip.point_at(u, v)
            })
            .collect();
        if samples
            .iter()
            .any(|p| surface_distance(&g_support, p, None) > pos_tol)
        {
            return Err(Reason::MarchUnresolved);
        }
        // Samples run along the canonical half-edge.
        if strip_side == 1 {
            samples.reverse();
        }
        let mut params = vec![0.0];
        for w in samples.windows(2) {
            params.push(params[params.len() - 1] + w[0].distance_to(&w[1]));
        }
        let total = params[params.len() - 1].max(1e-300);
        params.iter_mut().for_each(|p| *p /= total);
        let rail = marched::MarchedIntersection {
            samples,
            params,
            band: pos_tol,
        };
        let curve = marched::intersection_curve(&rail, false)?;
        let mut staged = Vec::new();
        for (i, (he, face)) in sides.iter().enumerate() {
            let g = geometry(brep, *face).ok_or(Reason::EulerRejected)?;
            let h = brep.half_edges.get(*he).ok_or(Reason::EulerRejected)?;
            let seed = if i == 0 {
                neighbour_chart_point(brep, *face, h.prev, true)
            } else {
                neighbour_chart_point(brep, *face, h.next, false)
            };
            staged.push((
                *face,
                *he,
                marched::intersection_pcurve(&g, &rail, i == 1, seed)?,
            ));
        }
        commit_edge_authoring(brep, edge, curve, staged);
        if brep.measure_edge_deviation(edge) > pos_tol {
            return Err(Reason::ValidationFailed);
        }
    }
    Ok(())
}
