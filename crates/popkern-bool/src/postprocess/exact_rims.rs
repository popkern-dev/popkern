//! Exact rims (PK-417): intersection edges that trace a circle leave
//! the boolean as the analytic `Arc` they are.
//!
//! The pipeline carries every intersection curve as a NURBS triple —
//! including the rational circles the plane × cylinder / cone solvers
//! author exactly — so a bore's rims reached the result as NURBS:
//! `circular()` could not select them, fillets and STEP consumers lost
//! the analytic path, and two rims of one bore differed in spelling
//! from the extruded cylinder's own caps. This pass runs at the tail
//! of output post-processing: every NURBS edge whose samples lie on
//! one circle within the body's tolerance is re-spelled as an
//! angle-parameterised `Arc` with pcurves co-parameterised on both
//! faces (the exact plane arm where it applies, the chart fit
//! elsewhere), verified sample for sample before the edge is
//! rewritten. An edge that fails any check keeps its NURBS.
//!
//! Scope: rims and seam pieces on analytic faces (plane, cylinder,
//! cone, sphere, torus) of every boolean result — bores, countersinks
//! and cut slots (subtraction), boss junctions (union), sphere and
//! torus caps. Sphere and torus rims and seams joined the scope with
//! PK-418 once the mesher put every hole ring of a periodic chart on
//! the outer ring's period sheet and the sphere/torus iso-circle
//! pcurve arms made those rims exact by construction; union and
//! intersection rims joined with PK-419 once the fillet's spine
//! walker seeded a closed rim at its own station (a midpoint seed had
//! converged to the antipodal ball, which an exact arc exposed and a
//! NURBS rim only escaped by ULP luck).

use popkern_core::geometry::Point3;
use popkern_core::geometry::circle_recognition::{arc_between, arc_on_edge, recognize_circular};
use popkern_core::geometry::curve::CurveGeometry;
use popkern_core::geometry::curve::{CurveEval as _, CurveOnEdge};
use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::edge_pcurve_reauthor::coparameterised_pcurves;
use popkern_core::topology::{BRep, EdgeHandle};

/// What the pass did.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct ExactRimReport {
    /// NURBS edges re-spelled as arcs.
    pub promoted: usize,
    /// Circular edges left as NURBS because a pcurve could not be
    /// authored or verified on one of their faces.
    pub declined: usize,
}

/// Promote every circular NURBS edge of `brep` to an exact `Arc`.
pub fn promote_exact_rims(brep: &mut BRep, tolerance: &Tolerance) -> ExactRimReport {
    let gate = tolerance.linear();
    let mut report = ExactRimReport::default();
    let edges: Vec<EdgeHandle> = brep.edges.keys().collect();
    for edge in edges {
        // Only NURBS edges are candidates, and a two-point degree-1
        // span is a straight segment the recognizer's sagitta test
        // rejects anyway — skipped before sampling (and before cloning
        // every edge's curve: PK-804).
        let Some(CurveGeometry::Nurbs(n)) = brep
            .edges
            .get(edge)
            .and_then(|e| e.curve.as_ref())
            .map(|c| &c.geometry)
        else {
            continue;
        };
        if n.degree == 1
            && n.control_points.len() == 2
            && n.control_points[0] != n.control_points[1]
        {
            continue;
        }
        let Some(curve) = brep.edges.get(edge).and_then(|e| e.curve.clone()) else {
            continue;
        };
        if !rim_between_analytic_faces(brep, edge) {
            continue;
        }
        let range = (curve.t_start, curve.t_end);
        let Some(fit) = recognize_circular(&curve.geometry, range, gate) else {
            continue;
        };
        // The arc runs the old curve's way, from the vertex at its
        // start to the vertex at its end: each half-edge keeps the
        // sense it had along its curve.
        let arc = match vertex_ends(brep, edge, &curve, gate) {
            Some((start, end)) => arc_between(&fit, fit.normal, start, end),
            None => arc_on_edge(&curve.geometry, range, &fit),
        };
        match coparameterised_pcurves(brep, edge, &arc, gate * 10.0, tolerance) {
            Some(pcurves) => {
                if popkern_core::trace::trace_enabled("exact-rims") {
                    for (he, pcurve) in &pcurves {
                        let before = brep.half_edges.get(*he).and_then(|h| h.pcurve.as_ref());
                        let show = |p: &Pcurve| {
                            let a = p.point_at(p.t_start);
                            let b = p.point_at(p.t_end);
                            let m = p.point_at(p.t_start + (p.t_end - p.t_start) * 0.25);
                            format!(
                                "t=[{:.3},{:.3}] uv ({:.3},{:.3})->({:.3},{:.3}) q=({:.3},{:.3})",
                                p.t_start, p.t_end, a.x, a.y, b.x, b.y, m.x, m.y
                            )
                        };
                        eprintln!(
                            "[exact-rims] edge {edge:?} he {he:?}: before {} | after {}",
                            before.map(show).unwrap_or_else(|| "none".into()),
                            show(pcurve)
                        );
                    }
                }
                if let Some(e) = brep.edges.get_mut(edge) {
                    e.curve = Some(arc);
                }
                for (he, pcurve) in pcurves {
                    if let Some(h) = brep.half_edges.get_mut(he) {
                        h.pcurve = Some(pcurve);
                    }
                }
                report.promoted += 1;
            }
            None => report.declined += 1,
        }
    }
    report
}

/// The edge's vertex positions at the curve's start and end, when
/// each curve end lies within `gate` of one of the edge's two
/// vertices: the arc is anchored at them (PK-944). The forward
/// half-edge does not always run its curve forward (a graph-merge
/// rim's does not), so the pairing is by distance, not by sense. The vertices are what the neighbouring edges share; a
/// sampled curve's end can sit 2e-9 along the circle from its vertex
/// (a sphere's rim split at its seam), and an arc started there gave
/// the rim's iso-circle pcurve a u the seam meridian's pcurve did not
/// reach, so the face's parameter-space boundary missed closing by
/// that much and its area was 3.4e-10 off. A vertex farther than the
/// gate is not this edge's end in any sense the arc can honour, and
/// the curve's own ends stand.
fn vertex_ends(
    brep: &BRep,
    edge: EdgeHandle,
    curve: &CurveOnEdge,
    gate: f64,
) -> Option<(Point3, Point3)> {
    let forward = brep.half_edges.get(brep.edges.get(edge)?.half_edge)?;
    let backward = brep.half_edges.get(forward.twin)?;
    let a = brep.vertices.get(forward.origin)?.position;
    let b = brep.vertices.get(backward.origin)?.position;
    let (first, last) = (
        curve.geometry.point_at(curve.t_start),
        curve.geometry.point_at(curve.t_end),
    );
    let near = |v: Point3, p: Point3| v.distance_to(&p) <= gate;
    if near(a, first) && near(b, last) {
        Some((a, b))
    } else if near(b, first) && near(a, last) {
        Some((b, a))
    } else {
        None
    }
}

/// The promotion's scope: an edge whose faces are analytic surfaces of
/// the bore, boss, countersink and cap class — planes, cylinders,
/// cones, spheres and tori — including a periodic seam (both
/// half-edges on one face): a boolean splits a sphere's or torus's
/// seam meridian into pieces it re-spells as NURBS, and those pieces
/// are the same exact circle as the primitive's own seam. Free-form
/// (NURBS, Gregory, polyhedral-spline) faces keep their NURBS rims:
/// their pcurves are fits either way and the arc spelling buys the
/// consumer nothing it can select on.
fn rim_between_analytic_faces(brep: &BRep, edge: EdgeHandle) -> bool {
    let Some(e) = brep.edges.get(edge) else {
        return false;
    };
    let Some(forward) = brep.half_edges.get(e.half_edge) else {
        return false;
    };
    let Some(backward) = brep.half_edges.get(forward.twin) else {
        return false;
    };
    let face_of = |loop_| brep.loops.get(loop_).map(|l| l.face);
    let (Some(fa), Some(fb)) = (face_of(forward.loop_), face_of(backward.loop_)) else {
        return false;
    };
    [fa, fb].iter().all(|f| {
        matches!(
            brep.faces
                .get(*f)
                .and_then(|f| f.surface.as_ref())
                .map(|s| &s.geometry),
            Some(SurfaceGeometry::Plane { .. })
                | Some(SurfaceGeometry::Cylinder { .. })
                | Some(SurfaceGeometry::Cone { .. })
                | Some(SurfaceGeometry::Sphere { .. })
                | Some(SurfaceGeometry::Torus { .. })
        )
    })
}
