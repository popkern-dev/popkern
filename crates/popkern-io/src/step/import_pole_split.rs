//! Edges through sphere poles (PK-646).
//!
//! A sphere's chart collapses its poles to points, so a boundary curve
//! passing THROUGH a pole (a great-circle arc over the top of a ball,
//! a dome cut by a plane through the pole) runs up one meridian and
//! down the opposite one: its chart image jumps half a period in u at
//! the pole, and no continuous pcurve represents it. The projected
//! pcurve smears the jump across the chart (measured: a radius-10 dome
//! through the pole carried a 0.251 pcurve deviation; jack body 29's
//! corner sphere 1.07e-2), the winding audit read the smeared loop's
//! wrong side and "reconciled" an honest flag, and the mesher charted
//! the wrong region.
//!
//! The repair splits such an edge AT the pole, by Euler split: each
//! half then ends at the pole, where its pcurve's free azimuth snaps to
//! its own meridian, and the chart discontinuity sits at a vertex — the
//! pole-touching loop every chart consumer already handles (the
//! winding audit defers it to its tangency tier; the mesher runs the
//! pole row between the halves). The halves' pcurves are re-projected
//! on both faces. The repair runs on a copy and commits only when the
//! body validates; every repair and decline is logged.

use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle};

use crate::step::diagnostics::{DiagnosticLog, Phase};

/// Split every edge that passes through a pole of an adjacent sphere
/// face's chart at that pole; returns the number of splits.
pub(crate) fn split_edges_at_poles(brep: &mut BRep, log: &mut DiagnosticLog) -> usize {
    let mut splits = 0;
    // Each split changes the edge set; rescan until none is left (an
    // edge over both poles splits twice).
    for _ in 0..2 * brep.edges.len() + 2 {
        let Some((edge, pole, faces)) = next_pole_crossing(brep) else {
            break;
        };
        let mut candidate = brep.clone();
        match split_at_pole(&mut candidate, edge, pole, &faces) {
            Ok(()) => match candidate.validate() {
                Ok(()) => {
                    *brep = candidate;
                    splits += 1;
                    continue;
                }
                Err(e) => log.warn(
                    Phase::MapTopology,
                    None,
                    format!("{edge:?} through a sphere pole not split: the split body failed validation: {e}"),
                ),
            },
            Err(reason) => log.warn(
                Phase::MapTopology,
                None,
                format!("{edge:?} through a sphere pole not split: {reason}"),
            ),
        }
        break;
    }
    if splits > 0 {
        log.warn(
            Phase::MapTopology,
            None,
            format!("split {splits} edge(s) at the sphere pole they pass through"),
        );
    }
    splits
}

/// The first edge passing through (not ending at) a pole of an
/// adjacent sphere face, with the pole and the edge's two faces.
fn next_pole_crossing(brep: &BRep) -> Option<(EdgeHandle, Point3, [FaceHandle; 2])> {
    let gate = 10.0 * brep.tolerance.linear();
    for (eh, edge) in brep.edges.iter() {
        let Some(curve) = edge.curve.as_ref() else {
            continue;
        };
        let he = brep.half_edges.get(edge.half_edge)?;
        let twin = brep.half_edges.get(he.twin)?;
        let faces = [
            brep.loops.get(he.loop_)?.face,
            brep.loops.get(twin.loop_)?.face,
        ];
        let ends = [
            brep.vertices.get(he.origin)?.position,
            brep.vertices.get(twin.origin)?.position,
        ];
        for face in faces {
            let Some(SurfaceGeometry::Sphere { .. }) =
                brep.faces.get(face)?.surface.as_ref().map(|s| &s.geometry)
            else {
                continue;
            };
            let sphere = &brep.faces.get(face)?.surface.as_ref()?.geometry;
            for v_pole in [0.0, std::f64::consts::PI] {
                let pole = sphere.point_at(0.0, v_pole);
                if ends.iter().any(|p| p.distance_to(&pole) <= gate) {
                    continue; // the pole is already a vertex
                }
                let (t, distance) = closest_parameter(curve, &pole);
                let span = curve.t_end - curve.t_start;
                let fraction = (t - curve.t_start) / span;
                if distance <= gate && fraction > 1e-6 && fraction < 1.0 - 1e-6 {
                    return Some((eh, curve.geometry.point_at(t), faces));
                }
            }
        }
    }
    None
}

/// The edge-curve parameter nearest `target` within the trim, and its
/// distance: a 256-sample scan refined by golden section.
fn closest_parameter(
    curve: &popkern_core::geometry::curve::CurveOnEdge,
    target: &Point3,
) -> (f64, f64) {
    let (a, b) = (curve.t_start, curve.t_end);
    let d = |t: f64| curve.geometry.point_at(t).distance_to(target);
    let n = 256;
    let k = (0..=n)
        .min_by(|&i, &j| {
            let (ti, tj) = (
                a + (b - a) * i as f64 / n as f64,
                a + (b - a) * j as f64 / n as f64,
            );
            d(ti).total_cmp(&d(tj))
        })
        .unwrap_or(0);
    let step = (b - a) / n as f64;
    let at = a + step * k as f64;
    let (mut lo, mut hi) = ((at - step).max(a.min(b)), (at + step).min(a.max(b)));
    for _ in 0..80 {
        let m1 = lo + (hi - lo) * 0.382;
        let m2 = lo + (hi - lo) * 0.618;
        if d(m1) < d(m2) {
            hi = m2;
        } else {
            lo = m1;
        }
    }
    let t = 0.5 * (lo + hi);
    (t, d(t))
}

fn split_at_pole(
    brep: &mut BRep,
    edge: EdgeHandle,
    at: Point3,
    faces: &[FaceHandle; 2],
) -> Result<(), String> {
    // Split by the curve fraction (re-inverting the point returns the
    // principal angle, outside a trim past one turn). `split_edge`
    // reads the fraction along the edge's own half-edge, whose
    // direction a closed curve leaves to a pcurve witness — and these
    // pcurves are the smeared ones — so the split vertex is checked
    // against the pole, mirrored when it missed.
    let gate = 10.0 * brep.tolerance.linear();
    let curve = brep.edges[edge].curve.clone().ok_or("edge has no curve")?;
    let (t, _) = closest_parameter(&curve, &at);
    let fraction = (t - curve.t_start) / (curve.t_end - curve.t_start);
    let original = brep.clone();
    let mut split = None;
    for f in [fraction, 1.0 - fraction] {
        *brep = original.clone();
        let s = brep
            .split_edge(edge, f)
            .map_err(|e| format!("split: {e}"))?;
        if brep.vertices[s.vertex].position.distance_to(&at) <= gate {
            split = Some(s);
            break;
        }
    }
    let split = split.ok_or("the split vertex missed the pole both ways")?;
    // Both halves' pcurves were cut from the smeared projection;
    // re-project them on both faces.
    for half in [split.edge_a, split.edge_b] {
        let first = brep.edges[half].half_edge;
        for he in [first, brep.half_edges[first].twin] {
            brep.half_edges[he].pcurve = None;
        }
    }
    for &face in faces {
        brep.recompute_face_pcurves(face)
            .map_err(|e| format!("pcurves: {e}"))?;
    }
    Ok(())
}
