//! Exact rims for a rolling-ball corner carved by the setback rung
//! (PK-864). The carve leaves the corner region bounded by the coarse
//! cubic arcs the topology splits (3.9e-3 off the ball on a rounded
//! tetrahedron, r 0.7); the KP20 commit replaces them with its cap rows.
//! A ball cap needs the exact curves instead: each boundary edge lies on
//! its round's cross-section through the ball's centre — the circle of
//! radius r about the centre in the plane normal to the round's axis,
//! which lies on both the round's cylinder and the ball.

use popkern_core::geometry::curve::{CurveEval, CurveGeometry, CurveOnEdge};
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::topology::{BRep, FaceHandle};

use super::{face_of, loop_half_edges};

/// The arc of the circle about `center` (radius `radius`, plane normal
/// `normal`) from `start` to `end` the short way, as an edge curve.
fn short_arc(
    center: Point3,
    normal: Vec3,
    radius: f64,
    start: Point3,
    end: Point3,
) -> Option<CurveOnEdge> {
    let tau = std::f64::consts::TAU;
    for n in [normal, -normal] {
        let geometry = CurveGeometry::Arc {
            center,
            normal: n,
            radius,
            start_angle: 0.0,
            end_angle: tau,
        };
        let (t0, t1) = (
            geometry.closest_point(&start)?,
            geometry.closest_point(&end)?,
        );
        let sweep = (t1 - t0).rem_euclid(tau);
        if sweep > 0.0 && sweep <= std::f64::consts::PI {
            // The arc's own angle range is the edge's window: a window
            // past 2π on a [0, 2π] circle is clipped by its sub-range
            // (the corner's pcurve sat 8.5e-3 off one rim per corner).
            return Some(CurveOnEdge {
                geometry: CurveGeometry::Arc {
                    center,
                    normal: n,
                    radius,
                    start_angle: t0,
                    end_angle: t0 + sweep,
                },
                t_start: t0,
                t_end: t0 + sweep,
            });
        }
    }
    None
}

/// Give every boundary edge of `face` (a carved ball corner) the exact
/// cross-section arc of the round across it, after checking its ends sit
/// on the ball within `band`. Returns the faces whose pcurves the new
/// curves stale (the corner and its rounds), or why it cannot.
pub(crate) fn exact_ball_rims(
    brep: &mut BRep,
    face: FaceHandle,
    center: Point3,
    radius: f64,
    band: f64,
) -> Result<Vec<FaceHandle>, String> {
    let mut faces = vec![face];
    for he in loop_half_edges(brep, face) {
        let h = brep.half_edges.get(he).ok_or("stale half-edge")?;
        let (edge, twin) = (h.edge, h.twin);
        let round = face_of(brep, twin).ok_or("stale twin")?;
        let Some(SurfaceGeometry::Cylinder { axis, .. }) = brep
            .faces
            .get(round)
            .and_then(|f| f.surface.as_ref())
            .map(|s| &s.geometry)
        else {
            return Err("a corner boundary edge borders a face that is not a round".to_string());
        };
        let normal = axis.normalized().ok_or("a round without an axis")?;
        let e = brep.edges.get(edge).ok_or("stale edge")?;
        let first = brep.half_edges.get(e.half_edge).ok_or("stale half-edge")?;
        let start = brep
            .vertices
            .get(first.origin)
            .ok_or("stale vertex")?
            .position;
        let end = brep
            .half_edges
            .get(first.twin)
            .and_then(|t| brep.vertices.get(t.origin))
            .ok_or("stale vertex")?
            .position;
        for p in [start, end] {
            let off = (p.distance_to(&center) - radius).abs();
            if off > band {
                return Err(format!(
                    "a corner vertex leaves the ball by {off:.3e} (band {band:.1e})"
                ));
            }
        }
        let curve = short_arc(center, normal, radius, start, end)
            .ok_or("no short arc between the rim's ends")?;
        if let Some(e) = brep.edges.get_mut(edge) {
            e.curve = Some(curve);
        }
        if !faces.contains(&round) {
            faces.push(round);
        }
    }
    for &f in &faces {
        for he in loop_half_edges(brep, f) {
            if let Some(h) = brep.half_edges.get_mut(he) {
                h.pcurve = None;
            }
        }
        brep.face_use_pcurves.retain(|key, _| key.face != f);
    }
    Ok(faces)
}
