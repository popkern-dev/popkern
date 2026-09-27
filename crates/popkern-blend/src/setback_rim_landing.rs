//! Setback rim landing (PK-584): how well each cap-boundary rim edge
//! lies ON its neighbour surface, as the served body records it.
//!
//! A rim that lies on the surface has its own tangent in the surface's
//! tangent plane. The cap normal at the rim is (rim tangent) × (cross
//! jet), and the rim reconciles only correct the cross jet, so a rim
//! tangent tilted off the neighbour is a G1 floor no later pass can
//! lower. That floor ruled variable-radius saddles until the
//! junction-continuous band refit: the rim tangent sat 1e-6..1.2e-5
//! off the strip, while the constant-radius rims read 1e-8. This record
//! makes the tilt a structured, per-edge fact on the result.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::uv_projection::project_point_to_uv;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle};

/// Samples per rim edge.
const SAMPLES: usize = 65;

/// One cap-boundary rim edge's landing on its neighbour face.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SetbackRimLanding {
    /// The rim edge (cap piece ↔ non-cap face).
    pub edge: EdgeHandle,
    /// The non-cap face it lands on.
    pub neighbour: FaceHandle,
    /// Worst sine between the edge curve's tangent and the neighbour's
    /// tangent plane (0 for a rim lying on the neighbour).
    pub worst_tangent_sine: f64,
    /// Edge-curve fraction where [`Self::worst_tangent_sine`] occurs.
    pub worst_fraction: f64,
    /// Worst distance from the edge curve to the neighbour surface.
    pub worst_distance: f64,
}

/// Landing records for every edge between a face in `cap_faces` and a
/// non-cap face. Edges whose curve or neighbour surface cannot be
/// evaluated are skipped.
pub fn setback_rim_landing(brep: &BRep, cap_faces: &[FaceHandle]) -> Vec<SetbackRimLanding> {
    let caps: std::collections::BTreeSet<FaceHandle> = cap_faces.iter().copied().collect();
    let mut out = Vec::new();
    for (eh, edge) in brep.edges.iter() {
        let Some(he) = brep.half_edges.get(edge.half_edge) else {
            continue;
        };
        let Some(tw) = brep.half_edges.get(he.twin) else {
            continue;
        };
        let (Some(la), Some(lb)) = (brep.loops.get(he.loop_), brep.loops.get(tw.loop_)) else {
            continue;
        };
        let neighbour = match (caps.contains(&la.face), caps.contains(&lb.face)) {
            (true, false) => lb.face,
            (false, true) => la.face,
            _ => continue,
        };
        let (Some(curve), Some(surface)) = (
            edge.curve.as_ref(),
            brep.faces
                .get(neighbour)
                .and_then(|f| f.surface.as_ref())
                .map(|s| &s.geometry),
        ) else {
            continue;
        };
        let eval: &dyn CurveEval = &curve.geometry;
        let (t0, t1) = (curve.t_start, curve.t_end);
        if !(t0.is_finite() && t1.is_finite()) || t0 == t1 {
            continue;
        }
        let mut record = SetbackRimLanding {
            edge: eh,
            neighbour,
            worst_tangent_sine: 0.0,
            worst_fraction: 0.0,
            worst_distance: 0.0,
        };
        let mut measured = false;
        for s in 0..SAMPLES {
            let f = s as f64 / (SAMPLES - 1) as f64;
            let t = t0 + (t1 - t0) * f;
            let p = eval.point_at(t);
            let Some(tangent) = eval.derivative_at(t, 1).normalized() else {
                continue;
            };
            let Some((q, n)) = foot_and_normal(surface, &p) else {
                continue;
            };
            measured = true;
            let sine = tangent.dot(&n).abs();
            if sine > record.worst_tangent_sine {
                record.worst_tangent_sine = sine;
                record.worst_fraction = f;
            }
            record.worst_distance = record.worst_distance.max(q.distance_to(&p));
        }
        if measured {
            out.push(record);
        }
    }
    out
}

/// Foot point and unit normal of `surface` at the projection of `p`
/// (central differences; `None` at a degenerate chart point).
fn foot_and_normal(surface: &SurfaceGeometry, p: &Point3) -> Option<(Point3, Vec3)> {
    let (u, v) = project_point_to_uv(surface, p)?;
    let h = 1e-6;
    let du = surface.point_at(u + h, v) - surface.point_at(u - h, v);
    let dv = surface.point_at(u, v + h) - surface.point_at(u, v - h);
    let (lu, lv) = (du.length(), dv.length());
    if lu.min(lv) < 1e-3 * lu.max(lv) {
        return None;
    }
    Some((surface.point_at(u, v), du.cross(&dv).normalized()?))
}
