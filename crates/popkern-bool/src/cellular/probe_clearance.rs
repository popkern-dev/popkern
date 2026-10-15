//! How far an offset-pair probe may reach from a fragment sample and
//! still read the sample's own cell: the distance, along the probe axis,
//! to the first face of the other body.
//!
//! The classifier probes a fragment at a fine offset and, to resolve
//! grazing, at coarser ones. A probe point past the other body's nearest
//! face along the axis classifies a different cell of space, so its
//! verdict says nothing about the fragment (an engraved glyph's
//! top-face fragment sat inside a cutter only 0.2 above and 0.5 below
//! it; the 0.89 probe read Outside on both sides and kept the fragment).

use popkern_core::geometry::curve::CurveGeometry;
use popkern_core::geometry::intersection::surface_curve::{
    SurfaceCurveIntersection, intersect_surface_curve_dense,
};
use popkern_core::geometry::intersection::surface_subdivision::Aabb;
use popkern_core::geometry::uv_projection::bounded_parameter_range;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};
use popkern_core::validation::FaceUvContainment;

use crate::error::PointClassification;

/// An offset pair's verdicts: the `+n` probe's, then the `-n` probe's.
pub(super) type PairVerdict = (PointClassification, PointClassification);

/// The verdict of a fragment's offset-pair probes, `pair_at(e)` probing
/// at offset `e` along the outward normal.
///
/// Probe at the tuned offset `eps` AND at the coarse one: 1 µm probes
/// graze corner-adjacent faces when the sample sits on the other body's
/// edge or corner fold, and the ray parity inverts the verdict
/// (measured: a flush cutter's flank fragment read anti-oriented at its
/// corner sample and the relief-mouth membrane was kept). When the two
/// scales disagree, a third pair between them breaks the tie by
/// majority: blindly trusting the coarser scale mis-served
/// coarse-tolerance bodies whose sample clearance to the other body
/// sits between the two offsets — at linear 1e-3 the 16× probe reached
/// 0.09, PAST the sphere face 0.053 away, and the straddling verdict
/// dropped every cube annulus (measured: centered-sphere subtraction,
/// six half-claimed rim edges). All-distinct keeps the coarse verdict.
///
/// The ladder is for grazing, which needs the sample `near_other_edge`.
/// Clear of the other body's edges, an agreeing fine pair already reads
/// the fragment's own neighbourhood. And a larger probe counts only
/// while it stays in the sample's own cell, nearer than half the
/// `clearance` to the other body's first face along the axis: past it
/// it reads a different region. Without that bound, a 0.56 mm probe
/// found an engraved ring's top-face fragment inside a 7 mm-thick
/// cutter, the 8.9 mm probe overshot the cutter on both sides, the
/// 2.2 mm tie-break straddled its top, and the coarse "both Outside"
/// kept the fragment (the cut declined on 121 curveless edges; at three
/// times the size the coarse probe stayed inside and it served).
pub(super) fn ladder_verdict(
    pair_at: &dyn Fn(f64) -> PairVerdict,
    eps: f64,
    coarse_eps: f64,
    near_other_edge: bool,
    clearance: impl FnOnce() -> f64,
) -> PairVerdict {
    let fine = pair_at(eps);
    let agrees = fine.0 == fine.1
        && matches!(
            fine.0,
            PointClassification::Inside | PointClassification::Outside
        );
    if agrees && !near_other_edge {
        return fine;
    }
    let reach = clearance();
    let local = |e: f64| e < 0.5 * reach;
    let mid_eps = coarse_eps.sqrt() * eps.sqrt();
    if !local(mid_eps) {
        // No larger probe stays in the cell.
        return fine;
    }
    if !local(coarse_eps) {
        // Only the tie-break scale stays in the cell: it is the
        // larger local reading.
        return pair_at(mid_eps);
    }
    let coarse = pair_at(coarse_eps);
    if fine == coarse {
        return fine;
    }
    if pair_at(mid_eps) == fine {
        fine
    } else {
        coarse
    }
}

/// The distance from `sample` along `+axis` or `-axis` (unit) to the
/// nearest point of `brep`'s faces, searched out to `reach`; `reach`
/// when nothing is nearer. Crossings within the linear tolerance of the
/// sample are the sample's own contact, not a wall ahead of the probe,
/// and are skipped.
pub(super) fn axis_clearance(
    brep: &BRep,
    face_aabbs: &[(FaceHandle, Aabb)],
    sample: &Point3,
    axis: &Vec3,
    reach: f64,
    tolerance: &Tolerance,
) -> f64 {
    let mut nearest = reach;
    for sign in [1.0, -1.0] {
        let direction = *axis * sign;
        let end = *sample + direction * reach;
        let (lo, hi) = (
            Point3::new(
                sample.x.min(end.x),
                sample.y.min(end.y),
                sample.z.min(end.z),
            ),
            Point3::new(
                sample.x.max(end.x),
                sample.y.max(end.y),
                sample.z.max(end.z),
            ),
        );
        let line = CurveGeometry::Line {
            origin: *sample,
            direction,
        };
        for (face, aabb) in face_aabbs {
            let pad = tolerance.linear();
            if aabb.max.x < lo.x - pad
                || aabb.min.x > hi.x + pad
                || aabb.max.y < lo.y - pad
                || aabb.min.y > hi.y + pad
                || aabb.max.z < lo.z - pad
                || aabb.min.z > hi.z + pad
            {
                continue;
            }
            let Some(surface) = brep.faces.get(*face).and_then(|f| f.surface.as_ref()) else {
                continue;
            };
            let positions = brep.face_outer_loop_positions(*face);
            if positions.is_empty() {
                continue;
            }
            let (u_range, v_range) = bounded_parameter_range(&positions, &surface.geometry);
            let Ok(hits) = intersect_surface_curve_dense(
                &surface.geometry,
                u_range,
                v_range,
                &line,
                (0.0, nearest),
                tolerance,
            ) else {
                continue;
            };
            for hit in &hits.intersections {
                let SurfaceCurveIntersection::Point { t, u, v, .. } = hit else {
                    continue;
                };
                if *t <= tolerance.linear() || *t >= nearest {
                    continue;
                }
                if brep.face_uv_trim_containment(*face, &surface.geometry, *u, *v)
                    != Some(FaceUvContainment::Outside)
                {
                    nearest = *t;
                }
            }
        }
    }
    nearest
}
