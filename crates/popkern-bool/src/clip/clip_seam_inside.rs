//! Containment of a point ON a face's own seam line (PK-931).
//!
//! A periodic face's seam pair is a chart cut, not a boundary: a point
//! on it is interior to the face as soon as the face covers either side
//! of the cut. A uv polygon test there is a coin flip between the two
//! branch copies, so a record lying exactly on the seam read inside or
//! outside by rotation (measured: a cone minus a box whose face plane
//! contains the cone's seam generator — kept unrotated, dropped about z
//! by 0.7, and the assembly stranded four one-sided record edges). The
//! test steps a hair off the seam to each side.

use popkern_core::geometry::Point3;
use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};

/// Offset off the seam, and the band that counts as on it, as fractions
/// of the period.
const SEAM_STEP: f64 = 1e-4;
const SEAM_BAND: f64 = 1e-6;

/// [`crate::split_inner_loop::point_inside_face_uv`], deciding a point
/// on the face's own seam line by the face's cover on either side.
pub(super) fn point_inside_face_uv_seam_aware(
    brep: &BRep,
    face: FaceHandle,
    q: Point3,
    tolerance: &Tolerance,
) -> bool {
    let inside =
        |p: Point3| crate::split_inner_loop::point_inside_face_uv(brep, face, p, tolerance);
    if inside(q) {
        return true;
    }
    let Some(period) = brep
        .faces
        .get(face)
        .and_then(|f| f.surface.as_ref())
        .and_then(|s| s.geometry.parameter_periodicity().0)
        .filter(|p| *p > 0.0)
    else {
        return false;
    };
    let on_seam = seam_us(brep, face).iter().any(|&s| {
        let d = q.x - s;
        (d - period * (d / period).round()).abs() < SEAM_BAND * period
    });
    on_seam
        && [1.0, -1.0]
            .iter()
            .any(|side| inside(Point3::new(q.x + side * SEAM_STEP * period, q.y, q.z)))
}

/// The u of every seam-pair half-edge of `face` (both halves on the
/// face, constant u along the pcurve).
fn seam_us(brep: &BRep, face: FaceHandle) -> Vec<f64> {
    let Some(f) = brep.faces.get(face) else {
        return Vec::new();
    };
    let mut out = Vec::new();
    for l in std::iter::once(f.outer_loop).chain(f.inner_loops.iter().copied()) {
        let Some(start) = brep.loops.get(l).map(|l| l.half_edge) else {
            continue;
        };
        let mut h = start;
        for _ in 0..=brep.half_edges.len() {
            let Some(he) = brep.half_edges.get(h) else {
                break;
            };
            let twin_here = brep
                .half_edges
                .get(he.twin)
                .and_then(|t| brep.loops.get(t.loop_))
                .is_some_and(|tl| tl.face == face);
            if twin_here && let Some(pc) = he.pcurve.as_ref() {
                let (t0, t1) = pc.parameter_range();
                let (a, b) = (pc.point_at(t0), pc.point_at(t1));
                if (a.x - b.x).abs() < 1e-9 && (a.y - b.y).abs() > 0.0 {
                    out.push(a.x);
                }
            }
            h = he.next;
            if h == start {
                break;
            }
        }
    }
    out
}
