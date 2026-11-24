//! A fragment sample beside the other body's skin but clear of it: its
//! side read off the nearest face, not by straddling probes.
//!
//! The offset-pair classifier probes a fragment's sample at ± a fixed
//! offset (a thousandth of the model, at least a thousand tolerances)
//! and reads straddling verdicts as a coincident contact. A sample
//! nearer the other body than that offset but past the linear reach is
//! not coincident: the material between them is real (a lens two
//! tolerances deep where a cylinder dips into a rounded bracket, YL-E
//! 1142), yet both probes straddled it, the lens fragments read
//! `OnBoundary` anti-oriented and the common served empty. The nearest
//! face of the other body, its foot inside its trim, decides: the
//! sample lies inside where it is behind that face's outward normal.

use popkern_core::geometry::Point3;
use popkern_core::geometry::intersection::surface_subdivision::Aabb;
use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};
use popkern_core::validation::FaceUvContainment;

use crate::error::PointClassification;

/// The side of `other` that `sample` lies on, read off the nearest face
/// whose foot is inside its trim, when the sample is past the linear
/// reach from `other` and nearer than `probe`; `None` otherwise (the
/// offset-pair classifier decides).
pub(super) fn clear_side(
    other: &BRep,
    aabbs: &[(FaceHandle, Aabb)],
    sample: &Point3,
    probe: f64,
    tolerance: &Tolerance,
) -> Option<PointClassification> {
    let near = |b: &Aabb| {
        [
            (sample.x, b.min.x, b.max.x),
            (sample.y, b.min.y, b.max.y),
            (sample.z, b.min.z, b.max.z),
        ]
        .iter()
        .all(|&(c, lo, hi)| c >= lo - probe && c <= hi + probe)
    };
    let mut best: Option<(f64, FaceHandle, (f64, f64))> = None;
    for (face, _) in aabbs.iter().filter(|(_, b)| near(b)) {
        let Some(place) = other.face_point_trim_containment_within(*face, sample, 1e-12) else {
            continue;
        };
        if place.containment != FaceUvContainment::Interior {
            continue;
        }
        if best.is_none_or(|(d, _, _)| place.distance < d) {
            best = Some((place.distance, *face, place.uv));
        }
    }
    let (distance, face, (u, v)) = best?;
    if distance <= tolerance.linear_reach() || distance >= probe {
        return None;
    }
    let sof = other.faces.get(face)?.surface.as_ref()?;
    let foot = sof.geometry.point_at(u, v);
    let outward = sof.geometry.normal_at(u, v) * if sof.orientation { 1.0 } else { -1.0 };
    let side = (*sample - foot).dot(&outward);
    Some(if side < 0.0 {
        PointClassification::Inside
    } else {
        PointClassification::Outside
    })
}
