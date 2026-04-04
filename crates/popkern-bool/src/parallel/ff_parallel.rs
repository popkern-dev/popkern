//! Parallel face-face intersection dispatch.
//!
//! Divides face pairs into work units and processes them. Currently
//! sequential; designed so that switching to rayon's `par_iter` requires
//! minimal changes.

use popkern_core::geometry::curve_on_surface::CurveOnSurface;
use popkern_core::geometry::intersection::intersect_surfaces;
use popkern_core::geometry::intersection::surface_subdivision::Aabb;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};

use crate::clip::{clip_curve_to_faces, face_vertex_positions};

/// A face pair candidate for intersection.
#[derive(Debug, Clone)]
pub struct FacePairCandidate {
    /// Face from body A.
    pub face_a: FaceHandle,
    /// Face from body B.
    pub face_b: FaceHandle,
}

/// Result of intersecting a single face pair.
#[derive(Debug, Clone)]
pub struct FacePairResult {
    /// Face from body A.
    pub face_a: FaceHandle,
    /// Face from body B.
    pub face_b: FaceHandle,
    /// Intersection curves found (may be empty).
    pub curves: Vec<CurveOnSurface>,
    /// Whether the pair was coincident.
    pub is_coincident: bool,
}

/// Collect face pair candidates using AABB overlap filtering.
pub fn collect_candidates(
    aabbs_a: &[(FaceHandle, Aabb)],
    aabbs_b: &[(FaceHandle, Aabb)],
    padding: f64,
) -> Vec<FacePairCandidate> {
    let mut candidates = Vec::new();
    for (fa, aa) in aabbs_a {
        for (fb, ab) in aabbs_b {
            if aa.overlaps(ab, padding) {
                candidates.push(FacePairCandidate {
                    face_a: *fa,
                    face_b: *fb,
                });
            }
        }
    }
    candidates
}

/// Process face pairs sequentially.
///
/// This is the current implementation. To parallelize with rayon:
/// replace `candidates.iter()` with `candidates.par_iter()` and
/// collect into a `Vec` via `par_iter().map().collect()`.
pub fn process_pairs_sequential(
    candidates: &[FacePairCandidate],
    a: &BRep,
    b: &BRep,
    tolerance: &Tolerance,
) -> Vec<FacePairResult> {
    candidates
        .iter()
        .filter_map(|candidate| process_one_pair(candidate, a, b, tolerance))
        .collect()
}

/// Process face pairs in parallel using rayon.
///
/// Requires the `parallel` feature flag. Each face pair is independent
/// and can be processed on a separate thread.
#[cfg(feature = "parallel")]
pub fn process_pairs_parallel(
    candidates: &[FacePairCandidate],
    a: &BRep,
    b: &BRep,
    tolerance: &Tolerance,
) -> Vec<FacePairResult> {
    use rayon::prelude::*;
    candidates
        .par_iter()
        .filter_map(|candidate| process_one_pair(candidate, a, b, tolerance))
        .collect()
}

/// Process a single face pair.
fn process_one_pair(
    candidate: &FacePairCandidate,
    a: &BRep,
    b: &BRep,
    tolerance: &Tolerance,
) -> Option<FacePairResult> {
    let surf_a = a.faces.get(candidate.face_a)?.surface.as_ref()?;
    let surf_b = b.faces.get(candidate.face_b)?.surface.as_ref()?;

    let pos_a = face_vertex_positions(a, candidate.face_a);
    let pos_b = face_vertex_positions(b, candidate.face_b);
    let ((u1_min, u1_max), (v1_min, v1_max)) =
        popkern_core::geometry::uv_projection::bounded_parameter_range(&pos_a, &surf_a.geometry);
    let ((u2_min, u2_max), (v2_min, v2_max)) =
        popkern_core::geometry::uv_projection::bounded_parameter_range(&pos_b, &surf_b.geometry);

    let result = intersect_surfaces(
        &surf_a.geometry,
        (u1_min, u1_max),
        (v1_min, v1_max),
        &surf_b.geometry,
        (u2_min, u2_max),
        (v2_min, v2_max),
        tolerance,
    )
    .ok()?;

    let mut curves = Vec::new();
    let mut is_coincident = false;

    for isect in result.intersections {
        match isect {
            popkern_core::geometry::intersection::SurfaceSurfaceIntersection::Curve(curve) => {
                if let Some(clipped) =
                    clip_curve_to_faces(a, candidate.face_a, b, candidate.face_b, &curve)
                {
                    curves.push(clipped);
                }
            }
            popkern_core::geometry::intersection::SurfaceSurfaceIntersection::CoincidentRegion {
                ..
            } => {
                is_coincident = true;
            }
        }
    }

    Some(FacePairResult {
        face_a: candidate.face_a,
        face_b: candidate.face_b,
        curves,
        is_coincident,
    })
}
