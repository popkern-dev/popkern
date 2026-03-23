//! Stage 1: Face-face intersection curve computation.
//!
//! Given two BRep solids, finds all intersection curves between their faces
//! using AABB filtering followed by surface-surface intersection.

use popkern_core::geometry::curve_on_surface::CurveOnSurface;
use popkern_core::{
    geometry::{
        intersection::{SurfaceSurfaceIntersection, intersect_surfaces, surface_subdivision::Aabb},
        surface::SurfaceEval,
    },
    tolerance::Tolerance,
    topology::{BRep, FaceHandle},
};

use crate::clip::{clip_curve_to_faces, face_vertex_positions};
use crate::error::{BoolDiagnostics, BoolError, Stage};

/// Compute bounded parameter ranges for a face's surface.
///
/// Delegates to the shared `popkern_core::geometry::uv_projection` helper,
/// which handles analytic UV projection for planes, cylinders, cones, and
/// spheres.
fn bounded_parameter_range(
    brep: &BRep,
    face_handle: FaceHandle,
    surface: &popkern_core::geometry::surface::SurfaceGeometry,
) -> ((f64, f64), (f64, f64)) {
    let positions = face_vertex_positions(brep, face_handle);
    popkern_core::geometry::uv_projection::bounded_parameter_range(&positions, surface)
}

/// A single face-face intersection record.
#[derive(Debug, Clone)]
pub(crate) struct IntersectionRecord {
    /// Face handle from solid A.
    pub face_a: FaceHandle,
    /// Face handle from solid B.
    pub face_b: FaceHandle,
    /// The intersection curve on both surfaces.
    pub curve: CurveOnSurface,
}

/// Complete result of Stage 1 intersection computation.
#[derive(Debug)]
pub(crate) struct IntersectionResult {
    /// All face-face intersection records found.
    pub records: Vec<IntersectionRecord>,
    /// AABBs computed for each face of solid A.
    pub aabbs_a: Vec<(FaceHandle, Aabb)>,
    /// AABBs computed for each face of solid B.
    pub aabbs_b: Vec<(FaceHandle, Aabb)>,
}

/// Sample a face to compute its AABB.
///
/// For faces with finite parameter ranges (e.g. NURBS, sphere, torus), samples
/// a 5×5 grid over the surface. For faces with infinite ranges (e.g. planes),
/// falls back to collecting vertex positions from the face's outer loop.
///
/// Returns `None` if the face has no surface geometry or no valid vertices.
fn face_aabb(brep: &BRep, face_handle: FaceHandle) -> Option<Aabb> {
    let surface_on_face = brep.faces[face_handle].surface.as_ref()?;
    let surf = &surface_on_face.geometry;
    let ((u_min, u_max), (v_min, v_max)) = surf.parameter_range();

    let mut points = Vec::new();

    // Always add vertex positions — accurate for planar faces and provides a
    // safe fallback when parameter range is infinite.
    points.extend(face_vertex_positions(brep, face_handle));

    // If the surface has a finite parameter range, also sample a 5×5 grid for
    // curved surfaces where vertices may not capture the full extent.
    if u_min.is_finite() && u_max.is_finite() && v_min.is_finite() && v_max.is_finite() {
        const N: usize = 5;
        for i in 0..N {
            let u = u_min + (u_max - u_min) * (i as f64) / ((N - 1) as f64);
            for j in 0..N {
                let v = v_min + (v_max - v_min) * (j as f64) / ((N - 1) as f64);
                points.push(surf.point_at(u, v));
            }
        }
    }

    Aabb::from_points(&points)
}

/// Compute AABBs for all faces in a BRep that have surface geometry.
///
/// Returns a list of `(FaceHandle, Aabb)` pairs.
pub(crate) fn compute_face_aabbs(brep: &BRep) -> Vec<(FaceHandle, Aabb)> {
    brep.faces
        .keys()
        .filter_map(|fh| face_aabb(brep, fh).map(|aabb| (fh, aabb)))
        .collect()
}

/// Find all intersection curves between faces of two BRep solids.
///
/// Uses AABB filtering to reduce the number of expensive surface-surface
/// intersection calls. Returns an [`IntersectionResult`] on success, or a
/// [`BoolDiagnostics`] containing the error and any partial results on failure.
///
/// # Arguments
///
/// * `a` — the first (left-hand) solid.
/// * `b` — the second (right-hand) solid.
/// * `tolerance` — geometric tolerance for all distance and coincidence tests.
pub(crate) fn find_intersections(
    a: &BRep,
    b: &BRep,
    tolerance: &Tolerance,
) -> Result<IntersectionResult, BoolDiagnostics> {
    let aabbs_a = compute_face_aabbs(a);
    let aabbs_b = compute_face_aabbs(b);
    let mut records: Vec<IntersectionRecord> = Vec::new();

    for (face_a, aabb_a) in &aabbs_a {
        for (face_b, aabb_b) in &aabbs_b {
            if !aabb_a.overlaps(aabb_b, tolerance.linear) {
                continue;
            }

            // Both faces are guaranteed to have surface geometry because
            // compute_face_aabbs only includes faces with surfaces.
            let surf_a = match a.faces[*face_a].surface.as_ref() {
                Some(s) => &s.geometry,
                None => continue,
            };
            let surf_b = match b.faces[*face_b].surface.as_ref() {
                Some(s) => &s.geometry,
                None => continue,
            };

            let ((u1_min, u1_max), (v1_min, v1_max)) = bounded_parameter_range(a, *face_a, surf_a);
            let ((u2_min, u2_max), (v2_min, v2_max)) = bounded_parameter_range(b, *face_b, surf_b);

            match intersect_surfaces(
                surf_a,
                (u1_min, u1_max),
                (v1_min, v1_max),
                surf_b,
                (u2_min, u2_max),
                (v2_min, v2_max),
                tolerance,
            ) {
                Ok(result) => {
                    for isect in result.intersections {
                        match isect {
                            SurfaceSurfaceIntersection::Curve(curve) => {
                                // Clip the curve to both face boundaries.
                                match clip_curve_to_faces(a, *face_a, b, *face_b, &curve) {
                                    Some(clipped) => {
                                        records.push(IntersectionRecord {
                                            face_a: *face_a,
                                            face_b: *face_b,
                                            curve: clipped,
                                        });
                                    }
                                    None => {
                                        // Curve doesn't intersect both face
                                        // polygons — skip it.
                                    }
                                }
                            }
                            SurfaceSurfaceIntersection::CoincidentRegion { .. } => {
                                return Err(BoolDiagnostics {
                                    error: Some(BoolError::CoincidentFaces {
                                        face_a: *face_a,
                                        face_b: *face_b,
                                    }),
                                    intersection_curves: records
                                        .into_iter()
                                        .map(|r| r.curve)
                                        .collect(),
                                    face_classifications: Vec::new(),
                                    failed_stage: Some(Stage::Intersection),
                                });
                            }
                        }
                    }
                }
                Err(e) => {
                    return Err(BoolDiagnostics {
                        error: Some(BoolError::IntersectionFailed(format!(
                            "face_a={face_a:?} face_b={face_b:?}: {e}"
                        ))),
                        intersection_curves: records.into_iter().map(|r| r.curve).collect(),
                        face_classifications: Vec::new(),
                        failed_stage: Some(Stage::Intersection),
                    });
                }
            }
        }
    }

    Ok(IntersectionResult {
        records,
        aabbs_a,
        aabbs_b,
    })
}

#[cfg(test)]
#[path = "intersect_tests.rs"]
mod tests;
