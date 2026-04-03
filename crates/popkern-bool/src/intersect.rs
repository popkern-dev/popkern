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
    /// Edge-face intersection hits (diagnostic only, not yet used by split).
    #[allow(dead_code)]
    pub edge_face_hits: Vec<EdgeFaceHit>,
    /// Sub-FF intersection hierarchy results (VV/VE/EE/VF/EF).
    #[allow(dead_code)]
    pub hierarchy: crate::intersect_hierarchy::HierarchyResult,
    /// Coincidence map: vertex merges, edge overlaps, and face overlaps.
    #[allow(dead_code)]
    pub coincidence: crate::coincidence::CoincidenceMap,
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

    // N^2 face-pair loop with AABB pre-filtering.
    // NOTE: BVH acceleration is available via Bvh::find_overlapping_pairs but
    // the graph-merge pipeline is currently order-sensitive, so we preserve the
    // deterministic iteration order until the cellular decomposition (Phase 6)
    // replaces it.
    for (face_a, aabb_a) in &aabbs_a {
        for (face_b, aabb_b) in &aabbs_b {
            if !aabb_a.overlaps(aabb_b, tolerance.linear) {
                continue;
            }

            let face_a = *face_a;
            let face_b = *face_b;

            let surf_a = match a.faces[face_a].surface.as_ref() {
                Some(s) => &s.geometry,
                None => continue,
            };
            let surf_b = match b.faces[face_b].surface.as_ref() {
                Some(s) => &s.geometry,
                None => continue,
            };

            let ((u1_min, u1_max), (v1_min, v1_max)) = bounded_parameter_range(a, face_a, surf_a);
            let ((u2_min, u2_max), (v2_min, v2_max)) = bounded_parameter_range(b, face_b, surf_b);

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
                                if let Some(clipped) =
                                    clip_curve_to_faces(a, face_a, b, face_b, &curve)
                                {
                                    records.push(IntersectionRecord {
                                        face_a,
                                        face_b,
                                        curve: clipped,
                                    });
                                }
                            }
                            SurfaceSurfaceIntersection::CoincidentRegion { .. } => {
                                // Coincident face pair — graph_merge cannot yet handle
                                // these without overlap boundary curves. Return error
                                // until Phase 3.3 overlap boundary computation is done.
                                return Err(BoolDiagnostics {
                                    error: Some(BoolError::CoincidentFaces { face_a, face_b }),
                                    intersection_curves: records
                                        .into_iter()
                                        .map(|r| r.curve)
                                        .collect(),
                                    face_classifications: Vec::new(),
                                    failed_stage: Some(Stage::Intersection),
                                    entries: Vec::new(),
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
                        entries: Vec::new(),
                    });
                }
            }
        }
    }

    // Diagnostic: detect edge-face intersections (legacy path).
    let edge_face_hits = detect_edge_face_intersections(a, b, &aabbs_b, &aabbs_a, tolerance);

    // Compute the sub-FF intersection hierarchy.
    let bool_tol = crate::tolerance::BoolTolerance::from(*tolerance);
    let hierarchy =
        crate::intersect_hierarchy::compute_hierarchy(a, b, &aabbs_a, &aabbs_b, &bool_tol);

    // Build the coincidence map from dedicated coincidence detection.
    let bool_tol2 = crate::tolerance::BoolTolerance::from(*tolerance);
    let coincidence =
        crate::coincidence::detect::detect_coincidences(a, b, &aabbs_a, &aabbs_b, &bool_tol2);

    Ok(IntersectionResult {
        records,
        aabbs_a,
        aabbs_b,
        edge_face_hits,
        hierarchy,
        coincidence,
    })
}

// ---------------------------------------------------------------------------
// Edge-face intersection detection (diagnostic only)
// ---------------------------------------------------------------------------

/// A detected edge-face intersection point.
#[derive(Debug, Clone)]
#[allow(dead_code)] // Fields used for diagnostic inspection and future integration.
pub(crate) struct EdgeFaceHit {
    /// The edge handle from the source solid.
    pub edge: popkern_core::topology::EdgeHandle,
    /// The face handle from the target solid.
    pub face: FaceHandle,
    /// The 3D intersection point.
    pub point: popkern_core::geometry::Point3,
    /// Parameter t on the edge curve at the intersection.
    pub t_edge: f64,
    /// UV parameters (u, v) on the face surface at the intersection.
    pub uv: (f64, f64),
}

/// Detect edge-face intersections between two solids (diagnostic only).
///
/// For each edge of solid A, checks if it intersects any face of solid B,
/// and vice versa. Uses `intersect_surface_curve` with the edge's curve
/// geometry and the face's surface geometry.
///
/// Returns a list of [`EdgeFaceHit`] records. These are informational and
/// are not yet integrated into the split pipeline.
pub(crate) fn detect_edge_face_intersections(
    a: &BRep,
    b: &BRep,
    aabbs_b: &[(FaceHandle, Aabb)],
    aabbs_a: &[(FaceHandle, Aabb)],
    tolerance: &Tolerance,
) -> Vec<EdgeFaceHit> {
    let mut hits = Vec::new();

    // Edges of A vs faces of B.
    detect_edges_vs_faces(a, b, aabbs_b, tolerance, &mut hits);

    // Edges of B vs faces of A.
    detect_edges_vs_faces(b, a, aabbs_a, tolerance, &mut hits);

    if !hits.is_empty() {
        #[cfg(debug_assertions)]
        eprintln!(
            "[intersect] detected {} edge-face intersection(s)",
            hits.len()
        );
    }

    hits
}

/// Check edges of `edge_brep` against faces of `face_brep`.
fn detect_edges_vs_faces(
    edge_brep: &BRep,
    face_brep: &BRep,
    face_aabbs: &[(FaceHandle, Aabb)],
    tolerance: &Tolerance,
    hits: &mut Vec<EdgeFaceHit>,
) {
    use popkern_core::geometry::CurveEval;
    use popkern_core::geometry::intersection::{SurfaceCurveIntersection, intersect_surface_curve};

    use popkern_core::geometry::curve::CurveGeometry;

    for (edge_h, edge) in &edge_brep.edges {
        // Get endpoint positions from vertices (always up-to-date, even if
        // the edge curve geometry is stale after vertex translation).
        let he = match edge_brep.half_edges.get(edge.half_edge) {
            Some(h) => h,
            None => continue,
        };
        let v_start = match edge_brep.vertices.get(he.origin) {
            Some(v) => v.position,
            None => continue,
        };
        let twin = match edge_brep.half_edges.get(he.twin) {
            Some(t) => t,
            None => continue,
        };
        let v_end = match edge_brep.vertices.get(twin.origin) {
            Some(v) => v.position,
            None => continue,
        };

        let edge_aabb = match Aabb::from_points(&[v_start, v_end]) {
            Some(a) => a,
            None => continue,
        };

        // Construct a line curve from vertex positions (handles cases where
        // stored CurveOnEdge geometry is out of sync with vertex positions).
        let direction = v_end - v_start;
        if direction.length_squared() < 1e-30 {
            continue;
        }
        let edge_curve = CurveGeometry::Line {
            origin: v_start,
            direction,
        };
        let t_range = (0.0, 1.0);

        for (face_h, face_aabb) in face_aabbs {
            if !edge_aabb.overlaps(face_aabb, tolerance.linear) {
                continue;
            }

            let face = &face_brep.faces[*face_h];
            let surface_on_face = match face.surface.as_ref() {
                Some(s) => &s.geometry,
                None => continue,
            };

            let ((u_min, u_max), (v_min, v_max)) =
                bounded_parameter_range(face_brep, *face_h, surface_on_face);

            let result = intersect_surface_curve(
                surface_on_face,
                (u_min, u_max),
                (v_min, v_max),
                &edge_curve,
                t_range,
                tolerance,
            );

            if let Ok(result) = result {
                for isect in &result.intersections {
                    if let SurfaceCurveIntersection::Point { t, u, v, .. } = *isect {
                        // Filter out intersections at edge endpoints (these are
                        // vertex-face touches, not edge-face piercings).
                        let eps = tolerance.linear;
                        if (t - t_range.0).abs() < eps || (t - t_range.1).abs() < eps {
                            continue;
                        }
                        hits.push(EdgeFaceHit {
                            edge: edge_h,
                            face: *face_h,
                            point: edge_curve.point_at(t),
                            t_edge: t,
                            uv: (u, v),
                        });
                    }
                }
            }
        }
    }
}

#[cfg(test)]
#[path = "intersect_tests.rs"]
mod tests;

#[cfg(test)]
#[path = "intersect_edge_face_tests.rs"]
mod edge_face_tests;
