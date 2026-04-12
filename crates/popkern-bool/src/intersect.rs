//! Stage 1: Face-face intersection curve computation.
//!
//! Given two BRep solids, finds all intersection curves between their faces
//! using AABB filtering followed by surface-surface intersection.

use popkern_core::geometry::curve_on_surface::CurveOnSurface;
use popkern_core::{
    geometry::{
        CurveEval,
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

/// Boundary type at an intersection curve endpoint (spec §2.7).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
#[allow(dead_code)]
pub(crate) enum EndpointBoundaryType {
    /// Exits through an edge of Face A.
    #[default]
    OnEdgeA,
    /// Exits through an edge of Face B.
    OnEdgeB,
    /// At a boundary of both faces.
    OnEdgeBoth,
    /// At an existing vertex.
    OnVertex,
    /// Meets another intersection curve (branch point).
    OnBranchPoint,
    /// In the face interior (closed loop endpoint).
    Interior,
}

/// A single face-face intersection record.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub(crate) struct IntersectionRecord {
    /// Face handle from solid A.
    pub face_a: FaceHandle,
    /// Face handle from solid B.
    pub face_b: FaceHandle,
    /// The intersection curve on both surfaces.
    pub curve: CurveOnSurface,
    /// True if surfaces are tangent along this curve.
    #[allow(dead_code)]
    pub is_tangent: bool,
    /// True if this bounds a coincident-face region.
    pub is_coincident: bool,
    /// Boundary type at the start endpoint.
    pub start_boundary: EndpointBoundaryType,
    /// Boundary type at the end endpoint.
    pub end_boundary: EndpointBoundaryType,
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
    find_intersections_with_aabbs(a, b, None, None, tolerance)
}

/// Find intersections with optional precomputed AABBs.
///
/// When `pre_aabbs_a` or `pre_aabbs_b` is Some, uses the provided AABBs
/// instead of recomputing. This avoids redundant work in multi-body
/// operations where the same body's AABBs are needed multiple times.
pub(crate) fn find_intersections_with_aabbs(
    a: &BRep,
    b: &BRep,
    pre_aabbs_a: Option<&[(FaceHandle, Aabb)]>,
    pre_aabbs_b: Option<&[(FaceHandle, Aabb)]>,
    tolerance: &Tolerance,
) -> Result<IntersectionResult, BoolDiagnostics> {
    let owned_aabbs_a;
    let owned_aabbs_b;
    let aabbs_a = match pre_aabbs_a {
        Some(a) => a,
        None => {
            owned_aabbs_a = compute_face_aabbs(a);
            &owned_aabbs_a
        }
    };
    let aabbs_b = match pre_aabbs_b {
        Some(b) => b,
        None => {
            owned_aabbs_b = compute_face_aabbs(b);
            &owned_aabbs_b
        }
    };
    let mut records: Vec<IntersectionRecord> = Vec::new();

    // BVH-accelerated face-pair filtering. The cellular assembly is
    // order-independent, so BVH pair ordering is safe.
    use popkern_core::validation::bvh::Bvh;
    let n_a = aabbs_a.len();
    let mut combined: Vec<(usize, Aabb)> = Vec::with_capacity(n_a + aabbs_b.len());
    for (i, (_, aabb)) in aabbs_a.iter().enumerate() {
        combined.push((i, aabb.clone()));
    }
    for (i, (_, aabb)) in aabbs_b.iter().enumerate() {
        combined.push((n_a + i, aabb.clone()));
    }
    let bvh = Bvh::build(&combined);
    let candidate_pairs = bvh.find_overlapping_pairs(tolerance.linear());

    for (id_x, id_y) in &candidate_pairs {
        let (a_idx, b_idx) = if *id_x < n_a && *id_y >= n_a {
            (*id_x, *id_y - n_a)
        } else if *id_y < n_a && *id_x >= n_a {
            (*id_y, *id_x - n_a)
        } else {
            continue;
        };

        let face_a = aabbs_a[a_idx].0;
        let face_b = aabbs_b[b_idx].0;

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
                            if let Some(clipped) = clip_curve_to_faces(a, face_a, b, face_b, &curve)
                            {
                                // Classify tangency at curve midpoint
                                let (t0, t1) = clipped.curve_3d.parameter_range();
                                let t_mid = (t0 + t1) * 0.5;
                                let mid_uv1 = clipped.curve_uv1.point_at(t_mid);
                                let mid_uv2 = clipped.curve_uv2.point_at(t_mid);
                                let n1 = surf_a.normal_at(mid_uv1.x, mid_uv1.y);
                                let n2 = surf_b.normal_at(mid_uv2.x, mid_uv2.y);
                                let cross_len = n1.cross(&n2).length();
                                let norms = n1.length() * n2.length();
                                let is_tangent = norms > tolerance.numerical_floor()
                                    && (cross_len / norms) < tolerance.angular();

                                records.push(IntersectionRecord {
                                    face_a,
                                    face_b,
                                    curve: clipped,
                                    is_tangent,
                                    is_coincident: false,
                                    start_boundary: Default::default(),
                                    end_boundary: Default::default(),
                                });
                            }
                        }
                        SurfaceSurfaceIntersection::CoincidentRegion { .. } => {
                            // Coincident face pair — skip. Coincidence engine
                            // records the overlap separately.
                        }
                    }
                }
            }
            Err(e) => {
                // §13.1: report a 3D location for the failure. The midpoint
                // of face_a's bounded UV range is the closest cheap proxy
                // for "where" the intersection solver gave up.
                let u_mid = 0.5 * (u1_min + u1_max);
                let v_mid = 0.5 * (v1_min + v1_max);
                let location = surf_a.point_at(u_mid, v_mid);
                return Err(BoolDiagnostics {
                    error: Some(BoolError::IntersectionFailed(format!(
                        "face_a={face_a:?} face_b={face_b:?}: {e}"
                    ))),
                    intersection_curves: records.into_iter().map(|r| r.curve).collect(),
                    face_classifications: Vec::new(),
                    failed_stage: Some(Stage::Intersection),
                    entries: vec![crate::error::DiagnosticEntry {
                        severity: crate::error::Severity::Error,
                        phase: Stage::Intersection,
                        message: format!("Surface intersection failed: {e}"),
                        entities: vec![crate::error::EntityRef::FacePair(face_a, face_b)],
                        location: Some(location),
                        measured_value: None,
                        remediation: Some(crate::error::Remediation {
                            action: crate::error::RemediationAction::RetryWithPerturbation,
                            parameter: None,
                            value: None,
                            rationale: "Surface intersection solver failed. \
                                Try boolean_with_perturbation() to slightly \
                                offset the geometry and retry."
                                .into(),
                        }),
                    }],
                    phase_metrics: Vec::new(),
                });
            }
        }
    }

    // §4.2 branch-point labeling: if two intersection-curve endpoints
    // are coincident within tolerance, mark both as `OnBranchPoint`.
    // The label is informational today (downstream consumers are
    // gated by `#[allow(dead_code)]`) but is required by future split
    // logic to preserve branch topology when subdividing faces.
    label_branch_point_endpoints(&mut records, tolerance);

    // Diagnostic: detect edge-face intersections (legacy path).
    let edge_face_hits = detect_edge_face_intersections(a, b, aabbs_b, aabbs_a, tolerance);

    // Compute the sub-FF intersection hierarchy.
    let bool_tol = crate::tolerance::BoolTolerance::from(*tolerance);
    let hierarchy =
        crate::intersect_hierarchy::compute_hierarchy(a, b, aabbs_a, aabbs_b, &bool_tol);

    // Build the coincidence map from dedicated coincidence detection.
    let bool_tol2 = crate::tolerance::BoolTolerance::from(*tolerance);
    let coincidence =
        crate::coincidence::detect::detect_coincidences(a, b, aabbs_a, aabbs_b, &bool_tol2);

    Ok(IntersectionResult {
        records,
        aabbs_a: aabbs_a.to_vec(),
        aabbs_b: aabbs_b.to_vec(),
        edge_face_hits,
        hierarchy,
        coincidence,
    })
}

// ---------------------------------------------------------------------------
// §4.2 branch-point labeling
// ---------------------------------------------------------------------------

/// Label endpoints of intersection records that coincide with another
/// curve's endpoint as `OnBranchPoint`.
///
/// Two endpoints "coincide" when their 3D positions are within
/// `tolerance.linear() * 10.0`. The 10× factor matches the consumed
/// proxy used by the SSI marching loop in `popkern-core`, so it
/// catches the same junctions the solver itself treats as coincident.
///
/// The check is O(n²) in the number of intersection records, which is
/// fine for typical boolean operations (well under 100 curves). The
/// label is consumed by future split logic that needs to preserve
/// branch topology when subdividing faces; today it is informational.
fn label_branch_point_endpoints(records: &mut [IntersectionRecord], tolerance: &Tolerance) {
    use popkern_core::geometry::CurveEval;
    if records.len() < 2 {
        return;
    }

    let snap = tolerance.linear() * 10.0;
    let snap_sq = snap * snap;

    // Snapshot every endpoint with its (record_idx, is_start) location.
    let endpoints: Vec<(usize, bool, popkern_core::geometry::Point3)> = records
        .iter()
        .enumerate()
        .flat_map(|(i, r)| {
            let (t0, t1) = r.curve.curve_3d.parameter_range();
            let p0 = r.curve.curve_3d.point_at(t0);
            let p1 = r.curve.curve_3d.point_at(t1);
            [(i, true, p0), (i, false, p1)]
        })
        .collect();

    let mut is_branch = vec![[false; 2]; records.len()];

    for i in 0..endpoints.len() {
        for j in (i + 1)..endpoints.len() {
            // Endpoints from the same record are not branches between
            // each other (they are the two ends of one curve).
            if endpoints[i].0 == endpoints[j].0 {
                continue;
            }
            let dx = endpoints[i].2.x - endpoints[j].2.x;
            let dy = endpoints[i].2.y - endpoints[j].2.y;
            let dz = endpoints[i].2.z - endpoints[j].2.z;
            if dx * dx + dy * dy + dz * dz < snap_sq {
                let (idx_i, start_i, _) = endpoints[i];
                let (idx_j, start_j, _) = endpoints[j];
                is_branch[idx_i][if start_i { 0 } else { 1 }] = true;
                is_branch[idx_j][if start_j { 0 } else { 1 }] = true;
            }
        }
    }

    for (i, record) in records.iter_mut().enumerate() {
        if is_branch[i][0] {
            record.start_boundary = EndpointBoundaryType::OnBranchPoint;
        }
        if is_branch[i][1] {
            record.end_boundary = EndpointBoundaryType::OnBranchPoint;
        }
    }
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
        if direction.length_squared() < tolerance.numerical_floor().powi(2) {
            continue;
        }
        let edge_curve = CurveGeometry::Line {
            origin: v_start,
            direction,
        };
        let t_range = (0.0, 1.0);

        for (face_h, face_aabb) in face_aabbs {
            if !edge_aabb.overlaps(face_aabb, tolerance.linear()) {
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
                        let eps = tolerance.linear();
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
