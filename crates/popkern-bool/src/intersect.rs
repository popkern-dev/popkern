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

use crate::clip::clip_curve_to_faces;
use crate::error::{BoolDiagnostics, BoolError, Stage};
use crate::intersect_edge_face::{EdgeFaceHit, detect_edge_face_intersections};

/// Compute bounded parameter ranges for a face's surface.
///
/// Delegates to the shared `popkern_core::geometry::uv_projection` helper,
/// which handles analytic UV projection for planes, cylinders, cones, and
/// spheres.
pub(super) fn bounded_parameter_range(
    brep: &BRep,
    face_handle: FaceHandle,
    surface: &popkern_core::geometry::surface::SurfaceGeometry,
) -> ((f64, f64), (f64, f64)) {
    // Boundary CURVE samples, not just vertices: a face bounded by
    // a closed rim (cylinder/boss cap) has one vertex, and a
    // vertex-only box collapses to a point at the seam — every
    // pocket cut into such a cap silently no-opped because the
    // intersection search never covered the face interior.
    let positions = brep.face_outer_loop_sample_positions(face_handle, 8);
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

    // Boundary CURVE samples, not just vertices: a planar face
    // bounded by a closed rim (cylinder cap) has ONE vertex, and a
    // vertex-only AABB collapses to the seam point — every face
    // pair with such a cap was pruned before intersection ran
    // (measured: pockets into swept-cylinder caps silently
    // no-opped). Curve samples span the true extent; for straight
    // polygonal loops they add interpolants on the chords, which
    // never widen the box.
    points.extend(brep.face_outer_loop_sample_positions(face_handle, 8));

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

    if std::env::var_os("PK_PAIR_TRACE").is_some() {
        eprintln!("  [pairs] {} candidate pairs", candidate_pairs.len());
    }
    for (id_x, id_y) in &candidate_pairs {
        if std::env::var_os("PK_PAIR_TRACE").is_some() {
            eprintln!("  [pair] {id_x:?} x {id_y:?}");
        }
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
                            let clipped_opt = clip_curve_to_faces(a, face_a, b, face_b, &curve);
                            if std::env::var_os("PK_CLIP_TRACE").is_some() {
                                use popkern_core::geometry::curve::CurveEval;
                                let (t0, t1) = curve.curve_3d.parameter_range();
                                let p0 = curve.curve_3d.point_at(t0);
                                let p1 = curve.curve_3d.point_at(t1);
                                eprintln!(
                                    "  [clip] A{face_a:?} x B{face_b:?}: t[{t0:.3},{t1:.3}] ({:.3},{:.3},{:.3})->({:.3},{:.3},{:.3}) kept={}",
                                    p0.x,
                                    p0.y,
                                    p0.z,
                                    p1.x,
                                    p1.y,
                                    p1.z,
                                    clipped_opt.is_some()
                                );
                            }
                            if let Some(clipped) = clipped_opt {
                                // Classify tangency at curve midpoint.
                                // Uses popkern_core::adaptive_precision::
                                // adaptive_cross_ratio_sq_lt so near-parallel
                                // normals produce deterministic decisions —
                                // the f64 fast path handles clear cases,
                                // and double-precision recomputation kicks
                                // in only when the comparison lands inside
                                // the ambiguity margin.
                                let (t0, t1) = clipped.curve_3d.parameter_range();
                                let t_mid = (t0 + t1) * 0.5;
                                let mid_uv1 = clipped.curve_uv1.point_at(t_mid);
                                let mid_uv2 = clipped.curve_uv2.point_at(t_mid);
                                let n1 = surf_a.normal_at(mid_uv1.x, mid_uv1.y);
                                let n2 = surf_b.normal_at(mid_uv2.x, mid_uv2.y);
                                let is_tangent =
                                    popkern_core::adaptive_precision::adaptive_cross_ratio_sq_lt(
                                        (n1.x, n1.y, n1.z),
                                        (n2.x, n2.y, n2.z),
                                        tolerance.angular(),
                                    );

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

    // Consolidate the marching output: the walker can deliver the
    // same intersection curve several times (different walk seeds,
    // split at different offsets) and fragment junctions carry up to
    // ~1e-3 of marching slop. Every downstream consumer (split-stage
    // face splitting, the cellular graph's endpoint snapping) needs
    // deduplicated fragments with welded junctions.
    consolidate_records(&mut records, tolerance);

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

#[cfg(test)]
#[path = "intersect_tests.rs"]
mod tests;

/// Consolidate marching output per face pair: drop geometric
/// duplicates of the same curve piece and weld junction endpoints of
/// pieces that chain (endpoint slop up to the marching class,
/// measured 3.5e-4..1.2e-3 — far above the linear tolerance, far
/// below any fragment's own extent).
fn consolidate_records(records: &mut Vec<IntersectionRecord>, tolerance: &Tolerance) {
    use popkern_core::geometry::curve::CurveEval as _;
    use std::collections::HashMap;
    if records.len() < 2 {
        return;
    }
    let mut groups: HashMap<(FaceHandle, FaceHandle), Vec<usize>> = HashMap::new();
    for (i, r) in records.iter().enumerate() {
        groups.entry((r.face_a, r.face_b)).or_default().push(i);
    }
    let mut drop = vec![false; records.len()];
    // Junction weld targets: (record index, is_start, new point set).
    for idxs in groups.values() {
        if idxs.len() < 2 {
            continue;
        }
        let info: Vec<(
            popkern_core::geometry::Point3,
            popkern_core::geometry::Point3,
            popkern_core::geometry::Point3,
            f64,
        )> = idxs
            .iter()
            .map(|&i| {
                let c = &records[i].curve.curve_3d;
                let (t0, t1) = c.parameter_range();
                let s = c.point_at(t0);
                let e = c.point_at(t1);
                let m = c.point_at(0.5 * (t0 + t1));
                let chord = s.distance_to(&e);
                (s, e, m, chord)
            })
            .collect();
        // Median chord: a single short fragment must not collapse
        // the junction gate below the marching slop class.
        let mut chords: Vec<f64> = info.iter().map(|x| x.3).collect();
        chords.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        let median_chord = chords[chords.len() / 2].max(tolerance.linear());
        let gate = (0.05 * median_chord).max(tolerance.linear() * 10.0);
        // Duplicates: same endpoints (either direction) + same midpoint.
        for a in 0..idxs.len() {
            if drop[idxs[a]] {
                continue;
            }
            for b in (a + 1)..idxs.len() {
                if drop[idxs[b]] {
                    continue;
                }
                let fwd = info[a].0.distance_to(&info[b].0) < gate
                    && info[a].1.distance_to(&info[b].1) < gate;
                let rev = info[a].0.distance_to(&info[b].1) < gate
                    && info[a].1.distance_to(&info[b].0) < gate;
                if (fwd || rev) && info[a].2.distance_to(&info[b].2) < gate {
                    drop[idxs[b]] = true;
                }
            }
        }
        // Loop-level dedupe: the walker can deliver the SAME closed
        // intersection loop several times, split at different
        // parameter offsets (piecewise dedupe can't see those).
        // Chain the survivors by endpoints; closed chains whose
        // sampled geometry matches an earlier chain are dropped
        // whole.
        {
            let alive: Vec<usize> = idxs.iter().copied().filter(|&i| !drop[i]).collect();
            let ends_of = |i: usize| -> (
                popkern_core::geometry::Point3,
                popkern_core::geometry::Point3,
            ) {
                let c = &records[i].curve.curve_3d;
                let (t0, t1) = c.parameter_range();
                (c.point_at(t0), c.point_at(t1))
            };
            let mut used = vec![false; alive.len()];
            let mut kept_loops: Vec<Vec<popkern_core::geometry::Point3>> = Vec::new();
            for seed in 0..alive.len() {
                if used[seed] {
                    continue;
                }
                let mut chain: Vec<(usize, bool)> = vec![(seed, false)];
                let mut chain_used = used.clone();
                chain_used[seed] = true;
                let (start, mut tip) = ends_of(alive[seed]);
                let closed = loop {
                    if tip.distance_to(&start) < gate && chain.len() >= 2 {
                        break true;
                    }
                    let mut ext = false;
                    for k in 0..alive.len() {
                        if chain_used[k] {
                            continue;
                        }
                        let (ks, ke) = ends_of(alive[k]);
                        if ks.distance_to(&tip) < gate {
                            chain.push((k, false));
                            chain_used[k] = true;
                            tip = ke;
                            ext = true;
                            break;
                        }
                        if ke.distance_to(&tip) < gate {
                            chain.push((k, true));
                            chain_used[k] = true;
                            tip = ks;
                            ext = true;
                            break;
                        }
                    }
                    if !ext {
                        break false;
                    }
                };
                if std::env::var_os("PK_BOOL_TRACE").is_some() {
                    eprintln!(
                        "  [consolidate] chain from seed {seed}: len {} closed={closed}",
                        chain.len()
                    );
                }
                if !closed {
                    continue;
                }
                // Sample in chain-traversal order (reversed pieces
                // sample backward — otherwise the comparison
                // polyline zigzags across the loop).
                let samples: Vec<popkern_core::geometry::Point3> = chain
                    .iter()
                    .flat_map(|&(k, rev)| {
                        let c = &records[alive[k]].curve.curve_3d;
                        let (t0, t1) = c.parameter_range();
                        (0..8).map(move |q| {
                            let f = q as f64 / 8.0;
                            let f = if rev { 1.0 - f } else { f };
                            c.point_at(t0 + (t1 - t0) * f)
                        })
                    })
                    .collect();
                // Distance to the kept loop's POLYLINE (segment
                // projection) — sample-to-sample comparison misses
                // duplicates whose samples interleave.
                let dist_to_polyline =
                    |p: &popkern_core::geometry::Point3,
                     poly: &[popkern_core::geometry::Point3]| {
                        let mut best = f64::INFINITY;
                        for w in 0..poly.len() {
                            let a = poly[w];
                            let b = poly[(w + 1) % poly.len()];
                            let ab = b - a;
                            let len2 = ab.dot(&ab);
                            let t = if len2 > 1e-300 {
                                ((*p - a).dot(&ab) / len2).clamp(0.0, 1.0)
                            } else {
                                0.0
                            };
                            let q = a + ab * t;
                            best = best.min(p.distance_to(&q));
                        }
                        best
                    };
                let duplicate = kept_loops
                    .iter()
                    .any(|kl| samples.iter().all(|p| dist_to_polyline(p, kl) < gate));
                if std::env::var_os("PK_BOOL_TRACE").is_some() {
                    let n = samples.len() as f64;
                    let cx = samples.iter().map(|p| p.x).sum::<f64>() / n;
                    let cy = samples.iter().map(|p| p.y).sum::<f64>() / n;
                    let cz = samples.iter().map(|p| p.z).sum::<f64>() / n;
                    let c = popkern_core::geometry::Point3::new(cx, cy, cz);
                    let rmin = samples
                        .iter()
                        .map(|p| p.distance_to(&c))
                        .fold(f64::INFINITY, f64::min);
                    let rmax = samples
                        .iter()
                        .map(|p| p.distance_to(&c))
                        .fold(0.0_f64, f64::max);
                    eprintln!(
                        "  [consolidate] duplicate={duplicate} centroid ({cx:.3},{cy:.3},{cz:.3}) r {rmin:.4}..{rmax:.4} gate {gate:.2e}"
                    );
                }
                used = chain_used;
                if duplicate {
                    for &(k, _) in &chain {
                        drop[alive[k]] = true;
                    }
                } else {
                    kept_loops.push(samples);
                }
            }
        }

        // Self-weld: a nearly-closed loop record (closure gap within
        // the marching slop) becomes exactly closed — 3D and both UV
        // charts — so the graph sees one node, not an open arc with
        // a 1e-3 gap.
        for &i in idxs.iter().filter(|&&i| !drop[i]) {
            let (s3, e3) = {
                let c = &records[i].curve.curve_3d;
                let (t0, t1) = c.parameter_range();
                (c.point_at(t0), c.point_at(t1))
            };
            let d = s3.distance_to(&e3);
            if d > tolerance.linear() && d < gate {
                let weld = popkern_core::geometry::Point3::new(
                    0.5 * (s3.x + e3.x),
                    0.5 * (s3.y + e3.y),
                    0.5 * (s3.z + e3.z),
                );
                set_endpoint(&mut records[i], false, weld);
                set_endpoint(&mut records[i], true, weld);
                for uv in [true, false] {
                    let c = if uv {
                        &mut records[i].curve.curve_uv1
                    } else {
                        &mut records[i].curve.curve_uv2
                    };
                    if let (Some(first), Some(last)) = (
                        c.control_points.first().copied(),
                        c.control_points.last().copied(),
                    ) {
                        let m = popkern_core::geometry::Point3::new(
                            0.5 * (first.x + last.x),
                            0.5 * (first.y + last.y),
                            0.5 * (first.z + last.z),
                        );
                        if let Some(cp) = c.control_points.first_mut() {
                            *cp = m;
                        }
                        if let Some(cp) = c.control_points.last_mut() {
                            *cp = m;
                        }
                    }
                }
            }
        }

        // Weld junctions among the survivors: cluster endpoints and
        // move each cluster's endpoints (3D + both UV charts stay
        // consistent because the curves are only ADJUSTED at their
        // endpoint control points, which both charts share).
        let survivors: Vec<usize> = idxs.iter().copied().filter(|&i| !drop[i]).collect();
        for a_pos in 0..survivors.len() {
            for b_pos in (a_pos + 1)..survivors.len() {
                let (ia, ib) = (survivors[a_pos], survivors[b_pos]);
                for (a_end, b_end) in [(false, true), (true, false), (false, false), (true, true)] {
                    let pa = endpoint_of(&records[ia], a_end);
                    let pb = endpoint_of(&records[ib], b_end);
                    let d = pa.distance_to(&pb);
                    if d > tolerance.linear() && d < gate {
                        let weld = popkern_core::geometry::Point3::new(
                            0.5 * (pa.x + pb.x),
                            0.5 * (pa.y + pb.y),
                            0.5 * (pa.z + pb.z),
                        );
                        set_endpoint(&mut records[ia], a_end, weld);
                        set_endpoint(&mut records[ib], b_end, weld);
                    }
                }
            }
        }
    }
    if std::env::var_os("PK_BOOL_TRACE").is_some() {
        let dropped = drop.iter().filter(|&&d| d).count();
        eprintln!(
            "  [consolidate] {} records in {} groups, dropped {}",
            records.len(),
            groups.len(),
            dropped
        );
        for (k, v) in &groups {
            eprintln!("    group {k:?}: {} records", v.len());
        }
    }
    let mut keep_iter = drop.into_iter();
    records.retain(|_| !keep_iter.next().unwrap_or(false));
}

/// Endpoint position of a record's 3D curve (`end = true` for t_end).
fn endpoint_of(record: &IntersectionRecord, end: bool) -> popkern_core::geometry::Point3 {
    use popkern_core::geometry::curve::CurveEval as _;
    let c = &record.curve.curve_3d;
    let (t0, t1) = c.parameter_range();
    c.point_at(if end { t1 } else { t0 })
}

/// Move a record's 3D endpoint control point to `p` (clamped NURBS:
/// the endpoint IS a control point).
fn set_endpoint(record: &mut IntersectionRecord, end: bool, p: popkern_core::geometry::Point3) {
    let c = &mut record.curve.curve_3d;
    if end {
        if let Some(cp) = c.control_points.last_mut() {
            *cp = p;
        }
    } else if let Some(cp) = c.control_points.first_mut() {
        *cp = p;
    }
}

#[cfg(test)]
#[path = "intersect_edge_face_tests.rs"]
mod edge_face_tests;
