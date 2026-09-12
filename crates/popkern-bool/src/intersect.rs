//! Stage 1: Face-face intersection curve computation.
//!
//! Given two BRep solids, finds all intersection curves between their faces
//! using AABB filtering followed by surface-surface intersection.

use popkern_core::geometry::curve_on_surface::CurveOnSurface;
use popkern_core::{
    geometry::{intersection::surface_subdivision::Aabb, surface::SurfaceEval},
    tolerance::Tolerance,
    topology::{BRep, FaceHandle},
};

use crate::error::BoolDiagnostics;
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
    /// TRIM-RESIDUE marker (PK-465): `Some(true)` = this piece lies
    /// inside face A's trim but OUTSIDE face B's — it partitions face
    /// A only (registered to A on both sides, like a seam) and never
    /// becomes a two-face result edge. `Some(false)` = the mirror
    /// (face B's residue). `None` = an ordinary shared record. The
    /// residue is the other operand's SURFACE crossing beyond the
    /// shared-edge window: without it, regions of one face that touch
    /// only across that crossing walk as ONE connected fragment and
    /// classify with one (necessarily wrong) verdict — measured on
    /// the angled tap's flank, whose above-plate and in-material
    /// parts merged around each z-plane record endpoint.
    pub residue_owner_a: Option<bool>,
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

    // CONSERVATIVE curved-edge extents. Finite sampling of a curved
    // boundary UNDER-approximates its extent — 8 samples of a
    // rational circle (whose parameter is not angle-uniform) missed
    // the 180° extreme by 2.2 mm on an r = 36 bore, the box stopped
    // at x = −58.8 instead of −61, the breached wall at x = −60 fell
    // outside it, and the BVH never paired the faces: the whole
    // plane×cylinder curve set went missing and the subtract failed
    // LowValence. An AABB used as a pruning gate must CONTAIN the
    // face; over-covering only costs candidate pairs. Analytic
    // curves get exact per-axis amplitude boxes (full-sweep
    // superset for arcs), NURBS curves their control-point hull
    // (positive weights ⇒ curve ⊂ hull), parabola/hyperbola their
    // endpoint + interior derivative-zero extremes.
    let mut extend_curve_extent =
        |c: &popkern_core::geometry::curve::CurveOnEdge| conservative_curve_extent(c, &mut points);
    let face = &brep.faces[face_handle];
    let mut loops = vec![face.outer_loop];
    loops.extend(face.inner_loops.iter().copied());
    for lh in loops {
        let Some(lp) = brep.loops.get(lh) else {
            continue;
        };
        let start = lp.half_edge;
        let mut he = start;
        for _ in 0..brep.half_edges.len() + 1 {
            let Some(h) = brep.half_edges.get(he) else {
                break;
            };
            if let Some(edge) = brep.edges.get(h.edge)
                && let Some(c) = edge.curve.as_ref()
            {
                extend_curve_extent(c);
            }
            he = h.next;
            if he == start {
                break;
            }
        }
    }

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

/// Push CONSERVATIVE extent points for a curve into `points` — a
/// superset of the curve's true axis-aligned extent. Finite sampling
/// UNDER-approximates curved extents (a rational circle's parameter
/// is not angle-uniform; a closed rim's two vertices coincide at the
/// seam), and AABBs used as pruning gates must CONTAIN the geometry.
/// Analytic curves contribute exact per-axis amplitude boxes
/// (full-sweep superset for arcs), NURBS their control-point hull
/// (positive weights ⇒ curve ⊂ hull), parabola/hyperbola their
/// interior derivative-zero extremes.
pub(crate) fn conservative_curve_extent(
    c: &popkern_core::geometry::curve::CurveOnEdge,
    points: &mut Vec<popkern_core::geometry::Point3>,
) {
    use popkern_core::geometry::curve::{CurveEval, CurveGeometry};
    match &c.geometry {
        CurveGeometry::Line { .. } => {} // endpoints suffice
        CurveGeometry::Circle {
            center,
            normal,
            radius,
        }
        | CurveGeometry::Arc {
            center,
            normal,
            radius,
            ..
        } => {
            let amp = |n: f64| radius * (1.0 - n * n).max(0.0).sqrt();
            let (ax, ay, az) = (amp(normal.x), amp(normal.y), amp(normal.z));
            points.push(popkern_core::geometry::Point3::new(
                center.x - ax,
                center.y - ay,
                center.z - az,
            ));
            points.push(popkern_core::geometry::Point3::new(
                center.x + ax,
                center.y + ay,
                center.z + az,
            ));
        }
        CurveGeometry::Ellipse {
            center,
            u_basis,
            v_basis,
            ..
        } => {
            let amp = |u: f64, v: f64| (u * u + v * v).sqrt();
            let (ax, ay, az) = (
                amp(u_basis.x, v_basis.x),
                amp(u_basis.y, v_basis.y),
                amp(u_basis.z, v_basis.z),
            );
            points.push(popkern_core::geometry::Point3::new(
                center.x - ax,
                center.y - ay,
                center.z - az,
            ));
            points.push(popkern_core::geometry::Point3::new(
                center.x + ax,
                center.y + ay,
                center.z + az,
            ));
        }
        CurveGeometry::Parabola {
            linear, quadratic, ..
        } => {
            // Interior extreme along axis w at linear·w + 2t·quadratic·w = 0.
            for (l, q) in [
                (linear.x, quadratic.x),
                (linear.y, quadratic.y),
                (linear.z, quadratic.z),
            ] {
                if q.abs() > 1e-30 {
                    let t = -l / (2.0 * q);
                    if t > c.t_start.min(c.t_end) && t < c.t_start.max(c.t_end) {
                        points.push(c.geometry.point_at(t));
                    }
                }
            }
        }
        CurveGeometry::Hyperbola {
            semi_major,
            semi_minor,
            ..
        } => {
            // Extreme along w where sinh(t)·(sM·w) + cosh(t)·(sm·w) = 0.
            for (a, b) in [
                (semi_major.x, semi_minor.x),
                (semi_major.y, semi_minor.y),
                (semi_major.z, semi_minor.z),
            ] {
                if a.abs() > 1e-30 && (b / a).abs() < 1.0 {
                    let t = (-b / a).atanh();
                    if t > c.t_start.min(c.t_end) && t < c.t_start.max(c.t_end) {
                        points.push(c.geometry.point_at(t));
                    }
                }
            }
        }
        CurveGeometry::Nurbs(n) => {
            points.extend(n.control_points.iter().copied());
        }
    }
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

    let trace = crate::parallel::ff_parallel::TraceFlags::from_env();
    if trace.pair {
        eprintln!("  [pairs] {} candidate pairs", candidate_pairs.len());
        for (side, aabbs) in [("A", aabbs_a), ("B", aabbs_b)] {
            for (fh, aabb) in aabbs.iter() {
                eprintln!(
                    "  [pair-aabb] {side} {fh:?} ({:.1},{:.1},{:.1})..({:.1},{:.1},{:.1})",
                    aabb.min.x, aabb.min.y, aabb.min.z, aabb.max.x, aabb.max.y, aabb.max.z
                );
            }
        }
    }
    // Edge-face pierce points, computed BEFORE the face-face solves:
    // each transversal pierce of an edge through the other operand's
    // face is a junction the face-face records must reach, and a
    // known-exact point ON the intersection curves of both face pairs
    // adjacent to the edge. Passing them to the SSI as junction seeds
    // rescues narrow near-tangent lens curves the patch-midpoint
    // seeding cannot find (measured, PK-465: the angled tap's
    // thread-runout crest band crossed the tilted bore wall in
    // 0.05–0.15-long lenses; the wall×crest pair returned EMPTY while
    // all lens endpoints existed as pierce hits, and the missing
    // records cascaded into 20 half-claimed strand declines).
    let edge_face_hits = detect_edge_face_intersections(a, b, aabbs_b, aabbs_a, tolerance);
    let mut pair_seeds: std::collections::HashMap<
        (FaceHandle, FaceHandle),
        Vec<popkern_core::geometry::Point3>,
    > = std::collections::HashMap::new();
    let faces_of_edge =
        |brep: &BRep, eh: popkern_core::topology::EdgeHandle| -> Option<(FaceHandle, FaceHandle)> {
            let he_h = brep.edges.get(eh)?.half_edge;
            let he = brep.half_edges.get(he_h)?;
            let f1 = brep.loops.get(he.loop_)?.face;
            let tw = brep.half_edges.get(he.twin)?;
            let f2 = brep.loops.get(tw.loop_)?.face;
            Some((f1, f2))
        };
    for h in &edge_face_hits {
        let owner = if h.edge_of_a { a } else { b };
        let Some((f1, f2)) = faces_of_edge(owner, h.edge) else {
            continue;
        };
        if popkern_core::trace::trace_enabled("pair-scan") {
            eprintln!(
                "  [hit-map] edge {:?} (of {}) faces ({f1:?},{f2:?}) pierces {:?} at ({:.4},{:.4},{:.4})",
                h.edge,
                if h.edge_of_a { "A" } else { "B" },
                h.face,
                h.point.x,
                h.point.y,
                h.point.z
            );
        }
        for fe in [f1, f2] {
            let key = if h.edge_of_a {
                (fe, h.face)
            } else {
                (h.face, fe)
            };
            let seeds = pair_seeds.entry(key).or_default();
            if !seeds
                .iter()
                .any(|q| q.distance_to(&h.point) < tolerance.linear())
            {
                seeds.push(h.point);
            }
        }
    }

    if popkern_core::trace::trace_enabled("pair-scan") {
        let mut keys: Vec<_> = pair_seeds.iter().collect();
        keys.sort_by_key(|(k, _)| *k);
        for ((fa, fb), pts) in keys {
            eprintln!(
                "  [pair-seeds] {fa:?} x {fb:?}: {} pierce seed(s)",
                pts.len()
            );
        }
    }

    // Resolve BVH index pairs into face-pair tasks, in BVH order. The
    // task list is a deterministic pure function of the two solids, and
    // ff_parallel::process_pairs returns outcomes in task order, so the
    // record sequence below is identical at any thread count — and
    // positionally identical to the retired sequential pair loop.
    let tasks: Vec<crate::parallel::ff_parallel::FacePairTask> = candidate_pairs
        .iter()
        .filter_map(|(id_x, id_y)| {
            let (a_idx, b_idx) = if *id_x < n_a && *id_y >= n_a {
                (*id_x, *id_y - n_a)
            } else if *id_y < n_a && *id_x >= n_a {
                (*id_y, *id_x - n_a)
            } else {
                return None;
            };
            let face_a = aabbs_a[a_idx].0;
            let face_b = aabbs_b[b_idx].0;
            Some(crate::parallel::ff_parallel::FacePairTask {
                face_a,
                face_b,
                seeds: pair_seeds
                    .get(&(face_a, face_b))
                    .cloned()
                    .unwrap_or_default(),
            })
        })
        .collect();

    // First-failure fold in task order: records from pairs before the
    // failing pair are preserved in the diagnostics, exactly as the
    // sequential loop reported them.
    for outcome in crate::parallel::ff_parallel::process_pairs(&tasks, a, b, tolerance, trace) {
        match outcome {
            Ok(mut recs) => records.append(&mut recs),
            Err(failure) => return Err(failure.into_diagnostics(records)),
        }
    }

    // §4.2 branch-point labeling: if two intersection-curve endpoints
    // are coincident within tolerance, mark both as `OnBranchPoint`.
    // The label is informational today (downstream consumers are
    // gated by `#[allow(dead_code)]`) but is required by future split
    // logic to preserve branch topology when subdividing faces.
    label_branch_point_endpoints(&mut records, tolerance);

    if popkern_core::trace::trace_enabled("pair-scan") {
        eprintln!("  [hits] {} edge-face hits", edge_face_hits.len());
        for h in &edge_face_hits {
            eprintln!(
                "    hit at ({:.2},{:.2},{:.2})",
                h.point.x, h.point.y, h.point.z
            );
        }
    }

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

    if popkern_core::trace::trace_enabled("uv-scramble") {
        use popkern_core::geometry::curve::CurveEval as _;
        use popkern_core::geometry::surface::SurfaceEval as _;
        for (i, r) in records.iter().enumerate() {
            for (label, brep, fh, uv) in [
                ("uv1", a, r.face_a, &r.curve.curve_uv1),
                ("uv2", b, r.face_b, &r.curve.curve_uv2),
            ] {
                let Some(surf) = brep.faces.get(fh).and_then(|f| f.surface.as_ref()) else {
                    continue;
                };
                let (c0, c1) = r.curve.curve_3d.parameter_range();
                let (u0, u1) = uv.parameter_range();
                let mut worst = 0.0_f64;
                for k in 0..=8 {
                    let f = k as f64 / 8.0;
                    let p3 = r.curve.curve_3d.point_at(c0 + (c1 - c0) * f);
                    let q = uv.point_at(u0 + (u1 - u0) * f);
                    let ps = surf.geometry.point_at(q.x, q.y);
                    worst = worst.max(ps.distance_to(&p3));
                }
                if worst > tolerance.linear() * 100.0 {
                    eprintln!(
                        "  [uv-scramble] rec {i} ({:?},{:?}) {label} images {worst:.3e} off its 3D curve",
                        r.face_a, r.face_b
                    );
                }
            }
        }
    }

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

#[path = "intersect_consolidate.rs"]
mod intersect_consolidate;
use intersect_consolidate::consolidate_records;

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
