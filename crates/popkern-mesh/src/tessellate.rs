//! Core tessellation: per-face meshing and top-level BRep tessellation.

use std::collections::HashMap;

use popkern_core::geometry::{Point3, SurfaceEval, SurfaceGeometry, Vec3};
use popkern_core::topology::{BRep, FaceHandle};

use crate::adaptive::adaptive_subdivide;

#[path = "tessellate_diagnostics.rs"]
mod tessellate_diagnostics;
#[path = "tessellate_pcurve_range.rs"]
mod tessellate_pcurve_range;
use crate::mesh::{FaceMesh, MeshVertex, TriMesh};
use crate::{MeshError, TessellationOptions};
use tessellate_diagnostics::deepchord_trace;
pub(crate) use tessellate_diagnostics::face_mesh_diagonal;
use tessellate_pcurve_range::pcurve_parameter_range;
#[cfg(test)]
use tessellate_pcurve_range::unwrap_periodic;

/// Tessellate a single BRep face into a triangle mesh.
pub fn tessellate_face(
    brep: &BRep,
    face_handle: FaceHandle,
    options: &TessellationOptions,
) -> Result<FaceMesh, MeshError> {
    tessellate_face_with_samples(brep, face_handle, options, None)
}

/// [`tessellate_face`] with canonical per-edge stations (PK-299) so a
/// legacy-fallback face's trim boundary matches its boundary-first
/// neighbours exactly.
pub(crate) fn tessellate_face_with_samples(
    brep: &BRep,
    face_handle: FaceHandle,
    options: &TessellationOptions,
    edge_samples: Option<
        &std::collections::HashMap<
            popkern_core::topology::EdgeHandle,
            crate::boundary_first::EdgeSamples,
        >,
    >,
) -> Result<FaceMesh, MeshError> {
    let face = brep
        .faces
        .get(face_handle)
        .ok_or_else(|| MeshError::TessellationFailed("invalid face handle".into()))?;

    let surface_on_face = face.surface.as_ref().ok_or(MeshError::NoSurfaceGeometry)?;

    let surface = &surface_on_face.geometry;
    let orientation = surface_on_face.orientation;

    // Get bounded parameter range (handles infinite planes).
    let (mut u_range, mut v_range) = bounded_parameter_range(brep, face_handle, surface);

    // Trim loops drive both the padding decision and the clip below.
    let mut trim_loops = crate::trim::build_trim_loops_with_samples(
        brep,
        face_handle,
        surface,
        options.chord_height,
        edge_samples,
    );

    // FULL-PERIODIC wall face: the outer loop is a seam-pair spine
    // (every outer half-edge's twin lives on this same face — e.g. a
    // full sphere wall surviving a boolean with imprint holes). The
    // spine's pcurves all lie ON the seam iso-line, so the pcurve
    // bbox collapses the u-range to a sliver (measured: cube∪sphere
    // wall meshed u ∈ [0, 0.05] of a 2π chart, mesh volume 0.002)
    // and the degenerate outer polygon suppresses hole trimming. The
    // honest domain is the surface's full natural rectangle; the
    // outer trim is that rectangle, and the imprint holes cut it.
    if outer_loop_is_seam_spine(brep, face_handle) {
        let ((nu0, nu1), (nv0, nv1)) = surface.parameter_range();
        if nu0.is_finite() && nu1.is_finite() && nv0.is_finite() && nv1.is_finite() {
            u_range = (nu0, nu1);
            v_range = (nv0, nv1);
            trim_loops.retain(|l| !l.is_outer);
            // Push the synthetic outer rectangle a hair OUTSIDE the
            // grid domain: a rectangle coinciding exactly with the
            // first/last grid rows classifies every pole-row triangle
            // as Boundary and the exact-overlap clip drops the pole
            // fans (measured: Ø0.16 holes at both sphere poles). The
            // outer rect exists only so the imprint holes cut — it
            // must never itself clip.
            let eu = (nu1 - nu0).max(1.0) * 1.0e-6;
            let ev = (nv1 - nv0).max(1.0) * 1.0e-6;
            let corners = [
                (nu0 - eu, nv0 - ev),
                (nu1 + eu, nv0 - ev),
                (nu1 + eu, nv1 + ev),
                (nu0 - eu, nv1 + ev),
            ];
            trim_loops.push(crate::trim::TrimLoop {
                points: corners.to_vec(),
                points_3d: corners
                    .iter()
                    .map(|&(u, v)| surface.point_at(u, v))
                    .collect(),
                is_outer: true,
                from_pcurves: false,
            });
        }
    }
    let trim_loops = trim_loops;
    // Pad by 2% unless the outer polygon is an axis-aligned rectangle
    // spanning the range (untrimmed faces keep the exact fast grid):
    // trim-polygon iso edges coinciding with grid lines degenerate
    // the boundary clip and the mesh stops half a cell short of the
    // trim (PK-095, measured 0.04-0.13 gaps on fillet strips).
    let rectangular = trim_loops
        .iter()
        .find(|l| l.is_outer)
        .map(|l| {
            let (mut bb_u, mut bb_v) = (
                (f64::INFINITY, f64::NEG_INFINITY),
                (f64::INFINITY, f64::NEG_INFINITY),
            );
            let mut area2 = 0.0;
            for i in 0..l.points.len() {
                let (x1, y1) = l.points[i];
                let (x2, y2) = l.points[(i + 1) % l.points.len()];
                area2 += x1 * y2 - x2 * y1;
                bb_u = (bb_u.0.min(x1), bb_u.1.max(x1));
                bb_v = (bb_v.0.min(y1), bb_v.1.max(y1));
            }
            let bbox_area = (bb_u.1 - bb_u.0) * (bb_v.1 - bb_v.0);
            (0.5 * area2).abs() >= bbox_area * (1.0 - 1.0e-9)
        })
        .unwrap_or(true);
    if !rectangular {
        let pad_u = 0.02 * (u_range.1 - u_range.0).max(1.0e-9);
        let pad_v = 0.02 * (v_range.1 - v_range.0).max(1.0e-9);
        u_range = (u_range.0 - pad_u, u_range.1 + pad_u);
        v_range = (v_range.0 - pad_v, v_range.1 + pad_v);
    }

    // Adaptive subdivision of UV space.
    let cells = adaptive_subdivide(
        surface,
        u_range,
        v_range,
        options.min_subdivisions,
        options.max_subdivisions,
        options.chord_height,
    );

    // Deduplicate vertices by UV coordinate using bit-exact keys.
    let mut vertex_map: HashMap<(u64, u64), usize> = HashMap::new();
    let mut vertices: Vec<MeshVertex> = Vec::new();
    let mut triangles: Vec<[usize; 3]> = Vec::new();

    for cell in &cells {
        // Four corners of the quad cell.
        let corners = [
            (cell.u_min, cell.v_min),
            (cell.u_max, cell.v_min),
            (cell.u_min, cell.v_max),
            (cell.u_max, cell.v_max),
        ];

        let mut indices = [0_usize; 4];
        for (i, &(u, v)) in corners.iter().enumerate() {
            let key = (u.to_bits(), v.to_bits());
            indices[i] = *vertex_map.entry(key).or_insert_with(|| {
                let pos = surface.point_at(u, v);
                let mut normal = surface.normal_at(u, v);
                // Normalize the normal vector.
                let len = normal.length();
                if len > brep.tolerance.numerical_floor() {
                    normal = Vec3::new(normal.x / len, normal.y / len, normal.z / len);
                }
                // Flip normal if face orientation is reversed.
                if !orientation {
                    normal = -normal;
                }
                let idx = vertices.len();
                vertices.push(MeshVertex {
                    position: pos,
                    normal,
                    uv: (u, v),
                });
                idx
            });
        }

        // Split quad into 2 triangles along the shorter diagonal.
        // Indices: 0=bottom-left, 1=bottom-right, 2=top-left, 3=top-right
        let diag_02_13 = vertices[indices[0]]
            .position
            .distance_to(&vertices[indices[3]].position);
        let diag_01_23 = vertices[indices[1]]
            .position
            .distance_to(&vertices[indices[2]].position);

        if diag_02_13 <= diag_01_23 {
            // Split along 0-3 diagonal.
            triangles.push([indices[0], indices[1], indices[3]]);
            triangles.push([indices[0], indices[3], indices[2]]);
        } else {
            // Split along 1-2 diagonal.
            triangles.push([indices[0], indices[1], indices[2]]);
            triangles.push([indices[1], indices[3], indices[2]]);
        }
    }

    // Trim: use pcurve-based adaptive polylines for precise boundary clipping.
    if !trim_loops.is_empty() {
        let poly_area = crate::trim::trim_loop_area(&trim_loops);
        let uv_area = (u_range.1 - u_range.0) * (v_range.1 - v_range.0);
        // Skip trimming for degenerate boundaries (sphere seams, etc.)
        if poly_area > uv_area * 0.01 {
            // All construction pipelines (constructors, sweep, boolean,
            // fillet, shell, draft) now generate pcurves via
            // recompute_all_pcurves. The full clipping pipeline runs
            // unconditionally. The 3D validation in apply_trim handles
            // any residual UV imprecision.
            let (trimmed_verts, trimmed_tris) = crate::trim_clip::apply_trim(
                &vertices,
                &triangles,
                &trim_loops,
                surface,
                orientation,
            );
            vertices = trimmed_verts;
            triangles = trimmed_tris;
        }
    }

    // When face orientation is reversed, flip triangle winding to match.
    // The adaptive subdivision always produces CCW winding from the surface
    // geometry's normal direction. If orientation=false, the effective face
    // normal is flipped, so triangles must be rewound to remain front-facing
    // under CCW back-face culling.
    if !orientation {
        for tri in &mut triangles {
            tri.swap(1, 2);
        }
    }

    Ok(FaceMesh {
        face: face_handle,
        vertices,
        triangles,
    })
}

/// Compute a bounded parameter range for a face's surface for tessellation.
///
/// Uses the pcurve-based trim loop bounding box when pcurves are available,
/// giving exact UV bounds that match the trim boundary. Falls back to
/// sample-and-project with 5% padding for faces without pcurves, where
/// the padding ensures coverage for degenerate topologies (e.g., sphere
/// seam edges where all vertices project to u=0).
///
/// On periodic surfaces (cylinder, torus), if the raw UV range spans
/// more than half the period, the boundary likely crosses the seam.
/// In that case, shift all UV samples by half-period and recompute
/// min/max to find the shorter arc that doesn't straddle the seam.
/// True when every half-edge of the face's outer loop has its twin on
/// the SAME face — the outer boundary is a periodic seam spine (two
/// meridian copies bounding a full-period chart), not a real border.
pub(crate) fn outer_loop_is_seam_spine(brep: &BRep, face_handle: FaceHandle) -> bool {
    let Some(face) = brep.faces.get(face_handle) else {
        return false;
    };
    let Some(loop_) = brep.loops.get(face.outer_loop) else {
        return false;
    };
    let start = loop_.half_edge;
    if !brep.half_edges.contains_key(start) {
        return false;
    }
    let mut he = start;
    for _ in 0..10_000 {
        let data = &brep.half_edges[he];
        let twin_on_face = brep
            .half_edges
            .get(data.twin)
            .and_then(|twin| brep.loops.get(twin.loop_))
            .is_some_and(|l| l.face == face_handle);
        if !twin_on_face {
            return false;
        }
        he = data.next;
        if he == start {
            return true;
        }
    }
    false
}

fn bounded_parameter_range(
    brep: &BRep,
    face_handle: FaceHandle,
    surface: &SurfaceGeometry,
) -> ((f64, f64), (f64, f64)) {
    let positions = brep.face_boundary_sample_points(face_handle, 16);
    let mut result =
        popkern_core::geometry::uv_projection::tight_bounded_parameter_range(&positions, surface);

    // If pcurve-based bounds are available, intersect with the fallback.
    if let Some(pc_range) = pcurve_parameter_range(brep, face_handle) {
        let u = (pc_range.0.0.max(result.0.0), pc_range.0.1.min(result.0.1));
        let v = (pc_range.1.0.max(result.1.0), pc_range.1.1.min(result.1.1));
        if u.0 < u.1 && v.0 < v.1 {
            result = (u, v);
        }
    }

    result
}

/// Honest stage ceiling for [`tessellate_brep`] progress totals
/// (PK-276): edge sampling, up to 8 encroachment rounds, face pass,
/// boundary sew, orientation, and assembly+zipper. Early-converging
/// rounds emit fewer events, never more.
pub(crate) const TESSELLATE_STAGE_CEILING: usize = 13;

/// Tessellate all faces of a BRep into a combined [`TriMesh`].
///
/// `tracker` threads PK-276 interactive control: a cancellation poll
/// per encroachment round and per epilogue stage
/// (`MeshError::Cancelled`, nothing produced), a poll at every
/// rayon-parallel per-face closure entry (a fired token short-circuits
/// remaining faces), and one progress event per round/stage.
pub(crate) fn tessellate_brep(
    brep: &BRep,
    options: &TessellationOptions,
    tracker: &mut popkern_core::progress::StageTracker<'_, '_>,
) -> Result<TriMesh, MeshError> {
    use std::time::Instant;
    macro_rules! stage_gate {
        () => {
            if tracker.cancelled() {
                return Err(MeshError::Cancelled);
            }
        };
    }
    let cancel = tracker.token();
    // Canonical per-edge samples: adjacent faces share bit-identical
    // boundary vertices (PK-095 closure).
    stage_gate!();
    let t0 = Instant::now();
    let mut edge_samples = crate::boundary_first::sample_all_edges(brep, options.chord_height);
    tracker.finish_stage("tessellate.edge_samples", t0, edge_samples.len());
    // Encroachment refinement (Ruppert-style): a face whose Delaunay
    // coverage fails because boundary points crowd a ring segment
    // (engraved pockets, thin annuli) requests SPLITS of the
    // canonical edge spans involved. Splits refine the SHARED sample
    // sets and all faces re-run, so watertightness is preserved by
    // construction. Bounded rounds; convergence is guaranteed by the
    // no-new-fraction check.
    const ENCROACHMENT_ROUNDS: usize = 8;
    let bf_trace = std::env::var_os("PK_BF_TRACE").is_some();
    // Incremental rounds (PK-277): round 0 checks every face; later
    // rounds only faces incident to an edge whose canonical samples
    // were refined last round. A face's requests are a pure function
    // of its own boundary samples, and every skip reason below is
    // monotone (fractions only grow, per-edge caps only fill), so an
    // unchanged-boundary face can only reproduce requests that were
    // already inserted or are permanently deduped.
    let mut active: Option<std::collections::HashSet<FaceHandle>> = None;
    for round in 0..ENCROACHMENT_ROUNDS {
        // §17 per-round poll (the per-face closures also short-circuit
        // on a fired token, so a cancelled round finishes fast).
        stage_gate!();
        let round_started = Instant::now();
        // Slotmap key order, filtered — request order matches the
        // retired all-faces loop with the no-op faces removed.
        let faces: Vec<FaceHandle> = brep
            .faces
            .keys()
            .filter(|&fh| brep.faces[fh].surface.is_some())
            .filter(|fh| active.as_ref().is_none_or(|set| set.contains(fh)))
            .collect();
        let requests = crate::face_pass::collect_split_requests(
            brep,
            &faces,
            &edge_samples,
            options,
            bf_trace,
            cancel,
        );
        stage_gate!();
        if requests.is_empty() {
            break;
        }
        let mut split_edges: std::collections::BTreeSet<popkern_core::topology::EdgeHandle> =
            std::collections::BTreeSet::new();
        let mut inserted = 0usize;
        for (eh, fraction) in requests {
            let Some(samples) = edge_samples.get_mut(&eh) else {
                continue;
            };
            if samples.fractions.len() >= 96 {
                continue; // refinement cap per edge
            }
            let idx = samples.fractions.partition_point(|&f| f < fraction);
            let near = |f: f64| (f - fraction).abs() < 1e-4;
            if (idx > 0 && near(samples.fractions[idx - 1]))
                || (idx < samples.fractions.len() && near(samples.fractions[idx]))
            {
                continue; // already refined here
            }
            let Some(curve) = brep.edges.get(eh).and_then(|e| e.curve.as_ref()) else {
                continue;
            };
            use popkern_core::geometry::curve::CurveEval;
            let ceval: &dyn CurveEval = &curve.geometry;
            let point = ceval.point_at(curve.t_start + (curve.t_end - curve.t_start) * fraction);
            samples.fractions.insert(idx, fraction);
            samples.points.insert(idx, point);
            split_edges.insert(eh);
            inserted += 1;
        }
        if bf_trace && inserted > 0 {
            eprintln!("  [bf] encroachment round {round}: split {inserted} edge spans");
        }
        tracker.finish_stage("tessellate.encroachment_round", round_started, inserted);
        if inserted == 0 {
            break;
        }
        active = Some(crate::face_pass::faces_incident_to_edges(
            brep,
            &split_edges,
        ));
    }
    stage_gate!();
    let t0 = Instant::now();
    let mut face_meshes = Vec::new();
    let mut legacy_faces = Vec::new();
    let faces: Vec<FaceHandle> = brep
        .faces
        .keys()
        .filter(|&fh| brep.faces[fh].surface.is_some())
        .collect();
    let outcomes = crate::face_pass::tessellate_face_pass(
        brep,
        &faces,
        &edge_samples,
        options,
        bf_trace,
        cancel,
    );
    // A token fired mid-pass leaves short-circuited outcomes; surface
    // the cancellation rather than a fabricated per-face error.
    stage_gate!();
    for outcome in outcomes {
        match outcome {
            crate::face_pass::FacePassOutcome::BoundaryFirst(fm) => face_meshes.push(fm),
            crate::face_pass::FacePassOutcome::Legacy { seam_spine, mesh } => {
                // Seam-spine wall faces (full-periodic charts) keep
                // their exact chart sampling: their boundary IS the
                // seam iso-line plus pcurve-sampled hole rims, and the
                // generous projection gate (20% of the face diagonal —
                // the WHOLE model for a full sphere wall) drags rim
                // vertices onto unrelated edge polylines (measured:
                // cube∪sphere wall area inflated 7.93 → 8.95 by
                // sewing flaps).
                if !seam_spine {
                    legacy_faces.push(face_meshes.len());
                }
                face_meshes.push(mesh?);
            }
        }
    }
    tracker.finish_stage("tessellate.face_pass", t0, face_meshes.len());
    // The projection-based sewing stopgap applies ONLY to legacy-path
    // faces: boundary-first faces carry exact canonical boundary
    // points, and re-projecting them can snap a point onto a
    // DIFFERENT nearby incident edge (measured: floor arc points
    // pulled onto the base-edge line, reopening the very seams the
    // canonical sampling closed).
    stage_gate!();
    let t0 = Instant::now();
    deepchord_trace(brep, &face_meshes, options, "pre-sew");
    crate::boundary_sew::sew_boundaries_to_edge_polylines(brep, &mut face_meshes, &legacy_faces);
    deepchord_trace(brep, &face_meshes, options, "post-sew");
    tracker.finish_stage("tessellate.sew_boundaries", t0, legacy_faces.len());
    stage_gate!();
    let t0 = Instant::now();
    orient_face_meshes_coherently(brep, &mut face_meshes);
    tracker.finish_stage("tessellate.orient", t0, face_meshes.len());
    stage_gate!();
    let t0 = Instant::now();
    let mut mesh = TriMesh::from_face_meshes(face_meshes);
    // Tolerance-scoped T-junction zippering (PK-208): differently
    // sampled boundary chains along RECORDED-tolerance edges are
    // subdivided to share their points, so the weld closes. Unrecorded
    // gaps are left open — see crack_zipper for why that scoping is
    // the whole design.
    mesh.zipped_t_junctions =
        crate::crack_zipper::zip_recorded_tolerance_cracks(brep, &mut mesh, options.chord_height);
    tracker.finish_stage("tessellate.assemble_zipper", t0, mesh.triangles.len());
    Ok(mesh)
}

/// Global winding coherence: no per-face orientation witness is
/// trustworthy on its own — `SurfaceOnFace.orientation` is measurably
/// dishonest on extrude bottom caps (the breached-bore void floor
/// meshed inside-out: −170516 wedge where +170516 is exact, 10% of
/// the solid volume), while chamfer cones carry honest flags but
/// loops wound AGAINST the CCW-about-outward convention (calibrating
/// facing from loop winding inverted them instead). Adjacent faces of
/// a closed solid must traverse a shared edge in OPPOSITE directions.
///
/// The traversal-direction constraint is measured against the BREP
/// TOPOLOGY, not against bit-identical mesh vertices: for each edge
/// with two half-edges on different faces, each face mesh's boundary
/// segments near the edge's canonical polyline vote with their
/// signed arc-length advance. (The earlier bit-identity segment
/// matching only paired the handful of boundary vertices two faces
/// happened to sample at the same parameter — 3-vote constraints on
/// 60-segment edges — and its first-vote-wins propagation over
/// HashMap-ordered adjacency was nondeterministic: a 2×3×1 prism's
/// mesh wedge volume measured {+6, +2, −2, +4} across identical
/// runs.) Flips propagate across the shared-edge graph per component,
/// then each component's global side is fixed by the sign of its
/// wedge volume. Faces with no measurable shared-edge constraint keep
/// their flag-derived facing, as before.
/// Position-weld key of a mesh vertex: the bit patterns of (x, y, z).
type PositionBitKey = (u64, u64, u64);

/// Winding votes for a face-mesh pair (mesh indices, lower first):
/// summed matched arc length voting that the two boundaries traverse
/// shared edges in the (same, opposite) direction.
type PairWindingVotes = std::collections::HashMap<(usize, usize), (f64, f64)>;

fn orient_face_meshes_coherently(brep: &BRep, face_meshes: &mut [crate::mesh::FaceMesh]) {
    use std::collections::HashMap;

    let n = face_meshes.len();
    if std::env::var_os("PK_MESH_ORIENT_TRACE").is_some() {
        for (fi, fm) in face_meshes.iter().enumerate() {
            let class = brep
                .faces
                .get(fm.face)
                .and_then(|f| f.surface.as_ref())
                .map(|s| match &s.geometry {
                    SurfaceGeometry::Plane { origin, .. } => format!("Plane z={:.3}", origin.z),
                    SurfaceGeometry::Cylinder { radius, .. } => format!("Cyl r={radius:.3}"),
                    SurfaceGeometry::Cone { apex, .. } => format!("Cone apex_z={:.3}", apex.z),
                    _ => "other".to_string(),
                });
            let mut v6 = 0.0;
            for t in &fm.triangles {
                let (a, b, c) = (
                    fm.vertices[t[0]].position,
                    fm.vertices[t[1]].position,
                    fm.vertices[t[2]].position,
                );
                v6 += a.x * (b.y * c.z - b.z * c.y) - a.y * (b.x * c.z - b.z * c.x)
                    + a.z * (b.x * c.y - b.y * c.x);
            }
            eprintln!(
                "  [orient] mesh {fi} = {:?} {class:?} tris={} verts={} wedge={:.5}",
                fm.face,
                fm.triangles.len(),
                fm.vertices.len(),
                v6 / 6.0
            );
        }
    }
    let mesh_index_of: HashMap<popkern_core::topology::FaceHandle, usize> = face_meshes
        .iter()
        .enumerate()
        .map(|(fi, fm)| (fm.face, fi))
        .collect();

    // Boundary directed segments per face mesh: triangle edges used
    // exactly once (position-welded so duplicated boundary vertices
    // at the same position don't hide the match).
    let boundary_segments: Vec<Vec<(Point3, Point3)>> = face_meshes
        .iter()
        .map(|fm| {
            let key =
                |p: &Point3| -> PositionBitKey { (p.x.to_bits(), p.y.to_bits(), p.z.to_bits()) };
            let mut counts: HashMap<(PositionBitKey, PositionBitKey), i64> = HashMap::new();
            for t in &fm.triangles {
                for (a, b) in [(t[0], t[1]), (t[1], t[2]), (t[2], t[0])] {
                    let (ka, kb) = (key(&fm.vertices[a].position), key(&fm.vertices[b].position));
                    if ka == kb {
                        continue;
                    }
                    if ka < kb {
                        *counts.entry((ka, kb)).or_default() += 1;
                    } else {
                        *counts.entry((kb, ka)).or_default() -= 1;
                    }
                }
            }
            let from_bits = |k: (u64, u64, u64)| -> Point3 {
                Point3::new(
                    f64::from_bits(k.0),
                    f64::from_bits(k.1),
                    f64::from_bits(k.2),
                )
            };
            counts
                .into_iter()
                .filter(|&(_, c)| c != 0)
                .map(|((ka, kb), c)| {
                    if c > 0 {
                        (from_bits(ka), from_bits(kb))
                    } else {
                        (from_bits(kb), from_bits(ka))
                    }
                })
                .collect()
        })
        .collect();

    // Signed arc-length advance of a face's boundary along a
    // canonical edge polyline. Positive = traverses the edge in
    // increasing parameter, negative = decreasing; magnitude = the
    // summed length of matched segments (the vote's weight).
    let advance_along = |segs: &[(Point3, Point3)], polyline: &[(f64, Point3)]| -> f64 {
        let mut length = 0.0_f64;
        for w in polyline.windows(2) {
            length += w[0].1.distance_to(&w[1].1);
        }
        if length <= 0.0 {
            return 0.0;
        }
        let tol = (length * 1.0e-3).max(1.0e-9);
        // Arc-length parameter of the nearest polyline point, or None
        // if the point is farther than tol from every segment.
        let project = |p: &Point3| -> Option<f64> {
            let mut best: Option<(f64, f64)> = None; // (dist, s)
            let mut s0 = 0.0_f64;
            for w in polyline.windows(2) {
                let (a, b) = (w[0].1, w[1].1);
                let ab = b - a;
                let len = ab.length();
                if len <= 0.0 {
                    continue;
                }
                let t = ((*p - a).dot(&ab) / (len * len)).clamp(0.0, 1.0);
                let q = a + ab * t;
                let d = p.distance_to(&q);
                if best.is_none_or(|(bd, _)| d < bd) {
                    best = Some((d, s0 + len * t));
                }
                s0 += len;
            }
            best.and_then(|(d, s)| (d <= tol).then_some(s))
        };
        // Closed polylines (rim circles): a raw arc-length delta
        // telescopes to zero around the loop — the wraparound
        // segment's jump (−circumference + step) cancels every other
        // step and the cap↔wall constraint silently vanishes
        // (measured: a cylinder cap stayed inside-out, mesh wedge
        // volume exactly ⅔ of 2πr²h). Wrap each per-segment delta to
        // the short way around.
        let closed = polyline
            .first()
            .zip(polyline.last())
            .is_some_and(|(a, b)| a.1.distance_to(&b.1) <= tol);
        let mut advance = 0.0_f64;
        for (p, q) in segs {
            if let (Some(sp), Some(sq)) = (project(p), project(q)) {
                let mut d = sq - sp;
                if closed {
                    d -= length * (d / length).round();
                }
                advance += d;
            }
        }
        advance
    };

    // Face-pair parity votes over shared BRep edges, accumulated by
    // weight (summed matched boundary length). Same-signed advance on
    // both sides means one of the pair is inside-out.
    let mut votes: PairWindingVotes = HashMap::new(); // (same, opposite)
    for (_, edge) in brep.edges.iter() {
        let Some(curve) = edge.curve.as_ref() else {
            continue;
        };
        if !(curve.t_start.is_finite() && curve.t_end.is_finite()) {
            continue;
        }
        let he = edge.half_edge;
        let (Some(h1), Some(h2)) = (
            brep.half_edges.get(he),
            brep.half_edges
                .get(he)
                .and_then(|h| brep.half_edges.get(h.twin)),
        ) else {
            continue;
        };
        let (Some(l1), Some(l2)) = (brep.loops.get(h1.loop_), brep.loops.get(h2.loop_)) else {
            continue;
        };
        if l1.face == l2.face {
            continue; // seam edge: no cross-face constraint
        }
        let (Some(&f1), Some(&f2)) = (mesh_index_of.get(&l1.face), mesh_index_of.get(&l2.face))
        else {
            continue;
        };
        const ORIENT_SAMPLES: usize = 33;
        let eval: &dyn popkern_core::geometry::curve::CurveEval = &curve.geometry;
        let polyline: Vec<(f64, Point3)> = (0..ORIENT_SAMPLES)
            .map(|i| {
                let t = curve.t_start
                    + (curve.t_end - curve.t_start) * i as f64 / (ORIENT_SAMPLES - 1) as f64;
                (t, eval.point_at(t))
            })
            .collect();
        let d1 = advance_along(&boundary_segments[f1], &polyline);
        let d2 = advance_along(&boundary_segments[f2], &polyline);
        // Require a non-trivial matched length on BOTH sides before
        // trusting the vote (a lone straggler segment is sew noise).
        let mut edge_len = 0.0;
        for w in polyline.windows(2) {
            edge_len += w[0].1.distance_to(&w[1].1);
        }
        let significant = edge_len * 1.0e-3;
        if d1.abs() < significant || d2.abs() < significant {
            if std::env::var_os("PK_MESH_ORIENT_TRACE").is_some() {
                eprintln!(
                    "  [orient] pair ({f1},{f2}) advance below floor: d1={d1:.6} d2={d2:.6} \
                     edge_len={edge_len:.4}"
                );
            }
            continue;
        }
        let weight = d1.abs().min(d2.abs());
        let entry = votes.entry((f1.min(f2), f1.max(f2))).or_default();
        if d1 * d2 > 0.0 {
            entry.0 += weight;
        } else {
            entry.1 += weight;
        }
    }

    // Majority-resolved adjacency, sorted for deterministic traversal.
    let mut adj: Vec<Vec<(usize, bool)>> = vec![Vec::new(); n]; // (other, same_direction)
    let mut pairs: Vec<_> = votes.iter().collect();
    pairs.sort_by_key(|(k, _)| **k);
    for (&(f1, f2), &(same, opposite)) in pairs {
        if same == opposite {
            // Tied evidence — no trustworthy constraint. Leave the
            // pair unlinked rather than guess.
            if std::env::var_os("PK_MESH_ORIENT_TRACE").is_some() {
                eprintln!("  [orient] pair ({f1},{f2}) tied votes {same}:{opposite} — unlinked");
            }
            continue;
        }
        let same_dir = same > opposite;
        if same.min(opposite) > 0.0 && std::env::var_os("PK_MESH_ORIENT_TRACE").is_some() {
            eprintln!(
                "  [orient] pair ({f1},{f2}) conflicted votes same={same:.3} opposite={opposite:.3} — majority {}",
                if same_dir { "same" } else { "opposite" }
            );
        }
        adj[f1].push((f2, same_dir));
        adj[f2].push((f1, same_dir));
    }
    // BFS parity propagation per component (deterministic: seeds in
    // index order, adjacency sorted above).
    let mut flip: Vec<Option<bool>> = vec![None; n];
    for seed in 0..n {
        if flip[seed].is_some() {
            continue;
        }
        flip[seed] = Some(false);
        let mut queue = std::collections::VecDeque::from([seed]);
        let mut component = vec![seed];
        while let Some(f) = queue.pop_front() {
            let f_flip = flip[f].unwrap_or(false);
            for &(g, same_dir) in &adj[f] {
                // Opposite traversal (same_dir=false) is CONSISTENT:
                // g keeps f's parity. Same-direction traversal means
                // one of the two is inside-out relative to the other.
                let want = if same_dir { !f_flip } else { f_flip };
                match flip[g] {
                    None => {
                        flip[g] = Some(want);
                        component.push(g);
                        queue.push_back(g);
                    }
                    Some(assigned) => {
                        if assigned != want && std::env::var_os("PK_MESH_ORIENT_TRACE").is_some() {
                            eprintln!(
                                "  [orient] parity cycle conflict at faces ({f},{g}) — keeping first assignment"
                            );
                        }
                    }
                }
            }
        }
        // Component-global side: an OUTER boundary component encloses
        // positive volume (outward facing); a VOID boundary — a
        // component enclosed inside another component of the same body
        // (ray-crossing parity, position-only, so it needs no facing
        // decisions) — encloses negative volume, facing into the void.
        // The former unconditional-positive rule meshed every hollow
        // body's inner skin inside-out (see winding_side).
        let mut vol6 = 0.0;
        for &f in &component {
            let s = if flip[f] == Some(true) { -1.0 } else { 1.0 };
            let fm = &face_meshes[f];
            for t in &fm.triangles {
                let (a, b, c) = (
                    fm.vertices[t[0]].position,
                    fm.vertices[t[1]].position,
                    fm.vertices[t[2]].position,
                );
                vol6 += s
                    * (a.x * (b.y * c.z - b.z * c.y) - a.y * (b.x * c.z - b.z * c.x)
                        + a.z * (b.x * c.y - b.y * c.x));
            }
        }
        let enclosed = crate::winding_side::component_is_enclosed(&component, face_meshes);
        let wrong_side = if enclosed { vol6 > 0.0 } else { vol6 < 0.0 };
        if std::env::var_os("PK_MESH_ORIENT_TRACE").is_some() {
            eprintln!(
                "  [orient] component {component:?} vol6 {vol6:.4} enclosed={enclosed} \
                 wrong_side={wrong_side}"
            );
        }
        if wrong_side {
            for &f in &component {
                flip[f] = flip[f].map(|x| !x);
            }
        }
    }
    for (fi, fm) in face_meshes.iter_mut().enumerate() {
        if flip[fi] == Some(true) {
            for t in &mut fm.triangles {
                t.swap(1, 2);
            }
            for v in &mut fm.vertices {
                v.normal = -v.normal;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tessellate;
    use popkern_core::tolerance::Tolerance;

    #[test]
    fn tessellate_unit_cube_correct_counts() {
        let cube = BRep::unit_cube(Tolerance::default()).unwrap();
        let mesh = tessellate(&cube, &TessellationOptions::default()).unwrap();
        // Boundary-first tessellation: counts depend on the interior
        // grid + boundary sampling; the structural contract is one
        // mesh per face, consistent maps, and full area coverage.
        assert!(!mesh.triangles.is_empty());
        assert_eq!(mesh.face_map.len(), mesh.triangles.len());
        assert_eq!(mesh.face_meshes.len(), 6);
        let mut area = 0.0;
        for t in &mesh.triangles {
            let a = mesh.vertices[t[0]].position;
            let b = mesh.vertices[t[1]].position;
            let c = mesh.vertices[t[2]].position;
            let ab = b - a;
            let ac = c - a;
            let cx = popkern_core::geometry::Vec3::new(
                ab.y * ac.z - ab.z * ac.y,
                ab.z * ac.x - ab.x * ac.z,
                ab.x * ac.y - ab.y * ac.x,
            );
            area += 0.5 * cx.length();
        }
        assert!((area - 6.0).abs() < 1e-9, "cube surface area, got {area}");
    }

    #[test]
    fn tessellate_unit_cube_min_subdiv_1() {
        let cube = BRep::unit_cube(Tolerance::default()).unwrap();
        let opts = TessellationOptions {
            min_subdivisions: 1,
            ..Default::default()
        };
        let mesh = tessellate(&cube, &opts).unwrap();
        // Coarsest setting still covers every face.
        assert!(!mesh.triangles.is_empty());
        assert_eq!(mesh.face_map.len(), mesh.triangles.len());
        assert_eq!(mesh.face_meshes.len(), 6);
    }

    #[test]
    fn tessellate_sphere_adaptive() {
        let sphere = BRep::sphere(1.0, Tolerance::default()).unwrap();
        let mesh = tessellate(&sphere, &TessellationOptions::default()).unwrap();
        assert!(mesh.triangles.len() > 12);
    }

    #[test]
    fn tessellate_cube_normals_axis_aligned() {
        let cube = BRep::unit_cube(Tolerance::default()).unwrap();
        let mesh = tessellate(&cube, &TessellationOptions::default()).unwrap();
        for v in &mesh.vertices {
            let abs_sum = v.normal.x.abs() + v.normal.y.abs() + v.normal.z.abs();
            assert!(
                (abs_sum - 1.0).abs() < 0.1,
                "Normal should be roughly unit axis-aligned, got ({}, {}, {})",
                v.normal.x,
                v.normal.y,
                v.normal.z
            );
        }
    }

    #[test]
    fn tessellate_sphere_vertices_on_surface() {
        let sphere = BRep::sphere(1.0, Tolerance::default()).unwrap();
        let mesh = tessellate(&sphere, &TessellationOptions::default()).unwrap();
        for v in &mesh.vertices {
            let dist = (v.position.x.powi(2) + v.position.y.powi(2) + v.position.z.powi(2)).sqrt();
            assert!((dist - 1.0).abs() < 1e-6);
        }
    }

    #[test]
    fn tighter_chord_height_more_triangles() {
        let sphere = BRep::sphere(1.0, Tolerance::default()).unwrap();
        let coarse = tessellate(
            &sphere,
            &TessellationOptions {
                chord_height: 0.1,
                ..Default::default()
            },
        )
        .unwrap();
        let fine = tessellate(
            &sphere,
            &TessellationOptions {
                chord_height: 0.001,
                ..Default::default()
            },
        )
        .unwrap();
        assert!(fine.triangles.len() > coarse.triangles.len());
    }

    #[test]
    fn all_indices_in_bounds() {
        let cube = BRep::unit_cube(Tolerance::default()).unwrap();
        let mesh = tessellate(&cube, &TessellationOptions::default()).unwrap();
        for tri in &mesh.triangles {
            for &idx in tri {
                assert!(idx < mesh.vertices.len());
            }
        }
    }
}

#[cfg(test)]
#[path = "tessellate_edge_tests.rs"]
mod tessellate_edge_tests;
