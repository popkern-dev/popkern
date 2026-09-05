//! Core tessellation: per-face meshing and top-level BRep tessellation.

use std::collections::HashMap;

use popkern_core::geometry::{SurfaceEval, SurfaceGeometry, Vec3};
use popkern_core::topology::{BRep, FaceHandle};

use crate::adaptive::adaptive_subdivide;

#[path = "tessellate_diagnostics.rs"]
mod tessellate_diagnostics;
#[path = "tessellate_orientation.rs"]
mod tessellate_orientation;
#[path = "tessellate_pcurve_range.rs"]
mod tessellate_pcurve_range;
use crate::mesh::{FaceMesh, MeshVertex, TriMesh};
use crate::{MeshError, TessellationOptions};
use tessellate_diagnostics::deepchord_trace;
pub(crate) use tessellate_diagnostics::face_mesh_diagonal;
use tessellate_orientation::orient_face_meshes_coherently;
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
/// uniform-grid-fallback face's trim boundary matches its boundary-first
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
        // Full-period axes are NOT padded: the chart wraps onto
        // itself there, so padding double-covers a seam slice (the
        // grid meshes the same surface strip twice — measured +3%
        // volume on an offset cone once the full-period grid landed).
        let (pp_u, pp_v) = surface.periodicity();
        let full_u = pp_u.is_some_and(|p| p > 0.0 && u_range.1 - u_range.0 >= 0.95 * p);
        let full_v = pp_v.is_some_and(|p| p > 0.0 && v_range.1 - v_range.0 >= 0.95 * p);
        if !full_u {
            let pad_u = 0.02 * (u_range.1 - u_range.0).max(1.0e-9);
            u_range = (u_range.0 - pad_u, u_range.1 + pad_u);
        }
        if !full_v {
            let pad_v = 0.02 * (v_range.1 - v_range.0).max(1.0e-9);
            v_range = (v_range.0 - pad_v, v_range.1 + pad_v);
        }
        // The pad must not cross the surface's NATURAL parameter
        // bounds: past them the chart is either undefined or a mirror
        // image (a cone's v < 0 is the second nappe — the padded grid
        // meshed a tiny beyond-apex disk whose ring read as open mesh
        // edges on every cone).
        let ((nu0, nu1), (nv0, nv1)) = surface.parameter_range();
        if nu0.is_finite() {
            u_range.0 = u_range.0.max(nu0);
        }
        if nu1.is_finite() {
            u_range.1 = u_range.1.min(nu1);
        }
        if nv0.is_finite() {
            v_range.0 = v_range.0.max(nv0);
        }
        if nv1.is_finite() {
            v_range.1 = v_range.1.min(nv1);
        }
    }
    if popkern_core::trace::trace_enabled("trim") {
        eprintln!(
            "  [trim] face {face_handle:?} rectangular={rectangular} grid u={u_range:?} v={v_range:?}"
        );
        for l in &trim_loops {
            let (mut bu, mut bv) = (
                (f64::INFINITY, f64::NEG_INFINITY),
                (f64::INFINITY, f64::NEG_INFINITY),
            );
            for &(u, v) in &l.points {
                bu = (bu.0.min(u), bu.1.max(u));
                bv = (bv.0.min(v), bv.1.max(v));
            }
            eprintln!(
                "  [trim]   loop outer={} from_pcurves={} n={} bbox u={bu:?} v={bv:?}",
                l.is_outer,
                l.from_pcurves,
                l.points.len()
            );
            if l.points.len() <= 8 {
                for &(u, v) in &l.points {
                    eprintln!("  [trim]     pt ({u:.4},{v:.4})");
                }
            }
        }
    }

    let max_subdivisions = crate::adaptive::effective_max_subdivisions(surface, options);
    let cells = adaptive_subdivide(
        surface,
        u_range,
        v_range,
        options.min_subdivisions,
        max_subdivisions,
        options.chord_height,
    );

    // Deduplicate vertices by UV coordinate using bit-exact keys.
    let mut vertex_map: HashMap<(u64, u64), usize> = HashMap::new();
    let mut vertices: Vec<MeshVertex> = Vec::new();
    let mut triangles: Vec<[usize; 3]> = Vec::new();

    // Conforming closure for the adaptive quadtree (PK-333): a coarse
    // cell adjacent to finer neighbours carries their corner points as
    // HANGING NODES on its edges — triangulating each quad
    // independently left an open crack along every refinement
    // transition (measured: an 80-edge open ring at the cone's
    // near-apex refinement boundary, invisible to the crack zipper,
    // which is scoped to BRep-edge welds). Iso-line registries of all
    // cell corner coordinates let each cell include its neighbours'
    // points; children share bit-exact parent arithmetic, so lookup is
    // exact. Cells with no hanging nodes keep the original
    // shorter-diagonal split bit-for-bit.
    let sortable = |f: f64| -> u64 {
        let b = f.to_bits();
        if b >> 63 == 1 { !b } else { b ^ (1u64 << 63) }
    };
    let mut corners_on_u: HashMap<u64, std::collections::BTreeMap<u64, f64>> = HashMap::new();
    let mut corners_on_v: HashMap<u64, std::collections::BTreeMap<u64, f64>> = HashMap::new();
    for cell in &cells {
        for (u, v) in [
            (cell.u_min, cell.v_min),
            (cell.u_max, cell.v_min),
            (cell.u_min, cell.v_max),
            (cell.u_max, cell.v_max),
        ] {
            corners_on_u
                .entry(u.to_bits())
                .or_default()
                .insert(sortable(v), v);
            corners_on_v
                .entry(v.to_bits())
                .or_default()
                .insert(sortable(u), u);
        }
    }
    let hanging = |registry: &HashMap<u64, std::collections::BTreeMap<u64, f64>>,
                   line: f64,
                   lo: f64,
                   hi: f64|
     -> Vec<f64> {
        use std::ops::Bound::Excluded;
        // A zero-extent cell side (a seam-touching sliver trim) has
        // lo == hi; an equal excluded range panics in BTreeMap.
        if sortable(lo) >= sortable(hi) {
            return Vec::new();
        }
        registry
            .get(&line.to_bits())
            .map(|m| {
                m.range((Excluded(sortable(lo)), Excluded(sortable(hi))))
                    .map(|(_, &x)| x)
                    .collect()
            })
            .unwrap_or_default()
    };

    fn intern_vertex(
        u: f64,
        v: f64,
        surface: &SurfaceGeometry,
        orientation: bool,
        floor: f64,
        vertex_map: &mut HashMap<(u64, u64), usize>,
        vertices: &mut Vec<MeshVertex>,
    ) -> usize {
        let key = (u.to_bits(), v.to_bits());
        *vertex_map.entry(key).or_insert_with(|| {
            let pos = surface.point_at(u, v);
            let mut normal = surface.normal_at(u, v);
            let len = normal.length();
            if len > floor {
                normal = Vec3::new(normal.x / len, normal.y / len, normal.z / len);
            }
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
        })
    }
    let floor = brep.tolerance.numerical_floor();

    for cell in &cells {
        let bottom = hanging(&corners_on_v, cell.v_min, cell.u_min, cell.u_max);
        let right = hanging(&corners_on_u, cell.u_max, cell.v_min, cell.v_max);
        let top = hanging(&corners_on_v, cell.v_max, cell.u_min, cell.u_max);
        let left = hanging(&corners_on_u, cell.u_min, cell.v_min, cell.v_max);

        if bottom.is_empty() && right.is_empty() && top.is_empty() && left.is_empty() {
            // Four corners of the quad cell.
            let corners = [
                (cell.u_min, cell.v_min),
                (cell.u_max, cell.v_min),
                (cell.u_min, cell.v_max),
                (cell.u_max, cell.v_max),
            ];
            let mut indices = [0_usize; 4];
            for (i, &(u, v)) in corners.iter().enumerate() {
                indices[i] = intern_vertex(
                    u,
                    v,
                    surface,
                    orientation,
                    floor,
                    &mut vertex_map,
                    &mut vertices,
                );
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
            continue;
        }

        // Transition cell: fan-triangulate around the cell centre over
        // the counter-clockwise boundary ring including hanging nodes.
        let mut ring: Vec<(f64, f64)> = Vec::new();
        ring.push((cell.u_min, cell.v_min));
        ring.extend(bottom.iter().map(|&u| (u, cell.v_min)));
        ring.push((cell.u_max, cell.v_min));
        ring.extend(right.iter().map(|&v| (cell.u_max, v)));
        ring.push((cell.u_max, cell.v_max));
        ring.extend(top.iter().rev().map(|&u| (u, cell.v_max)));
        ring.push((cell.u_min, cell.v_max));
        ring.extend(left.iter().rev().map(|&v| (cell.u_min, v)));

        let u_mid = 0.5 * (cell.u_min + cell.u_max);
        let v_mid = 0.5 * (cell.v_min + cell.v_max);
        let center = intern_vertex(
            u_mid,
            v_mid,
            surface,
            orientation,
            floor,
            &mut vertex_map,
            &mut vertices,
        );
        let ring_idx: Vec<usize> = ring
            .iter()
            .map(|&(u, v)| {
                intern_vertex(
                    u,
                    v,
                    surface,
                    orientation,
                    floor,
                    &mut vertex_map,
                    &mut vertices,
                )
            })
            .collect();
        for i in 0..ring_idx.len() {
            let a = ring_idx[i];
            let b = ring_idx[(i + 1) % ring_idx.len()];
            if a != b {
                triangles.push([center, a, b]);
            }
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
        // FULL-PERIOD axes take the pcurve range verbatim: the
        // position-sampled fallback tops out below the period (16
        // boundary samples of a rim circle end at 15/16·2π), and
        // intersecting truncated the grid — every revolution
        // primitive's lateral face meshed with a missing seam slice
        // (~1.3% of the angle, visible as a crack on cones and
        // cylinders in any render or STL).
        let (period_u, period_v) = surface.periodicity();
        if let Some(p) = period_u.filter(|p| *p > 0.0)
            && pc_range.0.1 - pc_range.0.0 >= 0.95 * p
        {
            result.0 = pc_range.0;
        }
        if let Some(p) = period_v.filter(|p| *p > 0.0)
            && pc_range.1.1 - pc_range.1.0 >= 0.95 * p
        {
            result.1 = pc_range.1;
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
    let t0 = popkern_core::time::Stopwatch::start();
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
    let bf_trace = popkern_core::trace::trace_enabled("boundary-first");
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
        let round_started = popkern_core::time::Stopwatch::start();
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
            if samples.fractions.len() >= 512 {
                // Refinement cap per edge — matches the adaptive
                // sampler's own cap (96 starved the xfine cut-curve
                // splits and the coverage check bounced the whole
                // wall to the uniform-grid path).
                continue;
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
    let t0 = popkern_core::time::Stopwatch::start();
    let mut face_meshes = Vec::new();
    let mut sew_candidate_faces = Vec::new();
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
            crate::face_pass::FacePassOutcome::UniformGrid { seam_spine, mesh } => {
                // Seam-spine wall faces (full-periodic charts) keep
                // their exact chart sampling: their boundary IS the
                // seam iso-line plus pcurve-sampled hole rims, and the
                // generous projection gate (20% of the face diagonal —
                // the WHOLE model for a full sphere wall) drags rim
                // vertices onto unrelated edge polylines (measured:
                // cube∪sphere wall area inflated 7.93 → 8.95 by
                // sewing flaps).
                if !seam_spine {
                    sew_candidate_faces.push(face_meshes.len());
                }
                face_meshes.push(mesh?);
            }
        }
    }
    tracker.finish_stage("tessellate.face_pass", t0, face_meshes.len());
    // The projection-based sewing stopgap applies ONLY to uniform-grid-path
    // faces: boundary-first faces carry exact canonical boundary
    // points, and re-projecting them can snap a point onto a
    // DIFFERENT nearby incident edge (measured: floor arc points
    // pulled onto the base-edge line, reopening the very seams the
    // canonical sampling closed).
    stage_gate!();
    let t0 = popkern_core::time::Stopwatch::start();
    deepchord_trace(brep, &face_meshes, options, "pre-sew");
    crate::boundary_sew::sew_boundaries_to_edge_polylines(
        brep,
        &mut face_meshes,
        &sew_candidate_faces,
    );
    deepchord_trace(brep, &face_meshes, options, "post-sew");
    tracker.finish_stage("tessellate.sew_boundaries", t0, sew_candidate_faces.len());
    stage_gate!();
    let t0 = popkern_core::time::Stopwatch::start();
    orient_face_meshes_coherently(brep, &mut face_meshes);
    tracker.finish_stage("tessellate.orient", t0, face_meshes.len());
    stage_gate!();
    let t0 = popkern_core::time::Stopwatch::start();
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
pub(crate) type PositionBitKey = (u64, u64, u64);

/// Winding votes for a face-mesh pair (mesh indices, lower first):
/// summed matched arc length voting that the two boundaries traverse
/// shared edges in the (same, opposite) direction.
pub(crate) type PairWindingVotes = std::collections::HashMap<(usize, usize), (f64, f64)>;

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
