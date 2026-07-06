//! Core tessellation: per-face meshing and top-level BRep tessellation.

use std::collections::HashMap;

use popkern_core::geometry::{Point3, SurfaceEval, SurfaceGeometry, Vec3};
use popkern_core::topology::{BRep, FaceHandle};

use crate::adaptive::adaptive_subdivide;
use crate::mesh::{FaceMesh, MeshVertex, TriMesh};
use crate::{MeshError, TessellationOptions};

/// Tessellate a single BRep face into a triangle mesh.
pub fn tessellate_face(
    brep: &BRep,
    face_handle: FaceHandle,
    options: &TessellationOptions,
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
    let trim_loops =
        crate::trim::build_trim_loops(brep, face_handle, surface, options.chord_height);
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

/// Compute UV parameter range from the outer loop's pcurves.
///
/// Returns `Some((u_range, v_range))` if every half-edge in the outer loop
/// has a pcurve; `None` otherwise, signaling the caller to fall back to
/// 3D projection.
///
/// On periodic surfaces, UV values are unwrapped relative to the first
/// sample to prevent seam crossings from inflating the range to the
/// full period.
fn pcurve_parameter_range(
    brep: &BRep,
    face_handle: FaceHandle,
) -> Option<((f64, f64), (f64, f64))> {
    let face = brep.faces.get(face_handle)?;
    let surface = face.surface.as_ref()?;
    let (period_u, period_v) = surface.geometry.periodicity();

    let loop_ = brep.loops.get(face.outer_loop)?;
    let start_he = loop_.half_edge;
    if !brep.half_edges.contains_key(start_he) {
        return None;
    }

    let samples = 16_usize;
    let mut all_u = Vec::new();
    let mut all_v = Vec::new();

    let mut current = start_he;
    loop {
        let he = &brep.half_edges[current];
        let pcurve = he.pcurve.as_ref()?;
        let (t0, t1) = pcurve.parameter_range();
        for i in 0..=samples {
            let t = t0 + (t1 - t0) * (i as f64 / samples as f64);
            let pt = pcurve.point_at(t);
            all_u.push(pt.x);
            all_v.push(pt.y);
        }
        current = he.next;
        if current == start_he {
            break;
        }
    }

    if all_u.is_empty() {
        return None;
    }

    // Unwrap periodic values relative to the first sample.
    unwrap_periodic(&mut all_u, period_u);
    unwrap_periodic(&mut all_v, period_v);

    let u_min = all_u.iter().cloned().fold(f64::INFINITY, f64::min);
    let u_max = all_u.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let v_min = all_v.iter().cloned().fold(f64::INFINITY, f64::min);
    let v_max = all_v.iter().cloned().fold(f64::NEG_INFINITY, f64::max);

    if u_min.is_finite() && u_max.is_finite() && v_min.is_finite() && v_max.is_finite() {
        Some(((u_min, u_max), (v_min, v_max)))
    } else {
        None
    }
}

/// Unwrap a sequence of values on a periodic domain so consecutive
/// values don't jump across the seam. Adjusts values in-place by
/// adding/subtracting the period to stay within half-period of the
/// previous value.
fn unwrap_periodic(values: &mut [f64], period: Option<f64>) {
    let Some(p) = period else { return };
    if p <= 0.0 || values.len() < 2 {
        return;
    }
    let half = p * 0.5;
    for i in 1..values.len() {
        let diff = values[i] - values[i - 1];
        if diff > half {
            values[i] -= p;
        } else if diff < -half {
            values[i] += p;
        }
    }
}

/// Tessellate all faces of a BRep into a combined [`TriMesh`].
pub(crate) fn tessellate_brep(
    brep: &BRep,
    options: &TessellationOptions,
) -> Result<TriMesh, MeshError> {
    // Canonical per-edge samples: adjacent faces share bit-identical
    // boundary vertices (PK-095 closure).
    let mut edge_samples = crate::boundary_first::sample_all_edges(brep, options.chord_height);
    // Encroachment refinement (Ruppert-style): a face whose Delaunay
    // coverage fails because boundary points crowd a ring segment
    // (engraved pockets, thin annuli) requests SPLITS of the
    // canonical edge spans involved. Splits refine the SHARED sample
    // sets and all faces re-run, so watertightness is preserved by
    // construction. Bounded rounds; convergence is guaranteed by the
    // no-new-fraction check.
    const ENCROACHMENT_ROUNDS: usize = 8;
    for round in 0..ENCROACHMENT_ROUNDS {
        let mut requests: Vec<crate::boundary_first::SplitRequest> = Vec::new();
        for face_handle in brep.faces.keys() {
            if brep.faces[face_handle].surface.is_none() {
                continue;
            }
            crate::boundary_first::tessellate_face_boundary_first(
                brep,
                face_handle,
                &edge_samples,
                options,
                &mut requests,
            );
        }
        if requests.is_empty() {
            break;
        }
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
            inserted += 1;
        }
        if std::env::var_os("PK_BF_TRACE").is_some() && inserted > 0 {
            eprintln!("  [bf] encroachment round {round}: split {inserted} edge spans");
        }
        if inserted == 0 {
            break;
        }
    }
    let mut face_meshes = Vec::new();
    let mut legacy_faces = Vec::new();
    for face_handle in brep.faces.keys() {
        if brep.faces[face_handle].surface.is_none() {
            continue;
        }
        let mut sink = Vec::new();
        if let Some(fm) = crate::boundary_first::tessellate_face_boundary_first(
            brep,
            face_handle,
            &edge_samples,
            options,
            &mut sink,
        ) {
            face_meshes.push(fm);
            continue;
        }
        if std::env::var_os("PK_BF_TRACE").is_some() {
            let kind = brep
                .faces
                .get(face_handle)
                .and_then(|f| f.surface.as_ref())
                .map(|sof| match &sof.geometry {
                    popkern_core::geometry::SurfaceGeometry::Plane { .. } => "Plane",
                    popkern_core::geometry::SurfaceGeometry::Cylinder { .. } => "Cylinder",
                    popkern_core::geometry::SurfaceGeometry::Nurbs(_) => "Nurbs",
                    popkern_core::geometry::SurfaceGeometry::Gregory(_) => "Gregory",
                    _ => "Other",
                })
                .unwrap_or("None");
            eprintln!("  [bf] face {face_handle:?} ({kind}) fell back");
        }
        legacy_faces.push(face_meshes.len());
        face_meshes.push(tessellate_face(brep, face_handle, options)?);
    }
    // The projection-based sewing stopgap applies ONLY to legacy-path
    // faces: boundary-first faces carry exact canonical boundary
    // points, and re-projecting them can snap a point onto a
    // DIFFERENT nearby incident edge (measured: floor arc points
    // pulled onto the base-edge line, reopening the very seams the
    // canonical sampling closed).
    sew_boundaries_to_edge_polylines(brep, &mut face_meshes, &legacy_faces);
    Ok(TriMesh::from_face_meshes(face_meshes))
}

/// Boundary sewing (PK-095): face meshes are sampled independently,
/// so their boundary vertices land at slightly different positions
/// along the shared edges — visible cracks between every pair of
/// adjacent faces. Each edge gets ONE canonical polyline (sampled
/// from its 3D curve, endpoints snapped to the topology vertices),
/// and every face-mesh vertex within tolerance of an incident edge's
/// polyline is projected exactly onto it. Both sides then lie on the
/// same 3D polyline: the seam is visually sealed (residual
/// T-junctions between differently-parameterized samples are the
/// recorded follow-up for exact topological closure).
fn sew_boundaries_to_edge_polylines(
    brep: &BRep,
    face_meshes: &mut [crate::mesh::FaceMesh],
    legacy_faces: &[usize],
) {
    use std::collections::HashMap;

    const SEW_SAMPLES: usize = 129;
    // Canonical (t, point) polylines per edge — the polyline seeds a
    // projection whose result is evaluated ON the curve (projecting
    // onto chords would pull on-surface vertices off the geometry by
    // the chord error).
    let mut polylines: HashMap<popkern_core::topology::EdgeHandle, Vec<(f64, Point3)>> =
        HashMap::new();
    for (eh, edge) in brep.edges.iter() {
        let Some(curve) = edge.curve.as_ref() else {
            continue;
        };
        if !(curve.t_start.is_finite() && curve.t_end.is_finite()) {
            continue;
        }
        let eval: &dyn popkern_core::geometry::curve::CurveEval = &curve.geometry;
        let pts: Vec<(f64, Point3)> = (0..SEW_SAMPLES)
            .map(|i| {
                let t = curve.t_start
                    + (curve.t_end - curve.t_start) * i as f64 / (SEW_SAMPLES - 1) as f64;
                (t, eval.point_at(t))
            })
            .collect();
        polylines.insert(eh, pts);
    }
    // Per face: incident edge set from its loops.
    for &fi in legacy_faces {
        let fm = &mut face_meshes[fi];
        let Some(face) = brep.faces.get(fm.face) else {
            continue;
        };
        let mut edges: Vec<popkern_core::topology::EdgeHandle> = Vec::new();
        let mut loops = vec![face.outer_loop];
        loops.extend(face.inner_loops.iter().copied());
        for lh in loops {
            let Some(loop_) = brep.loops.get(lh) else {
                continue;
            };
            let start = loop_.half_edge;
            let mut he = start;
            for _ in 0..10_000 {
                edges.push(brep.half_edges[he].edge);
                he = brep.half_edges[he].next;
                if he == start {
                    break;
                }
            }
        }
        // Only the face mesh's own BOUNDARY vertices are sewn (its
        // interior is a grid cell away from any edge). Boundary =
        // endpoints of position-welded edges used by exactly one
        // triangle. The gate is generous (20% of the face diagonal):
        // boundary chords at coarse tessellation sit several percent
        // off the neighbor's finer chords of the same 3D curve.
        let tol = 0.2 * face_mesh_diagonal(fm).max(1.0e-6);
        // Only trust edges whose 3D curve actually lies on THIS
        // face's surface (a sphere seam records a pole-to-pole chord
        // Line — sewing onto it pulls vertices inside the solid).
        // Agreement is checked through the half-edge's own pcurve:
        // surface(pcurve(t)) must match curve(t) at three params.
        let surface = match face.surface.as_ref() {
            Some(sf) => sf.geometry.clone(),
            None => continue,
        };
        let on_surface_gate = 1.0e-4 * face_mesh_diagonal(fm).max(1.0);
        let mut usable: Vec<popkern_core::topology::EdgeHandle> = Vec::new();
        for lh in [face.outer_loop]
            .into_iter()
            .chain(face.inner_loops.iter().copied())
        {
            let Some(loop_) = brep.loops.get(lh) else {
                continue;
            };
            let start = loop_.half_edge;
            let mut he = start;
            for _ in 0..10_000 {
                let (eh, pcurve) = {
                    let data = &brep.half_edges[he];
                    he = data.next;
                    (data.edge, data.pcurve.clone())
                };
                let agree = (|| {
                    // THIS face's half-edge pcurve (the twin maps into
                    // the neighboring face's surface).
                    let pc = pcurve.as_ref()?;
                    let curve = brep.edges.get(eh)?.curve.as_ref()?;
                    use popkern_core::geometry::SurfaceEval;
                    use popkern_core::geometry::curve::CurveEval;
                    let _ = curve;
                    // Param-free agreement: the surface image of the
                    // pcurve at three params must lie ON the edge's
                    // canonical polyline (nearest-distance — pcurve
                    // and curve parameterizations need not correspond).
                    let pts = polylines.get(&eh)?;
                    for f in [0.25_f64, 0.5, 0.75] {
                        let tp = pc.t_start + (pc.t_end - pc.t_start) * f;
                        let uv = pc.point_at(tp);
                        let p = surface.point_at(uv.x, uv.y);
                        let d = pts
                            .windows(2)
                            .map(|w| {
                                project_to_segment_param(p, w[0].1, w[1].1)
                                    .1
                                    .distance_to(&p)
                            })
                            .fold(f64::INFINITY, f64::min);
                        if d > on_surface_gate {
                            return None;
                        }
                    }
                    Some(())
                })()
                .is_some();
                if agree {
                    usable.push(eh);
                }
                if he == start {
                    break;
                }
            }
        }
        let edges = usable;
        let boundary = face_mesh_boundary_vertices(fm);
        for (vi, v) in fm.vertices.iter_mut().enumerate() {
            if !boundary.contains(&vi) {
                continue;
            }
            let mut best: Option<(f64, popkern_core::topology::EdgeHandle, f64)> = None;
            for eh in &edges {
                let Some(pts) = polylines.get(eh) else {
                    continue;
                };
                for w in pts.windows(2) {
                    let (frac, proj) = project_to_segment_param(v.position, w[0].1, w[1].1);
                    let d = proj.distance_to(&v.position);
                    if best.map(|(bd, _, _)| d < bd).unwrap_or(true) {
                        let t = w[0].0 + frac * (w[1].0 - w[0].0);
                        best = Some((d, *eh, t));
                    }
                }
            }
            if let Some((d, eh, t)) = best
                && d < tol
                && let Some(curve) = brep.edges.get(eh).and_then(|e| e.curve.as_ref())
            {
                let eval: &dyn popkern_core::geometry::curve::CurveEval = &curve.geometry;
                // Refine the chord-seeded parameter on the true curve
                // (ternary search on distance²) so vertices already on
                // the curve project back to themselves.
                let dt = (curve.t_end - curve.t_start) / (SEW_SAMPLES - 1) as f64;
                let (mut lo, mut hi) = (
                    (t - dt).max(curve.t_start.min(curve.t_end)),
                    (t + dt).min(curve.t_start.max(curve.t_end)),
                );
                for _ in 0..60 {
                    let m1 = lo + (hi - lo) / 3.0;
                    let m2 = hi - (hi - lo) / 3.0;
                    if eval.point_at(m1).distance_to(&v.position)
                        < eval.point_at(m2).distance_to(&v.position)
                    {
                        hi = m2;
                    } else {
                        lo = m1;
                    }
                }
                let t = 0.5 * (lo + hi);
                let on_curve = eval.point_at(t);
                v.position = on_curve;
            }
        }
    }
}

fn face_mesh_diagonal(fm: &crate::mesh::FaceMesh) -> f64 {
    let mut min = [f64::INFINITY; 3];
    let mut max = [f64::NEG_INFINITY; 3];
    for v in &fm.vertices {
        for (k, c) in [v.position.x, v.position.y, v.position.z]
            .into_iter()
            .enumerate()
        {
            min[k] = min[k].min(c);
            max[k] = max[k].max(c);
        }
    }
    ((max[0] - min[0]).powi(2) + (max[1] - min[1]).powi(2) + (max[2] - min[2]).powi(2)).sqrt()
}

/// Vertex indices on the face mesh's boundary: endpoints of edges
/// used by exactly one triangle after position welding (clipping can
/// duplicate vertices at identical positions).
fn face_mesh_boundary_vertices(fm: &crate::mesh::FaceMesh) -> std::collections::HashSet<usize> {
    use std::collections::HashMap;
    let mut canon: HashMap<(i64, i64, i64), usize> = HashMap::new();
    let mut remap = Vec::with_capacity(fm.vertices.len());
    for v in &fm.vertices {
        let k = (
            (v.position.x * 1.0e9).round() as i64,
            (v.position.y * 1.0e9).round() as i64,
            (v.position.z * 1.0e9).round() as i64,
        );
        let next = canon.len();
        remap.push(*canon.entry(k).or_insert(next));
    }
    let mut counts: HashMap<(usize, usize), usize> = HashMap::new();
    for t in &fm.triangles {
        let w = [remap[t[0]], remap[t[1]], remap[t[2]]];
        if w[0] == w[1] || w[1] == w[2] || w[0] == w[2] {
            continue;
        }
        for (a, b) in [(w[0], w[1]), (w[1], w[2]), (w[2], w[0])] {
            *counts.entry((a.min(b), a.max(b))).or_insert(0) += 1;
        }
    }
    let mut boundary_canon = std::collections::HashSet::new();
    for (&(a, b), &c) in &counts {
        if c == 1 {
            boundary_canon.insert(a);
            boundary_canon.insert(b);
        }
    }
    (0..fm.vertices.len())
        .filter(|&i| boundary_canon.contains(&remap[i]))
        .collect()
}

fn project_to_segment_param(p: Point3, a: Point3, b: Point3) -> (f64, Point3) {
    let ab = b - a;
    let len2 = ab.dot(&ab);
    if len2 < 1.0e-30 {
        return (0.0, a);
    }
    let t = ((p - a).dot(&ab) / len2).clamp(0.0, 1.0);
    (
        t,
        Point3::new(a.x + ab.x * t, a.y + ab.y * t, a.z + ab.z * t),
    )
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
