//! Projection-based boundary sewing for legacy-path face meshes.
//! Extracted from `tessellate.rs` (file-size cap, PK-276); applies
//! ONLY to legacy-path faces — boundary-first faces carry exact
//! canonical boundary points, and re-projecting them can snap a point
//! onto a DIFFERENT nearby incident edge (see `tessellate_brep`).

use popkern_core::geometry::Point3;
use popkern_core::topology::BRep;

use crate::tessellate::face_mesh_diagonal;

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
pub(crate) fn sew_boundaries_to_edge_polylines(
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
                let (eh, pcurve, twin_same_face) = {
                    let data = &brep.half_edges[he];
                    // Boundary sides have a dead twin: not same-face.
                    let same = brep
                        .half_edges
                        .get(data.twin)
                        .and_then(|twin| brep.loops.get(twin.loop_))
                        .map(|l| l.face == fm.face)
                        .unwrap_or(false);
                    he = data.next;
                    (data.edge, data.pcurve.clone(), same)
                };
                // SEAM edges (twin on this same face) have no neighbor
                // face to weld against — their boundary points already
                // come from this face's own pcurves, and re-snapping
                // them onto the canonical polyline nudged a
                // pole-adjacent seam vertex into a folded hairline fan
                // sliver (measured on the primitive sphere once its
                // seam carried true arc geometry; the chord era
                // excluded the seam by accident via the on-surface
                // gate).
                if twin_same_face {
                    if he == start {
                        break;
                    }
                    continue;
                }
                let agree = (|| {
                    // THIS face's half-edge pcurve (the twin maps into
                    // the neighboring face's surface).
                    let pc = pcurve.as_ref()?;
                    let curve = brep.edges.get(eh)?.curve.as_ref()?;
                    use popkern_core::geometry::SurfaceEval;
                    let _ = curve;
                    // The gate honors the edge's RECORDED tolerance: a
                    // construction that accepted a bounded boundary
                    // inexactness (the setback fan's fitted rim rows,
                    // 1.5e-3..2.3e-3 — a kinked deviation field a
                    // degree-9 row cannot represent below ~9e-4)
                    // stamps it on the edge; excluding those edges
                    // here left every rim boundary vertex unsewn and
                    // the blunt-pentagon fan cracked open (236 welded
                    // open edges). An UNRECORDED disagreement keeps
                    // the tight gate — the sphere-seam pole chord this
                    // gate exists for stays excluded.
                    let recorded = brep.edges.get(eh).and_then(|e| e.tolerance).unwrap_or(0.0);
                    let on_surface_gate = on_surface_gate.max(1.5 * recorded);
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
