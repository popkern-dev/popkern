//! Vertex position recomputation for inner skin.

use std::collections::{BTreeMap, HashMap, HashSet};

use crate::error::ShellError;
use crate::skin_map::SkinMap;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle, VertexHandle};

/// Recompute all vertex positions in the inner skin.
///
/// For each vertex, finds the adjacent faces, maps them back to the outer
/// BRep, retrieves the corresponding offset surfaces, and computes the
/// new vertex position as the intersection of those surfaces.
///
/// For vertices adjacent to 3+ planes, uses Cramer's rule for an exact
/// solution. Otherwise falls back to a Newton-iteration approach from
/// an initial guess.
///
/// # Errors
///
/// Returns [`ShellError::Topology`] if adjacency lookups fail.
pub fn recompute_vertices(
    inner: &mut BRep,
    offset_surfaces: &HashMap<FaceHandle, SurfaceGeometry>,
    skin_map: &SkinMap,
    outer: &BRep,
    tolerance: &Tolerance,
) -> Result<(), ShellError> {
    recompute_vertices_with_gap_corners(
        inner,
        offset_surfaces,
        skin_map,
        outer,
        tolerance,
        &BTreeMap::new(),
    )
}

/// [`recompute_vertices`] with gap-corner handling (PK-325/PK-326): a
/// vertex scheduled for arc-join corner surgery takes its RETAINED
/// sector's target position instead of an intersection solve — the
/// other sectors split off to their own offset targets later, so an
/// intersection here would place the sharp point the arc join removes.
/// Single-face sectors project the original vertex exactly onto that
/// face's offset surface; edge-pair sectors (mixed corners, E1) take
/// the nearest point on the two offset surfaces' intersection so the
/// retained resolvable edge recomputes against its final endpoint.
pub(crate) fn recompute_vertices_with_gap_corners(
    inner: &mut BRep,
    offset_surfaces: &HashMap<FaceHandle, SurfaceGeometry>,
    skin_map: &SkinMap,
    outer: &BRep,
    tolerance: &Tolerance,
    gap_corners: &BTreeMap<VertexHandle, crate::join_classification::CornerPlanEntry>,
) -> Result<(), ShellError> {
    let reverse_face = skin_map.reverse_face_map();
    let reverse_vertex = skin_map.reverse_vertex_map();

    let inner_vertex_handles: Vec<_> = inner.vertices.keys().collect();

    // Phase 1 (parallel, PK-329): each vertex's solve reads the skin
    // topology, the offset surfaces, and the OUTER body immutably, and
    // writes only its own position — so the per-vertex solves dispatch
    // in parallel with an ordered collect and the sequential apply
    // below is positionally identical to the single loop (first
    // failing vertex's error included).
    let solve_one =
        |inner_vh: VertexHandle| -> Result<Option<popkern_core::geometry::Point3>, ShellError> {
            if let Some(entry) = gap_corners.get(&inner_vh) {
                let p_orig = reverse_vertex
                    .get(&inner_vh)
                    .and_then(|ov| outer.vertices.get(*ov).map(|v| v.position))
                    .unwrap_or_else(|| inner.vertices[inner_vh].position);
                let target = match entry.retained {
                    crate::join_classification::CornerTarget::Face(fh) => offset_surfaces
                        .get(&fh)
                        .and_then(|s| crate::curved_recompute::project_to_surface(s, &p_orig)),
                    crate::join_classification::CornerTarget::EdgePair(f1, f2) => {
                        match (offset_surfaces.get(&f1), offset_surfaces.get(&f2)) {
                            (Some(s1), Some(s2)) => {
                                crate::curved_recompute::solve_vertex_on_surfaces(
                                    &[s1, s2],
                                    p_orig,
                                    tolerance,
                                )
                            }
                            _ => None,
                        }
                    }
                };
                if let Some(q) = target {
                    return Ok(Some(q));
                }
                // Projection ambiguity: fall through to the general solve;
                // the corner surgery repositions the vertex again anyway.
            }
            // Find adjacent faces by walking outgoing half-edges.
            let adjacent_outer_faces = find_adjacent_outer_faces(inner, inner_vh, &reverse_face);

            // Collect offset surfaces for adjacent faces, keyed by outer
            // face handle so curved solves can order deterministically
            // (the adjacency set iterates in hash order — a PK-277-class
            // hazard for any order-sensitive iteration).
            let mut planes: Vec<(Vec3, f64)> = Vec::new();
            let mut surfaces: Vec<&SurfaceGeometry> = Vec::new();
            let mut keyed: Vec<(FaceHandle, &SurfaceGeometry)> = Vec::new();

            for outer_fh in &adjacent_outer_faces {
                if let Some(surf) = offset_surfaces.get(outer_fh) {
                    surfaces.push(surf);
                    keyed.push((*outer_fh, surf));
                    if let SurfaceGeometry::Plane { origin, normal, .. } = surf {
                        let d = normal.dot(&Vec3::new(origin.x, origin.y, origin.z));
                        planes.push((*normal, d));
                    }
                }
            }

            // Opening-rim vertex (PK-307): faces adjacent in the OUTER
            // body but absent from the inner skin were removed — the rim
            // vertex must stay ON each removed face's ORIGINAL surface so
            // the opening remains flush. Without this constraint a rim
            // corner sees only two kept offsets and lands wherever the
            // averaging iteration stops (measured ~9e-7 off the exact
            // corner — a ShellGap right at the tolerance boundary); with
            // it, three planes solve exactly.
            if let Some(&outer_vh) = reverse_vertex.get(&inner_vh) {
                // "Kept" means the face contributed an offset constraint. A
                // face may be PRESENT in the skin yet unclassified (midplane
                // outer skins keep removed faces so assembly can discover
                // the rims) — its original surface is the rim constraint
                // either way.
                let kept: HashSet<FaceHandle> = adjacent_outer_faces
                    .iter()
                    .copied()
                    .filter(|fh| offset_surfaces.contains_key(fh))
                    .collect();
                for outer_fh in faces_at_vertex(outer, outer_vh) {
                    if kept.contains(&outer_fh) {
                        continue;
                    }
                    if let Some(sof) = outer.faces.get(outer_fh).and_then(|f| f.surface.as_ref()) {
                        surfaces.push(&sof.geometry);
                        keyed.push((outer_fh, &sof.geometry));
                        if let SurfaceGeometry::Plane { origin, normal, .. } = &sof.geometry {
                            let d = normal.dot(&Vec3::new(origin.x, origin.y, origin.z));
                            planes.push((*normal, d));
                        }
                    }
                }
            }

            if surfaces.is_empty() {
                return Ok(None);
            }

            let has_curved = surfaces
                .iter()
                .any(|s| !matches!(s, SurfaceGeometry::Plane { .. }));

            // PK-322: a vertex sitting at a cone face's parametric apex
            // tracks the OFFSET cone's apex exactly — projection cannot
            // find it (every direction from the axis is equidistant).
            if has_curved && let Some(&outer_vh) = reverse_vertex.get(&inner_vh) {
                let orig_pos = outer.vertices[outer_vh].position;
                let mut apex_target = None;
                let mut sorted_faces = adjacent_outer_faces.clone();
                sorted_faces.sort();
                for outer_fh in &sorted_faces {
                    if let (
                        Some(SurfaceGeometry::Cone { apex: a0, .. }),
                        Some(SurfaceGeometry::Cone { apex: a1, .. }),
                    ) = (
                        outer
                            .faces
                            .get(*outer_fh)
                            .and_then(|f| f.surface.as_ref())
                            .map(|s| &s.geometry),
                        offset_surfaces.get(outer_fh),
                    ) && orig_pos.distance_to(a0) <= tolerance.linear()
                    {
                        apex_target = Some(*a1);
                        break;
                    }
                }
                if let Some(pos) = apex_target {
                    return Ok(Some(pos));
                }
            }

            // Solve the new position. Pure-plane valence-3+ keeps its exact
            // Cramer arm; any curved participant routes through sequential
            // exact projection (deterministic: sorted by face handle,
            // seeded from the original vertex for branch selection); the
            // plane-averaging Newton covers the plane-only
            // under-determined cases.
            let new_pos = if !has_curved && planes.len() == 3 {
                solve_3_planes(&planes[0], &planes[1], &planes[2])
            } else if !has_curved && planes.len() > 3 {
                // Valence>3: least-squares over ALL incident plane
                // constraints in sorted-handle order (design §5 — the
                // first-3-of-a-hash-set pick both dropped constraints and
                // was nondeterministic). Exact when the system is
                // consistent; the best point otherwise, with validation as
                // the arbiter.
                solve_planes_least_squares(&planes)
            } else if has_curved {
                let outer_vh = reverse_vertex
                    .get(&inner_vh)
                    .ok_or_else(|| ShellError::Topology("vertex not in reverse map".into()))?;
                let seed = outer.vertices[*outer_vh].position;
                keyed.sort_by_key(|(fh, _)| *fh);
                let ordered: Vec<&SurfaceGeometry> = keyed.iter().map(|(_, s)| *s).collect();
                crate::curved_recompute::solve_vertex_on_surfaces(&ordered, seed, tolerance)
            } else {
                // Under-determined plane sets (free-rim vertices of open
                // sheets: 1 or 2 incident planes) solve EXACTLY as the
                // minimal displacement of the original vertex onto the
                // constraint set — the previous averaging iteration stopped
                // ~7e-7 off the intersection line, which the quarter-scale
                // ladder tolerance rightly rejected downstream (measured:
                // every open-box rim edge kept its stale clone).
                let outer_vh = reverse_vertex
                    .get(&inner_vh)
                    .ok_or_else(|| ShellError::Topology("vertex not in reverse map".into()))?;
                let seed = outer.vertices[*outer_vh].position;
                solve_underdetermined_planes(&planes, seed)
            };

            Ok(new_pos)
        };

    let solved = crate::parallel_dispatch::dispatch_ordered(&inner_vertex_handles, false, |&vh| {
        solve_one(vh)
    });

    // Phase 2 (sequential): apply position writes in vertex order.
    for (inner_vh, outcome) in inner_vertex_handles.iter().zip(solved) {
        if let Some(pos) = outcome? {
            inner.vertices[*inner_vh].position = pos;
        }
    }

    Ok(())
}

/// Faces adjacent to `vertex` in `brep`, by walking the outgoing fan.
fn faces_at_vertex(brep: &BRep, vertex: VertexHandle) -> Vec<FaceHandle> {
    let mut faces = HashSet::new();
    let start_he = brep.vertices[vertex].outgoing;
    if !brep.half_edges.contains_key(start_he) {
        return Vec::new();
    }
    let mut current = start_he;
    loop {
        let he = &brep.half_edges[current];
        if let Some(lp) = brep.loops.get(he.loop_) {
            faces.insert(lp.face);
        }
        let twin = he.twin;
        if !brep.half_edges.contains_key(twin) {
            break;
        }
        let next_outgoing = brep.half_edges[twin].next;
        if next_outgoing == start_he || !brep.half_edges.contains_key(next_outgoing) {
            break;
        }
        current = next_outgoing;
    }
    // Sorted handles: the fan set iterates in hash order, and every
    // downstream constraint list (planes, curved projection order)
    // inherits this order — a PK-277-class determinism hazard.
    let mut faces: Vec<FaceHandle> = faces.into_iter().collect();
    faces.sort();
    faces
}

/// Find outer face handles adjacent to an inner vertex.
fn find_adjacent_outer_faces(
    inner: &BRep,
    vertex: VertexHandle,
    reverse_face: &HashMap<FaceHandle, FaceHandle>,
) -> Vec<FaceHandle> {
    let mut faces = HashSet::new();
    let start_he = inner.vertices[vertex].outgoing;

    if !inner.half_edges.contains_key(start_he) {
        return Vec::new();
    }

    // Walk the half-edge fan around the vertex.
    let mut current = start_he;
    loop {
        let he = &inner.half_edges[current];
        let loop_h = he.loop_;
        if inner.loops.contains_key(loop_h) {
            let face_h = inner.loops[loop_h].face;
            if let Some(&outer_fh) = reverse_face.get(&face_h) {
                faces.insert(outer_fh);
            }
        }

        // Move to the next half-edge around the vertex:
        // twin's next gives the next outgoing half-edge from this vertex.
        let twin = he.twin;
        if !inner.half_edges.contains_key(twin) {
            break;
        }
        let next_outgoing = inner.half_edges[twin].next;
        if next_outgoing == start_he || !inner.half_edges.contains_key(next_outgoing) {
            break;
        }
        current = next_outgoing;
    }

    // Sorted handles for deterministic constraint order (PK-277 class).
    let mut faces: Vec<FaceHandle> = faces.into_iter().collect();
    faces.sort();
    faces
}

/// Solve for the intersection of 3 planes using Cramer's rule.
///
/// Given planes `n_i . x = d_i`, builds the 3x3 system and solves.
/// Returns `None` if the system is singular.
fn solve_3_planes(
    plane0: &(Vec3, f64),
    plane1: &(Vec3, f64),
    plane2: &(Vec3, f64),
) -> Option<Point3> {
    let (n0, d0) = plane0;
    let (n1, d1) = plane1;
    let (n2, d2) = plane2;

    // Compute determinant using triple product: n0 . (n1 x n2).
    let cross12 = n1.cross(n2);
    let det = n0.dot(&cross12);

    // Use a tight numerical threshold for the 3x3 determinant singularity
    // check.  This is a pure linear-algebra guard — not a geometric tolerance —
    // so 1e-12 is appropriate regardless of the tolerance object.
    if det.abs() < 1e-12 {
        return None;
    }

    let cross20 = n2.cross(n0);
    let cross01 = n0.cross(n1);

    let inv_det = 1.0 / det;
    let x = (*d0 * cross12.x + *d1 * cross20.x + *d2 * cross01.x) * inv_det;
    let y = (*d0 * cross12.y + *d1 * cross20.y + *d2 * cross01.y) * inv_det;
    let z = (*d0 * cross12.z + *d1 * cross20.z + *d2 * cross01.z) * inv_det;

    Some(Point3::new(x, y, z))
}

/// Least-squares point for N ≥ 3 plane constraints `n_i · x = d_i`:
/// solves the 3×3 normal equations `(Σ n_i n_iᵀ) x = Σ d_i n_i` by
/// Cramer's rule. Equals the exact intersection when one exists.
fn solve_planes_least_squares(planes: &[(Vec3, f64)]) -> Option<Point3> {
    let mut m = [[0.0_f64; 3]; 3];
    let mut b = [0.0_f64; 3];
    for (n, d) in planes {
        let nv = [n.x, n.y, n.z];
        for r in 0..3 {
            for c in 0..3 {
                m[r][c] += nv[r] * nv[c];
            }
            b[r] += d * nv[r];
        }
    }
    let det3 = |a: &[[f64; 3]; 3]| -> f64 {
        a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
            - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
            + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0])
    };
    let det = det3(&m);
    if det.abs() < 1e-12 {
        return None;
    }
    let mut out = [0.0_f64; 3];
    for (axis, slot) in out.iter_mut().enumerate() {
        let mut mm = m;
        for r in 0..3 {
            mm[r][axis] = b[r];
        }
        *slot = det3(&mm) / det;
    }
    Some(Point3::new(out[0], out[1], out[2]))
}

/// Exact minimal-displacement solve for 1 or 2 independent plane
/// constraints: project the seed onto the plane, or onto the two
/// planes' intersection line. Near-parallel pairs collapse to the
/// single-plane case (their constraint directions coincide).
fn solve_underdetermined_planes(planes: &[(Vec3, f64)], seed: Point3) -> Option<Point3> {
    // Distinct constraint directions only.
    let mut distinct: Vec<(Vec3, f64)> = Vec::new();
    for (n, d) in planes {
        let n_unit = n.normalized()?;
        let dup = distinct
            .iter()
            .any(|(m, _)| m.cross(&n_unit).length() < 1e-9);
        if !dup {
            let scale = n.length();
            distinct.push((n_unit, d / scale));
        }
    }
    match distinct.len() {
        1 => {
            let (n, d) = distinct[0];
            let dist = seed.x * n.x + seed.y * n.y + seed.z * n.z - d;
            Some(seed - n * dist)
        }
        2 => {
            // Line of intersection, then closest point on it to seed.
            let (n1, d1) = distinct[0];
            let (n2, d2) = distinct[1];
            let dir = n1.cross(&n2).normalized()?;
            let p = solve_point_on_two_planes(n1, d1, n2, d2, dir)?;
            let t = (seed - p).dot(&dir);
            Some(p + dir * t)
        }
        _ => solve_planes_least_squares(&distinct),
    }
}

/// A point on the intersection of two planes `n_i · x = d_i`, fixing
/// the coordinate with the largest |dir| component to zero.
fn solve_point_on_two_planes(n1: Vec3, d1: f64, n2: Vec3, d2: f64, dir: Vec3) -> Option<Point3> {
    let solve_2x2 = |a: f64, b: f64, c1: f64, d: f64, e: f64, c2: f64| -> Option<(f64, f64)> {
        let det = a * e - b * d;
        if det.abs() < 1e-14 {
            return None;
        }
        Some(((c1 * e - b * c2) / det, (a * c2 - c1 * d) / det))
    };
    let (ax, ay, az) = (dir.x.abs(), dir.y.abs(), dir.z.abs());
    if ax >= ay && ax >= az {
        let (y, z) = solve_2x2(n1.y, n1.z, d1, n2.y, n2.z, d2)?;
        Some(Point3::new(0.0, y, z))
    } else if ay >= ax && ay >= az {
        let (x, z) = solve_2x2(n1.x, n1.z, d1, n2.x, n2.z, d2)?;
        Some(Point3::new(x, 0.0, z))
    } else {
        let (x, y) = solve_2x2(n1.x, n1.y, d1, n2.x, n2.y, d2)?;
        Some(Point3::new(x, y, 0.0))
    }
}

#[cfg(test)]
#[path = "vertex_recompute_tests.rs"]
mod tests;
