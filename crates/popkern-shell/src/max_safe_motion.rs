//! Motion-model arms of [`crate::max_safe`] (PK-328): first-order
//! vertex velocities for the pre-flight edge-collapse bound, and the
//! thin-wall clearance sweep over trimmed patches. Split from
//! `max_safe.rs` (500-line ratchet).

use std::collections::{HashMap, HashSet};

use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::geometry::{Point3, Vec3, uv_projection};
use popkern_core::topology::{BRep, FaceHandle, VertexHandle};

use crate::body_distance::FaceTrimEntry;
use crate::max_safe::OffsetConstraint;

/// Vertex → adjacent classified faces, in first-encounter order.
pub(crate) fn vertex_face_adjacency(brep: &BRep) -> HashMap<VertexHandle, Vec<FaceHandle>> {
    let mut adjacency: HashMap<VertexHandle, Vec<FaceHandle>> = HashMap::new();
    for (_, he) in &brep.half_edges {
        let Some(loop_data) = brep.loops.get(he.loop_) else {
            continue;
        };
        let faces = adjacency.entry(he.origin).or_default();
        if !faces.contains(&loop_data.face) {
            faces.push(loop_data.face);
        }
    }
    adjacency
}

/// First-order vertex velocities per unit offset distance: the
/// least-squares solution of `nᵢ·v = mᵢ` over the vertex's adjacent
/// faces (exact for plane triple points). A small ridge keeps
/// under-determined seam vertices at the minimal-norm solution.
pub(crate) fn vertex_velocities(
    brep: &BRep,
    motions: &HashMap<FaceHandle, f64>,
) -> HashMap<VertexHandle, Vec3> {
    let adjacency = vertex_face_adjacency(brep);

    let mut velocities = HashMap::new();
    for (vh, faces) in adjacency {
        let Some(vertex) = brep.vertices.get(vh) else {
            continue;
        };
        let p = vertex.position;
        // Normal equations AᵀA v = Aᵀ b over (n_f, m_f) rows.
        let mut ata = [[0.0f64; 3]; 3];
        let mut atb = [0.0f64; 3];
        let mut rows = 0usize;
        for fh in faces {
            let Some(&m) = motions.get(&fh) else { continue };
            let Some(face) = brep.faces.get(fh) else {
                continue;
            };
            let Some(sof) = face.surface.as_ref() else {
                continue;
            };
            let Some((u, v)) = uv_projection::project_point_to_uv(&sof.geometry, &p) else {
                continue;
            };
            let n_raw = sof.geometry.normal_at(u, v);
            let len = n_raw.length();
            if len < 1e-12 {
                continue;
            }
            let n = n_raw * (1.0 / len);
            let row = [n.x, n.y, n.z];
            for (r, &row_r) in row.iter().enumerate() {
                for (c, &row_c) in row.iter().enumerate() {
                    ata[r][c] += row_r * row_c;
                }
                atb[r] += row_r * m;
            }
            rows += 1;
        }
        if rows == 0 {
            continue;
        }
        for (r, row) in ata.iter_mut().enumerate() {
            row[r] += 1.0e-9; // ridge: minimal-norm on under-determined vertices
        }
        if let Some(v) = solve3(&ata, &atb) {
            velocities.insert(vh, v);
        }
    }
    velocities
}

/// Curved-aware refinement of a linear edge-collapse bound (PK-328):
/// the linear model reads normals at the ORIGINAL corner, but on
/// curved faces the recomputed vertex lands where the surface has
/// risen or fallen (measured: the wavy-plate bump's tails slow its
/// corner-edge shrink from rate 2.0 to 1.89 — a 5.5% bound error,
/// outside the slider-consistency pin). Re-solve both endpoint
/// vertices against the exact offset tangent planes at trial
/// distances and secant the signed chord projection to zero.
pub(crate) fn refine_edge_collapse_bound(
    brep: &BRep,
    adjacency: &HashMap<VertexHandle, Vec<FaceHandle>>,
    motions: &HashMap<FaceHandle, f64>,
    v0: VertexHandle,
    v1: VertexHandle,
    d_linear: f64,
) -> Option<f64> {
    let p0 = brep.vertices.get(v0)?.position;
    let p1 = brep.vertices.get(v1)?.position;
    let chord = p1 - p0;
    let len0 = chord.length();
    if len0 <= 0.0 {
        return None;
    }
    let dir = chord * (1.0 / len0);
    let shrink = |d: f64| -> Option<f64> {
        let q0 = vertex_at_distance(brep, adjacency, motions, v0, d)?;
        let q1 = vertex_at_distance(brep, adjacency, motions, v1, d)?;
        Some((q1 - q0).dot(&dir))
    };
    // Secant from the linear bound (where the linear model predicts
    // projection zero) and a nearby second point.
    let mut a = d_linear;
    let mut fa = shrink(a)?;
    let mut b = d_linear * 0.9;
    let mut fb = shrink(b)?;
    for _ in 0..6 {
        if (fa - fb).abs() < 1e-14 {
            break;
        }
        let c = a - fa * (a - b) / (fa - fb);
        // Keep iterates in a sane positive window around the estimate.
        let c = c.clamp(d_linear * 0.25, d_linear * 4.0);
        if (c - a).abs() < d_linear * 1e-6 {
            return Some(c);
        }
        b = a;
        fb = fa;
        a = c;
        fa = shrink(a)?;
    }
    (a > 0.0).then_some(a)
}

/// Solve one vertex position at offset distance `d`: iterate the
/// tangent-plane intersection of the adjacent faces' EXACT offsets
/// (`nᵢ·p = nᵢ·qᵢ + d·mᵢ`, with `qᵢ` the foot of the previous iterate
/// on the source surface). Exact for plane corners in one step;
/// converges on curved corners where the linear model is first-order.
fn vertex_at_distance(
    brep: &BRep,
    adjacency: &HashMap<VertexHandle, Vec<FaceHandle>>,
    motions: &HashMap<FaceHandle, f64>,
    vh: VertexHandle,
    d: f64,
) -> Option<Point3> {
    let faces = adjacency.get(&vh)?;
    let start = brep.vertices.get(vh)?.position;
    let mut p = start;
    let mut seeds: HashMap<FaceHandle, (f64, f64)> = HashMap::new();
    for _ in 0..4 {
        let mut ata = [[0.0f64; 3]; 3];
        let mut atb = [0.0f64; 3];
        let mut rows = 0usize;
        for &fh in faces {
            let Some(&m) = motions.get(&fh) else { continue };
            let Some(sof) = brep.faces.get(fh).and_then(|f| f.surface.as_ref()) else {
                continue;
            };
            let uv = match seeds.get(&fh) {
                Some(&(su, sv)) => {
                    uv_projection::project_point_to_uv_near(&sof.geometry, &p, su, sv)
                }
                None => uv_projection::project_point_to_uv(&sof.geometry, &p),
            };
            let Some((u, v)) = uv else { continue };
            seeds.insert(fh, (u, v));
            let q = sof.geometry.point_at(u, v);
            let n_raw = sof.geometry.normal_at(u, v);
            let n_len = n_raw.length();
            if n_len < 1e-12 {
                continue;
            }
            let n = n_raw * (1.0 / n_len);
            let rhs = n.dot(&(q - Point3::new(0.0, 0.0, 0.0))) + d * m;
            let row = [n.x, n.y, n.z];
            for (r, &row_r) in row.iter().enumerate() {
                for (c, &row_c) in row.iter().enumerate() {
                    ata[r][c] += row_r * row_c;
                }
                atb[r] += row_r * rhs;
            }
            rows += 1;
        }
        if rows == 0 {
            return None;
        }
        // Ridge toward the CURRENT iterate: under-determined vertices
        // (two-face seams) keep their tangential position instead of
        // drifting to the origin.
        let lambda = 1.0e-9;
        for (r, row) in ata.iter_mut().enumerate() {
            row[r] += lambda;
        }
        atb[0] += lambda * p.x;
        atb[1] += lambda * p.y;
        atb[2] += lambda * p.z;
        let v = solve3(&ata, &atb)?;
        p = Point3::new(v.x, v.y, v.z);
    }
    Some(p)
}

/// Solve a 3×3 linear system by Cramer's rule (`None` on a singular
/// matrix — the ridge above makes this unreachable in practice).
fn solve3(a: &[[f64; 3]; 3], b: &[f64; 3]) -> Option<Vec3> {
    let det = |m: &[[f64; 3]; 3]| -> f64 {
        m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
            - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
    };
    let d = det(a);
    if d.abs() < 1e-30 {
        return None;
    }
    let mut out = [0.0f64; 3];
    for (k, o) in out.iter_mut().enumerate() {
        let mut m = *a;
        for r in 0..3 {
            m[r][k] = b[r];
        }
        *o = det(&m) / d;
    }
    Some(Vec3::new(out[0], out[1], out[2]))
}

/// Thin-wall clearance candidates: non-adjacent (no shared vertex)
/// face pairs whose trimmed patches approach under the per-face
/// motions. Boundary samples of each face are measured against the
/// other's trimmed patch; the closing speed per unit distance is
/// `s = (v_a − v_b)·û`, so the wall closes at `gap/s`.
pub(crate) fn thin_wall_candidates(
    brep: &BRep,
    motions: &HashMap<FaceHandle, f64>,
    lin_tol: f64,
    candidates: &mut Vec<(f64, OffsetConstraint)>,
) {
    let mut face_handles: Vec<FaceHandle> = motions.keys().copied().collect();
    face_handles.sort();

    struct WallFace {
        fh: FaceHandle,
        entry: FaceTrimEntry,
        verts: HashSet<VertexHandle>,
        aabb: (Point3, Point3),
        motion: f64,
    }
    // Face → boundary vertices (adjacency exemption).
    let mut face_verts: HashMap<FaceHandle, HashSet<VertexHandle>> = HashMap::new();
    for (_, he) in &brep.half_edges {
        if let Some(loop_data) = brep.loops.get(he.loop_) {
            face_verts
                .entry(loop_data.face)
                .or_default()
                .insert(he.origin);
        }
    }

    let mut walls: Vec<WallFace> = Vec::new();
    for &fh in &face_handles {
        let Some(entry) = FaceTrimEntry::build(brep, fh) else {
            continue;
        };
        let samples = entry.boundary_samples();
        if samples.is_empty() {
            continue;
        }
        let mut lo = samples[0];
        let mut hi = samples[0];
        for p in samples {
            lo = Point3::new(lo.x.min(p.x), lo.y.min(p.y), lo.z.min(p.z));
            hi = Point3::new(hi.x.max(p.x), hi.y.max(p.y), hi.z.max(p.z));
        }
        walls.push(WallFace {
            fh,
            entry,
            verts: face_verts.get(&fh).cloned().unwrap_or_default(),
            aabb: (lo, hi),
            motion: motions[&fh],
        });
    }

    let aabb_gap = |a: &(Point3, Point3), b: &(Point3, Point3)| -> f64 {
        let dx = (b.0.x - a.1.x).max(a.0.x - b.1.x).max(0.0);
        let dy = (b.0.y - a.1.y).max(a.0.y - b.1.y).max(0.0);
        let dz = (b.0.z - a.1.z).max(a.0.z - b.1.z).max(0.0);
        (dx * dx + dy * dy + dz * dz).sqrt()
    };

    let mut best_wall = f64::INFINITY;
    for i in 0..walls.len() {
        for j in (i + 1)..walls.len() {
            let (a, b) = (&walls[i], &walls[j]);
            if !a.verts.is_disjoint(&b.verts) {
                continue; // adjacent: edge/join machinery owns these
            }
            // Closing speed is at most 2, so a pair whose AABB gap
            // already exceeds twice the best bound cannot improve it.
            if aabb_gap(&a.aabb, &b.aabb) / 2.0 >= best_wall {
                continue;
            }
            for (from, to) in [(a, b), (b, a)] {
                for (idx, p) in from.entry.boundary_samples().iter().enumerate() {
                    let Some(near) = to.entry.nearest(p) else {
                        continue;
                    };
                    if near.distance <= lin_tol * 10.0 {
                        continue; // touching band: not a wall, a contact
                    }
                    let u_hat = (near.point - *p) * (1.0 / near.distance);
                    let n_from = from.entry.boundary_normal(idx);
                    let n_to = to.entry.surface_normal(near.uv);
                    let (Some(n_from), Some(n_to)) = (n_from, n_to) else {
                        continue;
                    };
                    // A real wall collision needs the near face to
                    // ADVANCE along the gap: a raw closing-speed test
                    // is material-blind and reported the two corner
                    // spheres of a convex filleted cube "closing"
                    // THROUGH the solid under outward growth (their
                    // facing flanks drift together at s≈0.74 while
                    // both retreat from the gap segment). An oblique
                    // far face may retreat slower than the near face
                    // advances — that closure is real and stays.
                    let advance = (n_from * from.motion).dot(&u_hat);
                    if advance <= 1.0e-3 {
                        continue;
                    }
                    let s = advance - (n_to * to.motion).dot(&u_hat);
                    if s > 1.0e-9 {
                        let d = near.distance / s;
                        if d >= best_wall {
                            continue; // only the tightest bound is reported
                        }
                        // SHIELD exclusion: when a third face passes
                        // near the gap segment, that face owns the
                        // local geometry and the skins chain through
                        // it instead of colliding (measured: the
                        // bracket's 45° countersink flank CONTAINS the
                        // whole top-rim→bore-junction segment, and the
                        // hollow skins stay consistent to 1.5× the
                        // naive gap/s estimate).
                        let mid = *p + (near.point - *p) * 0.5;
                        let shielded = walls.iter().enumerate().any(|(k, w)| {
                            k != i
                                && k != j
                                && w.entry
                                    .distance(&mid)
                                    .is_some_and(|sd| sd < near.distance * 0.25)
                        });
                        if shielded {
                            continue;
                        }
                        if std::env::var_os("PK_MAXSAFE_TRACE").is_some() {
                            eprintln!(
                                "[maxsafe-wall] {:?}x{:?} gap={:.4} s={:.3} d={:.4} \
                                 p=({:.3},{:.3},{:.3}) q=({:.3},{:.3},{:.3})",
                                from.fh,
                                to.fh,
                                near.distance,
                                s,
                                d,
                                p.x,
                                p.y,
                                p.z,
                                near.point.x,
                                near.point.y,
                                near.point.z
                            );
                        }
                        if d < best_wall {
                            best_wall = d;
                        }
                        candidates.push((
                            d,
                            OffsetConstraint::ThinWall {
                                face_a: a.fh,
                                face_b: b.fh,
                            },
                        ));
                    }
                }
            }
        }
    }
}
