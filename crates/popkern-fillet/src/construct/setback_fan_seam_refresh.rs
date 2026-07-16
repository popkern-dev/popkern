//! Post-reconcile refresh of the setback fan's internal seam curves.
//!
//! The fan wiring writes cap-cap edge curves from the ASSEMBLY's
//! boundary rows; the rim reconciliation then moves the committed
//! patches' depth-0..2 rows onto the live sweeps, leaving those
//! curves stale by the move size (measured 3.5e-5 curve-vs-surface
//! on the pyramid5_mixed sweep fan while the two patch surfaces
//! agree to 7.6e-7 across the same seams). Every internal seam of
//! the fan is a full bi-9 boundary row shared verbatim by its two
//! patches after CP re-unification, so the honest curve is that row
//! itself: this pass rewrites each cap-cap edge curve as the exact
//! degree-9 Bezier of the matched boundary row.
//!
//! A row only replaces a curve when both endpoints AND the mid-span
//! agree with the existing curve within a staleness band — an edge
//! that is not a full boundary row (nothing in the current wiring
//! produces one) is left untouched rather than guessed at.
//!
//! Nested inside `setback_vertex_topology` (shares its private
//! bezier9 helpers).

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::gregory::GregoryPatch;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{NurbsCurve, Point3};
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, VertexHandle};

use super::bezier9_eval;

/// Worst cross-seam quality of the committed fan, measured the same
/// way the acceptance tests do: per cap-cap edge, each side's pcurve
/// image against the shared curve at matched fractions (`g0`) and the
/// two sides' finite-difference normal misalignment (`g1`, cross
/// product norm). Samples whose chart is length-degenerate (collapsed
/// patch corners: one partial under 1e-3 of the other) are skipped —
/// their FD normals are direction noise, not geometry. This is the
/// surgery's COMMIT POSTCONDITION: a fan whose seams are not
/// commercial must decline (rollback) rather than ship — measured on
/// pentagon5_blunt, whose wide flat fan creases 2.9e-1 at one seam
/// corner while every sharp-apex fixture sits at g1 ≤ 2e-3.
pub(super) fn worst_fan_seam_quality(brep: &BRep, cap_faces: &[FaceHandle]) -> (f64, f64, f64) {
    use popkern_core::geometry::surface::SurfaceEval;
    let cap_set: std::collections::BTreeSet<FaceHandle> = cap_faces.iter().copied().collect();
    let (mut worst_g0, mut worst_g1, mut worst_g2) = (0.0_f64, 0.0_f64, 0.0_f64);
    for (_, edge) in brep.edges.iter() {
        let Some(he_data) = brep.half_edges.get(edge.half_edge) else {
            continue;
        };
        let Some(twin_data) = brep.half_edges.get(he_data.twin) else {
            continue;
        };
        let (Some(loop_a), Some(loop_b)) = (
            brep.loops.get(he_data.loop_),
            brep.loops.get(twin_data.loop_),
        ) else {
            continue;
        };
        if loop_a.face == loop_b.face
            || !cap_set.contains(&loop_a.face)
            || !cap_set.contains(&loop_b.face)
        {
            continue;
        }
        let Some(curve) = edge.curve.as_ref() else {
            continue;
        };
        let ce: &dyn CurveEval = &curve.geometry;
        let p_start = ce.point_at(curve.t_start);
        let side = |face: FaceHandle,
                    he: popkern_core::topology::HalfEdgeHandle|
         -> Option<(
            &popkern_core::geometry::SurfaceGeometry,
            &popkern_core::geometry::pcurve::Pcurve,
            bool,
        )> {
            let pc = brep.half_edges.get(he)?.pcurve.as_ref()?;
            let surf = &brep.faces.get(face)?.surface.as_ref()?.geometry;
            let uv0 = pc.point_at(pc.t_start);
            let uv1 = pc.point_at(pc.t_end);
            let rev = surf.point_at(uv1.x, uv1.y).distance_to(&p_start)
                < surf.point_at(uv0.x, uv0.y).distance_to(&p_start);
            Some((surf, pc, rev))
        };
        let (Some(sa), Some(sb)) = (
            side(loop_a.face, edge.half_edge),
            side(loop_b.face, he_data.twin),
        ) else {
            continue;
        };
        let normal_at = |surf: &popkern_core::geometry::SurfaceGeometry, u: f64, v: f64| {
            let h = 1e-6;
            let du = surf.point_at(u + h, v) - surf.point_at(u - h, v);
            let dv = surf.point_at(u, v + h) - surf.point_at(u, v - h);
            let (lu, lv) = (du.length(), dv.length());
            if lu.min(lv) < 1e-3 * lu.max(lv) {
                return None;
            }
            du.cross(&dv).normalized()
        };
        // Sampling density matches the acceptance metric (100 per
        // seam): the pentagon5_blunt corner crease is a narrow spike
        // (2.9e-1 at f≈0.95) that a 9-sample grid reads as 2.6e-2 —
        // it must not slip through the gate between samples.
        for k in 0..100 {
            let f = 0.05 + 0.90 * k as f64 / 99.0;
            let t = curve.t_start + f * (curve.t_end - curve.t_start);
            let c = ce.point_at(t);
            // Seam direction for the normal-curvature comparison
            // (matches the live-seam acceptance metric).
            let dt = 1e-4 * (curve.t_end - curve.t_start).abs().max(1e-12);
            let seam_dir = (ce.point_at(t + dt) - c).normalized();
            let mut normals = Vec::with_capacity(2);
            let mut curvatures = Vec::with_capacity(2);
            for (surf, pc, rev) in [sa, sb] {
                let g = if rev { 1.0 - f } else { f };
                let uv = pc.point_at(pc.t_start + g * (pc.t_end - pc.t_start));
                let s = surf.point_at(uv.x, uv.y);
                worst_g0 = worst_g0.max(s.distance_to(&c));
                let normal = normal_at(surf, uv.x, uv.y);
                if let (Some(normal), Some(dir)) = (normal, seam_dir) {
                    curvatures.push(normal_curvature_at(surf, uv.x, uv.y, normal, dir));
                } else {
                    curvatures.push(None);
                }
                normals.push(normal);
            }
            if let (Some(Some(na)), Some(Some(nb))) = (normals.first(), normals.get(1)) {
                worst_g1 = worst_g1.max(na.cross(nb).length());
            }
            if let (Some(Some(ka)), Some(Some(kb))) = (curvatures.first(), curvatures.get(1)) {
                let kmax = ka.abs().max(kb.abs());
                worst_g2 = worst_g2.max((ka - kb).abs() / (1.0 + kmax));
            }
        }
    }
    (worst_g0, worst_g1, worst_g2)
}

/// Normal curvature of `surf` at `(u, v)` along the 3-D direction
/// `dir` (κ_n = S''·n / |S'|², matching the live-seam acceptance
/// metric): solve `dir ≈ a·Su + b·Sv` in the tangent plane, then take
/// a second difference along `(a, b)` in UV.
fn normal_curvature_at(
    surf: &popkern_core::geometry::SurfaceGeometry,
    u: f64,
    v: f64,
    normal: popkern_core::geometry::Vec3,
    dir: popkern_core::geometry::Vec3,
) -> Option<f64> {
    use popkern_core::geometry::surface::SurfaceEval;
    let h = 1e-6;
    let su = (surf.point_at(u + h, v) - surf.point_at(u - h, v)) * (0.5 / h);
    let sv = (surf.point_at(u, v + h) - surf.point_at(u, v - h)) * (0.5 / h);
    let (e, f_, g) = (su.dot(&su), su.dot(&sv), sv.dot(&sv));
    let (r1, r2) = (dir.dot(&su), dir.dot(&sv));
    let det = e * g - f_ * f_;
    if det.abs() < 1e-30 {
        return None;
    }
    let a = (r1 * g - r2 * f_) / det;
    let b = (r2 * e - r1 * f_) / det;
    let step = 1e-4;
    let p0 = surf.point_at(u, v);
    let p1 = surf.point_at(u + a * step, v + b * step);
    let p2 = surf.point_at(u - a * step, v - b * step);
    let second = (p1 - p0) + (p2 - p0);
    let first_len = (p1 - p2).length() * 0.5;
    if first_len < 1e-30 {
        return None;
    }
    Some(second.dot(&normal) / (first_len * first_len))
}

/// Max distance between the stale curve and the candidate row,
/// relative to the edge chord. Reconcile moves are a few 1e-4 on
/// ~1e-2 chords; a WRONG row sits at least a patch width (≈ one
/// chord) away. An absolute band was measured unsafe: 5e-2 exceeded
/// octagon8_mixed's whole patch size, and wrong-row rewrites tore
/// the fan to g0 = 2.8.
fn staleness_band(chord: f64) -> f64 {
    (0.2 * chord).max(1.0e-4)
}

/// Rewrite every cap-cap edge curve from the shared boundary row of
/// its two (already re-unified) patches. Returns the number of edges
/// rewritten and the worst curve movement applied.
pub(super) fn refresh_fan_seam_curves(brep: &mut BRep, cap_faces: &[FaceHandle]) -> (usize, f64) {
    let cap_set: std::collections::BTreeSet<FaceHandle> = cap_faces.iter().copied().collect();
    let mut rewrites: Vec<(EdgeHandle, [Point3; 10], f64)> = Vec::new();
    let mut vertex_moves: std::collections::BTreeMap<VertexHandle, Point3> =
        std::collections::BTreeMap::new();
    for (eh, edge) in brep.edges.iter() {
        let he = edge.half_edge;
        let Some(he_data) = brep.half_edges.get(he) else {
            continue;
        };
        let Some(twin_data) = brep.half_edges.get(he_data.twin) else {
            continue;
        };
        let (Some(loop_a), Some(loop_b)) = (
            brep.loops.get(he_data.loop_),
            brep.loops.get(twin_data.loop_),
        ) else {
            continue;
        };
        if loop_a.face == loop_b.face
            || !cap_set.contains(&loop_a.face)
            || !cap_set.contains(&loop_b.face)
        {
            continue;
        }
        let Some(curve) = edge.curve.as_ref() else {
            continue;
        };
        let ce: &dyn CurveEval = &curve.geometry;
        let p0 = ce.point_at(curve.t_start);
        let p1 = ce.point_at(curve.t_end);
        let pm = ce.point_at(0.5 * (curve.t_start + curve.t_end));
        let Some(face) = brep.faces.get(loop_a.face) else {
            continue;
        };
        let Some(SurfaceGeometry::Gregory(GregoryPatch::RectangularBiNonic(bp))) =
            face.surface.as_ref().map(|s| &s.geometry)
        else {
            continue;
        };
        // Candidate rows: the four boundary rows of the bi-9 grid, in
        // both orientations.
        let mut best: Option<([Point3; 10], f64)> = None;
        for side in 0..4usize {
            let row: [Point3; 10] = core::array::from_fn(|k| match side {
                0 => bp.grid[0][k],
                1 => bp.grid[9][k],
                2 => bp.grid[k][0],
                _ => bp.grid[k][9],
            });
            for reversed in [false, true] {
                let oriented: [Point3; 10] = if reversed {
                    core::array::from_fn(|k| row[9 - k])
                } else {
                    row
                };
                let score = oriented[0]
                    .distance_to(&p0)
                    .max(oriented[9].distance_to(&p1))
                    .max(bezier9_eval(&oriented, 0.5).distance_to(&pm));
                if best.as_ref().map(|(_, s)| score < *s).unwrap_or(true) {
                    best = Some((oriented, score));
                }
            }
        }
        let Some((row, score)) = best else {
            continue;
        };
        if score > staleness_band(p0.distance_to(&p1)) {
            continue;
        }
        // Endpoint policy: EVERY end vertex follows its row end. The
        // rim vertices were split on the extracted arc cubics, but
        // the rim reconciliation then landed the cap rows on the
        // sweep SURFACE — measured 1e-5..3.2e-4 vertex-to-row-end
        // gaps on the pyramid5_mixed sweep fan, and snapping the row
        // to the stale vertex bends the seam curve off both patch
        // surfaces by exactly the residual this pass exists to close.
        // Incident non-seam curves are re-snapped to the moved
        // vertex below.
        let (Some(pa), Some(pb)) = (
            brep.vertices.get(he_data.origin).map(|v| v.position),
            brep.vertices.get(twin_data.origin).map(|v| v.position),
        ) else {
            continue;
        };
        let v_start = if pa.distance_to(&p0) + pb.distance_to(&p1)
            <= pa.distance_to(&p1) + pb.distance_to(&p0)
        {
            (he_data.origin, twin_data.origin)
        } else {
            (twin_data.origin, he_data.origin)
        };
        for (slot, vertex) in [(0usize, v_start.0), (9usize, v_start.1)] {
            vertex_moves.entry(vertex).or_insert(row[slot]);
        }
        // Movement actually applied (how stale the old curve was).
        let moved = (0..=8)
            .map(|k| {
                let f = 0.02 + 0.96 * k as f64 / 8.0;
                let t = curve.t_start + f * (curve.t_end - curve.t_start);
                ce.point_at(t).distance_to(&bezier9_eval(&row, f))
            })
            .fold(0.0_f64, f64::max);
        rewrites.push((eh, row, moved));
    }
    // Apply vertex moves, then re-snap the endpoint CPs of every
    // incident curve that is NOT being rewritten from a row (the arc
    // sub-cubics on the sweeps and the corner segments on the lateral
    // faces). Both the old and new endpoint lie on the neighbouring
    // surface, so the snap keeps those curves on-surface to second
    // order in the (≤ 3.2e-4) move.
    let rewritten: std::collections::BTreeSet<EdgeHandle> =
        rewrites.iter().map(|(eh, _, _)| *eh).collect();
    let mut incident_snaps: Vec<(EdgeHandle, bool, Point3)> = Vec::new();
    for (eh, edge) in brep.edges.iter() {
        if rewritten.contains(&eh) {
            continue;
        }
        let Some(he_data) = brep.half_edges.get(edge.half_edge) else {
            continue;
        };
        let Some(twin_data) = brep.half_edges.get(he_data.twin) else {
            continue;
        };
        let Some(curve) = edge.curve.as_ref() else {
            continue;
        };
        let ce: &dyn CurveEval = &curve.geometry;
        for vertex in [he_data.origin, twin_data.origin] {
            let (Some(&target), Some(old)) = (
                vertex_moves.get(&vertex),
                brep.vertices.get(vertex).map(|v| v.position),
            ) else {
                continue;
            };
            let d_start = ce.point_at(curve.t_start).distance_to(&old);
            let d_end = ce.point_at(curve.t_end).distance_to(&old);
            let at_start = d_start <= d_end;
            if d_start.min(d_end) < 1.0e-3 {
                incident_snaps.push((eh, at_start, target));
            }
        }
    }
    for (vertex, position) in &vertex_moves {
        if let Some(v) = brep.vertices.get_mut(*vertex) {
            v.position = *position;
        }
    }
    for (eh, at_start, target) in incident_snaps {
        let Some(edge) = brep.edges.get_mut(eh) else {
            continue;
        };
        let Some(curve) = edge.curve.as_mut() else {
            continue;
        };
        if let popkern_core::geometry::CurveGeometry::Nurbs(nc) = &mut curve.geometry
            && !nc.control_points.is_empty()
        {
            // Clamped curves: first/last CP are the t_start/t_end
            // images.
            let idx = if at_start {
                0
            } else {
                nc.control_points.len() - 1
            };
            nc.control_points[idx] = target;
        }
    }
    let mut worst_moved = 0.0_f64;
    let count = rewrites.len();
    for (eh, row, moved) in rewrites {
        worst_moved = worst_moved.max(moved);
        if let Some(edge) = brep.edges.get_mut(eh) {
            edge.curve = Some(popkern_core::geometry::CurveOnEdge {
                geometry: popkern_core::geometry::CurveGeometry::Nurbs(NurbsCurve {
                    degree: 9,
                    knots: [vec![0.0; 10], vec![1.0; 10]].concat(),
                    control_points: row.to_vec(),
                    weights: vec![1.0; 10],
                }),
                t_start: 0.0,
                t_end: 1.0,
            });
        }
    }
    (count, worst_moved)
}
