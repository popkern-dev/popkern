//! Tolerance-scoped T-junction zippering (PK-208).
//!
//! When a face falls back to the legacy tessellation path, its
//! boundary is sampled at its own density; the neighbouring face —
//! boundary-first or legacy — samples the same edge at a different
//! density. The sewing pass lands both chains ON the shared curve,
//! but at STAGGERED positions, so position-welding can never pair
//! them: the mesh reads open along a seam whose geometry is closed
//! (measured: 236 welded open edges on the blunt-pentagon fan whose
//! seams agree to 2e-13).
//!
//! This pass closes exactly that class: an open mesh edge is split at
//! opposing open-boundary vertices that lie on it, so both chains end
//! up sharing every subdivision point and the weld closes.
//!
//! # Why the scoping is the whole design
//!
//! Zippering unconditionally would blind the watertightness gates: a
//! REAL geometric hole (a flush-boolean gap, a dropped face) would be
//! silently stitched shut, and "watertight" would stop being
//! evidence. The zipper therefore only acts where the BRep TOPOLOGY
//! says the mesh must be closed: on open mesh edges lying along a
//! TWO-SIDED BRep edge (distinct faces on both half-edges). There,
//! openness can only be a sampling artifact. The matching band is
//! tight (2e-4) unless the edge carries a RECORDED entity tolerance —
//! the check-9 contract — in which case the construction explicitly
//! accepted that inexactness and the band widens to 1.5x it. A hole
//! whose opposing chain is absent or beyond the band stays open —
//! honest.
//!
//! The pass reports how many junctions it closed via
//! [`crate::mesh::TriMesh::zipped_t_junctions`]; a consumer that
//! wants to gate on "watertight without repair" can assert that
//! count is zero.

use std::collections::HashMap;

use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::CurveEval;
use popkern_core::topology::BRep;

use crate::mesh::TriMesh;

/// Weld-grid resolution — must match the welded-open-edge census the
/// example gates use.
const WELD: f64 = 1.0e7;
/// Maximum zip rounds; each round can expose junctions the previous
/// round's splits created.
const MAX_ROUNDS: usize = 8;
/// Fractional interior window for a split along an open edge: a
/// candidate at (or numerically at) an endpoint is already welded.
const END_FRACTION: f64 = 1.0e-6;

/// Default matching band along a topologically two-sided edge with
/// no recorded tolerance: covers sampling/sewing offsets while
/// staying far under any real geometric gap worth catching.
const DEFAULT_BAND: f64 = 2.0e-4;

/// Close T-junction cracks along two-sided BRep edges. Returns the
/// number of junction splits performed.
pub(crate) fn zip_recorded_tolerance_cracks(brep: &BRep, mesh: &mut TriMesh) -> usize {
    // Sampled polylines of every TWO-SIDED edge (distinct faces on
    // the two half-edges — the topology says the mesh must close
    // here), with the per-edge matching band.
    let mut closed_edges: Vec<(Vec<Point3>, f64)> = Vec::new();
    for (_eh, edge) in brep.edges.iter() {
        let two_sided = (|| {
            let he = brep.half_edges.get(edge.half_edge)?;
            let twin = brep.half_edges.get(he.twin)?;
            let fa = brep.loops.get(he.loop_)?.face;
            let fb = brep.loops.get(twin.loop_)?.face;
            Some(fa != fb)
        })()
        .unwrap_or(false);
        if !two_sided {
            continue;
        }
        let Some(curve) = edge.curve.as_ref() else {
            continue;
        };
        if !(curve.t_start.is_finite() && curve.t_end.is_finite()) {
            continue;
        }
        let ce: &dyn CurveEval = &curve.geometry;
        let pts: Vec<Point3> = (0..=64)
            .map(|i| ce.point_at(curve.t_start + (curve.t_end - curve.t_start) * i as f64 / 64.0))
            .collect();
        let band = DEFAULT_BAND.max(edge.tolerance.map(|t| 1.5 * t).unwrap_or(0.0));
        closed_edges.push((pts, band));
    }
    if closed_edges.is_empty() {
        return 0;
    }
    let near_tolerant_edge = |p: &Point3| -> Option<f64> {
        let mut best: Option<f64> = None;
        for (pts, band) in &closed_edges {
            let d = pts
                .windows(2)
                .map(|w| dist_point_segment(p, &w[0], &w[1]))
                .fold(f64::MAX, f64::min);
            if d <= *band && best.map(|b| *band > b).unwrap_or(true) {
                best = Some(*band);
            }
        }
        best
    };

    let mut zipped = 0usize;
    for _round in 0..MAX_ROUNDS {
        let splits = zip_round(mesh, &near_tolerant_edge);
        if splits == 0 {
            break;
        }
        zipped += splits;
    }
    zipped
}

fn zip_round(mesh: &mut TriMesh, near_tolerant_edge: &dyn Fn(&Point3) -> Option<f64>) -> usize {
    // Weld census: open edges are used by exactly one triangle.
    let mut weld: HashMap<(i64, i64, i64), usize> = HashMap::new();
    let mut weld_id = |p: Point3| -> usize {
        let k = (
            (p.x * WELD).round() as i64,
            (p.y * WELD).round() as i64,
            (p.z * WELD).round() as i64,
        );
        let next = weld.len();
        *weld.entry(k).or_insert(next)
    };
    let mut ids: Vec<[usize; 3]> = Vec::with_capacity(mesh.triangles.len());
    for t in &mesh.triangles {
        ids.push([
            weld_id(mesh.vertices[t[0]].position),
            weld_id(mesh.vertices[t[1]].position),
            weld_id(mesh.vertices[t[2]].position),
        ]);
    }
    let mut edge_use: HashMap<(usize, usize), usize> = HashMap::new();
    for w in &ids {
        if w[0] == w[1] || w[1] == w[2] || w[0] == w[2] {
            continue;
        }
        for (a, b) in [(w[0], w[1]), (w[1], w[2]), (w[2], w[0])] {
            *edge_use.entry((a.min(b), a.max(b))).or_insert(0) += 1;
        }
    }
    // Open boundary vertices (weld-id -> a representative position).
    let mut open_vertex_pos: HashMap<usize, Point3> = HashMap::new();
    let mut open_edges: Vec<(usize, usize, usize)> = Vec::new(); // (tri, corner, corner+1)
    for (ti, t) in mesh.triangles.iter().enumerate() {
        let w = ids[ti];
        if w[0] == w[1] || w[1] == w[2] || w[0] == w[2] {
            continue;
        }
        for k in 0..3 {
            let (a, b) = (w[k], w[(k + 1) % 3]);
            if edge_use.get(&(a.min(b), a.max(b))) == Some(&1) {
                open_edges.push((ti, k, (k + 1) % 3));
                open_vertex_pos
                    .entry(a)
                    .or_insert(mesh.vertices[t[k]].position);
                open_vertex_pos
                    .entry(b)
                    .or_insert(mesh.vertices[t[(k + 1) % 3]].position);
            }
        }
    }
    if open_edges.is_empty() {
        return 0;
    }
    let candidates: Vec<(usize, Point3)> = open_vertex_pos.into_iter().collect();

    // One split per triangle per round keeps the bookkeeping simple;
    // subsequent rounds finish the job.
    let mut split_tris: std::collections::BTreeSet<usize> = std::collections::BTreeSet::new();
    let mut new_triangles: Vec<([usize; 3], usize, popkern_core::topology::FaceHandle)> =
        Vec::new();
    let mut splits = 0usize;
    let mut stats = (0usize, 0usize, 0usize);
    for (ti, ka, kb) in open_edges {
        stats.0 += 1;
        if split_tris.contains(&ti) {
            continue;
        }
        let tri = mesh.triangles[ti];
        let (ia, ib) = (tri[ka], tri[kb]);
        let (pa, pb) = (mesh.vertices[ia].position, mesh.vertices[ib].position);
        let (wa, wb) = (ids[ti][ka], ids[ti][kb]);
        let seg_len = pa.distance_to(&pb);
        if seg_len < 1e-12 {
            continue;
        }
        // The zip band comes from the recorded tolerance of the BRep
        // edge this junction sits on; no recorded edge, no zip.
        let mid = Point3::new(
            0.5 * (pa.x + pb.x),
            0.5 * (pa.y + pb.y),
            0.5 * (pa.z + pb.z),
        );
        let Some(band) = near_tolerant_edge(&mid) else {
            stats.1 += 1;
            continue;
        };
        // Nearest opposing open vertex strictly interior to the edge.
        let mut best: Option<(f64, usize)> = None;
        for (ci, (wv, pv)) in candidates.iter().enumerate() {
            if *wv == wa || *wv == wb {
                continue;
            }
            let (f, foot) = project_point_segment(pv, &pa, &pb);
            if !(END_FRACTION..=1.0 - END_FRACTION).contains(&f) {
                continue;
            }
            // The candidate must sit within the recorded band of the
            // edge AND split it at a genuinely interior arc position.
            let d = foot.distance_to(pv);
            if d <= band
                && f * seg_len > 1e-9
                && (1.0 - f) * seg_len > 1e-9
                && best.map(|(bd, _)| d < bd).unwrap_or(true)
            {
                best = Some((d, ci));
            }
        }
        let Some((_, ci)) = best else {
            stats.2 += 1;
            continue;
        };
        let (_, pv) = candidates[ci];
        // Split: replace (a, b, c) with (a, v, c) + (v, b, c). The
        // split vertex takes the OPPOSING side's exact position, so
        // the weld pairs them; the displacement is bounded by the
        // recorded tolerance band.
        let kc = 3 - ka - kb;
        let ic = tri[kc];
        let vi = mesh.vertices.len();
        let mut nv = mesh.vertices[ia].clone();
        nv.position = pv;
        // Blend the normal from the split edge's endpoints.
        let nb = mesh.vertices[ib].normal;
        nv.normal = popkern_core::geometry::Vec3::new(
            0.5 * (nv.normal.x + nb.x),
            0.5 * (nv.normal.y + nb.y),
            0.5 * (nv.normal.z + nb.z),
        );
        mesh.vertices.push(nv);
        let face = mesh.face_map[ti];
        // Preserve winding: a->v->c and v->b->c keep the original
        // a->b->c orientation.
        mesh.triangles[ti] = [ia, vi, ic];
        new_triangles.push(([vi, ib, ic], ti, face));
        split_tris.insert(ti);
        splits += 1;
    }
    for (tri, _src, face) in new_triangles {
        mesh.triangles.push(tri);
        mesh.face_map.push(face);
    }
    if std::env::var_os("PK_ZIP_TRACE").is_some() {
        eprintln!(
            "  [zip] open={} no_band={} no_candidate={} splits={splits}",
            stats.0, stats.1, stats.2
        );
    }
    splits
}

fn project_point_segment(p: &Point3, a: &Point3, b: &Point3) -> (f64, Point3) {
    let ab = *b - *a;
    let len2 = ab.dot(&ab);
    if len2 < 1e-300 {
        return (0.0, *a);
    }
    let f = ((*p - *a).dot(&ab) / len2).clamp(0.0, 1.0);
    (
        f,
        Point3::new(a.x + ab.x * f, a.y + ab.y * f, a.z + ab.z * f),
    )
}

fn dist_point_segment(p: &Point3, a: &Point3, b: &Point3) -> f64 {
    project_point_segment(p, a, b).1.distance_to(p)
}
