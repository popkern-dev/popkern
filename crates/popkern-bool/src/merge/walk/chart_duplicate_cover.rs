//! I5 — duplicate-cover dedup (split from `chart_subdivision` per
//! the file-size rule).

use popkern_core::topology::FaceHandle;

use super::chart_subdivision::{FaceChartSubdivision, FilteredRecord};
use crate::merge::graph::IntersectionGraph;

/// I5 — duplicate-cover dedup: coplanar overlap records trace
/// chart segments that COINCIDE with boundary edges (or with each
/// other). The duplicate cover bounds no region and unbalances
/// every junction it touches (measured on the flush-gusset weld:
/// ±1 in/out at each overlap endpoint, whole faces falling back).
/// Boundary edges outrank records; between records, the lower
/// index wins. Comparison is geometric (endpoint + midpoint UVs,
/// either direction) so it is authoring-order independent.
pub(super) fn dedup_duplicate_covers(
    graph: &IntersectionGraph,
    sub: &mut FaceChartSubdivision,
    overlap_pairs: &[(FaceHandle, FaceHandle)],
) {
    let hes: Vec<usize> = sub.embeddings.keys().copied().collect();
    let mut scale = 0.0_f64;
    for emb in sub.embeddings.values() {
        let (s, e) = (emb.uv[0], emb.uv[emb.uv.len() - 1]);
        scale = scale.max((s.0 - e.0).abs() + (s.1 - e.1).abs());
    }
    let tol = 1e-7 * (1.0 + scale);
    // On a u-periodic chart, runs are not yet anchored to one branch
    // here: a closed rim and the record tracing it can sit a whole
    // period apart in u, so u compares modulo the period (measured,
    // PK-925: a cone's base rim on a coaxial cylinder wall and the
    // record along it were never deduplicated, the face fell back to
    // the 3D-Newell walk and served the rim as its own face).
    let u_period = sub.u_period;
    let near = |a: (f64, f64), b: (f64, f64)| -> bool {
        let du = a.0 - b.0;
        let du = match u_period {
            Some(p) if p > 0.0 => du - p * (du / p).round(),
            _ => du,
        };
        du.abs() < tol && (a.1 - b.1).abs() < tol
    };
    let seg = |h: usize| -> ((f64, f64), (f64, f64), (f64, f64)) {
        let uv = &sub.embeddings[&h].uv;
        (uv[0], uv[uv.len() / 2], uv[uv.len() - 1])
    };
    // The middle SAMPLE is a parameter statement: a record clipped
    // from a closed circle samples its span non-uniformly, and its
    // middle sample sat 5e-3 along the boundary edge it traces from
    // the edge's own (measured, PK-942: a torus halved on its equator
    // under rotation_z(0.75) kept all eight record duplicates and every
    // face fell back to the 3D-Newell walk). Either middle lying on the
    // other's polyline is the same statement made geometrically.
    let on_polyline = |p: (f64, f64), h: usize| -> bool {
        sub.embeddings[&h].uv.windows(2).any(|w| {
            let (a, b) = (w[0], w[1]);
            let pu = match u_period {
                Some(per) if per > 0.0 => p.0 - per * ((p.0 - a.0) / per).round(),
                _ => p.0,
            };
            let (dx, dy) = (b.0 - a.0, b.1 - a.1);
            let l2 = (dx * dx + dy * dy).max(1e-300);
            let t = (((pu - a.0) * dx + (p.1 - a.1) * dy) / l2).clamp(0.0, 1.0);
            (pu - a.0 - dx * t).powi(2) + (p.1 - a.1 - dy * t).powi(2) < tol * tol
        })
    };
    let is_record = |h: usize| -> bool {
        matches!(
            graph.half_edges[h].source,
            crate::merge::graph::EdgeSource::Shared { .. }
        )
    };
    let mut drop: Vec<(usize, usize)> = Vec::new();
    for (i, &a) in hes.iter().enumerate() {
        for &b in hes.iter().skip(i + 1) {
            if graph.half_edges[a].twin == b {
                continue; // twins are the two sides of one edge
            }
            let (sa, ma, ea) = seg(a);
            let (sb, mb, eb) = seg(b);
            let middles = near(ma, mb) || (on_polyline(ma, b) && on_polyline(mb, a));
            let same = near(sa, sb) && near(ea, eb) && middles;
            let reversed = near(sa, eb) && near(ea, sb) && middles;
            if !(same || reversed) {
                continue;
            }
            match (is_record(a), is_record(b)) {
                (true, false) => drop.push((a, b)),
                (false, true) => drop.push((b, a)),
                (true, true) => drop.push((b, a)),
                (false, false) => {} // two boundary edges: not ours to judge
            }
        }
    }
    for (d, k) in drop {
        if sub.embeddings.remove(&d).is_some() {
            sub.report.filtered.push(FilteredRecord::DuplicateCover {
                dropped: d,
                kept: k,
            });
        }
    }
    // PARTIAL covers: overlap records are windowed at different
    // nodes than the boundary edges they trace (a weld trace can
    // span half of one boundary edge and all of the next), so
    // segment-identity misses them. SCOPE: only records authored
    // by a known coplanar overlap pair — unconditional
    // sample-on-polyline dropping was measured to also remove
    // records that legitimately bound flush regions (two bool
    // fixtures regressed). Within the overlap pair the semantics
    // are exact: the record traces IN the shared plane, so where
    // it lies on this face's boundary the boundary already
    // bounds, and the record covers nothing new.
    let from_overlap_pair = |h: usize| -> bool {
        match &graph.half_edges[h].source {
            crate::merge::graph::EdgeSource::Shared { face_a, face_b } => overlap_pairs
                .iter()
                .any(|&(a, b)| *face_a == a && *face_b == b),
            _ => false,
        }
    };
    let boundary_polys: Vec<Vec<(f64, f64)>> = sub
        .embeddings
        .iter()
        .filter(|(h, _)| !is_record(**h))
        .map(|(_, e)| e.uv.clone())
        .collect();
    let on_boundary = |p: (f64, f64)| -> bool {
        for poly in &boundary_polys {
            for w in poly.windows(2) {
                let (a, b) = (w[0], w[1]);
                let (dx, dy) = (b.0 - a.0, b.1 - a.1);
                let l2 = (dx * dx + dy * dy).max(1e-300);
                let t = (((p.0 - a.0) * dx + (p.1 - a.1) * dy) / l2).clamp(0.0, 1.0);
                let d2 = (p.0 - a.0 - dx * t).powi(2) + (p.1 - a.1 - dy * t).powi(2);
                if d2 < tol * tol {
                    return true;
                }
            }
        }
        false
    };
    let overlap_hes: Vec<usize> = sub
        .embeddings
        .keys()
        .copied()
        .filter(|&h| from_overlap_pair(h))
        .collect();
    for h in overlap_hes {
        let covered = sub.embeddings[&h].uv.iter().all(|&p| on_boundary(p));
        if covered && sub.embeddings.remove(&h).is_some() {
            sub.report.filtered.push(FilteredRecord::DuplicateCover {
                dropped: h,
                kept: usize::MAX,
            });
        }
    }
}
