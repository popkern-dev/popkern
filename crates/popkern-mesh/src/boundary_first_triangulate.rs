//! Metric-scaled Delaunay of a boundary-first face chart, filtered to
//! the trim region (split from `boundary_first.rs`; PK-650).

use popkern_math::PortableMath;

use super::boundary_first_cdt::recover_segments;
use super::boundary_first_delaunay::bowyer_watson;
use super::boundary_first_edge_samples;
use super::boundary_first_ring_index::RingIndex;

/// Delaunay of `uv_pts` in the metric-scaled chart (`(u, v) ->
/// (u ku + v kuv, v kv)`: the chart metric's Cholesky factor), constrained to
/// the boundary rings (`ring_ranges`: each ring's index range, closed),
/// filtered to the trim region and wound counter-clockwise in the
/// chart. A ring segment the unconstrained triangulation dropped is
/// recovered by flips when it can be (`boundary_first_cdt`); one that
/// cannot stays missing for the caller's coverage check to split.
///
/// Ring points that COINCIDE — a pinch, where a lobe cut from a ring
/// (`dissolve_pinches`) shares its point with the ring it left — are
/// one triangulation vertex: the triangles and segments use the first
/// index of each coincident group, given back as `alias[i]` for the
/// caller's coverage check. Two copies of one point left the
/// triangulation to keep one, every segment at the other could not be
/// recovered, and the face fell back to a mesh off the canonical edge
/// samples (measured, PK-732: the equal-radius tee's main wall, whose
/// hole lobes touch the lens at the singular points — 893 open mesh
/// edges).
pub(super) fn triangulate_region(
    face_handle: popkern_core::topology::FaceHandle,
    uv_pts: &[(f64, f64)],
    (ku, kuv, kv): (f64, f64, f64),
    ring_index: &RingIndex,
    ring_ranges: &[(usize, usize)],
) -> Option<(Vec<[usize; 3]>, Vec<usize>)> {
    // Delaunay on METRIC-SCALED coordinates: the chart mixes units
    // (radians × millimetres on a cylinder), and the incircle
    // predicate degrades continuously with anisotropy (measured on a
    // full-period wall: aspect ~312 lost 30% of the wedge volume,
    // ~1250 went inside-out). Stored vertices stay in raw UV.
    let alias = coincident_ring_points(uv_pts, ring_ranges);
    let canon: Vec<usize> = (0..uv_pts.len()).filter(|&i| alias[i] == i).collect();
    let mut local = vec![usize::MAX; uv_pts.len()];
    for (k, &i) in canon.iter().enumerate() {
        local[i] = k;
    }
    let bw_pts: Vec<(f64, f64)> = canon
        .iter()
        .map(|&i| (uv_pts[i].0 * ku + uv_pts[i].1 * kuv, uv_pts[i].1 * kv))
        .collect();
    let mut tris = bowyer_watson(&bw_pts)?;
    boundary_first_edge_samples::delaunay_audit(face_handle, &bw_pts, &tris);
    let segments: Vec<(usize, usize)> = ring_ranges
        .iter()
        .flat_map(|&(start, end)| {
            (start..end).map(move |a| (a, start + (a + 1 - start) % (end - start)))
        })
        .map(|(a, b)| (local[alias[a]], local[alias[b]]))
        .filter(|(a, b)| a != b)
        .collect();
    let present: popkern_math::hash::FastHashSet<(usize, usize)> = tris
        .iter()
        .flat_map(|t| (0..3).map(move |k| (t[k].min(t[(k + 1) % 3]), t[k].max(t[(k + 1) % 3]))))
        .collect();
    let missing = segments
        .iter()
        .filter(|&&(a, b)| !present.contains(&(a.min(b), a.max(b))))
        .count();
    if missing > 0 {
        let recovered = recover_segments(&bw_pts, &tris, &segments);
        if popkern_core::trace::trace_enabled("boundary-first") {
            eprintln!(
                "    [bf] face {face_handle:?} constrained recovery of {missing} ring segment(s): {}",
                if recovered.is_some() {
                    "recovered"
                } else {
                    "declined"
                }
            );
        }
        if let Some(recovered) = recovered {
            tris = recovered;
        }
    }
    let mut triangles = Vec::with_capacity(tris.len());
    for t in tris {
        let t = t.map(|k| canon[k]);
        let (a, b, c) = (uv_pts[t[0]], uv_pts[t[1]], uv_pts[t[2]]);
        let centroid = ((a.0 + b.0 + c.0) / 3.0, (a.1 + b.1 + c.1) / 3.0);
        if !ring_index.inside(centroid) {
            continue;
        }
        let signed = (b.0 - a.0) * (c.1 - a.1) - (c.0 - a.0) * (b.1 - a.1);
        if signed.abs() < 1e-16 {
            continue;
        }
        triangles.push(if signed > 0.0 { t } else { [t[0], t[2], t[1]] });
    }
    Some((triangles, alias))
}

/// For each point, the smallest index of the ring points it coincides
/// with (itself when none): ring points closer than a billionth of the
/// chart's extent are one point, closed transitively (a chain of points
/// each within the bound of the next is one point). Interior points
/// alias themselves. Every alias is a root (`alias[alias[i]] ==
/// alias[i]`): a chain merged pairwise once left a point aliased to a
/// point itself re-aliased later, which the triangulation's local index
/// read as unmapped (the generated corpus's NURBS ball union panicked
/// indexing with `usize::MAX`).
fn coincident_ring_points(uv_pts: &[(f64, f64)], ring_ranges: &[(usize, usize)]) -> Vec<usize> {
    let mut alias: Vec<usize> = (0..uv_pts.len()).collect();
    let ring: Vec<usize> = ring_ranges.iter().flat_map(|&(s, e)| s..e).collect();
    let (mut lo, mut hi) = (
        (f64::INFINITY, f64::INFINITY),
        (f64::NEG_INFINITY, f64::NEG_INFINITY),
    );
    for &i in &ring {
        lo = (lo.0.min(uv_pts[i].0), lo.1.min(uv_pts[i].1));
        hi = (hi.0.max(uv_pts[i].0), hi.1.max(uv_pts[i].1));
    }
    let eps = 1e-9 * (hi.0 - lo.0).phypot(hi.1 - lo.1).max(1e-12);
    let mut order = ring;
    order.sort_by(|&a, &b| uv_pts[a].0.total_cmp(&uv_pts[b].0).then(a.cmp(&b)));
    fn root(alias: &mut [usize], mut i: usize) -> usize {
        while alias[i] != i {
            alias[i] = alias[alias[i]];
            i = alias[i];
        }
        i
    }
    for (k, &a) in order.iter().enumerate() {
        for &b in &order[k + 1..] {
            if uv_pts[b].0 - uv_pts[a].0 > eps {
                break;
            }
            if (uv_pts[b].1 - uv_pts[a].1).abs() <= eps {
                let (ra, rb) = (root(&mut alias, a), root(&mut alias, b));
                let (keep, drop) = (ra.min(rb), ra.max(rb));
                alias[drop] = keep;
            }
        }
    }
    for i in 0..alias.len() {
        alias[i] = root(&mut alias, i);
    }
    alias
}

#[cfg(test)]
mod tests {
    use super::coincident_ring_points;

    /// A chain a ~ b ~ c with a and c past the bound of each other, in
    /// an index order that merged b into c after a into b: every alias
    /// must be a root.
    #[test]
    fn coincident_chains_alias_to_one_root() {
        let pts = [(1.6e-9, 0.0), (0.8e-9, 0.0), (0.0, 0.0), (1.0, 1.0)];
        let alias = coincident_ring_points(&pts, &[(0, 4)]);
        assert_eq!(alias, vec![0, 0, 0, 3]);
        for &a in &alias {
            assert_eq!(alias[a], a);
        }
    }
}
