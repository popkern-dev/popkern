//! Metric-scaled Delaunay of a boundary-first face chart, filtered to
//! the trim region (split from `boundary_first.rs`; PK-650).

use super::boundary_first_cdt::recover_segments;
use super::boundary_first_delaunay::bowyer_watson;
use super::boundary_first_edge_samples;
use super::boundary_first_ring_index::RingIndex;

/// Delaunay of `uv_pts` in the metric-scaled chart, constrained to
/// the boundary rings (`ring_ranges`: each ring's index range, closed),
/// filtered to the trim region and wound counter-clockwise in the
/// chart. A ring segment the unconstrained triangulation dropped is
/// recovered by flips when it can be (`boundary_first_cdt`); one that
/// cannot stays missing for the caller's coverage check to split.
pub(super) fn triangulate_region(
    face_handle: popkern_core::topology::FaceHandle,
    uv_pts: &[(f64, f64)],
    (su, sv): (f64, f64),
    ring_index: &RingIndex,
    ring_ranges: &[(usize, usize)],
) -> Option<Vec<[usize; 3]>> {
    // Delaunay on METRIC-SCALED coordinates: the chart mixes units
    // (radians × millimetres on a cylinder), and the incircle
    // predicate degrades continuously with anisotropy (measured on a
    // full-period wall: aspect ~312 lost 30% of the wedge volume,
    // ~1250 went inside-out). Stored vertices stay in raw UV.
    let bw_pts: Vec<(f64, f64)> = uv_pts.iter().map(|&(u, v)| (u * su, v * sv)).collect();
    let mut tris = bowyer_watson(&bw_pts)?;
    boundary_first_edge_samples::delaunay_audit(face_handle, &bw_pts, &tris);
    let segments: Vec<(usize, usize)> = ring_ranges
        .iter()
        .flat_map(|&(start, end)| {
            (start..end).map(move |a| (a, start + (a + 1 - start) % (end - start)))
        })
        .collect();
    let present: std::collections::HashSet<(usize, usize)> = tris
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
    Some(triangles)
}
