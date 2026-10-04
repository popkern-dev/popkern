//! Metric-scaled Delaunay of a boundary-first face chart, filtered to
//! the trim region (split from `boundary_first.rs`; PK-650).

use super::boundary_first_delaunay::bowyer_watson;
use super::boundary_first_edge_samples;
use super::boundary_first_ring_index::RingIndex;

/// Delaunay of `uv_pts` in the metric-scaled chart, filtered to the
/// trim region and wound counter-clockwise in the chart.
pub(super) fn triangulate_region(
    face_handle: popkern_core::topology::FaceHandle,
    uv_pts: &[(f64, f64)],
    (su, sv): (f64, f64),
    ring_index: &RingIndex,
) -> Option<Vec<[usize; 3]>> {
    // Delaunay on METRIC-SCALED coordinates: the chart mixes units
    // (radians × millimetres on a cylinder), and the incircle
    // predicate degrades continuously with anisotropy (measured on a
    // full-period wall: aspect ~312 lost 30% of the wedge volume,
    // ~1250 went inside-out). Stored vertices stay in raw UV.
    let bw_pts: Vec<(f64, f64)> = uv_pts.iter().map(|&(u, v)| (u * su, v * sv)).collect();
    let tris = bowyer_watson(&bw_pts)?;
    boundary_first_edge_samples::delaunay_audit(face_handle, &bw_pts, &tris);
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
