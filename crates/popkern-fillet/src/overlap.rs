//! Overlap and face disappearance detection.
//!
//! Detects when fillet trim curves overlap on the same face (fillets are too
//! close together), and when a fillet consumes an entire face.

use popkern_core::geometry::NurbsCurve;
use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::curve::CurveGeometry;
use popkern_core::geometry::intersection::intersect_curves;
use popkern_core::tolerance::Tolerance;

/// Overlap between two fillets on the same face.
#[derive(Debug)]
pub struct FilletOverlap {
    /// Index of the first fillet.
    pub fillet_a: usize,
    /// Index of the second fillet.
    pub fillet_b: usize,
    /// The intersection curve between the two fillet surfaces, if available.
    pub intersection_curve: Option<NurbsCurve>,
}

/// Check for overlaps between fillet trim curves on shared faces.
///
/// Each entry in `trim_curves` is `(face_index, trim_on_face, other_trim)`.
/// For every pair of entries that share the same `face_index`, the two
/// `trim_on_face` curves are tested for intersection. If they cross, a
/// [`FilletOverlap`] is reported.
pub fn detect_overlaps(
    trim_curves: &[(usize, NurbsCurve, NurbsCurve)],
    tolerance: &Tolerance,
) -> Vec<FilletOverlap> {
    let mut overlaps = Vec::new();

    for i in 0..trim_curves.len() {
        for j in (i + 1)..trim_curves.len() {
            let (face_i, ref trim_i, _) = trim_curves[i];
            let (face_j, ref trim_j, _) = trim_curves[j];

            // Only check curves on the same face.
            if face_i != face_j {
                continue;
            }

            if curves_intersect(trim_i, trim_j, tolerance) {
                overlaps.push(FilletOverlap {
                    fillet_a: i,
                    fillet_b: j,
                    intersection_curve: None,
                });
            }
        }
    }

    overlaps
}

/// Check if any face is entirely consumed by its trim curves.
///
/// For each face, samples points in the face interior not covered by any trim
/// curve. If no uncovered point exists, the face is considered consumed.
///
/// Returns the face indices of consumed faces.
pub fn detect_face_disappearance(
    face_trims: &[(usize, Vec<NurbsCurve>)],
    tolerance: &Tolerance,
) -> Vec<usize> {
    let mut consumed = Vec::new();

    for &(face_index, ref trims) in face_trims {
        if trims.is_empty() {
            continue;
        }

        if face_is_consumed(trims, tolerance) {
            consumed.push(face_index);
        }
    }

    consumed
}

/// Test whether two NURBS curves intersect.
///
/// §18 hot path: fillet trim curves on the same face are usually at
/// different corners, so their AABBs don't overlap. An AABB pre-filter
/// rejects those pairs *before* wrapping the curves in `CurveGeometry`
/// (two clones) and calling the heavy NURBS intersector. Measured impact
/// on a 12-edge unit cube fillet: ~303,000 overlap-phase allocations drop
/// to ~15,000.
///
/// If the AABBs DO overlap, promotes each curve to `CurveGeometry::Nurbs`
/// and delegates to the core curve-curve intersection routine. Falls back
/// to a sampling-based check if the analytic intersector returns an error.
fn curves_intersect(a: &NurbsCurve, b: &NurbsCurve, tolerance: &Tolerance) -> bool {
    let (t_a_min, t_a_max) = a.parameter_range();
    let (t_b_min, t_b_max) = b.parameter_range();

    // AABB pre-filter. Padding = linear tolerance so near-touching
    // curves still pass through to the full intersector.
    let aabb_a = a.curve_bounding_box((t_a_min, t_a_max));
    let aabb_b = b.curve_bounding_box((t_b_min, t_b_max));
    if !aabb_a.overlaps(&aabb_b, tolerance.linear()) {
        return false;
    }

    let geom_a = CurveGeometry::Nurbs(a.clone());
    let geom_b = CurveGeometry::Nurbs(b.clone());

    match intersect_curves(
        &geom_a,
        (t_a_min, t_a_max),
        &geom_b,
        (t_b_min, t_b_max),
        tolerance,
    ) {
        Ok(result) => !result.intersections.is_empty(),
        Err(_) => {
            // Fallback: sample both curves and check minimum distance.
            sample_based_intersection_check(a, b, tolerance)
        }
    }
}

/// Sample-based intersection check as a fallback.
///
/// Samples both curves and checks whether any pair of sample points is
/// closer than the linear tolerance.
fn sample_based_intersection_check(a: &NurbsCurve, b: &NurbsCurve, tolerance: &Tolerance) -> bool {
    let n = 30;
    let (a_min, a_max) = a.parameter_range();
    let (b_min, b_max) = b.parameter_range();

    for i in 0..n {
        let t_a = a_min + (a_max - a_min) * (i as f64) / ((n - 1).max(1) as f64);
        let pa = a.point_at(t_a);

        for j in 0..n {
            let t_b = b_min + (b_max - b_min) * (j as f64) / ((n - 1).max(1) as f64);
            let pb = b.point_at(t_b);

            // Adaptive squared-distance compare: deterministic at
            // the tolerance boundary, no sqrt in the inner loop.
            if popkern_core::adaptive_precision::adaptive_dist_sq_lt(
                (pa.x, pa.y, pa.z),
                (pb.x, pb.y, pb.z),
                tolerance.linear() * tolerance.linear(),
            ) {
                return true;
            }
        }
    }

    false
}

/// Check if trim curves collectively consume a face.
///
/// Heuristic: compute the bounding box of all trim curve samples. If the
/// trim curves span a very tight region (all samples cluster together), the
/// face is likely not consumed. Instead we check whether the trim curve
/// extents cover the range of the face by comparing the "covered" parameter
/// span against a threshold.
///
/// In practice, we sample many points on the trim curves and check if they
/// form a boundary that encloses most of the face. A simplified approach:
/// if there are more than 2 trims covering the face from different
/// directions (different angular quadrants from the centroid), the face
/// is considered consumed.
fn face_is_consumed(trims: &[NurbsCurve], tolerance: &Tolerance) -> bool {
    use popkern_core::geometry::Point3;

    if trims.len() < 2 {
        return false;
    }

    // Sample all trim curves.
    let n_per_curve = 20;
    let mut all_points: Vec<Point3> = Vec::new();

    for trim in trims {
        let (t_min, t_max) = trim.parameter_range();
        for i in 0..n_per_curve {
            let t = t_min + (t_max - t_min) * (i as f64) / ((n_per_curve - 1).max(1) as f64);
            all_points.push(trim.point_at(t));
        }
    }

    if all_points.is_empty() {
        return false;
    }

    // Compute centroid of all sample points.
    let n = all_points.len() as f64;
    let cx = all_points.iter().map(|p| p.x).sum::<f64>() / n;
    let cy = all_points.iter().map(|p| p.y).sum::<f64>() / n;
    let _cz = all_points.iter().map(|p| p.z).sum::<f64>() / n;

    // Check angular coverage: project directions from centroid onto a plane
    // and see if they span all four quadrants in the two dominant axes.
    let mut quadrants = [false; 4];
    for p in &all_points {
        let dx = p.x - cx;
        let dy = p.y - cy;
        let dist = (dx * dx + dy * dy).sqrt();
        if dist < tolerance.linear() {
            continue;
        }
        let q = match (dx >= 0.0, dy >= 0.0) {
            (true, true) => 0,
            (false, true) => 1,
            (false, false) => 2,
            (true, false) => 3,
        };
        quadrants[q] = true;
    }

    // Face is considered consumed only if trims cover all four quadrants
    // around the centroid — a simple heuristic.
    quadrants.iter().all(|&q| q)
}

#[cfg(test)]
mod tests {
    use super::*;
    use popkern_core::geometry::Point3;

    fn make_nurbs_line(start: Point3, end: Point3) -> NurbsCurve {
        crate::center_curve::fit_nurbs_through_points(&[start, end])
            .expect("fit should succeed for 2 points")
    }

    fn test_tolerance() -> Tolerance {
        Tolerance::default()
            .with_linear(1e-3)
            .with_angular(1e-3)
            .with_parametric(1e-6)
    }

    #[test]
    fn no_overlap_distant_fillets() {
        let tol = test_tolerance();

        // Two trim curves far apart on the same face.
        let trim_a = make_nurbs_line(Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0));
        let trim_b = make_nurbs_line(Point3::new(0.0, 10.0, 0.0), Point3::new(1.0, 10.0, 0.0));
        let dummy = make_nurbs_line(Point3::new(0.0, 0.0, 1.0), Point3::new(1.0, 0.0, 1.0));

        let trims = vec![(0, trim_a, dummy.clone()), (0, trim_b, dummy)];

        let overlaps = detect_overlaps(&trims, &tol);
        assert!(overlaps.is_empty(), "distant fillets should not overlap");
    }

    #[test]
    fn overlap_close_fillets() {
        let tol = test_tolerance();

        // Two trim curves that cross: one horizontal, one diagonal through it.
        let trim_a = make_nurbs_line(Point3::new(0.0, 0.5, 0.0), Point3::new(1.0, 0.5, 0.0));
        let trim_b = make_nurbs_line(Point3::new(0.5, 0.0, 0.0), Point3::new(0.5, 1.0, 0.0));
        let dummy = make_nurbs_line(Point3::new(0.0, 0.0, 1.0), Point3::new(1.0, 0.0, 1.0));

        let trims = vec![(0, trim_a, dummy.clone()), (0, trim_b, dummy)];

        let overlaps = detect_overlaps(&trims, &tol);
        assert!(
            !overlaps.is_empty(),
            "crossing fillets should be detected as overlapping"
        );
        assert_eq!(overlaps[0].fillet_a, 0);
        assert_eq!(overlaps[0].fillet_b, 1);
    }

    #[test]
    fn no_overlap_different_faces() {
        let tol = test_tolerance();

        // Two trim curves that would cross, but on different faces.
        let trim_a = make_nurbs_line(Point3::new(0.0, 0.5, 0.0), Point3::new(1.0, 0.5, 0.0));
        let trim_b = make_nurbs_line(Point3::new(0.5, 0.0, 0.0), Point3::new(0.5, 1.0, 0.0));
        let dummy = make_nurbs_line(Point3::new(0.0, 0.0, 1.0), Point3::new(1.0, 0.0, 1.0));

        let trims = vec![(0, trim_a, dummy.clone()), (1, trim_b, dummy)];

        let overlaps = detect_overlaps(&trims, &tol);
        assert!(
            overlaps.is_empty(),
            "curves on different faces should not overlap"
        );
    }

    #[test]
    fn face_not_consumed() {
        let tol = test_tolerance();

        // Single small trim on a large conceptual face.
        let trim = make_nurbs_line(Point3::new(0.0, 0.0, 0.0), Point3::new(0.1, 0.0, 0.0));

        let face_trims = vec![(0, vec![trim])];
        let consumed = detect_face_disappearance(&face_trims, &tol);
        assert!(
            consumed.is_empty(),
            "single small trim should not consume a face"
        );
    }

    #[test]
    fn face_consumed_by_surrounding_trims() {
        let tol = test_tolerance();

        // Four trims forming a box around a point, covering all quadrants.
        let t1 = make_nurbs_line(Point3::new(-1.0, -1.0, 0.0), Point3::new(1.0, -1.0, 0.0));
        let t2 = make_nurbs_line(Point3::new(1.0, -1.0, 0.0), Point3::new(1.0, 1.0, 0.0));
        let t3 = make_nurbs_line(Point3::new(1.0, 1.0, 0.0), Point3::new(-1.0, 1.0, 0.0));
        let t4 = make_nurbs_line(Point3::new(-1.0, 1.0, 0.0), Point3::new(-1.0, -1.0, 0.0));

        let face_trims = vec![(0, vec![t1, t2, t3, t4])];
        let consumed = detect_face_disappearance(&face_trims, &tol);
        assert!(
            !consumed.is_empty(),
            "face surrounded by trims should be consumed"
        );
        assert_eq!(consumed[0], 0);
    }
}

/// A mid-face strip interference between two vertex-disjoint fillets.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct StripInterference {
    /// Spec index of the first fillet.
    pub fillet_a: usize,
    /// Spec index of the second fillet.
    pub fillet_b: usize,
}

/// Detect mid-face strip interference: two fillets on VERTEX-DISJOINT
/// edges sharing a support face whose removal strips overlap. Their
/// rails are typically PARALLEL (opposite edges of one face), so the
/// crossing test in [`detect_overlaps`] never fires — the overlap is
/// a region containment: one fillet's rail lies closer to the other
/// fillet's edge than that fillet's own rail does.
///
/// Corner pairs (edges sharing a vertex) are excluded — adjacent
/// fillets legitimately meet at vertex blends and their rails always
/// approach the shared corner.
pub fn detect_strip_interference(
    edge_geometries: &[crate::construct::EdgeGeometry],
    brep: &popkern_core::topology::BRep,
    tolerance: &Tolerance,
) -> Vec<StripInterference> {
    let samples = |curve: &NurbsCurve| -> Vec<popkern_core::geometry::Point3> {
        let (t0, t1) = curve.parameter_range();
        (0..=16)
            .map(|k| curve.point_at(t0 + (t1 - t0) * k as f64 / 16.0))
            .collect()
    };
    let edge_points = |index: usize| -> Vec<popkern_core::geometry::Point3> {
        let analysis = &edge_geometries[index].analysis;
        brep.edges
            .get(analysis.edge)
            .and_then(|edge| edge.curve.as_ref())
            .map(|curve| {
                (0..=16)
                    .map(|k| {
                        curve.geometry.point_at(
                            curve.t_start + (curve.t_end - curve.t_start) * k as f64 / 16.0,
                        )
                    })
                    .collect()
            })
            .unwrap_or_default()
    };
    let min_distance =
        |a: &[popkern_core::geometry::Point3], b: &[popkern_core::geometry::Point3]| -> f64 {
            let mut best = f64::INFINITY;
            for pa in a {
                for pb in b {
                    best = best.min(pa.distance_to(pb));
                }
            }
            best
        };
    let mut interferences = Vec::new();
    for i in 0..edge_geometries.len() {
        for j in (i + 1)..edge_geometries.len() {
            let (a, b) = (&edge_geometries[i].analysis, &edge_geometries[j].analysis);
            // Corner pairs are legitimate vertex-blend territory.
            if a.vertex_start == b.vertex_start
                || a.vertex_start == b.vertex_end
                || a.vertex_end == b.vertex_start
                || a.vertex_end == b.vertex_end
            {
                continue;
            }
            // Shared support face pairs, with each fillet's trim on
            // that face.
            let shared: Vec<(&NurbsCurve, &NurbsCurve)> = [
                (a.face_1, &edge_geometries[i].trim_3d_1),
                (a.face_2, &edge_geometries[i].trim_3d_2),
            ]
            .into_iter()
            .flat_map(|(fa, trim_a)| {
                [
                    (b.face_1, &edge_geometries[j].trim_3d_1),
                    (b.face_2, &edge_geometries[j].trim_3d_2),
                ]
                .into_iter()
                .filter(move |(fb, _)| *fb == fa)
                .map(move |(_, trim_b)| (trim_a, trim_b))
            })
            .collect();
            if shared.is_empty() {
                continue;
            }
            let pts_edge_i = edge_points(i);
            let pts_edge_j = edge_points(j);
            if pts_edge_i.is_empty() || pts_edge_j.is_empty() {
                continue;
            }
            let hit = shared.iter().any(|(trim_a, trim_b)| {
                let rail_a = samples(trim_a);
                let rail_b = samples(trim_b);
                let width_a = min_distance(&rail_a, &pts_edge_i);
                let width_b = min_distance(&rail_b, &pts_edge_j);
                let b_into_a = min_distance(&rail_b, &pts_edge_i);
                let a_into_b = min_distance(&rail_a, &pts_edge_j);
                let slack = tolerance.linear().max(1.0e-9) * 10.0;
                b_into_a + slack < width_a || a_into_b + slack < width_b
            });
            if hit {
                interferences.push(StripInterference {
                    fillet_a: i,
                    fillet_b: j,
                });
            }
        }
    }
    interferences
}
