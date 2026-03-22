//! Starting point detection for surface-surface intersection.

use crate::geometry::intersection::marching::{march_project, MarchSample};
use crate::geometry::intersection::surface_subdivision::{decompose_to_patches, refine_patch, SurfacePatch};
use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};
use crate::geometry::Vec3;
use crate::tolerance::Tolerance;

/// A classified starting point for surface-surface intersection marching.
pub(crate) struct StartingPoint {
    pub sample: MarchSample,
    pub is_tangent: bool,
}

/// Find starting points for surface-surface intersection.
///
/// Uses AABB-based subdivision to narrow down candidate region pairs, then applies
/// Newton 4×4 iteration (via `march_project`) to refine intersection points.
/// Results are deduplicated and classified as transversal or tangent.
pub(crate) fn find_starting_points(
    surface1: &SurfaceGeometry,
    u1_range: (f64, f64),
    v1_range: (f64, f64),
    surface2: &SurfaceGeometry,
    u2_range: (f64, f64),
    v2_range: (f64, f64),
    tolerance: &Tolerance,
) -> Vec<StartingPoint> {
    // Step 1: Promote both surfaces to NURBS over the given parameter ranges.
    let nurbs1 = surface1.to_nurbs_in_range(u1_range, v1_range);
    let nurbs2 = surface2.to_nurbs_in_range(u2_range, v2_range);

    // Step 2: Decompose both into patches at internal knots.
    let patches1 = decompose_to_patches(&nurbs1, u1_range, v1_range);
    let patches2 = decompose_to_patches(&nurbs2, u2_range, v2_range);

    // Step 3 & 4: AABB culling with adaptive refinement.
    // We process a queue of (patch1, patch2, depth) pairs.
    let mut candidate_pairs: Vec<(SurfacePatch, SurfacePatch)> = Vec::new();

    let initial_pairs: Vec<(SurfacePatch, SurfacePatch, u32)> = patches1
        .iter()
        .flat_map(|p1| {
            patches2.iter().filter_map(move |p2| {
                if p1.bbox.overlaps(&p2.bbox, tolerance.linear) {
                    Some((p1.clone(), p2.clone(), 0u32))
                } else {
                    None
                }
            })
        })
        .collect();

    let mut queue = initial_pairs;

    while let Some((patch1, patch2, depth)) = queue.pop() {
        let diag = patch1.bbox.diagonal().max(patch2.bbox.diagonal());
        if depth < 6 && diag > 100.0 * tolerance.linear {
            // Refine both patches and re-cull sub-pairs.
            let sub1 = refine_patch(&patch1);
            let sub2 = refine_patch(&patch2);
            for sp1 in &sub1 {
                for sp2 in &sub2 {
                    if sp1.bbox.overlaps(&sp2.bbox, tolerance.linear) {
                        queue.push((sp1.clone(), sp2.clone(), depth + 1));
                    }
                }
            }
        } else {
            candidate_pairs.push((patch1, patch2));
        }
    }

    // Step 5: For each small overlapping pair, try Newton convergence.
    let mut raw_samples: Vec<MarchSample> = Vec::new();

    for (patch1, patch2) in &candidate_pairs {
        // Use the midpoint of each patch as the initial guess.
        let u1_mid = (patch1.u_range.0 + patch1.u_range.1) / 2.0;
        let v1_mid = (patch1.v_range.0 + patch1.v_range.1) / 2.0;
        let u2_mid = (patch2.u_range.0 + patch2.u_range.1) / 2.0;
        let v2_mid = (patch2.v_range.0 + patch2.v_range.1) / 2.0;

        let target = surface1.point_at(u1_mid, v1_mid);

        // Compute a march direction: cross-product of the two normals gives the
        // tangent direction of the intersection curve. If normals are nearly
        // parallel (tangent case), use an arbitrary non-degenerate direction.
        let n1 = surface1.normal_at(u1_mid, v1_mid);
        let n2 = surface2.normal_at(u2_mid, v2_mid);
        let march_dir = pick_march_dir(n1, n2, surface1, u1_mid, v1_mid);

        if let Some(sample) = march_project(
            surface1, u1_mid, v1_mid,
            surface2, u2_mid, v2_mid,
            &target, &march_dir,
            tolerance,
        ) {
            raw_samples.push(sample);
        }
    }

    // Step 6: Deduplicate points within linear tolerance.
    let unique_samples = deduplicate(raw_samples, tolerance);

    // Step 7: Classify each point as transversal or tangent.
    unique_samples
        .into_iter()
        .map(|sample| {
            let is_tangent = classify_tangent(surface1, surface2, &sample, tolerance);
            StartingPoint { sample, is_tangent }
        })
        .collect()
}

/// Choose a march direction for Newton convergence.
///
/// The natural march direction is N1 × N2 (intersection curve tangent).
/// If the normals are nearly parallel (tangent intersection or coincident surfaces),
/// fall back to the first non-degenerate partial derivative of surface 1.
fn pick_march_dir(
    n1: Vec3,
    n2: Vec3,
    surface1: &SurfaceGeometry,
    u1: f64,
    v1: f64,
) -> Vec3 {
    let cross = n1.cross(&n2);
    if cross.length() > 1e-8 {
        return cross;
    }

    // Fallback: use a surface partial (guaranteed non-degenerate for regular surfaces).
    let du = surface1.partial_u(u1, v1);
    if du.length() > 1e-8 {
        return du;
    }
    let dv = surface1.partial_v(u1, v1);
    if dv.length() > 1e-8 {
        return dv;
    }

    // Last resort: cardinal X axis.
    Vec3::new(1.0, 0.0, 0.0)
}

/// Remove duplicate samples: keep first occurrence, skip any within `tolerance.linear`.
fn deduplicate(samples: Vec<MarchSample>, tolerance: &Tolerance) -> Vec<MarchSample> {
    let mut unique: Vec<MarchSample> = Vec::new();
    'outer: for s in samples {
        for u in &unique {
            if s.point.distance_to(&u.point) < tolerance.linear {
                continue 'outer;
            }
        }
        unique.push(s);
    }
    unique
}

/// Returns true if the two surface normals at the sample point are nearly parallel,
/// indicating a tangent intersection.
fn classify_tangent(
    surface1: &SurfaceGeometry,
    surface2: &SurfaceGeometry,
    sample: &MarchSample,
    tolerance: &Tolerance,
) -> bool {
    let n1 = surface1.normal_at(sample.u1, sample.v1);
    let n2 = surface2.normal_at(sample.u2, sample.v2);

    let len1 = n1.length();
    let len2 = n2.length();

    // Degenerate normals can't be classified; treat as non-tangent.
    if len1 < 1e-14 || len2 < 1e-14 {
        return false;
    }

    let cos_angle = n1.dot(&n2).abs() / (len1 * len2);
    cos_angle > 1.0 - tolerance.angular
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::Point3;

    fn xy_plane() -> SurfaceGeometry {
        SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        }
    }

    fn xz_plane() -> SurfaceGeometry {
        SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 1.0, 0.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 0.0, 1.0),
        }
    }

    #[test]
    fn find_starting_points_two_planes() {
        let s1 = xy_plane();
        let s2 = xz_plane();
        let points = find_starting_points(
            &s1, (-1.0, 1.0), (-1.0, 1.0),
            &s2, (-1.0, 1.0), (-1.0, 1.0),
            &Tolerance::default(),
        );
        assert!(!points.is_empty(), "Should find starting points");
        assert!(!points[0].is_tangent, "XY/XZ intersection should be transversal");
    }

    #[test]
    fn no_starting_points_disjoint() {
        let s1 = SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        };
        let s2 = SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 10.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        };
        let points = find_starting_points(
            &s1, (0.0, 1.0), (0.0, 1.0),
            &s2, (0.0, 1.0), (0.0, 1.0),
            &Tolerance::default(),
        );
        assert!(points.is_empty(), "Disjoint parallel planes should yield no starting points");
    }

    #[test]
    fn tangent_classification() {
        // Two identical planes should produce tangent starting points (or not panic).
        let s = SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        };
        // For identical surfaces Newton may not converge (over-determined / under-constrained).
        // The key requirement is that this does not panic.
        let _ = find_starting_points(
            &s, (0.0, 1.0), (0.0, 1.0),
            &s, (0.0, 1.0), (0.0, 1.0),
            &Tolerance::default(),
        );
    }

    #[test]
    fn classify_tangent_parallel_planes() {
        // Two parallel planes at the same height — normals are identical.
        let s = xy_plane();
        let sample = MarchSample {
            u1: 0.5,
            v1: 0.5,
            u2: 0.5,
            v2: 0.5,
            point: Point3::new(0.5, 0.5, 0.0),
        };
        let tol = Tolerance::default();
        assert!(
            classify_tangent(&s, &s, &sample, &tol),
            "Identical-normal planes should be classified tangent"
        );
    }

    #[test]
    fn classify_transversal_perpendicular_planes() {
        let s1 = xy_plane();
        let s2 = xz_plane();
        let sample = MarchSample {
            u1: 0.0,
            v1: 0.0,
            u2: 0.0,
            v2: 0.0,
            point: Point3::new(0.0, 0.0, 0.0),
        };
        let tol = Tolerance::default();
        assert!(
            !classify_tangent(&s1, &s2, &sample, &tol),
            "Perpendicular planes should be classified transversal"
        );
    }
}
