//! UV-centroid sampling for face classification.
//!
//! Provides a robust method for picking a sample point on a face by computing
//! the UV centroid of the outer loop's pcurves and evaluating the surface at
//! that parameter point.

use popkern_core::{
    geometry::{Point3, Vec3, surface::SurfaceEval},
    topology::{BRep, FaceHandle},
};

use super::face_vertex_positions;

/// Try to pick a sample point using the UV centroid of the face's outer loop.
///
/// Collects UV sample points from the face's outer-loop pcurves, computes their
/// centroid, verifies it lies inside the face boundary using ray-casting, and
/// evaluates the surface at that UV point to get a 3D sample. Returns `None` if
/// the face lacks surface geometry, pcurves, or the UV centroid falls outside
/// the face boundary.
pub(super) fn uv_centroid_sample(brep: &BRep, face_handle: FaceHandle) -> Option<Point3> {
    use popkern_core::geometry::curve::CurveEval;

    let face = &brep.faces[face_handle];
    let surface_on_face = face.surface.as_ref()?;
    let boundary = crate::clip::face_uv_boundary(brep, face_handle)?;

    // Collect UV sample points from outer loop pcurves.
    let mut uv_sum = (0.0f64, 0.0f64);
    let mut n = 0usize;
    let samples_per_edge = 4;
    for (_, pcurve) in &boundary {
        let (t0, t1) = pcurve.parameter_range();
        for i in 0..samples_per_edge {
            let t = t0 + (t1 - t0) * i as f64 / samples_per_edge as f64;
            let p = pcurve.point_at(t);
            uv_sum.0 += p.x;
            uv_sum.1 += p.y;
            n += 1;
        }
    }
    if n == 0 {
        return None;
    }

    let u_mid = uv_sum.0 / n as f64;
    let v_mid = uv_sum.1 / n as f64;
    let uv_pt = Point3::new(u_mid, v_mid, 0.0);

    // Verify the UV centroid is inside the face's outer boundary.
    if !crate::split_inner_loop::point_inside_face_uv(brep, face_handle, uv_pt, &brep.tolerance) {
        return None;
    }

    // For faces with inner loops, check that the centroid doesn't fall inside a hole.
    for &inner_loop_h in &face.inner_loops {
        if inner_loop_uv_contains(brep, inner_loop_h, u_mid, v_mid) {
            return None;
        }
    }

    let sample_3d = surface_on_face.geometry.point_at(u_mid, v_mid);

    // Sanity check: the evaluated 3D point should be close to the face's
    // vertex centroid. If not, the surface parameterization may be stale
    // (e.g. vertices were moved without updating the surface). In that case,
    // fall back to the vertex-based heuristics.
    let positions = face_vertex_positions(brep, face_handle);
    if !positions.is_empty() {
        let cnt = positions.len() as f64;
        let centroid = positions.iter().fold(Vec3::new(0.0, 0.0, 0.0), |acc, p| {
            Vec3::new(acc.x + p.x, acc.y + p.y, acc.z + p.z)
        });
        let centroid_pt = Point3::new(centroid.x / cnt, centroid.y / cnt, centroid.z / cnt);
        // The UV sample should be reasonably close to the vertex centroid.
        // Use a threshold proportional to the face size but not overly generous.
        let max_dist = positions
            .iter()
            .map(|p| p.distance_to(&centroid_pt))
            .fold(0.0f64, f64::max);
        let threshold = (max_dist * 0.5).max(brep.tolerance.linear() * 10.0);
        if sample_3d.distance_to(&centroid_pt) > threshold {
            return None;
        }
    }

    Some(sample_3d)
}

/// Check whether a UV point lies inside an inner loop polygon.
///
/// Samples the inner loop's pcurves to build a UV polygon and uses
/// ray-casting for containment. Returns `false` if pcurves are missing.
fn inner_loop_uv_contains(
    brep: &BRep,
    loop_handle: popkern_core::topology::LoopHandle,
    u: f64,
    v: f64,
) -> bool {
    use popkern_core::geometry::curve::CurveEval;

    let loop_data = match brep.loops.get(loop_handle) {
        Some(l) => l,
        None => return false,
    };
    let start = loop_data.half_edge;
    if !brep.half_edges.contains_key(start) {
        return false;
    }

    let samples_per_edge = 8;
    let mut polygon = Vec::new();
    let mut cursor = start;
    let limit = brep.half_edges.len() + 1;
    for _ in 0..limit {
        let he = match brep.half_edges.get(cursor) {
            Some(h) => h,
            None => return false,
        };
        if let Some(pcurve) = &he.pcurve {
            let (t0, t1) = pcurve.parameter_range();
            let range = t1 - t0;
            if range.abs() > 1e-15 {
                for i in 0..samples_per_edge {
                    let t = t0 + range * (i as f64) / (samples_per_edge as f64);
                    polygon.push(pcurve.point_at(t));
                }
            }
        }
        cursor = he.next;
        if cursor == start {
            break;
        }
    }

    crate::split_inner_loop::ray_cast_contains(&polygon, u, v)
}
