//! Inner-loop detection helpers for the split stage.
//!
//! When a closed intersection curve has 0 boundary crossings and lies
//! entirely inside a face, it should become an inner loop (hole) rather
//! than being dropped.

use popkern_core::{
    geometry::Point3,
    tolerance::Tolerance,
    topology::{BRep, FaceHandle},
};

use crate::clip::face_uv_boundary;

/// Test whether a 2D point (using x,y of a `Point3`) is inside a polygon
/// via the even-odd ray-casting rule.
///
/// Casts a horizontal ray from `(px, py)` in the +x direction and counts
/// crossings with the polygon edges.
pub(crate) fn ray_cast_contains(polygon: &[Point3], px: f64, py: f64) -> bool {
    let n = polygon.len();
    if n < 3 {
        return false;
    }
    let mut inside = false;
    let mut j = n - 1;
    for i in 0..n {
        let yi = polygon[i].y;
        let yj = polygon[j].y;
        if (yi > py) != (yj > py) {
            let xi = polygon[i].x;
            let xj = polygon[j].x;
            let x_intersect = xi + (py - yi) / (yj - yi) * (xj - xi);
            if px < x_intersect {
                inside = !inside;
            }
        }
        j = i;
    }
    inside
}

/// Check whether a UV-space point lies inside a face's outer boundary.
///
/// Samples the face boundary pcurves to build a polygon, then applies
/// even-odd ray casting for containment.
pub(crate) fn point_inside_face_uv(
    brep: &BRep,
    face: FaceHandle,
    uv_point: Point3,
    tolerance: &Tolerance,
) -> bool {
    let boundary = match face_uv_boundary(brep, face) {
        Some(b) => b,
        None => return false,
    };
    // Sample each boundary pcurve to build a polygon.
    let samples_per_edge = 8;
    let mut polygon = Vec::new();
    for (_he, pcurve) in &boundary {
        let (t0, t1) = pcurve.parameter_range();
        let range = t1 - t0;
        if range.abs() < tolerance.linear() {
            continue;
        }
        for i in 0..samples_per_edge {
            let t = t0 + range * (i as f64) / (samples_per_edge as f64);
            polygon.push(pcurve.point_at(t));
        }
    }
    ray_cast_contains(&polygon, uv_point.x, uv_point.y)
}

#[cfg(test)]
#[path = "split_inner_loop_tests.rs"]
mod tests;
