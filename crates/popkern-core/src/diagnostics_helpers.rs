//! Geometric helper functions shared by `diagnostics_detectors.rs`.
//!
//! Split out via `#[path]` so the detector file stays under the
//! 500-line ceiling.

use crate::geometry::Point3;
use crate::geometry::curve::CurveEval;
use crate::geometry::surface::SurfaceGeometry;
use crate::topology::{BRep, EdgeHandle, FaceHandle};

pub(super) fn surface_kind_label(sg: &SurfaceGeometry) -> &'static str {
    match sg {
        SurfaceGeometry::Plane { .. } => "Plane",
        SurfaceGeometry::Cylinder { .. } => "Cyl",
        SurfaceGeometry::Sphere { .. } => "Sph",
        SurfaceGeometry::Cone { .. } => "Cone",
        SurfaceGeometry::Torus { .. } => "Torus",
        SurfaceGeometry::Nurbs(_) => "Nurbs",
        SurfaceGeometry::Gregory(_) => "Gregory",
    }
}

pub(super) fn edge_length(brep: &BRep, eh: EdgeHandle) -> f64 {
    let Some(edge) = brep.edges.get(eh) else {
        return 0.0;
    };
    if let Some(curve_on_edge) = edge.curve.as_ref() {
        let len = curve_on_edge
            .geometry
            .length_between(curve_on_edge.t_start, curve_on_edge.t_end);
        if len.is_finite() && len > 0.0 {
            return len;
        }
    }
    let Some(he) = brep.half_edges.get(edge.half_edge) else {
        return 0.0;
    };
    let Some(twin) = brep.half_edges.get(he.twin) else {
        return 0.0;
    };
    let pa = brep.vertices[he.origin].position;
    let pb = brep.vertices[twin.origin].position;
    pa.distance_to(&pb)
}

pub(super) fn edge_midpoint(brep: &BRep, eh: EdgeHandle) -> Point3 {
    let Some(edge) = brep.edges.get(eh) else {
        return Point3::new(0.0, 0.0, 0.0);
    };
    if let Some(curve_on_edge) = edge.curve.as_ref() {
        let t_mid = 0.5 * (curve_on_edge.t_start + curve_on_edge.t_end);
        return curve_on_edge.geometry.point_at(t_mid);
    }
    let Some(he) = brep.half_edges.get(edge.half_edge) else {
        return Point3::new(0.0, 0.0, 0.0);
    };
    let Some(twin) = brep.half_edges.get(he.twin) else {
        return Point3::new(0.0, 0.0, 0.0);
    };
    let pa = brep.vertices[he.origin].position;
    let pb = brep.vertices[twin.origin].position;
    Point3::new(
        0.5 * (pa.x + pb.x),
        0.5 * (pa.y + pb.y),
        0.5 * (pa.z + pb.z),
    )
}

/// Approximate face area by fan-triangulating outer-loop vertex
/// positions and summing unsigned triangle areas. Cheap defect proxy
/// — not exact for curved surfaces.
pub(super) fn face_area(brep: &BRep, fh: FaceHandle) -> f64 {
    let positions = brep.face_outer_loop_positions(fh);
    if positions.len() < 3 {
        return 0.0;
    }
    let p0 = positions[0];
    let mut area = 0.0;
    for i in 1..positions.len() - 1 {
        let p1 = positions[i];
        let p2 = positions[i + 1];
        let ax = p1.x - p0.x;
        let ay = p1.y - p0.y;
        let az = p1.z - p0.z;
        let bx = p2.x - p0.x;
        let by = p2.y - p0.y;
        let bz = p2.z - p0.z;
        let cx = ay * bz - az * by;
        let cy = az * bx - ax * bz;
        let cz = ax * by - ay * bx;
        area += 0.5 * (cx * cx + cy * cy + cz * cz).sqrt();
    }
    area
}

pub(super) fn face_edge_aspect_ratio(brep: &BRep, fh: FaceHandle) -> f64 {
    let positions = brep.face_outer_loop_positions(fh);
    if positions.len() < 2 {
        return 0.0;
    }
    let mut min_len = f64::INFINITY;
    let mut max_len = 0.0_f64;
    for i in 0..positions.len() {
        let pa = positions[i];
        let pb = positions[(i + 1) % positions.len()];
        let len = pa.distance_to(&pb);
        if len > 0.0 && len < min_len {
            min_len = len;
        }
        if len > max_len {
            max_len = len;
        }
    }
    if min_len <= 0.0 || !min_len.is_finite() {
        return f64::INFINITY;
    }
    max_len / min_len
}

pub(super) fn face_centroid(brep: &BRep, fh: FaceHandle) -> Point3 {
    let positions = brep.face_outer_loop_positions(fh);
    if positions.is_empty() {
        return Point3::new(0.0, 0.0, 0.0);
    }
    let n = positions.len() as f64;
    let (mut x, mut y, mut z) = (0.0_f64, 0.0_f64, 0.0_f64);
    for p in &positions {
        x += p.x;
        y += p.y;
        z += p.z;
    }
    Point3::new(x / n, y / n, z / n)
}

pub(super) fn min_median_max(values: &[f64]) -> (f64, f64, f64) {
    let mut sorted: Vec<f64> = values.iter().copied().filter(|v| v.is_finite()).collect();
    if sorted.is_empty() {
        return (0.0, 0.0, 0.0);
    }
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let min = sorted[0];
    let max = sorted[sorted.len() - 1];
    let median = sorted[sorted.len() / 2];
    (min, median, max)
}
