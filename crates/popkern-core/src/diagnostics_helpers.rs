//! Geometric helper functions shared by `diagnostics_detectors.rs`.
//!
//! Split out via `#[path]` so the detector file stays under the
//! 500-line ceiling.

use crate::geometry::Point3;
use crate::geometry::curve::CurveEval;
use crate::geometry::nurbs_curve::NurbsCurve;
use crate::geometry::surface::SurfaceGeometry;
use crate::topology::{BRep, EdgeHandle, FaceHandle};

use super::NurbsKind;

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

pub(super) fn face_bbox(brep: &BRep, fh: FaceHandle) -> (Point3, Point3) {
    let positions = brep.face_outer_loop_positions(fh);
    if positions.is_empty() {
        let z = Point3::new(0.0, 0.0, 0.0);
        return (z, z);
    }
    let mut min = positions[0];
    let mut max = positions[0];
    for p in &positions[1..] {
        if p.x < min.x {
            min.x = p.x;
        }
        if p.y < min.y {
            min.y = p.y;
        }
        if p.z < min.z {
            min.z = p.z;
        }
        if p.x > max.x {
            max.x = p.x;
        }
        if p.y > max.y {
            max.y = p.y;
        }
        if p.z > max.z {
            max.z = p.z;
        }
    }
    (min, max)
}

pub(super) fn aabb_overlap(a: &(Point3, Point3), b: &(Point3, Point3), margin: f64) -> bool {
    !(a.1.x + margin < b.0.x
        || b.1.x + margin < a.0.x
        || a.1.y + margin < b.0.y
        || b.1.y + margin < a.0.y
        || a.1.z + margin < b.0.z
        || b.1.z + margin < a.0.z)
}

/// True when two surface-geometries are the same surface within
/// tolerance. Conservative: only `Plane` and `Cylinder` checked
/// directly; other surface kinds compared by surface-kind discriminant
/// only (so two NURBS faces always count as "potentially equal").
pub(super) fn surfaces_match(a: &SurfaceGeometry, b: &SurfaceGeometry, tol: f64) -> bool {
    use SurfaceGeometry::*;
    match (a, b) {
        (
            Plane {
                origin: o1,
                normal: n1,
                ..
            },
            Plane {
                origin: o2,
                normal: n2,
                ..
            },
        ) => {
            let n_aligned = (n1.x * n2.x + n1.y * n2.y + n1.z * n2.z).abs();
            let n_norm = (n1.x.powi(2) + n1.y.powi(2) + n1.z.powi(2)).sqrt()
                * (n2.x.powi(2) + n2.y.powi(2) + n2.z.powi(2)).sqrt();
            if (n_aligned - n_norm).abs() > tol {
                return false;
            }
            let dx = o2.x - o1.x;
            let dy = o2.y - o1.y;
            let dz = o2.z - o1.z;
            (dx * n1.x + dy * n1.y + dz * n1.z).abs() < tol
        }
        (
            Cylinder {
                origin: o1,
                axis: a1,
                radius: r1,
            },
            Cylinder {
                origin: o2,
                axis: a2,
                radius: r2,
            },
        ) => {
            (r1 - r2).abs() < tol
                && o1.distance_to(o2) < tol
                && (a1.x - a2.x).abs() < tol
                && (a1.y - a2.y).abs() < tol
                && (a1.z - a2.z).abs() < tol
        }
        _ => std::mem::discriminant(a) == std::mem::discriminant(b),
    }
}

/// Validate a `NurbsCurve`'s knot/cp/weight invariants and finite-ness.
pub(super) fn check_nurbs_curve(nurbs: &NurbsCurve) -> Vec<NurbsKind> {
    let mut out = Vec::new();
    let n = nurbs.control_points.len();
    let degree = nurbs.degree;
    if nurbs.weights.len() != n {
        out.push(NurbsKind::CpWeightMismatch {
            cps: n,
            weights: nurbs.weights.len(),
        });
    }
    if nurbs.knots.len() != n + degree + 1 {
        out.push(NurbsKind::CpKnotMismatch {
            cps: n,
            degree,
            knots: nurbs.knots.len(),
        });
    }
    let mut prev = f64::NEG_INFINITY;
    for (i, &k) in nurbs.knots.iter().enumerate() {
        if !k.is_finite() {
            out.push(NurbsKind::NonFiniteKnot { index: i });
            continue;
        }
        if k < prev - 1e-12 {
            out.push(NurbsKind::NonMonotonicKnots { index: i });
            break;
        }
        prev = k;
    }
    let mut run = 1;
    for w in nurbs.knots.windows(2) {
        if (w[1] - w[0]).abs() < 1e-12 {
            run += 1;
            if run > degree + 1 {
                out.push(NurbsKind::HighKnotMultiplicity { degree, run });
                break;
            }
        } else {
            run = 1;
        }
    }
    for (i, &w) in nurbs.weights.iter().enumerate() {
        if !w.is_finite() || w <= 0.0 {
            out.push(NurbsKind::InvalidWeight {
                index: i,
                weight: w,
            });
        }
    }
    for (i, p) in nurbs.control_points.iter().enumerate() {
        if !p.x.is_finite() || !p.y.is_finite() || !p.z.is_finite() {
            out.push(NurbsKind::NonFiniteControlPoint { index: i });
        }
    }
    out
}

/// Validate a `NurbsSurface` CP/knot/weight grid.
pub(super) fn check_nurbs_surface(
    n: &crate::geometry::nurbs_surface::NurbsSurface,
) -> Vec<NurbsKind> {
    let mut out = Vec::new();
    let expected_cps = n.num_u * n.num_v;
    if n.control_points.len() != expected_cps {
        out.push(NurbsKind::CpKnotMismatch {
            cps: n.control_points.len(),
            degree: n.degree_u,
            knots: expected_cps,
        });
    }
    if n.weights.len() != expected_cps {
        out.push(NurbsKind::CpWeightMismatch {
            cps: n.control_points.len(),
            weights: n.weights.len(),
        });
    }
    if n.knots_u.len() != n.num_u + n.degree_u + 1 {
        out.push(NurbsKind::CpKnotMismatch {
            cps: n.num_u,
            degree: n.degree_u,
            knots: n.knots_u.len(),
        });
    }
    if n.knots_v.len() != n.num_v + n.degree_v + 1 {
        out.push(NurbsKind::CpKnotMismatch {
            cps: n.num_v,
            degree: n.degree_v,
            knots: n.knots_v.len(),
        });
    }
    for (i, &k) in n.knots_u.iter().chain(n.knots_v.iter()).enumerate() {
        if !k.is_finite() {
            out.push(NurbsKind::NonFiniteKnot { index: i });
        }
    }
    for (i, &w) in n.weights.iter().enumerate() {
        if !w.is_finite() || w <= 0.0 {
            out.push(NurbsKind::InvalidWeight {
                index: i,
                weight: w,
            });
        }
    }
    for (i, p) in n.control_points.iter().enumerate() {
        if !p.x.is_finite() || !p.y.is_finite() || !p.z.is_finite() {
            out.push(NurbsKind::NonFiniteControlPoint { index: i });
        }
    }
    out
}
