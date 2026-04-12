//! Geometric helper functions for constraint checking and projection.
//!
//! Provides point reflection, best-fit plane fitting (PCA), entity
//! representative points, and related linear algebra utilities used
//! by the constraint checking and repair modules.

use crate::geometry::curve::CurveEval;
use crate::geometry::{Point3, Vec3};
use crate::healing::diagnostics::EntityRef;
use crate::topology::BRep;

// ---------------------------------------------------------------------------
// Projection helpers (pub(crate) for use by constraints and constraints_repair)
// ---------------------------------------------------------------------------

/// Project a point across a symmetry plane.
///
/// Given a point and a plane defined by a point and normal, returns
/// the mirror image of the point across the plane.
pub(crate) fn reflect_point(point: &Point3, plane_point: &Point3, plane_normal: &Vec3) -> Point3 {
    let diff = *point - *plane_point;
    let dist = diff.dot(plane_normal);
    let normal = *plane_normal;
    *point - normal * (2.0 * dist)
}

/// Fit a best-fit plane through a set of points (centroid + PCA normal).
///
/// Returns `(centroid, unit_normal)`, or `None` if fewer than 3 points
/// or the points are collinear (degenerate normal).
///
/// The normal is the eigenvector of the covariance matrix with the
/// smallest eigenvalue. We find it by computing the dominant eigenvector
/// and the second eigenvector via deflation, then taking their cross
/// product. Three different initial guesses ensure convergence on
/// axis-aligned data.
pub(crate) fn best_fit_plane(points: &[Point3]) -> Option<(Point3, Vec3)> {
    if points.len() < 3 {
        return None;
    }

    let centroid = compute_centroid(points);

    // Build the 3x3 covariance matrix.
    let mut cov = [[0.0f64; 3]; 3];
    for pt in points {
        let d = [pt.x - centroid.x, pt.y - centroid.y, pt.z - centroid.z];
        for r in 0..3 {
            for c in 0..3 {
                cov[r][c] += d[r] * d[c];
            }
        }
    }

    // Find the dominant eigenvector (largest eigenvalue).
    let v1 = power_iteration_robust(&cov)?;

    // Deflate and find the second eigenvector.
    let deflated = deflate_matrix(&cov, &v1);
    let v2 = power_iteration_robust(&deflated)?;

    // The plane normal is perpendicular to both dominant directions.
    let normal = v1.cross(&v2).normalized()?;
    Some((centroid, normal))
}

/// Get a representative point for an entity.
///
/// For vertices, returns the position. For edges, returns the midpoint
/// of the edge curve (or the origin vertex position). For faces,
/// returns the centroid of the outer loop positions.
pub(crate) fn entity_representative_point(brep: &BRep, entity: &EntityRef) -> Option<Point3> {
    match entity {
        EntityRef::Vertex(vh) => Some(brep.vertices.get(*vh)?.position),
        EntityRef::Edge(eh) => {
            let edge = brep.edges.get(*eh)?;
            if let Some(ref curve) = edge.curve {
                let t_mid = (curve.t_start + curve.t_end) * 0.5;
                Some(curve.geometry.point_at(t_mid))
            } else {
                let he = brep.half_edges.get(edge.half_edge)?;
                Some(brep.vertices.get(he.origin)?.position)
            }
        }
        EntityRef::Face(fh) => {
            let pts = brep.face_outer_loop_positions(*fh);
            if pts.is_empty() {
                None
            } else {
                Some(compute_centroid(&pts))
            }
        }
        EntityRef::HalfEdge(heh) => {
            let he = brep.half_edges.get(*heh)?;
            Some(brep.vertices.get(he.origin)?.position)
        }
        EntityRef::Loop(lh) => {
            let loop_ = brep.loops.get(*lh)?;
            let he = brep.half_edges.get(loop_.half_edge)?;
            Some(brep.vertices.get(he.origin)?.position)
        }
        EntityRef::Shell(_) | EntityRef::EdgePair(_, _) | EntityRef::FacePair(_, _) => None,
    }
}

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

/// Compute the centroid of a set of points.
pub(crate) fn compute_centroid(points: &[Point3]) -> Point3 {
    let n = points.len() as f64;
    if n < 1.0 {
        return Point3::new(0.0, 0.0, 0.0);
    }
    let mut cx = 0.0;
    let mut cy = 0.0;
    let mut cz = 0.0;
    for pt in points {
        cx += pt.x;
        cy += pt.y;
        cz += pt.z;
    }
    Point3::new(cx / n, cy / n, cz / n)
}

/// Clamp infinite parameter range to a bounded midpoint.
pub(crate) fn clamp_param(lo: f64, hi: f64) -> f64 {
    let clamped_lo = lo.max(-100.0);
    let clamped_hi = hi.min(100.0);
    (clamped_lo + clamped_hi) * 0.5
}

/// Power iteration with multiple initial guesses for robustness.
///
/// Tries three axis-aligned initial vectors and returns the eigenvector
/// that produces the largest eigenvalue. This handles axis-aligned data
/// where a single initial guess may be orthogonal to the dominant direction.
fn power_iteration_robust(mat: &[[f64; 3]; 3]) -> Option<Vec3> {
    let starts = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    let mut best_vec: Option<Vec3> = None;
    let mut best_eigenvalue: f64 = -1.0;

    for start in &starts {
        if let Some((vec, eigenvalue)) = power_iteration_single(mat, start)
            && eigenvalue > best_eigenvalue
        {
            best_eigenvalue = eigenvalue;
            best_vec = Some(vec);
        }
    }
    best_vec
}

/// Single power iteration run from a given starting vector.
///
/// Returns the converged eigenvector and its eigenvalue, or `None` if
/// the iteration collapses to zero.
fn power_iteration_single(mat: &[[f64; 3]; 3], start: &[f64; 3]) -> Option<(Vec3, f64)> {
    let mut v = *start;
    for _ in 0..100 {
        let mut nv = [0.0f64; 3];
        for r in 0..3 {
            for c in 0..3 {
                nv[r] += mat[r][c] * v[c];
            }
        }
        let len = (nv[0] * nv[0] + nv[1] * nv[1] + nv[2] * nv[2]).sqrt();
        if len < f64::EPSILON {
            return None;
        }
        v = [nv[0] / len, nv[1] / len, nv[2] / len];
    }
    // Compute eigenvalue: v^T * M * v.
    let mut mv = [0.0f64; 3];
    for r in 0..3 {
        for c in 0..3 {
            mv[r] += mat[r][c] * v[c];
        }
    }
    let eigenvalue = v[0] * mv[0] + v[1] * mv[1] + v[2] * mv[2];
    let vec = Vec3::new(v[0], v[1], v[2]).normalized()?;
    Some((vec, eigenvalue))
}

/// Deflate a 3x3 matrix by removing the component along a vector.
fn deflate_matrix(mat: &[[f64; 3]; 3], v: &Vec3) -> [[f64; 3]; 3] {
    let va = [v.x, v.y, v.z];
    // Compute eigenvalue lambda = v^T * M * v.
    let mut mv = [0.0f64; 3];
    for r in 0..3 {
        for c in 0..3 {
            mv[r] += mat[r][c] * va[c];
        }
    }
    let lambda = va[0] * mv[0] + va[1] * mv[1] + va[2] * mv[2];

    // Deflated matrix: M - lambda * v * v^T.
    let mut result = *mat;
    for r in 0..3 {
        for c in 0..3 {
            result[r][c] -= lambda * va[r] * va[c];
        }
    }
    result
}
