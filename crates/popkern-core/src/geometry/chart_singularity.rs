//! Whether a surface's chart has a singular row: a row of the chart
//! that maps to one point (a sphere's poles, a cone's apex, a spindle
//! torus's axis points, a NURBS boundary row of coincident control
//! points). Beside such a row the chart's speed across it falls to
//! zero, so a curve passing within `d` of the point turns the other
//! coordinate by up to half a period within a stretch of length `d`:
//! uniform samples along the curve cannot follow that sweep however
//! often they double, while on a chart with no singular row the
//! metric's stretch is bounded and doubling converges (PK-1185).

use super::Point3;
use super::nurbs_surface::NurbsSurface;
use super::surface::SurfaceGeometry;

/// Whether `surface`'s chart has a row mapping to one point, rows of a
/// NURBS net counted collapsed when their control points lie within
/// `linear` of one another.
pub fn has_singular_row(surface: &SurfaceGeometry, linear: f64) -> bool {
    match surface {
        SurfaceGeometry::Sphere { .. } | SurfaceGeometry::Cone { .. } => true,
        SurfaceGeometry::Torus {
            major_radius,
            minor_radius,
            ..
        } => minor_radius >= major_radius,
        SurfaceGeometry::Nurbs(n) => nurbs_has_collapsed_row(n, linear),
        SurfaceGeometry::Plane { .. }
        | SurfaceGeometry::Cylinder { .. }
        | SurfaceGeometry::Gregory(_)
        | SurfaceGeometry::PolyhedralSpline(_) => false,
    }
}

/// Whether one of the net's four boundary rows has all its control
/// points within `linear` of its first.
fn nurbs_has_collapsed_row(n: &NurbsSurface, linear: f64) -> bool {
    let (nu, nv) = (n.num_u, n.num_v);
    if nu == 0 || nv == 0 || n.control_points.len() != nu * nv {
        return false;
    }
    let at = |(i, j): (usize, usize)| n.control_points[i * nv + j];
    collapsed((0..nv).map(|j| at((0, j))), linear)
        || collapsed((0..nv).map(|j| at((nu - 1, j))), linear)
        || collapsed((0..nu).map(|i| at((i, 0))), linear)
        || collapsed((0..nu).map(|i| at((i, nv - 1))), linear)
}

/// Whether every point of `row` lies within `linear` of its first.
fn collapsed(mut row: impl Iterator<Item = Point3>, linear: f64) -> bool {
    let Some(first) = row.next() else {
        return false;
    };
    row.all(|p| p.distance_to(&first) <= linear)
}

#[cfg(test)]
#[path = "chart_singularity_tests.rs"]
mod tests;
