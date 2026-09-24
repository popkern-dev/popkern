//! Global tensor-product B-spline surface interpolation (Piegl &
//! Tiller §9.2.5): interpolate every grid row in `u`, then every
//! column of the resulting control rows in `v`. Because each pass
//! uses one parameter vector, every row shares one knot vector (and
//! every column the other), so the two passes compose into a single
//! non-rational tensor surface that passes through every grid point.

use crate::error::CoreError;
use crate::geometry::Point3;
use crate::geometry::nurbs_surface::NurbsSurface;

/// **Interpolate** a non-rational tensor B-spline surface through a
/// grid of points: `grid[j][i]` is the point at `(u_params[i],
/// v_params[j])`. Parameter vectors must be strictly increasing and
/// carry at least `degree + 1` entries each.
///
/// # Errors
/// `CoreError::PreconditionFailed` for a ragged or undersized grid,
/// mismatched parameter counts, or a singular interpolation system.
pub fn interpolate_surface_grid(
    grid: &[Vec<Point3>],
    u_params: &[f64],
    v_params: &[f64],
    degree_u: usize,
    degree_v: usize,
) -> Result<NurbsSurface, CoreError> {
    if grid.len() != v_params.len() || grid.iter().any(|row| row.len() != u_params.len()) {
        return Err(CoreError::PreconditionFailed(
            "interpolate_surface_grid: grid shape must match the parameter vectors".to_string(),
        ));
    }
    let rows: Vec<crate::geometry::nurbs_curve::NurbsCurve> = grid
        .iter()
        .map(|row| super::interpolate_curve(row, u_params, degree_u))
        .collect::<Result<_, _>>()?;
    let knots_u = rows.first().map(|row| row.knots.clone()).ok_or_else(|| {
        CoreError::PreconditionFailed("interpolate_surface_grid: empty grid".to_string())
    })?;
    let num_u = u_params.len();
    let num_v = v_params.len();
    let mut control_points = vec![Point3::new(0.0, 0.0, 0.0); num_u * num_v];
    let mut knots_v = Vec::new();
    for i in 0..num_u {
        let column: Vec<Point3> = rows.iter().map(|row| row.control_points[i]).collect();
        let curve = super::interpolate_curve(&column, v_params, degree_v)?;
        for (j, point) in curve.control_points.iter().enumerate() {
            control_points[i * num_v + j] = *point;
        }
        knots_v = curve.knots;
    }
    Ok(NurbsSurface {
        degree_u,
        degree_v,
        knots_u,
        knots_v,
        control_points,
        num_u,
        num_v,
        weights: vec![1.0; num_u * num_v],
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::surface::SurfaceEval;

    /// A polynomial surface of degree <= the interpolation degree is
    /// reproduced exactly (interpolation is a projection onto the
    /// spline space, which contains it), including off-grid.
    #[test]
    fn reproduces_a_polynomial_surface_off_grid() {
        let f = |u: f64, v: f64| Point3::new(u, v, u * u * v - 0.5 * v * v * v + u * v);
        let params: Vec<f64> = (0..9).map(|k| k as f64 / 8.0).collect();
        let grid: Vec<Vec<Point3>> = params
            .iter()
            .map(|&v| params.iter().map(|&u| f(u, v)).collect())
            .collect();
        let surface = interpolate_surface_grid(&grid, &params, &params, 3, 3).expect("fit");
        for &(u, v) in &[(0.13, 0.71), (0.5, 0.5), (0.97, 0.03), (0.0, 1.0)] {
            let got = surface.point_at(u, v);
            assert!(
                got.distance_to(&f(u, v)) < 1.0e-12,
                "at ({u}, {v}): {got:?}"
            );
        }
    }

    #[test]
    fn rejects_a_ragged_grid() {
        let params = [0.0, 0.5, 1.0, 1.5];
        let grid = vec![
            vec![Point3::new(0.0, 0.0, 0.0); 4],
            vec![Point3::new(0.0, 0.0, 0.0); 3],
        ];
        assert!(interpolate_surface_grid(&grid, &params, &params[..2], 3, 3).is_err());
    }
}
