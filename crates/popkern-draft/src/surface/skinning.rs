//! NURBS surface skinning from rows of 3D points.
//!
//! Given a grid of 3D points (rows of cross-section profiles), builds a
//! degree-1 NurbsSurface with clamped uniform knots. This is the same
//! approach used in fillet and sweep operations.

use crate::error::DraftError;
use popkern_core::geometry::{NurbsSurface, Point3};

/// Build a NurbsSurface by skinning rows of 3D points.
///
/// Each row is a cross-section profile (same number of points per row).
/// The resulting surface is degree 1 in both u and v directions with
/// clamped uniform knots.
///
/// # Errors
///
/// Returns [`DraftError::Topology`] if the input grid is degenerate
/// (empty rows, inconsistent row lengths, or fewer than 2 rows/columns).
pub fn skin_points(rows: &[Vec<Point3>]) -> Result<NurbsSurface, DraftError> {
    if rows.len() < 2 {
        return Err(DraftError::Topology(
            "skinning requires at least 2 rows".into(),
        ));
    }

    let num_v = rows[0].len();
    if num_v < 2 {
        return Err(DraftError::Topology(
            "skinning requires at least 2 points per row".into(),
        ));
    }

    for (i, row) in rows.iter().enumerate() {
        if row.len() != num_v {
            return Err(DraftError::Topology(format!(
                "row {} has {} points, expected {}",
                i,
                row.len(),
                num_v
            )));
        }
    }

    let num_u = rows.len();
    let degree = 1;

    // Build clamped uniform knots.
    let knots_u = clamped_uniform_knots(num_u, degree);
    let knots_v = clamped_uniform_knots(num_v, degree);

    // Flatten control points row-major.
    let control_points: Vec<Point3> = rows.iter().flat_map(|row| row.iter().copied()).collect();
    let weights = vec![1.0_f64; num_u * num_v];

    Ok(NurbsSurface {
        degree_u: degree,
        degree_v: degree,
        knots_u,
        knots_v,
        control_points,
        num_u,
        num_v,
        weights,
    })
}

/// Build a clamped uniform knot vector for `n` control points and given `degree`.
fn clamped_uniform_knots(n: usize, degree: usize) -> Vec<f64> {
    let num_knots = n + degree + 1;
    let mut knots = Vec::with_capacity(num_knots);
    let n_interior = n - degree; // number of interior spans

    knots.extend(std::iter::repeat_n(0.0, degree + 1));
    for i in 1..n_interior {
        knots.push(i as f64 / n_interior as f64);
    }
    knots.extend(std::iter::repeat_n(1.0, degree + 1));

    knots
}

#[cfg(test)]
mod tests {
    use super::*;
    use popkern_core::geometry::surface::SurfaceEval;

    #[test]
    fn skin_2x2_grid() {
        let rows = vec![
            vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)],
            vec![Point3::new(0.0, 1.0, 0.0), Point3::new(1.0, 1.0, 0.0)],
        ];
        let surf = skin_points(&rows).unwrap();
        assert_eq!(surf.num_u, 2);
        assert_eq!(surf.num_v, 2);
        assert_eq!(surf.degree_u, 1);
        assert_eq!(surf.degree_v, 1);

        // Check corner points.
        let p00 = surf.point_at(0.0, 0.0);
        assert!(p00.distance_to(&Point3::new(0.0, 0.0, 0.0)) < 1e-10);
        let p11 = surf.point_at(1.0, 1.0);
        assert!(p11.distance_to(&Point3::new(1.0, 1.0, 0.0)) < 1e-10);
    }

    #[test]
    fn skin_3x3_grid() {
        let rows = vec![
            vec![
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(0.5, 0.0, 0.0),
                Point3::new(1.0, 0.0, 0.0),
            ],
            vec![
                Point3::new(0.0, 0.5, 0.0),
                Point3::new(0.5, 0.5, 0.1),
                Point3::new(1.0, 0.5, 0.0),
            ],
            vec![
                Point3::new(0.0, 1.0, 0.0),
                Point3::new(0.5, 1.0, 0.0),
                Point3::new(1.0, 1.0, 0.0),
            ],
        ];
        let surf = skin_points(&rows).unwrap();
        assert_eq!(surf.num_u, 3);
        assert_eq!(surf.num_v, 3);
    }

    #[test]
    fn skin_too_few_rows() {
        let rows = vec![vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)]];
        assert!(skin_points(&rows).is_err());
    }

    #[test]
    fn skin_too_few_columns() {
        let rows = vec![
            vec![Point3::new(0.0, 0.0, 0.0)],
            vec![Point3::new(0.0, 1.0, 0.0)],
        ];
        assert!(skin_points(&rows).is_err());
    }

    #[test]
    fn skin_inconsistent_row_lengths() {
        let rows = vec![
            vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)],
            vec![Point3::new(0.0, 1.0, 0.0)],
        ];
        assert!(skin_points(&rows).is_err());
    }

    #[test]
    fn clamped_knots_correct_length() {
        let knots = clamped_uniform_knots(4, 1);
        // n=4, degree=1 => 4+1+1 = 6 knots
        assert_eq!(knots.len(), 6);
        assert_eq!(knots[0], 0.0);
        assert_eq!(knots[1], 0.0);
        assert_eq!(knots[knots.len() - 1], 1.0);
        assert_eq!(knots[knots.len() - 2], 1.0);
    }
}
