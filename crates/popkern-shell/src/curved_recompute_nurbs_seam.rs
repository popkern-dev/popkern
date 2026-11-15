//! A NURBS face's seam on its offset surface (PK-1215).
//!
//! A face closed round its own chart (a NURBS ball, a revolved NURBS
//! wall) has a seam edge on a boundary row of the chart: the row where
//! the chart closes on itself. Its offset surface keeps the chart, so
//! the offset seam is that surface's own iso-curve on the same row,
//! exact. The parallel map declined on it (a ball's seam joins two poles
//! and maps to no similar curve), and the seam stayed at the old radius,
//! the offset ball open along it by the distance.
//!
//! The row is found from the edge's own ends, already recomputed onto
//! the offset surface: the boundary iso-curve whose two ends are the
//! edge's, either way round, within the linear tolerance.

use popkern_core::geometry::curve::{CurveEval, CurveGeometry, CurveOnEdge};
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::{NurbsCurve, Point3};
use popkern_core::tolerance::Tolerance;

/// The seam on the offset NURBS `surface` from `start` to `end`, when
/// they are the two ends of one of its boundary iso-curves.
pub(crate) fn nurbs_seam_curve(
    surface: &SurfaceGeometry,
    start: Point3,
    end: Point3,
    tolerance: &Tolerance,
) -> Option<CurveOnEdge> {
    let SurfaceGeometry::Nurbs(n) = surface else {
        return None;
    };
    let lin = tolerance.linear();
    let ((u0, u1), (v0, v1)) = n.parameter_range();
    let rows: [NurbsCurve; 4] = [
        n.iso_curve_at_u(u0),
        n.iso_curve_at_u(u1),
        n.iso_curve_at_v(v0),
        n.iso_curve_at_v(v1),
    ];
    rows.into_iter().find_map(|curve| {
        let (t0, t1) = curve.parameter_range();
        let (a, b) = (curve.point_at(t0), curve.point_at(t1));
        // A collapsed row (a pole) is a point, no seam.
        if a.distance_to(&b) <= lin && curve.point_at(0.5 * (t0 + t1)).distance_to(&a) <= lin {
            return None;
        }
        let (t_start, t_end) = if a.distance_to(&start) <= lin && b.distance_to(&end) <= lin {
            (t0, t1)
        } else if b.distance_to(&start) <= lin && a.distance_to(&end) <= lin {
            (t1, t0)
        } else {
            return None;
        };
        Some(CurveOnEdge {
            geometry: CurveGeometry::Nurbs(curve),
            t_start,
            t_end,
        })
    })
}
