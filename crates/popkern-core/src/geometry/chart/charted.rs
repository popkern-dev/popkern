//! A surface evaluated through a regular chart (design §6a): its
//! parameters are the chart's local coordinates, so a consumer written
//! against `SurfaceEval` (the mesher) works in a chart with no cut and no
//! singular row, unchanged.

use crate::geometry::intersection::surface_subdivision::Aabb;
use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};
use crate::geometry::surface_helpers::bounding_box_by_sampling;
use crate::geometry::{Point3, Vec3};

use super::regular::RegularChart;

/// A remembered inverse: local coordinates and the face parameters they
/// map to (`None` where the inverse failed).
type LastInverse = ((f64, f64), Option<(f64, f64)>);

/// `surface` seen through `chart`, over the local box `domain`.
#[derive(Debug, Clone)]
pub struct ChartedSurface<'a> {
    /// The face's surface.
    pub surface: &'a SurfaceGeometry,
    /// The regular chart.
    pub chart: RegularChart,
    /// The local coordinates' box the chart is valid over (a pole
    /// projection's disc's box; a shifted chart's lifted window).
    pub domain: ((f64, f64), (f64, f64)),
    /// The last inverse, `(x, y) → (u, v)`: a consumer evaluates a point
    /// and its normal (or partials) at one place in turn, each needing
    /// the chart's inverse.
    last: std::cell::Cell<Option<LastInverse>>,
}

impl<'a> ChartedSurface<'a> {
    /// `surface` through `chart` over `domain`.
    pub fn new(
        surface: &'a SurfaceGeometry,
        chart: RegularChart,
        domain: ((f64, f64), (f64, f64)),
    ) -> Self {
        ChartedSurface {
            surface,
            chart,
            domain,
            last: std::cell::Cell::new(None),
        }
    }

    /// Face parameters of local coordinates `(x, y)`; `None` where the
    /// chart's inverse fails (a typed `ChartError` from `to_surface`).
    pub fn face_uv(&self, x: f64, y: f64) -> Option<(f64, f64)> {
        if let Some((at, uv)) = self.last.get()
            && at.0.to_bits() == x.to_bits()
            && at.1.to_bits() == y.to_bits()
        {
            return uv;
        }
        let uv = self.chart.to_surface(self.surface, (x, y)).ok();
        self.last.set(Some(((x, y), uv)));
        uv
    }

    /// Both partials with respect to the local coordinates. In a pole
    /// projection the surface is the graph `p = c + x e1 + y e2 + h(x, y) n`
    /// over its tangent plane, so `p_x = e1 + h_x n` with `h_x = −(N·e1) /
    /// (N·n)` for the surface normal `N` at the point: exact, and defined at
    /// the singular point itself (where the face partials vanish and no
    /// Jacobian inverse exists). A shifted chart is the face chart.
    fn local_partials(&self, x: f64, y: f64) -> (Vec3, Vec3) {
        let zero = Vec3::new(0.0, 0.0, 0.0);
        let Some((u, v)) = self.face_uv(x, y) else {
            return (zero, zero);
        };
        match &self.chart {
            RegularChart::PoleProjection { e1, e2, .. } => {
                let n = e1.cross(e2);
                let normal = self.surface.normal_at(u, v);
                let nn = normal.dot(&n);
                // A degenerate normal at the row itself: the graph is
                // tangent to its plane there (h's gradient vanishes).
                if !nn.is_finite() || nn.abs() < 1e-12 {
                    return (*e1, *e2);
                }
                (
                    *e1 - n * (normal.dot(e1) / nn),
                    *e2 - n * (normal.dot(e2) / nn),
                )
            }
            _ => self.surface.first_partials(u, v),
        }
    }
}

impl SurfaceEval for ChartedSurface<'_> {
    fn point_at(&self, x: f64, y: f64) -> Point3 {
        match self.face_uv(x, y) {
            Some((u, v)) => self.surface.point_at(u, v),
            None => match &self.chart {
                // The planar point: off the surface by the inverse's
                // residual, which the caller's deviation checks see.
                RegularChart::PoleProjection { centre, e1, e2, .. } => *centre + *e1 * x + *e2 * y,
                _ => self.surface.point_at(x, y),
            },
        }
    }

    fn partial_u(&self, x: f64, y: f64) -> Vec3 {
        self.local_partials(x, y).0
    }

    fn partial_v(&self, x: f64, y: f64) -> Vec3 {
        self.local_partials(x, y).1
    }

    fn first_partials(&self, x: f64, y: f64) -> (Vec3, Vec3) {
        self.local_partials(x, y)
    }

    fn normal_at(&self, x: f64, y: f64) -> Vec3 {
        match self.face_uv(x, y) {
            Some((u, v)) => self.surface.normal_at(u, v),
            None => match &self.chart {
                RegularChart::PoleProjection { e1, e2, .. } => e1.cross(e2),
                _ => Vec3::new(0.0, 0.0, 1.0),
            },
        }
    }

    fn parameter_range(&self) -> ((f64, f64), (f64, f64)) {
        self.domain
    }

    fn closest_point(&self, point: &Point3) -> Option<(f64, f64)> {
        let uv = self.surface.closest_point(point)?;
        Some(self.chart.from_surface(self.surface, uv))
    }

    fn bounding_box(&self, u_range: (f64, f64), v_range: (f64, f64)) -> Aabb {
        bounding_box_by_sampling(self, u_range, v_range)
    }
}
