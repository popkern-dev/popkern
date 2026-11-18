//! A surface evaluated through a regular chart (design §6a): its
//! parameters are the chart's local coordinates, so a consumer written
//! against `SurfaceEval` (the mesher) works in a chart with no cut and no
//! singular row, unchanged.

use crate::geometry::intersection::surface_subdivision::Aabb;
use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};
use crate::geometry::surface_helpers::bounding_box_by_sampling;
use crate::geometry::{Point3, PortableMath, Vec3};

use super::regular::RegularChart;

/// A remembered inverse: local coordinates and the face parameters they
/// map to (`None` where the inverse failed).
type LastInverse = ((f64, f64), Option<(f64, f64)>);

/// A seed for the inverse: chart `(x, y)` and the face parameters there.
type SeedEntry = ((f64, f64), (f64, f64));

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
    /// Seeds for the inverse, `(x, y)` with the face parameters there, on
    /// a pole projection of a surface with no closed-form projection.
    seeds: std::cell::OnceCell<Vec<SeedEntry>>,
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
            seeds: std::cell::OnceCell::new(),
        }
    }

    /// The inverse's seed at `(x, y)`: the face parameters of the nearest
    /// of a polar table's samples round the singular point (a fixed
    /// table, so the seed and the converged parameters depend on the
    /// point alone). None for an analytic surface, whose projection is
    /// closed-form: unseeded, a NURBS cap's every evaluation was a
    /// multi-start closest point (a paraboloid cup's side meshed 18k
    /// triangles in 10 s).
    fn seed(&self, x: f64, y: f64) -> Option<(f64, f64)> {
        let table = self.seeds.get_or_init(|| self.seed_table());
        table
            .iter()
            .min_by(|a, b| {
                let d = |q: (f64, f64)| (q.0 - x).phypot(q.1 - y);
                d(a.0).total_cmp(&d(b.0))
            })
            .map(|&(_, uv)| uv)
    }

    /// Rings of samples about the singular row, at distances doubling from
    /// 1/4096 of the face's range across it until a ring clears the
    /// domain.
    fn seed_table(&self) -> Vec<SeedEntry> {
        let RegularChart::PoleProjection { row, .. } = &self.chart else {
            return Vec::new();
        };
        if matches!(
            self.surface,
            SurfaceGeometry::Plane { .. }
                | SurfaceGeometry::Cylinder { .. }
                | SurfaceGeometry::Cone { .. }
                | SurfaceGeometry::Sphere { .. }
                | SurfaceGeometry::Torus { .. }
        ) {
            return Vec::new();
        }
        let ((u0, u1), (v0, v1)) = self.surface.parameter_range();
        let ((across_lo, across_hi), (along_lo, along_hi)) = match row.axis {
            super::Axis::V => ((v0, v1), (u0, u1)),
            super::Axis::U => ((u0, u1), (v0, v1)),
        };
        if ![across_lo, across_hi, along_lo, along_hi]
            .iter()
            .all(|t| t.is_finite())
        {
            return Vec::new();
        }
        let span = across_hi - across_lo;
        let inward = if row.at - across_lo <= across_hi - row.at {
            1.0
        } else {
            -1.0
        };
        let reach = (self.domain.0.0.abs())
            .max(self.domain.0.1.abs())
            .max(self.domain.1.0.abs())
            .max(self.domain.1.1.abs());
        const AROUND: usize = 32;
        let mut out = Vec::new();
        let mut step = span / 4096.0;
        while step <= span {
            let across = (row.at + inward * step).clamp(across_lo, across_hi);
            let mut nearest = f64::INFINITY;
            for k in 0..AROUND {
                let along = along_lo + (along_hi - along_lo) * (k as f64 + 0.5) / AROUND as f64;
                let uv = match row.axis {
                    super::Axis::V => (along, across),
                    super::Axis::U => (across, along),
                };
                let xy = self.chart.from_surface(self.surface, uv);
                nearest = nearest.min(xy.0.phypot(xy.1));
                out.push((xy, uv));
            }
            if nearest > reach {
                break;
            }
            step *= 2.0;
        }
        out
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
        let uv = self
            .chart
            .to_surface_near(self.surface, (x, y), self.seed(x, y))
            .ok();
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
