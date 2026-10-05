use crate::geometry::curve::{CurveEval, closest_point_on_curve, curve_bounding_box_by_sampling};
use crate::geometry::intersection::surface_subdivision::Aabb;
use crate::geometry::{Point3, Vec3};
use popkern_math::PortableMath;

/// Returns an arbitrary vector perpendicular to `normal` (assumed non-zero).
pub(crate) fn arbitrary_perpendicular(normal: Vec3) -> Vec3 {
    // Pick an axis not parallel to normal, then cross to get a perpendicular.
    let axis = if normal.x.abs() < 0.9 {
        Vec3::new(1.0, 0.0, 0.0)
    } else {
        Vec3::new(0.0, 1.0, 0.0)
    };
    normal.cross(&axis)
}

/// Compute the two orthonormal basis vectors (u, v) for a circle/arc in the plane
/// defined by `normal`. `u` and `v` satisfy:
///   - u and v are unit vectors
///   - u, v are perpendicular to each other and to normal
///   - normal ≈ u × v (right-hand rule)
///
/// This is the in-plane frame the [`CurveGeometry::Arc`] / Circle
/// evaluators measure angles in. Distinct from
/// `analytic_surfaces::orthonormal_basis` — callers constructing
/// angle-parameterized arcs MUST use this frame or their angles
/// land in a reflected/rotated chart.
///
/// [`CurveGeometry::Arc`]: crate::geometry::curve::CurveGeometry::Arc
pub fn circle_basis(normal: Vec3) -> (Vec3, Vec3) {
    let perp = arbitrary_perpendicular(normal);
    // u = normalize(perp)
    // SAFETY: arbitrary_perpendicular guarantees a non-zero vector,
    // and normal × u is non-zero when both are non-zero and not parallel.
    // unwrap_or fallback to unit X/Y is defensive but should never trigger.
    let u = perp.normalized().unwrap_or(Vec3::new(1.0, 0.0, 0.0));
    let v = normal
        .cross(&u)
        .normalized()
        .unwrap_or(Vec3::new(0.0, 1.0, 0.0));
    (u, v)
}

// ---------------------------------------------------------------------------
// Line
// ---------------------------------------------------------------------------

pub(super) struct Line {
    pub(super) origin: Point3,
    pub(super) direction: Vec3,
}

impl CurveEval for Line {
    fn point_at(&self, t: f64) -> Point3 {
        self.origin + self.direction * t
    }

    fn derivative_at(&self, _t: f64, order: usize) -> Vec3 {
        if order == 1 {
            self.direction
        } else {
            Vec3::zero()
        }
    }

    fn tangent_at(&self, _t: f64) -> Vec3 {
        self.direction.normalized().unwrap_or(Vec3::zero())
    }

    fn parameter_range(&self) -> (f64, f64) {
        (f64::NEG_INFINITY, f64::INFINITY)
    }

    fn length_between(&self, t0: f64, t1: f64) -> f64 {
        self.direction.length() * (t1 - t0).abs()
    }

    fn closest_point(&self, point: &Point3) -> Option<f64> {
        // Analytic: project point onto line
        let d = self.direction;
        let len_sq = d.dot(&d);
        if len_sq < 1e-30 {
            return None;
        }
        let diff = *point - self.origin;
        Some(diff.dot(&d) / len_sq)
    }

    fn curve_bounding_box(&self, t_range: (f64, f64)) -> Aabb {
        curve_bounding_box_by_sampling(self, t_range)
    }
}

// ---------------------------------------------------------------------------
// Circle
// ---------------------------------------------------------------------------

pub(super) struct Circle {
    pub(super) center: Point3,
    pub(super) radius: f64,
    pub(super) u: Vec3,
    pub(super) v: Vec3,
}

impl Circle {
    pub(super) fn new(center: Point3, normal: Vec3, radius: f64) -> Self {
        let (u, v) = circle_basis(normal);
        Self {
            center,
            radius,
            u,
            v,
        }
    }
}

impl CurveEval for Circle {
    fn point_at(&self, t: f64) -> Point3 {
        self.center + self.u * (self.radius * t.pcos()) + self.v * (self.radius * t.psin())
    }

    fn derivative_at(&self, t: f64, order: usize) -> Vec3 {
        match order % 4 {
            0 => self.u * (self.radius * t.pcos()) + self.v * (self.radius * t.psin()),
            1 => self.u * (-self.radius * t.psin()) + self.v * (self.radius * t.pcos()),
            2 => self.u * (-self.radius * t.pcos()) + self.v * (-self.radius * t.psin()),
            3 => self.u * (self.radius * t.psin()) + self.v * (-self.radius * t.pcos()),
            _ => unreachable!(),
        }
    }

    fn tangent_at(&self, t: f64) -> Vec3 {
        let d = self.derivative_at(t, 1);
        d.normalized().unwrap_or(Vec3::zero())
    }

    fn parameter_range(&self) -> (f64, f64) {
        (0.0, std::f64::consts::TAU)
    }

    fn length_between(&self, t0: f64, t1: f64) -> f64 {
        self.radius * (t1 - t0).abs()
    }

    fn closest_point(&self, point: &Point3) -> Option<f64> {
        closest_point_on_curve(self, point)
    }

    fn curve_bounding_box(&self, t_range: (f64, f64)) -> Aabb {
        curve_bounding_box_by_sampling(self, t_range)
    }
}

// ---------------------------------------------------------------------------
// Arc
// ---------------------------------------------------------------------------

pub(super) struct Arc {
    pub(super) center: Point3,
    pub(super) radius: f64,
    pub(super) start_angle: f64,
    pub(super) end_angle: f64,
    pub(super) u: Vec3,
    pub(super) v: Vec3,
}

impl Arc {
    pub(super) fn new(
        center: Point3,
        normal: Vec3,
        radius: f64,
        start_angle: f64,
        end_angle: f64,
    ) -> Self {
        let (u, v) = circle_basis(normal);
        Self {
            center,
            radius,
            start_angle,
            end_angle,
            u,
            v,
        }
    }
}

impl CurveEval for Arc {
    fn point_at(&self, t: f64) -> Point3 {
        self.center + self.u * (self.radius * t.pcos()) + self.v * (self.radius * t.psin())
    }

    fn derivative_at(&self, t: f64, order: usize) -> Vec3 {
        match order % 4 {
            0 => self.u * (self.radius * t.pcos()) + self.v * (self.radius * t.psin()),
            1 => self.u * (-self.radius * t.psin()) + self.v * (self.radius * t.pcos()),
            2 => self.u * (-self.radius * t.pcos()) + self.v * (-self.radius * t.psin()),
            3 => self.u * (self.radius * t.psin()) + self.v * (-self.radius * t.pcos()),
            _ => unreachable!(),
        }
    }

    fn tangent_at(&self, t: f64) -> Vec3 {
        let d = self.derivative_at(t, 1);
        d.normalized().unwrap_or(Vec3::zero())
    }

    fn parameter_range(&self) -> (f64, f64) {
        (self.start_angle, self.end_angle)
    }

    fn length_between(&self, t0: f64, t1: f64) -> f64 {
        self.radius * (t1 - t0).abs()
    }

    fn closest_point(&self, point: &Point3) -> Option<f64> {
        closest_point_on_curve(self, point)
    }

    fn curve_bounding_box(&self, t_range: (f64, f64)) -> Aabb {
        curve_bounding_box_by_sampling(self, t_range)
    }
}

// ---------------------------------------------------------------------------
// Ellipse — `center + cos(t) · u_basis + sin(t) · v_basis`.
// Used directly by `CurveGeometry::Ellipse` dispatch; the (u_basis,
// v_basis) form keeps both the axis directions and the semi-radii
// in one struct, so non-orthogonal bases work too (the trace is
// still an ellipse for any non-degenerate basis pair).
// ---------------------------------------------------------------------------

pub(super) struct Ellipse {
    pub(super) center: Point3,
    pub(super) u_basis: Vec3,
    pub(super) v_basis: Vec3,
    pub(super) start_angle: f64,
    pub(super) end_angle: f64,
}

impl Ellipse {
    pub(super) fn new(
        center: Point3,
        u_basis: Vec3,
        v_basis: Vec3,
        start_angle: f64,
        end_angle: f64,
    ) -> Self {
        Self {
            center,
            u_basis,
            v_basis,
            start_angle,
            end_angle,
        }
    }
}

impl CurveEval for Ellipse {
    fn point_at(&self, t: f64) -> Point3 {
        self.center + self.u_basis * t.pcos() + self.v_basis * t.psin()
    }

    fn derivative_at(&self, t: f64, order: usize) -> Vec3 {
        match order % 4 {
            0 => self.u_basis * t.pcos() + self.v_basis * t.psin(),
            1 => self.u_basis * (-t.psin()) + self.v_basis * t.pcos(),
            2 => self.u_basis * (-t.pcos()) + self.v_basis * (-t.psin()),
            3 => self.u_basis * t.psin() + self.v_basis * (-t.pcos()),
            _ => unreachable!(),
        }
    }

    fn tangent_at(&self, t: f64) -> Vec3 {
        self.derivative_at(t, 1)
            .normalized()
            .unwrap_or(Vec3::zero())
    }

    fn parameter_range(&self) -> (f64, f64) {
        (self.start_angle, self.end_angle)
    }

    fn length_between(&self, t0: f64, t1: f64) -> f64 {
        // No closed form — sample-and-sum. 32 samples are enough
        // for tessellation-grade length to ~1e-9 over a half-arc.
        let n = 32;
        let mut sum = 0.0;
        let mut prev = self.point_at(t0);
        for i in 1..=n {
            let t = t0 + (t1 - t0) * (i as f64) / (n as f64);
            let p = self.point_at(t);
            sum += prev.distance_to(&p);
            prev = p;
        }
        sum
    }

    fn closest_point(&self, point: &Point3) -> Option<f64> {
        closest_point_on_curve(self, point)
    }

    fn curve_bounding_box(&self, t_range: (f64, f64)) -> Aabb {
        curve_bounding_box_by_sampling(self, t_range)
    }
}

// ---------------------------------------------------------------------------
// Parabola — `vertex + t · linear + t² · quadratic`.
// ---------------------------------------------------------------------------

pub(super) struct Parabola {
    pub(super) vertex: Point3,
    pub(super) linear: Vec3,
    pub(super) quadratic: Vec3,
    pub(super) t_start: f64,
    pub(super) t_end: f64,
}

impl Parabola {
    pub(super) fn new(
        vertex: Point3,
        linear: Vec3,
        quadratic: Vec3,
        t_start: f64,
        t_end: f64,
    ) -> Self {
        Self {
            vertex,
            linear,
            quadratic,
            t_start,
            t_end,
        }
    }
}

impl CurveEval for Parabola {
    fn point_at(&self, t: f64) -> Point3 {
        self.vertex + self.linear * t + self.quadratic * (t * t)
    }

    fn derivative_at(&self, t: f64, order: usize) -> Vec3 {
        // Convention here matches `Ellipse::derivative_at`: order 0 is
        // the curve offset from its base point (i.e. excluding `vertex`).
        match order {
            0 => self.linear * t + self.quadratic * (t * t),
            1 => self.linear + self.quadratic * (2.0 * t),
            2 => self.quadratic * 2.0,
            _ => Vec3::zero(),
        }
    }

    fn tangent_at(&self, t: f64) -> Vec3 {
        self.derivative_at(t, 1)
            .normalized()
            .unwrap_or(Vec3::zero())
    }

    fn parameter_range(&self) -> (f64, f64) {
        (self.t_start, self.t_end)
    }

    fn length_between(&self, t0: f64, t1: f64) -> f64 {
        // Closed-form arc length involves an arctangent, but the
        // 32-sample Riemann sum used elsewhere is sufficient for
        // tessellation-grade accuracy here too.
        let n = 32;
        let mut sum = 0.0;
        let mut prev = self.point_at(t0);
        for i in 1..=n {
            let t = t0 + (t1 - t0) * (i as f64) / (n as f64);
            let p = self.point_at(t);
            sum += prev.distance_to(&p);
            prev = p;
        }
        sum
    }

    fn closest_point(&self, point: &Point3) -> Option<f64> {
        closest_point_on_curve(self, point)
    }

    fn curve_bounding_box(&self, t_range: (f64, f64)) -> Aabb {
        curve_bounding_box_by_sampling(self, t_range)
    }
}

// ---------------------------------------------------------------------------
// Hyperbola — `center + cosh(t) · semi_major + sinh(t) · semi_minor`.
// ---------------------------------------------------------------------------

pub(super) struct Hyperbola {
    pub(super) center: Point3,
    pub(super) semi_major: Vec3,
    pub(super) semi_minor: Vec3,
    pub(super) t_start: f64,
    pub(super) t_end: f64,
}

impl Hyperbola {
    pub(super) fn new(
        center: Point3,
        semi_major: Vec3,
        semi_minor: Vec3,
        t_start: f64,
        t_end: f64,
    ) -> Self {
        Self {
            center,
            semi_major,
            semi_minor,
            t_start,
            t_end,
        }
    }
}

impl CurveEval for Hyperbola {
    fn point_at(&self, t: f64) -> Point3 {
        self.center + self.semi_major * t.pcosh() + self.semi_minor * t.psinh()
    }

    fn derivative_at(&self, t: f64, order: usize) -> Vec3 {
        // Both `cosh` and `sinh` are their own second derivatives;
        // derivatives cycle with period 2: even orders give back the
        // (cosh, sinh) pair, odd orders swap them.
        if order == 0 {
            return self.semi_major * t.pcosh() + self.semi_minor * t.psinh();
        }
        if order % 2 == 1 {
            self.semi_major * t.psinh() + self.semi_minor * t.pcosh()
        } else {
            self.semi_major * t.pcosh() + self.semi_minor * t.psinh()
        }
    }

    fn tangent_at(&self, t: f64) -> Vec3 {
        self.derivative_at(t, 1)
            .normalized()
            .unwrap_or(Vec3::zero())
    }

    fn parameter_range(&self) -> (f64, f64) {
        (self.t_start, self.t_end)
    }

    fn length_between(&self, t0: f64, t1: f64) -> f64 {
        let n = 32;
        let mut sum = 0.0;
        let mut prev = self.point_at(t0);
        for i in 1..=n {
            let t = t0 + (t1 - t0) * (i as f64) / (n as f64);
            let p = self.point_at(t);
            sum += prev.distance_to(&p);
            prev = p;
        }
        sum
    }

    fn closest_point(&self, point: &Point3) -> Option<f64> {
        closest_point_on_curve(self, point)
    }

    fn curve_bounding_box(&self, t_range: (f64, f64)) -> Aabb {
        curve_bounding_box_by_sampling(self, t_range)
    }
}
