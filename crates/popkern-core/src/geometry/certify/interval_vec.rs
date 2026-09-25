//! Interval vectors in 3D (surface jets) and 2D (pcurve jets).

use super::interval::Iv;
use crate::geometry::{Point3, Vec3};
use std::ops::{Add, Neg, Sub};

/// A 3-vector of intervals: a box enclosing an unknown vector.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct IvVec3 {
    /// x component.
    pub x: Iv,
    /// y component.
    pub y: Iv,
    /// z component.
    pub z: Iv,
}

impl IvVec3 {
    /// The zero vector (exact).
    pub const ZERO: IvVec3 = IvVec3 {
        x: Iv { lo: 0.0, hi: 0.0 },
        y: Iv { lo: 0.0, hi: 0.0 },
        z: Iv { lo: 0.0, hi: 0.0 },
    };

    /// Unbounded in every component.
    pub const ENTIRE: IvVec3 = IvVec3 {
        x: Iv::ENTIRE,
        y: Iv::ENTIRE,
        z: Iv::ENTIRE,
    };

    /// From three intervals.
    pub fn new(x: Iv, y: Iv, z: Iv) -> IvVec3 {
        IvVec3 { x, y, z }
    }

    /// The exact vector with these float components.
    pub fn from_vec(v: Vec3) -> IvVec3 {
        IvVec3::new(Iv::point(v.x), Iv::point(v.y), Iv::point(v.z))
    }

    /// The exact position vector of a point.
    pub fn from_point(p: Point3) -> IvVec3 {
        IvVec3::new(Iv::point(p.x), Iv::point(p.y), Iv::point(p.z))
    }

    /// Scale every component by an interval.
    pub fn scale(self, k: Iv) -> IvVec3 {
        IvVec3::new(self.x * k, self.y * k, self.z * k)
    }

    /// Dot product.
    pub fn dot(self, o: IvVec3) -> Iv {
        self.x * o.x + self.y * o.y + self.z * o.z
    }

    /// Cross product.
    pub fn cross(self, o: IvVec3) -> IvVec3 {
        IvVec3::new(
            self.y * o.z - self.z * o.y,
            self.z * o.x - self.x * o.z,
            self.x * o.y - self.y * o.x,
        )
    }

    /// Squared length, from squares (never negative).
    pub fn norm_sqr(self) -> Iv {
        self.x.sqr() + self.y.sqr() + self.z.sqr()
    }

    /// Length.
    pub fn norm(self) -> Iv {
        self.norm_sqr().sqrt()
    }

    /// Largest length of any vector in the box.
    pub fn mag(self) -> f64 {
        let s = self.x.mag().powi(2) + self.y.mag().powi(2) + self.z.mag().powi(2);
        s.sqrt().next_up().next_up()
    }

    /// Component-wise hull.
    pub fn hull(self, o: IvVec3) -> IvVec3 {
        IvVec3::new(self.x.hull(o.x), self.y.hull(o.y), self.z.hull(o.z))
    }

    /// True when `v` lies in the box.
    pub fn contains(self, v: Vec3) -> bool {
        self.x.contains(v.x) && self.y.contains(v.y) && self.z.contains(v.z)
    }

    /// True when every component is finite.
    pub fn is_finite(self) -> bool {
        self.x.is_finite() && self.y.is_finite() && self.z.is_finite()
    }

    /// The midpoint vector.
    pub fn mid(self) -> Vec3 {
        Vec3::new(self.x.mid(), self.y.mid(), self.z.mid())
    }
}

impl Add for IvVec3 {
    type Output = IvVec3;
    fn add(self, o: IvVec3) -> IvVec3 {
        IvVec3::new(self.x + o.x, self.y + o.y, self.z + o.z)
    }
}

impl Sub for IvVec3 {
    type Output = IvVec3;
    fn sub(self, o: IvVec3) -> IvVec3 {
        IvVec3::new(self.x - o.x, self.y - o.y, self.z - o.z)
    }
}

impl Neg for IvVec3 {
    type Output = IvVec3;
    fn neg(self) -> IvVec3 {
        IvVec3::new(-self.x, -self.y, -self.z)
    }
}

/// A 2-vector of intervals: a `(u, v)` parameter box.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct IvVec2 {
    /// u component.
    pub u: Iv,
    /// v component.
    pub v: Iv,
}

impl IvVec2 {
    /// From two intervals.
    pub fn new(u: Iv, v: Iv) -> IvVec2 {
        IvVec2 { u, v }
    }

    /// Scale both components by an interval.
    pub fn scale(self, k: Iv) -> IvVec2 {
        IvVec2::new(self.u * k, self.v * k)
    }

    /// Component-wise hull.
    pub fn hull(self, o: IvVec2) -> IvVec2 {
        IvVec2::new(self.u.hull(o.u), self.v.hull(o.v))
    }

    /// True when `(u, v)` lies in the box.
    pub fn contains(self, u: f64, v: f64) -> bool {
        self.u.contains(u) && self.v.contains(v)
    }

    /// True when both components are finite.
    pub fn is_finite(self) -> bool {
        self.u.is_finite() && self.v.is_finite()
    }
}

impl Add for IvVec2 {
    type Output = IvVec2;
    fn add(self, o: IvVec2) -> IvVec2 {
        IvVec2::new(self.u + o.u, self.v + o.v)
    }
}

impl Sub for IvVec2 {
    type Output = IvVec2;
    fn sub(self, o: IvVec2) -> IvVec2 {
        IvVec2::new(self.u - o.u, self.v - o.v)
    }
}
