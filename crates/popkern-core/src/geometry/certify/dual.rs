//! Forward-mode interval differentiation along the seam parameter.
//!
//! An [`IvDual`] carries an enclosure of a quantity `q(s)` and of its
//! derivative `q′(s)` over the same s-set. Evaluated over an
//! s-interval `T`, the derivative part encloses `q′(T)`, which is
//! exactly what the mean-value form `q(T) ⊆ q(m) + q′(T)·(T − m)`
//! needs.

use super::interval::Iv;
use super::interval_vec::IvVec3;
use std::ops::{Add, Div, Mul, Neg, Sub};

/// `(q, dq/ds)` as intervals.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct IvDual {
    /// Enclosure of the value.
    pub v: Iv,
    /// Enclosure of the derivative.
    pub d: Iv,
}

impl IvDual {
    /// A value with the given derivative.
    pub fn new(v: Iv, d: Iv) -> IvDual {
        IvDual { v, d }
    }

    /// A constant (zero derivative).
    pub fn constant(v: Iv) -> IvDual {
        IvDual::new(v, Iv::point(0.0))
    }

    /// `√q`; the derivative is [`Iv::ENTIRE`] when the value can
    /// reach 0.
    pub fn sqrt(self) -> IvDual {
        let r = self.v.sqrt();
        IvDual::new(r, self.d.div(r.scale(2.0)))
    }

    /// `q²`.
    pub fn sqr(self) -> IvDual {
        IvDual::new(self.v.sqr(), (self.v * self.d).scale(2.0))
    }
}

impl Add for IvDual {
    type Output = IvDual;
    fn add(self, o: IvDual) -> IvDual {
        IvDual::new(self.v + o.v, self.d + o.d)
    }
}

impl Sub for IvDual {
    type Output = IvDual;
    fn sub(self, o: IvDual) -> IvDual {
        IvDual::new(self.v - o.v, self.d - o.d)
    }
}

impl Neg for IvDual {
    type Output = IvDual;
    fn neg(self) -> IvDual {
        IvDual::new(-self.v, -self.d)
    }
}

impl Mul for IvDual {
    type Output = IvDual;
    fn mul(self, o: IvDual) -> IvDual {
        IvDual::new(self.v * o.v, self.d * o.v + self.v * o.d)
    }
}

/// A 3-vector of [`IvDual`]s.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DualVec3 {
    /// x component.
    pub x: IvDual,
    /// y component.
    pub y: IvDual,
    /// z component.
    pub z: IvDual,
}

impl DualVec3 {
    /// From a value box and a derivative box.
    pub fn new(v: IvVec3, d: IvVec3) -> DualVec3 {
        DualVec3 {
            x: IvDual::new(v.x, d.x),
            y: IvDual::new(v.y, d.y),
            z: IvDual::new(v.z, d.z),
        }
    }

    /// The value box.
    pub fn value(self) -> IvVec3 {
        IvVec3::new(self.x.v, self.y.v, self.z.v)
    }

    /// The derivative box.
    pub fn deriv(self) -> IvVec3 {
        IvVec3::new(self.x.d, self.y.d, self.z.d)
    }

    /// Scale by a dual scalar.
    pub fn scale(self, k: IvDual) -> DualVec3 {
        DualVec3 {
            x: self.x * k,
            y: self.y * k,
            z: self.z * k,
        }
    }

    /// Dot product.
    pub fn dot(self, o: DualVec3) -> IvDual {
        self.x * o.x + self.y * o.y + self.z * o.z
    }

    /// Cross product.
    pub fn cross(self, o: DualVec3) -> DualVec3 {
        DualVec3 {
            x: self.y * o.z - self.z * o.y,
            y: self.z * o.x - self.x * o.z,
            z: self.x * o.y - self.y * o.x,
        }
    }

    /// Squared length (the value part never negative).
    pub fn norm_sqr(self) -> IvDual {
        self.x.sqr() + self.y.sqr() + self.z.sqr()
    }

    /// Length.
    pub fn norm(self) -> IvDual {
        self.norm_sqr().sqrt()
    }
}

impl Add for DualVec3 {
    type Output = DualVec3;
    fn add(self, o: DualVec3) -> DualVec3 {
        DualVec3 {
            x: self.x + o.x,
            y: self.y + o.y,
            z: self.z + o.z,
        }
    }
}

impl Sub for DualVec3 {
    type Output = DualVec3;
    fn sub(self, o: DualVec3) -> DualVec3 {
        DualVec3 {
            x: self.x - o.x,
            y: self.y - o.y,
            z: self.z - o.z,
        }
    }
}

impl Div for IvDual {
    type Output = IvDual;
    /// Quotient `self / o`; the value and derivative are
    /// [`Iv::ENTIRE`] when `o`'s value contains 0.
    fn div(self, o: IvDual) -> IvDual {
        let q = self.v.div(o.v);
        IvDual::new(q, (self.d - q * o.d).div(o.v))
    }
}
