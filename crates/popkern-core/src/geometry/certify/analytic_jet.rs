//! Closed-form jet enclosures of the analytic surface kinds.
//!
//! Each formula replicates the kernel's evaluator
//! (`analytic_surfaces.rs`) with the same float frame vectors from
//! [`orthonormal_basis`], so the enclosed surface is the one the
//! kernel evaluates. Only `cos`/`sin`/`tan` are transcendental and
//! carry the libm widening of [`Iv::cos`].

use super::interval::Iv;
use super::interval_vec::IvVec3;
use super::surface_jet::{JET_ORDERS, SurfaceJet};
use crate::geometry::analytic_surfaces::orthonormal_basis;
use crate::geometry::{Point3, Vec3};

/// An analytic surface with its evaluation frame resolved.
#[derive(Debug, Clone, Copy)]
pub enum AnalyticJetSource {
    /// `O + u·U + v·V`.
    Plane {
        /// Origin.
        origin: Point3,
        /// u axis.
        u_axis: Vec3,
        /// v axis.
        v_axis: Vec3,
    },
    /// `O + r(cos u·e1 + sin u·e2) + v·a`.
    Cylinder {
        /// Axis origin.
        origin: Point3,
        /// Axis.
        axis: Vec3,
        /// Frame from [`orthonormal_basis`].
        frame: (Vec3, Vec3),
        /// Radius.
        radius: f64,
    },
    /// `apex + v·a + v·tan(h)(cos u·e1 + sin u·e2)`.
    Cone {
        /// Apex.
        apex: Point3,
        /// Axis.
        axis: Vec3,
        /// Frame from [`orthonormal_basis`].
        frame: (Vec3, Vec3),
        /// Half-angle.
        half_angle: f64,
    },
    /// `C + r(sin v cos u, sin v sin u, cos v)`.
    Sphere {
        /// Center.
        center: Point3,
        /// Radius.
        radius: f64,
    },
    /// `C + (R + r cos v)(cos u·e1 + sin u·e2) + r sin v·a`.
    Torus {
        /// Center.
        center: Point3,
        /// Axis.
        axis: Vec3,
        /// Frame from [`orthonormal_basis`].
        frame: (Vec3, Vec3),
        /// Major radius.
        major: f64,
        /// Minor radius.
        minor: f64,
    },
}

/// Frame of an axis, as the kernel's evaluator builds it.
pub fn frame_of(axis: Vec3) -> (Vec3, Vec3) {
    orthonormal_basis(axis)
}

/// `d^k/dx^k cos x` over `x`, from the enclosures of `cos x`, `sin x`.
fn cos_deriv(k: usize, c: Iv, s: Iv) -> Iv {
    match k % 4 {
        0 => c,
        1 => -s,
        2 => -c,
        _ => s,
    }
}

/// `d^k/dx^k sin x`.
fn sin_deriv(k: usize, c: Iv, s: Iv) -> Iv {
    match k % 4 {
        0 => s,
        1 => c,
        2 => -s,
        _ => -c,
    }
}

fn v3(v: Vec3) -> IvVec3 {
    IvVec3::from_vec(v)
}

fn p3(p: Point3) -> IvVec3 {
    IvVec3::from_point(p)
}

/// `a·e1 + b·e2`.
fn comb(frame: (Vec3, Vec3), a: Iv, b: Iv) -> IvVec3 {
    v3(frame.0).scale(a) + v3(frame.1).scale(b)
}

impl AnalyticJetSource {
    /// Enclosure of every partial up to order 6 over the box.
    pub fn jet(&self, u: Iv, v: Iv) -> SurfaceJet {
        let mut jet = SurfaceJet::zero();
        let (cu, su) = (Iv::cos(u), Iv::sin(u));
        let (cv, sv) = (Iv::cos(v), Iv::sin(v));
        for &(k, l) in JET_ORDERS.iter() {
            let value = match *self {
                AnalyticJetSource::Plane {
                    origin,
                    u_axis,
                    v_axis,
                } => match (k, l) {
                    (0, 0) => p3(origin) + v3(u_axis).scale(u) + v3(v_axis).scale(v),
                    (1, 0) => v3(u_axis),
                    (0, 1) => v3(v_axis),
                    _ => IvVec3::ZERO,
                },
                AnalyticJetSource::Cylinder {
                    origin,
                    axis,
                    frame,
                    radius,
                } => {
                    let r = Iv::point(radius);
                    let ring = comb(frame, r * cos_deriv(k, cu, su), r * sin_deriv(k, cu, su));
                    match (k, l) {
                        (0, 0) => p3(origin) + ring + v3(axis).scale(v),
                        (_, 0) => ring,
                        (0, 1) => v3(axis),
                        _ => IvVec3::ZERO,
                    }
                }
                AnalyticJetSource::Cone {
                    apex,
                    axis,
                    frame,
                    half_angle,
                } => {
                    let t = Iv::tan_point(half_angle);
                    let ring = comb(frame, t * cos_deriv(k, cu, su), t * sin_deriv(k, cu, su));
                    match (k, l) {
                        (0, 0) => p3(apex) + v3(axis).scale(v) + ring.scale(v),
                        (_, 0) => ring.scale(v),
                        (0, 1) => v3(axis) + ring,
                        (_, 1) => ring,
                        _ => IvVec3::ZERO,
                    }
                }
                AnalyticJetSource::Sphere { center, radius } => {
                    let r = Iv::point(radius);
                    let s_l = sin_deriv(l, cv, sv);
                    let z = if k == 0 {
                        r * cos_deriv(l, cv, sv)
                    } else {
                        Iv::point(0.0)
                    };
                    let local = IvVec3::new(
                        r * s_l * cos_deriv(k, cu, su),
                        r * s_l * sin_deriv(k, cu, su),
                        z,
                    );
                    if (k, l) == (0, 0) {
                        p3(center) + local
                    } else {
                        local
                    }
                }
                AnalyticJetSource::Torus {
                    center,
                    axis,
                    frame,
                    major,
                    minor,
                } => {
                    let r = Iv::point(minor);
                    let f = if l == 0 {
                        Iv::point(major) + r * cv
                    } else {
                        r * cos_deriv(l, cv, sv)
                    };
                    let ring = comb(frame, f * cos_deriv(k, cu, su), f * sin_deriv(k, cu, su));
                    let axial = if k == 0 {
                        v3(axis).scale(r * sin_deriv(l, cv, sv))
                    } else {
                        IvVec3::ZERO
                    };
                    let local = ring + axial;
                    if (k, l) == (0, 0) {
                        p3(center) + local
                    } else {
                        local
                    }
                }
            };
            jet.set(k, l, value);
        }
        jet
    }
}
