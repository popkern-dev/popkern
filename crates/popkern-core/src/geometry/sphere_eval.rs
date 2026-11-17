//! A sphere's chart evaluated with each trigonometric value computed
//! once (PK-1223).
//!
//! The chart is `S(u, v) = c + r·(sin v cos u, sin v sin u, cos v)` in
//! the frame of [`super::analytic_surfaces::sphere_to_world`]. Its point
//! and partials need only `sin u, cos u, sin v, cos v`, but each was
//! written with its own calls (the point five, the partials nine), and
//! the portable trigonometry held ~20% of the tangency corpus. Every
//! expression keeps its operation order, so values are bit-identical to
//! the per-call form.

use popkern_math::PortableMath;

use super::Vec3;
use super::analytic_surfaces::orthonormal_basis;

/// `sin u, cos u, sin v, cos v`.
#[derive(Clone, Copy)]
pub(crate) struct SphereTrig {
    su: f64,
    cu: f64,
    sv: f64,
    cv: f64,
}

impl SphereTrig {
    /// The four values at `(u, v)`.
    pub(crate) fn at(u: f64, v: f64) -> SphereTrig {
        SphereTrig {
            su: u.psin(),
            cu: u.pcos(),
            sv: v.psin(),
            cv: v.pcos(),
        }
    }

    /// The point relative to the centre, in the chart frame.
    pub(crate) fn point(self, radius: f64) -> Vec3 {
        Vec3::new(
            radius * self.sv * self.cu,
            radius * self.sv * self.su,
            radius * self.cv,
        )
    }

    /// `∂S/∂u` in the chart frame.
    pub(crate) fn partial_u(self, radius: f64) -> Vec3 {
        Vec3::new(-radius * self.sv * self.su, radius * self.sv * self.cu, 0.0)
    }

    /// `∂S/∂v` in the chart frame.
    pub(crate) fn partial_v(self, radius: f64) -> Vec3 {
        Vec3::new(
            radius * self.cu * self.cv,
            radius * self.su * self.cv,
            -radius * self.sv,
        )
    }
}

/// Both first partials in world coordinates, the trigonometry and the
/// chart frame computed once.
pub(crate) fn sphere_first_partials(radius: f64, axis: &Vec3, u: f64, v: f64) -> (Vec3, Vec3) {
    let trig = SphereTrig::at(u, v);
    let (du, dv) = (trig.partial_u(radius), trig.partial_v(radius));
    if *axis == Vec3::Z {
        return (du, dv);
    }
    let (e1, e2) = orthonormal_basis(*axis);
    let world = |local: Vec3| e1 * local.x + e2 * local.y + *axis * local.z;
    (world(du), world(dv))
}

#[cfg(test)]
mod tests {
    use super::super::analytic_surfaces::sphere_to_world;
    use super::super::surface::{SurfaceEval, SurfaceGeometry};
    use super::super::{Point3, Vec3};
    use popkern_math::PortableMath;

    fn bits(v: Vec3) -> [u64; 3] {
        [v.x.to_bits(), v.y.to_bits(), v.z.to_bits()]
    }

    /// The per-call formulas the evaluation used before sharing.
    fn reference(r: f64, axis: &Vec3, u: f64, v: f64) -> (Vec3, Vec3, Vec3) {
        let p = Vec3::new(
            r * v.psin() * u.pcos(),
            r * v.psin() * u.psin(),
            r * v.pcos(),
        );
        let du = Vec3::new(-r * v.psin() * u.psin(), r * v.psin() * u.pcos(), 0.0);
        let dv = Vec3::new(
            r * u.pcos() * v.pcos(),
            r * u.psin() * v.pcos(),
            -r * v.psin(),
        );
        (
            sphere_to_world(axis, p),
            sphere_to_world(axis, du),
            sphere_to_world(axis, dv),
        )
    }

    #[test]
    fn shared_trigonometry_is_the_per_call_evaluation_bit_for_bit() {
        let center = Point3::new(0.3, -1.7, 2.2);
        for axis in [Vec3::Z, Vec3::new(0.48, -0.6, 0.64)] {
            for radius in [1.0, 37.25, 1e-3] {
                let s = SurfaceGeometry::Sphere {
                    center,
                    radius,
                    axis,
                };
                for i in 0..=24 {
                    for j in 0..=12 {
                        let u = -1.0 + 0.37 * i as f64;
                        let v = std::f64::consts::PI * j as f64 / 12.0;
                        let (p, du, dv) = reference(radius, &axis, u, v);
                        let (got, want) = (s.point_at(u, v), center + p);
                        assert_eq!(
                            [got.x.to_bits(), got.y.to_bits(), got.z.to_bits()],
                            [want.x.to_bits(), want.y.to_bits(), want.z.to_bits()]
                        );
                        assert_eq!(bits(s.partial_u(u, v)), bits(du));
                        assert_eq!(bits(s.partial_v(u, v)), bits(dv));
                        let (fu, fv) = s.first_partials(u, v);
                        assert_eq!((bits(fu), bits(fv)), (bits(du), bits(dv)));
                        assert_eq!(bits(s.normal_at(u, v)), bits(du.cross(&dv)));
                    }
                }
            }
        }
    }
}
