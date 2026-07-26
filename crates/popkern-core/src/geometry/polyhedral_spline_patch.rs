//! A polyhedral-net spline patch as first-class surface geometry.
//!
//! # Why this is a variant and not just a NURBS surface
//!
//! Be precise about what this buys, because the obvious justification
//! is the wrong one.
//!
//! A PnS bi-4 patch is a bi-quartic Bézier surface, and every Bézier
//! surface is a NURBS surface. Converting one to the other is
//! **geometrically lossless**, and `SurfaceGeometry::Nurbs` already
//! has an analytic second derivative, so curvature analytics work
//! perfectly well on the converted form. **This variant is not an
//! analytics enabler, and must not be justified as one.**
//!
//! What the conversion destroys is *identity*. Once a cap ships as
//! anonymous NURBS:
//!
//! - `unfillet` and `modify` cannot recognise the faces they are
//!   asked to remove or re-radius as cap sectors of one construction;
//!   they see 5 unrelated spline faces.
//! - the sector index and valence — which say *which* piece of *which*
//!   cap a face is — are gone, so a stable handle cannot be rebuilt
//!   after a topology edit.
//! - STEP export writes a B-spline with no record that it was a cap,
//!   and import cannot classify it back.
//! - every consumer that wants the patch pays a conversion round to
//!   get at control points the kernel already had.
//!
//! So the variant carries the BB net plus the provenance, and
//! evaluates the net directly rather than converting first.
//!
//! # Relationship to `Gregory`
//!
//! `SurfaceGeometry::Gregory` is the precedent: it too could be
//! stored as its rational image and is not, for the same reason.

use serde::{Deserialize, Serialize};

use super::NurbsSurface;
use super::surface::SurfaceEval;
use super::{Point3, Vec3};
use crate::geometry::intersection::surface_subdivision::Aabb;

/// Degree of the bi-4 patch in each parameter.
const DEGREE: usize = 4;
/// Control points per row/column (`DEGREE + 1`).
const ORDER: usize = DEGREE + 1;

/// Where a patch came from, so a face can still say what it is after
/// the construction that made it has returned.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub struct PolyhedralSplineProvenance {
    /// Valence of the cap or vertex this patch belongs to — the
    /// number of sectors meeting at its pole.
    pub valence: usize,
    /// Which sector of that cap this patch is, in `0..valence`.
    pub sector: usize,
}

/// A bi-quartic polyhedral-net spline patch.
///
/// The control net is stored as `net[j][i]` matching the cap
/// assembler's own indexing, where `net[0][0]` is the pole corner.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct PolyhedralSplinePatch {
    /// Bézier control net, `net[j][i]`.
    pub net: [[Point3; ORDER]; ORDER],
    /// Which cap and sector this patch is.
    pub provenance: PolyhedralSplineProvenance,
}

impl PolyhedralSplinePatch {
    /// Wrap a cap-assembler patch with its provenance.
    pub fn new(net: [[Point3; ORDER]; ORDER], provenance: PolyhedralSplineProvenance) -> Self {
        Self { net, provenance }
    }

    /// The equivalent NURBS surface.
    ///
    /// Lossless: a bi-4 Bézier patch *is* a NURBS patch with clamped
    /// knots and unit weights. Provided for consumers (STEP export,
    /// interop) that want the generic form; nothing inside the kernel
    /// needs to convert in order to evaluate.
    pub fn to_nurbs(&self) -> NurbsSurface {
        let knots: Vec<f64> = std::iter::repeat_n(0.0, ORDER)
            .chain(std::iter::repeat_n(1.0, ORDER))
            .collect();
        // NurbsSurface indexes `control_points[i_u * num_v + i_v]`,
        // while the cap assembler's net is `net[j][i]` with `i` the
        // u index. Emit u-major so the two agree — flattening the net
        // rows directly transposes the patch, which the
        // direct-vs-converted evaluation gate catches immediately.
        let mut control_points = Vec::with_capacity(ORDER * ORDER);
        for i in 0..ORDER {
            for row in &self.net {
                control_points.push(row[i]);
            }
        }
        NurbsSurface {
            degree_u: DEGREE,
            degree_v: DEGREE,
            knots_u: knots.clone(),
            knots_v: knots,
            control_points,
            weights: vec![1.0; ORDER * ORDER],
            num_u: ORDER,
            num_v: ORDER,
        }
    }

    /// The pole corner — the point all sectors of a cap share.
    pub fn pole(&self) -> Point3 {
        self.net[0][0]
    }
}

/// Bernstein basis of degree 4 at `t`.
fn bernstein4(t: f64) -> [f64; ORDER] {
    let s = 1.0 - t;
    [
        s * s * s * s,
        4.0 * t * s * s * s,
        6.0 * t * t * s * s,
        4.0 * t * t * t * s,
        t * t * t * t,
    ]
}

/// Derivative of the degree-4 Bernstein basis at `t`.
fn bernstein4_deriv(t: f64) -> [f64; ORDER] {
    let s = 1.0 - t;
    // d/dt B_i^4 = 4 (B_{i-1}^3 - B_i^3)
    let b3 = [s * s * s, 3.0 * t * s * s, 3.0 * t * t * s, t * t * t];
    [
        -4.0 * b3[0],
        4.0 * (b3[0] - b3[1]),
        4.0 * (b3[1] - b3[2]),
        4.0 * (b3[2] - b3[3]),
        4.0 * b3[3],
    ]
}

impl SurfaceEval for PolyhedralSplinePatch {
    fn point_at(&self, u: f64, v: f64) -> Point3 {
        let (bu, bv) = (bernstein4(u), bernstein4(v));
        let mut acc = [0.0_f64; 3];
        for (j, row) in self.net.iter().enumerate() {
            for (i, cp) in row.iter().enumerate() {
                let w = bu[i] * bv[j];
                acc[0] += w * cp.x;
                acc[1] += w * cp.y;
                acc[2] += w * cp.z;
            }
        }
        Point3::new(acc[0], acc[1], acc[2])
    }

    fn partial_u(&self, u: f64, v: f64) -> Vec3 {
        let (du, bv) = (bernstein4_deriv(u), bernstein4(v));
        let mut acc = [0.0_f64; 3];
        for (j, row) in self.net.iter().enumerate() {
            for (i, cp) in row.iter().enumerate() {
                let w = du[i] * bv[j];
                acc[0] += w * cp.x;
                acc[1] += w * cp.y;
                acc[2] += w * cp.z;
            }
        }
        Vec3::new(acc[0], acc[1], acc[2])
    }

    fn partial_v(&self, u: f64, v: f64) -> Vec3 {
        let (bu, dv) = (bernstein4(u), bernstein4_deriv(v));
        let mut acc = [0.0_f64; 3];
        for (j, row) in self.net.iter().enumerate() {
            for (i, cp) in row.iter().enumerate() {
                let w = bu[i] * dv[j];
                acc[0] += w * cp.x;
                acc[1] += w * cp.y;
                acc[2] += w * cp.z;
            }
        }
        Vec3::new(acc[0], acc[1], acc[2])
    }

    fn normal_at(&self, u: f64, v: f64) -> Vec3 {
        self.partial_u(u, v).cross(&self.partial_v(u, v))
    }

    fn parameter_range(&self) -> ((f64, f64), (f64, f64)) {
        ((0.0, 1.0), (0.0, 1.0))
    }

    fn closest_point(&self, point: &Point3) -> Option<(f64, f64)> {
        super::surface_helpers::closest_point_gauss_newton(self, point)
    }

    fn bounding_box(&self, u_range: (f64, f64), v_range: (f64, f64)) -> Aabb {
        super::surface_helpers::bounding_box_by_sampling(self, u_range, v_range)
    }
}

#[cfg(test)]
#[path = "polyhedral_spline_patch_tests.rs"]
mod tests;
