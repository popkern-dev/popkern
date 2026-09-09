//! Sample and branch-point records produced by the SSI march.

use crate::geometry::Point3;

/// A detected branch point where the intersection curve splits.
#[derive(Debug, Clone)]
pub struct BranchPoint {
    /// The 3D location of the branch point.
    pub point: Point3,
    /// Parameters on surface 1.
    pub u1: f64,
    /// Parameters on surface 1.
    pub v1: f64,
    /// The Jacobian condition indicator (sin of normal angle) at this point.
    pub condition: f64,
}

/// A sample point on a surface-surface intersection curve.
#[derive(Debug, Clone)]
pub struct MarchSample {
    /// Parameter on surface 1, u direction.
    pub u1: f64,
    /// Parameter on surface 1, v direction.
    pub v1: f64,
    /// Parameter on surface 2, u direction.
    pub u2: f64,
    /// Parameter on surface 2, v direction.
    pub v2: f64,
    /// The 3D point (equal to surface1.point_at(u1, v1) at convergence).
    pub point: Point3,
}
