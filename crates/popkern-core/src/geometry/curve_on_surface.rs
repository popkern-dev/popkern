//! CurveOnSurface — an intersection curve lying on two surfaces.

use crate::geometry::nurbs_curve::NurbsCurve;
use serde::{Deserialize, Serialize};

/// An intersection curve lying on two surfaces.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CurveOnSurface {
    /// The intersection curve in 3D space.
    pub curve_3d: NurbsCurve,
    /// The curve in surface 1's (u,v) parameter space.
    /// Uses Point3(u, v, 0.0) convention for 2D points.
    pub curve_uv1: NurbsCurve,
    /// The curve in surface 2's (u,v) parameter space.
    pub curve_uv2: NurbsCurve,
    /// Whether the curve forms a closed loop.
    pub is_closed: bool,
}
