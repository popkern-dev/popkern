//! Input descriptors: anchors, hints, surface-side identifiers.

use serde::{Deserialize, Serialize};

use crate::geometry::nurbs_curve::NurbsCurve;
use crate::geometry::{Point3, Vec3};

/// Endpoint anchors for the edge. The factory hits these exactly
/// (G1) via CP equality on the analytic dispatch paths; marching
/// dispatch snaps endpoints to anchors as a final step.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub enum EdgeAnchors {
    /// Open edge with distinct start/end points.
    Open {
        /// Start anchor point — `curve_3d.point_at(t_start)`.
        start: Point3,
        /// End anchor point — `curve_3d.point_at(t_end)`.
        end: Point3,
    },
    /// Closed-loop edge. `sample` lies on the loop; `direction`
    /// disambiguates which traversal direction is canonical.
    Closed {
        /// A point on the loop used to disambiguate the start.
        sample: Point3,
        /// Tangent direction at `sample` defining traversal sense.
        direction: Vec3,
    },
    /// Self-loop / seam: surf_a == surf_b on a degenerate seam.
    Seam {
        /// A point on the seam.
        sample: Point3,
        /// Which UV axis is degenerate.
        seam_axis: SeamAxis,
    },
}

/// Which UV axis is the degenerate seam coordinate (i.e. which
/// dimension is identified at u=0 with u=2π or analogous).
#[derive(Copy, Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum SeamAxis {
    /// Seam runs along constant `u` (e.g. cylinder).
    U,
    /// Seam runs along constant `v` (e.g. some torus
    /// parameterisations).
    V,
}

/// Optional caller-side hints to drive multi-branch disambiguation
/// or to pick the Gregory-boundary specialisation.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub enum EdgeHint {
    /// Caller has a 3-D curve they want preserved. Used to
    /// disambiguate when intersection returns multiple branches
    /// (sphere×cyl Viviani, cyl×cyl parallel ellipses).
    Curve3D(NurbsCurve),
    /// Caller knows the edge corresponds to a canonical surface
    /// boundary (e.g. v=0 on a Gregory cap). The factory builds
    /// from the boundary, not from intersection.
    SurfaceBoundary {
        /// Which surface owns the canonical boundary curve.
        side: SurfaceSide,
        /// Which of the four (or N) boundaries.
        boundary: BoundarySide,
    },
    /// Caller knows the desired branch by index in the canonical
    /// enumeration (preview / UI workflows).
    BranchIndex(u32),
}

/// Identifies one of the two surfaces in a pair.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum SurfaceSide {
    /// First surface argument to `compute_edge_geometry`.
    A,
    /// Second surface argument.
    B,
}

/// One of the four boundaries of a rectangular UV patch
/// (or one of N for an N-sided Gregory).
#[derive(Copy, Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum BoundarySide {
    /// `v = v_min`.
    VMin,
    /// `v = v_max`.
    VMax,
    /// `u = u_min`.
    UMin,
    /// `u = u_max`.
    UMax,
    /// For N-sided patches, the i-th boundary in the patch's own
    /// enumeration (Gregory triangular, polar K-P, …).
    Indexed(u8),
}
