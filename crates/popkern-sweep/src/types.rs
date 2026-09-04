//! API types for sweep and loft operations.

use popkern_core::geometry::{NurbsCurve, Point3, Vec3};
use popkern_core::topology::{BRep, FaceHandle};

use crate::error::SweepWarning;

/// A profile to sweep along a path.
#[derive(Debug, Clone)]
pub enum SweepProfile {
    /// A single wire (open or closed NURBS curve) in the PATH FRAME's
    /// 2-D coordinates: `x` maps to the frame normal, `y` to the
    /// binormal, and `z` must be constant (it is dropped after the
    /// wire is recentred). A world-space profile on a plane other than
    /// XY has to be expressed in that plane's frame first and the
    /// swept body placed afterwards; the sweep rejects a wire with a
    /// `z` extent instead of silently collapsing its spans (PK-416).
    Wire(NurbsCurve),
    /// A face from an existing BRep.
    Face {
        /// The BRep containing the face.
        brep: Box<BRep>,
        /// Handle to the face to sweep.
        face: FaceHandle,
    },
    /// Multiple wires forming a compound profile.
    MultiWire(Vec<NurbsCurve>),
}

/// A cross-section for lofting.
#[derive(Debug, Clone)]
pub struct LoftSection {
    /// The section curve.
    pub curve: NurbsCurve,
    /// Optional explicit position (overrides curve start point).
    pub position: Option<Point3>,
}

/// Pipe radius specification.
#[derive(Debug, Clone)]
pub enum PipeRadius {
    /// Constant radius along the path.
    Constant(f64),
    /// Variable radius specified as (arc_length_fraction, radius) pairs.
    Variable(Vec<(f64, f64)>),
}

/// Helix specification.
#[derive(Debug, Clone)]
pub struct HelixSpec {
    /// Origin point of the helix axis.
    pub axis_origin: Point3,
    /// Direction of the helix axis (will be normalized).
    pub axis_dir: Vec3,
    /// Radius at the start of the helix.
    pub start_radius: f64,
    /// Radius at the end of the helix.
    pub end_radius: f64,
    /// Pitch specification.
    pub pitch: HelixPitch,
    /// Number of turns.
    pub turns: f64,
    /// Phase: angle in radians about the axis at which the helix
    /// starts, measured from the axis frame's reference direction.
    /// Lets multiple helical bodies about one axis be phased apart
    /// (propeller blades, multi-start threads).
    pub start_angle: f64,
    /// Winding sense. `false` = right-hand (counter-clockwise about
    /// `axis_dir` while advancing along it); `true` = left-hand.
    pub left_handed: bool,
}

/// Helix pitch specification.
#[derive(Debug, Clone)]
pub enum HelixPitch {
    /// Constant pitch (axial distance per turn).
    Constant(f64),
    /// Variable pitch specified as (turn_fraction, pitch) pairs.
    Variable(Vec<(f64, f64)>),
}

/// Options controlling sweep behavior.
#[derive(Debug, Clone)]
pub struct SweepOptions {
    /// Frame computation method.
    pub frame: FrameMethod,
    /// Total twist angle (radians) applied along the path.
    pub twist: f64,
    /// Scaling method along the path.
    pub scale: ScaleMethod,
    /// Whether to allow self-intersecting results.
    pub allow_self_intersection: bool,
}

impl Default for SweepOptions {
    fn default() -> Self {
        Self {
            frame: FrameMethod::Bishop,
            twist: 0.0,
            scale: ScaleMethod::None,
            allow_self_intersection: false,
        }
    }
}

/// Method for computing orientation frames along a path.
#[derive(Debug, Clone)]
pub enum FrameMethod {
    /// Rotation-minimizing (Bishop) frame.
    Bishop,
    /// Frenet-Serret frame (T, N, B from curvature).
    Frenet,
    /// Fixed up-direction projected onto the normal plane.
    Fixed(Vec3),
    /// Auxiliary curve defining the up-direction at each point.
    Auxiliary(NurbsCurve),
}

/// Scaling method along a sweep path.
#[derive(Debug, Clone)]
pub enum ScaleMethod {
    /// No scaling.
    None,
    /// Uniform scaling specified as (arc_length_fraction, scale_factor) pairs.
    Uniform(Vec<(f64, f64)>),
}

/// Options controlling loft behavior.
#[derive(Debug, Clone)]
pub struct LoftOptions {
    /// Continuity constraint at the start section.
    pub start_continuity: Continuity,
    /// Continuity constraint at the end section.
    pub end_continuity: Continuity,
    /// Whether the loft is closed (last section connects to first).
    pub closed: bool,
}

impl Default for LoftOptions {
    fn default() -> Self {
        Self {
            start_continuity: Continuity::Free,
            end_continuity: Continuity::Free,
            closed: false,
        }
    }
}

/// Continuity constraint for loft ends.
#[derive(Debug, Clone)]
pub enum Continuity {
    /// No constraint (natural end).
    Free,
    /// Tangent continuity with an adjacent face.
    G1 {
        /// The BRep containing the adjacent face.
        face: BRep,
        /// Handle to the adjacent face.
        face_handle: FaceHandle,
    },
    /// Curvature continuity with an adjacent face.
    G2 {
        /// The BRep containing the adjacent face.
        face: BRep,
        /// Handle to the adjacent face.
        face_handle: FaceHandle,
    },
}

/// Result of a sweep or loft operation.
#[derive(Debug)]
pub struct SweepResult {
    /// The resulting BRep solid.
    pub brep: BRep,
    /// Non-fatal warnings emitted during the operation.
    pub warnings: Vec<SweepWarning>,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sweep_options_default() {
        let opts = SweepOptions::default();
        assert!(matches!(opts.frame, FrameMethod::Bishop));
        assert!((opts.twist - 0.0).abs() < f64::EPSILON);
        assert!(matches!(opts.scale, ScaleMethod::None));
        assert!(!opts.allow_self_intersection);
    }

    #[test]
    fn loft_options_default() {
        let opts = LoftOptions::default();
        assert!(matches!(opts.start_continuity, Continuity::Free));
        assert!(matches!(opts.end_continuity, Continuity::Free));
        assert!(!opts.closed);
    }
}
