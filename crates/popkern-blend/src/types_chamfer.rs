use popkern_core::topology::{EdgeHandle, FaceHandle};

/// Specification for a single chamfer edge.
#[derive(Debug, Clone)]
pub struct ChamferSpec {
    /// The edge to chamfer.
    pub edge: EdgeHandle,
    /// The chamfer geometry.
    pub kind: ChamferKind,
}

/// Chamfer geometry definition.
#[derive(Debug, Clone)]
pub enum ChamferKind {
    /// Equal-distance chamfer (45-degree bevel).
    Equal(f64),
    /// Two-distance chamfer with specified face orientation.
    TwoDistance {
        /// Reference face for the first distance.
        face_1: FaceHandle,
        /// Distance on the first face.
        dist_1: f64,
        /// Distance on the second face.
        dist_2: f64,
    },
    /// Distance-and-angle chamfer.
    DistanceAngle {
        /// Reference face for the distance measurement.
        face: FaceHandle,
        /// Distance from the edge on the reference face.
        distance: f64,
        /// Angle in radians from the reference face.
        angle: f64,
    },
}

/// Options controlling chamfer behaviour.
#[derive(Debug, Clone)]
pub struct ChamferOptions {
    /// Whether to automatically propagate chamfers along tangent-continuous edge chains.
    pub propagate_tangent: bool,
}

impl Default for ChamferOptions {
    fn default() -> Self {
        Self {
            propagate_tangent: true,
        }
    }
}
