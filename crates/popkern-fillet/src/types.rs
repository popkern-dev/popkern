//! API types for fillet and chamfer specifications.

use popkern_core::topology::{EdgeHandle, FaceHandle};

use crate::error::BlendWarning;

/// Result of a successful blend operation.
#[derive(Debug, Clone)]
pub struct BlendResult {
    /// The modified BRep with blends applied.
    pub brep: popkern_core::topology::BRep,
    /// Per-edge warnings emitted during blending.
    pub edge_warnings: Vec<(EdgeHandle, BlendWarning)>,
}

/// Specification for a single fillet edge.
#[derive(Debug, Clone)]
pub struct FilletSpec {
    /// The edge to fillet.
    pub edge: EdgeHandle,
    /// The fillet radius profile.
    pub radius: FilletRadius,
}

/// Fillet radius profile along an edge.
#[derive(Debug, Clone)]
pub enum FilletRadius {
    /// Constant radius along the entire edge.
    Constant(f64),
    /// Linearly varying radius from start to end of the edge.
    Variable {
        /// Radius at the start vertex.
        start: f64,
        /// Radius at the end vertex.
        end: f64,
    },
    /// Radius defined at specific parameter positions along the edge.
    ByPoints(Vec<(f64, f64)>),
}

/// Options controlling fillet behaviour.
#[derive(Debug, Clone)]
pub struct FilletOptions {
    /// Whether to automatically propagate fillets along tangent-continuous edge chains.
    pub propagate_tangent: bool,
    /// How to handle vertex blends where three or more fillets meet.
    pub vertex_mode: VertexMode,
}

impl Default for FilletOptions {
    fn default() -> Self {
        Self {
            propagate_tangent: true,
            vertex_mode: VertexMode::Smooth,
        }
    }
}

/// Vertex blend strategy at multi-fillet junctions.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum VertexMode {
    /// Smooth G1-continuous vertex blend.
    Smooth,
    /// Mitred (flat) vertex blend.
    Mitred,
}

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
        /// Distance on the second (opposite) face.
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fillet_options_default() {
        let opts = FilletOptions::default();
        assert!(opts.propagate_tangent);
        assert_eq!(opts.vertex_mode, VertexMode::Smooth);
    }

    #[test]
    fn chamfer_options_default() {
        let opts = ChamferOptions::default();
        assert!(opts.propagate_tangent);
    }

    #[test]
    fn fillet_radius_constant() {
        let r = FilletRadius::Constant(2.5);
        match r {
            FilletRadius::Constant(v) => assert!((v - 2.5).abs() < f64::EPSILON),
            _ => panic!("expected Constant"),
        }
    }
}
