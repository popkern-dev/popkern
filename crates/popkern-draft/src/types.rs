//! Public types for draft/taper operations.

use crate::error::DraftWarning;
use popkern_core::geometry::{NurbsCurve, Point3, Vec3};
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle};

/// Specification for drafting a single face.
#[derive(Debug, Clone)]
pub struct DraftSpec {
    /// The face to draft.
    pub face: FaceHandle,
    /// Draft angle in radians (positive tilts away from pull direction).
    pub angle: f64,
    /// Pull direction (direction of mold withdrawal).
    pub pull_direction: Vec3,
    /// How the neutral line is determined.
    pub neutral: NeutralReference,
    /// How the angle varies along the neutral line.
    pub angle_variation: AngleVariation,
}

/// How the neutral (stationary) line on the face is determined.
#[derive(Debug, Clone)]
pub enum NeutralReference {
    /// Neutral line is the intersection of the face with a plane.
    Plane {
        /// A point on the neutral plane.
        origin: Point3,
        /// Normal of the neutral plane.
        normal: Vec3,
    },
    /// Neutral line is an existing edge of the face.
    Edge(EdgeHandle),
    /// Automatically compute the neutral line perpendicular to the pull direction.
    Auto,
}

/// How the draft angle varies along the neutral line.
#[derive(Debug, Clone)]
pub enum AngleVariation {
    /// Constant angle along the entire neutral line.
    Constant,
    /// Variable angle: list of (parameter, angle) pairs along the neutral line.
    Variable(Vec<(f64, f64)>),
}

/// How the parting line (boundary between drafted and undrafted regions) is defined.
#[derive(Debug, Clone)]
pub enum PartingLine {
    /// Parting line is the intersection with a plane.
    Plane {
        /// A point on the parting plane.
        origin: Point3,
        /// Normal of the parting plane.
        normal: Vec3,
    },
    /// Parting line is an explicit 3D curve.
    Curve(NurbsCurve),
    /// Parting line follows existing edges.
    Edges(Vec<EdgeHandle>),
}

/// Options controlling draft operation behavior.
#[derive(Debug, Clone)]
pub struct DraftOptions {
    /// How adjacent faces are handled at draft boundaries.
    pub adjacent_mode: AdjacentMode,
    /// Maximum allowed draft angle in radians.
    pub max_angle: f64,
    /// Number of sample points for NURBS surface skinning.
    pub nurbs_samples: usize,
}

impl Default for DraftOptions {
    fn default() -> Self {
        Self {
            adjacent_mode: AdjacentMode::ExtendTrimWithFallback,
            max_angle: std::f64::consts::FRAC_PI_4,
            nurbs_samples: 20,
        }
    }
}

/// How adjacent (non-drafted) faces are handled at draft boundaries.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AdjacentMode {
    /// Extend and trim adjacent faces to meet the drafted surface.
    ExtendTrim,
    /// Insert transition faces between drafted and adjacent faces.
    TransitionFaces,
    /// Try extend-trim first; fall back to transition faces on failure.
    ExtendTrimWithFallback,
}

/// Result of a successful draft operation.
#[derive(Debug, Clone)]
pub struct DraftResult {
    /// The resulting BRep solid with drafted faces.
    pub brep: BRep,
    /// Any warnings generated during the operation.
    pub warnings: Vec<DraftWarning>,
}

#[cfg(test)]
mod tests {
    use super::*;
    use popkern_core::tolerance::Tolerance;

    #[test]
    fn draft_options_default() {
        let opts = DraftOptions::default();
        assert!((opts.max_angle - std::f64::consts::FRAC_PI_4).abs() < 1e-12);
        assert_eq!(opts.nurbs_samples, 20);
        assert_eq!(opts.adjacent_mode, AdjacentMode::ExtendTrimWithFallback);
    }

    #[test]
    fn draft_spec_construction() {
        let spec = DraftSpec {
            face: Default::default(),
            angle: 0.1,
            pull_direction: Vec3::new(0.0, 0.0, 1.0),
            neutral: NeutralReference::Auto,
            angle_variation: AngleVariation::Constant,
        };
        assert!((spec.angle - 0.1).abs() < 1e-12);
    }

    #[test]
    fn neutral_reference_variants() {
        let _plane = NeutralReference::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
        };
        let _edge = NeutralReference::Edge(Default::default());
        let _auto = NeutralReference::Auto;
    }

    #[test]
    fn angle_variation_variants() {
        let _constant = AngleVariation::Constant;
        let _variable = AngleVariation::Variable(vec![(0.0, 0.05), (1.0, 0.1)]);
    }

    #[test]
    fn parting_line_variants() {
        let _plane = PartingLine::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
        };
        let _edges = PartingLine::Edges(vec![]);
    }

    #[test]
    fn adjacent_mode_variants() {
        assert_ne!(AdjacentMode::ExtendTrim, AdjacentMode::TransitionFaces);
        let _ = format!("{:?}", AdjacentMode::ExtendTrimWithFallback);
    }

    #[test]
    fn draft_result_construction() {
        let result = DraftResult {
            brep: BRep::new(Tolerance::default()),
            warnings: vec![],
        };
        let _ = format!("{:?}", result);
    }
}
