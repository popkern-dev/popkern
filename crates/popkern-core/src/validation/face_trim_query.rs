//! Public UV trim-containment query (PK-264).
//!
//! Exposes the self-intersection validator's even-odd point-in-face
//! classification as a stable query: given a face and a UV point on
//! its surface chart, is the point inside the face's trimmed region?
//!
//! The classification is **winding-independent** (even-odd crossing
//! parity — loop direction never enters), **hole-correct** (a point
//! inside an inner loop classifies `Outside`), and **periodic-aware**
//! (u-periodic charts test every period image against the
//! branch-unwrapped loop polygon). That combination is exactly what
//! consumers need when stored loop windings cannot be trusted — the
//! popkern-blend convexity classifier's interior-side witness is the
//! first client (primitive constructors encode outwardness in the
//! orientation flag alone and leave loop windings unenforced, so any
//! interior-on-the-left rule silently mirrors on those faces).

use crate::geometry::surface::SurfaceGeometry;
use crate::topology::{BRep, FaceHandle};

use super::self_intersection_trim::{TrimContainmentForSense, UvPointForSense};

/// Where a UV point sits relative to a face's trimmed region.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FaceUvContainment {
    /// Outside the outer loop, or inside a hole.
    Outside,
    /// Within the parametric boundary band of a loop.
    Boundary,
    /// Strictly inside the trimmed region.
    Interior,
}

impl BRep {
    /// Classify `(u, v)` on `surface` against `face`'s trim loops by
    /// even-odd containment.
    ///
    /// `surface` must be the face's surface geometry (passed
    /// explicitly so callers holding a resolved [`SurfaceGeometry`]
    /// avoid a second lookup; the loops are sampled through their own
    /// pcurves, with 3D-projection fallback for pcurve-less
    /// half-edges).
    ///
    /// Returns `None` when the face or its loop data is missing or
    /// unsampleable — callers must treat that as "no answer", not
    /// "outside".
    pub fn face_uv_trim_containment(
        &self,
        face: FaceHandle,
        surface: &SurfaceGeometry,
        u: f64,
        v: f64,
    ) -> Option<FaceUvContainment> {
        let containment = self.face_trim_classify_uv(face, surface, UvPointForSense { u, v })?;
        Some(match containment {
            TrimContainmentForSense::Outside => FaceUvContainment::Outside,
            TrimContainmentForSense::Boundary => FaceUvContainment::Boundary,
            TrimContainmentForSense::Interior => FaceUvContainment::Interior,
        })
    }
}
