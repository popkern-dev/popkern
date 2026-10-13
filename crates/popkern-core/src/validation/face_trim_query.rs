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

/// A point's place on a face: its chart parameters, its distance from
/// the surface there, and its trim containment.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FacePointContainment {
    /// Chart parameters of the point's closest preimage.
    pub uv: (f64, f64),
    /// Distance from the point to the surface at `uv`.
    pub distance: f64,
    /// Trim containment of `uv`.
    pub containment: FaceUvContainment,
}

impl BRep {
    /// Classify a 3D point against `face`: project it onto the face's
    /// surface with the multi-sheet projection (a point ON a multi-turn
    /// helical face must not land on a neighbouring turn — the
    /// construction projector's tuned basins did for 1551..2806 display
    /// triangle centroids per threaded body, all read Outside, PK-821)
    /// and classify the preimage. The caller judges `distance` against
    /// its own tolerance (a display facet's centroid sits its sagitta
    /// off the surface). `None` when the face has no surface or the
    /// projection or trim data is unusable.
    pub fn face_point_trim_containment(
        &self,
        face: FaceHandle,
        point: &crate::geometry::Point3,
    ) -> Option<FacePointContainment> {
        use crate::geometry::surface::SurfaceEval;
        let surface = &self.faces.get(face)?.surface.as_ref()?.geometry;
        let (u, v) =
            crate::geometry::uv_projection::project_point_to_uv_multisheet(surface, point)?;
        Some(FacePointContainment {
            uv: (u, v),
            distance: surface.point_at(u, v).distance_to(point),
            containment: self.face_uv_trim_containment(face, surface, u, v)?,
        })
    }
}

impl BRep {
    /// The face's trim window in its chart: the bounding box of its
    /// loops' pcurve samples, `((u_min, u_max), (v_min, v_max))`
    /// (vertex-projected when no loop carries a pcurve). Sampled, so a
    /// caller that must contain the face pads it.
    pub fn face_trim_uv_window(&self, face: FaceHandle) -> Option<((f64, f64), (f64, f64))> {
        let surface = &self.faces.get(face)?.surface.as_ref()?.geometry;
        Some(super::face_uv_window::face_trim_uv_window(
            self, face, surface,
        ))
    }
}
