//! A solid prepared for repeated ray classification (PK-814): the face
//! windows and hole boxes the crossing count reads, computed once per
//! face and kept across a boolean's classifications.

use popkern_core::geometry::Point3;
use popkern_core::geometry::intersection::surface_subdivision::Aabb;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};

use crate::error::{BoolError, PointClassification};

/// A face's UV window for ray crossings.
pub(crate) type UvWindow = ((f64, f64), (f64, f64));

/// A solid prepared for repeated ray classification: each face's
/// parameter window, computed once on first need and kept across every
/// classification against it (PK-814: a cut of a 64-slot disc classified
/// ~23 points against it and recomputed the windows of the faces each
/// ray met, ~0.2 ms a cut).
pub(crate) struct RaySolid<'a> {
    pub(crate) solid: &'a BRep,
    pub(crate) aabbs: &'a [(FaceHandle, Aabb)],
    windows: Vec<std::cell::OnceCell<UvWindow>>,
    /// Each face's holes with their UV boxes: a hit on a disc's top face
    /// tested its 63 holes, each rebuilding its box (~16 µs a hit).
    holes: Vec<std::cell::OnceCell<HoleBoxes>>,
}

type HoleBoxes = Vec<(
    popkern_core::topology::LoopHandle,
    Option<crate::pcurve_box::UvBox>,
)>;

impl<'a> RaySolid<'a> {
    pub(crate) fn new(solid: &'a BRep, aabbs: &'a [(FaceHandle, Aabb)]) -> Self {
        RaySolid {
            solid,
            aabbs,
            windows: (0..aabbs.len())
                .map(|_| std::cell::OnceCell::new())
                .collect(),
            holes: (0..aabbs.len())
                .map(|_| std::cell::OnceCell::new())
                .collect(),
        }
    }

    /// Whether a hit at `(u, v)` on `aabbs[k]`'s face lies in one of its
    /// holes (`hit_inside_inner_loop`, each hole's box tested first).
    pub(crate) fn hits_hole(&self, k: usize, u: f64, v: f64) -> bool {
        crate::classify_uv_loops::hit_inside_inner_loop_boxed(self.solid, self.hole_boxes(k), u, v)
    }

    /// The holes of `aabbs[k]`'s face with their UV boxes.
    fn hole_boxes(&self, k: usize) -> &HoleBoxes {
        self.holes[k].get_or_init(|| {
            let face = self.aabbs[k].0;
            self.solid.faces.get(face).map_or_else(Vec::new, |f| {
                f.inner_loops
                    .iter()
                    .map(|&lh| (lh, crate::classify_uv_loops::loop_uv_box(self.solid, lh)))
                    .collect()
            })
        })
    }

    /// The window of `aabbs[k]`'s face.
    pub(crate) fn window(&self, k: usize) -> UvWindow {
        *self.windows[k].get_or_init(|| super::face_parameter_range(self.solid, self.aabbs[k].0))
    }
}

/// [`super::point_in_solid_final`] against a prepared solid.
pub(crate) fn point_in_solid_final_on(
    point: &Point3,
    rays: &RaySolid<'_>,
    tolerance: &Tolerance,
) -> Result<PointClassification, BoolError> {
    super::point_in_solid_impl(point, rays.solid, rays.aabbs, tolerance, true, Some(rays))
        .map(|(c, _)| c)
}
