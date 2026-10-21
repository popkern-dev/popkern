//! The shared context of one cellular fragment classification, and the
//! point classification through its prepared operands (PK-814).

use std::collections::HashSet;

use popkern_core::geometry::Point3;
use popkern_core::geometry::intersection::surface_subdivision::Aabb;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};

use crate::adaptive_precision::{robust_point_in_solid_final, robust_point_in_solid_final_on};
use crate::error::PointClassification;
use crate::merge::graph::IntersectionGraph;

/// Shared context for one `classify_fragments_cellular` run: the
/// operands, their face AABBs, the coincident-overlap face sets, and
/// the calibrated probe offsets.
pub(in crate::cellular) struct ClassifyCtx<'a> {
    pub graph: &'a IntersectionGraph,
    pub brep_a: &'a BRep,
    pub brep_b: &'a BRep,
    pub aabbs_a: &'a [(FaceHandle, Aabb)],
    pub aabbs_b: &'a [(FaceHandle, Aabb)],
    pub overlap_a: &'a HashSet<FaceHandle>,
    pub overlap_b: &'a HashSet<FaceHandle>,
    pub probe_eps: f64,
    pub probe_eps_coarse: f64,
    pub tolerance: &'a Tolerance,
    /// Each operand prepared for ray classification (face windows kept
    /// across the classification, PK-814).
    pub rays_a: &'a crate::classify::RaySolid<'a>,
    pub rays_b: &'a crate::classify::RaySolid<'a>,
}

/// [`robust_point_in_solid_final`] through `ctx`'s prepared operand when
/// `solid` and `aabbs` are one of the two (the uncached path otherwise).
pub(super) fn classify_against(
    ctx: &ClassifyCtx<'_>,
    point: &Point3,
    solid: &BRep,
    aabbs: &[(FaceHandle, Aabb)],
    tolerance: &Tolerance,
) -> PointClassification {
    let prepared = [ctx.rays_a, ctx.rays_b]
        .into_iter()
        .find(|r| std::ptr::eq(r.solid, solid) && std::ptr::eq(r.aabbs, aabbs));
    match prepared {
        Some(rays) => robust_point_in_solid_final_on(point, rays, tolerance),
        None => robust_point_in_solid_final(point, solid, aabbs, tolerance),
    }
}
