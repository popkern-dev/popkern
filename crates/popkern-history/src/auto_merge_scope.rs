//! The body a cut subtracts from when nothing names one. A cut's host
//! is normally the body owning the face its sketch stands on, or the
//! face it ends on (`auto_merge`); a cut sketched on a free plane with
//! a blind or through end names neither, and used to be kept as a
//! body of its own, a positive solid where material should have gone.
//!
//! Its scope is the bodies its tool can reach: those whose bounding
//! box overlaps the tool's (a necessary condition for removing any
//! material). One such body is the host; the cut's own check
//! (`CutRemovedNothing`) still judges whether material went. None
//! fails with `CutRemovedNothing`; several fail with
//! `CutTargetAmbiguous` instead of choosing among them.

use popkern_core::geometry::Point3;
use popkern_core::geometry::intersection::surface_subdivision::Aabb;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use crate::error::ReplayError;
use crate::failure::FeatureFailure;
use crate::node::BodyId;
use crate::working_bodies::WorkingBodies;

/// A body's box: the union of its faces' prefilter boxes, which
/// contain their faces (vertices alone miss a sphere's bulk).
fn body_box(brep: &BRep) -> Option<Aabb> {
    let corners: Vec<Point3> = popkern_bool::face_bounding_boxes(brep)
        .into_iter()
        .flat_map(|(_, b)| [b.min, b.max])
        .collect();
    Aabb::from_points(&corners)
}

/// The host of a cut whose tool is `tool_body`, chosen by reach.
pub(crate) fn cut_scope_host(
    tool_body: BodyId,
    body_breps: &WorkingBodies,
    tolerance: &Tolerance,
) -> Result<BodyId, ReplayError> {
    let fail = |error: FeatureFailure| ReplayError::FeatureFailed {
        label: "auto_merge".into(),
        error: Box::new(error),
    };
    let tool = body_breps
        .get(&tool_body)
        .ok_or_else(|| ReplayError::FeatureFailed {
            label: "auto_merge".into(),
            error: format!("tool body {tool_body:?} has no BRep").into(),
        })?;
    let tool_faces = tool.faces.len();
    let Some(reach) = body_box(tool) else {
        return Err(fail(FeatureFailure::CutRemovedNothing { tool_faces }));
    };
    let mut bodies: Vec<&BodyId> = body_breps.keys().collect();
    bodies.sort_by_key(|b| b.0);
    let candidates: Vec<BodyId> = bodies
        .into_iter()
        .copied()
        .filter(|b| *b != tool_body)
        .filter(|b| {
            body_breps
                .get(b)
                .and_then(body_box)
                .is_some_and(|other| other.overlaps(&reach, tolerance.linear()))
        })
        .collect();
    match candidates.as_slice() {
        [host] => Ok(*host),
        [] => Err(fail(FeatureFailure::CutRemovedNothing { tool_faces })),
        several => Err(fail(FeatureFailure::CutTargetAmbiguous {
            candidates: several.len(),
        })),
    }
}
