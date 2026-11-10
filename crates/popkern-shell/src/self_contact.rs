//! A shell whose skin touches itself: zero-thickness self-contact
//! (the tangency foundations' shell consumer, design §8.4, F3).
//!
//! At the thickness where a thin wall closes exactly, the offset skins
//! of the wall's two sides meet without crossing: a box with a through
//! hole whose lowest line runs 2t above the floor, shelled by t, puts
//! the hole's offset cylinder tangent to the floor's offset plane along
//! a line. Nothing crosses, so neither the self-intersection trim nor
//! `validate_all` sees it, and the shell served the touching skins as
//! one manifold solid (the material meeting itself along the line with
//! no shared edge, the cavity split in two). Under the kernel's default
//! policy a result that touches itself is refused with its contact, as
//! a boolean's is (`NonManifoldResult`).
//!
//! Detection reads the served body: every pair of faces that share no
//! vertex and whose boxes come within the tolerance is searched for
//! tangent contacts (`find_surface_contacts`, the critical points of
//! the signed gap) over their trims' parameter boxes, and a contact
//! counts when it lies within the tolerance of both trimmed patches.
//! Two parallel planes within the tolerance of each other (a closing
//! wall between flat sides) are a contact where a boundary sample of
//! one lies on the other's patch.

use std::collections::{HashMap, HashSet};

use popkern_core::geometry::Point3;
use popkern_core::geometry::intersection::contact_classify::ContactScale;
use popkern_core::geometry::intersection::contact_discovery::find_surface_contacts;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle, VertexHandle};

use crate::body_distance::FaceTrimEntry;

/// Where two faces of a body touch without sharing topology.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct SelfContact {
    /// A point of the contact.
    pub location: Point3,
    /// The two faces that touch there.
    pub faces: (FaceHandle, FaceHandle),
}

/// Shared implementation for [`crate::shell`] and
/// [`crate::shell_variable`]: the shell as built
/// ([`crate::shell_impl_unchecked`]), declined when it touches itself
/// ([`crate::ShellError::SelfContact`]).
#[allow(clippy::too_many_arguments)]
pub(crate) fn shell_impl(
    brep: &BRep,
    thickness: f64,
    direction: crate::ShellDirection,
    removed_faces: &[FaceHandle],
    face_thicknesses: &[(FaceHandle, f64)],
    options: &crate::ShellOptions,
    tolerance: Tolerance,
    tracker: &mut popkern_core::progress::StageTracker<'_, '_>,
    naming_sink: Option<&mut crate::ShellProvenance>,
) -> Result<crate::ShellResult, crate::ShellError> {
    let result = crate::shell_impl_unchecked(
        brep,
        thickness,
        direction,
        removed_faces,
        face_thicknesses,
        options,
        tolerance,
        tracker,
        naming_sink,
    )?;
    let Some(c) = first_self_contact(&result.brep, &tolerance) else {
        return Ok(result);
    };
    let declined = crate::ShellError::SelfContact {
        location: c.location,
        faces: c.faces,
    };
    if options.non_manifold != popkern_bool::NonManifoldPolicy::General {
        return Err(declined);
    }
    // One general body: the contact imprinted on both faces and glued
    // (`BRep::glue_touching`, a piece touching itself), served only
    // when the glued body validates and no longer touches itself.
    let mut result = result;
    if result.brep.glue_touching().is_err()
        || !result.brep.validate_all().is_empty()
        || first_self_contact(&result.brep, &tolerance).is_some()
    {
        return Err(declined);
    }
    Ok(result)
}

/// The first place, in face handle order, where `brep` touches itself.
pub(crate) fn first_self_contact(brep: &BRep, tolerance: &Tolerance) -> Option<SelfContact> {
    let tol = tolerance.linear();
    let boxes: HashMap<FaceHandle, _> = popkern_bool::face_bounding_boxes(brep)
        .into_iter()
        .collect();
    let mut verts: HashMap<FaceHandle, HashSet<VertexHandle>> = HashMap::new();
    for (_, he) in &brep.half_edges {
        if let Some(l) = brep.loops.get(he.loop_) {
            verts.entry(l.face).or_default().insert(he.origin);
        }
    }
    let mut faces: Vec<FaceHandle> = brep.faces.keys().collect();
    faces.sort();
    let mut entries: HashMap<FaceHandle, Option<FaceTrimEntry>> = HashMap::new();
    for (i, &fa) in faces.iter().enumerate() {
        for &fb in &faces[i + 1..] {
            let (Some(ba), Some(bb)) = (boxes.get(&fa), boxes.get(&fb)) else {
                continue;
            };
            let apart = (0..3).any(|k| {
                let (lo_a, hi_a) = (coord(&ba.min, k), coord(&ba.max, k));
                let (lo_b, hi_b) = (coord(&bb.min, k), coord(&bb.max, k));
                lo_b - hi_a > tol || lo_a - hi_b > tol
            });
            let share = match (verts.get(&fa), verts.get(&fb)) {
                (Some(a), Some(b)) => !a.is_disjoint(b),
                _ => true,
            };
            if apart || share {
                continue;
            }
            for f in [fa, fb] {
                entries
                    .entry(f)
                    .or_insert_with(|| FaceTrimEntry::build(brep, f));
            }
            let (Some(Some(ea)), Some(Some(eb))) = (entries.get(&fa), entries.get(&fb)) else {
                continue;
            };
            let extent = ba.min.distance_to(&ba.max).max(bb.min.distance_to(&bb.max));
            if let Some(location) = pair_contact(ea, eb, tol, extent) {
                return Some(SelfContact {
                    location,
                    faces: (fa, fb),
                });
            }
        }
    }
    None
}

fn coord(p: &Point3, k: usize) -> f64 {
    match k {
        0 => p.x,
        1 => p.y,
        _ => p.z,
    }
}

/// A point where the trimmed patches `a` and `b` touch, if any.
fn pair_contact(a: &FaceTrimEntry, b: &FaceTrimEntry, tol: f64, extent: f64) -> Option<Point3> {
    let on_both = |p: &Point3| {
        a.distance(p).is_some_and(|d| d <= tol) && b.distance(p).is_some_and(|d| d <= tol)
    };
    if let (SurfaceGeometry::Plane { .. }, SurfaceGeometry::Plane { .. }) = (&a.surface, &b.surface)
    {
        // Planes have no isolated critical points of their gap: a
        // contact is an area, found from either boundary.
        return a
            .boundary_samples()
            .iter()
            .chain(b.boundary_samples())
            .find(|p| on_both(p))
            .copied();
    }
    let window = |e: &FaceTrimEntry| e.uv_bbox().map(|(u0, v0, u1, v1)| ((u0, u1), (v0, v1)));
    let (wa, wb) = (window(a)?, window(b)?);
    let found = find_surface_contacts(
        &a.surface,
        wa,
        &b.surface,
        wb,
        ContactScale {
            linear: tol,
            extent,
        },
    );
    found.contacts.iter().map(|c| c.point).find(|p| on_both(p))
}

#[cfg(test)]
#[path = "self_contact_tests.rs"]
mod tests;
