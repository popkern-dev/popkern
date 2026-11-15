//! `Feature::DeleteFaces`: take blend faces out of a solid and close the
//! gap with their neighbours run on along their own surfaces until they
//! meet (PK-1199).
//!
//! A blend face (a round, a bevel, a countersink, a corner ball) lies
//! across the corner its neighbours would make if run on. That corner
//! region is a solid of its own: the face's cross-section closed by its
//! neighbours' sections to where they meet, extruded along a straight
//! spine or revolved about a circular one, taken whole to the face's
//! chord (what lies across the face from the corner is body where the
//! region is added and air where it is taken away); for a corner ball,
//! the cube its three tangent planes span, less the ball. Corner balls
//! close first, so the edges' regions meet closed corners. Deleting the face adds
//! the region where it lies outside the body (a convex blend) and takes
//! it away where it lies inside (a concave one). Each region is built
//! from the body as given, so several faces delete together. A face of
//! any other form, or whose neighbours do not meet in a line of the
//! section (a plane along a straight spine; a plane square to, or a
//! cylinder or cone about, a circular one), declines with the reason.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::transform::Transform3;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};
use popkern_sketch::profile::{AnalyticProfile, PathSegment, ProfileShape};

use crate::error::ReplayError;

#[path = "replay_delete_faces_regions.rs"]
mod regions;

/// `brep` with `faces` deleted and their gaps closed.
pub(super) fn delete_faces(
    brep: &BRep,
    faces: &[FaceHandle],
    tolerance: Tolerance,
) -> Result<BRep, ReplayError> {
    let mut fills = Vec::with_capacity(faces.len());
    let ball_first = |f: &FaceHandle| {
        !matches!(
            brep.faces
                .get(*f)
                .and_then(|x| x.surface.as_ref())
                .map(|s| &s.geometry),
            Some(SurfaceGeometry::Sphere { .. })
        )
    };
    let mut ordered = faces.to_vec();
    ordered.sort_by_key(ball_first);
    for face in ordered {
        fills.push(corner_fill(brep, face, tolerance)?);
    }
    let mut out = brep.clone();
    for Fill { region, adds } in fills {
        out = if adds {
            popkern_bool::union(out, region, tolerance, None)
        } else {
            popkern_bool::subtraction(out, region, tolerance, None)
        }
        .map_err(|e| failed(format!("closing the gap: {e}")))?;
    }
    Ok(out)
}

/// `brep` with each of `faces` moved `distance` along its outward
/// normal, its neighbours running on to meet it (PK-1199): the region
/// the face sweeps is added when the face moves out of the body and taken away when it
/// moves in. A flat face sweeps its prism (the face swept along its
/// normal), a hole's wall the ring between its radius and the new one. Neighbours run on along their
/// own surfaces only where they stand square to the face (an extruded
/// block's top, a hole's ends); a face whose neighbours lean onto it
/// would need them re-solved, and the region would serve the wrong
/// solid, so such a face declines.
pub(super) fn offset_faces(
    brep: &BRep,
    faces: &[FaceHandle],
    distance: f64,
    tolerance: Tolerance,
) -> Result<BRep, ReplayError> {
    if distance.abs() <= tolerance.linear() {
        return Ok(brep.clone());
    }
    let mut regions = Vec::with_capacity(faces.len());
    for &face in faces {
        regions::require_square_neighbours(brep, face, tolerance)?;
        regions.push(regions::face_sweep(brep, face, distance, tolerance)?);
    }
    let mut out = brep.clone();
    for region in regions {
        out = if distance > 0.0 {
            popkern_bool::union(out, region, tolerance, None)
        } else {
            popkern_bool::subtraction(out, region, tolerance, None)
        }
        .map_err(|e| failed(format!("closing the face's sweep: {e}")))?;
    }
    Ok(out)
}

/// `Feature::DeleteFaces`: the faces named taken out and their gaps
/// closed (`replay_delete_faces`).
pub(in crate::replay) fn delete_faces_feature(
    brep: &mut BRep,
    faces: &[popkern_core::naming::EntityName],
    tolerance: &Tolerance,
) -> Result<crate::replay_side_effects::ExecuteSideEffects, ReplayError> {
    let handles = crate::resolve::resolve_faces(faces, brep)?;
    *brep = delete_faces(brep, &handles, *tolerance)?;
    Ok(crate::replay_side_effects::ExecuteSideEffects::empty())
}

/// `Feature::RemoveFaces`: the faces named taken away, the body left
/// open along their edges (a sheet, thicken's input).
pub(in crate::replay) fn remove_faces_feature(
    brep: &mut BRep,
    faces: &[popkern_core::naming::EntityName],
) -> Result<crate::replay_side_effects::ExecuteSideEffects, ReplayError> {
    let handles = crate::resolve::resolve_faces(faces, brep)?;
    brep.retire_faces_to_boundary(&handles)
        .map_err(|e| ReplayError::FeatureFailed {
            label: "remove_faces".into(),
            error: crate::failure::FeatureFailure::other(e.to_string()).into(),
        })?;
    // The faces kept are the feature's result, a sheet.
    brep.settle_open_shells();
    Ok(crate::replay_side_effects::ExecuteSideEffects::empty())
}

/// `Feature::OffsetFaces`: each face named moved `distance` along its
/// outward normal, its neighbours following (`replay_delete_faces`).
pub(in crate::replay) fn offset_faces_feature(
    brep: &mut BRep,
    faces: &[popkern_core::naming::EntityName],
    distance: f64,
    tolerance: &Tolerance,
) -> Result<crate::replay_side_effects::ExecuteSideEffects, ReplayError> {
    let handles = crate::resolve::resolve_faces(faces, brep)?;
    *brep = offset_faces(brep, &handles, distance, *tolerance)?;
    Ok(crate::replay_side_effects::ExecuteSideEffects::empty())
}

/// A deleted face's corner region and whether it is added.
struct Fill {
    region: BRep,
    adds: bool,
}

fn failed(reason: impl Into<String>) -> ReplayError {
    ReplayError::FeatureFailed {
        label: "delete_faces".into(),
        error: crate::failure::FeatureFailure::other(reason.into()).into(),
    }
}

/// The corner region of `face`.
fn corner_fill(brep: &BRep, face: FaceHandle, tolerance: Tolerance) -> Result<Fill, ReplayError> {
    let surface = brep
        .faces
        .get(face)
        .and_then(|f| f.surface.as_ref())
        .map(|s| s.geometry.clone())
        .ok_or_else(|| failed("the face has no surface"))?;
    let (region, corner, face_mid) = match &surface {
        SurfaceGeometry::Cylinder {
            origin,
            axis,
            radius,
        } => regions::straight(brep, face, *axis, Some((*origin, *radius)), tolerance)?,
        SurfaceGeometry::Plane { .. } => {
            let d = regions::longest_straight_edge(brep, face)
                .ok_or_else(|| failed("a flat face deletes as a bevel along its straight edges"))?;
            regions::straight(brep, face, d, None, tolerance)?
        }
        SurfaceGeometry::Torus {
            center,
            axis,
            major_radius,
            minor_radius,
        } => regions::circular(
            brep,
            face,
            *center,
            *axis,
            Some((*major_radius, *minor_radius)),
            tolerance,
        )?,
        SurfaceGeometry::Cone { apex, axis, .. } => {
            regions::circular(brep, face, *apex, *axis, None, tolerance)?
        }
        SurfaceGeometry::Sphere { center, radius, .. } => {
            regions::corner_ball(brep, face, *center, *radius, tolerance)?
        }
        _ => return Err(failed("only rounds, bevels and corner balls delete")),
    };
    // The face's outward normal at its section's middle points toward
    // the corner when the corner lies outside the body (a convex blend:
    // its region is added); away from it for a concave one.
    let outward = {
        let sof = brep.faces.get(face).and_then(|f| f.surface.as_ref());
        let (u, v) = surface
            .closest_point(&face_mid)
            .ok_or_else(|| failed("the face's middle has no chart point"))?;
        let n = surface.normal_at(u, v);
        if sof.is_some_and(|s| s.orientation) {
            n
        } else {
            -n
        }
    };
    let adds = outward.dot(&(corner - face_mid)) > 0.0;
    Ok(Fill { region, adds })
}

/// The profile `start` → `segments` (closed) in a local frame, extruded
/// along local z from `z0` by `h`, placed by `frame`.
fn extrude_section(
    start: [f64; 2],
    segments: Vec<PathSegment<f64>>,
    (z0, h): (f64, f64),
    frame: &Transform3,
    tolerance: Tolerance,
) -> Result<BRep, ReplayError> {
    let profile = AnalyticProfile::new(ProfileShape::Path { start, segments });
    let local = super::extruded_profile("delete_faces", &profile, z0, h, tolerance)?;
    local.transform(frame).map_err(|e| failed(e.to_string()))
}

/// The profile in (ρ, h) on the local xz plane revolved a whole turn
/// about local z, placed by `frame`.
fn revolve_section(
    start: [f64; 2],
    segments: Vec<PathSegment<f64>>,
    frame: &Transform3,
    tolerance: Tolerance,
) -> Result<BRep, ReplayError> {
    let profile = AnalyticProfile::new(ProfileShape::Path { start, segments });
    let local = super::revolved_profile("delete_faces", &profile, tolerance)?;
    local.transform(frame).map_err(|e| failed(e.to_string()))
}

/// Points along an edge's curve: its two ends and its middle.
fn edge_points(brep: &BRep, edge: popkern_core::topology::EdgeHandle) -> Option<[Point3; 3]> {
    let c = brep.edges.get(edge)?.curve.as_ref()?;
    let (a, b) = (c.t_start, c.t_end);
    Some([
        c.geometry.point_at(a),
        c.geometry.point_at(0.5 * (a + b)),
        c.geometry.point_at(b),
    ])
}

/// The unit vector `v` and one square to it.
fn frame_of(origin: Point3, z: Vec3) -> Option<(Vec3, Vec3, Transform3)> {
    let z = z.normalized()?;
    let seed = if z.x.abs() < 0.9 {
        Vec3::new(1.0, 0.0, 0.0)
    } else {
        Vec3::new(0.0, 1.0, 0.0)
    };
    let x = (seed - z * seed.dot(&z)).normalized()?;
    let t = Transform3::from_axis_placement(origin, z, x).ok()?;
    Some((z, x, t))
}
