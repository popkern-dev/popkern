//! A face moved with the rounds that blend it into its neighbours
//! (PK-1221).
//!
//! `offset_faces` moves a flat face whose neighbours stand square to it.
//! A round between the face and a neighbour is tangent to the face, so
//! the move declined; the edit the user means carries the round: it
//! stays a round of its radius tangent to the moved face and to the
//! neighbour, moved with the face, not stretched (a block's top lifted 6
//! takes its back round from y = 27, z = 7 to z = 13). That is the round
//! deleted (the neighbour and the face run on until they meet), the face
//! moved (its neighbours now square), and the new edge between the
//! moved face and the neighbour rounded again at the round's radius.
//! Faces are found again on the healed and moved bodies by their planes,
//! which the edit keeps (the moved face's shifted by the distance).

use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle};

use super::regions::boundary;
use super::{edge_points, failed};
use crate::error::ReplayError;

/// A round tangent to the moved face: its face, radius, and the plane of
/// the neighbour it blends the face into.
struct Round {
    face: FaceHandle,
    radius: f64,
    partner: (Point3, Vec3),
}

/// The plane of `face`, when it is flat.
fn plane_of(brep: &BRep, face: FaceHandle) -> Option<(Point3, Vec3)> {
    match brep
        .faces
        .get(face)?
        .surface
        .as_ref()
        .map(|s| &s.geometry)?
    {
        SurfaceGeometry::Plane { origin, normal, .. } => Some((*origin, normal.normalized()?)),
        _ => None,
    }
}

/// The unit normal of `face` nearest `p`.
fn normal_near(brep: &BRep, face: FaceHandle, p: &Point3) -> Option<Vec3> {
    let s = &brep.faces.get(face)?.surface.as_ref()?.geometry;
    let (u, v) = s.closest_point(p)?;
    s.normal_at(u, v).normalized()
}

/// The rounds tangent to `face` along its edges, each with the flat
/// neighbour it blends into; `None` when a neighbour is neither square
/// to the face nor such a round.
fn rounds_at(brep: &BRep, face: FaceHandle, tolerance: Tolerance) -> Option<Vec<Round>> {
    let ang = tolerance.angular();
    let mut rounds = Vec::new();
    for (edge, across) in boundary(brep, face) {
        let across = across?;
        let [_, mid, _] = edge_points(brep, edge)?;
        let (n_face, n_across) = (
            normal_near(brep, face, &mid)?,
            normal_near(brep, across, &mid)?,
        );
        if n_face.dot(&n_across).abs() <= ang {
            continue;
        }
        let Some(SurfaceGeometry::Cylinder { radius, .. }) = brep
            .faces
            .get(across)
            .and_then(|f| f.surface.as_ref())
            .map(|s| &s.geometry)
        else {
            return None;
        };
        if (1.0 - n_face.dot(&n_across)).abs() > ang {
            return None;
        }
        // The round's other support: its neighbour it is tangent to,
        // other than the face.
        let partner = boundary(brep, across).into_iter().find_map(|(e, f)| {
            let f = f.filter(|&f| f != face)?;
            let [_, m, _] = edge_points(brep, e)?;
            let tangent =
                (1.0 - normal_near(brep, f, &m)?.dot(&normal_near(brep, across, &m)?)).abs();
            (tangent <= ang).then(|| plane_of(brep, f)).flatten()
        })?;
        rounds.push(Round {
            face: across,
            radius: *radius,
            partner,
        });
    }
    Some(rounds)
}

/// The face of `brep` on the plane through `origin` with `normal`.
fn face_on(brep: &BRep, (origin, normal): (Point3, Vec3), lin: f64) -> Option<FaceHandle> {
    brep.faces.keys().find(|&f| {
        plane_of(brep, f).is_some_and(|(o, n)| {
            n.dot(&normal) > 1.0 - 1e-9 && (o - origin).dot(&normal).abs() <= lin
        })
    })
}

/// The edge `a` and `b` share.
fn shared_edge(brep: &BRep, a: FaceHandle, b: FaceHandle) -> Option<EdgeHandle> {
    boundary(brep, a)
        .into_iter()
        .find_map(|(e, f)| (f == Some(b)).then_some(e))
}

/// `face` of `brep` moved `distance` along its normal with the rounds
/// that blend it carried; `None` when it has no such rounds (the plain
/// move decides).
pub(super) fn offset_carrying_rounds(
    brep: &BRep,
    face: FaceHandle,
    distance: f64,
    tolerance: Tolerance,
) -> Option<Result<BRep, ReplayError>> {
    let rounds = rounds_at(brep, face, tolerance)?;
    if rounds.is_empty() {
        return None;
    }
    let lin = tolerance.linear();
    let moved_plane = plane_of(brep, face)?;
    Some((|| {
        let round_faces: Vec<FaceHandle> = rounds.iter().map(|r| r.face).collect();
        let healed = super::delete_faces(brep, &round_faces, tolerance)?;
        let on_healed = face_on(&healed, moved_plane, lin)
            .ok_or_else(|| failed("the moved face is gone once its rounds are deleted"))?;
        let moved = super::offset_faces(&healed, &[on_healed], distance, tolerance)?;
        let shifted = (moved_plane.0 + moved_plane.1 * distance, moved_plane.1);
        let top = face_on(&moved, shifted, lin)
            .ok_or_else(|| failed("the moved face cannot be found after the move"))?;
        let mut out = moved;
        for round in &rounds {
            let partner = face_on(&out, round.partner, lin)
                .ok_or_else(|| failed("a round's other support cannot be found after the move"))?;
            let edge = shared_edge(&out, top, partner)
                .ok_or_else(|| failed("the moved face no longer meets a round's support"))?;
            let spec = popkern_blend::FilletSpec {
                edge,
                radius: popkern_blend::FilletRadius::Constant(round.radius),
            };
            out = popkern_blend::fillet_edges(
                &out,
                &[spec],
                &popkern_blend::FilletOptions::default(),
                tolerance,
                None,
            )
            .map_err(|e| failed(format!("rounding the moved edge again: {e}")))?
            .brep;
        }
        Ok(out)
    })())
}
