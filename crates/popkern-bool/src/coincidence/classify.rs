//! Orientation classification for coincident face pairs.
//!
//! Determines whether two coincident faces have the same or opposite
//! normal orientation. This is critical for correct boolean results:
//! the spec (SS3.3) prescribes different rules for each case.

use popkern_core::geometry::intersection::surface_subdivision::Aabb;
use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};

use crate::error::PointClassification;

use super::types::CoincidenceOrientation;

/// Classify the orientation relationship between two coincident faces.
///
/// Computes each face's material-aware outward normal (see
/// [`face_outward_normal`]) and compares them: normals agreeing (dot ≥ 0)
/// are `Same`-oriented, opposing (dot < 0) are `Opposite`-oriented. The
/// outward direction is resolved by which side the operand's solid
/// material lies on, so it is correct regardless of how the cap plane's
/// normal and `orientation` flag were stored at construction.
///
/// The two are compared at ONE point, the middle of face A: each face's
/// outward normal is known at its own middle, which on a curved surface
/// says nothing about the other's (two half balls of one sphere about
/// perpendicular axes have outward normals at their middles a right
/// angle apart, and read `Opposite`; PK-1147). A face's outward side is
/// one side of its surface throughout, so each outward normal is kept
/// as a sign against its surface's normal and both surfaces' normals
/// are taken at the common point.
pub(crate) fn classify_face_orientation(
    a: &BRep,
    face_a: FaceHandle,
    aabbs_a: &[(FaceHandle, Aabb)],
    b: &BRep,
    face_b: FaceHandle,
    aabbs_b: &[(FaceHandle, Aabb)],
    tolerance: &Tolerance,
) -> CoincidenceOrientation {
    let outward_a = face_outward_normal(a, face_a, aabbs_a, tolerance);
    let outward_b = face_outward_normal(b, face_b, aabbs_b, tolerance);
    let (normal_a, normal_b) = match (outward_a.at, outward_b.at) {
        (Some(pa), Some(pb)) => {
            let (Some(sa), Some(sb)) = (surface_of(a, face_a), surface_of(b, face_b)) else {
                return orientation_of(outward_a.normal, outward_b.normal);
            };
            // Each outward normal as a sign on its surface's normal.
            let sign = |s: &popkern_core::geometry::SurfaceGeometry, p: &Point3, n: Vec3| {
                s.closest_point(p)
                    .map(|(u, v)| s.normal_at(u, v).dot(&n).signum())
            };
            let (Some(ka), Some(kb)) = (
                sign(sa, &pa, outward_a.normal),
                sign(sb, &pb, outward_b.normal),
            ) else {
                return orientation_of(outward_a.normal, outward_b.normal);
            };
            let at = |s: &popkern_core::geometry::SurfaceGeometry| {
                s.closest_point(&pa).map(|(u, v)| s.normal_at(u, v))
            };
            match (at(sa), at(sb)) {
                (Some(na), Some(nb)) => (na * ka, nb * kb),
                _ => (outward_a.normal, outward_b.normal),
            }
        }
        _ => (outward_a.normal, outward_b.normal),
    };
    orientation_of(normal_a, normal_b)
}

fn surface_of(brep: &BRep, face: FaceHandle) -> Option<&popkern_core::geometry::SurfaceGeometry> {
    brep.faces.get(face)?.surface.as_ref().map(|s| &s.geometry)
}

fn orientation_of(normal_a: Vec3, normal_b: Vec3) -> CoincidenceOrientation {
    if normal_a.dot(&normal_b) >= 0.0 {
        CoincidenceOrientation::Same
    } else {
        CoincidenceOrientation::Opposite
    }
}

/// Compute the outward-facing normal of a face at its centroid.
///
/// The stored surface normal and `SurfaceOnFace::orientation` flag do
/// not reliably encode the OUTWARD (away-from-solid) direction for
/// every face: a circle-extrude stores both of a cylinder's cap planes
/// with the same `+n` normal and `orientation = true`, so the bottom
/// cap — whose true outward normal is `-n` — would read `+n` (PK-435:
/// this flipped a bottom disc/rectangle coincident pair to `Opposite`
/// and dropped the kept union skin). Instead, take the plane normal's
/// axis and orient it by MATERIAL SIDE: the outward normal is the one
/// under which a point just off the face centroid lies OUTSIDE the
/// solid. This is construction-independent and correct for any closed
/// operand. Falls back to the stored (flag-adjusted) normal only when
/// the material probe is inconclusive (e.g. an open-sheet operand).
fn face_outward_normal(
    brep: &BRep,
    face: FaceHandle,
    aabbs: &[(FaceHandle, Aabb)],
    tolerance: &Tolerance,
) -> Outward {
    let fallback = |normal| Outward { normal, at: None };
    let sof = match brep.faces[face].surface.as_ref() {
        Some(s) => s,
        None => return fallback(Vec3::new(0.0, 0.0, 1.0)),
    };

    // Boundary positions sampled along edge CURVES, not just vertices:
    // a disc's rim carries exactly ONE vertex (the seam), so a
    // vertex-only "centroid" sits ON the rim — radially on the
    // adjacent wall surface — and both material probes below read
    // boundary noise (measured: the hub's bottom cap fell back to the
    // flag-adjusted normal and its coincident pair mis-typed
    // Opposite, while the top cap's identical fallback was right by
    // luck). Curve samples put the centroid at the true face middle.
    let positions = brep.face_outer_loop_sample_positions(face, 8);
    if positions.is_empty() {
        return fallback(Vec3::new(0.0, 0.0, 1.0));
    }

    // Find UV centroid
    let ((u_min, u_max), (v_min, v_max)) =
        popkern_core::geometry::uv_projection::bounded_parameter_range(&positions, &sof.geometry);
    let u_mid = (u_min + u_max) * 0.5;
    let v_mid = (v_min + v_max) * 0.5;

    let mut normal = sof.geometry.normal_at(u_mid, v_mid);
    if !sof.orientation {
        normal = normal * (-1.0);
    }
    let axis_len = normal.length();
    if axis_len < 1.0e-12 {
        return fallback(normal);
    }
    let axis = normal * (1.0 / axis_len);

    // Material-side probe from the 3D centroid of the face boundary.
    let n = positions.len() as f64;
    let (mut cx, mut cy, mut cz) = (0.0, 0.0, 0.0);
    for p in &positions {
        cx += p.x;
        cy += p.y;
        cz += p.z;
    }
    let centroid = Point3::new(cx / n, cy / n, cz / n);
    let eps = (tolerance.linear() * 100.0).max(1.0e-5);
    let probe_plus = Point3::new(
        centroid.x + axis.x * eps,
        centroid.y + axis.y * eps,
        centroid.z + axis.z * eps,
    );
    let probe_minus = Point3::new(
        centroid.x - axis.x * eps,
        centroid.y - axis.y * eps,
        centroid.z - axis.z * eps,
    );
    let plus_out = matches!(
        crate::classify::point_in_solid(&probe_plus, brep, aabbs, tolerance),
        Ok(PointClassification::Outside)
    );
    let minus_out = matches!(
        crate::classify::point_in_solid(&probe_minus, brep, aabbs, tolerance),
        Ok(PointClassification::Outside)
    );
    let at = Some(sof.geometry.point_at(u_mid, v_mid));
    match (plus_out, minus_out) {
        // Material on the -axis side: outward is +axis.
        (true, false) => Outward { normal: axis, at },
        // Material on the +axis side: outward is -axis.
        (false, true) => Outward { normal: -axis, at },
        // Inconclusive: keep the flag-adjusted normal.
        _ => Outward { normal, at },
    }
}

/// A face's outward normal and the surface point it was taken at
/// (`None` when the face gave no point: no surface or no boundary).
struct Outward {
    normal: Vec3,
    at: Option<Point3>,
}
