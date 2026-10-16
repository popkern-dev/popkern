//! §8.3 blend cap geometry constructors.
//!
//! When a filleted edge terminates at a **free boundary** — a place
//! where the blend surface would otherwise end with no neighbouring
//! blend to meet — the engine closes the end with a dedicated cap
//! surface. The four cap kinds defined in §8.3 each address a
//! different continuity requirement:
//!
//! | Kind | Continuity with blend | Continuity with supports | Surface |
//! |---|---|---|---|
//! | [`CapKind::Open`] | — (no cap) | — | none |
//! | [`CapKind::Flat`] | G0 | G0 | `Plane` |
//! | [`CapKind::Rounded`] | G1 | G0 | `Sphere` |
//! | [`CapKind::TangentMatch`] | G2 on cylinder strips, G1 otherwise | G1 at the contacts | `Nurbs` (revolved quartic) |
//! | [`CapKind::UserDefined`] | caller's choice | caller's choice | caller |
//!
//! This module ships the **geometry constructors** — pure functions
//! that take the blend-endpoint contact data and return a
//! [`CapGeometry`] ready for topology insertion. The companion
//! **topology insertion** pass that wires the cap face into the
//! BRep is deferred to a follow-up commit; it needs detection logic
//! for "this blend endpoint terminates at a free boundary" plus
//! non-trivial topology surgery on the blend-face outer loop.
//!
//! # Inputs
//!
//! Every constructor takes a [`CapEndpoint`] struct built from the
//! blend's compute output at the terminating spine parameter. It
//! carries:
//!
//! - the 3D contact points on each support face (where the rolling
//!   ball touches)
//! - the unit normal of each support face at those points
//! - the rolling-ball centre at the endpoint
//! - the unit tangent to the centre curve (i.e. the spine direction
//!   leaving the blend toward the open end)
//! - the local minor radius `r`
//!
//! The blend surface's own end trim curve is not part of the input;
//! the cap boundary is built directly from the two contact points
//! and the arc on the sphere of radius `r` around the centre point.

use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{Point3, Vec3};
use popkern_math::PortableMath;

use crate::types::CapKind;

/// Blend-endpoint data required to build any cap kind.
///
/// Every field is in world coordinates. Unit-length vectors are
/// documented as unit-length but not runtime-checked — pass
/// normalised inputs.
#[derive(Debug, Clone, Copy)]
pub struct CapEndpoint {
    /// 3D contact point on support face 1 (where the rolling ball
    /// touches when it's at the spine endpoint).
    pub contact_1: Point3,
    /// 3D contact point on support face 2.
    pub contact_2: Point3,
    /// Unit outward normal of support face 1 at `contact_1`.
    pub normal_1: Vec3,
    /// Unit outward normal of support face 2 at `contact_2`.
    pub normal_2: Vec3,
    /// Rolling-ball centre at the endpoint (midpoint of the two
    /// contact points offset by `r` along each support's inward
    /// normal).
    pub center: Point3,
    /// Unit tangent to the spine/centre curve at the endpoint,
    /// pointing *away* from the blend interior (toward the open
    /// end). For a cap at `t = 0` this is the negative of the
    /// spine's start tangent; for a cap at `t = 1` it is the
    /// spine's end tangent.
    pub spine_tangent_outward: Vec3,
    /// Minor radius of the blend at the endpoint — the rolling-
    /// ball radius. Constant specs pass their single radius;
    /// variable specs pass the radius at the terminating parameter.
    pub radius: f64,
}

/// A constructed cap ready for topology insertion.
///
/// The `boundary_points` vector traces the 3D boundary of the cap
/// face in walking order, starting at `contact_1` and ending at
/// `contact_2`. The first element equals `contact_1` and the last
/// equals `contact_2`; intermediate points follow the circular arc
/// from `contact_1` to `contact_2` along the sphere of radius `r`
/// around `center`. Callers that need finer boundary resolution
/// call [`resample_cap_boundary`] with a higher sample count.
#[derive(Debug, Clone)]
pub struct CapGeometry {
    /// Which cap kind was produced.
    pub kind: CapKind,
    /// The cap's surface geometry. `None` for [`CapKind::Open`].
    pub surface: Option<SurfaceGeometry>,
    /// 3D polyline of the cap face's outer boundary, from
    /// `contact_1` to `contact_2`. At least 2 points (the two
    /// contact endpoints).
    pub boundary_points: Vec<Point3>,
    /// Unit spine tangent at the endpoint, pointing away from the
    /// blend interior: the normal of the strip's end section, which
    /// the inserted cap's boundary with the strip follows (PK-532).
    pub spine_tangent_outward: Vec3,
}

/// Dispatch to the cap constructor matching the requested kind.
///
/// Returns an empty [`CapGeometry`] for [`CapKind::Open`]. Returns
/// `None` for [`CapKind::UserDefined`] (not implemented) — the
/// caller is expected to handle the fallback.
pub fn build_cap(kind: CapKind, endpoint: &CapEndpoint) -> Option<CapGeometry> {
    match kind {
        CapKind::Open => Some(CapGeometry {
            kind,
            surface: None,
            boundary_points: vec![endpoint.contact_1, endpoint.contact_2],
            spine_tangent_outward: endpoint.spine_tangent_outward,
        }),
        CapKind::Flat => Some(build_flat_cap(endpoint)),
        CapKind::Rounded => Some(build_rounded_cap(endpoint)),
        CapKind::TangentMatch => Some(build_tangent_match_cap(endpoint)),
        CapKind::UserDefined => None,
    }
}

/// Build a planar cap perpendicular to the spine tangent at the
/// blend endpoint.
///
/// The cap's plane passes through the midpoint of `contact_1` and
/// `contact_2` with normal aligned to `spine_tangent_outward`. The
/// boundary is the straight segment from `contact_1` through the
/// centre to `contact_2` plus the arc from `contact_2` back to
/// `contact_1` along the sphere-on-plane intersection — in practice,
/// for the cap's topology role, we emit just the two contact points
/// and the centre as three collinear-in-plane points; the real
/// boundary is closed by the blend's own trim curve on the other
/// side.
pub fn build_flat_cap(endpoint: &CapEndpoint) -> CapGeometry {
    let midpoint = midpoint_of(&endpoint.contact_1, &endpoint.contact_2);
    // Plane normal = spine tangent (by definition of "perpendicular
    // to the spine"). The u_axis is arbitrary but chosen to point
    // from contact_1 toward contact_2 for deterministic
    // parameterisation.
    let u_axis = (endpoint.contact_2 - endpoint.contact_1)
        .normalized()
        .unwrap_or(Vec3::new(1.0, 0.0, 0.0));
    let v_axis = endpoint
        .spine_tangent_outward
        .cross(&u_axis)
        .normalized()
        .unwrap_or(Vec3::new(0.0, 1.0, 0.0));

    let plane = SurfaceGeometry::Plane {
        origin: midpoint,
        normal: endpoint.spine_tangent_outward,
        u_axis,
        v_axis,
    };

    CapGeometry {
        kind: CapKind::Flat,
        surface: Some(plane),
        boundary_points: vec![endpoint.contact_1, midpoint, endpoint.contact_2],
        spine_tangent_outward: endpoint.spine_tangent_outward,
    }
}

/// Build a spherical cap of radius equal to the blend's local
/// minor radius, centred at the rolling-ball centre.
///
/// The sphere surface continues the blend's G1 curvature at the
/// endpoint — a rolling ball of radius `r` sitting at `center`
/// touches both supports at exactly the input contact points, so
/// the sphere and the blend surface share their normal at those
/// contacts.
///
/// The boundary points traces the arc from `contact_1` to
/// `contact_2` along the sphere surface, sampled at 8 intermediate
/// points.
pub fn build_rounded_cap(endpoint: &CapEndpoint) -> CapGeometry {
    let sphere = SurfaceGeometry::Sphere {
        center: endpoint.center,
        radius: endpoint.radius,
        axis: popkern_core::geometry::Vec3::Z,
    };

    let boundary_points = sample_spherical_arc(
        &endpoint.contact_1,
        &endpoint.contact_2,
        &endpoint.center,
        endpoint.radius,
        8,
    );

    CapGeometry {
        kind: CapKind::Rounded,
        surface: Some(sphere),
        boundary_points,
        spine_tangent_outward: endpoint.spine_tangent_outward,
    }
}

/// Build the tangent-match cap: a surface of revolution about the end
/// tangent line through the ball centre, whose quartic meridian
/// leaves the end section along the strip's ruling with zero
/// curvature and closes flat at a pole r beyond the station (see
/// [`tangent_match`] for the construction and continuity argument).
///
/// Its boundary with the strip is the section arc, the same great
/// circle the Rounded sphere shares with the strip. A degenerate
/// endpoint yields no surface, and the cap is then not inserted.
pub fn build_tangent_match_cap(endpoint: &CapEndpoint) -> CapGeometry {
    CapGeometry {
        kind: CapKind::TangentMatch,
        surface: tangent_match::tangent_match_surface(endpoint).map(SurfaceGeometry::Nurbs),
        boundary_points: sample_spherical_arc(
            &endpoint.contact_1,
            &endpoint.contact_2,
            &endpoint.center,
            endpoint.radius,
            8,
        ),
        spine_tangent_outward: endpoint.spine_tangent_outward,
    }
}

/// Resample the boundary polyline of a rounded or tangent-match cap
/// (its section arc) at a different density. Flat caps return the
/// contacts and their midpoint; the rest return the contacts.
pub fn resample_cap_boundary(
    kind: CapKind,
    endpoint: &CapEndpoint,
    n_samples: usize,
) -> Vec<Point3> {
    let n = n_samples.max(2);
    match kind {
        CapKind::Rounded | CapKind::TangentMatch => sample_spherical_arc(
            &endpoint.contact_1,
            &endpoint.contact_2,
            &endpoint.center,
            endpoint.radius,
            n,
        ),
        CapKind::Flat | CapKind::Open | CapKind::UserDefined => {
            // Flat / Open caps have no meaningful arc to sample;
            // return the two endpoints (plus midpoint for Flat).
            if kind == CapKind::Flat {
                vec![
                    endpoint.contact_1,
                    midpoint_of(&endpoint.contact_1, &endpoint.contact_2),
                    endpoint.contact_2,
                ]
            } else {
                vec![endpoint.contact_1, endpoint.contact_2]
            }
        }
    }
}

fn midpoint_of(a: &Point3, b: &Point3) -> Point3 {
    Point3::new(0.5 * (a.x + b.x), 0.5 * (a.y + b.y), 0.5 * (a.z + b.z))
}

/// Sample `n` points along the circular arc from `a` to `b` on a
/// sphere of radius `r` centred at `center`. The arc plane is
/// defined by the vectors `a - center` and `b - center`; the shorter
/// arc (< π) is chosen.
fn sample_spherical_arc(a: &Point3, b: &Point3, center: &Point3, r: f64, n: usize) -> Vec<Point3> {
    let n = n.max(2);
    let va = (*a - *center)
        .normalized()
        .unwrap_or(Vec3::new(1.0, 0.0, 0.0));
    let vb = (*b - *center)
        .normalized()
        .unwrap_or(Vec3::new(0.0, 1.0, 0.0));
    let cos_angle = va.dot(&vb).clamp(-1.0, 1.0);
    let angle = cos_angle.pacos();

    // Orthonormal basis {va, perp} in the arc plane.
    let perp = (vb - va * cos_angle).normalized().unwrap_or_else(|| {
        // `a` and `b` coincide or are antipodal — fall back to
        // any perpendicular to va.
        if va.x.abs() < 0.9 {
            va.cross(&Vec3::new(1.0, 0.0, 0.0))
                .normalized()
                .unwrap_or(Vec3::new(0.0, 1.0, 0.0))
        } else {
            va.cross(&Vec3::new(0.0, 1.0, 0.0))
                .normalized()
                .unwrap_or(Vec3::new(0.0, 0.0, 1.0))
        }
    });

    let mut out = Vec::with_capacity(n);
    for i in 0..n {
        let t = (i as f64) / ((n - 1) as f64);
        let theta = angle * t;
        let cos_t = theta.pcos();
        let sin_t = theta.psin();
        let dir = va * cos_t + perp * sin_t;
        out.push(Point3::new(
            center.x + dir.x * r,
            center.y + dir.y * r,
            center.z + dir.z * r,
        ));
    }
    out
}

#[path = "cap_tangent_match.rs"]
pub(crate) mod tangent_match;

#[cfg(test)]
#[path = "cap_tests.rs"]
mod tests;
