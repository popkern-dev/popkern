//! Tilt-aware entry-ring cut tools (PK-480): the 45° revolved cones
//! that chamfer or countersink a hole mouth, shaped so that on a
//! tilted entry face only the cone and a flank-clear outer cylinder
//! ever meet material.
//!
//! The naive ring — a revolved triangle closed by a flat disk just
//! above the origin plane — is only correct for perpendicular
//! entries. On a tilted face its closing disk lies inside material on
//! the uphill side, where it grazes internal-thread flank tips at the
//! flank half-angle and amplifies the NURBS flank fit slop by
//! `1/sin α` (measured 1.6e-4..3.0e-4 curve-surface residuals against
//! the 8e-5 axis-aligned class — a torn mesh when tolerated), while
//! on the downhill side the cone stops short of the face and leaves
//! no chamfer at all.
//!
//! The tilt-aware ring is correct by construction instead:
//!
//! 1. the entry-face tilt is probed from the body (the planar face
//!    the hole axis pierces nearest the origin);
//! 2. the cone is plunged deeper by exactly the amount that restores
//!    the design rim radius at the steepest (downhill) azimuth, so
//!    the chamfer is nominal-or-wider at every azimuth;
//! 3. the tool is closed by an outer cylinder at the cone-top radius
//!    — strictly outside the thread-tip radius, so it never touches a
//!    flank — and a lid above the body's whole vertex hull, so the
//!    closure provably never meets material (curved or absent entry
//!    faces included).
//!
//! Every intersection class the tilted ring produces is either the
//! rigid image of a serving axis-aligned class (cone × flank, cone ×
//! bore circle — the ring is coaxial with the thread) or an analytic
//! quadric × plane section (cone/cylinder × tilted plane, the serving
//! bore-entry class). Perpendicular entries produce the identical cut
//! region: the plunge extension is zero and the closure faces lie
//! strictly outside the material half-space.

use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::topology::BRep;

use crate::FastenerError;

/// Entry-face tilt beyond which the ring declines: past ~40° the
/// 45° cone and the entry plane approach tangency along the uphill
/// rim (a near-tangent surface-intersection class), and a lead-in
/// chamfer cut coaxially into so steep a face is no longer a chamfer.
/// Kept as a tangent so the comparison is exact.
const MAX_ENTRY_TILT_TAN: f64 = 0.85;

/// How far the axis ray may pierce a candidate plane from the origin
/// (as a multiple of the caller's entry headroom) and still count as
/// the entry face.
const ENTRY_PIERCE_BAND: f64 = 1.0;

/// What the body looks like around the hole mouth: the measured
/// entry-face tilt and the height of the body above the origin, both
/// resolved against the hole axis.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct EntryProbe {
    /// Tangent of the angle between the hole axis and the entry-face
    /// normal (0 = perpendicular entry). Probed from the planar face
    /// the axis pierces nearest the origin; bodies with no such face
    /// (curved or absent entry) probe as perpendicular — the hull lid
    /// keeps the closure safe, only the plunge compensation is
    /// skipped.
    pub tilt_tan: f64,
    /// Height of the body's vertex hull above the origin measured
    /// against the axis direction (0 when the origin is at the top of
    /// the body). Places the closing lid provably outside material.
    pub hull_height: f64,
}

impl EntryProbe {
    /// A perpendicular entry at the top of the body (the legacy
    /// assumption; also the fallback when no entry face is found).
    pub fn perpendicular() -> EntryProbe {
        EntryProbe {
            tilt_tan: 0.0,
            hull_height: 0.0,
        }
    }
}

/// Probe the entry-face tilt and hull height for a hole cut into
/// `body` at `origin` along `direction` (into the material).
///
/// The entry face is the planar face whose plane the axis ray pierces
/// within `headroom` of the origin, nearest the origin; ties within
/// the model tolerance take the steeper plane (the conservative,
/// deeper-plunging read). Planes nearly parallel to the axis are
/// never candidates — a side wall close to the mouth must not
/// masquerade as the entry face.
pub(crate) fn probe_entry(
    body: &BRep,
    origin: Point3,
    direction: Vec3,
    headroom: f64,
) -> EntryProbe {
    let Some(axis) = direction.normalized() else {
        return EntryProbe::perpendicular();
    };
    let band = ENTRY_PIERCE_BAND * headroom;
    let tie = body.tolerance.linear();
    // (|pierce parameter|, tilt tangent) of the best candidate so far.
    let mut best: Option<(f64, f64)> = None;
    for (_, face) in body.faces.iter() {
        let Some(surface) = &face.surface else {
            continue;
        };
        let SurfaceGeometry::Plane {
            origin: po, normal, ..
        } = surface.geometry
        else {
            continue;
        };
        let cos = normal.dot(&axis);
        // Reject planes within ~8.6° of parallel to the axis: the
        // pierce parameter degenerates and no such plane is an entry.
        if cos.abs() < 0.15 {
            continue;
        }
        let t = (po - origin).dot(&normal) / cos;
        if !t.is_finite() || t.abs() > band {
            continue;
        }
        let tilt_tan = (1.0 - cos * cos).max(0.0).sqrt() / cos.abs();
        let closer = match best {
            None => true,
            Some((bt, btan)) => t.abs() + tie < bt || (t.abs() < bt + tie && tilt_tan > btan),
        };
        if closer {
            best = Some((t.abs(), tilt_tan));
        }
    }
    let tilt_tan = best.map_or(0.0, |(_, tan)| tan);
    let hull_height = body
        .vertices
        .iter()
        .map(|(_, v)| -(v.position - origin).dot(&axis))
        .fold(0.0_f64, f64::max);
    EntryProbe {
        tilt_tan,
        hull_height,
    }
}

/// Meridian polygon of a tilt-aware 45° entry ring, in the tool's
/// canonical frame (+Z into the material, origin at the hole mouth).
///
/// `r_in` is the inner cylinder radius (strictly inside the bore
/// void), `r_rim` the design rim radius — where the cone crosses the
/// origin plane on a perpendicular entry — and `eps` the transversal
/// overshoot margin. The cone keeps its 45° slope; the plunge is
/// deepened by `r_rim · tan θ` so the rim reaches `r_rim` at the
/// steepest azimuth, the outer cylinder sits one `eps` outside the
/// deepened cone top, and the lid sits `2·eps` above the probed hull.
///
/// # Errors
///
/// [`FastenerError::EntryTooSteep`] when the probed tilt exceeds
/// [`MAX_ENTRY_TILT_TAN`] — the cone/face intersection would go
/// near-tangent along the uphill rim.
pub(crate) fn entry_ring_profile(
    r_in: f64,
    r_rim: f64,
    eps: f64,
    probe: &EntryProbe,
    stage: &str,
) -> Result<Vec<(f64, f64)>, FastenerError> {
    if probe.tilt_tan > MAX_ENTRY_TILT_TAN {
        return Err(FastenerError::EntryTooSteep {
            stage: stage.to_string(),
            tilt_deg: probe.tilt_tan.atan().to_degrees(),
            limit_deg: MAX_ENTRY_TILT_TAN.atan().to_degrees(),
        });
    }
    let plunge = r_rim * probe.tilt_tan;
    let apex_sum = r_rim + plunge;
    let r_outer = apex_sum + eps;
    let z_lid = probe.hull_height + 2.0 * eps;
    Ok(vec![
        (r_in, -z_lid),
        (r_outer, -z_lid),
        (r_outer, -eps),
        (r_in, apex_sum - r_in),
    ])
}

#[cfg(test)]
#[path = "entry_ring_tests.rs"]
mod entry_ring_tests;
