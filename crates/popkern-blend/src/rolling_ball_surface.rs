//! The rolling ball's surface on a corner face (PK-861, PK-864): its
//! chart axis, the sphere oriented like the cap it replaces, and the
//! sense read from a neighbouring round.

use popkern_core::geometry::PortableMath;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry, SurfaceOnFace};
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::topology::{BRep, FaceHandle};

use super::{Ball, face_of, loop_half_edges};

/// The ball's chart axis for a corner whose cap faces `toward` from the
/// centre and reaches `reach` (radians) from that direction: +Z when
/// both of its poles stay clear of the cap by a margin, else
/// perpendicular to `toward`, so both poles sit 90° from the cap's
/// middle (a +Z chart put a pole ON one rim of each tetrahedron corner:
/// that rim's chart image jumps at the pole and its pcurve sat 8.5e-3
/// off), turned about `toward` so the chart's u seam (its `e1`
/// direction) lies behind the cap.
fn ball_chart_axis(toward: Vec3, reach: f64) -> Vec3 {
    let Some(d) = toward.normalized() else {
        return Vec3::Z;
    };
    let margin = 10f64.to_radians();
    let to_pole =
        d.z.clamp(-1.0, 1.0)
            .pacos()
            .min((-d.z).clamp(-1.0, 1.0).pacos());
    if to_pole > reach + margin {
        return Vec3::Z;
    }
    let (p, q) = popkern_core::geometry::analytic_surfaces::orthonormal_basis(d);
    let mut best = (f64::INFINITY, Vec3::Z);
    for k in 0..64 {
        let angle = std::f64::consts::TAU * f64::from(k) / 64.0;
        let axis = p * angle.pcos() + q * angle.psin();
        let (e1, _) = popkern_core::geometry::analytic_surfaces::orthonormal_basis(axis);
        let facing = e1.dot(&d);
        if facing < best.0 {
            best = (facing, axis);
        }
    }
    best.1
}

/// Give `face` the ball's sphere, oriented like the face it replaces
/// (its normal at `probe` away from the centre or toward it, as
/// `outward_from_center` says), and fresh pcurves.
pub(crate) fn set_ball_surface(
    brep: &mut BRep,
    face: FaceHandle,
    ball: &Ball,
    probe: Point3,
    outward_from_center: bool,
) {
    // The cap's middle and reach, from its boundary vertices.
    let dirs: Vec<Vec3> = loop_half_edges(brep, face)
        .into_iter()
        .filter_map(|he| {
            brep.half_edges
                .get(he)
                .and_then(|h| brep.vertices.get(h.origin))
        })
        .filter_map(|v| (v.position - ball.center).normalized())
        .collect();
    let toward = dirs.iter().fold(Vec3::new(0.0, 0.0, 0.0), |a, &d| a + d);
    let reach = toward.normalized().map_or(std::f64::consts::PI, |m| {
        dirs.iter()
            .map(|d| d.dot(&m).clamp(-1.0, 1.0).pacos())
            .fold(0.0, f64::max)
    });
    let sphere = SurfaceGeometry::Sphere {
        center: ball.center,
        radius: ball.radius,
        axis: ball_chart_axis(toward, reach),
    };
    let chart_outward = sphere
        .closest_point(&probe)
        .map(|(u, v)| sphere.normal_at(u, v).dot(&(probe - ball.center)) > 0.0)
        .unwrap_or(true);
    brep.faces[face].surface = Some(SurfaceOnFace {
        geometry: sphere,
        orientation: chart_outward == outward_from_center,
    });
    for he in loop_half_edges(brep, face) {
        if let Some(h) = brep.half_edges.get_mut(he) {
            h.pcurve = None;
        }
    }
    brep.face_use_pcurves.retain(|key, _| key.face != face);
    brep.recompute_pcurves_for_lenient(&[face]);
}

/// A point of `face`'s boundary and whether a ball cap there faces away
/// from the ball's centre, read from a round across the boundary: on a
/// ball corner's boundary the round's oriented normal is the ball's
/// (G1 by construction) and the arc lies in the round's cross-section
/// through the centre, so the round pointing away from its axis means
/// the cap points away from the centre.
pub(crate) fn round_sense(brep: &BRep, face: FaceHandle) -> Option<(Point3, bool)> {
    for he in loop_half_edges(brep, face) {
        let h = brep.half_edges.get(he)?;
        let other = face_of(brep, h.twin)?;
        let sof = brep.faces.get(other)?.surface.as_ref()?;
        let SurfaceGeometry::Cylinder { origin, axis, .. } = &sof.geometry else {
            continue;
        };
        let probe = brep.vertices.get(h.origin)?.position;
        let (u, v) = sof.geometry.closest_point(&probe)?;
        let n = sof.geometry.normal_at(u, v) * if sof.orientation { 1.0 } else { -1.0 };
        let a = axis.normalized()?;
        let d = probe - *origin;
        let radial = d - a * d.dot(&a);
        return Some((probe, n.dot(&radial) > 0.0));
    }
    None
}
