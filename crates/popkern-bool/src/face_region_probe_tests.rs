//! The probe answers exactly as the per-query containment it caches.

use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use super::FaceRegionProbe;
use crate::{BoolOperation, boolean};

/// A disc with eight radial slots: its caps are planar with holes.
fn slotted_disc() -> BRep {
    let tol = Tolerance::default();
    let mut body = BRep::cylinder(30.0, 1.5, tol).expect("disc");
    for k in 0..8 {
        let a = std::f64::consts::TAU * k as f64 / 8.0;
        let (ca, sa, w) = (a.cos(), a.sin(), 0.6);
        let pts = [
            Point3::new(22.0 * ca + w * sa, 22.0 * sa - w * ca, -1.0),
            Point3::new(27.0 * ca + w * sa, 27.0 * sa - w * ca, -1.0),
            Point3::new(27.0 * ca - w * sa, 27.0 * sa + w * ca, -1.0),
            Point3::new(22.0 * ca - w * sa, 22.0 * sa + w * ca, -1.0),
        ];
        let slot = BRep::prism(&pts, Vec3::new(0.0, 0.0, 1.0), 4.0, tol).expect("slot");
        body = boolean(BoolOperation::Subtraction, body, slot, tol, None).expect("cut");
    }
    body
}

#[test]
fn the_probe_agrees_with_the_per_query_containment_on_every_face() {
    let body = slotted_disc();
    let tol = body.tolerance;
    let mut planar_holed = 0;
    for (face, data) in body.faces.iter() {
        if !data.inner_loops.is_empty() {
            planar_holed += 1;
        }
        let probe = FaceRegionProbe::new(&body, face, &tol, 64);
        let bounds = crate::clip::face_uv_boundary(&body, face)
            .into_iter()
            .flatten()
            .flat_map(|(_, pc)| {
                let (t0, t1) = pc.parameter_range();
                (0..=8).map(move |i| pc.point_at(t0 + (t1 - t0) * i as f64 / 8.0))
            })
            .fold(None, |acc: Option<(f64, f64, f64, f64)>, q| {
                Some(acc.map_or((q.x, q.x, q.y, q.y), |(a, b, c, d)| {
                    (a.min(q.x), b.max(q.x), c.min(q.y), d.max(q.y))
                }))
            });
        let Some((u0, u1, v0, v1)) = bounds else {
            continue;
        };
        for i in 0..=40 {
            for j in 0..=40 {
                let q = Point3::new(
                    u0 + (u1 - u0) * (i as f64 / 40.0 * 1.2 - 0.1),
                    v0 + (v1 - v0) * (j as f64 / 40.0 * 1.2 - 0.1),
                    0.0,
                );
                assert_eq!(
                    probe.contains(q),
                    crate::split_inner_loop::point_inside_face_region_uv_with_density(
                        &body, face, q, &tol, 64
                    ),
                    "{face:?} at {q:?}"
                );
            }
        }
    }
    assert!(planar_holed >= 2, "both caps carry the slot holes");
}
