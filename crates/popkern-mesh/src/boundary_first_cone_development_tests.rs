use std::f64::consts::{PI, TAU};

use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::topology::EdgeHandle;

use super::{BoundaryRing, Development, development_for, ring_azimuths};

fn cone() -> SurfaceGeometry {
    SurfaceGeometry::Cone {
        apex: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        half_angle: 0.4,
    }
}

fn ring(surface: &SurfaceGeometry, uv: &[(f64, f64)]) -> BoundaryRing {
    BoundaryRing {
        uv: uv.to_vec(),
        pos: uv.iter().map(|&(u, v)| surface.point_at(u, v)).collect(),
        sources: uv
            .iter()
            .enumerate()
            .map(|(k, _)| (EdgeHandle::default(), k as f64, k as f64 + 1.0))
            .collect(),
    }
}

/// A face reaching a millimetre from the apex: a ruling down to
/// `v = 1e-3`, a conic turning round the apex there, a ruling back up and
/// the rim (the shape of a plane missing a cone's apex by a few
/// tolerances, PC-K).
fn sweeping_face(surface: &SurfaceGeometry) -> BoundaryRing {
    let mut uv = vec![(0.2, 1.0), (0.2, 0.5), (0.2, 1e-3)];
    uv.extend((1..=4).map(|k| (0.2 - 0.55 * k as f64, 1e-3)));
    uv.extend([(-2.0, 0.5), (-2.0, 1.0), (-1.0, 1.0)]);
    ring(surface, &uv)
}

fn development(surface: &SurfaceGeometry) -> Development<'_> {
    let rings = [sweeping_face(surface)];
    development_for(surface, &rings).map(|(d, _)| d).unwrap()
}

#[test]
fn the_development_is_an_isometry_keeping_the_face_charts_handedness() {
    let s = cone();
    let dev = development(&s);
    for &(u, v) in &[(0.1, 0.7), (-1.3, 0.02), (-0.4, 1e-4)] {
        let (x, y) = dev.develop((u, v));
        let (u2, v2) = dev.face_uv(x, y);
        assert!(
            (u2 - u).abs() < 1e-12 && (v2 - v).abs() < 1e-15,
            "round trip"
        );
        // Short steps along both parameters keep their length.
        for (du, dv) in [(1e-3, 0.0), (0.0, 1e-4)] {
            let (x2, y2) = dev.develop((u + du, v + dv));
            let planar = (x2 - x).hypot(y2 - y);
            let space = s.point_at(u, v).distance_to(&s.point_at(u + du, v + dv));
            assert!(
                ((planar - space) / space).abs() < 1e-3,
                "{planar} vs {space}"
            );
        }
        // The development's frame turns as the face chart's.
        let (px, py) = dev.first_partials(x, y);
        let (su, sv) = s.first_partials(u, v);
        assert!(px.cross(&py).dot(&su.cross(&sv)) > 0.0, "handedness");
    }
}

#[test]
fn a_ring_sweeping_round_the_apex_beside_it_develops() {
    let s = cone();
    let rings = [sweeping_face(&s)];
    let (dev, az) = development_for(&s, &rings).unwrap();
    // The azimuths are the face chart's, continuous along the ring.
    for (a, &(u, _)) in az[0].iter().zip(&rings[0].uv) {
        assert!((a - u).abs() < 1e-9, "{a} vs {u}");
    }
    assert!((dev.u_mid - (-0.9)).abs() < 1e-9);
}

#[test]
fn a_face_without_a_sweep_beside_the_apex_keeps_the_face_chart() {
    let s = cone();
    // The same face held back from the apex: its conic at v = 0.3.
    let mut uv = vec![(0.2, 1.0), (0.2, 0.3)];
    uv.extend((1..=4).map(|k| (0.2 - 0.55 * k as f64, 0.3)));
    uv.extend([(-2.0, 1.0), (-1.0, 1.0)]);
    assert!(development_for(&s, &[ring(&s, &uv)]).is_none());
    // A ruling sliver running onto the apex: no sweep either.
    let sliver = ring(&s, &[(0.0, 1.0), (0.0, 0.0), (0.01, 1.0)]);
    assert!(development_for(&s, &[sliver]).is_none());
}

#[test]
fn one_position_spread_along_the_apex_row_is_a_sweep() {
    let s = cone();
    // A ruling ends 1e-3 up, where the face chart spread its end along
    // the row in steps of 0.05 (each under the sweep step): one 3D point.
    let mut r = ring(&s, &[(0.2, 1.0), (0.2, 1e-3)]);
    let end = r.pos[1];
    for k in 1..=10 {
        r.uv.push((0.2 - 0.05 * k as f64, 1e-3));
        r.pos.push(end);
        r.sources.push((EdgeHandle::default(), 0.0, 0.0));
    }
    r.uv.push((0.3, 1.0));
    r.pos.push(s.point_at(0.3, 1.0));
    r.sources.push((EdgeHandle::default(), 0.0, 1.0));
    let (_, az) = development_for(&s, &[r]).unwrap();
    // Read from the position, the spread samples keep one azimuth.
    assert!(az[0][2..12].iter().all(|a| (a - 0.2).abs() < 1e-9));
}

fn frame() -> (Point3, Vec3, Vec3, Vec3) {
    let (e1, e2) =
        popkern_core::geometry::analytic_surfaces::orthonormal_basis(Vec3::new(0.0, 0.0, 1.0));
    (Point3::new(0.0, 0.0, 0.0), e1, e2, Vec3::new(0.0, 0.0, 1.0))
}

#[test]
fn a_ring_closing_round_the_axis_through_the_apex_turns_there() {
    let s = cone();
    // A face spanning a whole period between its seam's two uses: up the
    // seam, round the rim, down the seam to the apex.
    let mut uv = vec![(0.0, 0.0), (0.0, 1.0)];
    uv.extend((1..8).map(|k| (TAU * k as f64 / 8.0, 1.0)));
    uv.push((TAU, 1.0));
    let az = ring_azimuths(&ring(&s, &uv), frame(), 1e-9).unwrap();
    assert!(
        (az[1] - 0.0).abs() < 1e-9 && (az[9] - TAU).abs() < 1e-9,
        "{az:?}"
    );
    // Round the rim alone, never through the apex: it winds.
    let rim: Vec<_> = (0..8).map(|k| (TAU * k as f64 / 8.0, 1.0)).collect();
    assert!(ring_azimuths(&ring(&s, &rim), frame(), 1e-9).is_none());
}

#[test]
fn a_step_past_half_a_period_takes_the_pcurves_way() {
    let s = cone();
    // A conic turning 3.7 rad round the apex between two samples: from the
    // positions alone the nearer way is −2.58.
    let uv = [(0.0, 1.0), (0.0, 1e-3), (3.7, 1e-3), (3.7, 1.0), (1.8, 1.0)];
    let az = ring_azimuths(&ring(&s, &uv), frame(), 1e-9).unwrap();
    assert!((az[2] - 3.7).abs() < 1e-9, "{az:?}");
    assert!(az[2] - az[1] > PI);
}

#[test]
fn coincident_samples_merge_into_one_with_their_spans() {
    let s = cone();
    let p = s.point_at(0.0, 0.5);
    let q = s.point_at(0.5, 0.5);
    let e = EdgeHandle::default();
    let r = BoundaryRing {
        uv: vec![(0.0, 0.5), (0.1, 0.5), (0.5, 0.5), (0.9, 0.9)],
        pos: vec![p, p, q, s.point_at(0.9, 0.9)],
        sources: vec![
            (e, 0.0, 0.25),
            (e, 0.25, 0.5),
            (e, 0.5, 0.75),
            (e, 0.75, 1.0),
        ],
    };
    let m = super::mesh::merge_coincident(r);
    assert_eq!(m.pos.len(), 3);
    assert_eq!(m.sources[0], (e, 0.0, 0.5));
}
