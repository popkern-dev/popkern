use super::*;
use crate::geometry::intersection::contact_discovery::find_surface_contacts;

const SCALE: ContactScale = ContactScale {
    linear: 1e-6,
    extent: 10.0,
};

fn plane_z(z: f64) -> SurfaceGeometry {
    SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, z),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    }
}

/// The first tangent-curve contact `s1` over `w1` has with `s2`.
fn tangent_seed(
    s1: &SurfaceGeometry,
    w1: Window,
    s2: &SurfaceGeometry,
    w2: Window,
) -> SurfaceContact {
    find_surface_contacts(s1, w1, s2, w2, SCALE)
        .contacts
        .into_iter()
        .find(|c| {
            matches!(
                c.classification.kind,
                SurfaceContactKind::TangentCurve { .. }
            )
        })
        .expect("a tangent-curve contact")
}

/// A rod lying on a plane along its ruling, turned about its axis so its
/// chart's seam is off the line: the trace runs the ruling's length in
/// the plane's window, every point on the line.
#[test]
fn a_rod_on_a_plane_traces_its_ruling() {
    let plane = plane_z(0.0);
    let rod = SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 2.0, 0.5),
        axis: Vec3::new(1.0, 0.0, 0.0),
        radius: 0.5,
    };
    let (wp, wr) = (
        ((0.0, 4.0), (0.0, 4.0)),
        ((0.0, std::f64::consts::TAU), (0.5, 3.5)),
    );
    let seed = tangent_seed(&plane, wp, &rod, wr);
    let trace = trace_tangent_curve((&plane, wp), (&rod, wr), &seed, SCALE).expect("traced");
    assert!(!trace.closed);
    for p in &trace.points {
        assert!((p.y - 2.0).abs() < 1e-6 && p.z.abs() < 1e-9, "{p:?}");
    }
    let xs: Vec<f64> = trace.points.iter().map(|p| p.x).collect();
    let (lo, hi) = xs
        .iter()
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(a, b), &x| {
            (a.min(x), b.max(x))
        });
    assert!(lo < 0.5 + 0.2 && hi > 3.5 - 0.2, "the run {lo}..{hi}");
}

/// A torus lying on a plane touches it along a circle: the trace closes.
#[test]
fn a_torus_on_a_plane_traces_a_closed_circle() {
    let plane = plane_z(0.0);
    let (major, minor) = (2.0, 0.5);
    let torus = SurfaceGeometry::Torus {
        center: Point3::new(0.0, 0.0, minor),
        axis: Vec3::new(0.0, 0.0, 1.0),
        major_radius: major,
        minor_radius: minor,
    };
    let wp = ((-4.0, 4.0), (-4.0, 4.0));
    let seed = tangent_seed(&plane, wp, &torus, torus.parameter_range());
    let trace = trace_tangent_curve(
        (&plane, wp),
        (&torus, torus.parameter_range()),
        &seed,
        SCALE,
    )
    .expect("traced");
    assert!(trace.closed, "{} points", trace.points.len());
    for p in &trace.points {
        let r = (p.x * p.x + p.y * p.y).sqrt();
        assert!((r - major).abs() < 1e-5 && p.z.abs() < 1e-9, "{p:?}");
    }
}
