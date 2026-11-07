use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};
use crate::geometry::{Point3, Vec3};

use super::{cross_section_curvature, runs};

fn at(surface: &SurfaceGeometry, p: Point3) -> (f64, f64) {
    surface.closest_point(&p).expect("on the surface")
}

#[test]
fn a_plane_has_no_cross_section_curvature() {
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    let uv = at(&plane, Point3::new(0.3, 0.0, 0.0));
    let k = cross_section_curvature(
        &plane,
        uv,
        Vec3::new(0.0, 1.0, 0.0),
        Vec3::new(1.0, 0.0, 0.0),
    );
    assert!(k.is_some_and(|k| k.abs() < 1e-12), "{k:?}");
}

/// A cylinder of radius r resting on the x axis bends toward +z at
/// 1/r on either side; the sign follows `t × d`.
#[test]
fn a_resting_cylinder_bends_at_one_over_its_radius() {
    let r = 0.25;
    let cyl = SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, r),
        axis: Vec3::new(1.0, 0.0, 0.0),
        radius: r,
    };
    let uv = at(&cyl, Point3::new(0.5, 0.0, 0.0));
    let t = Vec3::new(1.0, 0.0, 0.0);
    // d = +y: t × d = +z, the bend: positive.
    let up = cross_section_curvature(&cyl, uv, Vec3::new(0.0, 1.0, 0.0), t).expect("k");
    // d = −y: t × d = −z, against the bend: negative.
    let down = cross_section_curvature(&cyl, uv, Vec3::new(0.0, -1.0, 0.0), t).expect("k");
    assert!((up - 1.0 / r).abs() < 1e-9, "{up}");
    assert!((down + 1.0 / r).abs() < 1e-9, "{down}");
    // Along the rulings (the axis direction) it does not bend.
    let ruling =
        cross_section_curvature(&cyl, uv, Vec3::new(1.0, 0.0, 0.0), Vec3::new(0.0, 1.0, 0.0))
            .expect("k");
    assert!(ruling.abs() < 1e-9, "{ruling}");
}

#[test]
fn a_sphere_bends_at_one_over_its_radius_in_every_direction() {
    let r = 3.0;
    let sphere = SurfaceGeometry::Sphere {
        center: Point3::new(0.0, 0.0, r),
        radius: r,
        axis: Vec3::new(1.0, 0.0, 0.0),
    };
    let uv = at(&sphere, Point3::new(0.0, 0.0, 0.0));
    for d in [Vec3::new(1.0, 0.0, 0.0), Vec3::new(0.0, 1.0, 0.0)] {
        let t = Vec3::new(0.0, 0.0, 1.0).cross(&d);
        let k = cross_section_curvature(&sphere, uv, d, t).expect("k");
        assert!((k.abs() - 1.0 / r).abs() < 1e-9, "{k}");
    }
}

#[test]
fn runs_split_where_consecutive_keys_part() {
    let got = runs(vec![0.0, 0.05, 0.1, 1.0, 1.02, 3.0], |&x| x, 0.06);
    assert_eq!(got, vec![vec![0.0, 0.05, 0.1], vec![1.0, 1.02], vec![3.0]]);
}
