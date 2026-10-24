//! face_area of trimmed cylinder, sphere and cone faces.
//!
//! PK-934: a face without inner loops read the closed form of its
//! pcurves' uv BOUNDING rectangle: the notched cone wall of a cone minus
//! a corner box read 8.781 for a lateral area of at most 7.025, the
//! octant-cut sphere half of its true 7/8. The closed form now applies
//! only when the outer loop bounds the rectangle.
//!
//! PK-935: any other trim, holes included, integrates its area exactly
//! over its boundary (it read the measuring mesh, 1.4e-4 relative off),
//! in every orientation. PK-939: torus faces too.

use popkern_bool::{BoolOperation, boolean};
use popkern_core::geometry::Transform3;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_measure::face_area;

const PI: f64 = std::f64::consts::PI;

fn curved_area(body: &BRep, tol: &Tolerance, pick: fn(&SurfaceGeometry) -> bool) -> f64 {
    body.faces
        .iter()
        .filter(|(_, f)| f.surface.as_ref().is_some_and(|s| pick(&s.geometry)))
        .map(|(h, _)| face_area(body, h, tol).expect("area"))
        .sum()
}

fn is_cone(g: &SurfaceGeometry) -> bool {
    matches!(g, SurfaceGeometry::Cone { .. })
}
fn is_sphere(g: &SurfaceGeometry) -> bool {
    matches!(g, SurfaceGeometry::Sphere { .. })
}
fn is_cylinder(g: &SurfaceGeometry) -> bool {
    matches!(g, SurfaceGeometry::Cylinder { .. })
}
fn is_torus(g: &SurfaceGeometry) -> bool {
    matches!(g, SurfaceGeometry::Torus { .. })
}

fn rotations() -> Vec<(&'static str, Transform3)> {
    vec![
        ("unrotated", Transform3::identity()),
        ("z 0.7", Transform3::rotation_z(0.7)),
        ("x 1.1", Transform3::rotation_x(1.1)),
        ("y -0.4", Transform3::rotation_y(-0.4)),
        (
            "mixed",
            Transform3::rotation_x(0.3)
                .compose(&Transform3::rotation_y(1.9))
                .compose(&Transform3::rotation_z(-2.2)),
        ),
    ]
}

fn op(kind: BoolOperation, a: BRep, b: BRep, r: &Transform3, tol: Tolerance) -> BRep {
    let a = a.transform(r).expect("rotate");
    let b = b.transform(r).expect("rotate");
    boolean(kind, a, b, tol, None).expect("boolean")
}

fn assert_close(what: &str, got: f64, want: f64) {
    assert!((got - want).abs() < 1e-10 * want, "{what}: {got} vs {want}");
}

#[test]
fn notched_and_octant_faces_measure_their_region_exactly() {
    let tol = Tolerance::default();
    for (name, r) in rotations() {
        // Cone (radius 1, height 2) less the corner box [0,1]^2 x
        // [0.8,1.8]: the lateral area less a quarter of the frustum band
        // between radii 0.6 and 0.1 (slant √1.25).
        let cube = BRep::unit_cube(tol)
            .and_then(|c| c.transform(&Transform3::translation(0.0, 0.0, 0.8)))
            .expect("box");
        let cone = BRep::cone(1.0, 2.0, tol).expect("cone");
        let notched = op(BoolOperation::Subtraction, cone, cube, &r, tol);
        let want = PI * 5.0_f64.sqrt() - PI * 0.7 * 1.25_f64.sqrt() / 4.0;
        assert_close(
            &format!("{name}: cone wall"),
            curved_area(&notched, &tol, is_cone),
            want,
        );

        // A box and a sphere on its corner: 7/8 of the sphere shows on
        // the union, the inward octant on the difference.
        let ball = || {
            BRep::sphere(0.5, tol)
                .and_then(|s| s.transform(&Transform3::translation(1.0, 1.0, 1.0)))
                .expect("sphere")
        };
        let cube = || BRep::unit_cube(tol).expect("box");
        let united = op(BoolOperation::Union, cube(), ball(), &r, tol);
        assert_close(
            &format!("{name}: 7/8 sphere"),
            curved_area(&united, &tol, is_sphere),
            7.0 * PI / 8.0,
        );
        let cut = op(BoolOperation::Subtraction, cube(), ball(), &r, tol);
        assert_close(
            &format!("{name}: inward octant"),
            curved_area(&cut, &tol, is_sphere),
            PI / 8.0,
        );
    }
}

#[test]
fn holed_cylinder_and_sphere_faces_measure_exactly() {
    let tol = Tolerance::default();
    for (name, r) in rotations() {
        // A box pocketed into a cylinder wall (radius 0.5) at u = π/2:
        // the hole spans |x| < 0.1, i.e. 2·asin(0.2) of azimuth, over
        // height 0.2.
        let cylinder = BRep::cylinder(0.5, 2.0, tol).expect("cylinder");
        let pocket = BRep::unit_cube(tol)
            .and_then(|c| c.transform(&Transform3::scale_uniform(0.2)))
            .and_then(|c| c.transform(&Transform3::translation(-0.1, 0.4, 0.5)))
            .expect("pocket");
        let holed = op(BoolOperation::Subtraction, cylinder, pocket, &r, tol);
        let want = 2.0 * PI * 0.5 * 2.0 - 0.5 * 2.0 * 0.2_f64.asin() * 0.2;
        assert_close(
            &format!("{name}: holed wall"),
            curved_area(&holed, &tol, is_cylinder),
            want,
        );

        // A sphere (radius 0.5) with a cap cut off by the plane x = 0.3:
        // its hole is a circle off the poles, and the cap has height 0.2.
        let ball = BRep::sphere(0.5, tol).expect("sphere");
        let slab = BRep::unit_cube(tol)
            .and_then(|c| c.transform(&Transform3::translation(0.3, -0.5, -0.5)))
            .expect("slab");
        let capped = op(BoolOperation::Subtraction, ball, slab, &r, tol);
        let want = PI - 2.0 * PI * 0.5 * 0.2;
        assert_close(
            &format!("{name}: capped sphere"),
            curved_area(&capped, &tol, is_sphere),
            want,
        );
    }
}

/// A sphere (radius 0.5) drilled through its centre by a radius-0.2
/// bore: two caps of height 0.5 − √0.21 leave it (PK-937). Across the
/// sphere's pole axis the hole circles were marched, 7e-8 off both
/// surfaces, and the area read 2e-7 high; along it they were iso-rows.
/// A hole straddling the sphere's seam meridian (turn 0 meets the bore
/// seam there too, turn 0.2 does not) also left the cellular path:
/// the kept face sampled its classification point inside its hole.
#[test]
fn a_drilled_sphere_measures_exactly() {
    let tol = Tolerance::default();
    let want = PI - 2.0 * 2.0 * PI * 0.5 * (0.5 - 0.21_f64.sqrt());
    let across =
        |turn: f64| Transform3::rotation_z(turn).compose(&Transform3::rotation_y(PI / 2.0));
    for (bore_axis, tilt) in [
        ("along the pole", Transform3::identity()),
        ("across, on the seam", across(0.0)),
        ("across, straddling the seam", across(0.2)),
        ("across, clear of the seam", across(1.0)),
    ] {
        for (name, r) in rotations() {
            let ball = BRep::sphere(0.5, tol).expect("sphere");
            let bore = BRep::cylinder(0.2, 2.0, tol)
                .and_then(|c| c.transform(&Transform3::translation(0.0, 0.0, -1.0)))
                .and_then(|c| c.transform(&tilt))
                .expect("bore");
            let _ = popkern_bool::merge_fallback::take_merge_fallback_events();
            let drilled = op(BoolOperation::Subtraction, ball, bore, &r, tol);
            let fallbacks = popkern_bool::merge_fallback::take_merge_fallback_events();
            assert!(
                fallbacks.is_empty(),
                "{bore_axis}, {name}: left the cellular path: {fallbacks:?}"
            );
            let got = curved_area(&drilled, &tol, is_sphere);
            assert!(
                (got - want).abs() < 1e-10 * want,
                "{bore_axis}, {name}: drilled sphere {got} vs {want}"
            );
        }
    }
}

#[test]
fn a_rectangle_bounded_face_keeps_its_closed_form() {
    let tol = Tolerance::default();
    let cylinder = BRep::cylinder(0.5, 2.0, tol).expect("cylinder");
    let got = curved_area(&cylinder, &tol, is_cylinder);
    let want = 2.0 * PI * 0.5 * 2.0;
    assert!(
        (got - want).abs() < 1e-12 * want,
        "cylinder wall {got} vs {want}"
    );
}

/// The torus of [`trimmed_torus_faces_measure_exactly`]: major and
/// minor radius.
const TORUS: (f64, f64) = (1.0, 0.3);

/// The area of the torus region `cos u > k(v)` over `v ∈ [a, b]` (a
/// hole about u = 0), `∫ r(R + r cos v)·2 acos k dv`. `k` reaches 1 at
/// both ends, where the integrand has a square-root edge; the
/// substitution `v = a + (b − a)(1 − cos θ)/2` smooths it, and
/// Gauss–Legendre on 400 panels resolves it to roundoff.
fn torus_hole(a: f64, b: f64, k: &dyn Fn(f64) -> f64) -> f64 {
    const GL5: [(f64, f64); 5] = [
        (-0.906_179_845_938_664, 0.236_926_885_056_189_1),
        (-0.538_469_310_105_683, 0.478_628_670_499_366_5),
        (0.0, 0.568_888_888_888_888_9),
        (0.538_469_310_105_683, 0.478_628_670_499_366_5),
        (0.906_179_845_938_664, 0.236_926_885_056_189_1),
    ];
    let (big_r, r) = TORUS;
    let panels = 400;
    let mut sum = 0.0;
    for i in 0..panels {
        let (t0, t1) = (
            PI * i as f64 / panels as f64,
            PI * (i + 1) as f64 / panels as f64,
        );
        for &(x, w) in &GL5 {
            let theta = 0.5 * (t0 + t1) + 0.5 * (t1 - t0) * x;
            let v = a + (b - a) * (1.0 - theta.cos()) / 2.0;
            let dv = (b - a) * theta.sin() / 2.0;
            let width = 2.0 * k(v).clamp(-1.0, 1.0).acos();
            sum += 0.5 * (t1 - t0) * w * r * (big_r + r * v.cos()) * width * dv;
        }
    }
    sum
}

/// PK-939: a trimmed torus face read its measuring mesh (~1e-4). It
/// now integrates its boundary, its pcurves lifted in the periodic v.
/// A plane through the centre, tilted 0.3 off the axis, halves the
/// torus by point symmetry: exact to roundoff on its exact trims. A
/// radius-0.15 bore parallel to the axis through the tube, and a
/// radius-0.2 sphere centred on the outer equator, leave holes whose
/// trims are marched (samples within 9e-7 of the surfaces, inside the
/// linear tolerance): the area holds their accuracy, 3e-10 relative
/// at worst, against a quadrature reference. Under rotation_y(-0.4)
/// and the mixed rotation the bore and bite read 1.5e-7 and 1.1e-6 off
/// until the hooked trims of PK-941 were fixed.
#[test]
fn trimmed_torus_faces_measure_exactly() {
    let tol = Tolerance::default();
    let (big_r, r) = TORUS;
    let whole = 4.0 * PI * PI * r * big_r;
    let rho = 0.15;
    let bore_k = |v: f64| {
        let p = big_r + r * v.cos();
        (p * p + big_r * big_r - rho * rho) / (2.0 * p * big_r)
    };
    let bore = 2.0 * torus_hole((rho / r).acos(), (-rho / r).acos(), &bore_k);
    let bite_r = 0.2;
    let bite_k = |v: f64| {
        let p = big_r + r * v.cos();
        let reach = (big_r + r).powi(2) + big_r * big_r + 2.0 * big_r * r * v.cos() + r * r;
        (reach - bite_r * bite_r) / (2.0 * p * (big_r + r))
    };
    let edge = 2.0 * (bite_r / (2.0 * r)).asin();
    let bite = torus_hole(-edge, edge, &bite_k);
    for (name, rot) in rotations() {
        let torus = || BRep::torus(big_r, r, tol).expect("torus");
        let half_space = BRep::unit_cube(tol)
            .and_then(|c| c.transform(&Transform3::scale_uniform(4.0)))
            .and_then(|c| c.transform(&Transform3::translation(0.0, -2.0, -2.0)))
            .and_then(|c| c.transform(&Transform3::rotation_y(0.3)))
            .expect("half space");
        let halved = op(BoolOperation::Subtraction, torus(), half_space, &rot, tol);
        let got = curved_area(&halved, &tol, is_torus);
        assert!(
            (got - 0.5 * whole).abs() < 1e-12 * whole,
            "{name}: halved torus {got} vs {}",
            0.5 * whole
        );
        let drill = BRep::cylinder(rho, 2.0, tol)
            .and_then(|c| c.transform(&Transform3::translation(big_r, 0.0, -1.0)))
            .expect("drill");
        let bitten_by = BRep::sphere(bite_r, tol)
            .and_then(|s| s.transform(&Transform3::translation(big_r + r, 0.0, 0.0)))
            .expect("sphere");
        for (case, tool, want) in [
            ("drilled", drill, whole - bore),
            ("bitten", bitten_by, whole - bite),
        ] {
            let cut = op(BoolOperation::Subtraction, torus(), tool, &rot, tol);
            let got = curved_area(&cut, &tol, is_torus);
            assert!(
                (got - want).abs() < 1e-9 * want,
                "{name}: {case} torus {got} vs {want}"
            );
        }
    }
}
