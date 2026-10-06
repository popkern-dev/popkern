//! Booleans of two overlapping cylinders with FLUSH caps: every cap
//! pair is a partial coplanar overlap (a lens) and every record of a
//! cap circle runs along the other wall's rim. Gated on VOLUME — the
//! face-count smoke test this family used to ride let a 0.77-volume
//! union (true 2.527) and a zero-volume subtraction serve as valid.

use popkern_bool::{intersection, subtraction, union};
use popkern_core::geometry::Transform3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

const R: f64 = 0.5;
const H: f64 = 2.0;
const OFFSET: f64 = 0.5;

/// Lens area of two radius-R disks OFFSET apart.
fn lens_area() -> f64 {
    let d = OFFSET;
    2.0 * R * R * (d / (2.0 * R)).acos() - 0.5 * d * (4.0 * R * R - d * d).sqrt()
}

fn operands(tol: Tolerance) -> (BRep, BRep) {
    let Ok(a) = BRep::cylinder(R, H, tol) else {
        panic!("cylinder A");
    };
    let Ok(b) = BRep::cylinder(R, H, tol)
        .and_then(|b| b.transform(&Transform3::translation(OFFSET, 0.0, 0.0)))
    else {
        panic!("cylinder B");
    };
    (a, b)
}

fn served_volume(body: &BRep, tol: &Tolerance) -> f64 {
    assert!(body.validate().is_ok(), "served body must validate");
    match popkern_measure::body_volume(body, tol) {
        Ok(v) => v,
        Err(e) => panic!("volume: {e}"),
    }
}

/// Tessellated volume agrees with the exact value to this fraction.
const VOLUME_FRACTION: f64 = 5e-3;

fn assert_volume(label: &str, got: f64, exact: f64) {
    assert!(
        (got - exact).abs() <= VOLUME_FRACTION * exact,
        "{label}: served volume {got} vs exact {exact}"
    );
}

#[test]
fn flush_cylinder_subtraction_serves_the_crescent_prism() {
    let disk = std::f64::consts::PI * R * R;
    for linear in [1e-4, 1e-6] {
        let tol = Tolerance::default().with_linear(linear);
        let (a, b) = operands(tol);
        match subtraction(a, b, tol, None) {
            Ok(r) => assert_volume(
                "subtraction",
                served_volume(&r, &tol),
                (disk - lens_area()) * H,
            ),
            Err(e) => panic!("subtraction declined at {linear:e}: {e}"),
        }
    }
}

#[test]
fn flush_cylinder_intersection_serves_the_lens_prism() {
    for linear in [1e-4, 1e-6] {
        let tol = Tolerance::default().with_linear(linear);
        let (a, b) = operands(tol);
        match intersection(a, b, tol, None) {
            Ok(r) => assert_volume("intersection", served_volume(&r, &tol), lens_area() * H),
            Err(e) => panic!("intersection declined at {linear:e}: {e}"),
        }
    }
}

#[test]
fn flush_cylinder_union_serves_both_disks_minus_one_lens() {
    let disk = std::f64::consts::PI * R * R;
    for linear in [1e-4, 1e-6] {
        let tol = Tolerance::default().with_linear(linear);
        let (a, b) = operands(tol);
        match union(a, b, tol, None) {
            Ok(r) => assert_volume(
                "union",
                served_volume(&r, &tol),
                (2.0 * disk - lens_area()) * H,
            ),
            Err(e) => panic!("union declined at {linear:e}: {e}"),
        }
    }
}
