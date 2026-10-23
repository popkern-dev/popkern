//! PK-926: a tube whose outer wall coincides with a coaxial D-bore's
//! wall (equal radius) serves every operation in both orders. D's top
//! rim lies on the tube wall; the clip judged the record along it by
//! containment samples on the cap's trim, its flush-span guard tested
//! hugging against a 25-sample polyline whose chords sit ~2.6e-3 off
//! the arc, and kept only the 5 % of the rim across the seam. The
//! wall's record ring dead-ended, the graze filter dropped the flat's
//! cut, and the intersection and subtraction stranded record sides.

use popkern_bool::{BoolOperation, boolean};
use popkern_core::geometry::{Point3, Transform3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_mesh::{TessellationOptions, tessellate};

fn volume(brep: &BRep) -> f64 {
    let options = TessellationOptions {
        chord_height: 1.0e-4,
        ..Default::default()
    };
    let mesh = tessellate(brep, &options).expect("mesh");
    let mut v6 = 0.0;
    for t in &mesh.triangles {
        let [a, b, c] = t.map(|i| mesh.vertices[i].position);
        v6 += a.x * (b.y * c.z - b.z * c.y) - a.y * (b.x * c.z - b.z * c.x)
            + a.z * (b.x * c.y - b.y * c.x);
    }
    (v6 / 6.0).abs()
}

/// Area of the radius-`r` disc above the chord y = -0.3.
fn d_area(r: f64) -> f64 {
    let c = 0.3_f64;
    std::f64::consts::PI * r * r - (r * r * (c / r).acos() - c * (r * r - c * c).sqrt())
}

#[test]
fn a_tube_flush_in_a_d_bore_serves_every_operation() {
    let tol = Tolerance::default();
    let outer = BRep::cylinder(1.0, 1.0, tol).expect("outer");
    let bore = BRep::cylinder(0.9, 1.6, tol)
        .expect("bore")
        .transform(&Transform3::translation(0.0, 0.0, -0.3))
        .expect("placed");
    let tube = boolean(BoolOperation::Subtraction, outer, bore, tol, None).expect("tube");
    let disc = BRep::cylinder(1.0, 1.0, tol)
        .expect("disc")
        .transform(&Transform3::translation(0.0, 0.0, -0.5))
        .expect("placed");
    let square = [
        Point3::new(-2.0, -0.3, -1.0),
        Point3::new(2.0, -0.3, -1.0),
        Point3::new(2.0, 2.0, -1.0),
        Point3::new(-2.0, 2.0, -1.0),
    ];
    let mut slab = BRep::prism(&square, Vec3::new(0.0, 0.0, 1.0), 2.0, tol).expect("slab");
    slab.reverse_body_orientation();
    let d = boolean(BoolOperation::Intersection, disc, slab, tol, None).expect("D");

    let v_tube = std::f64::consts::PI * (1.0 - 0.81);
    let v_d = d_area(1.0);
    let v_both = 0.5 * (d_area(1.0) - d_area(0.9));
    let cases = [
        ("tube ∩ D", BoolOperation::Intersection, &tube, &d, v_both),
        ("D ∩ tube", BoolOperation::Intersection, &d, &tube, v_both),
        (
            "tube ∪ D",
            BoolOperation::Union,
            &tube,
            &d,
            v_tube + v_d - v_both,
        ),
        (
            "D ∪ tube",
            BoolOperation::Union,
            &d,
            &tube,
            v_tube + v_d - v_both,
        ),
        (
            "tube − D",
            BoolOperation::Subtraction,
            &tube,
            &d,
            v_tube - v_both,
        ),
        (
            "D − tube",
            BoolOperation::Subtraction,
            &d,
            &tube,
            v_d - v_both,
        ),
    ];
    for (name, op, a, b, want) in cases {
        let body =
            boolean(op, a.clone(), b.clone(), tol, None).unwrap_or_else(|e| panic!("{name}: {e}"));
        let errors = body.validate_all();
        assert!(errors.is_empty(), "{name}: {errors:?}");
        let v = volume(&body);
        assert!(
            (v - want).abs() < 3e-3 * want,
            "{name}: volume {v} vs exact {want}"
        );
    }
}
