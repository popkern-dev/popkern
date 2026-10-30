//! A record piece outside the unused face's trim is a phantom split,
//! not a leak: the keyway slot cut into a bored block serves through
//! the complex.

use popkern_core::geometry::transform::Transform3;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use crate::BoolOperation;

/// The axis-aligned box over `[x0, x1] × [y0, y1] × [z0, z1]`.
fn cuboid(x: (f64, f64), y: (f64, f64), z: (f64, f64), tol: Tolerance) -> BRep {
    let profile = [
        Point3::new(x.0, y.0, z.0),
        Point3::new(x.1, y.0, z.0),
        Point3::new(x.1, y.1, z.0),
        Point3::new(x.0, y.1, z.0),
    ];
    BRep::prism(&profile, Vec3::new(0.0, 0.0, 1.0), z.1 - z.0, tol).expect("cuboid")
}

/// A 1 × 1 × 0.3 block centred on the z axis, bored through by a
/// 0.15-radius hole.
fn bored_block(tol: Tolerance) -> BRep {
    let block = cuboid((-0.5, 0.5), (-0.5, 0.5), (0.0, 0.3), tol);
    let bore = BRep::cylinder(0.15, 0.5, tol)
        .expect("bore")
        .transform(&Transform3::translation(0.0, 0.0, -0.1))
        .expect("placed bore");
    crate::boolean(BoolOperation::Subtraction, block, bore, tol, None).expect("bored block")
}

/// A 0.08 × 0.06 slot over y ∈ [0.14, 0.20], through the block: its
/// inner wall y = 0.14 crosses the bore, so the records where that wall
/// meets the block's top and bottom faces run, between x = ±0.04,
/// across those faces' holes.
fn key(tol: Tolerance) -> BRep {
    cuboid((-0.04, 0.04), (0.14, 0.2), (-0.14, 0.44), tol)
}

#[test]
fn a_record_across_a_hole_does_not_decline_the_complex() {
    let tol = Tolerance::default();
    let block = bored_block(tol);
    let key = key(tol);
    let ix = crate::intersect::find_intersections(&block, &key, &tol).expect("intersections");
    let served = crate::complex::serve::rescue(&block, &key, &ix, BoolOperation::Subtraction, &tol)
        .expect("manifold");
    let (body, _) = served.expect("the complex serves the keyway cut");
    body.validate().expect("valid");
    // Material where the block stays, none in the bore or the slot.
    let boxes = crate::intersect::compute_face_aabbs(&body);
    for (p, inside) in [
        ((0.3, 0.0, 0.15), true),
        ((0.0, -0.3, 0.15), true),
        ((0.1, 0.17, 0.15), true),
        ((0.0, 0.0, 0.15), false),
        ((0.0, 0.17, 0.15), false),
        ((0.0, 0.17, 0.02), false),
        ((0.0, 0.17, -0.05), false),
    ] {
        let point = Point3::new(p.0, p.1, p.2);
        let got =
            crate::classify::point_in_solid_final(&point, &body, &boxes, &tol).expect("classifies");
        let want = if inside {
            crate::PointClassification::Inside
        } else {
            crate::PointClassification::Outside
        };
        assert_eq!(got, want, "at {p:?}");
    }
}

/// `shape` at `center` sketched on `plane`, extruded along its normal
/// by `height`: the operands as a script builds them (seams and edge
/// layout where the sketch puts them).
fn sketched_prism(
    shape: popkern_sketch::ProfileShape<f64>,
    center: [f64; 2],
    plane: popkern_sketch::SketchPlane,
    height: f64,
) -> BRep {
    let normal = plane.normal;
    let sheet = popkern_sketch::AnalyticProfile { shape, center }
        .to_brep(plane, Tolerance::default(), None)
        .expect("sheet");
    let face = sheet.regions[0].face;
    let mut b = sheet.brep;
    b.extrude(face, normal, height, None).expect("extrude");
    b
}

/// A solid cylinder of `radius` and `height` along z, its seam turned
/// `seam_turn` about the axis from +x, then `placed`: one periodic wall
/// face, as a script's exact-wire extrude of a circle builds it.
fn cylinder(radius: f64, height: f64, seam_turn: f64, placed: [Transform3; 2]) -> BRep {
    let tol = Tolerance::default();
    BRep::cylinder(radius, height, tol)
        .and_then(|c| c.transform(&Transform3::rotation_z(seam_turn)))
        .and_then(|c| c.transform(&placed[0]))
        .and_then(|c| c.transform(&placed[1]))
        .expect("cylinder")
}

/// A ring disc (r 49, axis y, y ∈ [-26, -16]) with a web through its
/// rim at the +x extreme, and a boss (r 33, axis x) whose top generator
/// is tangent to the web's top face y = -17 and touches the ring's
/// extreme there. The boss's records with the ring wall run inside the
/// web, through the web's notch in the ring wall, and end on the notch
/// rim at that tangent point. An end on the rim is where any record
/// leaving the face lies, no evidence the record splits it (measured:
/// the six-axis arm's end casting declined in every step of its
/// tangency sweep). The boss stops short of the web's end face, whose
/// coplanar cap is a separate class.
#[test]
fn a_record_through_a_notch_ending_on_its_rim_does_not_decline_the_complex() {
    use std::f64::consts::{FRAC_PI_2, PI};
    let tol = Tolerance::default();
    // Axis -y from y = -16, seam on -x (the sketch plane's u axis).
    let ring = cylinder(
        49.0,
        10.0,
        PI,
        [
            Transform3::rotation_x(FRAC_PI_2),
            Transform3::translation(-817.0, -16.0, 162.5),
        ],
    );
    let web = sketched_prism(
        popkern_sketch::ProfileShape::Rect {
            width: 42.95,
            height: 8.0,
        },
        [-761.475, -21.0],
        popkern_sketch::SketchPlane {
            origin: Point3::new(0.0, 0.0, 137.5),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
        },
        50.0,
    );
    let a = crate::boolean(BoolOperation::Union, ring, web, tol, None).expect("ring and web");
    // Axis +x from x = -770 to -741, seam on -y.
    let boss = cylinder(
        33.0,
        29.0,
        -FRAC_PI_2,
        [
            Transform3::rotation_y(FRAC_PI_2),
            Transform3::translation(-770.0, -50.0, 162.5),
        ],
    );
    let op = BoolOperation::Union;
    let ix = crate::intersect::find_intersections(&a, &boss, &tol).expect("intersections");
    let (body, _) = crate::cellular::boolean::cellular_graph_merge(
        &a,
        &boss,
        &ix,
        op,
        &tol,
        None,
        crate::complex::serve::SelectionSource::Complex,
    )
    .unwrap_or_else(|d| panic!("the complex serves the casting: {}", d.into_error()));
    body.validate().expect("valid");
    let boxes = crate::intersect::compute_face_aabbs(&body);
    for (p, inside) in [
        ((-817.0, -21.0, 162.5), true),
        ((-760.0, -21.0, 160.0), true),
        ((-750.0, -50.0, 162.5), true),
        ((-817.0, -21.0, 220.0), false),
        ((-817.0, -10.0, 162.5), false),
        ((-750.0, -50.0, 200.0), false),
    ] {
        let point = Point3::new(p.0, p.1, p.2);
        let got =
            crate::classify::point_in_solid_final(&point, &body, &boxes, &tol).expect("classifies");
        let want = if inside {
            crate::PointClassification::Inside
        } else {
            crate::PointClassification::Outside
        };
        assert_eq!(got, want, "at {p:?}");
    }
}
