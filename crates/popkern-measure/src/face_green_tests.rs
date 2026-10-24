use super::*;
use popkern_core::geometry::Transform3;
use popkern_core::tolerance::Tolerance;
use std::f64::consts::PI;

/// `body`'s curved faces' boundary-integral areas, summed.
fn curved_area(body: &BRep) -> f64 {
    body.faces
        .keys()
        .filter_map(|f| boundary_integral_area(body, f))
        .sum()
}

#[test]
fn primitive_walls_match_their_closed_forms() {
    let tol = Tolerance::default();
    let tilt = Transform3::rotation_x(1.1).compose(&Transform3::rotation_z(0.7));
    let cases = [
        (
            "cylinder",
            BRep::cylinder(0.7, 1.3, tol),
            2.0 * PI * 0.7 * 1.3,
        ),
        ("cone", BRep::cone(1.0, 2.0, tol), PI * 5.0_f64.sqrt()),
        ("sphere", BRep::sphere(0.5, tol), PI),
    ];
    for (name, body, want) in cases {
        let body = body.expect("primitive");
        for (frame, b) in [
            ("home", body.clone()),
            ("tilted", body.transform(&tilt).expect("tilt")),
        ] {
            let got = curved_area(&b);
            assert!(
                (got - want).abs() < 1e-12 * want,
                "{name} {frame}: {got} vs {want}"
            );
        }
    }
}

#[test]
fn a_reversed_face_flag_reads_the_same_area() {
    let tol = Tolerance::default();
    let mut ball = BRep::sphere(0.5, tol).expect("sphere");
    // Turning the body inside out flips every flag and loop together.
    ball.reverse_body_orientation();
    let got = curved_area(&ball);
    assert!((got - PI).abs() < 1e-12 * PI, "{got}");
}

#[test]
fn planes_have_no_boundary_integral_here() {
    let cube = BRep::unit_cube(Tolerance::default()).expect("cube");
    assert!(
        cube.faces
            .keys()
            .all(|f| boundary_integral_area(&cube, f).is_none())
    );
}

/// A torus face (R 1, r 0.3) rewired to the band its two loops
/// bound, each loop one half-edge carrying `pcurves[i]` and no 3D
/// curve.
fn torus_band(pcurves: [Pcurve; 2]) -> (BRep, FaceHandle) {
    let mut body = BRep::torus(1.0, 0.3, Tolerance::default()).expect("torus");
    let face = body.faces.keys().next().expect("a face");
    let template = body.loops[body.faces[face].outer_loop].half_edge;
    let mut loops = Vec::new();
    for pc in pcurves {
        let mut edge = body.edges[body.half_edges[template].edge].clone();
        edge.curve = None;
        let e = body.edges.insert(edge);
        let mut he = body.half_edges[template].clone();
        (he.edge, he.pcurve) = (e, Some(pc));
        let h = body.half_edges.insert(he);
        let l = body
            .loops
            .insert(popkern_core::topology::Loop { half_edge: h, face });
        let he = &mut body.half_edges[h];
        (he.next, he.prev, he.twin, he.loop_) = (h, h, h, l);
        body.edges[e].half_edge = h;
        loops.push(l);
    }
    body.faces[face].outer_loop = loops[0];
    body.faces[face].inner_loops = vec![loops[1]];
    (body, face)
}

/// PK-939: a band across the tube is bounded by two loops winding v,
/// `du = 0` along both; `−∮ G du` alone reads 0 (modulo the torus).
/// The winding term restores `2πrR·Δu`, on any branch of either loop.
/// A band around the axis winds u instead and needs no term.
#[test]
fn torus_bands_measure_their_winding_loops() {
    use std::f64::consts::TAU;
    let (r, big_r) = (0.3, 1.0);
    let across = |a: f64, b: f64, lift: f64| {
        [
            Pcurve::linear((b, lift), (b, lift + TAU), 0.0, 1.0),
            Pcurve::linear((a, TAU - lift), (a, -lift), 0.0, 1.0),
        ]
    };
    for (a, b, lift) in [
        (1.0, 2.5, 0.0),
        (4.0, 7.0, 0.0),
        (1.0, 2.5, TAU),
        (-0.5, 0.8, -2.0),
    ] {
        let (body, face) = torus_band(across(a, b, lift));
        let got = boundary_integral_area(&body, face).expect("area");
        let want = TAU * r * big_r * (b - a);
        assert!(
            (got - want).abs() < 1e-12 * want,
            "u in [{a}, {b}] lift {lift}: {got} vs {want}"
        );
    }
    let (c, d) = (0.4, 2.9);
    let around = [
        Pcurve::linear((0.0, c), (TAU, c), 0.0, 1.0),
        Pcurve::linear((TAU, d), (0.0, d), 0.0, 1.0),
    ];
    let (body, face) = torus_band(around);
    let got = boundary_integral_area(&body, face).expect("area");
    let want = TAU * r * (big_r * (d - c) + r * (d.psin() - c.psin()));
    assert!(
        (got - want).abs() < 1e-12 * want,
        "v in [{c}, {d}]: {got} vs {want}"
    );
}
