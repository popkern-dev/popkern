//! PK-142 stage 2 — TrimBoth face–face blend acceptance on CLOSED
//! (blind) pockets, where the strip end arcs land on real end-wall
//! faces. Two fixtures: a rectangular pocket (parallel walls — the
//! line-contact closed-form march) and a trapezoid pocket (tilted
//! walls — the general contact-solver march with an analytically
//! pinned tangency height).

use popkern_blend::{FaceBlendAttach, FaceBlendSpec, face_blend};
use popkern_bool::subtraction;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};

fn tol() -> Tolerance {
    Tolerance::default()
}

fn block() -> BRep {
    BRep::prism(
        &[
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(10.0, 0.0, 0.0),
            Point3::new(10.0, 8.0, 0.0),
            Point3::new(0.0, 8.0, 0.0),
        ],
        Vec3::new(0.0, 0.0, 1.0),
        6.0,
        tol(),
    )
    .expect("block builds")
}

/// Blind rectangular pocket: 4 wide (x ∈ [3,7]) × 6 long (y ∈ [1,7])
/// sunk from the top to a floor at z = 1.
fn rect_pocket_body() -> BRep {
    let tool = BRep::prism(
        &[
            Point3::new(3.0, 1.0, 1.0),
            Point3::new(7.0, 1.0, 1.0),
            Point3::new(7.0, 7.0, 1.0),
            Point3::new(3.0, 7.0, 1.0),
        ],
        Vec3::new(0.0, 0.0, 1.0),
        6.0,
        tol(),
    )
    .expect("pocket tool builds");
    subtraction(block(), tool, tol(), None).expect("pocket subtract")
}

/// Blind trapezoid pocket: walls tilted 1:5 in x-z, converging from
/// a 4-wide mouth to a 2-wide floor at z = 2, running y ∈ [1,7].
fn trapezoid_pocket_body() -> BRep {
    let tool = BRep::prism(
        &[
            Point3::new(3.0, 1.0, 7.0),
            Point3::new(4.0, 1.0, 2.0),
            Point3::new(6.0, 1.0, 2.0),
            Point3::new(7.0, 1.0, 7.0),
        ],
        Vec3::new(0.0, 1.0, 0.0),
        6.0,
        tol(),
    )
    .expect("trapezoid tool builds");
    subtraction(block(), tool, tol(), None).expect("trapezoid subtract")
}

/// Pocket walls by plane signature: dominant ±x normal, plus the
/// required tilt component for the trapezoid case.
fn pocket_walls(brep: &BRep, min_tilt: f64, max_tilt: f64) -> (FaceHandle, FaceHandle) {
    let mut left = None;
    let mut right = None;
    for (fh, f) in brep.faces.iter() {
        let Some(sof) = f.surface.as_ref() else {
            continue;
        };
        let SurfaceGeometry::Plane { origin, normal, .. } = &sof.geometry else {
            continue;
        };
        let tilt = normal.z.abs();
        if normal.y.abs() > 0.01 || normal.x.abs() < 0.5 || tilt < min_tilt || tilt > max_tilt {
            continue;
        }
        // Exclude the block's outer x faces (x = 0 / x = 10).
        if origin.x < 2.0 || origin.x > 8.0 {
            continue;
        }
        if normal.x > 0.0 {
            left = Some(fh);
        } else {
            right = Some(fh);
        }
    }
    (
        left.expect("left pocket wall"),
        right.expect("right pocket wall"),
    )
}

fn strip_loop_edge_count(body: &BRep, strip: FaceHandle) -> usize {
    let start = body.loops[body.faces[strip].outer_loop].half_edge;
    let mut he = start;
    let mut count = 0;
    loop {
        count += 1;
        he = body.half_edges[he].next;
        if he == start || count > 32 {
            break;
        }
    }
    count
}

#[test]
fn trim_both_rect_pocket_attaches_watertight() {
    let body = rect_pocket_body();
    let faces_before = body.faces.len();
    let (wall_a, wall_b) = pocket_walls(&body, 0.0, 0.05);
    let result = face_blend(
        &body,
        &FaceBlendSpec {
            face_a: wall_a,
            face_b: wall_b,
            radius: 2.0,
            radius_end: None,
            attach: FaceBlendAttach::TrimBoth,
        },
        tol(),
    )
    .expect("rect-pocket TrimBoth must complete");
    let attached = result.attached.as_ref().expect("TrimBoth returns a body");
    let strip = result.strip_face.expect("strip face named");

    attached.validate().expect("attached body must validate");
    assert_eq!(
        attached.faces.len(),
        faces_before,
        "4 splits then 5 doomed faces merged into the strip — net face count preserved"
    );
    assert_eq!(
        strip_loop_edge_count(attached, strip),
        4,
        "strip boundary = 2 contact lines + 2 end arcs"
    );
    // Canonical line contact at the wall mid-extent height z = 3.5.
    for p in result.contacts_a.iter().chain(result.contacts_b.iter()) {
        assert!(
            (p.z - 3.5).abs() < 1e-6,
            "parallel line contact must sit at the wall mid-extent: {p:?}"
        );
    }
}

#[test]
fn trim_both_trapezoid_pocket_attaches_watertight() {
    let body = trapezoid_pocket_body();
    let faces_before = body.faces.len();
    let (wall_a, wall_b) = pocket_walls(&body, 0.1, 0.4);
    let radius = 1.3;
    let result = face_blend(
        &body,
        &FaceBlendSpec {
            face_a: wall_a,
            face_b: wall_b,
            radius,
            radius_end: None,
            attach: FaceBlendAttach::TrimBoth,
        },
        tol(),
    )
    .expect("trapezoid-pocket TrimBoth must complete");
    let attached = result.attached.as_ref().expect("TrimBoth returns a body");
    let strip = result.strip_face.expect("strip face named");

    attached.validate().expect("attached body must validate");
    assert_eq!(
        attached.faces.len(),
        faces_before,
        "net face count preserved"
    );
    assert_eq!(strip_loop_edge_count(attached, strip), 4);

    // Analytic tangency: walls tilt 1:5, so the ball centre sits at
    // z = r·sqrt(26) − 3 and contacts at z_c − r/sqrt(26).
    let s26 = 26.0_f64.sqrt();
    let z_expect = radius * s26 - 3.0 - radius / s26;
    for p in result.contacts_a.iter().chain(result.contacts_b.iter()) {
        assert!(
            (p.z - z_expect).abs() < 1e-6,
            "contact height off the analytic tangency: {p:?} vs z={z_expect}"
        );
    }
}
