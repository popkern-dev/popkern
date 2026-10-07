//! PK-702: a hinge grown on two existing bodies along their seam — a
//! box and the lid sitting on it. The knuckles join the bodies (one
//! shell each, watertight, the original faces kept), the wall check
//! refuses a wall too thin for the web, and the lid swings about the
//! pin axis until it meets the box exactly where the cross-sections say.

use popkern_core::geometry::transform::Transform3;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_hinge::{
    HingeComponent, HingeError, HingeKind, HingeParams, Seam, seam_hinge_component,
};

fn tol() -> Tolerance {
    Tolerance::default()
}

/// An axis-aligned box: a rectangle sketched at the base, extruded.
fn block(x: [f64; 2], y: [f64; 2], z: [f64; 2]) -> BRep {
    let plane = popkern_sketch::SketchPlane {
        origin: Point3::new(0.0, 0.0, z[0]),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
    };
    let sheet = popkern_sketch::AnalyticProfile {
        shape: popkern_sketch::ProfileShape::Rect {
            width: x[1] - x[0],
            height: y[1] - y[0],
        },
        center: [0.5 * (x[0] + x[1]), 0.5 * (y[0] + y[1])],
    }
    .to_brep(plane, tol(), None)
    .expect("sheet");
    let face = sheet.regions[0].face;
    let mut b = sheet.brep;
    b.extrude(face, Vec3::new(0.0, 0.0, 1.0), z[1] - z[0], None)
        .expect("extrude");
    b
}

fn volume(b: &BRep) -> f64 {
    popkern_measure::body_volume(b, &tol()).expect("volume")
}

fn watertight(b: &BRep) -> bool {
    let mesh =
        popkern_mesh::tessellate(b, &popkern_mesh::TessellationOptions::display()).expect("mesh");
    popkern_mesh::welded_open_edges(&mesh) == 0
}

/// Box 80 × 50 × 30, lid 80 × 50 × 10 on top; the seam is the top back
/// edge (y = 50, z = 30), the hinge 64 long starting at x = 8.
fn setup() -> (BRep, BRep, Seam, HingeParams) {
    let box_ = block([0.0, 80.0], [0.0, 50.0], [0.0, 30.0]);
    let lid = block([0.0, 80.0], [0.0, 50.0], [30.0, 40.0]);
    let seam = Seam {
        origin: Point3::new(8.0, 50.0, 30.0),
        axis: Vec3::new(1.0, 0.0, 0.0),
        outward: Vec3::new(0.0, 1.0, 0.0),
        up: Vec3::new(0.0, 0.0, 1.0),
    };
    let params = HingeParams {
        pin_diameter: 3.0,
        thickness: 1.6,
        ..HingeParams::print_in_place(64.0)
    };
    (box_, lid, seam, params)
}

#[test]
fn knuckles_join_the_box_and_the_lid() {
    let (box_, lid, seam, params) = setup();
    let a = seam_hinge_component(
        HingeKind::PrintInPlace,
        HingeComponent::LeafA,
        &seam,
        &params,
        Some(&box_),
        tol(),
    )
    .expect("box knuckles");
    let b = seam_hinge_component(
        HingeKind::PrintInPlace,
        HingeComponent::LeafB,
        &seam,
        &params,
        Some(&lid),
        tol(),
    )
    .expect("lid knuckles");
    for (name, before, after) in [("box", &box_, &a), ("lid", &lid, &b)] {
        assert_eq!(after.shells.len(), 1, "{name}: one shell");
        assert!(watertight(after), "{name} watertight");
        let grown = volume(after) - volume(before);
        assert!(grown > 0.0, "{name} gained the knuckles: {grown}");
    }
    let bad = seam_hinge_component(
        HingeKind::PrintInPlace,
        HingeComponent::Pin,
        &seam,
        &params,
        None,
        tol(),
    );
    assert!(bad.is_err(), "a print-in-place seam hinge has no pin");
}

#[test]
fn a_thin_wall_is_refused() {
    let (_, lid, seam, params) = setup();
    // A box shelled to a 0.5 wall at the back: the web would poke
    // through into the cavity.
    let outer = block([0.0, 80.0], [0.0, 50.0], [0.0, 30.0]);
    let cavity = block([0.5, 79.5], [0.5, 49.5], [0.5, 31.0]);
    let shelled = popkern_bool::boolean(
        popkern_bool::BoolOperation::Subtraction,
        outer,
        cavity,
        tol(),
        None,
    )
    .expect("shell");
    let err = seam_hinge_component(
        HingeKind::PrintInPlace,
        HingeComponent::LeafA,
        &seam,
        &params,
        Some(&shelled),
        tol(),
    )
    .expect_err("thin wall");
    assert!(
        matches!(err, HingeError::SeamWall { body: "A", .. }),
        "{err}"
    );
    let _ = lid;
}

/// Whether two convex polygons overlap (separating-axis theorem).
fn overlap(a: &[[f64; 2]], b: &[[f64; 2]]) -> bool {
    for poly in [a, b] {
        for i in 0..poly.len() {
            let (p, q) = (poly[i], poly[(i + 1) % poly.len()]);
            let n = [q[1] - p[1], p[0] - q[0]];
            let proj = |s: &[[f64; 2]]| {
                s.iter()
                    .map(|v| v[0] * n[0] + v[1] * n[1])
                    .fold((f64::MAX, f64::MIN), |(lo, hi), x| (lo.min(x), hi.max(x)))
            };
            let ((alo, ahi), (blo, bhi)) = (proj(a), proj(b));
            if ahi <= blo + 1e-9 || bhi <= alo + 1e-9 {
                return false;
            }
        }
    }
    true
}

#[test]
fn the_lid_swings_back_until_it_meets_the_box_where_the_sections_say() {
    let (box_, lid, seam, params) = setup();
    let a = seam_hinge_component(
        HingeKind::PrintInPlace,
        HingeComponent::LeafA,
        &seam,
        &params,
        Some(&box_),
        tol(),
    )
    .expect("box knuckles");
    let b = seam_hinge_component(
        HingeKind::PrintInPlace,
        HingeComponent::LeafB,
        &seam,
        &params,
        Some(&lid),
        tol(),
    )
    .expect("lid knuckles");
    // The pin axis: outward of the seam by the knuckle radius plus the
    // clearance.
    let r = params.pin_diameter / 2.0 + params.clearance + params.thickness;
    let a_off = r + params.clearance;
    let axis_y = 50.0 + a_off;
    // Sections relative to the axis (n along +y, v along +z): the box
    // below the seam, the lid above it; opening turns the lid so its
    // front rises (clockwise in (n, v)).
    let box_rect = [
        [-50.0 - a_off, -30.0],
        [-a_off, -30.0],
        [-a_off, 0.0],
        [-50.0 - a_off, 0.0],
    ];
    let lid_rect = [
        [-50.0 - a_off, 0.0],
        [-a_off, 0.0],
        [-a_off, 10.0],
        [-50.0 - a_off, 10.0],
    ];
    let hits = |deg: f64| {
        let (s, c) = deg.to_radians().sin_cos();
        let turned: Vec<[f64; 2]> = lid_rect
            .iter()
            .map(|p| [p[0] * c + p[1] * s, -p[0] * s + p[1] * c])
            .collect();
        overlap(&box_rect, &turned)
    };
    let mut deg = 1.0;
    while !hits(deg) {
        deg += 0.01;
        assert!(deg < 359.0);
    }
    let expected = deg;
    // The kernel: turn the knuckled lid about the pin axis, intersect
    // with the knuckled box, bisect the first overlap.
    let overlaps = |deg: f64| {
        let turn = Transform3::rotation_around_axis(
            Point3::new(0.0, axis_y, 30.0),
            Vec3::new(1.0, 0.0, 0.0),
            -deg.to_radians(),
        )
        .expect("rotation");
        let moved = b.transform(&turn).expect("move");
        match popkern_bool::boolean(
            popkern_bool::BoolOperation::Intersection,
            a.clone(),
            moved,
            tol(),
            None,
        ) {
            Ok(x) => !x.faces.is_empty() && volume(&x) > 1e-6,
            Err(_) => true,
        }
    };
    let (mut clear, mut hit) = (expected - 3.0, expected + 3.0);
    assert!(!overlaps(clear), "apart 3° before the analytic contact");
    assert!(overlaps(hit), "overlapping 3° past it");
    while hit - clear > 0.1 {
        let mid = 0.5 * (clear + hit);
        if overlaps(mid) {
            hit = mid;
        } else {
            clear = mid;
        }
    }
    eprintln!("seam hinge: lid meets the box at {hit:.2}° (sections: {expected:.2}°)");
    assert!((hit - expected).abs() < 0.5, "{hit} vs {expected}");
}
