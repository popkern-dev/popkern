//! PK-704: a partial cylinder face whose arc straddles the cylinder
//! chart's u = 0 seam, cut by a box. The hinge leaf profile — a line,
//! a knuckle arc round the far side, a line back — extruded, its arc in
//! pieces whose boundaries sit at a rotated phase, so a piece crosses
//! the seam; a box tooth bites across the knuckle and the plate.

use popkern_bool::merge_fallback::take_merge_fallback_events;
use popkern_bool::{BoolOperation, boolean};
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_sketch::{AnalyticProfile, PathSegment, ProfileShape, SketchPlane};

fn tol() -> Tolerance {
    Tolerance::default()
}

fn extrude(shape: ProfileShape<f64>, z0: f64, depth: f64) -> BRep {
    let plane = SketchPlane {
        origin: Point3::new(0.0, 0.0, z0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
    };
    let sheet = AnalyticProfile {
        shape,
        center: [0.0, 0.0],
    }
    .to_brep(plane, tol(), None)
    .expect("sheet");
    let face = sheet.regions[0].face;
    let mut b = sheet.brep;
    b.extrude(face, Vec3::new(0.0, 0.0, 1.0), depth, None)
        .expect("extrude");
    b
}

/// The leaf: plate x in [-15, 0] at y in [-4, -2], knuckle radius 4 round
/// the origin; the arc from -90° to its end split into `pieces` equal
/// parts starting `phase` degrees in (so piece boundaries rotate).
fn leaf(pieces: usize, phase: f64) -> BRep {
    let (r, t, w): (f64, f64, f64) = (4.0, 2.0, 15.0);
    let yf = -r + t;
    let xf = -(r * r - yf * yf).sqrt();
    let start = -90.0_f64.to_radians();
    let end = yf.atan2(xf).rem_euclid(std::f64::consts::TAU);
    let mut cuts = vec![start];
    let first = start + phase.to_radians();
    let step = (end - first) / pieces as f64;
    for k in 0..pieces {
        cuts.push(first + step * k as f64);
    }
    cuts.retain(|c| *c >= start);
    cuts.dedup_by(|a, b| (*a - *b).abs() < 1e-9);
    cuts.push(end);
    let point = |a: f64| [r * a.cos(), r * a.sin()];
    let mut segments = vec![PathSegment::Line { to: [0.0, -r] }];
    for (k, w2) in cuts.windows(2).enumerate() {
        let to = if k + 2 == cuts.len() {
            [xf, yf]
        } else {
            point(w2[1])
        };
        segments.push(PathSegment::Arc {
            via: point(0.5 * (w2[0] + w2[1])),
            to,
        });
    }
    segments.push(PathSegment::Line { to: [-w, yf] });
    extrude(
        ProfileShape::Path {
            start: [-w, -r],
            segments,
        },
        0.0,
        60.0,
    )
}

fn tooth() -> BRep {
    extrude(
        ProfileShape::Polygon {
            points: vec![[-4.5, -6.0], [6.0, -6.0], [6.0, 6.0], [-4.5, 6.0]],
        },
        10.0,
        10.0,
    )
}

#[test]
fn seam_straddling_partial_cylinders_cut_by_a_box_serve() {
    for (pieces, phase) in [(4, 0.0), (4, 30.0), (3, 45.0), (5, 70.0), (4, 85.0)] {
        let _ = take_merge_fallback_events();
        let body = leaf(pieces, phase);
        let cut = boolean(BoolOperation::Subtraction, body, tooth(), tol(), None)
            .unwrap_or_else(|e| panic!("{pieces} pieces at phase {phase}: {e}"));
        cut.validate().expect("valid");
        let events = take_merge_fallback_events();
        assert!(events.is_empty(), "{pieces} at {phase}: {events:?}");
    }
}
