//! Two blocks set a gap apart (PK-1131): a gap at or under the linear
//! tolerance is closed, one over it stays open.
//!
//! Blocks [0, 10]³ and [10 + g, 20 + g] × [y, 10 + y] × [0, 10], flush
//! (y = 0) or offset (y = 3), at identity and turned. At the
//! tolerance's resolution a gap at or under it is no gap: the union is
//! one block across it (its volume larger by at most the facing area
//! times the gap), reported as a closed slot; the cut leaves the first
//! block and the common part is empty. A gap over it is two bodies.
//! The union declined four ways in turn: the plane × plane lines a
//! tolerance off a face's edges were clipped away (an absolute 1e-9
//! for a line parallel to an edge), so no record joined the blocks and
//! they went side by side touching; the contact search read a
//! coplanar neighbour ten tolerances past a trim as touching it (the
//! trim test's boundary band), so a gap of 1e-5 declined as
//! non-manifold, turned; a gap of exactly the tolerance on the turned
//! blocks measured 1.8e-9 over it, past the rounding allowance; and
//! the validator's pair scan, run at twice the tolerance, read faces
//! 1.5e-6 apart as one lying on the other.

use popkern_bool::sub_tolerance::{SubToleranceKind, take_sub_tolerance_events};
use popkern_bool::{BoolOperation, boolean};
use popkern_core::geometry::{Point3, Transform3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

fn block(lo: (f64, f64, f64), size: f64, tol: Tolerance, place: &Transform3) -> BRep {
    let (x0, y0, z0) = lo;
    let rect = |z: f64| {
        vec![
            Point3::new(x0, y0, z),
            Point3::new(x0 + size, y0, z),
            Point3::new(x0 + size, y0 + size, z),
            Point3::new(x0, y0 + size, z),
        ]
    };
    BRep::polygon_frustum(&rect(z0), &rect(z0 + size), tol)
        .and_then(|b| b.transform(place))
        .expect("block")
}

/// What a served body must be.
struct Want {
    faces: usize,
    volume: f64,
    /// How far the volume may stray (the closed gap's slab).
    slack: f64,
}

fn check(
    name: &str,
    result: Result<BRep, popkern_bool::BoolError>,
    want: &Want,
    tol: &Tolerance,
) -> Option<String> {
    let body = match result {
        Ok(body) => body,
        Err(e) => return Some(format!("{name}: declined {e}")),
    };
    if body.faces.is_empty() && want.faces == 0 {
        return None;
    }
    if let Some(first) = body.validate_all().first() {
        return Some(format!("{name}: {first}"));
    }
    let v = popkern_measure::body_volume(&body, tol).expect("volume");
    if body.faces.len() != want.faces || (v - want.volume).abs() > want.slack + 1e-9 {
        return Some(format!(
            "{name}: {} faces, volume {v} (want {} faces, {} within {})",
            body.faces.len(),
            want.faces,
            want.volume,
            want.slack
        ));
    }
    None
}

#[test]
fn a_gap_at_or_under_the_tolerance_is_closed_and_one_over_it_stays_open() {
    let tol = Tolerance::default();
    let linear = tol.linear();
    let turned =
        Transform3::rotation_around_axis(Point3::new(0.0, 0.0, 0.0), Vec3::new(1.0, 0.7, 0.3), 0.6)
            .expect("axis");
    let mut problems = Vec::new();
    for (place_name, place) in [("identity", Transform3::identity()), ("turned", turned)] {
        for (layout, y, fused_faces) in [("flush", 0.0, 6), ("offset", 3.0, 10)] {
            for gap in [-1e-6, 0.0, 5e-7, 1e-6, 1.5e-6, 1e-5] {
                let name = |op: &str| format!("{place_name} {layout} gap {gap:e} {op}");
                let a = || block((0.0, 0.0, 0.0), 10.0, tol, &place);
                let b = || block((10.0 + gap, y, 0.0), 10.0, tol, &place);
                let closed = gap.abs() <= linear;
                // A body known to the tolerance has its volume known
                // to its area times the tolerance (1200 for the two
                // blocks); the slot's slab, the facing area times the
                // gap, lies inside that.
                let slack = 1200.0 * linear;
                let union = if closed {
                    Want {
                        faces: fused_faces,
                        volume: 2000.0,
                        slack,
                    }
                } else {
                    Want {
                        faces: 12,
                        volume: 2000.0,
                        slack: 0.0,
                    }
                };
                let _ = take_sub_tolerance_events();
                problems.extend(check(
                    &name("union"),
                    boolean(BoolOperation::Union, a(), b(), tol, None),
                    &union,
                    &tol,
                ));
                // A slot set at a gap is reported closed; blocks set
                // flush, or apart, or overlapping, have none.
                let slots = take_sub_tolerance_events()
                    .into_iter()
                    .filter(|e| e.kind == SubToleranceKind::SlotClosed)
                    .collect::<Vec<_>>();
                let reported = slots.first().map(|e| e.measured_width);
                let expect_slot = gap > 0.0 && closed;
                match reported {
                    Some(w) if expect_slot && (w - gap).abs() <= 1e-6 * gap => {}
                    None if !expect_slot => {}
                    _ => problems.push(format!("{}: slot reports {slots:?}", name("union"))),
                }
                let first = Want {
                    faces: 6,
                    volume: 1000.0,
                    slack,
                };
                problems.extend(check(
                    &name("cut"),
                    boolean(BoolOperation::Subtraction, a(), b(), tol, None),
                    &first,
                    &tol,
                ));
                let empty = Want {
                    faces: 0,
                    volume: 0.0,
                    slack: 0.0,
                };
                problems.extend(check(
                    &name("common"),
                    boolean(BoolOperation::Intersection, a(), b(), tol, None),
                    &empty,
                    &tol,
                ));
            }
        }
    }
    assert!(problems.is_empty(), "{problems:#?}");
}
