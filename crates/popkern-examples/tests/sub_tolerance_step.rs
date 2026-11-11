//! Two blocks united with their tops a step apart (PK-1163, PK-1165).
//!
//! Blocks [0, 10]³ and [x₀, x₀ + 10] × [0, 10] × [0, 10 + h], the
//! second flush against the first (x₀ = 10) or overlapping it (x₀ = 8),
//! at identity and turned, steps h of 0, 0.5, 1, 1.5, 3 and 10
//! tolerances. At the tolerance's resolution a step at or under it is
//! no step: the union is one block, its tops one face, and a step set
//! above rounding is reported (`StepFlattened`). A step over it stays:
//! the tops are two faces with the step face between, every face flat.
//!
//! The vertex table linked two precise candidates within the sum of
//! their radii, each a whole tolerance, so the blocks' corners 1.5e-6
//! apart were one vertex, and the higher top was served tilted down to
//! the lower's edge, 1.5e-6 off its own plane (valid by every check); a
//! record lying along an operand edge read as a crossing at angle zero
//! took a zone of twenty tolerances along it; the assembler's span
//! registry unified edges whose midpoints lay two tolerances apart; and
//! the face merge measured a boundary strictly under the tolerance, so
//! a step of exactly one merged in one layout and not another.

use popkern_bool::sub_tolerance::{SubToleranceKind, take_sub_tolerance_events};
use popkern_bool::{BoolOperation, boolean};
use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::geometry::{Point3, Transform3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

fn block(x0: f64, top: f64, tol: Tolerance, place: &Transform3) -> BRep {
    let rect = |z: f64| {
        vec![
            Point3::new(x0, 0.0, z),
            Point3::new(x0 + 10.0, 0.0, z),
            Point3::new(x0 + 10.0, 10.0, z),
            Point3::new(x0, 10.0, z),
        ]
    };
    BRep::polygon_frustum(&rect(0.0), &rect(top), tol)
        .and_then(|b| b.transform(place))
        .expect("block")
}

/// The farthest any vertex of `body` lies off a face it bounds.
fn worst_vertex_off_face(body: &BRep) -> f64 {
    let mut worst: f64 = 0.0;
    for (fh, face) in body.faces.iter() {
        let Some(surface) = face.surface.as_ref() else {
            continue;
        };
        for p in body.face_outer_loop_positions(fh) {
            if let Some((u, v)) = surface.geometry.closest_point(&p) {
                worst = worst.max(surface.geometry.point_at(u, v).distance_to(&p));
            }
        }
    }
    worst
}

#[test]
fn a_step_at_or_under_the_tolerance_is_flattened_and_one_over_it_stays() {
    let tol = Tolerance::default();
    let linear = tol.linear();
    let turned =
        Transform3::rotation_around_axis(Point3::new(0.0, 0.0, 0.0), Vec3::new(1.0, 0.7, 0.3), 0.6)
            .expect("axis");
    let mut problems = Vec::new();
    for (place_name, place) in [("identity", Transform3::identity()), ("turned", turned)] {
        for (layout, x0) in [("flush", 10.0), ("overlapping", 8.0)] {
            for k in [0.0, 0.5, 1.0, 1.5, 3.0, 10.0] {
                let h = k * linear;
                let name = format!("{place_name} {layout} step {k} tolerances");
                let _ = take_sub_tolerance_events();
                let result = boolean(
                    BoolOperation::Union,
                    block(0.0, 10.0, tol, &place),
                    block(x0, 10.0 + h, tol, &place),
                    tol,
                    None,
                );
                let steps: Vec<f64> = take_sub_tolerance_events()
                    .into_iter()
                    .filter(|e| e.kind == SubToleranceKind::StepFlattened)
                    .map(|e| e.measured_width)
                    .collect();
                let body = match result {
                    Ok(body) => body,
                    Err(e) => {
                        problems.push(format!("{name}: declined {e}"));
                        continue;
                    }
                };
                if let Some(first) = body.validate_all().first() {
                    problems.push(format!("{name}: {first}"));
                }
                let flattened = tol.within_linear(h);
                let faces = if flattened { 6 } else { 8 };
                if body.faces.len() != faces {
                    problems.push(format!("{name}: {} faces, want {faces}", body.faces.len()));
                }
                // A vertex lies within the tolerance of every face it
                // bounds, the flattened step's included.
                let off = worst_vertex_off_face(&body);
                if !tol.within_linear(off) {
                    problems.push(format!("{name}: a vertex {off:.2e} off its face"));
                }
                // The volume, to the tolerance times the area (1200).
                let want = 2000.0 - (10.0 - x0) * 100.0 + 100.0 * h;
                let v = popkern_measure::body_volume(&body, &tol).expect("volume");
                if (v - want).abs() > 1200.0 * linear {
                    problems.push(format!("{name}: volume {v}, want {want}"));
                }
                let reported = h > 0.0 && flattened;
                let ok = match steps.as_slice() {
                    [] => !reported,
                    [w] => reported && (w - h).abs() <= 1e-3 * h,
                    _ => false,
                };
                if !ok {
                    problems.push(format!("{name}: steps reported {steps:?}"));
                }
            }
        }
    }
    assert!(problems.is_empty(), "{problems:#?}");
}
