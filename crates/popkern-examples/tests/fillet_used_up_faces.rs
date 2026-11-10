//! Fillets that use a support face up exactly: the tangency
//! foundations' blend consumer (design §8.4, F2).
//!
//! A plate rounded fully (both front edges at half its thickness), a
//! slot's floor rounded fully, and a plate filleted at its thickness:
//! in the first two the blends meet tangent to each other where the
//! face between them had its width, in the third the rail runs along
//! the face's far edge. The face used up is gone, not kept with no
//! width: the full round and the slot kept their front face and floor
//! as faces of two coincident edges, and the plate filleted at its
//! thickness stopped its rail a thousandth of the edge short, leaving
//! a strip 2e-3 wide (the anchor walk kept to the open edge).

use popkern_blend::{FilletOptions, FilletRadius, FilletSpec, fillet_edges};
use popkern_core::geometry::curve::CurveEval as _;
use popkern_core::geometry::{Point3, Transform3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use std::f64::consts::PI;

fn block(lo: (f64, f64, f64), hi: (f64, f64, f64), tol: Tolerance) -> BRep {
    let rect = |z: f64| {
        vec![
            Point3::new(lo.0, lo.1, z),
            Point3::new(hi.0, lo.1, z),
            Point3::new(hi.0, hi.1, z),
            Point3::new(lo.0, hi.1, z),
        ]
    };
    BRep::polygon_frustum(&rect(lo.2), &rect(hi.2), tol).expect("block")
}

/// Fillet the edges whose midpoints are `mids` (before `place`) at
/// `radius`; the body placed by `place`. `Err` says why it failed.
fn check(
    body: &BRep,
    place: &Transform3,
    mids: &[Point3],
    radius: f64,
    (want_volume, want_faces): (f64, usize),
) -> Result<(), String> {
    let tol = Tolerance::default();
    let body = body.transform(place).map_err(|e| format!("{e:?}"))?;
    let specs = mids
        .iter()
        .map(|m| {
            let at = place.apply_point(m);
            body.edges
                .iter()
                .find(|(_, e)| {
                    e.curve.as_ref().is_some_and(|c| {
                        let mid = c.geometry.point_at(0.5 * (c.t_start + c.t_end));
                        mid.distance_to(&at) < 1e-6
                    })
                })
                .map(|(edge, _)| FilletSpec {
                    edge,
                    radius: FilletRadius::Constant(radius),
                })
                .ok_or(format!("no edge at {m:?}"))
        })
        .collect::<Result<Vec<_>, _>>()?;
    let served = fillet_edges(&body, &specs, &FilletOptions::default(), tol, None)
        .map_err(|e| format!("declined: {e}"))?
        .brep;
    if let Some(first) = served.validate_all().first() {
        return Err(format!("validate_all: {first}"));
    }
    let v = popkern_measure::body_volume(&served, &tol).map_err(|e| format!("{e}"))?;
    if (v - want_volume).abs() > 1e-9 * want_volume {
        return Err(format!("volume {v} for {want_volume}"));
    }
    if served.faces.len() != want_faces {
        return Err(format!("{} faces for {want_faces}", served.faces.len()));
    }
    Ok(())
}

/// (name, body, edge midpoints, radius, (volume, faces)).
type Case<'a> = (&'a str, &'a BRep, &'a [Point3], f64, (f64, usize));

fn placements() -> [(&'static str, Transform3); 2] {
    let turned =
        Transform3::rotation_around_axis(Point3::new(0.0, 0.0, 0.0), Vec3::new(1.0, 0.7, 0.3), 0.6)
            .expect("axis");
    [("identity", Transform3::identity()), ("turned", turned)]
}

#[test]
fn a_fillet_using_a_face_up_removes_it() {
    let tol = Tolerance::default();
    let plate = block((0.0, 0.0, 0.0), (10.0, 10.0, 2.0), tol);
    let slotted = popkern_bool::boolean(
        popkern_bool::BoolOperation::Subtraction,
        block((0.0, 0.0, 0.0), (10.0, 10.0, 5.0), tol),
        block((4.0, -1.0, 3.0), (6.0, 11.0, 6.0), tol),
        tol,
        None,
    )
    .expect("slot");
    // A quarter disc of radius r leaves r² (1 − π/4) per unit length.
    let corner = |r: f64| r * r * (1.0 - PI / 4.0) * 10.0;
    let front = [Point3::new(5.0, 0.0, 0.0), Point3::new(5.0, 0.0, 2.0)];
    let floor = [Point3::new(4.0, 5.0, 3.0), Point3::new(6.0, 5.0, 3.0)];
    let mut problems = Vec::new();
    for (name, place) in placements() {
        // The full round: the front face goes, the two blends (one
        // cylinder) remain beside top, bottom, back and the ends.
        let cases: [Case; 4] = [
            (
                "full round",
                &plate,
                &front,
                1.0,
                (200.0 - 2.0 * corner(1.0), 6),
            ),
            // The rail along the bottom edge: top, bottom, back, the
            // ends and the blend.
            (
                "at its thickness",
                &plate,
                &front[1..],
                2.0,
                (200.0 - corner(2.0), 6),
            ),
            (
                "half a tolerance past its thickness",
                &plate,
                &front[1..],
                2.0 + 0.5e-6,
                (200.0 - corner(2.0 + 0.5e-6), 6),
            ),
            // The slot's floor goes; its two blends (one cylinder)
            // join the walls.
            (
                "slot floor",
                &slotted,
                &floor,
                1.0,
                (460.0 + 2.0 * corner(1.0), 10),
            ),
        ];
        for (case, body, mids, r, want) in cases {
            if let Err(why) = check(body, &place, mids, r, want) {
                problems.push(format!("{case} ({name}): {why}"));
            }
        }
    }
    assert!(problems.is_empty(), "{problems:#?}");
}
