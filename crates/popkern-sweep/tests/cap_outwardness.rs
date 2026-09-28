//! PK-565: every planar cap of a swept or lofted solid faces OUT of
//! the solid, and every self-loop cap pcurve starts at its loop
//! vertex and runs in its half-edge's direction.
//!
//! The winding-honesty audit cannot see the defect this pins: a
//! full-circle start cap used to take the END cap's normal with its
//! flag and UV winding inverted together — honest winding, inward
//! face. Outwardness is measured directly: both caps of a two-cap
//! body must fold their normals away from the midpoint of the caps.

use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::{CurveEval, NurbsCurve, Point3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_sweep::{LoftOptions, LoftSection, SweepOptions, SweepProfile};

fn circle(r: f64, z: f64) -> NurbsCurve {
    let w = std::f64::consts::FRAC_1_SQRT_2;
    let unit = [
        (1.0, 0.0),
        (1.0, 1.0),
        (0.0, 1.0),
        (-1.0, 1.0),
        (-1.0, 0.0),
        (-1.0, -1.0),
        (0.0, -1.0),
        (1.0, -1.0),
        (1.0, 0.0),
    ];
    NurbsCurve {
        degree: 2,
        knots: vec![
            0.0, 0.0, 0.0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1.0, 1.0, 1.0,
        ],
        control_points: unit
            .iter()
            .map(|&(x, y)| Point3::new(r * x, r * y, z))
            .collect(),
        weights: vec![1.0, w, 1.0, w, 1.0, w, 1.0, w, 1.0],
    }
}

fn line(a: Point3, b: Point3) -> NurbsCurve {
    NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![a, b],
        weights: vec![1.0, 1.0],
    }
}

fn sweep(profile: SweepProfile, end: Point3) -> BRep {
    popkern_sweep::sweep(
        &profile,
        &line(Point3::new(0.0, 0.0, 0.0), end),
        &SweepOptions::default(),
        Tolerance::default(),
        None,
    )
    .expect("sweep")
    .brep
}

/// Planar caps as (plane origin, flag-folded unit normal).
fn planar_caps(brep: &BRep) -> Vec<(Point3, popkern_core::geometry::Vec3)> {
    brep.faces
        .iter()
        .filter_map(|(_, f)| {
            let sof = f.surface.as_ref()?;
            let SurfaceGeometry::Plane { origin, .. } = &sof.geometry else {
                return None;
            };
            let n = sof.geometry.normal_at(0.0, 0.0);
            Some((*origin, if sof.orientation { n } else { n * -1.0 }))
        })
        .collect()
}

fn assert_caps_outward(brep: &BRep, fixture: &str) {
    let caps = planar_caps(brep);
    assert_eq!(caps.len(), 2, "{fixture}: expected two planar caps");
    let mid = Point3::new(
        0.5 * (caps[0].0.x + caps[1].0.x),
        0.5 * (caps[0].0.y + caps[1].0.y),
        0.5 * (caps[0].0.z + caps[1].0.z),
    );
    for (i, (origin, n)) in caps.iter().enumerate() {
        let out = *origin - mid;
        let cos = n.dot(&out) / out.length().max(1e-30);
        assert!(
            cos > 0.5,
            "{fixture}: cap {i} at {origin:?} folds {n:?} INTO the solid (cos {cos:+.3})"
        );
    }
}

/// Every self-loop half-edge on a planar cap: the pcurve starts and
/// ends at the loop vertex and leaves it in the half-edge's 3-D
/// direction (the canonical half runs the edge curve forward).
fn assert_self_loop_pcurves_anchored(brep: &BRep, fixture: &str) {
    let mut checked = 0;
    for (he, h) in brep.half_edges.iter() {
        let twin = &brep.half_edges[h.twin];
        if h.origin != twin.origin {
            continue;
        }
        let face = &brep.faces[brep.loops[h.loop_].face];
        let Some(sof) = face.surface.as_ref() else {
            continue;
        };
        if !matches!(sof.geometry, SurfaceGeometry::Plane { .. }) {
            continue;
        }
        let pcurve = h.pcurve.as_ref().expect("self-loop cap pcurve");
        let (s0, s1) = pcurve.parameter_range();
        let lift = |s: f64| {
            let uv = pcurve.point_at(s);
            sof.geometry.point_at(uv.x, uv.y)
        };
        let vertex = brep.vertices[h.origin].position;
        assert!(
            lift(s0).distance_to(&vertex) < 1e-9 && lift(s1).distance_to(&vertex) < 1e-9,
            "{fixture}: self-loop pcurve of {he:?} does not start/end at its vertex"
        );
        let edge = &brep.edges[h.edge];
        let curve = edge.curve.as_ref().expect("edge curve");
        let span = curve.t_end - curve.t_start;
        let ahead_fwd = curve.geometry.point_at(curve.t_start + 0.05 * span);
        let ahead_rev = curve.geometry.point_at(curve.t_end - 0.05 * span);
        let step = lift(s0 + 0.05 * (s1 - s0));
        let runs_forward = step.distance_to(&ahead_fwd) < step.distance_to(&ahead_rev);
        assert_eq!(
            runs_forward,
            edge.half_edge == he,
            "{fixture}: self-loop pcurve of {he:?} runs against its half-edge"
        );
        checked += 1;
    }
    assert!(checked > 0, "{fixture}: no self-loop cap half-edges found");
}

fn assert_cap_contract(brep: &BRep, fixture: &str) {
    assert_caps_outward(brep, fixture);
    assert_self_loop_pcurves_anchored(brep, fixture);
    let dishonest: Vec<_> = brep
        .audit_winding_honesty()
        .into_iter()
        .filter(|v| v.honest == Some(false))
        .collect();
    assert!(dishonest.is_empty(), "{fixture}: dishonest {dishonest:?}");
}

#[test]
fn circle_extrude_up_caps_face_out() {
    let b = sweep(
        SweepProfile::Wire(circle(1.0, 0.0)),
        Point3::new(0.0, 0.0, 3.0),
    );
    assert_cap_contract(&b, "extrude +z");
}

#[test]
fn circle_extrude_tilted_caps_face_out() {
    let b = sweep(
        SweepProfile::Wire(circle(0.8, 0.0)),
        Point3::new(1.0, 0.5, 3.0),
    );
    assert_cap_contract(&b, "extrude tilted");
}

#[test]
fn circle_extrude_down_caps_face_out() {
    let b = sweep(
        SweepProfile::Wire(circle(0.5, 0.0)),
        Point3::new(0.0, 0.0, -2.0),
    );
    assert_cap_contract(&b, "extrude -z");
}

#[test]
fn annulus_extrude_caps_face_out() {
    let b = sweep(
        SweepProfile::MultiWire(vec![circle(1.0, 0.0), circle(0.4, 0.0)]),
        Point3::new(0.0, 0.0, 2.0),
    );
    assert_cap_contract(&b, "annulus extrude");
}

#[test]
fn circular_loft_caps_face_out() {
    let sections = vec![
        LoftSection {
            curve: circle(1.0, 0.0),
            position: None,
        },
        LoftSection {
            curve: circle(0.6, 2.0),
            position: None,
        },
    ];
    let b = popkern_sweep::loft(
        &sections,
        &LoftOptions::default(),
        Tolerance::default(),
        None,
    )
    .expect("loft")
    .brep;
    assert_cap_contract(&b, "circular loft");
}
