//! PK-527: the saddle setback serves notched-cube saddles G1.
//!
//! At a saddle vertex, one support face has a reflex corner; on the
//! notch that is x = 1, where the two convex strips' rails diverge.
//! The setback border there is hybrid:
//! - strip arcs meet at the convex faces' rail crossings;
//! - a rail-flat bridge on the reflex face continues both rails;
//! - the ring carve lays that border over the pre-carve's mitre
//!   complex.
//!
//! The KP20 fan then serves the corner. The boundary-ring N-quad it
//! replaces met at 90° across its spokes (live g1 1.0). Design:
//! `docs/architecture/saddle-corner-setback.md`.

#[path = "kp_gate_helpers.rs"]
mod helpers;

use helpers::{G1_GATE, N_SAMPLES};
use popkern_blend::types::{FilletOptions, FilletRadius, FilletSpec};
use popkern_blend::{BlendResult, fillet_edges};
use popkern_core::geometry::Point3;
use popkern_core::geometry::transform::Transform3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, VertexHandle};

/// A unit cube with its (+,+,+) corner notched out, and the saddle
/// vertex at `corner`.
fn notched_cube(corner: Point3) -> (BRep, VertexHandle) {
    let tol = Tolerance::default();
    let notch = Transform3::translation(0.6, 0.6, 0.6).compose(&Transform3::scale_uniform(0.6));
    let tool = BRep::unit_cube(tol)
        .expect("cube")
        .transform(&notch)
        .expect("notch transform");
    let body = popkern_bool::subtraction(BRep::unit_cube(tol).expect("cube"), tool, tol, None)
        .expect("notched cube builds");
    let vertex = body
        .vertices
        .iter()
        .find(|(_, v)| v.position.distance_to(&corner) < 1e-6)
        .map(|(h, _)| h)
        .expect("saddle vertex");
    (body, vertex)
}

fn fillet_saddle(corner: Point3, radii: [FilletRadius; 3]) -> (BlendResult, VertexHandle) {
    let (body, vertex) = notched_cube(corner);
    let mut edges: Vec<_> = body
        .half_edges
        .values()
        .filter(|he| he.origin == vertex)
        .map(|he| he.edge)
        .collect();
    edges.sort();
    edges.dedup();
    let specs: Vec<FilletSpec> = edges
        .into_iter()
        .zip(radii)
        .map(|(edge, radius)| FilletSpec { edge, radius })
        .collect();
    let result = fillet_edges(
        &body,
        &specs,
        &FilletOptions::default(),
        Tolerance::default(),
        None,
    )
    .expect("saddle fillet serves");
    (result, vertex)
}

/// Cap-cap bar for the saddle fan (PK-538): a fifth of the commercial
/// G1 gate. The rail-flat bridge measures 8e-8 (constant radius) and
/// 4.2e-3 (variable radius); the ring N-quad measured 1.0.
const CAP_CAP_G1: f64 = 1.0e-2;

/// Serve evidence: the setback fan (4 sectors x 4 pieces), a clean
/// `validate_all`, cap-cap seams under [`CAP_CAP_G1`], and cap-strip
/// and cap-support seams under the commercial G1 gate.
fn assert_saddle_serves(name: &str, result: &BlendResult, vertex: VertexHandle) {
    let setback = result
        .corner_patches
        .iter()
        .filter(|p| p.vertex == vertex && p.setback_cap)
        .count();
    assert_eq!(setback, 16, "{name}: the saddle serves the setback fan");
    let errors = result.brep.validate_all();
    assert!(errors.is_empty(), "{name}: validate_all: {errors:?}");
    let caps = helpers::cap_faces_from_result(result);
    let seams = helpers::measure_live_cap_seams(&result.brep, &caps, N_SAMPLES);
    assert!(seams.seam_count > 0, "{name}: no cap seams measured");
    assert!(
        seams.worst_g0_gap < 1e-9,
        "{name}: cap-cap g0 {:.3e}",
        seams.worst_g0_gap
    );
    assert!(
        seams.worst_g1_normal < CAP_CAP_G1,
        "{name}: cap-cap g1 {:.3e}",
        seams.worst_g1_normal
    );
    let outer = helpers::measure_cap_outer_boundary(&result.brep, &caps, N_SAMPLES);
    assert!(outer.edges_measured > 0, "{name}: no cap boundary measured");
    assert!(
        outer.worst_g0_gap < 1e-6,
        "{name}: cap-outer g0 {:.3e}",
        outer.worst_g0_gap
    );
    assert!(
        outer.worst_g1_coplanarity < G1_GATE,
        "{name}: cap-outer g1 {:.3e}",
        outer.worst_g1_coplanarity
    );
}

#[test]
fn constant_radius_saddle_serves_the_setback_fan_g1() {
    let (result, vertex) = fillet_saddle(
        Point3::new(1.0, 0.6, 0.6),
        [0.06, 0.08, 0.10].map(FilletRadius::Constant),
    );
    assert_saddle_serves("constant", &result, vertex);
}

#[test]
fn variable_radius_saddle_serves_the_setback_fan_g1() {
    let (result, vertex) = fillet_saddle(
        Point3::new(1.0, 0.6, 0.6),
        [(0.05, 0.09), (0.07, 0.11), (0.08, 0.06)]
            .map(|(start, end)| FilletRadius::Variable { start, end }),
    );
    assert_saddle_serves("variable", &result, vertex);
}

/// The notch's other two saddle vertices serve the same way: the
/// reflex face and the rail frames are read from geometry, not from a
/// fixed axis.
#[test]
fn every_notch_saddle_serves_the_setback_fan_g1() {
    for corner in [Point3::new(0.6, 1.0, 0.6), Point3::new(0.6, 0.6, 1.0)] {
        let (result, vertex) =
            fillet_saddle(corner, [0.06, 0.08, 0.10].map(FilletRadius::Constant));
        assert_saddle_serves(&format!("{corner:?}"), &result, vertex);
    }
}

/// A G2 request cannot be met at a saddle. At each bridge foot the
/// direction across the bridge is the strip's section direction, so
/// the plane asks for curvature 0 and the strip for 1/r at the same
/// point. The fan's seam-G2 postcondition declines, and the rung's G1
/// retry serves. Every corner patch reports achieved G1, and the
/// downgrade is a typed warning, not a G2 claim.
#[test]
fn a_g2_request_on_a_saddle_serves_an_honest_g1_downgrade() {
    let (body, vertex) = notched_cube(Point3::new(1.0, 0.6, 0.6));
    let mut edges: Vec<_> = body
        .half_edges
        .values()
        .filter(|he| he.origin == vertex)
        .map(|he| he.edge)
        .collect();
    edges.sort();
    edges.dedup();
    let specs: Vec<FilletSpec> = edges
        .into_iter()
        .zip([0.06, 0.08, 0.10])
        .map(|(edge, r)| FilletSpec {
            edge,
            radius: FilletRadius::Constant(r),
        })
        .collect();
    let options = FilletOptions {
        continuity: popkern_blend::ContinuityLevel::G2,
        ..FilletOptions::default()
    };
    let result = fillet_edges(&body, &specs, &options, Tolerance::default(), None)
        .expect("the saddle serves under a G2 request");
    let corner: Vec<_> = result
        .corner_patches
        .iter()
        .filter(|p| p.vertex == vertex)
        .collect();
    assert_eq!(corner.len(), 16, "the setback fan serves");
    assert!(
        corner
            .iter()
            .all(|p| p.achieved_continuity == popkern_blend::ContinuityLevel::G1),
        "the saddle cap reports achieved G1, never G2"
    );
    assert!(
        result.edge_warnings.iter().any(|(_, w)| matches!(
            w,
            popkern_blend::error::BlendWarning::CornerContinuityDowngraded { vertex: v, .. } if *v == vertex
        )),
        "the downgrade is reported: {:?}",
        result.edge_warnings
    );
}
