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

use helpers::N_SAMPLES;
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

/// G1 bars for the saddle fan, as the sine of the normal angle
/// (PK-538, PK-498), per seam class.
///
/// Cap-cap (the fan's spokes and sector seams): this instrument's
/// finite-difference normals read 8e-8 constant / 1.4e-7 variable; the
/// exact-jet seam profile reads ~1e-11 on every spoke in both radius
/// modes. Cap-outer (strip and reflex-support rims): 3.1e-6 constant
/// (this instrument's floor) and 1.6e-6 variable. On the exact profile
/// cap-strip reads 1.4e-9 constant and 1.2-1.5e-7 variable. Before
/// PK-576 every class sat at 4e-3..1e-2: the rim reconciles targeted
/// the bridge rim, which lies ON the reflex plane, against the nearest
/// strip. The variable-radius rims then read 1.2e-5 until the PK-584
/// band refit, because the G0 landing's translation-only end bands left
/// the rim tangent tilted off the strip. The ring N-quad this fan
/// replaced measured 1.0.
const CAP_CAP_G1: f64 = 1.0e-6;
/// Cap-outer bar, both radius modes.
const CAP_OUTER_G1: f64 = 1.0e-5;
/// Rim landing bar: worst sine of a rim's own tangent against its
/// neighbour's tangent plane (see the assertion for the measurements).
const RIM_TANGENT_SINE: f64 = 2.5e-7;

/// Serve evidence: the setback fan (4 sectors x 4 pieces), a clean
/// `validate_all`, cap-cap seams under [`CAP_CAP_G1`] and cap-boundary
/// (strip and support) seams under [`CAP_OUTER_G1`].
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
    // PK-584 rim landing: every rim edge recorded, its own tangent in
    // the neighbour's tangent plane. Measured worst sine 1.0e-9
    // (constant) and 6.6e-8 (variable, the degree-9 boundary row's end
    // layer against the variable strip); before the band refit the
    // variable rims tilted 1.2e-5.
    let rims = &result.repair_report.setback_rims;
    assert_eq!(rims.len(), 16, "{name}: one landing record per rim edge");
    let tilt = rims
        .iter()
        .map(|r| r.worst_tangent_sine)
        .fold(0.0, f64::max);
    assert!(
        tilt < RIM_TANGENT_SINE,
        "{name}: rim tangent tilt {tilt:.3e}"
    );
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
        outer.worst_g1_normal < CAP_OUTER_G1,
        "{name}: cap-outer g1 {:.3e}",
        outer.worst_g1_normal
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

/// `(body, vertex, G2 result)` for the saddle at (1, 0.6, 0.6) with
/// radii 0.06/0.08/0.10 under a G2 request and corner-cap preference
/// `corner_cap`.
fn g2_saddle(corner_cap: popkern_blend::CornerCapPreference) -> (BRep, VertexHandle, BlendResult) {
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
        corner_cap,
        ..FilletOptions::default()
    };
    let result = fillet_edges(&body, &specs, &options, Tolerance::default(), None)
        .expect("the saddle serves under a G2 request");
    (body, vertex, result)
}

/// The KP20 fan cannot meet a G2 request at a saddle. At each bridge
/// foot the direction across the bridge is the strip's section
/// direction, so the plane asks for curvature 0 and the strip for 1/r
/// at the same point. With the fit ruled out (strict setback
/// preference) the fan serves as a G1 downgrade: every corner patch
/// reports achieved G1, and the downgrade is a typed warning, not a G2
/// claim.
#[test]
fn a_g2_request_on_a_saddle_fan_serves_an_honest_g1_downgrade() {
    let (_, vertex, result) = g2_saddle(popkern_blend::CornerCapPreference {
        construction: popkern_blend::CornerConstruction::Setback,
        strict: true,
    });
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

/// PK-596: under the default ladder the fan's G2 failure hands the
/// saddle to the curvature-weighted fit, whose cap need not honour the
/// bridge foot's pointwise 0-vs-1/r conflict the fan's layout
/// inherits. Measured 2026-09-30: the fit serves and the whole-blend
/// evidence certifies G2 on every seam (g1 3.3e-6, g2 1.6e-3).
#[test]
fn a_g2_request_on_a_saddle_serves_the_certified_g2_fit() {
    use popkern_blend::blend_continuity_evidence::{
        ContinuityGates, measure_blend_continuity_evidence,
    };
    let (body, vertex, result) = g2_saddle(Default::default());
    let corner: Vec<_> = result
        .corner_patches
        .iter()
        .filter(|p| p.vertex == vertex)
        .collect();
    assert_eq!(corner.len(), 4, "the fitted cap serves: {corner:?}");
    assert!(
        corner.iter().all(|p| p.setback_cap
            && p.method == popkern_blend::CornerMethod::PlateMin
            && p.achieved_continuity == popkern_blend::ContinuityLevel::G2),
        "{corner:?}"
    );
    let pre: std::collections::BTreeSet<_> = body.faces.keys().collect();
    let blend: Vec<_> = result
        .brep
        .faces
        .keys()
        .filter(|f| !pre.contains(f))
        .collect();
    let evidence = measure_blend_continuity_evidence(
        &result.brep,
        &blend,
        popkern_blend::ContinuityLevel::G2,
        &ContinuityGates::commercial(),
        33,
    );
    assert_eq!(
        evidence.certified_achieved,
        Some(popkern_blend::ContinuityLevel::G2),
        "{evidence:#?}"
    );
}
