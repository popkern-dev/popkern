//! Validation of the blend continuity evidence instrument.
//!
//! Same discipline as `pcurve_fidelity_tests`: an instrument that has
//! only run on the thing it investigates is not evidence. These pin a
//! known-sharp reading, a known-G1 one (a real fillet), and the
//! honesty rules for unmeasured seams.

use super::*;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

fn cube() -> BRep {
    BRep::unit_cube(Tolerance::default()).expect("unit cube fixture")
}

/// A single cube fillet and its new-face delta.
fn filleted_cube(radius: f64) -> (BRep, Vec<FaceHandle>) {
    let cube = cube();
    let pre: std::collections::BTreeSet<FaceHandle> = cube.faces.keys().collect();
    let edge = cube.edges.keys().next().expect("cube has edges");
    let result = crate::fillet_edges(
        &cube,
        &[crate::FilletSpec {
            edge,
            radius: crate::FilletRadius::Constant(radius),
        }],
        &crate::FilletOptions::default(),
        Tolerance::default(),
        None,
    )
    .expect("single cube-edge fillet is the most-tested path in the kernel");
    let new_faces: Vec<FaceHandle> = result
        .brep
        .faces
        .keys()
        .filter(|f| !pre.contains(f))
        .collect();
    (result.brep, new_faces)
}

/// KNOWN SHARP: a bare cube face against its neighbours is 90°
/// creases on every side. The instrument must classify every seam
/// intended-sharp, claim no smooth continuity, and still measure the
/// positional closure.
#[test]
fn a_bare_cube_face_reads_all_sharp_and_claims_nothing() {
    let brep = cube();
    let face = brep.faces.keys().next().expect("cube face");
    let evidence = measure_blend_continuity_evidence(
        &brep,
        &[face],
        ContinuityLevel::G1,
        &ContinuityGates::commercial(),
        17,
    );
    assert_eq!(evidence.seams.len(), 4, "a cube face has four seams");
    assert_eq!(
        evidence.sharp_seams, 4,
        "every 90° crease must classify intended-sharp"
    );
    assert_eq!(
        evidence.achieved,
        ContinuityLevel::G0,
        "no gated seam means no smooth claim"
    );
    assert!(evidence.downgraded(), "G1 was requested, G0 measured");
    assert!(
        evidence.worst_g0_gap < 1e-9,
        "cube seams are positionally exact, got {:.3e}",
        evidence.worst_g0_gap
    );
}

/// KNOWN G1-NOT-G2: a rolling-ball fillet on a cube edge. The rails
/// (cylinder tangent to the two support planes) must gate as smooth
/// and pass G1; the curvature jump across a rail is 1/r by
/// construction, so G2 must NOT be claimed — and a G2 request must
/// therefore read as downgraded.
#[test]
fn a_cube_edge_fillet_reads_g1_and_honestly_refuses_g2() {
    let (brep, new_faces) = filleted_cube(0.2);
    assert!(!new_faces.is_empty(), "fillet must create faces");
    let gates = ContinuityGates::commercial();
    let evidence =
        measure_blend_continuity_evidence(&brep, &new_faces, ContinuityLevel::G1, &gates, 17);
    assert!(evidence.measured(), "a committed fillet has live seams");
    assert_eq!(
        evidence.unmeasured_seams, 0,
        "every seam of a served fillet must be measurable"
    );
    assert!(
        evidence.seams.iter().any(|s| s.gated),
        "the tangent rails must classify intended-smooth"
    );
    assert_eq!(
        evidence.achieved,
        ContinuityLevel::G1,
        "rolling-ball delivers G1 (worst gated g1 {:.3e}, g2 {:.3e})",
        evidence.worst_g1_normal,
        evidence.worst_g2_residual
    );
    assert!(!evidence.downgraded(), "G1 requested, G1 measured");
    // The curvature jump across a plane↔cylinder rail is exactly 1/r.
    assert!(
        (evidence.worst_g2_residual - 1.0 / 0.2).abs() < 0.5,
        "rail curvature jump should be ~1/r = 5, got {:.3}",
        evidence.worst_g2_residual
    );

    // Certified: every rail's G1 holds over the whole edge, and the
    // certified curvature jump brackets 1/r.
    assert_eq!(
        evidence.certified_achieved,
        Some(ContinuityLevel::G1),
        "{evidence:?}"
    );
    let rail_jumps: Vec<(f64, f64)> = evidence
        .seams
        .iter()
        .filter_map(|s| s.certified())
        .filter(|c| c.lower[2] > 1.0)
        .map(|c| (c.lower[2], c.g2_upper))
        .collect();
    assert!(!rail_jumps.is_empty(), "{evidence:?}");
    for (lower, upper) in rail_jumps {
        assert!(
            lower <= 5.0 && upper >= 5.0 && upper - lower < 0.1,
            "certified rail jump [{lower}, {upper}] must bracket 1/r = 5"
        );
    }
    let g2_request =
        measure_blend_continuity_evidence(&brep, &new_faces, ContinuityLevel::G2, &gates, 17);
    assert!(
        g2_request.downgraded(),
        "a G2 request against G1 geometry must read as a downgrade"
    );
}

/// HONESTY: strip one rail pcurve and the same body must stop
/// claiming G1 — an unmeasured seam is never a passing one.
#[test]
fn an_unmeasurable_seam_caps_the_claim_at_g0() {
    let (mut brep, new_faces) = filleted_cube(0.2);
    let set: std::collections::BTreeSet<FaceHandle> = new_faces.iter().copied().collect();
    let he = brep
        .half_edges
        .iter()
        .find_map(|(he, hd)| {
            let face = brep.loops.get(hd.loop_)?.face;
            (set.contains(&face) && hd.pcurve.is_some()).then_some(he)
        })
        .expect("a blend face half-edge with a pcurve");
    brep.half_edges[he].pcurve = None;
    let evidence = measure_blend_continuity_evidence(
        &brep,
        &new_faces,
        ContinuityLevel::G1,
        &ContinuityGates::commercial(),
        17,
    );
    assert!(
        evidence.unmeasured_seams >= 1,
        "the stripped seam must be counted, not skipped"
    );
    assert_eq!(
        evidence.achieved,
        ContinuityLevel::G0,
        "an unmeasured seam caps the claim at G0"
    );
}

/// A positional tear beyond the G0 gate floors the claim at G0 even
/// when the tangency still measures clean.
#[test]
fn a_positional_tear_floors_the_claim() {
    let (brep, new_faces) = filleted_cube(0.2);
    let strict = ContinuityGates {
        g0_gap: 1.0e-16,
        ..ContinuityGates::commercial()
    };
    let evidence =
        measure_blend_continuity_evidence(&brep, &new_faces, ContinuityLevel::G1, &strict, 17);
    // With an impossible G0 gate the classification must fall to G0
    // regardless of the (clean) tangency — proving g0 participates.
    assert_eq!(evidence.achieved, ContinuityLevel::G0);
}

/// No faces in, no claim out.
#[test]
fn an_empty_face_set_claims_nothing() {
    let brep = cube();
    let evidence = measure_blend_continuity_evidence(
        &brep,
        &[],
        ContinuityLevel::G2,
        &ContinuityGates::commercial(),
        17,
    );
    assert!(!evidence.measured());
    assert_eq!(evidence.achieved, ContinuityLevel::G0);
}
