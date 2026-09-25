//! PK-498 negative controls for [`measure_blend_continuity_evidence`]:
//! proof that the evidence can see each defect class it claims to
//! classify.
//!
//! The positive fixture is the PK-389 G2 rim fillet, an exact surface
//! of revolution whose section runs in v. It measures `achieved = G2`.
//! Each control perturbs ONE control row next to the v = 0 rail and
//! first self-checks, by direct evaluation, that the lower orders it
//! claims to preserve really are preserved. A control that silently
//! broke G0 would prove nothing about G1.
//! - **Row 2 moves along the normal:** G0 and G1 exact, transverse
//!   curvature jumps. The evidence must read G1. The along-seam normal
//!   curvature is unchanged here (Meusnier: it is fixed by the shared
//!   curve and normal), so only the transverse term can catch it.
//! - **Row 1 moves along the normal, a small tilt:** G0 exact, the
//!   normal tilts. The evidence must read G0.
//! - **Row 1, a tilt past the smooth-intent bound:** the seam
//!   reclassifies as intended-sharp, so G1 is no longer gated there.
//!   The sharp-seam census must record it; that census is what
//!   baselines pin.

use super::{ContinuityGates, measure_blend_continuity_evidence};
use crate::types::{ContinuityLevel, FilletOptions, FilletRadius, FilletSpec};
use crate::{BlendResult, fillet_edges};
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::{NurbsSurface, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};

fn g2_rim() -> BlendResult {
    let tol = Tolerance::default();
    let brep = BRep::cylinder(1.0, 2.0, tol).expect("cylinder");
    let specs: Vec<FilletSpec> = brep
        .edges
        .iter()
        .filter(|(_, d)| {
            let he = &brep.half_edges[d.half_edge];
            let a = brep.vertices[he.origin].position;
            let b = brep.vertices[brep.half_edges[he.twin].origin].position;
            brep.loops[he.loop_].face != brep.loops[brep.half_edges[he.twin].loop_].face
                && (a.z - 2.0).abs() < 1e-9
                && (b.z - 2.0).abs() < 1e-9
        })
        .map(|(edge, _)| FilletSpec {
            edge,
            radius: FilletRadius::Constant(0.2),
        })
        .collect();
    let opts = FilletOptions {
        continuity: ContinuityLevel::G2,
        ..FilletOptions::default()
    };
    fillet_edges(&brep, &specs, &opts, tol, None).expect("G2 rim fillet")
}

fn strip(result: &BlendResult) -> (FaceHandle, NurbsSurface) {
    result
        .brep
        .faces
        .iter()
        .find_map(|(h, f)| match f.surface.as_ref().map(|s| &s.geometry) {
            Some(SurfaceGeometry::Nurbs(n)) => Some((h, n.clone())),
            _ => None,
        })
        .expect("the G2 strip")
}

/// `surface` with control row `row` (counted from v = 0) of every
/// u-column moved by `offset` along the surface normal at that
/// column's rail point.
fn perturb_rail_row(surface: &NurbsSurface, row: usize, offset: f64) -> NurbsSurface {
    let mut out = surface.clone();
    let (u0, u1) = (
        surface.knots_u[0],
        surface.knots_u[surface.knots_u.len() - 1],
    );
    let v0 = surface.knots_v[0];
    for i in 0..surface.num_u {
        // Greville abscissa of column i: where its row sits.
        let greville = surface.knots_u[i + 1..=i + surface.degree_u]
            .iter()
            .sum::<f64>()
            / surface.degree_u as f64;
        let u = greville.clamp(u0, u1);
        let normal = surface
            .normal_at(u, v0)
            .normalized()
            .unwrap_or(Vec3::new(0.0, 0.0, 1.0));
        let k = i * surface.num_v + row;
        out.control_points[k] = out.control_points[k] + normal * offset;
    }
    out
}

/// Worst position and normal change along the v = 0 rail between two
/// versions of the strip, as the self-check that a perturbation kept
/// the orders it claims to keep.
fn rail_change(before: &NurbsSurface, after: &NurbsSurface) -> (f64, f64) {
    let (u0, u1) = (before.knots_u[0], before.knots_u[before.knots_u.len() - 1]);
    let v0 = before.knots_v[0];
    let (mut position, mut normal) = (0.0_f64, 0.0_f64);
    for k in 0..=64 {
        let u = u0 + (u1 - u0) * k as f64 / 64.0;
        position = position.max(before.point_at(u, v0).distance_to(&after.point_at(u, v0)));
        if let (Some(a), Some(b)) = (
            before.normal_at(u, v0).normalized(),
            after.normal_at(u, v0).normalized(),
        ) {
            normal = normal.max(a.cross(&b).length());
        }
    }
    (position, normal)
}

/// Measure the strip's evidence with `surface` swapped in.
fn evidence_with(
    result: &BlendResult,
    face: FaceHandle,
    surface: NurbsSurface,
) -> super::BlendContinuityEvidence {
    let mut brep = result.brep.clone();
    if let Some(f) = brep.faces.get_mut(face)
        && let Some(on_face) = f.surface.as_mut()
    {
        on_face.geometry = SurfaceGeometry::Nurbs(surface);
    }
    measure_blend_continuity_evidence(
        &brep,
        &[face],
        ContinuityLevel::G2,
        &ContinuityGates::commercial(),
        33,
    )
}

/// The positive fixture reads G2, the precondition of every control.
#[test]
fn the_g2_rim_fixture_reads_g2() {
    let result = g2_rim();
    let (face, surface) = strip(&result);
    let evidence = evidence_with(&result, face, surface);
    assert_eq!(evidence.achieved, ContinuityLevel::G2, "{evidence:?}");
    // The claim is certified over every gated seam, not only sampled.
    assert_eq!(
        evidence.certified_achieved,
        Some(ContinuityLevel::G2),
        "{evidence:?}"
    );
}

/// The worst certified lower bound (a value attained on the seam) of
/// quantity `k` (0 gap, 1 sine, 2 curvature) across gated seams.
fn certified_witness(evidence: &crate::BlendContinuityEvidence, k: usize) -> f64 {
    evidence
        .seams
        .iter()
        .filter_map(|s| s.certificate.as_ref())
        .map(|c| c.lower[k])
        .fold(0.0, f64::max)
}

/// Curvature only: row 2 moves, so the rail keeps its position and
/// normal exactly, and the transverse curvature jumps. The evidence
/// must fall to G1. The along-seam curvature cannot see this
/// (Meusnier), so the transverse term is what fails.
#[test]
fn a_curvature_only_perturbation_fails_g2_and_passes_g1() {
    let result = g2_rim();
    let (face, surface) = strip(&result);
    let perturbed = perturb_rail_row(&surface, 2, 5.0e-3);
    let (position, normal) = rail_change(&surface, &perturbed);
    assert!(
        position < 1e-12 && normal < 1e-12,
        "G0/G1 kept: {position:.3e} {normal:.3e}"
    );
    let evidence = evidence_with(&result, face, perturbed);
    assert_eq!(evidence.achieved, ContinuityLevel::G1, "{evidence:?}");
    assert!(
        evidence.worst_g0_gap <= evidence.gates.g0_gap,
        "{evidence:?}"
    );
    assert!(
        evidence.worst_g1_normal <= evidence.gates.g1_normal,
        "{evidence:?}"
    );
    assert!(
        evidence.worst_g2_residual > 10.0 * evidence.gates.g2_curvature,
        "the jump reads far past the gate: {evidence:?}"
    );
    // Certified failure: a curvature jump past the gate is attained.
    assert!(
        certified_witness(&evidence, 2) > evidence.gates.g2_curvature,
        "{evidence:?}"
    );
}

/// Normal only: row 1 moves by a small tilt (about 10°), so the rail
/// keeps its position but the normal turns. The evidence must fall to
/// G0, with the seam still gated as intended-smooth.
#[test]
fn a_normal_only_perturbation_fails_g1() {
    let result = g2_rim();
    let (face, surface) = strip(&result);
    let perturbed = perturb_rail_row(&surface, 1, 1.0e-2);
    let (position, normal) = rail_change(&surface, &perturbed);
    assert!(position < 1e-12, "G0 kept: {position:.3e}");
    assert!(normal > 5.0e-2, "the normal really turned: {normal:.3e}");
    let evidence = evidence_with(&result, face, perturbed);
    assert_eq!(evidence.achieved, ContinuityLevel::G0, "{evidence:?}");
    assert!(
        evidence.worst_g0_gap <= evidence.gates.g0_gap,
        "{evidence:?}"
    );
    assert_eq!(evidence.sharp_seams, 0, "still gated: {evidence:?}");
    // Certified failure: a normal turn past the gate is attained.
    assert!(
        certified_witness(&evidence, 1) > evidence.gates.g1_normal,
        "{evidence:?}"
    );
}

/// A tilt past the 30° smooth-intent bound reclassifies the rail as an
/// intended-sharp crease, which is G0-gated only. The classification
/// alone would then pass G1 on the remaining seams. The sharp-seam
/// census is the guard, and it must record the change.
#[test]
fn a_tilt_past_the_intent_bound_is_caught_by_the_sharp_seam_census() {
    let result = g2_rim();
    let (face, surface) = strip(&result);
    let before = evidence_with(&result, face, surface.clone());
    let perturbed = perturb_rail_row(&surface, 1, 5.0e-2);
    let after = evidence_with(&result, face, perturbed);
    assert_eq!(before.sharp_seams, 0, "{before:?}");
    assert!(
        after.sharp_seams > before.sharp_seams,
        "a rail collapsing into 'sharp' must show in the census: {after:?}"
    );
}

/// PK-499 exploration E2 acceptance: the G2 cylinder rim's rails pair
/// a rational-arc strip parameter with an angle (fraction
/// correspondence fails by ~1e-3), and must certify under projected
/// correspondence at the tight default targets: G1 below 1e-6 and G2
/// below 1e-3 over the whole seam.
#[test]
fn the_g2_rim_rails_certify_under_projected_correspondence() {
    use popkern_core::geometry::certify::SeamCorrespondence;
    use popkern_core::topology::edge_seam_certificate::{certify_edge_seam, seam_certify_options};
    let result = g2_rim();
    let (face, _) = strip(&result);
    let brep = &result.brep;
    let options = seam_certify_options(brep);
    let mut rails = 0;
    for (edge, data) in brep.edges.iter() {
        let he = &brep.half_edges[data.half_edge];
        let faces = [
            brep.loops[he.loop_].face,
            brep.loops[brep.half_edges[he.twin].loop_].face,
        ];
        if !faces.contains(&face) || faces[0] == faces[1] {
            continue;
        }
        let cert = certify_edge_seam(brep, edge, &options);
        if cert.lower[1] > 0.5 {
            continue; // the strip's designed end crease
        }
        rails += 1;
        assert!(cert.is_certified(), "{cert:?}");
        assert_eq!(cert.correspondence, SeamCorrespondence::Projected);
        assert!(cert.g1_upper < 1e-6 && cert.g2_upper < 1e-3, "{cert:?}");
    }
    assert_eq!(rails, 2, "the strip has two rails");
}
