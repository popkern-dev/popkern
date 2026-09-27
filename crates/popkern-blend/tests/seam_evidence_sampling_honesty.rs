//! PK-579: the sampled seam evidence must not be blind where the
//! residual lives. Before, the instrument skipped curvature samples
//! whose finite-difference estimates did not converge and inset every
//! seam by 5%. On pentagon5_blunt at the G2 request it read cap-rim g2
//! 0.059 where the whole-edge certificate PROVES >= 0.29 (first 3% of
//! the seam). Now curvature comes from exact jets and the inset is
//! 1e-3.
//!
//! Invariant: the certificate's lower bound is a value attained on the
//! seam, so a sampled worst far below it is a blind sampler. Every
//! gated seam's sampled g2 must reach 80% of its certified lower bound,
//! with no skipped curvature samples. The blind sampler read 20%. At 33
//! samples a narrow peak between samples reads ~90%
//! (pentagon5_variable's cap-cap seams), which is sampling phase, not
//! blindness.

use popkern_blend::blend_continuity_evidence::{
    ContinuityGates, measure_blend_continuity_evidence,
};
use popkern_blend::{ContinuityLevel, FilletOptions, FilletRadius, FilletSpec, fillet_edges};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};

fn apex_g2_radii(n: usize, height: f64, radii: &[f64]) -> (BRep, Vec<FaceHandle>) {
    let tol = Tolerance::default();
    let pyr = BRep::n_gonal_pyramid(n, 1.0, height, tol).expect("pyramid");
    let apex = pyr
        .vertices
        .iter()
        .max_by(|a, b| a.1.position.z.total_cmp(&b.1.position.z))
        .map(|(h, _)| h)
        .expect("apex");
    let mut edges: Vec<_> = pyr
        .half_edges
        .values()
        .filter(|he| he.origin == apex)
        .map(|he| he.edge)
        .collect();
    edges.sort();
    edges.dedup();
    let specs: Vec<FilletSpec> = edges
        .into_iter()
        .zip(radii.iter())
        .map(|(edge, &r)| FilletSpec {
            edge,
            radius: FilletRadius::Constant(r),
        })
        .collect();
    let options = FilletOptions {
        continuity: ContinuityLevel::G2,
        ..FilletOptions::default()
    };
    let pre: std::collections::BTreeSet<_> = pyr.faces.keys().collect();
    let result = fillet_edges(&pyr, &specs, &options, tol, None).expect("fillet");
    let blend = result
        .brep
        .faces
        .keys()
        .filter(|f| !pre.contains(f))
        .collect();
    (result.brep, blend)
}

#[test]
fn sampled_seam_evidence_reaches_the_certified_lower_bound() {
    for (name, n, height, radii) in [
        ("pentagon_flagship", 5, 1.2, vec![0.16; 5]),
        ("pentagon5_blunt", 5, 3.0, vec![0.05; 5]),
        ("hexagon6_blunt", 6, 3.0, vec![0.05; 6]),
        (
            "hexagon6_mixed",
            6,
            0.8,
            vec![0.04, 0.06, 0.05, 0.07, 0.05, 0.06],
        ),
        (
            "pentagon5_variable",
            5,
            0.8,
            vec![0.02, 0.04, 0.06, 0.04, 0.02],
        ),
    ] {
        let (brep, blend) = apex_g2_radii(n, height, &radii);
        let evidence = measure_blend_continuity_evidence(
            &brep,
            &blend,
            ContinuityLevel::G2,
            &ContinuityGates::commercial(),
            101,
        );
        let mut worst_ratio = f64::INFINITY;
        let worst_lower = evidence
            .seams
            .iter()
            .filter_map(|s| s.certificate.as_ref().map(|c| c.lower[2]))
            .fold(0.0, f64::max);
        eprintln!(
            "{name}: achieved {:?}, worst sampled g2 {:.3e}, worst certified lower {worst_lower:.3e}, unmeasured {}",
            evidence.achieved, evidence.worst_g2_residual, evidence.unmeasured_seams
        );
        for seam in evidence.seams.iter().filter(|s| s.gated) {
            assert_eq!(
                seam.skipped_curvature_samples, 0,
                "{name}: seam {:?} skipped curvature samples",
                seam.edge
            );
            let Some(cert) = seam.certificate.as_ref() else {
                continue;
            };
            let lower = cert.lower[2];
            if lower < 1e-3 {
                continue; // below any gate: sampling noise dominates
            }
            let ratio = seam.worst_g2_residual / lower;
            worst_ratio = worst_ratio.min(ratio);
            assert!(
                ratio >= 0.8,
                "{name}: seam {:?} sampled g2 {:.3e} under the certified lower bound {lower:.3e}",
                seam.edge,
                seam.worst_g2_residual
            );
        }
        eprintln!("{name}: worst sampled/certified-lower ratio {worst_ratio:.3}");
    }
}
