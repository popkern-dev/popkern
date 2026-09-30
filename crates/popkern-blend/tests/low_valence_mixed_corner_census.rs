//! PK-591 / PK-593: the low-valence mixed-radius apex census gate.
//!
//! Pyramid apexes n = 3, 4 with mixed, variable and wide radii at G2.
//! Before PK-591 the setback carve declined most of them ("expected 2
//! arc/boundary crossings, found 3": the S-shaped mixed-radius arc
//! grazed a strip boundary mid-arc) and the corner fell to a creased
//! vertex-blend Gregory cap (body G0). PK-593 found a second cause
//! beneath: direct-apex promotion claimed these G2 corners although
//! DirectApex v1 builds only on cylinder strips, and its fallback built
//! the strips on direct-apex anchors (every rail creased 13°). Each row
//! pins the whole-blend verdict (PK-499 certificates, commercial gates)
//! and that the setback rung serves.

use std::collections::BTreeSet;

use popkern_blend::blend_continuity_evidence::{
    ContinuityGates, measure_blend_continuity_evidence,
};
use popkern_blend::{ContinuityLevel, FilletOptions, FilletRadius, FilletSpec, fillet_edges};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};

const MIXED: [f64; 4] = [0.04, 0.06, 0.05, 0.07];
const VARIABLE: [f64; 4] = [0.02, 0.04, 0.06, 0.04];
const WIDE: [f64; 4] = [0.05, 0.15, 0.08, 0.2];

fn verdict(n: usize, h: f64, radii: &[f64]) -> ContinuityLevel {
    let tol = Tolerance::default();
    let body = BRep::n_gonal_pyramid(n, 1.0, h, tol).expect("pyramid");
    let apex = body
        .vertices
        .iter()
        .max_by(|a, b| a.1.position.z.total_cmp(&b.1.position.z))
        .map(|(v, _)| v)
        .expect("apex");
    let mut edges: Vec<_> = body
        .half_edges
        .values()
        .filter(|he| he.origin == apex)
        .map(|he| he.edge)
        .collect();
    edges.sort();
    edges.dedup();
    let specs: Vec<FilletSpec> = edges
        .iter()
        .zip(radii.iter().cycle())
        .map(|(&edge, &r)| FilletSpec {
            edge,
            radius: FilletRadius::Constant(r),
        })
        .collect();
    let pre: BTreeSet<FaceHandle> = body.faces.keys().collect();
    let options = FilletOptions {
        continuity: ContinuityLevel::G2,
        ..Default::default()
    };
    let label = format!("pyr{n} h{h} {radii:?}");
    let result = fillet_edges(&body, &specs, &options, tol, None)
        .unwrap_or_else(|e| panic!("{label}: {e:?}"));
    assert!(
        !result.corner_patches.is_empty() && result.corner_patches.iter().all(|p| p.setback_cap),
        "{label}: the setback rung serves: {:?}",
        result.corner_patches
    );
    let blend: Vec<FaceHandle> = result
        .brep
        .faces
        .keys()
        .filter(|f| !pre.contains(f))
        .collect();
    let evidence = measure_blend_continuity_evidence(
        &result.brep,
        &blend,
        ContinuityLevel::G2,
        &ContinuityGates::commercial(),
        33,
    );
    assert!(result.brep.validate().is_ok(), "{label}");
    evidence.achieved
}

#[test]
fn triangular_apexes_serve_g2_on_the_setback_rung() {
    // Measured 2026-09-30. h 0.4 variable/wide serve the fitted cap at
    // certified G1 (KP20 declines there, the fit is G1 on these).
    for (h, radii, level) in [
        (0.4, &MIXED[..3], ContinuityLevel::G2),
        (0.4, &VARIABLE[..3], ContinuityLevel::G1),
        (0.4, &WIDE[..3], ContinuityLevel::G1),
        (0.8, &MIXED[..3], ContinuityLevel::G2),
        (0.8, &VARIABLE[..3], ContinuityLevel::G2),
        (0.8, &WIDE[..3], ContinuityLevel::G2),
        (1.2, &VARIABLE[..3], ContinuityLevel::G2),
        (1.2, &WIDE[..3], ContinuityLevel::G2),
    ] {
        assert_eq!(verdict(3, h, radii), level, "pyr3 h{h} {radii:?}");
    }
}

#[test]
fn quad_apexes_certify_at_least_g1_on_the_setback_rung() {
    // Measured 2026-09-30: every row certifies G1 (KP20's measured G1
    // downgrade, C3 class, or the fitted cap at h 0.4 wide); none G0.
    for h in [0.4, 0.8, 1.2] {
        for radii in [&MIXED[..], &VARIABLE[..], &WIDE[..]] {
            assert!(
                verdict(4, h, radii) >= ContinuityLevel::G1,
                "pyr4 h{h} {radii:?}"
            );
        }
    }
}
