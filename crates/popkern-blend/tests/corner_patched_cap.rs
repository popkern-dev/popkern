//! PK-619: steep apexes whose star caps bend past the fairness bound
//! serve the corner-patched cap. A thin strip's arc between long
//! bridges gives the star's two sectors there rim halves 3.5–6x apart,
//! and the cap bends in them (stage 2, honest reading: 4.53–50 / r_min).
//! A four-sided corner patch over each thin side, joined to the main cap
//! along a seam ribbon, brings every row under the bound with every
//! seam, the patch seams included, certified G1.

use std::collections::BTreeSet;

use popkern_blend::blend_continuity_evidence::{
    ContinuityGates, measure_blend_continuity_evidence,
};
use popkern_blend::{
    ApexConstructionDiagnostic, BlendResult, ContinuityLevel, CornerSubPatchRole, FilletOptions,
    FilletRadius, FilletSpec, SetbackFitBorder, SetbackFitMode, SetbackFitOutcome, fillet_edges,
};
use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::healing::curvature::principal_curvatures;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};

fn pyramid_apex(n: usize, h: f64, radii: &[f64]) -> (BlendResult, BTreeSet<FaceHandle>) {
    let body = BRep::n_gonal_pyramid(n, 1.0, h, Tolerance::default()).expect("pyramid");
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
        .into_iter()
        .zip(radii.iter().cycle())
        .map(|(edge, &r)| FilletSpec {
            edge,
            radius: FilletRadius::Constant(r),
        })
        .collect();
    let options = FilletOptions {
        continuity: ContinuityLevel::G1,
        ..Default::default()
    };
    let pre = body.faces.keys().collect();
    let result = fillet_edges(&body, &specs, &options, Tolerance::default(), None).expect("fillet");
    (result, pre)
}

/// Largest principal curvature over the corner caps times `r_min`, on
/// the rung's graded grid (PK-621).
fn cap_fairness(result: &BlendResult, r_min: f64) -> f64 {
    let mut fractions: Vec<f64> = (1..32).map(|i| i as f64 / 32.0).collect();
    for d in [1.0 / 64.0, 1.0 / 256.0, 1.0 / 1024.0] {
        fractions.extend([d, 1.0 - d]);
    }
    let mut worst = 0.0_f64;
    for patch in &result.corner_patches {
        let Some(surface) = result
            .brep
            .faces
            .get(patch.face)
            .and_then(|f| f.surface.as_ref())
        else {
            continue;
        };
        let geometry: &dyn SurfaceEval = &surface.geometry;
        let ((u0, u1), (v0, v1)) = geometry.parameter_range();
        for &fu in &fractions {
            for &fv in &fractions {
                let (u, v) = (u0 + (u1 - u0) * fu, v0 + (v1 - v0) * fv);
                if let Some((k1, k2)) = principal_curvatures(geometry, u, v) {
                    worst = worst.max(k1.abs().max(k2.abs()));
                }
            }
        }
    }
    worst * r_min
}

/// The row serves the corner-patched cap on a bowed border (the fairest
/// of the patched attempts, PK-625):
/// certified G1 over the whole blend, under the bound, fairer than the
/// stage-2 star (`star`), its faces naming their piece.
fn assert_corner_patched(n: usize, h: f64, radii: &[f64], star: f64, at_most: f64) {
    let name = format!("pyramid n={n} h={h} {radii:?}");
    let (result, pre) = pyramid_apex(n, h, radii);
    let blend: Vec<FaceHandle> = result
        .brep
        .faces
        .keys()
        .filter(|f| !pre.contains(f))
        .collect();
    let evidence = measure_blend_continuity_evidence(
        &result.brep,
        &blend,
        ContinuityLevel::G1,
        &ContinuityGates::commercial(),
        33,
    );
    assert_eq!(
        evidence.achieved,
        ContinuityLevel::G1,
        "{name}: {evidence:?}"
    );
    assert!(result.brep.validate().is_ok(), "{name}");
    let record = result
        .repair_report
        .apex_diagnostics
        .iter()
        .find_map(|d| match d {
            ApexConstructionDiagnostic::SetbackFit(record) => Some(record),
            _ => None,
        })
        .unwrap_or_else(|| panic!("{name}: the rung's record"));
    let served = record
        .attempts
        .iter()
        .find(|a| matches!(a.outcome, SetbackFitOutcome::Served(_)))
        .unwrap_or_else(|| panic!("{name}: {record:#?}"));
    assert!(
        served.mode == SetbackFitMode::TangentCornerPatched
            && matches!(served.border, SetbackFitBorder::BowedBridges(_)),
        "{name}: {:?} on {:?}",
        served.mode,
        served.border
    );
    let r_min = radii.iter().take(n).copied().fold(f64::INFINITY, f64::min);
    let fairness = cap_fairness(&result, r_min);
    assert!(
        fairness <= at_most && fairness < star,
        "{name}: cap {fairness:.2} / r_min (star {star})"
    );
    let recorded = served.interior.unwrap_or(f64::NAN);
    assert!((recorded - fairness).abs() < 1e-9, "{name}: {recorded}");
    // Every face names its piece: the main cap's sectors and at least one
    // corner patch's.
    let roles: Vec<_> = result
        .corner_patches
        .iter()
        .filter_map(|p| p.sub_patch_provenance)
        .collect();
    assert_eq!(roles.len(), result.corner_patches.len(), "{name}");
    assert!(roles.iter().any(|r| r.role == CornerSubPatchRole::Sector));
    assert!(
        roles
            .iter()
            .any(|r| r.role == CornerSubPatchRole::CornerPatchSector && r.piece_index > 0)
    );
}

#[test]
fn steep_trihedral_variable_apex_serves_the_corner_patched_cap() {
    // Stage 2 served the round-bow star at 6.29 / r_min; stage 3 measured
    // 2.51, and the fairest patched attempt (PK-625) 1.20.
    assert_corner_patched(3, 3.0, &[0.02, 0.04, 0.06], 6.29, 1.5);
}

#[test]
fn steep_trihedral_wide_apex_serves_the_corner_patched_cap() {
    // Stage 2: 8.27; stage 3 measured 2.49; PK-625: 1.24.
    assert_corner_patched(3, 3.0, &[0.05, 0.15, 0.08], 8.27, 1.5);
}

#[test]
fn steep_pentagonal_wide_apex_serves_the_corner_patched_cap() {
    // Stage 2: 22.3, peaking on the bridge between the two thinnest
    // arcs; measured 1.47.
    assert_corner_patched(5, 3.0, &[0.05, 0.15, 0.08, 0.2, 0.06], 22.3, 2.0);
}

#[test]
fn steep_pentagonal_mixed_apex_serves_the_corner_patched_cap() {
    // Stage 2: 4.53; stage 3 measured 2.56; PK-625: 2.14.
    assert_corner_patched(5, 3.0, &[0.04, 0.06, 0.05, 0.07, 0.04], 4.53, 2.5);
}

#[test]
fn steep_hexagonal_variable_apex_serves_the_corner_patched_cap() {
    // Stage 2: 50.4; measured 1.69 (PK-625: 1.67).
    assert_corner_patched(6, 3.0, &[0.02, 0.04, 0.06, 0.04, 0.02, 0.05], 50.4, 2.0);
}
