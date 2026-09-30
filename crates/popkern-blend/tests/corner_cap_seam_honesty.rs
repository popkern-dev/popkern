//! PK-609 + PK-610: corner-cap seams are judged as smooth by
//! construction, and the fitted cap that then serves the creased cube
//! corners carries sound strip pcurves.
//!
//! The G1 cube corners 0.05/0.1/0.2 and 0.05/0.2/0.3 served a
//! three-sector direct-apex cap whose sector seams certify normal-sine
//! lower bounds 0.48-0.9998; the evidence recorded those seams as
//! designed creases (past the 30° blend-to-blend bound) and read the
//! blend G1. With cap-internal seams gated at any angle (PK-609) the
//! direct-apex cap reads G0, the serving crease gate takes the
//! exclusion re-run, and the tangent fit serves. That fit had left two
//! kinds of broken strip pcurve (PK-610): a pcurve still spanning the
//! pre-carve edge, snapped at its ends into a fold, and a rim pcurve
//! authored backwards (an out-voted vertex sense) and snapped the same
//! way. Both validated; the evidence read up to 7.2e-2 of gap.

use std::collections::BTreeSet;

use popkern_blend::blend_continuity_evidence::{
    BlendContinuityEvidence, ContinuityGates, measure_blend_continuity_evidence,
};
use popkern_blend::{
    BlendResult, ContinuityLevel, CornerCapPreference, CornerConstruction, CornerMethod,
    FilletOptions, FilletRadius, FilletSpec, fillet_edges,
};
use popkern_core::geometry::Point3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};

fn cube_corner_g1(
    radii: [f64; 3],
    options: FilletOptions,
) -> (BlendResult, BlendContinuityEvidence) {
    let tol = Tolerance::default();
    let cube = BRep::unit_cube(tol).expect("cube");
    let corner = cube
        .vertices
        .iter()
        .find(|(_, v)| v.position.distance_to(&Point3::new(1.0, 1.0, 1.0)) < 1e-9)
        .map(|(h, _)| h)
        .expect("corner");
    let mut edges: Vec<_> = cube
        .half_edges
        .values()
        .filter(|he| he.origin == corner)
        .map(|he| he.edge)
        .collect();
    edges.sort();
    edges.dedup();
    let specs: Vec<FilletSpec> = edges
        .iter()
        .zip(radii)
        .map(|(&edge, r)| FilletSpec {
            edge,
            radius: FilletRadius::Constant(r),
        })
        .collect();
    let options = FilletOptions {
        continuity: ContinuityLevel::G1,
        ..options
    };
    let result = fillet_edges(&cube, &specs, &options, tol, None).expect("fillet");
    let pre: BTreeSet<FaceHandle> = cube.faces.keys().collect();
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
    (result, evidence)
}

const CREASED: [[f64; 3]; 2] = [[0.05, 0.1, 0.2], [0.05, 0.2, 0.3]];

#[test]
fn creased_direct_apex_caps_read_g0() {
    // Construction observed with the serving retry off: the direct-apex
    // cap and its creased sector seams, now gated.
    for radii in CREASED {
        let (result, evidence) = cube_corner_g1(
            radii,
            FilletOptions {
                serving_self_intersection_retry: false,
                ..Default::default()
            },
        );
        assert!(
            result.corner_patches.iter().all(|p| !p.setback_cap),
            "{radii:?}: the direct-apex cap is observed"
        );
        assert_eq!(evidence.achieved, ContinuityLevel::G0, "{radii:?}");
        assert!(
            evidence.worst_g1_normal > 0.7,
            "{radii:?}: {}",
            evidence.worst_g1_normal
        );
    }
}

#[test]
fn creased_corners_serve_a_certified_g1_fit_with_sound_pcurves() {
    for radii in CREASED {
        for (label, corner_cap) in [
            ("default", CornerCapPreference::default()),
            (
                "setback",
                CornerCapPreference {
                    construction: CornerConstruction::Setback,
                    strict: false,
                },
            ),
        ] {
            let (result, evidence) = cube_corner_g1(
                radii,
                FilletOptions {
                    corner_cap,
                    ..Default::default()
                },
            );
            assert!(
                result
                    .corner_patches
                    .iter()
                    .all(|p| p.setback_cap && p.method == CornerMethod::PlateMin),
                "{radii:?} {label}: the tangent fit serves: {:?}",
                result.corner_patches
            );
            assert_eq!(
                evidence.certified_achieved,
                Some(ContinuityLevel::G1),
                "{radii:?} {label}"
            );
            assert!(
                evidence.worst_g0_gap <= 1e-5,
                "{radii:?} {label}: pcurve gap {:.3e}",
                evidence.worst_g0_gap
            );
            result.brep.validate().expect("valid");
        }
    }
}
