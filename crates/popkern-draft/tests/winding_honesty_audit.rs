//! PK-245: winding-flag honesty audit at the draft serving boundary.
//!
//! The stored orientation flag must satisfy the interior-left
//! invariant (outer loop CCW in UV ⇔ flag true) at every operator
//! serving boundary. Draft tilts face surfaces in place and
//! recomputes trims — a rewritten pcurve set with an untouched flag
//! is exactly the class that ships dishonest flags. A failure lists
//! the dishonest faces with their surface class.

use popkern_core::geometry::Vec3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_draft::{DraftOptions, DraftSpec, NeutralReference, draft_face, draft_faces};

fn dishonest_report(brep: &BRep) -> Vec<String> {
    brep.audit_winding_honesty()
        .into_iter()
        .filter(|v| v.honest == Some(false))
        .map(|v| {
            let surface = brep
                .faces
                .get(v.face)
                .and_then(|f| f.surface.as_ref())
                .map(|s| format!("{:?}", std::mem::discriminant(&s.geometry)))
                .unwrap_or_else(|| "<no-surface>".into());
            format!(
                "{:?} [{surface}] area {:.3e} flag {}",
                v.face, v.signed_uv_area, v.orientation
            )
        })
        .collect()
}

fn assert_honest(brep: &BRep, fixture: &str) {
    let report = dishonest_report(brep);
    assert!(
        report.is_empty(),
        "{fixture}: {} dishonest face(s):\n  {}",
        report.len(),
        report.join("\n  ")
    );
}

fn tol() -> Tolerance {
    Tolerance::default()
}

/// Single drafted cube face (auto neutral line).
#[test]
fn draft_single_face_is_winding_honest() {
    let brep = BRep::unit_cube(tol()).unwrap();
    let face = brep.faces.keys().next().unwrap();
    let result = draft_face(
        &brep,
        face,
        5.0_f64.to_radians(),
        Vec3::new(0.0, 0.0, 1.0),
        NeutralReference::Auto,
        tol(),
        None,
    )
    .expect("draft");
    assert_honest(&result.brep, "draft_single_face");
}

/// All four cube side faces drafted together (mold-taper shape).
#[test]
fn draft_four_sides_is_winding_honest() {
    let brep = BRep::unit_cube(tol()).unwrap();
    let pull = Vec3::new(0.0, 0.0, 1.0);
    let specs: Vec<DraftSpec> = brep
        .faces
        .iter()
        .filter(|(fh, _)| {
            // Side faces: normal perpendicular to the pull direction.
            let positions = brep.face_outer_loop_positions(*fh);
            let zs: Vec<f64> = positions.iter().map(|p| p.z).collect();
            zs.iter().any(|z| (z - 0.0).abs() < 1e-9) && zs.iter().any(|z| (z - 1.0).abs() < 1e-9)
        })
        .map(|(fh, _)| DraftSpec {
            face: fh,
            angle: 4.0_f64.to_radians(),
            pull_direction: pull,
            neutral: NeutralReference::Auto,
            angle_variation: popkern_draft::AngleVariation::Constant,
        })
        .collect();
    assert_eq!(specs.len(), 4, "cube has four side faces");
    let result =
        draft_faces(&brep, &specs, &DraftOptions::default(), tol(), None).expect("draft sides");
    assert_honest(&result.brep, "draft_four_sides");
}
