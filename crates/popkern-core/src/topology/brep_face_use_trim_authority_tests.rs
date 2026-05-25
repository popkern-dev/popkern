use super::*;
use crate::tolerance::Tolerance;

fn options() -> BrepFaceUseTrimAuthorityOptions {
    BrepFaceUseTrimAuthorityOptions {
        tolerance: 1e-6,
        samples: 9,
    }
}

#[test]
fn clean_unit_cube_has_coherent_face_use_trim_authority() {
    let brep = BRep::unit_cube(Tolerance::default()).unwrap();

    let report = analyze_brep_face_use_trim_authority(&brep, options());

    assert!(report.is_clean(), "{report:#?}");
    assert_eq!(report.inspected_face_uses, 24);
}

#[test]
fn detects_loop_reused_by_wrong_face_as_trim_authority_incoherence() {
    let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let faces = brep.faces.keys().take(2).collect::<Vec<_>>();
    let source_loop = brep.faces[faces[0]].outer_loop;
    brep.faces[faces[1]].outer_loop = source_loop;

    let report = analyze_brep_face_use_trim_authority(&brep, options());

    assert!(
        report.loop_ownership_mismatches > 0,
        "wrong-face loop reuse must be explicit trim-authority evidence: {report:#?}"
    );
    assert!(
        report.findings.iter().any(|finding| {
            finding.face == faces[1]
                && finding.loop_ == source_loop
                && finding.kind == BrepFaceUseTrimAuthorityFindingKind::LoopOwnershipMismatch
        }),
        "finding should identify the face-use that borrowed another face's loop: {report:#?}"
    );
}
