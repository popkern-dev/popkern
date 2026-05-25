use super::*;
use crate::tolerance::Tolerance;
use crate::topology::{BrepFaceUseTrimAuthorityFindingKind, analyze_brep_face_use_trim_authority};

fn options() -> BrepFaceUseTrimAuthorityOptions {
    BrepFaceUseTrimAuthorityOptions {
        tolerance: 1e-6,
        samples: 9,
    }
}

#[test]
fn repair_materializes_explicit_authority_for_reused_loop_on_same_surface() {
    let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let faces = brep.faces.keys().take(2).collect::<Vec<_>>();
    let source_loop = brep.faces[faces[0]].outer_loop;
    brep.faces[faces[1]].outer_loop = source_loop;
    brep.faces[faces[1]].surface = brep.faces[faces[0]].surface.clone();

    let before = analyze_brep_face_use_trim_authority(&brep, options());
    assert!(before.findings.iter().any(|finding| {
        finding.kind == BrepFaceUseTrimAuthorityFindingKind::LoopOwnershipMismatch
    }));

    let repaired = repair_brep_face_use_trim_authority_transactionally(&brep, options())
        .expect("same-surface loop reuse should be materializable as face-use authority");

    assert!(
        repaired.report.after.is_clean(),
        "{:#?}",
        repaired.report.after
    );
    assert!(repaired.report.materialized_face_uses > 0);
    assert_eq!(
        repaired.report.after.explicit_face_use_authorities,
        repaired.report.materialized_face_uses
    );
    assert!(
        analyze_brep_face_use_trim_authority(&brep, options()).loop_ownership_mismatches > 0,
        "transaction must not mutate the input BRep"
    );
}

#[test]
fn repair_declines_when_face_surface_cannot_support_loop_vertices() {
    let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let faces = brep.faces.keys().take(2).collect::<Vec<_>>();
    brep.faces[faces[1]].outer_loop = brep.faces[faces[0]].outer_loop;

    let failure = repair_brep_face_use_trim_authority_transactionally(&brep, options())
        .expect_err("different-surface loop reuse must not be silently accepted");

    assert!(failure.report.before.loop_ownership_mismatches > 0);
    assert!(failure.report.after.has_blocking_findings());
}
