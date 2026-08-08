//! PK-245: winding-flag honesty audit at the shell serving boundary.
//!
//! The stored orientation flag must satisfy the interior-left
//! invariant (outer loop CCW in UV ⇔ flag true) at every operator
//! serving boundary. Shell authors an inner skin whose material side
//! is the MIRROR of the outer skin's, plus side walls — both are
//! high-risk authors for dishonest flags. A failure lists the
//! dishonest faces with their surface class.

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_shell::offset::offset_solid;
use popkern_shell::thicken::thicken;
use popkern_shell::{ShellDirection, ShellOptions, shell};

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

/// Hollow cube (no removed faces): outer + mirrored inner skin.
#[test]
fn shell_hollow_cube_is_winding_honest() {
    let brep = BRep::unit_cube(tol()).unwrap();
    let result = shell(
        &brep,
        0.1,
        ShellDirection::Inward,
        &[],
        &[],
        &ShellOptions::default(),
        tol(),
        None,
    )
    .expect("hollow shell");
    assert_honest(&result.brep, "shell_hollow_cube");
}

/// Whole-body offset: same-orientation clone, no mirror expected.
#[test]
fn offset_solid_is_winding_honest() {
    let brep = BRep::unit_cube(tol()).unwrap();
    let result = offset_solid(&brep, 0.2, &ShellOptions::default(), tol()).expect("offset solid");
    assert_honest(&result.brep, "offset_solid");
}

/// Thicken: outer skin same-orientation + mirrored inner skin.
#[test]
fn thicken_is_winding_honest() {
    let brep = BRep::unit_cube(tol()).unwrap();
    let result = thicken(&brep, 0.2, &ShellOptions::default(), tol()).expect("thicken");
    assert_honest(&result.brep, "thicken");
}

/// Open box (one face removed): inner skin + side-wall authors.
#[test]
fn shell_open_box_is_winding_honest() {
    let brep = BRep::unit_cube(tol()).unwrap();
    let removed = brep.faces.keys().next().unwrap();
    let result = shell(
        &brep,
        0.1,
        ShellDirection::Inward,
        &[removed],
        &[],
        &ShellOptions::default(),
        tol(),
        None,
    )
    .expect("open-box shell");
    assert_honest(&result.brep, "shell_open_box");
}
