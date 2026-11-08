//! Shell and offset take solids, thicken solids or sheets; each refuses
//! a general body by its kind (PK-1101 S5).

use popkern_core::error::CoreError;
use popkern_core::geometry::transform::Transform3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, BodyKind};

fn tol() -> Tolerance {
    Tolerance::default()
}

use popkern_shell::offset::offset_solid;
use popkern_shell::thicken::thicken;
use popkern_shell::{ShellDirection, ShellError, ShellOptions, shell};

/// A cube beside a sheet: lumps of two kinds, a general body.
fn general_body() -> BRep {
    let mut body = BRep::unit_cube(tol()).expect("cube");
    let mut sheet = BRep::unit_cube(tol())
        .and_then(|c| c.transform(&Transform3::translation(3.0, 0.0, 0.0)))
        .expect("cube");
    let face = sheet.faces.keys().next().expect("a face");
    sheet.retire_faces_to_boundary(&[face]).expect("retires");
    sheet.settle_open_shells();
    body.append_body(&sheet);
    assert_eq!(body.body_kind(), BodyKind::General);
    body
}

fn assert_refused<T: std::fmt::Debug>(result: Result<T, ShellError>, operation: &str) {
    match result {
        Err(ShellError::Core(CoreError::UnsupportedBodyKind {
            operation: op,
            found: BodyKind::General,
            ..
        })) => assert_eq!(op, operation),
        other => panic!("{operation}: expected the kind gate, got {other:?}"),
    }
}

#[test]
fn shells_offsets_and_thickens_refuse_a_general_body() {
    let body = general_body();
    let options = ShellOptions::default();
    assert_refused(
        shell(
            &body,
            0.1,
            ShellDirection::Inward,
            &[],
            &[],
            &options,
            tol(),
            None,
        ),
        "shell",
    );
    assert_refused(offset_solid(&body, 0.1, &options, tol()), "offset");
    assert_refused(thicken(&body, 0.1, &options, tol()), "thicken");
}
