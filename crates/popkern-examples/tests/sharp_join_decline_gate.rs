//! PK-325 row 473: gap edges the sharp (Intersection) join cannot
//! close decline per-edge with a typed `JoinFailed` naming the arc
//! alternative — and no partial-join body is served. The wavy plate is
//! the acceptance body: its fitted NURBS top offset is domain-bounded,
//! so the window extensions never meet the side-plane offsets.

use popkern_core::tolerance::Tolerance;
use popkern_examples::offset_corpus::corpus;
use popkern_shell::{JoinStyle, ShellError, ShellOptions};

#[test]
fn intersection_join_declines_wavy_plate_with_typed_error() {
    let fixtures = corpus(Tolerance::default()).unwrap();
    let wavy = fixtures
        .iter()
        .find(|f| f.name == "wavy_plate")
        .expect("corpus has wavy_plate");

    let options = ShellOptions::default().with_join_style(JoinStyle::Intersection);
    let result = popkern_shell::thicken::thicken(&wavy.brep, 0.02, &options, Tolerance::default());
    match result {
        Err(ShellError::JoinFailed { style, declines }) => {
            assert_eq!(style, JoinStyle::Intersection);
            assert!(!declines.is_empty(), "per-edge declines must be named");
            assert!(
                declines.iter().all(|d| d.arc_alternative),
                "every decline names the arc join as the working alternative: {declines:?}"
            );
            assert!(
                format!(
                    "{}",
                    ShellError::JoinFailed {
                        style,
                        declines: declines.clone()
                    }
                )
                .contains("JoinStyle::Arc serves this edge"),
                "the Display text points the caller at the alternative"
            );
        }
        other => panic!("expected typed JoinFailed decline, got {other:?}"),
    }

    // The named alternative is real: the arc default serves this body
    // (the corpus campaign gate holds it to validate_all + the SDF
    // oracle; here the decline/serve pairing is the pin).
    let arc = popkern_shell::thicken::thicken(
        &wavy.brep,
        0.02,
        &ShellOptions::default(),
        Tolerance::default(),
    )
    .expect("arc join serves the wavy plate");
    assert!(arc.brep.validate_all().is_empty());
}
