//! The serve post-check on a body known to validate (PK-815).

use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use super::{PostCheckReport, post_check_gates};

/// A cylinder whose side surface sits 0.01 off its rim circles: it
/// passes `validate()` (no check reads curve-surface agreement) and
/// fails the curve-surface gate.
fn cylinder_with_lying_side() -> BRep {
    let mut b = BRep::cylinder(1.0, 2.0, Tolerance::default()).unwrap();
    for face in b.faces.values_mut() {
        if let Some(s) = face.surface.as_mut()
            && let SurfaceGeometry::Cylinder { radius, .. } = &mut s.geometry
        {
            *radius = 1.01;
        }
    }
    b
}

fn messages(r: &PostCheckReport) -> Vec<String> {
    r.findings.iter().map(|f| f.message.clone()).collect()
}

/// A body known to pass `validate()` skips the gates `validate()`
/// contains and keeps the curve-surface check: the same findings as
/// the full gates, and a curve-surface lie is still reported.
#[test]
fn a_validated_body_keeps_the_curve_surface_gate() {
    let tol = Tolerance::default();
    for body in [
        BRep::unit_cube(tol).unwrap(),
        BRep::torus(3.0, 1.0, tol).unwrap(),
        cylinder_with_lying_side(),
    ] {
        assert!(body.validate().is_ok());
        let full = post_check_gates(&body, &tol, &Default::default(), false);
        let given = post_check_gates(&body, &tol, &Default::default(), true);
        assert_eq!(messages(&given), messages(&full));
    }
    let given = post_check_gates(&cylinder_with_lying_side(), &tol, &Default::default(), true);
    assert!(
        !given.is_ok(),
        "the curve-surface lie must still be reported"
    );
}
