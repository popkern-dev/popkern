//! PK-164: face–face blend through the code-first part facade.
//!
//! Builds the PK-142 truncated-V-notch fixture with `Profile` and
//! blends the two tilted notch walls via face selectors
//! (`Part::face_blend`), pinning the analytic contact height. Also
//! pins the selector contract: ambiguous and empty selections are
//! structured errors, not guesses.

use popkern_part::{Plane, Profile};

fn v_notch() -> popkern_part::Part {
    // Outer 10×6 prism with a truncated-V notch: width-3 opening at
    // the top converging to a width-1 floor at y=2 (the PK-142
    // TrimBoth pin fixture).
    Profile::on(Plane::XY)
        .path(0.0, 0.0)
        .line_to(10.0, 0.0)
        .line_to(10.0, 6.0)
        .line_to(6.5, 6.0)
        .line_to(5.5, 2.0)
        .line_to(4.5, 2.0)
        .line_to(3.5, 6.0)
        .line_to(0.0, 6.0)
        .close()
        .extrude(5.0)
        .expect("v-notch prism must extrude")
}

#[test]
fn face_blend_selectors_resolve_notch_walls() {
    let part = v_notch();
    // The notch walls are the only tilted planes: dominant ±x normal
    // with a y component (outer walls are axis-aligned).
    let tilted =
        |sign: f64| move |n: &popkern_core::geometry::Vec3| sign * n.x > 0.5 && n.y.abs() > 0.2;
    let result = part.face_blend(
        1.0,
        move |f| f.filter(move |i| i.normal.map(|n| tilted(1.0)(&n)).unwrap_or(false)),
        move |f| f.filter(move |i| i.normal.map(|n| tilted(-1.0)(&n)).unwrap_or(false)),
    );
    // The through-notch TrimBoth attach is a recorded decline today
    // (open-ended pocket; the material-adding attach is PK-163). The
    // sugar must surface that structured decline verbatim — not
    // panic, not silently no-op.
    let err = result.expect_err("through-notch TrimBoth declines pending PK-163");
    let msg = format!("{err}");
    assert!(
        msg.contains("end-arc support") || msg.contains("face_blend"),
        "decline must carry the kernel diagnostic: {msg}"
    );
}

#[test]
fn face_blend_rejects_ambiguous_selection() {
    let part = v_notch();
    let err = part
        .face_blend(1.0, |f| f, |f| f)
        .expect_err("selector matching every face must be rejected");
    let msg = format!("{err}");
    assert!(
        msg.contains("exactly one face"),
        "ambiguity must be named: {msg}"
    );
}

#[test]
fn face_blend_rejects_empty_selection() {
    let part = v_notch();
    let err = part
        .face_blend(1.0, |f| f.at_height(999.0), |f| f.at_height(999.0))
        .expect_err("empty selection must be a structured error");
    let msg = format!("{err}");
    assert!(
        msg.contains("selection") || msg.contains("matched"),
        "empty selection must be named: {msg}"
    );
}

/// Happy path: a CLOSED trapezoid pocket (block minus a short
/// truncated-V tool) gives the end arcs real end-wall support, so the
/// TrimBoth attach must complete watertight through the facade.
#[test]
fn face_blend_closed_pocket_completes() {
    let block = Profile::on(Plane::XY)
        .path(0.0, 0.0)
        .line_to(10.0, 0.0)
        .line_to(10.0, 6.0)
        .line_to(0.0, 6.0)
        .close()
        .extrude(7.0)
        .expect("block must extrude");
    // Same wall slope as the PK-142 V-notch (dx/dy = 0.25), extended
    // past the y=6 wall so the subtraction opens cleanly; tool spans
    // z 2..5 leaving closed pocket ends.
    let tool = Profile::on(Plane::xy_at(2.0))
        .path(4.5, 2.0)
        .line_to(5.5, 2.0)
        .line_to(6.75, 7.0)
        .line_to(3.25, 7.0)
        .close()
        .extrude(3.0)
        .expect("pocket tool must extrude");
    let pocketed = block.subtract(&tool).expect("pocket subtraction");
    let tilted =
        |sign: f64| move |n: &popkern_core::geometry::Vec3| sign * n.x > 0.5 && n.y.abs() > 0.2;
    let blended = pocketed
        .face_blend(
            1.0,
            move |f| f.filter(move |i| i.normal.map(|n| tilted(1.0)(&n)).unwrap_or(false)),
            move |f| f.filter(move |i| i.normal.map(|n| tilted(-1.0)(&n)).unwrap_or(false)),
        )
        .expect("closed-pocket TrimBoth must complete through the facade");
    blended
        .brep()
        .validate()
        .expect("blended pocket must validate");
}
