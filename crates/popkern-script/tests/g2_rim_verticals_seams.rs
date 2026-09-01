//! PK-398 — the rim+verticals G2 class measures ACHIEVED G2 on the
//! whole-body instrument: every gated seam inside the commercial
//! G0/G1/G2 gates, with the cap↔strip seams orders of magnitude inside.
//!
//! History: the served cap↔strip seam class read 6.85e-2 (3.9°) on
//! every fixture of this class (PK-392). Its source was the setback
//! fan's rim row — the 4-point cubic interpolant of the strip's
//! section, whose sagitta (7.775e-4 at r = 0.1) the rim reconcile
//! could only land at four free columns; the along-seam slope of the
//! residual was the measured tangent jump. Beneath it sat a
//! transverse-curvature class of 1.4 (marches started from off-surface
//! control points) and a corner-block class of 1.4e-1 (the two pieces'
//! independent corner claims averaged into normal curvature). PK-398
//! fits the section as one C²-at-m degree-6 spline per arc, samples
//! the cross fields at true rim points, and pins both pieces at every
//! junction to shared Hermite data.

use popkern_part::prelude::*;
use popkern_script::Plane;
use popkern_script::prelude::*;

fn rim_verticals(w: f64, d: f64, h: f64, r: f64) -> Model {
    let mut m = Model::new();
    let plate = m
        .profile(Plane::XY)
        .rect(w, d)
        .sketch("plate")
        .expect("sketch")
        .extrude(&mut m, End::Blind(h.into()))
        .expect("extrude");
    plate
        .fillet_continuity(
            &mut m,
            edges()
                .straight()
                .at_height(h)
                .or(edges().straight().direction(Dir::Z)),
            r,
            ContinuityKind::G2,
            "rounds",
        )
        .expect("G2 fillet serves");
    m
}

/// Acceptance 593 (exceeded): cube and plate measure achieved G2; the
/// cap↔strip class that held them at G1 is pinned at its new floor.
#[test]
fn rim_verticals_g2_measures_achieved_g2_on_every_gated_seam() {
    for (w, d, h, r) in [
        (1.0, 1.0, 1.0, 0.1),
        (2.4, 1.6, 0.4, 0.06),
        (2.4, 1.6, 0.4, 0.08),
        (2.4, 1.6, 0.4, 0.1),
    ] {
        let mut m = rim_verticals(w, d, h, r);
        m.assert_watertight().expect("watertight");
        let report = m
            .evidence(&m.feature("rounds").expect("label"))
            .expect("evidence");
        let blend = report.blend.as_ref().expect("blend");
        let c = blend.continuity.as_ref().expect("continuity");
        let worst_gated = c
            .seams
            .iter()
            .filter(|s| s.gated)
            .map(|s| {
                (
                    s.worst_g1_normal,
                    s.kind,
                    s.worst_g0_gap,
                    s.worst_g2_residual,
                )
            })
            .fold(None, |acc: Option<(f64, _, f64, f64)>, s| match acc {
                Some(a) if a.0 >= s.0 => Some(a),
                _ => Some(s),
            });
        println!(
            "rim+verticals {w}x{d}x{h} r={r}: achieved {:?} g0={:.3e} g1={:.3e} g2={:.3e} worst gated seam {:?}",
            c.achieved, c.worst_g0_gap, c.worst_g1_normal, c.worst_g2_residual, worst_gated
        );
        assert!(c.measured(), "{w}x{d}x{h} r={r}: seams measured");
        assert_eq!(
            c.achieved,
            Continuity::G2,
            "{w}x{d}x{h} r={r}: achieved level (worst gated seam {:?})",
            worst_gated
        );
        // Floors of the retired classes, pinned an order of magnitude
        // above the measured values (g1 1.6e-7, g2 1.4e-4, g0 5e-11).
        assert!(
            c.worst_g1_normal < 1e-5,
            "{w}x{d}x{h} r={r}: worst gated tangent residual {:.3e} (the cubic-sagitta class was 6.85e-2)",
            c.worst_g1_normal,
        );
        assert!(
            c.worst_g2_residual < 2e-3,
            "{w}x{d}x{h} r={r}: worst gated curvature jump {:.3e} (the off-surface-march class was 1.4)",
            c.worst_g2_residual,
        );
        assert!(
            c.worst_g0_gap < 1e-9,
            "{w}x{d}x{h} r={r}: worst seam gap {:.3e} (the cubic-sagitta class was 4.2e-4)",
            c.worst_g0_gap
        );
    }
}
