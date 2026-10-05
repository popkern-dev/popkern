//! Closed-source window pins (PK-673).

use super::SourceWindow;
use popkern_core::geometry::curve::{CurveEval, CurveGeometry};
use popkern_core::geometry::curve_on_surface::CurveOnSurface;
use popkern_core::geometry::{Point3, Vec3};

fn circle_source(closed: bool) -> CurveOnSurface {
    let circle = CurveGeometry::Circle {
        center: Point3::new(0.0, 0.0, 2.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        radius: 1.25,
    }
    .to_nurbs();
    CurveOnSurface {
        curve_uv1: circle.clone(),
        curve_uv2: circle.clone(),
        curve_3d: circle,
        is_closed: closed,
    }
}

#[test]
fn open_sources_pass_through_unchanged() {
    let source = circle_source(false);
    let w = SourceWindow::of(&source);
    assert_eq!(w.wrap(1.7), 1.7);
    assert_eq!(w.lift_near(0.1, 5.0), 0.1);
    assert_eq!(w.unwrap_piece(0.76, 0.01, 0.9), (0.76, 0.01));
    let (t0, _) = source.curve_3d.parameter_range();
    assert!(
        w.seam_crossing_subcurve(&source.curve_3d, t0 + 0.1, t0 + 0.3)
            .is_none()
    );
}

/// The M3 nut's crest arc: kept window 0.761 → 1.009 on a [0, 1]
/// loop, so its ends project to 0.761 and 0.009 and its midpoint to
/// ~0.885. Unwrapped, the traversal runs through the seam.
#[test]
fn seam_crossing_piece_unwraps_through_its_midpoint() {
    let source = circle_source(true);
    let (s0, s1) = source.curve_3d.parameter_range();
    let p = s1 - s0;
    let at = |f: f64| s0 + f * p;
    let w = SourceWindow::of(&source);
    // This end at 0.761, the other end at 0.009 through the seam.
    let (a, b) = w.unwrap_piece(at(0.761), at(0.009), at(0.885));
    assert_eq!(a, at(0.761));
    assert!((b - at(1.009)).abs() < 1e-12, "far end lifted to {b}");
    // Same ends, midpoint on the short arc: no seam, raw.
    assert_eq!(
        w.unwrap_piece(at(0.761), at(0.009), at(0.4)),
        (at(0.761), at(0.009))
    );
    // The junction parameter lifts next to the end it corrects.
    assert!((w.lift_near(at(0.0099), at(1.009)) - at(1.0099)).abs() < 1e-12);
    assert!((w.wrap(at(1.25)) - at(0.25)).abs() < 1e-12);
}

#[test]
fn seam_crossing_subcurve_is_the_source_arc() {
    let source = circle_source(true);
    let (s0, s1) = source.curve_3d.parameter_range();
    let p = s1 - s0;
    let w = SourceWindow::of(&source);
    let (lo, hi) = (s0 + 0.761 * p, s0 + 1.009 * p);
    let arc = w
        .seam_crossing_subcurve(&source.curve_3d, lo, hi)
        .expect("window crosses the seam");
    let (a0, a1) = arc.parameter_range();
    assert!((a0 - lo).abs() < 1e-12 && (a1 - hi).abs() < 1e-12);
    for k in 0..=64 {
        let t = lo + (hi - lo) * k as f64 / 64.0;
        let on_arc = arc.point_at(t);
        let on_source = source.curve_3d.point_at(w.wrap(t));
        assert!(
            on_arc.distance_to(&on_source) < 1e-12,
            "arc leaves its source at t = {t}"
        );
    }
    // A window that stays inside the domain is not a seam crossing.
    assert!(
        w.seam_crossing_subcurve(&source.curve_3d, s0 + 0.2 * p, s0 + 0.4 * p)
            .is_none()
    );
}
