//! PK-781: conics crossing rods, against the analytic crossing count.

use std::f64::consts::TAU;

use super::super::intersect_surface_curve_general;
use super::{curve_period, lift_into_window};
use crate::geometry::curve::{CurveEval, CurveGeometry};
use crate::geometry::intersection::SurfaceCurveIntersection;
use crate::geometry::surface::SurfaceGeometry;
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;

const R: f64 = 5.0;

fn z() -> Vec3 {
    Vec3::new(0.0, 0.0, 1.0)
}

fn circle() -> CurveGeometry {
    CurveGeometry::Circle {
        center: Point3::new(0.0, 0.0, 0.0),
        normal: z(),
        radius: R,
    }
}

fn arc(stored: (f64, f64)) -> CurveGeometry {
    CurveGeometry::Arc {
        center: Point3::new(0.0, 0.0, 0.0),
        normal: z(),
        radius: R,
        start_angle: stored.0,
        end_angle: stored.1,
    }
}

/// A rod of radius `r` standing on the circle's radial line at angle
/// `phi`, its axis `R + d` from the centre. The circle's angle-0 point
/// comes from the kernel's own basis, not +x.
fn rod(phi: f64, d: f64, r: f64) -> SurfaceGeometry {
    let p = circle().point_at(phi);
    let s = (R + d) / R;
    SurfaceGeometry::Cylinder {
        origin: Point3::new(p.x * s, p.y * s, -1.0),
        axis: z(),
        radius: r,
    }
}

/// The crossing angles of the circle with the rod's section circle.
fn analytic_roots(phi: f64, d: f64, r: f64) -> Vec<f64> {
    let dist = R + d;
    let c = (R * R + dist * dist - r * r) / (2.0 * R * dist);
    if c.abs() >= 1.0 {
        return Vec::new();
    }
    let alpha = c.acos();
    vec![phi - alpha, phi + alpha]
}

fn solve(surface: &SurfaceGeometry, curve: &CurveGeometry, window: (f64, f64)) -> Vec<f64> {
    intersect_surface_curve_general(
        surface,
        (0.0, TAU),
        (-2.0, 2.0),
        curve,
        window,
        &Tolerance::default(),
        false,
    )
    .iter()
    .filter_map(|i| match i {
        SurfaceCurveIntersection::Point { t, .. } => Some(*t),
        SurfaceCurveIntersection::Overlap { .. } => None,
    })
    .collect()
}

/// Every analytic root inside `window` (lifted by whole turns) is found
/// at its own parameter, and nothing else is.
fn matches(found: &[f64], phi: f64, d: f64, r: f64, window: (f64, f64)) -> bool {
    let want: Vec<f64> = analytic_roots(phi, d, r)
        .into_iter()
        .filter_map(|t| lift_into_window(t, window, 0.0, Some(TAU)))
        .collect();
    found.len() == want.len()
        && want
            .iter()
            .all(|w| found.iter().any(|f| (f - w).abs() < 1e-7))
}

#[test]
fn an_arc_finds_both_crossings_of_a_rod_on_its_own_and_a_shifted_window() {
    for (stored, window) in [
        ((-0.5, 0.5), (-0.5, 0.5)),
        // Stored one turn above the query window.
        ((TAU - 0.5, TAU + 0.5), (-0.5, 0.5)),
    ] {
        let found = solve(&rod(0.0, 0.0, 1.0), &arc(stored), window);
        assert!(
            matches(&found, 0.0, 0.0, 1.0, window),
            "stored {stored:?} window {window:?}: {found:?}"
        );
    }
}

#[test]
fn a_circle_finds_a_crossing_just_below_the_rods_seam() {
    // Measured: this pair's second crossing polished to cylinder
    // u = −0.096 and was range-filtered out (the u-period lift).
    let found = solve(&rod(0.3, 0.2997, 1.0), &circle(), (0.0, TAU));
    assert!(matches(&found, 0.3, 0.2997, 1.0, (0.0, TAU)), "{found:?}");
}

#[test]
fn a_circle_isolates_two_crossings_of_one_segment_cell_pair() {
    // Measured: with the period lifts but no pair isolation, these
    // offsets returned one crossing of two.
    let mut missed = Vec::new();
    for d in [
        -0.2997, -0.2498, -0.0499, 0.05, 0.0999, 0.2498, 0.3496, 0.5495,
    ] {
        let found = solve(&rod(0.3, d, 1.0), &circle(), (0.0, TAU));
        if !matches(&found, 0.3, d, 1.0, (0.0, TAU)) {
            missed.push((d, found));
        }
    }
    assert!(missed.is_empty(), "{missed:?}");
}

/// The PK-781 sweep: arcs (own and one-turn-shifted windows) and
/// circles (full and partial windows) through rods of three radii,
/// offset through tangency to 0.999 r. At 0.99 r and beyond, a 0.02 rod
/// is pierced only 2e-4..2e-5 deep, and a 16-sample arc-segment box
/// missed the arc's apex by more than that: all four arc misses were
/// that cull (PK-785; the segment box is now exact).
#[test]
fn conics_through_rods_find_every_analytic_crossing() {
    let mut missed = Vec::new();
    let mut cases = 0;
    for r in [1.0, 0.2, 0.02] {
        let offsets = (0..=20)
            .map(|k| -0.95 + 1.9 * k as f64 / 20.0)
            .chain([-0.999, -0.99, 0.99, 0.999]);
        for f in offsets {
            let d = r * f;
            for phi in [0.0, 0.3, 1.0, 2.5, 4.0] {
                let near = (phi - 0.5, phi + 0.5);
                for (name, curve, window) in [
                    ("arc", arc(near), near),
                    ("arc-shift", arc((near.0 + TAU, near.1 + TAU)), near),
                    ("circle", circle(), (0.0, TAU)),
                    ("circle-window", circle(), near),
                ] {
                    cases += 1;
                    let found = solve(&rod(phi, d, r), &curve, window);
                    if !matches(&found, phi, d, r, window) {
                        missed.push(format!("{name} r={r} d={d:.4} phi={phi}: {found:?}"));
                    }
                }
            }
        }
    }
    assert!(
        missed.is_empty(),
        "{} of {cases}: {missed:#?}",
        missed.len()
    );
}

#[test]
fn lifts_move_by_whole_periods_only_when_periodic() {
    assert_eq!(lift_into_window(-0.1, (0.0, TAU), 0.0, None), None);
    let lifted = lift_into_window(-0.1, (0.0, TAU), 0.0, Some(TAU));
    assert!(lifted.is_some_and(|t| (t - (TAU - 0.1)).abs() < 1e-12));
    let down = lift_into_window(TAU + 0.2, (-0.5, 0.5), 0.0, Some(TAU));
    assert!(down.is_some_and(|t| (t - 0.2).abs() < 1e-12));
    // Reversed windows are windows too.
    assert_eq!(lift_into_window(0.2, (0.5, -0.5), 0.0, None), Some(0.2));
    assert_eq!(lift_into_window(1.0, (-0.5, 0.5), 0.0, Some(TAU)), None);
    assert_eq!(curve_period(&circle()), Some(TAU));
    let line = CurveGeometry::Line {
        origin: Point3::new(0.0, 0.0, 0.0),
        direction: z(),
    };
    assert_eq!(curve_period(&line), None);
}
