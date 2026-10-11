//! `closest_on_piece` on open and closed pieces (PK-777).

use std::f64::consts::TAU;

use popkern_core::geometry::curve::{CurveEval, CurveGeometry, CurveOnEdge};
use popkern_core::geometry::{Point3, Vec3};

use super::closest_on_piece;

fn rim(window: (f64, f64)) -> CurveOnEdge {
    CurveOnEdge {
        geometry: CurveGeometry::Circle {
            center: Point3::new(-6.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            radius: 4.5,
        },
        t_start: window.0,
        t_end: window.1,
    }
}

#[test]
fn a_point_just_short_of_a_closed_rims_closing_end_lies_on_it() {
    // Measured: a point just short of a [0, 2π] rim's closing end read
    // 0.30 away — the opening end's sample tied the closing end's (the
    // tie breaks on the last ulp of the two evaluations) and the polish
    // clamped at the opening end. Windows at many phases and points on
    // both sides of the joint hit both tie breaks.
    let mut off = Vec::new();
    for k in 0..48 {
        let a = -3.0 + 0.137 * k as f64;
        for window in [(a, a + TAU), (a + TAU, a)] {
            let piece = rim(window);
            for dt in [-0.13, -0.07, -0.01, 0.01, 0.07, 0.13] {
                let p = piece.geometry.point_at(a + dt);
                let (_, d) = closest_on_piece(&piece, &p);
                if d > 1e-9 {
                    off.push((window, dt, d));
                }
            }
        }
    }
    assert!(
        off.is_empty(),
        "{} off: {:?}",
        off.len(),
        &off[..off.len().min(4)]
    );
}

#[test]
fn an_open_piece_still_measures_to_its_own_window() {
    // A point on the carrier past the window's end is NOT on the piece.
    let piece = rim((0.5, 2.0));
    let p = piece.geometry.point_at(3.0);
    let (t, d) = closest_on_piece(&piece, &p);
    assert!((t - 2.0).abs() < 1e-9 && d > 1.0, "t {t} d {d}");
}
