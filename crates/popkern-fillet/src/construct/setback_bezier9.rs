//! Degree-9 Bézier row utilities for the setback fan (de Casteljau
//! evaluation, nearest-point, subdivision) — split from
//! `setback_vertex_topology.rs` for the 500-line/frozen-ratchet
//! structural gate.

use popkern_core::geometry::Point3;

pub(crate) fn bezier9_eval(cps: &[Point3; 10], t: f64) -> Point3 {
    let mut work = *cps;
    for level in (1..10).rev() {
        for i in 0..level {
            work[i] = Point3::new(
                work[i].x + (work[i + 1].x - work[i].x) * t,
                work[i].y + (work[i + 1].y - work[i].y) * t,
                work[i].z + (work[i + 1].z - work[i].z) * t,
            );
        }
    }
    work[0]
}

/// Nearest parameter on a degree-9 Bezier (dense scan + ternary).
pub(crate) fn bezier9_nearest(cps: &[Point3; 10], p: Point3) -> (f64, f64) {
    let mut best = (f64::INFINITY, 0.0_f64);
    for i in 0..=128 {
        let t = i as f64 / 128.0;
        let d = bezier9_eval(cps, t).distance_to(&p);
        if d < best.0 {
            best = (d, t);
        }
    }
    let (mut lo, mut hi) = (
        (best.1 - 1.0 / 128.0).max(0.0),
        (best.1 + 1.0 / 128.0).min(1.0),
    );
    for _ in 0..60 {
        let m1 = lo + (hi - lo) / 3.0;
        let m2 = hi - (hi - lo) / 3.0;
        if bezier9_eval(cps, m1).distance_to(&p) < bezier9_eval(cps, m2).distance_to(&p) {
            hi = m2;
        } else {
            lo = m1;
        }
    }
    let t = 0.5 * (lo + hi);
    (t, bezier9_eval(cps, t).distance_to(&p))
}

/// De Casteljau extraction of the sub-segment [ta, tb].
pub(crate) fn bezier9_extract(cps: &[Point3; 10], ta: f64, tb: f64) -> [Point3; 10] {
    // Split at ta, keep right; split that at (tb-ta)/(1-ta), keep left.
    let right = bezier9_split_right(cps, ta);
    let f = if (1.0 - ta).abs() < 1e-15 {
        0.0
    } else {
        (tb - ta) / (1.0 - ta)
    };
    bezier9_split_left(&right, f)
}

fn bezier9_split_right(cps: &[Point3; 10], t: f64) -> [Point3; 10] {
    // Right-segment control points are the terminal point of each de
    // Casteljau level: [C(t), b_1^8, ..., b_8^1, b_9].
    let mut work = *cps;
    let mut right = [cps[9]; 10];
    for level in (1..10).rev() {
        for i in 0..level {
            work[i] = Point3::new(
                work[i].x + (work[i + 1].x - work[i].x) * t,
                work[i].y + (work[i + 1].y - work[i].y) * t,
                work[i].z + (work[i + 1].z - work[i].z) * t,
            );
        }
        right[level - 1] = work[level - 1];
    }
    right[0] = work[0];
    right
}

fn bezier9_split_left(cps: &[Point3; 10], t: f64) -> [Point3; 10] {
    let mut work = *cps;
    let mut left = [cps[0]; 10];
    for (level, slot) in (1..10).rev().zip(1..10) {
        for i in 0..level {
            work[i] = Point3::new(
                work[i].x + (work[i + 1].x - work[i].x) * t,
                work[i].y + (work[i + 1].y - work[i].y) * t,
                work[i].z + (work[i + 1].z - work[i].z) * t,
            );
        }
        left[slot] = work[0];
    }
    left
}
