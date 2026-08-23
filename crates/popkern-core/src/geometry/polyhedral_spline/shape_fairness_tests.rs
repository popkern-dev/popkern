//! Interior-fairness pins for the committed macro cap tables (PK-318).
//!
//! The 2x2 bi-3 macro caps (star/ngon, n >= 6) are C1 across each sector's
//! internal junctions by construction, so no continuity gate sees their
//! interior quality. What isophotes DO see is the jump of the normal's
//! derivative across those junctions: the committed tables use the guided
//! (extension-agreement) interiors, which keep that jump at reference-package
//! level, while the paper-literal F3 interiors measured 5-10x larger
//! (docs/reports/pk318-macro-shape-quality.md). These tests pin the guided
//! level on a curved canonical net so a table regression (e.g. regenerating
//! without the shape gate) fails here instead of shipping silently.
//!
//! The nets are the canonical dome battery from
//! `tools/reference/pns-tables/shape_quality.py`; the metric is the maximal
//! cross-seam jump of the finite-differenced normal derivative along both
//! internal junction lines of every sector. Bounds are ~2x the measured
//! committed-table values (see the report for the measured battery).

use super::*;
use crate::geometry::{Point3, Vec3};

fn dome_z(x: f64, y: f64) -> f64 {
    0.35 * (x * x + y * y)
}

/// CC-net stencil for the star cap: columns q1^0, q2^0, ..., q1^{n-1},
/// q2^{n-1}, E — the layout of the committed `eopSct{n}` tables.
fn star_dome_net(n: usize) -> Vec<Point3> {
    let alpha = 2.0 * std::f64::consts::PI / n as f64;
    let mut pts = Vec::with_capacity(2 * n + 1);
    for s in 0..n {
        let a1 = s as f64 * alpha;
        let a2 = a1 + alpha / 2.0;
        let (x1, y1) = (a1.cos(), a1.sin());
        let (x2, y2) = (1.4 * a2.cos(), 1.4 * a2.sin());
        pts.push(Point3::new(x1, y1, dome_z(x1, y1)));
        pts.push(Point3::new(x2, y2, dome_z(x2, y2)));
    }
    pts.push(Point3::new(0.0, 0.0, 0.0));
    pts
}

/// DS-net stencil for the ngon cap: sector-major groups d0 (diagonal),
/// d1, d2 (edge-adjacent), d3 (central n-gon vertex) — the layout of the
/// committed `ngonSct{n}` tables.
fn ngon_dome_net(n: usize) -> Vec<Point3> {
    let alpha = 2.0 * std::f64::consts::PI / n as f64;
    let delta = alpha / 3.0;
    let mut pts = Vec::with_capacity(4 * n);
    for k in 0..n {
        let th = k as f64 * alpha;
        for (r, a) in [(2.4, th), (1.9, th + delta), (1.9, th - delta), (1.0, th)] {
            let (x, y) = (r * a.cos(), r * a.sin());
            pts.push(Point3::new(x, y, dome_z(x, y)));
        }
    }
    pts
}

fn unit_normal(patch: &PolyhedralPatch, u: f64, v: f64) -> Vec3 {
    let n = patch.partial_u(u, v).cross(&patch.partial_v(u, v));
    let len = n.length();
    assert!(len > 1e-12, "degenerate normal at ({u},{v})");
    n * (1.0 / len)
}

/// Max cross-seam jump of the FD normal derivative over the internal C1
/// junctions of a macro cap (patch index 4k + I + 2J; u-seams I=0|1,
/// v-seams J=0|1).
fn max_internal_normal_derivative_jump(spline: &PolyhedralSpline, n: usize) -> f64 {
    assert_eq!(spline.patches.len(), 4 * n, "macro cap expected");
    let h = 1e-4;
    let samples = 17;
    let mut worst = 0.0f64;
    let mut record = |jump: f64| {
        if jump > worst {
            worst = jump;
        }
    };
    for k in 0..n {
        for half in 0..2 {
            // u-seam: patch (I=0) at u=1 against patch (I=1) at u=0.
            let left = &spline.patches[4 * k + 2 * half];
            let right = &spline.patches[4 * k + 1 + 2 * half];
            for s in 0..samples {
                let t = s as f64 / (samples - 1) as f64;
                let dl = (unit_normal(left, 1.0, t) - unit_normal(left, 1.0 - h, t)) * (1.0 / h);
                let dr = (unit_normal(right, h, t) - unit_normal(right, 0.0, t)) * (1.0 / h);
                record((dl - dr).length());
            }
            // v-seam: patch (J=0) at v=1 against patch (J=1) at v=0.
            let bottom = &spline.patches[4 * k + half];
            let top = &spline.patches[4 * k + 2 + half];
            for s in 0..samples {
                let t = s as f64 / (samples - 1) as f64;
                let dl =
                    (unit_normal(bottom, t, 1.0) - unit_normal(bottom, t, 1.0 - h)) * (1.0 / h);
                let dr = (unit_normal(top, t, h) - unit_normal(top, t, 0.0)) * (1.0 / h);
                record((dl - dr).length());
            }
        }
    }
    worst
}

#[test]
fn macro_star_caps_keep_guided_interior_fairness() {
    // Measured committed-table values (guided interiors): n=6 0.0101,
    // n=7 0.0211, n=8 0.0292; F3-class interiors measure 5-10x larger.
    for (n, bound) in [(6usize, 0.025), (7, 0.05), (8, 0.07)] {
        let spline = construct_extraordinary_star(n, &star_dome_net(n)).unwrap();
        let jump = max_internal_normal_derivative_jump(&spline, n);
        assert!(
            jump < bound,
            "star n={n}: internal-junction normal-derivative jump {jump} \
             exceeds the guided-interior bound {bound} — did the committed \
             macro tables regress to F3-class interiors?"
        );
    }
}

#[test]
fn macro_ngon_caps_keep_guided_interior_fairness() {
    // Measured committed-table values (guided interiors): n=6 0.0208,
    // n=7 0.0441, n=8 0.0611; F3-class interiors measure 5-10x larger.
    for (n, bound) in [(6usize, 0.05), (7, 0.09), (8, 0.13)] {
        let spline = construct_ngon(n, &ngon_dome_net(n)).unwrap();
        let jump = max_internal_normal_derivative_jump(&spline, n);
        assert!(
            jump < bound,
            "ngon n={n}: internal-junction normal-derivative jump {jump} \
             exceeds the guided-interior bound {bound} — did the committed \
             macro tables regress to F3-class interiors?"
        );
    }
}
