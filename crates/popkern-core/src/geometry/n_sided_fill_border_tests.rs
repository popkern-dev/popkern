//! Does one bi-2 border per (arc + bridge) PAIR, or one per EDGE,
//! track a corner rim better?
//!
//! The PnS adapter builds one bi-2 position row spanning a sector's
//! arc *concatenated with* its bridge. A degree-2 polynomial cannot
//! represent a circular arc at all, let alone an arc joined to a
//! line, and the committed caps measure 0.0796 from their rim against
//! a 5.000e-3 acceptance bar.
//!
//! The proposed alternative is one PnS sector per rim EDGE (valence
//! 2N instead of N), so each border spans a single arc or a single
//! line. That change needs a different topology poke, so it is worth
//! knowing whether the geometry actually improves before paying for
//! it. This measures exactly that, on exact geometry, with no
//! pipeline involved.

use crate::geometry::{Point3, Vec3};

/// Sample a quarter arc of `radius` in the xz-plane, then a straight
/// run of `line_len` continuing from its end.
fn arc_then_line(radius: f64, line_len: f64, samples: usize) -> Vec<Point3> {
    let mut pts = Vec::new();
    for k in 0..=samples {
        let a = std::f64::consts::FRAC_PI_2 * k as f64 / samples as f64;
        pts.push(Point3::new(
            radius * f64::sin(a),
            0.0,
            radius * (1.0 - f64::cos(a)),
        ));
    }
    let end = *pts.last().expect("arc sampled");
    for k in 1..=samples {
        let t = line_len * k as f64 / samples as f64;
        pts.push(Point3::new(end.x, 0.0, end.z + t));
    }
    pts
}

/// Max distance from `pts` to the bi-2 Bezier through the first,
/// middle and last of them, at matched fractions.
fn bi2_deviation(pts: &[Point3]) -> f64 {
    let (a, c) = (pts[0], pts[pts.len() - 1]);
    let m = pts[pts.len() / 2];
    let b = Point3::new(
        2.0 * m.x - 0.5 * (a.x + c.x),
        2.0 * m.y - 0.5 * (a.y + c.y),
        2.0 * m.z - 0.5 * (a.z + c.z),
    );
    let mut worst = 0.0_f64;
    for (i, p) in pts.iter().enumerate() {
        let t = i as f64 / (pts.len() - 1) as f64;
        let (u, v, w) = ((1.0 - t) * (1.0 - t), 2.0 * t * (1.0 - t), t * t);
        let q = Point3::new(
            u * a.x + v * b.x + w * c.x,
            u * a.y + v * b.y + w * c.y,
            u * a.z + v * b.z + w * c.z,
        );
        worst = worst.max(q.distance_to(p));
    }
    worst
}

/// FALSIFICATION. Splitting the border per rim edge is WORSE, not
/// better — measured before anyone paid for the topology rework it
/// would have required.
///
/// The reasoning that motivated it was that no quadratic represents
/// an arc concatenated with a line, so one border per edge should
/// track the rim more closely. Measured, the pair spans the rim at
/// 1.119e-2 while the worst per-edge border is 1.715e-2: splitting
/// makes it 0.7x as good.
///
/// The reason is the real constraint: a degree-2 polynomial cannot
/// represent a 90-degree circular ARC either. Isolated, the arc is
/// the worse half. Pairing it with an exactly-representable straight
/// run lets the shared middle coefficient sit somewhere that serves
/// both. So the border problem is DEGREE, not partitioning — the fix
/// has to raise the border's representational power (higher degree or
/// rational), and re-partitioning cannot help.
#[test]
fn splitting_the_border_per_edge_is_worse_not_better() {
    // A 0.05 fillet radius with a bridge of comparable length — the
    // shape of a real low-valence corner sector.
    let radius = 0.05;
    let line_len = 0.05;
    let arc = arc_then_line(radius, 0.0, 24);
    let line = {
        let end = *arc.last().expect("arc");
        (0..=24)
            .map(|k| {
                let t = line_len * f64::from(k) / 24.0;
                Point3::new(end.x, 0.0, end.z + t)
            })
            .collect::<Vec<_>>()
    };
    let paired = arc_then_line(radius, line_len, 24);

    let paired_dev = bi2_deviation(&paired);
    let per_edge_dev = bi2_deviation(&arc).max(bi2_deviation(&line));

    println!(
        "[pk192-border] one-per-pair={paired_dev:.4e}  one-per-edge={per_edge_dev:.4e}  \
         ratio={:.1}x",
        paired_dev / per_edge_dev.max(f64::MIN_POSITIVE)
    );

    assert!(
        per_edge_dev > paired_dev,
        "measured falsification: per-edge borders were WORSE \
         (pair={paired_dev:.4e}, edge={per_edge_dev:.4e}). If this ever \
         reverses, the partitioning idea is back on the table and this \
         gate should say so."
    );
    // And the deeper point: even the arc ALONE is far over the
    // 5.000e-3 acceptance bar, so no partitioning of degree-2 borders
    // can reach it.
    assert!(
        bi2_deviation(&arc) > 5.0e-3,
        "a degree-2 border over a quarter arc must be shown to miss the \
         acceptance bar on its own, or the conclusion above is unfounded"
    );
}

/// A straight run is exactly representable, so the whole error of the
/// paired form comes from forcing one quadratic across the arc/line
/// junction.
#[test]
fn a_line_alone_is_exact() {
    let line: Vec<Point3> = (0..=24)
        .map(|k| Point3::new(0.1 * f64::from(k) / 24.0, 0.0, 0.0))
        .collect();
    assert!(bi2_deviation(&line) < 1e-15);
}

/// Guard the fixture: the arc must actually be curved, or the
/// comparison above would be measuring nothing.
#[test]
fn the_arc_fixture_is_genuinely_curved() {
    let arc = arc_then_line(0.05, 0.0, 24);
    let chord = arc[0].distance_to(&arc[arc.len() - 1]);
    let mid = arc[arc.len() / 2];
    let straight = Point3::new(
        0.5 * (arc[0].x + arc[arc.len() - 1].x),
        0.5 * (arc[0].y + arc[arc.len() - 1].y),
        0.5 * (arc[0].z + arc[arc.len() - 1].z),
    );
    assert!(
        mid.distance_to(&straight) > 0.1 * chord,
        "the arc fixture must bow well away from its chord"
    );
    let _ = Vec3::new(0.0, 0.0, 0.0);
}

/// Max deviation from a degree-4 Bezier least-squares fitted to
/// `pts`, ends interpolated.
///
/// The cap's boundary row is ALREADY degree 4 — it is just derived by
/// elevating a degree-2 row, which throws away the freedom the
/// degree-4 form has. This measures what that freedom is worth,
/// because if it does not clear the acceptance bar then raising the
/// border's representational power is not the answer either.
fn bi4_fit_deviation(pts: &[Point3]) -> f64 {
    let (p0, p4) = (pts[0], pts[pts.len() - 1]);
    // Bernstein degree 4; solve for the three interior coefficients by
    // normal equations on each axis independently.
    let basis = |t: f64| {
        let s = 1.0 - t;
        [
            s * s * s * s,
            4.0 * t * s * s * s,
            6.0 * t * t * s * s,
            4.0 * t * t * t * s,
            t * t * t * t,
        ]
    };
    let mut ata = [[0.0_f64; 3]; 3];
    let mut atb = [[0.0_f64; 3]; 3]; // [coef][axis]
    for (i, p) in pts.iter().enumerate() {
        let t = i as f64 / (pts.len() - 1) as f64;
        let b = basis(t);
        let residual = [
            p.x - b[0] * p0.x - b[4] * p4.x,
            p.y - b[0] * p0.y - b[4] * p4.y,
            p.z - b[0] * p0.z - b[4] * p4.z,
        ];
        for r in 0..3 {
            for c in 0..3 {
                ata[r][c] += b[r + 1] * b[c + 1];
            }
            for axis in 0..3 {
                atb[r][axis] += b[r + 1] * residual[axis];
            }
        }
    }
    // Gaussian elimination on the 3x3, three right-hand sides.
    let mut m = ata;
    let mut rhs = atb;
    for col in 0..3 {
        let pivot = (col..3)
            .max_by(|a, b| m[*a][col].abs().partial_cmp(&m[*b][col].abs()).unwrap())
            .expect("pivot exists");
        m.swap(col, pivot);
        rhs.swap(col, pivot);
        let d = m[col][col];
        if d.abs() < 1e-18 {
            return f64::INFINITY;
        }
        for r in 0..3 {
            if r == col {
                continue;
            }
            let f = m[r][col] / d;
            let pivot_row = m[col];
            for (target, &pv) in m[r].iter_mut().zip(&pivot_row) {
                *target -= f * pv;
            }
            let pivot_rhs = rhs[col];
            for (target, &pv) in rhs[r].iter_mut().zip(&pivot_rhs) {
                *target -= f * pv;
            }
        }
    }
    let mut interior = [[0.0_f64; 3]; 3];
    for r in 0..3 {
        for axis in 0..3 {
            interior[r][axis] = rhs[r][axis] / m[r][r];
        }
    }
    let cps = [
        [p0.x, p0.y, p0.z],
        interior[0],
        interior[1],
        interior[2],
        [p4.x, p4.y, p4.z],
    ];
    let mut worst = 0.0_f64;
    for (i, p) in pts.iter().enumerate() {
        let t = i as f64 / (pts.len() - 1) as f64;
        let b = basis(t);
        let mut q = [0.0_f64; 3];
        for k in 0..5 {
            for axis in 0..3 {
                q[axis] += b[k] * cps[k][axis];
            }
        }
        worst = worst.max(Point3::new(q[0], q[1], q[2]).distance_to(p));
    }
    worst
}

/// Is raising the border's degree enough to clear the acceptance bar?
///
/// The cap boundary is already degree 4; today it is an ELEVATED
/// degree-2 row, which wastes that. If fitting degree 4 directly
/// clears 5.000e-3 then the fix is contained — no new topology, no
/// rational forms, just stop throwing the freedom away.
#[test]
fn fitting_the_border_at_its_true_degree_clears_the_acceptance_bar() {
    let paired = arc_then_line(0.05, 0.05, 48);
    let elevated = bi2_deviation(&paired);
    let fitted = bi4_fit_deviation(&paired);
    println!(
        "[pk192-degree] bi2-elevated={elevated:.4e}  bi4-fitted={fitted:.4e}  \
         improvement={:.0}x",
        elevated / fitted.max(f64::MIN_POSITIVE)
    );
    assert!(
        fitted < elevated,
        "degree-4 fitting must beat the elevated degree-2 row"
    );
    assert!(
        fitted < 5.0e-3,
        "MEASURED VERDICT: degree-4 fitting reaches {fitted:.4e} against the \
         5.000e-3 acceptance bar. If this fails, raising the border degree is \
         NOT sufficient either and the next idea must come from elsewhere."
    );
}
