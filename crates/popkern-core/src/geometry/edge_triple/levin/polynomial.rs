//! Real-root solvers for low-degree polynomials.
//!
//! Used by:
//! - Stage 2b.32: characteristic-polynomial roots of the quadric pencil.
//! - Stage 2b.34: quartic from the parametric reduction.
//!
//! All roots are returned in ascending order with near-duplicates
//! merged. Tolerance constants are deliberately conservative (1e-12
//! to 1e-9) — Levin's algorithm is forgiving of small root errors
//! since the quartic is just an *intermediate* — the final curve
//! parameters are projected back onto the original quadrics, which
//! polishes residuals.
//!
//! These are **f64** implementations; stage 2b.36 lifts to `Q[√d]`
//! with exact discriminants.

const EPS_LEAD: f64 = 1e-14;
const EPS_DISC: f64 = 1e-12;
const EPS_DEDUP: f64 = 1e-9;

/// Real roots of `a·x² + b·x + c = 0`, sorted ascending.
pub(crate) fn solve_quadratic(coef: [f64; 3]) -> Vec<f64> {
    let [a, b, c] = coef;
    if a.abs() < EPS_LEAD {
        if b.abs() < EPS_LEAD {
            return Vec::new();
        }
        return vec![-c / b];
    }
    let disc = b * b - 4.0 * a * c;
    if disc < -EPS_DISC {
        return Vec::new();
    }
    let s = disc.max(0.0).sqrt();
    // Numerically stable form: avoid catastrophic cancellation.
    let q = -0.5 * (b + b.signum() * s);
    let mut r = if q.abs() < EPS_LEAD {
        vec![0.0, -b / a]
    } else {
        vec![q / a, c / q]
    };
    r.sort_by(|x, y| x.partial_cmp(y).unwrap_or(std::cmp::Ordering::Equal));
    r
}

/// Real roots of `a·x³ + b·x² + c·x + d = 0`, sorted ascending.
/// Uses the trigonometric form when the discriminant is positive
/// (three real roots) and Cardano otherwise.
pub(crate) fn solve_cubic(coef: [f64; 4]) -> Vec<f64> {
    let [a, b, c, d] = coef;
    if a.abs() < EPS_LEAD {
        return solve_quadratic([b, c, d]);
    }
    // Normalise + depress: x = y − b/(3a).
    let bn = b / a;
    let cn = c / a;
    let dn = d / a;
    let p = cn - bn * bn / 3.0;
    let q = 2.0 * bn.powi(3) / 27.0 - bn * cn / 3.0 + dn;
    // y³ + p·y + q = 0.
    let disc = -4.0 * p.powi(3) - 27.0 * q * q;
    let mut roots = Vec::new();
    if disc > EPS_DISC {
        // Three distinct real roots — trigonometric form.
        let r = (-p / 3.0).sqrt();
        let arg = ((-q / 2.0) / r.powi(3)).clamp(-1.0, 1.0);
        let theta = arg.acos();
        for k in 0..3 {
            let angle = (theta + std::f64::consts::TAU * (k as f64)) / 3.0;
            let y = 2.0 * r * angle.cos();
            roots.push(y - bn / 3.0);
        }
    } else if disc.abs() <= EPS_DISC {
        // Multiple root.
        if p.abs() < EPS_LEAD {
            roots.push(-bn / 3.0);
        } else {
            let y1 = 3.0 * q / p;
            let y2 = -y1 / 2.0;
            roots.push(y1 - bn / 3.0);
            roots.push(y2 - bn / 3.0);
            roots.push(y2 - bn / 3.0);
        }
    } else {
        // One real, two complex conjugates — Cardano.
        let inner = q * q / 4.0 + p.powi(3) / 27.0;
        let s = inner.sqrt();
        let u = cbrt_signed(-q / 2.0 + s);
        let v = cbrt_signed(-q / 2.0 - s);
        roots.push(u + v - bn / 3.0);
    }
    roots.sort_by(|x, y| x.partial_cmp(y).unwrap_or(std::cmp::Ordering::Equal));
    roots.dedup_by(|x, y| (*x - *y).abs() < EPS_DEDUP);
    roots
}

/// Real roots of `c[0]·x⁴ + c[1]·x³ + c[2]·x² + c[3]·x + c[4] = 0`,
/// sorted ascending. Uses Ferrari's method via the resolvent cubic.
pub(crate) fn solve_quartic(coef: [f64; 5]) -> Vec<f64> {
    let [a4, a3, a2, a1, a0] = coef;
    if a4.abs() < EPS_LEAD {
        return solve_cubic([a3, a2, a1, a0]);
    }
    let b = a3 / a4;
    let c = a2 / a4;
    let d = a1 / a4;
    let e = a0 / a4;
    // Depress: x = y − b/4.
    let p = c - 3.0 * b * b / 8.0;
    let q = d - b * c / 2.0 + b.powi(3) / 8.0;
    let r = e - b * d / 4.0 + b * b * c / 16.0 - 3.0 * b.powi(4) / 256.0;
    let mut y_roots = Vec::new();
    if q.abs() < EPS_DISC {
        // Biquadratic: y⁴ + p·y² + r = 0.
        for z in solve_quadratic([1.0, p, r]) {
            if z >= -EPS_DISC {
                let s = z.max(0.0).sqrt();
                y_roots.push(s);
                if s > EPS_DEDUP {
                    y_roots.push(-s);
                }
            }
        }
    } else {
        // Resolvent cubic: 8z³ + 8p·z² + (2p² − 8r)·z − q² = 0.
        let cubic = solve_cubic([8.0, 8.0 * p, 2.0 * p * p - 8.0 * r, -q * q]);
        // Pick a real root with z > 0 (so √(2z) is real). Largest is
        // most numerically stable.
        let z_opt = cubic
            .iter()
            .rev()
            .find(|&&z| z > EPS_LEAD)
            .copied()
            .or_else(|| cubic.iter().copied().find(|&z| z.abs() > EPS_LEAD));
        let Some(z) = z_opt else {
            return Vec::new();
        };
        if z <= 0.0 {
            return Vec::new();
        }
        let sq2z = (2.0 * z).sqrt();
        let qf = q / (2.0 * sq2z);
        // Ferrari completion:
        //   y⁴ + p·y² + q·y + r
        //     = (y² + p/2 + z)² − (sq2z·y − qf)²
        // → y² + p/2 + z = ±(sq2z·y − qf)
        // → two quadratics:
        //   y² − sq2z·y + (p/2 + z + qf) = 0
        //   y² + sq2z·y + (p/2 + z − qf) = 0
        for y in solve_quadratic([1.0, -sq2z, p / 2.0 + z + qf]) {
            y_roots.push(y);
        }
        for y in solve_quadratic([1.0, sq2z, p / 2.0 + z - qf]) {
            y_roots.push(y);
        }
    }
    let mut roots: Vec<f64> = y_roots.into_iter().map(|y| y - b / 4.0).collect();
    roots.sort_by(|x, y| x.partial_cmp(y).unwrap_or(std::cmp::Ordering::Equal));
    roots.dedup_by(|x, y| (*x - *y).abs() < EPS_DEDUP);
    roots
}

/// Sign-preserving cube root (`f64::cbrt` already handles this on
/// stable Rust, but we wrap for clarity at call sites).
fn cbrt_signed(x: f64) -> f64 {
    x.cbrt()
}
