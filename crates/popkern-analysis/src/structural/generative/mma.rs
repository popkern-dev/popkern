//! The method of moving asymptotes (Svanberg 1987, MMA) for one
//! inequality constraint: minimise f₀(x) subject to f₁(x) ≤ 0 and
//! x̲ ≤ x ≤ x̄.
//!
//! Each update replaces f₀ and f₁ by convex separable approximations
//! (§3, eqs 2–5) in 1/(U − x) and 1/(x − L) about the current point,
//! with p_ij = (U − x)² max(∂f_i/∂x_j, 0) and
//! q_ij = (x − L)² max(−∂f_i/∂x_j, 0), and solves that subproblem
//! exactly through its dual (§4). With one constraint the dual is a
//! concave function of a single multiplier y, whose derivative (the
//! approximated constraint at x(y), eq 19) falls monotonically, so y is
//! found by bisection to round-off. When no point of the move box meets
//! the approximated constraint, the artificial variable of §5 with an
//! unbounded cost takes the point that comes closest. A constraint that
//! is explicit and cheap (a design's volume) can be met exactly instead:
//! the multiplier is then the one whose x(y) meets f₁ itself, as the
//! optimality criteria update meets its volume.
//!
//! The asymptotes move by §3's rule (eqs 11–13): one unit of the
//! variable's range away at the first two updates, then pulled in by
//! s = 0.7 where the variable oscillates and pushed out by 1/√s where
//! it moves steadily (§3's more conservative variant: two relaxations
//! undo one tightening). Their distance from the point is held within
//! [0.01, 10] ranges, so a long run of one or the other neither stalls
//! the variable nor flattens the approximation. The move box is eq 8's
//! (a tenth of the way to each asymptote) intersected with a move limit.

/// Bisections of the multiplier: halving 200 times takes any bracket to
/// the last bit.
const BISECTIONS: usize = 200;

/// §3's s: the asymptotes' contraction where a variable oscillates.
const CONTRACTION: f64 = 0.7;

/// Asymptote distance at the first two updates, in ranges (eq 11).
const INITIAL_DISTANCE: f64 = 1.0;

/// Asymptote distance bounds, in ranges.
const DISTANCE_BOUNDS: (f64, f64) = (0.01, 10.0);

/// An explicit constraint function f₁(x).
pub type Explicit<'a> = dyn Fn(&[f64]) -> f64 + 'a;

/// One convex subproblem P⁽ᵏ⁾ about the point `x`.
pub struct Subproblem<'a> {
    /// The current point x⁽ᵏ⁾.
    pub x: &'a [f64],
    /// Lower asymptotes L (< α).
    pub lower: &'a [f64],
    /// Upper asymptotes U (> β).
    pub upper: &'a [f64],
    /// Lower end of the move box α.
    pub alpha: &'a [f64],
    /// Upper end of the move box β (≥ α).
    pub beta: &'a [f64],
    /// ∂f₀/∂x at x.
    pub df0: &'a [f64],
    /// f₁(x): the constraint holds where it is ≤ 0.
    pub f1: f64,
    /// ∂f₁/∂x at x.
    pub df1: &'a [f64],
    /// f₁ itself, when it is explicit and cheap: the multiplier is then
    /// chosen so the point meets it, not its approximation. It must not
    /// increase with any variable (so f₁(x(y)) falls as y grows).
    pub exact: Option<&'a Explicit<'a>>,
}

impl Subproblem<'_> {
    fn pq(&self, j: usize, d: f64) -> (f64, f64) {
        let (u, l) = (self.upper[j] - self.x[j], self.x[j] - self.lower[j]);
        (u * u * d.max(0.0), l * l * (-d).max(0.0))
    }

    /// The minimiser of the Lagrangian for multiplier `y` (eq 18,
    /// clamped to the move box: §4's cases 1–3). Above y = 1 the
    /// coefficients are divided by y, which leaves the minimiser alone
    /// and keeps them finite however large y grows on the way to an
    /// unreachable constraint.
    fn primal(&self, y: f64) -> Vec<f64> {
        let (w0, w1) = if y > 1.0 { (1.0 / y, 1.0) } else { (1.0, y) };
        (0..self.x.len())
            .map(|j| {
                let (p0, q0) = self.pq(j, self.df0[j]);
                let (p1, q1) = self.pq(j, self.df1[j]);
                let (p, q) = (w0 * p0 + w1 * p1, w0 * q0 + w1 * q1);
                let (a, b) = (self.alpha[j], self.beta[j]);
                let x = if p > 0.0 && q > 0.0 {
                    let (sp, sq) = (p.sqrt(), q.sqrt());
                    (sp * self.lower[j] + sq * self.upper[j]) / (sp + sq)
                } else if p > 0.0 {
                    a
                } else if q > 0.0 {
                    b
                } else {
                    // The rare case: the Lagrangian does not depend on
                    // x_j; it stays.
                    self.x[j]
                };
                x.clamp(a, b)
            })
            .collect()
    }

    /// The approximated constraint f̃₁ at `z`.
    fn constraint(&self, z: &[f64]) -> f64 {
        self.f1
            + (0..self.x.len())
                .map(|j| {
                    let (p, q) = self.pq(j, self.df1[j]);
                    let (u, l) = (self.upper[j], self.lower[j]);
                    let mut s = 0.0;
                    if p > 0.0 {
                        s += p * (1.0 / (u - z[j]) - 1.0 / (u - self.x[j]));
                    }
                    if q > 0.0 {
                        s += q * (1.0 / (z[j] - l) - 1.0 / (self.x[j] - l));
                    }
                    s
                })
                .sum::<f64>()
    }

    /// The subproblem's solution and its multiplier.
    pub fn solve(&self) -> (Vec<f64>, f64) {
        let constraint = |z: &[f64]| match self.exact {
            Some(f) => f(z),
            None => self.constraint(z),
        };
        let free = self.primal(0.0);
        if constraint(&free) <= 0.0 {
            return (free, 0.0);
        }
        // Bracket y: grow the upper end until the constraint holds. If it
        // never does, the largest y gives the point closest to feasible.
        let (mut y1, mut y2) = (0.0, 1.0);
        while constraint(&self.primal(y2)) > 0.0 {
            if y2 > 1e300 {
                return (self.primal(y2), y2);
            }
            y1 = y2;
            y2 *= 2.0;
        }
        for _ in 0..BISECTIONS {
            let mid = 0.5 * (y1 + y2);
            if mid <= y1 || mid >= y2 {
                break;
            }
            if constraint(&self.primal(mid)) > 0.0 {
                y1 = mid;
            } else {
                y2 = mid;
            }
        }
        (self.primal(y2), y2)
    }
}

/// MMA's state across updates: bounds, the last two points and the
/// asymptotes.
#[derive(Debug, Clone)]
pub struct Mma {
    min: Vec<f64>,
    max: Vec<f64>,
    move_limit: f64,
    previous: Vec<Vec<f64>>,
    lower: Vec<f64>,
    upper: Vec<f64>,
}

impl Mma {
    /// MMA over the box [`min`, `max`] with at most `move_limit` change
    /// per variable per update.
    pub fn new(min: Vec<f64>, max: Vec<f64>, move_limit: f64) -> Mma {
        Mma {
            min,
            max,
            move_limit,
            previous: Vec::new(),
            lower: Vec::new(),
            upper: Vec::new(),
        }
    }

    /// The next point from `x`, the objective's gradient `df0`, and the
    /// constraint's value `f1` (≤ 0 where it holds) and gradient `df1`;
    /// `exact` is the constraint itself when it is explicit
    /// ([`Subproblem::exact`]).
    pub fn update(
        &mut self,
        x: &[f64],
        df0: &[f64],
        (f1, df1): (f64, &[f64]),
        exact: Option<&Explicit<'_>>,
    ) -> Vec<f64> {
        let n = x.len();
        let x: Vec<f64> = (0..n)
            .map(|j| x[j].clamp(self.min[j], self.max[j].max(self.min[j])))
            .collect();
        let x = x.as_slice();
        let ranges: Vec<f64> = (0..n).map(|j| self.max[j] - self.min[j]).collect();
        let mut lower = vec![0.0; n];
        let mut upper = vec![0.0; n];
        for j in 0..n {
            let range = ranges[j].max(f64::MIN_POSITIVE);
            let (below, above) = match self.previous.as_slice() {
                [.., older, last] if self.lower.len() == n => {
                    let trend = (x[j] - last[j]) * (last[j] - older[j]);
                    let factor = if trend < 0.0 {
                        CONTRACTION
                    } else if trend > 0.0 {
                        1.0 / CONTRACTION.sqrt()
                    } else {
                        1.0
                    };
                    (
                        factor * (last[j] - self.lower[j]),
                        factor * (self.upper[j] - last[j]),
                    )
                }
                _ => (INITIAL_DISTANCE * range, INITIAL_DISTANCE * range),
            };
            let (lo, hi) = (DISTANCE_BOUNDS.0 * range, DISTANCE_BOUNDS.1 * range);
            lower[j] = x[j] - below.clamp(lo, hi);
            upper[j] = x[j] + above.clamp(lo, hi);
        }
        let alpha: Vec<f64> = (0..n)
            .map(|j| {
                (0.9 * lower[j] + 0.1 * x[j])
                    .max(x[j] - self.move_limit * ranges[j])
                    .max(self.min[j])
                    .min(x[j].max(self.min[j]))
            })
            .collect();
        let beta: Vec<f64> = (0..n)
            .map(|j| {
                (0.9 * upper[j] + 0.1 * x[j])
                    .min(x[j] + self.move_limit * ranges[j])
                    .min(self.max[j])
                    .max(alpha[j])
            })
            .collect();
        let (next, _) = Subproblem {
            x,
            lower: &lower,
            upper: &upper,
            alpha: &alpha,
            beta: &beta,
            df0,
            f1,
            df1,
            exact,
        }
        .solve();
        self.lower = lower;
        self.upper = upper;
        self.previous.push(x.to_vec());
        if self.previous.len() > 2 {
            self.previous.remove(0);
        }
        next
    }
}

#[cfg(test)]
#[path = "mma_tests.rs"]
mod tests;
