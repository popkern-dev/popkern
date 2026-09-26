//! Assignment search for the periodic loop chain (split from
//! `export_periodic_branches.rs` under the 500-line budget).
//!
//! Each loop use gets a traversal direction and a whole-period uv
//! shift so the uses chain into one closed contour. Two freedoms are
//! searched:
//! - the direction of CLOSED uses (both ends at one 3D point), which
//!   geometry alone cannot pin;
//! - the branch of the free coordinate at a POLE junction: where two
//!   uses meet at a collapsed chart point (a sphere pole, where every
//!   `u` maps to one point), the next use may start on any `u` branch
//!   (PK-558: a sphere's `[seam, seam⁻¹]` loop chained both seam uses
//!   onto u = 2π, a zero-area face OpenCASCADE rejects as
//!   unorientable once it reads popkern's pcurves).
//!
//! Candidates are ranked lexicographically: chain gap (closure and
//! per-junction drift, free pole coordinates ignored), then SEAM
//! BRANCH violations (the two uses of one edge in one loop must lie a
//! period apart — the uv contour of a SEAM_CURVE), then reversals of
//! the stored pcurve directions, then the required winding.

/// One loop use as the chain sees it.
pub(super) struct ChainUse {
    /// The use's edge (seam uses share it).
    pub edge: popkern_core::topology::EdgeHandle,
    /// Both pcurve endpoints, in pcurve parameter order.
    pub end_a: (f64, f64),
    pub end_b: (f64, f64),
    /// `end_a` is the walk start when known from 3D (open edges).
    pub a_is_start: Option<bool>,
    /// Sampled uv polyline in pcurve parameter order.
    pub samples: Vec<(f64, f64)>,
}

/// Per-use `(shift, forward)` in loop-walk order.
pub(super) type ChainAssignment = Vec<((f64, f64), bool)>;

/// The winning assignment: `(gap, seam violations, flips, winding_ok,
/// assignment)`.
pub(super) type ChainChoice = (f64, usize, u32, bool, ChainAssignment);

/// Chain gap treated as exact closure.
const EXACT_CLOSURE: f64 = 1e-9;

/// Most candidate assignments enumerated per loop.
const MAX_CANDIDATES: usize = 1 << 14;

/// The best chain assignment of `uses` (see the module docs).
/// `pole_at(uv)` reports which coordinates are free at `uv`
/// (`(u_free, v_free)`); `periods` are the surface's periods.
pub(super) fn best_assignment(
    uses: &[ChainUse],
    periods: (Option<f64>, Option<f64>),
    want_ccw: bool,
    pole_at: &dyn Fn((f64, f64)) -> (bool, bool),
) -> Option<ChainChoice> {
    let search = Search {
        uses,
        periods,
        want_ccw,
        pole_at,
    };
    let mut best: Option<ChainChoice> = None;
    let mut budget = MAX_CANDIDATES;
    let first_dirs: &[bool] = match uses.first()?.a_is_start {
        Some(true) => &[true],
        Some(false) => &[false],
        None => &[true, false],
    };
    for &fwd0 in first_dirs {
        let (_, e0) = oriented(&uses[0], fwd0);
        let flips = u32::from(uses[0].a_is_start.is_none() && !fwd0);
        let mut assignment = vec![((0.0, 0.0), fwd0)];
        search.extend(&mut assignment, e0, 0.0, flips, &mut best, &mut budget);
    }
    best
}

struct Search<'a> {
    uses: &'a [ChainUse],
    periods: (Option<f64>, Option<f64>),
    want_ccw: bool,
    pole_at: &'a dyn Fn((f64, f64)) -> (bool, bool),
}

/// `(start, end)` of `use_` walked forward or reversed.
fn oriented(use_: &ChainUse, forward: bool) -> ((f64, f64), (f64, f64)) {
    if forward {
        (use_.end_a, use_.end_b)
    } else {
        (use_.end_b, use_.end_a)
    }
}

impl Search<'_> {
    /// Whole-period shift moving `from` nearest `to`.
    fn shift_for(&self, from: (f64, f64), to: (f64, f64)) -> (f64, f64) {
        let whole = |p: Option<f64>, d: f64| p.map(|p| (d / p).round() * p).unwrap_or(0.0);
        (
            whole(self.periods.0, to.0 - from.0),
            whole(self.periods.1, to.1 - from.1),
        )
    }

    /// Distance between `a` and `b`, ignoring the coordinates that are
    /// free at the junction `at`.
    fn junction_gap(&self, a: (f64, f64), b: (f64, f64), at: (f64, f64)) -> f64 {
        let (u_free, v_free) = (self.pole_at)(at);
        let du = if u_free { 0.0 } else { a.0 - b.0 };
        let dv = if v_free { 0.0 } else { a.1 - b.1 };
        (du * du + dv * dv).sqrt()
    }

    /// Branch alternatives for a use starting at pole junction `at`.
    fn pole_offsets(&self, at: (f64, f64)) -> Vec<(f64, f64)> {
        let (u_free, v_free) = (self.pole_at)(at);
        let options = |free: bool, period: Option<f64>| match (free, period) {
            (true, Some(p)) => vec![0.0, p, -p],
            _ => vec![0.0],
        };
        let mut out = Vec::new();
        for du in options(u_free, self.periods.0) {
            for dv in options(v_free, self.periods.1) {
                out.push((du, dv));
            }
        }
        out
    }

    fn extend(
        &self,
        assignment: &mut ChainAssignment,
        cursor: (f64, f64),
        drift: f64,
        flips: u32,
        best: &mut Option<ChainChoice>,
        budget: &mut usize,
    ) {
        if *budget == 0 {
            return;
        }
        let idx = assignment.len();
        if idx == self.uses.len() {
            *budget -= 1;
            self.score(assignment, cursor, drift, flips, best);
            return;
        }
        let use_ = &self.uses[idx];
        let dirs: &[bool] = match use_.a_is_start {
            Some(true) => &[true],
            Some(false) => &[false],
            None => &[true, false],
        };
        for &fwd in dirs {
            let (start, end) = oriented(use_, fwd);
            let base = self.shift_for(start, cursor);
            for (pu, pv) in self.pole_offsets(cursor) {
                let shift = (base.0 + pu, base.1 + pv);
                let shifted_start = (start.0 + shift.0, start.1 + shift.1);
                let gap = self.junction_gap(shifted_start, cursor, cursor);
                assignment.push((shift, fwd));
                self.extend(
                    assignment,
                    (end.0 + shift.0, end.1 + shift.1),
                    drift.max(gap),
                    flips + u32::from(use_.a_is_start.is_none() && !fwd),
                    best,
                    budget,
                );
                assignment.pop();
            }
        }
    }

    fn score(
        &self,
        assignment: &ChainAssignment,
        cursor: (f64, f64),
        drift: f64,
        flips: u32,
        best: &mut Option<ChainChoice>,
    ) {
        let (s0, _) = oriented(&self.uses[0], assignment[0].1);
        // A chain within round-off of closing is exactly closed: flips
        // and winding, not float noise, decide between such chains.
        let gap = match drift.max(self.junction_gap(cursor, s0, s0)) {
            g if g <= EXACT_CLOSURE => 0.0,
            g => g,
        };
        let violations = self.seam_violations(assignment);
        let winding_ok = (self.signed_area(assignment) > 0.0) == self.want_ccw;
        let better = best.as_ref().is_none_or(|(g, v, f, w, _)| {
            (gap, violations, flips, !winding_ok) < (*g, *v, *f, !*w)
        });
        if better {
            *best = Some((gap, violations, flips, winding_ok, assignment.clone()));
        }
    }

    /// Pairs of uses of one edge whose shifted traces coincide (the
    /// same branch) instead of lying a period apart.
    fn seam_violations(&self, assignment: &ChainAssignment) -> usize {
        let traced = |i: usize| -> Vec<(f64, f64)> {
            let ((du, dv), _) = assignment[i];
            self.uses[i]
                .samples
                .iter()
                .map(|p| (p.0 + du, p.1 + dv))
                .collect()
        };
        let mut violations = 0;
        for i in 0..self.uses.len() {
            for j in (i + 1)..self.uses.len() {
                if self.uses[i].edge != self.uses[j].edge {
                    continue;
                }
                let (a, mut b) = (traced(i), traced(j));
                let coincide = |b: &[(f64, f64)]| {
                    a.iter().zip(b).all(|(p, q)| {
                        let scale = 1.0 + p.0.abs().max(p.1.abs());
                        (p.0 - q.0).abs() <= 1e-9 * scale && (p.1 - q.1).abs() <= 1e-9 * scale
                    })
                };
                // The two uses may store their pcurves in either
                // parameter direction.
                let same = coincide(&b) || {
                    b.reverse();
                    coincide(&b)
                };
                violations += usize::from(same);
            }
        }
        violations
    }

    /// Twice the signed area of the traversed contour (shoelace over
    /// the shifted samples; pole traverses are the implicit closing
    /// segments between consecutive uses).
    fn signed_area(&self, assignment: &ChainAssignment) -> f64 {
        let mut area2 = 0.0_f64;
        let mut prev: Option<(f64, f64)> = None;
        let mut first: Option<(f64, f64)> = None;
        for (use_, ((du, dv), fwd)) in self.uses.iter().zip(assignment) {
            let mut pts = use_.samples.clone();
            if !*fwd {
                pts.reverse();
            }
            for p in pts {
                let q = (p.0 + du, p.1 + dv);
                if let Some(pr) = prev {
                    area2 += pr.0 * q.1 - q.0 * pr.1;
                }
                first.get_or_insert(q);
                prev = Some(q);
            }
        }
        if let (Some(pr), Some(f0)) = (prev, first) {
            area2 += pr.0 * f0.1 - f0.0 * pr.1;
        }
        area2
    }
}

#[cfg(test)]
#[path = "export_periodic_chain_tests.rs"]
mod tests;
