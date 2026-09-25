//! Projected seam correspondence (PK-499 design §Correspondence):
//! when the two sides are not fraction-co-parameterized (a closed rail
//! whose pcurves start at different phases, a rational-arc parameter
//! against an angle), side b is compared at the point of its seam
//! curve lying in the normal plane of side a's seam curve.
//!
//! In side b's own seam fraction σ that point solves
//! `f(s, σ) = (P_b(σ) − P_a(s)) · T_a(s) = 0`. The condition uses side
//! b only through its position, which is continuous across every knot
//! and (within the closure tolerance) across the closure of a closed
//! side, so no tangent continuity of side b is assumed.
//!
//! # Branches
//!
//! Side b is split into smooth branches over a σ-set Σ: each knot
//! span's polynomial (extended over all of Σ) and, for a closed side,
//! each shift by a whole period that meets `[0, 1]`. Per branch and
//! per s-interval `T` the Krawczyk operator
//! `K = σ̃ − C f(m, σ̃) + (1 − C f_σ(T, Σ))(Σ − σ̃) − C f_s(T, Σ)(T − m)`
//! either certifies no root in Σ (`K ∩ Σ = ∅`: the branch drops out),
//! or a unique root in `K ⊂ int Σ` for every `s ∈ T`, or fails (the
//! interval is refined). The true correspondence lies in the branch
//! whose span contains it, so the union of the surviving branches'
//! enclosures encloses it.
//!
//! The Taylor coefficients of `σ(s)` follow order by order from
//! `F(s) = f(s, σ(s)) ≡ 0`: coefficient `k` of `F` is `f_σ σ_k` plus
//! terms in lower coefficients, so `σ_k = −[F]_k / f_σ` evaluated with
//! `σ_k = 0`. The branch's σ-series state is composed with that
//! series. `f_σ` and `f_s` over `T × Σ` are enclosed by Taylor forms
//! in σ: a natural enclosure of a turning rail's nearly constant
//! `|T|²` overestimates it many times over and stalls the contraction
//! (measured on the PK-389 rim: `[7.6, 78]` around 39.4 at
//! `|Σ| = 0.06`).

use super::CertifyDecline;
use super::interval::Iv;
use super::seam_state::{SeamBranch, SeamSideSource, Series, SeriesVec, SideState};
use super::taylor::SEAM_TERMS;

/// Why a projected correspondence could not be certified on an
/// interval.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ProjectionFailure {
    /// A side could not be enclosed.
    Side(CertifyDecline),
    /// Krawczyk could neither certify a unique root nor exclude one.
    NotContracting,
}

/// Most smooth branches of side b tried per σ-set.
const MAX_BRANCHES: usize = 32;

/// Seeding resolution of side b's seam.
const SEED_SAMPLES: usize = 256;

/// Side b with its closure resolved.
#[derive(Debug, Clone)]
pub struct ProjectedSide {
    source: SeamSideSource,
    closed: bool,
    /// Coarse float images of side b's seam (σ, P) for seeding.
    table: Vec<(f64, [f64; 3])>,
}

/// A smooth branch of side b: a knot branch at a period shift.
#[derive(Debug, Clone, Copy, PartialEq)]
struct Branch {
    shift: f64,
    knots: SeamBranch,
}

/// Side b corresponded to side a over an s-interval.
#[derive(Debug, Clone)]
pub struct ProjectedStates {
    /// `(over T, at m)` per surviving branch, as series in s.
    pub pieces: Vec<(SideState, SideState)>,
    /// Side b's true state (in σ; value coefficients meaningful) at
    /// each candidate for the midpoint's corresponding point: every
    /// root the true, hulled side b also puts on the normal plane.
    pub witnesses: Vec<SideState>,
}

impl ProjectedSide {
    /// Prepare side b. `closure_tolerance` is the largest seam-image
    /// gap at σ = 0 ↔ 1 accepted as a closed side.
    pub fn new(mut source: SeamSideSource, closure_tolerance: f64) -> Result<Self, CertifyDecline> {
        let mut table = Vec::with_capacity(SEED_SAMPLES + 1);
        for k in 0..=SEED_SAMPLES {
            let sigma = k as f64 / SEED_SAMPLES as f64;
            let p = source.state(Iv::point(sigma))?.0.p.value().mid();
            table.push((sigma, [p.x, p.y, p.z]));
        }
        let (start, end) = (table[0].1, table[SEED_SAMPLES].1);
        let gap = ((start[0] - end[0]).powi(2)
            + (start[1] - end[1]).powi(2)
            + (start[2] - end[2]).powi(2))
        .sqrt();
        Ok(Self {
            source,
            closed: gap <= closure_tolerance,
            table,
        })
    }

    /// Smooth branches meeting the σ-set, with their states over it
    /// (in each branch's own, shifted, fraction).
    fn branches(&mut self, sigma: Iv) -> Result<Vec<(Branch, SideState)>, ProjectionFailure> {
        let shifts: Vec<f64> = if self.closed {
            [-1.0, 0.0, 1.0]
                .into_iter()
                .filter(|k| sigma.hi - k >= 0.0 && sigma.lo - k <= 1.0)
                .collect()
        } else {
            vec![0.0]
        };
        let mut out = Vec::new();
        for shift in shifts {
            let local = sigma - Iv::point(shift);
            for (knots, state) in self
                .source
                .state_branches(local)
                .map_err(ProjectionFailure::Side)?
            {
                out.push((Branch { shift, knots }, state));
            }
        }
        if out.len() > MAX_BRANCHES {
            return Err(ProjectionFailure::NotContracting);
        }
        Ok(out)
    }

    fn branch_state(&mut self, branch: Branch, sigma: Iv) -> Result<SideState, ProjectionFailure> {
        self.source
            .state_in(sigma - Iv::point(branch.shift), branch.knots)
            .map_err(ProjectionFailure::Side)
    }

    /// Float guess of the corresponding σ for the point `p` with seam
    /// direction `t`: nearest seed, then Newton on
    /// `(P_b(σ) − p) · t = 0`.
    fn guess(&mut self, p: [f64; 3], t: [f64; 3]) -> f64 {
        let dist = |q: &[f64; 3]| (0..3).map(|k| (q[k] - p[k]).powi(2)).sum::<f64>();
        let mut sigma = self
            .table
            .iter()
            .min_by(|a, b| dist(&a.1).total_cmp(&dist(&b.1)))
            .map_or(0.5, |(s, _)| *s);
        for _ in 0..8 {
            let wrapped = if self.closed {
                sigma - sigma.floor()
            } else {
                sigma
            };
            let Ok((state, _)) = self.source.state(Iv::point(wrapped)) else {
                break;
            };
            let (pb, tb) = (state.p.value().mid(), state.tangent.value().mid());
            let f = (pb.x - p[0]) * t[0] + (pb.y - p[1]) * t[1] + (pb.z - p[2]) * t[2];
            let fs = tb.x * t[0] + tb.y * t[1] + tb.z * t[2];
            if fs.abs() <= 0.0 || !fs.is_finite() {
                break;
            }
            let step = f / fs;
            sigma -= step;
            if step.abs() < 1e-15 {
                break;
            }
        }
        sigma
    }

    /// A rough `|dP_b/dσ|` from the seed table (radius scale only).
    fn table_speed(&self) -> f64 {
        self.table
            .windows(2)
            .map(|w| {
                (0..3)
                    .map(|k| (w[1].1[k] - w[0].1[k]).powi(2))
                    .sum::<f64>()
                    .sqrt()
            })
            .sum()
    }

    /// Whether side b's seam closes (its images at σ = 0 and 1
    /// coincide within the closure tolerance).
    pub fn is_closed(&self) -> bool {
        self.closed
    }

    /// Float estimate of a phase-shifted correspondence
    /// `σ(s) = φ + d·s (mod 1)`, `d = ±1`, from the foot points of side
    /// a's seam at eight fractions — `None` unless side b is closed and
    /// the foot points are affine in `s` to `1e-6`. Both sides of a
    /// closed rail parameterized by angle, starting at different
    /// phases, have exactly this correspondence (PK-542).
    pub fn phase_estimate(&mut self, a: &mut SeamSideSource) -> Option<(f64, f64)> {
        if !self.closed {
            return None;
        }
        const SAMPLES: usize = 8;
        let mut sigmas = Vec::with_capacity(SAMPLES + 1);
        for k in 0..=SAMPLES {
            let (state, _) = a.state(Iv::point(k as f64 / SAMPLES as f64)).ok()?;
            let (p, t) = (state.p.value().mid(), state.tangent.value().mid());
            let sigma = self.guess([p.x, p.y, p.z], [t.x, t.y, t.z]);
            sigmas.push(sigma - sigma.floor());
        }
        let wrap = |x: f64| x - x.round();
        let steps: Vec<f64> = sigmas.windows(2).map(|w| wrap(w[1] - w[0])).collect();
        let total: f64 = steps.iter().sum();
        let dir = if total >= 0.0 { 1.0 } else { -1.0 };
        let expected = dir / SAMPLES as f64;
        if (total - dir).abs() > 1e-6 || steps.iter().any(|d| (d - expected).abs() > 1e-6) {
            return None;
        }
        Some((sigmas[0], dir))
    }

    /// Side b corresponded to side a by the fixed phase shift
    /// `σ(s) = phase + dir·s` over the s-interval `s` (midpoint `m`):
    /// per smooth branch meeting `σ(s)`, its state composed with the
    /// linear series `Δσ = dir·δ`, and the true state at `σ(m)` as the
    /// witness. The correspondence is certified by the gap bound
    /// itself, as fraction correspondence is.
    pub fn shifted_states(
        &mut self,
        s: Iv,
        m: f64,
        phase: f64,
        dir: f64,
    ) -> Result<ProjectedStates, ProjectionFailure> {
        let map = |x: Iv| Iv::point(phase) + Iv::point(dir) * x;
        let (sigma, sigma_m) = (map(s), map(Iv::point(m)));
        let mut delta = Series::constant(Iv::point(0.0));
        delta.c[1] = Iv::point(dir);
        let mut out = ProjectedStates {
            pieces: Vec::new(),
            witnesses: Vec::new(),
        };
        for (branch, b_over) in self.branches(sigma)? {
            let b_m = self.branch_state(branch, sigma_m)?;
            out.pieces.push((b_over.compose(delta), b_m.compose(delta)));
        }
        let wrapped = sigma_m - Iv::point(sigma_m.mid().floor());
        if let Ok((truth, _)) = self.source.state(wrapped) {
            out.witnesses.push(truth);
        }
        Ok(out)
    }

    /// Side b corresponded to one smooth branch of side a (`a` over
    /// the s-interval `s`, `a_at_m` at its midpoint `m`).
    pub fn states_along(
        &mut self,
        a: &SideState,
        s: Iv,
        a_at_m: &SideState,
        m: f64,
    ) -> Result<ProjectedStates, ProjectionFailure> {
        let (pa, ta) = (a_at_m.p.value().mid(), a_at_m.tangent.value().mid());
        let guess = self.guess([pa.x, pa.y, pa.z], [ta.x, ta.y, ta.z]);
        let rate = a_at_m.tangent.value().mag() / self.table_speed().max(1e-300);
        let mut radius = (1.5 * rate * s.width()).max(1e-12);
        for _ in 0..4 {
            let sigma = Iv::new(guess - radius, guess + radius);
            match self.contract_all(a, s, a_at_m, m, guess, sigma) {
                Ok(states) => return Ok(states),
                Err(ProjectionFailure::NotContracting) => radius *= 2.0,
                Err(other) => return Err(other),
            }
        }
        Err(ProjectionFailure::NotContracting)
    }

    fn contract_all(
        &mut self,
        a: &SideState,
        s: Iv,
        a_at_m: &SideState,
        m: f64,
        guess: f64,
        sigma: Iv,
    ) -> Result<ProjectedStates, ProjectionFailure> {
        let mut out = ProjectedStates {
            pieces: Vec::new(),
            witnesses: Vec::new(),
        };
        let point = Iv::point(guess);
        let offset_s = s - Iv::point(m);
        for (branch, b_over) in self.branches(sigma)? {
            let b_point = self.branch_state(branch, point)?;
            let Some(root) = krawczyk(a, a_at_m, &b_point, &b_over, point, sigma, offset_s)? else {
                continue; // certified: no root on this branch
            };
            // Series of σ(s) on the (tighter) root set.
            let b_root = self.branch_state(branch, root)?;
            let (f_sigma, _) = derivatives(a, &b_point, &b_root, root - point);
            let over = b_root.compose(sigma_series(a, &b_root, f_sigma));
            // The same branch's root for s = m alone, contracted to
            // roundoff: one step from the interval's root set leaves an
            // O(h²) enclosure, which loosens every midpoint coefficient
            // and drops the Taylor forms to O(h³) (measured).
            // The midpoint's root: a fresh Krawczyk box of the Newton
            // step's size around σ̃ (inside `root`, where uniqueness is
            // already certified), then one refining image: roundoff
            // tight without contracting down from the interval's root
            // set. Any Krawczyk image contains every root of its box.
            let f_point = (b_point.p.value() - a_at_m.p.value()).dot(a_at_m.tangent.value());
            let newton = (f_point / b_point.tangent.value().dot(a_at_m.tangent.value())).mag();
            let radius = 4.0 * newton + 1e-15 * (1.0 + guess.abs());
            let seed = Iv::new(guess - radius, guess + radius);
            let seed = Iv::new(seed.lo.max(root.lo), seed.hi.min(root.hi));
            let mut root_m = if seed.lo <= seed.hi { seed } else { root };
            for _ in 0..2 {
                let b_iter = self.branch_state(branch, root_m)?;
                let next = krawczyk_image(
                    a_at_m,
                    a_at_m,
                    &b_point,
                    &b_iter,
                    point,
                    root_m,
                    Iv::point(0.0),
                )?;
                let (lo, hi) = (next.lo.max(root_m.lo), next.hi.min(root_m.hi));
                if lo > hi {
                    // The seed missed the root: fall back to the root set.
                    root_m = root;
                    continue;
                }
                root_m = Iv::new(lo, hi);
            }
            let b_m = self.branch_state(branch, root_m)?;
            let (f_sigma_m, _) = derivatives(a_at_m, &b_point, &b_m, root_m - point);
            let at_m = b_m.compose(sigma_series(a_at_m, &b_m, f_sigma_m));
            out.pieces.push((over, at_m));
            // A witness only where the true side b also puts the root
            // on the normal plane.
            if let Ok((truth, _)) = self.source.state(root_m - Iv::point(branch.shift)) {
                let f_true = (truth.p.value() - a_at_m.p.value()).dot(a_at_m.tangent.value());
                if f_true.contains_zero() {
                    out.witnesses.push(truth);
                }
            }
        }
        if out.pieces.is_empty() {
            // No branch has a root in Σ: Σ was too small or misplaced.
            return Err(ProjectionFailure::NotContracting);
        }
        Ok(out)
    }
}

/// Enclosures of `f_σ = T_b · T_a` and
/// `f_s = (P_b − P_a) · T_a′ − T_a · T_a` over `T × Σ`: side a over
/// `T` naturally, side b by Taylor forms in σ about `σ̃` (`b_point`
/// at σ̃, `b_over` over Σ, `offset = Σ − σ̃`).
fn derivatives(a: &SideState, b_point: &SideState, b_over: &SideState, offset: Iv) -> (Iv, Iv) {
    let (pa, ta, ta1) = (a.p.value(), a.tangent.value(), a.tangent.coefficient(1));
    let series = |b: &SideState| {
        let f_sigma = b.tangent.dot(SeriesVec::constant(ta));
        let f_s = (b.p - SeriesVec::constant(pa)).dot(SeriesVec::constant(ta1))
            - Series::constant(ta.dot(ta));
        (f_sigma, f_s)
    };
    let ((sigma_point, s_point), (sigma_over, s_over)) = (series(b_point), series(b_over));
    (
        sigma_point.enclose_over(&sigma_over, offset),
        s_point.enclose_over(&s_over, offset),
    )
}

/// The Krawczyk image `K` over `T × Σ` (it contains every root in Σ
/// of every `s ∈ T`).
fn krawczyk_image(
    a: &SideState,
    a_at_m: &SideState,
    b_point: &SideState,
    b_over: &SideState,
    point: Iv,
    sigma: Iv,
    offset_s: Iv,
) -> Result<Iv, ProjectionFailure> {
    let ta_m = a_at_m.tangent.value();
    let f_point = (b_point.p.value() - a_at_m.p.value()).dot(ta_m);
    let c = 1.0 / b_point.tangent.value().dot(ta_m).mid();
    if !c.is_finite() {
        return Err(ProjectionFailure::NotContracting);
    }
    let c_iv = Iv::point(c);
    let (f_sigma, f_s) = derivatives(a, b_point, b_over, sigma - point);
    Ok(
        point - c_iv * f_point + (Iv::point(1.0) - c_iv * f_sigma) * (sigma - point)
            - c_iv * f_s * offset_s,
    )
}

/// Krawczyk on one branch: `Ok(None)` certifies no root in Σ,
/// `Ok(Some(K))` a unique root in `K ⊂ int Σ` for every `s ∈ T`.
fn krawczyk(
    a: &SideState,
    a_at_m: &SideState,
    b_point: &SideState,
    b_over: &SideState,
    point: Iv,
    sigma: Iv,
    offset_s: Iv,
) -> Result<Option<Iv>, ProjectionFailure> {
    let k = krawczyk_image(a, a_at_m, b_point, b_over, point, sigma, offset_s)?;
    if k.hi < sigma.lo || k.lo > sigma.hi {
        return Ok(None);
    }
    if k.lo > sigma.lo && k.hi < sigma.hi {
        return Ok(Some(k));
    }
    Err(ProjectionFailure::NotContracting)
}

/// Taylor series of `σ(s) − σ₀` from `F(s) = (P_b(σ(s)) − P_a(s)) ·
/// T_a(s) ≡ 0`, solved order by order with the enclosure `f_sigma`
/// of `∂f/∂σ` (see module docs).
fn sigma_series(a: &SideState, b: &SideState, f_sigma: Iv) -> Series {
    let mut delta = Series::constant(Iv::point(0.0));
    for k in 1..SEAM_TERMS {
        delta.c[k] = Iv::point(0.0);
        let residual = (b.p.compose(delta) - a.p).dot(a.tangent);
        delta.c[k] = -(residual.c[k] / f_sigma);
    }
    delta
}
