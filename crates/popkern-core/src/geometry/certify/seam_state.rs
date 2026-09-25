//! Per-side seam state along the seam fraction `s`, and the
//! continuity quantities compared across a seam, all as truncated
//! Taylor series in `s` ([`Series`]).
//!
//! A side is a surface with the pcurve of its half-edge, walked in
//! the edge-curve direction: pcurve parameter
//! `t(s) = t0 + g(s)·(t1 − t0)` with `g = s` or `1 − s`. With the
//! pcurve increments `Δu(δ) = u(s+δ) − u(s)` (and `Δv`) as series,
//! every surface partial `X = ∂^{a,b}S` expands by Faà di Bruno as
//! `X(s+δ) = Σ_{i,j} ∂^{a+i,b+j}S · (Δu^i/i!)(Δv^j/j!)`, truncated;
//! the seam tangent is the derivative of the (one order longer)
//! position series.
#![allow(clippy::needless_range_loop)] // Faà di Bruno sums pair power tables by index (`Δu^i Δv^j`).

use super::CertifyDecline;
use super::interval::Iv;
use super::pcurve_jet::{PCURVE_JET_TERMS, PcurveJetSource};
use super::surface_jet::{SurfaceJet, SurfaceJetSource};
use super::taylor::{SEAM_TERMS, Taylor, TaylorVec3};
use crate::geometry::pcurve::Pcurve;
use crate::geometry::surface::SurfaceGeometry;
use std::ops::Div;

/// A seam quantity's Taylor series in `s`.
pub type Series = Taylor<SEAM_TERMS>;

/// A seam vector quantity's Taylor series in `s`.
pub type SeriesVec = TaylorVec3<SEAM_TERMS>;

/// One side of a seam, prepared for enclosure.
#[derive(Debug, Clone)]
pub struct SeamSideSource {
    surface: SurfaceJetSource,
    pcurve: PcurveJetSource,
    t_range: (f64, f64),
    reversed: bool,
}

/// A side's seam state over an s-set, as series in `s`.
#[derive(Debug, Clone, Copy)]
pub struct SideState {
    /// Position.
    pub p: SeriesVec,
    /// `S_u`.
    pub su: SeriesVec,
    /// `S_v`.
    pub sv: SeriesVec,
    /// `S_uu`.
    pub suu: SeriesVec,
    /// `S_uv`.
    pub suv: SeriesVec,
    /// `S_vv`.
    pub svv: SeriesVec,
    /// Seam tangent `dP/ds`.
    pub tangent: SeriesVec,
}

impl SideState {
    /// Component-wise hull (used across a closure or knot).
    pub fn hull(&self, o: &SideState) -> SideState {
        SideState {
            p: self.p.hull(o.p),
            su: self.su.hull(o.su),
            sv: self.sv.hull(o.sv),
            suu: self.suu.hull(o.suu),
            suv: self.suv.hull(o.suv),
            svv: self.svv.hull(o.svv),
            tangent: self.tangent.hull(o.tangent),
        }
    }

    /// Every field composed with `σ(s) = σ₀ + Δ(s)` (see
    /// [`Taylor::compose`]).
    pub fn compose(&self, delta: Series) -> SideState {
        SideState {
            p: self.p.compose(delta),
            su: self.su.compose(delta),
            sv: self.sv.compose(delta),
            suu: self.suu.compose(delta),
            suv: self.suv.compose(delta),
            svv: self.svv.compose(delta),
            tangent: self.tangent.compose(delta),
        }
    }
}

/// `Σ_{i+j<M} ∂^{a+i,b+j}S · (Δu^i/i!)(Δv^j/j!)`.
fn expand<const M: usize>(
    jet: &SurfaceJet,
    (a, b): (usize, usize),
    du_powers: &[Taylor<M>],
    dv_powers: &[Taylor<M>],
) -> TaylorVec3<M> {
    let zero = Taylor::<M>::constant(Iv::point(0.0));
    let mut out = TaylorVec3 {
        x: zero,
        y: zero,
        z: zero,
    };
    for i in 0..M {
        for j in 0..M - i {
            let weight = du_powers[i] * dv_powers[j];
            let partial = jet.partial(a + i, b + j);
            out = out
                + TaylorVec3 {
                    x: weight.scale(partial.x),
                    y: weight.scale(partial.y),
                    z: weight.scale(partial.z),
                };
        }
    }
    out
}

/// `Δ^i / i!` for `i < M`.
fn scaled_powers<const M: usize>(delta: Taylor<M>) -> Vec<Taylor<M>> {
    let mut out = Vec::with_capacity(M);
    let mut power = Taylor::<M>::constant(Iv::point(1.0));
    for i in 0..M {
        if i > 0 {
            power = (power * delta).scale(Iv::point(1.0) / Iv::point(i as f64));
        }
        out.push(power);
    }
    out
}

impl SeamSideSource {
    /// Prepare a side. `reversed` means the pcurve runs against the
    /// edge curve.
    pub fn new(
        surface: &SurfaceGeometry,
        pcurve: &Pcurve,
        reversed: bool,
    ) -> Result<Self, CertifyDecline> {
        Ok(Self {
            surface: SurfaceJetSource::new(surface)?,
            pcurve: PcurveJetSource::new(pcurve)?,
            t_range: (pcurve.t_start, pcurve.t_end),
            reversed,
        })
    }

    /// Seam fractions of the pcurve's interior knots.
    pub fn fraction_breaks(&self) -> Vec<f64> {
        let (t0, t1) = self.t_range;
        self.pcurve
            .breaks()
            .into_iter()
            .map(|t| {
                let g = (t - t0) / (t1 - t0);
                if self.reversed { 1.0 - g } else { g }
            })
            .filter(|s| *s > 0.0 && *s < 1.0)
            .collect()
    }

    /// `t(s)` and `dt/ds`.
    fn pcurve_parameter(&self, s: Iv) -> (Iv, Iv) {
        let (t0, t1) = (Iv::point(self.t_range.0), Iv::point(self.t_range.1));
        let span = t1 - t0;
        let (g, dg) = if self.reversed {
            (Iv::point(1.0) - s, Iv::point(-1.0))
        } else {
            (s, Iv::point(1.0))
        };
        (t0 + g * span, dg * span)
    }

    /// State over the s-interval `s`, hulled across every knot branch
    /// it meets, and whether it met more than one (the Taylor forms
    /// are then unsound for quantities that jump at a knot).
    pub fn state(&mut self, s: Iv) -> Result<(SideState, bool), CertifyDecline> {
        let (t, dt) = self.pcurve_parameter(s);
        let (pc, pc_straddle) = self.pcurve.jet_piecewise(t)?;
        let (du, dv) = increments(&pc, dt);
        let (jet, surface_straddle) = self.surface.jet_piecewise(pc.d[0].u, pc.d[0].v)?;
        Ok((assemble(&jet, du, dv), pc_straddle || surface_straddle))
    }

    /// One state per smooth branch the s-interval meets: each pairs a
    /// pcurve span with a surface span, both polynomials extended over
    /// the whole interval. Every `s ∈ T` is represented exactly by the
    /// branch whose spans contain its images, and every branch is
    /// smooth on `T`, so a quantity's Taylor forms may be taken per
    /// branch and united.
    pub fn state_branches(
        &mut self,
        s: Iv,
    ) -> Result<Vec<(SeamBranch, SideState)>, CertifyDecline> {
        let (t, dt) = self.pcurve_parameter(s);
        let mut out = Vec::new();
        for pcurve in self.pcurve.spans_touching(t) {
            let pc = self.pcurve.jet_in(t, pcurve)?;
            let (du, dv) = increments(&pc, dt);
            for surface in self.surface.spans_touching(pc.d[0].u, pc.d[0].v)? {
                let jet = self.surface.jet_in(pc.d[0].u, pc.d[0].v, surface)?;
                out.push((SeamBranch { pcurve, surface }, assemble(&jet, du, dv)));
            }
        }
        Ok(out)
    }

    /// State over `s` of one branch (see [`Self::state_branches`]).
    pub fn state_in(&mut self, s: Iv, branch: SeamBranch) -> Result<SideState, CertifyDecline> {
        let (t, dt) = self.pcurve_parameter(s);
        let pc = self.pcurve.jet_in(t, branch.pcurve)?;
        let (du, dv) = increments(&pc, dt);
        let jet = self.surface.jet_in(pc.d[0].u, pc.d[0].v, branch.surface)?;
        Ok(assemble(&jet, du, dv))
    }
}

/// A smooth branch of a side: which pcurve span and which surface
/// span (`None` for kinds with a single smooth branch).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SeamBranch {
    /// Pcurve knot span.
    pub pcurve: Option<usize>,
    /// Surface knot span `(u, v)`.
    pub surface: Option<(usize, usize)>,
}

/// Increment series `Δu`, `Δv` in s: coefficient k is
/// `dᵏu/dtᵏ · (dt/ds)ᵏ / k!`.
fn increments(
    pc: &super::pcurve_jet::PcurveJet,
    dt: Iv,
) -> ([Iv; PCURVE_JET_TERMS], [Iv; PCURVE_JET_TERMS]) {
    let mut du = [Iv::point(0.0); PCURVE_JET_TERMS];
    let mut dv = [Iv::point(0.0); PCURVE_JET_TERMS];
    let (mut rate, mut factorial) = (Iv::point(1.0), 1.0);
    for k in 1..PCURVE_JET_TERMS {
        rate = rate * dt;
        factorial *= k as f64;
        let scale = rate / Iv::point(factorial);
        du[k] = pc.d[k].u * scale;
        dv[k] = pc.d[k].v * scale;
    }
    (du, dv)
}

/// The side state from a surface jet and pcurve increments.
fn assemble(jet: &SurfaceJet, du: [Iv; PCURVE_JET_TERMS], dv: [Iv; PCURVE_JET_TERMS]) -> SideState {
    // Position one order longer: its derivative is the tangent.
    let (du_long, dv_long) = (Taylor::new(du), Taylor::new(dv));
    let position = expand(
        jet,
        (0, 0),
        &scaled_powers(du_long),
        &scaled_powers(dv_long),
    );
    let (du, dv) = (
        du_long.truncate::<SEAM_TERMS>(),
        dv_long.truncate::<SEAM_TERMS>(),
    );
    let (pu, pv) = (scaled_powers(du), scaled_powers(dv));
    let at = |a: usize, b: usize| expand(jet, (a, b), &pu, &pv);
    SideState {
        p: position.truncate(),
        su: at(1, 0),
        sv: at(0, 1),
        suu: at(2, 0),
        suv: at(1, 1),
        svv: at(0, 2),
        tangent: position.derivative_series().truncate(),
    }
}

/// The compared quantities at an s-set. Each is a scalar component
/// whose magnitude the certificate bounds.
#[derive(Debug, Clone, Copy)]
pub struct SeamQuantities {
    /// `P_a − P_b`, per component.
    pub gap: [Series; 3],
    /// `n̂_a × n̂_b`, per component (its length is the normal-angle
    /// sine).
    pub normal_cross: [Series; 3],
    /// `II_a − II_b` on the shared frame `{e1, e2}`: components
    /// `(11, 12, 22)`.
    pub curvature: [Series; 3],
    /// `|S_u × S_v|` on each side (Jacobian magnitude).
    pub jacobian: [Iv; 2],
    /// `n_a · n_b` (hemisphere of side b's normal).
    pub normal_dot: Iv,
}

fn components(v: SeriesVec) -> [Series; 3] {
    [v.x, v.y, v.z]
}

/// The shape operator in a given tangent direction pair, by least
/// squares chart coordinates `w = (JᵀJ)⁻¹Jᵀe` and `II(e, e′) =
/// wᵀ H w′`, `H` the Hessian dotted with the unit normal.
fn second_form(side: &SideState, unit_normal: SeriesVec, e: SeriesVec, f: SeriesVec) -> Series {
    let (su, sv) = (side.su, side.sv);
    let (ee, ff, gg) = (su.dot(su), su.dot(sv), sv.dot(sv));
    let det = ee * gg - ff * ff;
    let chart = |d: SeriesVec| {
        let (r1, r2) = (d.dot(su), d.dot(sv));
        ((r1 * gg - r2 * ff).div(det), (r2 * ee - r1 * ff).div(det))
    };
    let (w1, w2) = (chart(e), chart(f));
    let (l, m, n) = (
        side.suu.dot(unit_normal),
        side.suv.dot(unit_normal),
        side.svv.dot(unit_normal),
    );
    w1.0 * l * w2.0 + (w1.0 * m * w2.1 + w1.1 * m * w2.0) + w1.1 * n * w2.1
}

/// Quantities comparing side `a` with side `b` at the same s-set.
/// `flip_b` negates side b's normal into side a's hemisphere (chosen
/// once per seam interval from its midpoint).
pub fn seam_quantities(a: &SideState, b: &SideState, flip_b: bool) -> SeamQuantities {
    let n_a = a.su.cross(a.sv);
    let n_b = b.su.cross(b.sv);
    let (len_a, len_b) = (n_a.norm(), n_b.norm());
    let sign = Series::constant(Iv::point(if flip_b { -1.0 } else { 1.0 }));
    let unit_a = n_a.scale(Series::constant(Iv::point(1.0)).div(len_a));
    let unit_b = n_b.scale(sign.div(len_b));
    let e1 = a
        .tangent
        .scale(Series::constant(Iv::point(1.0)).div(a.tangent.norm()));
    let e2 = unit_a.cross(e1);
    let form = |side: &SideState, unit: SeriesVec| {
        [
            second_form(side, unit, e1, e1),
            second_form(side, unit, e1, e2),
            second_form(side, unit, e2, e2),
        ]
    };
    let (ii_a, ii_b) = (form(a, unit_a), form(b, unit_b));
    SeamQuantities {
        gap: components(a.p - b.p),
        normal_cross: components(unit_a.cross(unit_b)),
        curvature: [ii_a[0] - ii_b[0], ii_a[1] - ii_b[1], ii_a[2] - ii_b[2]],
        jacobian: [len_a.value(), len_b.value()],
        normal_dot: n_a.dot(n_b).value(),
    }
}
