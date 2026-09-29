//! Sampled residuals of an iterate: rim against the neighbours, spokes
//! between sectors, and regularity. These drive refinement and fill the
//! report; the caller's certifier is the acceptance instrument.

use super::jets::{NeighbourJet, form_from_partials, form_jump};
use super::space::{SectorSpace, apply3};
use super::{RimHalf, RimResidual};
use crate::geometry::Vec3;

/// The cap's local jet at one parameter of one sector.
pub(super) struct CapJet {
    pub p: Vec3,
    pub su: Vec3,
    pub sv: Vec3,
    pub suu: Vec3,
    pub suv: Vec3,
    pub svv: Vec3,
}

/// Evaluate the cap jet of `sector` at `(u, v)`.
pub(super) fn cap_jet(
    space: &SectorSpace,
    coords: &[Vec<f64>; 3],
    sector: usize,
    u: f64,
    v: f64,
) -> CapJet {
    let bu = space.basis(u);
    let bv = space.basis(v);
    let d = |a, b| apply3(&space.functional_from(sector, &bu, &bv, a, b), coords);
    CapJet {
        p: d(0, 0),
        su: d(1, 0),
        sv: d(0, 1),
        suu: d(2, 0),
        suv: d(1, 1),
        svv: d(0, 2),
    }
}

/// One rim sample: where it sits and what it is compared with.
pub(super) struct RimSample {
    pub side: usize,
    pub half: RimHalf,
    /// Sector owning this half.
    pub sector: usize,
    /// Rim parameter in [0, 1] (`v` on a trailing half, `u` on a
    /// leading half), measured from `m_s`.
    pub x: f64,
    /// Side-curve parameter.
    pub t: f64,
    /// The neighbour at the rim point (model-to-normalised scaled).
    pub jet: NeighbourJet,
}

/// Parameters `(u, v)` of a rim sample and whether the cross direction
/// is u (trailing) or v (leading).
pub(super) fn rim_uv(sample: &RimSample) -> ((f64, f64), bool) {
    match sample.half {
        RimHalf::Trailing => ((1.0, sample.x), true),
        RimHalf::Leading => ((sample.x, 1.0), false),
    }
}

/// `(g0, g1, g2)` of one rim sample.
pub(super) fn rim_residual(
    space: &SectorSpace,
    coords: &[Vec<f64>; 3],
    sample: &RimSample,
) -> (f64, f64, f64) {
    let ((u, v), _) = rim_uv(sample);
    let j = cap_jet(space, coords, sample.sector, u, v);
    let n_nb = sample.jet.normal;
    let foot = Vec3::new(sample.jet.foot.x, sample.jet.foot.y, sample.jet.foot.z);
    let g0 = (j.p - foot).length();
    let Some(n_cap) = j.su.cross(&j.sv).normalized() else {
        return (g0, 1.0, f64::INFINITY);
    };
    if n_cap.dot(&n_nb) <= 0.0 {
        return (g0, 1.0, f64::INFINITY);
    }
    let g1 = n_cap.cross(&n_nb).length();
    let tangent = if sample.half == RimHalf::Trailing {
        j.sv
    } else {
        j.su
    };
    let Some(e1) = tangent.normalized() else {
        return (g0, g1, f64::INFINITY);
    };
    let e2 = n_nb.cross(&e1);
    let g2 = match form_from_partials(j.su, j.sv, j.suu, j.suv, j.svv, n_cap) {
        Some(form) => form_jump(&form, &sample.jet.form, e1, e2),
        None => f64::INFINITY,
    };
    (g0, g1, g2)
}

/// Worst residual per side half.
pub(super) fn rim_report(
    space: &SectorSpace,
    coords: &[Vec<f64>; 3],
    samples: &[RimSample],
    n: usize,
) -> Vec<RimResidual> {
    let mut out: Vec<RimResidual> = (0..n)
        .flat_map(|side| {
            [RimHalf::Leading, RimHalf::Trailing].map(|half| RimResidual {
                side,
                half,
                parameter: f64::NAN,
                g0: 0.0,
                g1: 0.0,
                g2: 0.0,
            })
        })
        .collect();
    for s in samples {
        let (g0, g1, g2) = rim_residual(space, coords, s);
        let slot = &mut out[2 * s.side + usize::from(s.half == RimHalf::Trailing)];
        let worse = g1 > slot.g1 || g2 > slot.g2 || g0 > slot.g0 || slot.parameter.is_nan();
        slot.g0 = slot.g0.max(g0);
        slot.g1 = slot.g1.max(g1);
        slot.g2 = slot.g2.max(g2);
        if worse {
            slot.parameter = s.t;
        }
    }
    out
}

/// Worst `(g0, g1, g2)` across every spoke, `samples` per spoke.
pub(super) fn spoke_report(
    space: &SectorSpace,
    coords: &[Vec<f64>; 3],
    samples: usize,
) -> (f64, f64, f64) {
    let (mut w0, mut w1, mut w2) = (0.0f64, 0.0f64, 0.0f64);
    for s in 0..space.n {
        let prev = (s + space.n - 1) % space.n;
        for k in 0..=samples {
            // Stay off the centre, where both sides' charts meet at a
            // corner and the tangent direction is the spoke itself.
            let t = 1e-3 + (1.0 - 2e-3) * k as f64 / samples as f64;
            let a = cap_jet(space, coords, s, t, 0.0);
            let b = cap_jet(space, coords, prev, 0.0, t);
            w0 = w0.max((a.p - b.p).length());
            let (Some(na), Some(nb)) = (
                a.su.cross(&a.sv).normalized(),
                b.su.cross(&b.sv).normalized(),
            ) else {
                return (w0, 1.0, f64::INFINITY);
            };
            w1 = w1.max(if na.dot(&nb) <= 0.0 {
                1.0
            } else {
                na.cross(&nb).length()
            });
            let Some(e1) = a.su.normalized() else {
                continue;
            };
            let e2 = na.cross(&e1);
            if let (Some(fa), Some(fb)) = (
                form_from_partials(a.su, a.sv, a.suu, a.suv, a.svv, na),
                form_from_partials(b.su, b.sv, b.suu, b.suv, b.svv, na),
            ) {
                w2 = w2.max(form_jump(&fa, &fb, e1, e2));
            }
        }
    }
    (w0, w1, w2)
}

/// The worst regularity sample of a grid: `min |S_u × S_v| / mean`
/// (0 on a normal flip between neighbouring grid points), its sector
/// and parameters.
#[derive(Debug, Clone, Copy)]
pub(super) struct Regularity {
    pub ratio: f64,
    pub sector: usize,
    pub uv: (f64, f64),
}

/// Regularity over a per-span grid: `per_span` evenly spaced points in
/// every knot span of each direction, plus the ends. A uniform grid is
/// blind to the thin spans graded refinement puts beside the corners
/// (measured: a fold at v = 0.958 inside the span [0.9375, 0.96875]
/// passed a 16-point uniform grid).
pub(super) fn regularity(
    space: &SectorSpace,
    coords: &[Vec<f64>; 3],
    per_span: usize,
) -> Regularity {
    let mut worst = Regularity {
        ratio: f64::INFINITY,
        sector: 0,
        uv: (0.0, 0.0),
    };
    let mut params: Vec<f64> = space
        .spans()
        .iter()
        .flat_map(|&(a, b)| (0..per_span).map(move |i| a + (b - a) * i as f64 / per_span as f64))
        .collect();
    params.push(1.0);
    let grid = params.len() - 1;
    let at = |i: usize| params[i];
    for s in 0..space.n {
        let mut values = Vec::with_capacity((grid + 1) * (grid + 1));
        let mut normals = Vec::with_capacity((grid + 1) * (grid + 1));
        for i in 0..=grid {
            for j in 0..=grid {
                let jt = cap_jet(space, coords, s, at(i), at(j));
                let cross = jt.su.cross(&jt.sv);
                values.push(cross.length());
                normals.push(cross.normalized());
            }
        }
        let mean = values.iter().sum::<f64>() / values.len() as f64;
        let idx = |i: usize, j: usize| i * (grid + 1) + j;
        for i in 0..=grid {
            for j in 0..=grid {
                // A fold reverses the normal between neighbouring grid
                // points even when no sample lands on the zero of the
                // Jacobian.
                let flipped = match normals[idx(i, j)] {
                    None => true,
                    Some(h) => [(i + 1, j), (i, j + 1)]
                        .into_iter()
                        .filter(|&(a, b)| a <= grid && b <= grid)
                        .any(|(a, b)| normals[idx(a, b)].is_none_or(|o| h.dot(&o) <= 0.0)),
                };
                let ratio = if flipped || mean.is_nan() || mean <= 0.0 {
                    0.0
                } else {
                    values[idx(i, j)] / mean
                };
                if ratio < worst.ratio {
                    worst = Regularity {
                        ratio,
                        sector: s,
                        uv: (at(i), at(j)),
                    };
                }
            }
        }
    }
    worst
}
