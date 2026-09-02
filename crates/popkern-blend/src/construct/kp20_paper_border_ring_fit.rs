//! Ring-wide cross-field fit for the setback border (PK-400).
//!
//! The setback fan's cross rows (cap-side jet `row1 − row0` and
//! cross-curvature `row2 − 2·row1 + row0` of every degree-6 border
//! piece) must carry data SHARED between neighbouring pieces: one
//! value / step / second difference at each spoke point `m` (the ring
//! is C² across the spoke by construction) and one KP20 corner block
//! at each sector corner `c` (both pieces reproduce the same 3×3
//! corner block, so the assembly's overlap average is a no-op).
//! Fitting each piece separately and PINNING the shared data to the
//! mean of the two pieces' independent claims left one free
//! coefficient per row — the c-quarter of every half then followed the
//! pins, not the sampled field, and the cap left the strips' tangent
//! planes by 2–4e-3 mid-quarter on every n-gonal pyramid fixture (the
//! cube's straight rims hide it: a constant tangent plane makes any
//! polynomial through in-plane pins in-plane).
//!
//! This module fits the whole ring as ONE linear least-squares system
//! in which the shared junction data are the unknowns: the spoke
//! junctions, the corner blocks' four interior slots, and each row's
//! middle coefficient. Every row coefficient is affine in those
//! unknowns, so the C²-at-`m` and shared-block invariants hold exactly
//! whatever the solve returns. The objective is weighted: tangency
//! rows (`n·d1 = 0` at every sampled rim point, `n` the strip's unit
//! normal there) dominate, then the jet field's in-plane shape, then
//! the curvature field — on G¹ rolling-ball input the two strips'
//! curvature claims at `c` are incompatible (each section carries
//! `1/r` where the other is straight) and the solve parks that
//! disagreement in the curvature rows instead of tilting the jets.

use popkern_core::geometry::gregory::bernstein6;
use popkern_core::geometry::{Point3, Vec3};

use super::ring_solve::System;
use super::{CornerJet, CrossFieldSamples};

/// Weight of the tangency rows (normal tilt `n·d1 / scale`).
const TANGENCY_WEIGHT: f64 = 1.0;
/// Weight of the jet field's vector fidelity (its in-plane shape),
/// relative to the sampled jet magnitude.
const JET_FIDELITY_WEIGHT: f64 = 0.1;
/// Weight of the cross-curvature field's fidelity (relative to the
/// jet magnitude): the honest home of a corner's G² incompatibility,
/// so it yields first. A curvature disagreement of relative size ε
/// between the two pieces at a corner leaks ~ε·w² into the jets
/// through the shared block; at 0.01 a 20 % disagreement costs 2e-5
/// of tilt (unit-tested), and G² sweeps carry no disagreement.
const CURVATURE_FIDELITY_WEIGHT: f64 = 0.01;

/// One half-arc's input to the ring fit.
pub(super) struct RingPiece<'a> {
    /// The rim row (exact strip section) of this half.
    pub row0: &'a [Point3; 7],
    /// Sampled cross fields at true rim points.
    pub samples: &'a CrossFieldSamples,
    /// Tangency stations `(t, unit strip normal, jet scale)`: the
    /// scale is the sampled jet's component PERPENDICULAR to the rim
    /// there, so the residual `n·d1 / scale` is the cap's normal tilt
    /// (a jet skewed toward the rim — the blunt pyramids' 159° cap
    /// corners — amplifies its out-of-plane fraction by 1/sin of the
    /// skew; measured 7.2e-4 → 2.0e-3 next to c). Denser than the
    /// march samples: normals are one surface evaluation each.
    pub tangency: Vec<(f64, Vec3, f64)>,
}

/// Corner data at sector corner `c_k` between arc `k`'s second half
/// (`A`, ending at `c`) and arc `k+1`'s first half (`B`, starting
/// there): `Some((jet_a, jet_b))` shares one corner block; `None`
/// (gapped or reflex corner) leaves both ends' coefficients free with
/// the end values taken from the sampled field.
pub(super) type RingCorner = Option<(CornerJet, CornerJet)>;

/// Fitted cross rows per arc, per half: `(d1, d2)`.
pub(super) type RingRows = [([Vec3; 7], [Vec3; 7]); 2];

/// A row coefficient as an affine function of the Vec3 unknowns.
#[derive(Clone, Debug)]
pub(super) struct Affine {
    pub(super) constant: Vec3,
    pub(super) terms: Vec<(usize, f64)>,
}

impl Affine {
    fn constant(v: Vec3) -> Self {
        Self {
            constant: v,
            terms: Vec::new(),
        }
    }
    fn unknown(index: usize) -> Self {
        Self {
            constant: Vec3::zero(),
            terms: vec![(index, 1.0)],
        }
    }
    fn scaled(&self, s: f64) -> Self {
        Self {
            constant: self.constant * s,
            terms: self.terms.iter().map(|&(i, w)| (i, w * s)).collect(),
        }
    }
    fn plus(&self, other: &Self) -> Self {
        let mut terms = self.terms.clone();
        terms.extend(other.terms.iter().copied());
        Self {
            constant: self.constant + other.constant,
            terms,
        }
    }
    fn plus_const(&self, v: Vec3) -> Self {
        Self {
            constant: self.constant + v,
            terms: self.terms.clone(),
        }
    }
    fn evaluate(&self, x: &[Vec3]) -> Vec3 {
        let mut v = self.constant;
        for &(i, w) in &self.terms {
            v = v + x[i] * w;
        }
        v
    }
}

fn combine(coefficients: &[Affine; 7], weights: [f64; 7]) -> Affine {
    let mut acc = Affine::constant(Vec3::zero());
    for (c, w) in coefficients.iter().zip(weights) {
        if w != 0.0 {
            acc = acc.plus(&c.scaled(w));
        }
    }
    acc
}

fn vec_of(p: Point3) -> Vec3 {
    Vec3::new(p.x, p.y, p.z)
}

/// De Casteljau c-quarter of a rim row: the last four coefficients
/// of the half at `t = ½` for a row ending at `c` (`at_end`), the
/// first four otherwise. Returns the quarter's 7 coefficients.
fn quarter_of(row0: &[Point3; 7], at_end: bool) -> [Point3; 7] {
    let (first, second) = popkern_core::geometry::gregory::split_sextic_row(row0);
    if at_end { second } else { first }
}

/// Coefficient slots of one piece's two cross rows, as affine
/// expressions in the ring unknowns.
struct PieceRows {
    d1: [Affine; 7],
    d2: [Affine; 7],
}

fn sample_at(samples: &[(f64, Vec3)], target: f64) -> Vec3 {
    samples
        .iter()
        .min_by(|a, b| {
            (a.0 - target)
                .abs()
                .partial_cmp(&(b.0 - target).abs())
                .unwrap_or(core::cmp::Ordering::Equal)
        })
        .map(|s| s.1)
        .unwrap_or(Vec3::zero())
}

/// Fit the ring's cross rows. `arcs[k]` holds arc `k`'s two halves
/// (`[c_{k-1} → m_k]`, `[m_k → c_k]`); `corners[k]` the data at `c_k`.
/// Returns the fitted `(d1, d2)` rows per arc and half, or `None` when
/// the normal equations are singular (degenerate sampling).
pub(super) fn fit_ring_cross_fields(
    arcs: &[[RingPiece<'_>; 2]],
    corners: &[RingCorner],
) -> Option<Vec<RingRows>> {
    let n = arcs.len();
    if n == 0 || corners.len() != n {
        return None;
    }
    let mut next_unknown = 0_usize;
    let mut alloc = |count: usize| -> usize {
        let base = next_unknown;
        next_unknown += count;
        base
    };
    // Spoke junctions per arc: jet (value, step, second) and curvature.
    let junction: Vec<(usize, usize)> = (0..n).map(|_| (alloc(3), alloc(3))).collect();
    // Corner blocks per corner (four interior slots), when shared.
    let block: Vec<Option<usize>> = corners
        .iter()
        .map(|c| c.as_ref().map(|_| alloc(4)))
        .collect();
    // Free end triples at unshared corners: (A: 4, 5 jet + 4, 5 curv),
    // (B: 1, 2 jet + 1, 2 curv).
    let free_end: Vec<Option<usize>> = corners
        .iter()
        .map(|c| c.is_none().then(|| alloc(8)))
        .collect();
    // Middle coefficient (index 3) per piece, jet and curvature.
    let middle: Vec<[usize; 2]> = (0..n).map(|_| [alloc(2), alloc(2)]).collect();

    let unknowns = next_unknown;
    let mut pieces: Vec<[PieceRows; 2]> = Vec::with_capacity(n);
    for (k, arc) in arcs.iter().enumerate() {
        let prev = (k + n - 1) % n;
        let mut rows = [
            PieceRows {
                d1: core::array::from_fn(|_| Affine::constant(Vec3::zero())),
                d2: core::array::from_fn(|_| Affine::constant(Vec3::zero())),
            },
            PieceRows {
                d1: core::array::from_fn(|_| Affine::constant(Vec3::zero())),
                d2: core::array::from_fn(|_| Affine::constant(Vec3::zero())),
            },
        ];
        // Spoke junction at m_k: half 0 ends with it, half 1 starts.
        let (j1, j2) = junction[k];
        for (row, base) in [(&mut rows[0].d1, j1), (&mut rows[0].d2, j2)] {
            let (value, step, second) = (
                Affine::unknown(base),
                Affine::unknown(base + 1),
                Affine::unknown(base + 2),
            );
            row[6] = value.clone();
            row[5] = value.plus(&step.scaled(-1.0));
            row[4] = value.plus(&step.scaled(-2.0)).plus(&second);
        }
        for (row, base) in [(&mut rows[1].d1, j1), (&mut rows[1].d2, j2)] {
            let (value, step, second) = (
                Affine::unknown(base),
                Affine::unknown(base + 1),
                Affine::unknown(base + 2),
            );
            row[0] = value.clone();
            row[1] = value.plus(&step);
            row[2] = value.plus(&step.scaled(2.0)).plus(&second);
        }
        // Middle coefficients.
        for (half, row_set) in rows.iter_mut().enumerate() {
            row_set.d1[3] = Affine::unknown(middle[k][half]);
            row_set.d2[3] = Affine::unknown(middle[k][half] + 1);
        }
        // Corner end of half 1 (piece A at c_k).
        match (&corners[k], block[k], free_end[k]) {
            (Some((jet_a, _)), Some(base), _) => {
                let q = quarter_of(arc[1].row0, true);
                let (m55, m45, m54, m44) = (
                    Affine::unknown(base),
                    Affine::unknown(base + 1),
                    Affine::unknown(base + 2),
                    Affine::unknown(base + 3),
                );
                // Quarter coefficients from the block slots.
                let a1_5 = m55.scaled(2.0).plus_const(vec_of(q[5]) * -2.0);
                let a1_4 = m45.scaled(2.0).plus_const(vec_of(q[4]) * -2.0);
                let a2_5 = m54
                    .scaled(5.0)
                    .plus(&a1_5.scaled(-5.0))
                    .plus_const(vec_of(q[5]) * -5.0);
                let a2_4 = m44
                    .scaled(5.0)
                    .plus(&a1_4.scaled(-5.0))
                    .plus_const(vec_of(q[4]) * -5.0);
                lift_end(&mut rows[1].d1, Affine::constant(jet_a.d1), &a1_5, &a1_4);
                lift_end(&mut rows[1].d2, Affine::constant(jet_a.d2), &a2_5, &a2_4);
            }
            (None, _, Some(base)) => {
                rows[1].d1[6] = Affine::constant(sample_at(&arc[1].samples.jets, 1.0));
                rows[1].d2[6] = Affine::constant(sample_at(&arc[1].samples.curvatures, 1.0));
                rows[1].d1[5] = Affine::unknown(base);
                rows[1].d1[4] = Affine::unknown(base + 1);
                rows[1].d2[5] = Affine::unknown(base + 2);
                rows[1].d2[4] = Affine::unknown(base + 3);
            }
            _ => return None,
        }
        // Corner start of half 0 (piece B at c_{k-1}).
        match (&corners[prev], block[prev], free_end[prev]) {
            (Some((_, jet_b)), Some(base), _) => {
                let q = quarter_of(arc[0].row0, false);
                let (m55, m45, m54, m44) = (
                    Affine::unknown(base),
                    Affine::unknown(base + 1),
                    Affine::unknown(base + 2),
                    Affine::unknown(base + 3),
                );
                let b1_1 = m55.scaled(2.0).plus_const(vec_of(q[1]) * -2.0);
                let b1_2 = m54.scaled(2.0).plus_const(vec_of(q[2]) * -2.0);
                let b2_1 = m45
                    .scaled(5.0)
                    .plus(&b1_1.scaled(-5.0))
                    .plus_const(vec_of(q[1]) * -5.0);
                let b2_2 = m44
                    .scaled(5.0)
                    .plus(&b1_2.scaled(-5.0))
                    .plus_const(vec_of(q[2]) * -5.0);
                lift_start(&mut rows[0].d1, Affine::constant(jet_b.d1), &b1_1, &b1_2);
                lift_start(&mut rows[0].d2, Affine::constant(jet_b.d2), &b2_1, &b2_2);
            }
            (None, _, Some(base)) => {
                rows[0].d1[0] = Affine::constant(sample_at(&arc[0].samples.jets, 0.0));
                rows[0].d2[0] = Affine::constant(sample_at(&arc[0].samples.curvatures, 0.0));
                rows[0].d1[1] = Affine::unknown(base + 4);
                rows[0].d1[2] = Affine::unknown(base + 5);
                rows[0].d2[1] = Affine::unknown(base + 6);
                rows[0].d2[2] = Affine::unknown(base + 7);
            }
            _ => return None,
        }
        pieces.push(rows);
    }

    // Residuals.
    let mut system = System::new(unknowns);
    for (arc, rows) in arcs.iter().zip(&pieces) {
        for (piece, row_set) in arc.iter().zip(rows) {
            // Field rows are scale-free: relative to the sampled jet
            // magnitude, so the weights compare like with like.
            let jet_scale = piece
                .samples
                .jets
                .iter()
                .map(|(_, j)| j.length())
                .fold(0.0_f64, f64::max)
                .max(f64::MIN_POSITIVE);
            for &(t, normal, scale) in &piece.tangency {
                let value = combine(&row_set.d1, bernstein6(t));
                system.add_normal(
                    &value,
                    normal,
                    TANGENCY_WEIGHT / scale.max(1.0e-3 * jet_scale),
                );
            }
            for &(t, jet) in &piece.samples.jets {
                let value = combine(&row_set.d1, bernstein6(t));
                system.add_vector(&value, jet, JET_FIDELITY_WEIGHT / jet_scale);
            }
            for &(t, curvature) in &piece.samples.curvatures {
                let value = combine(&row_set.d2, bernstein6(t));
                system.add_vector(&value, curvature, CURVATURE_FIDELITY_WEIGHT / jet_scale);
            }
        }
    }
    let x = system.solve()?;
    Some(
        pieces
            .iter()
            .map(|rows| {
                core::array::from_fn(|half| {
                    (
                        core::array::from_fn(|i| rows[half].d1[i].evaluate(&x)),
                        core::array::from_fn(|i| rows[half].d2[i].evaluate(&x)),
                    )
                })
            })
            .collect(),
    )
}

/// Full-row end coefficients reproducing given c-quarter coefficients
/// for a row ENDING at `c` (inverse de Casteljau half: `q6 = b6`,
/// `q5 = (b5 + b6)/2`, `q4 = (b4 + 2b5 + b6)/4`).
fn lift_end(row: &mut [Affine; 7], b6: Affine, q5: &Affine, q4: &Affine) {
    let b5 = q5.scaled(2.0).plus(&b6.scaled(-1.0));
    let b4 = q4.scaled(4.0).plus(&b5.scaled(-2.0)).plus(&b6.scaled(-1.0));
    row[4] = b4;
    row[5] = b5;
    row[6] = b6;
}

/// Full-row start coefficients reproducing given c-quarter
/// coefficients for a row STARTING at `c`.
fn lift_start(row: &mut [Affine; 7], b0: Affine, q1: &Affine, q2: &Affine) {
    let b1 = q1.scaled(2.0).plus(&b0.scaled(-1.0));
    let b2 = q2.scaled(4.0).plus(&b1.scaled(-2.0)).plus(&b0.scaled(-1.0));
    row[0] = b0;
    row[1] = b1;
    row[2] = b2;
}

/// Tangency stations for one half: the strip's unit normal at the rim
/// point from the rim's ANALYTIC uv (`uv_at(span_t)` is the arc's
/// section path on the strip from corner to corner, this half
/// covering `[half/2, (half+1)/2]` — one surface evaluation, no
/// projection; the projection serves only where the fitted half is
/// not that path, i.e. the cubic fallback), scaled by the sampled
/// jet's component perpendicular to the rim, interpolated between the
/// march samples.
pub(super) fn tangency_stations(
    surface: &popkern_core::geometry::SurfaceGeometry,
    uv_at: &dyn Fn(f64) -> (f64, f64),
    half: usize,
    row0: &[Point3; 7],
    samples: &CrossFieldSamples,
) -> Vec<(f64, Vec3, f64)> {
    use popkern_core::geometry::gregory::eval_sextic_row;
    use popkern_core::geometry::surface::SurfaceEval;
    use popkern_core::geometry::uv_projection::project_point_to_uv;
    const TANGENCY_STATIONS: usize = 97;
    let perp: Vec<(f64, f64)> = samples
        .jets
        .iter()
        .map(|&(t, jet)| {
            let h = 1.0e-4;
            let tangent =
                eval_sextic_row(row0, (t + h).min(1.0)) - eval_sextic_row(row0, (t - h).max(0.0));
            let scale = match tangent.normalized() {
                Some(tau) => (jet - tau * jet.dot(&tau)).length(),
                None => jet.length(),
            };
            (t, scale)
        })
        .collect();
    let scale_at = |t: f64| -> f64 {
        let upper = perp
            .iter()
            .position(|(pt, _)| *pt >= t)
            .unwrap_or(perp.len().saturating_sub(1));
        if upper == 0 {
            return perp.first().map(|p| p.1).unwrap_or(1.0);
        }
        let (t0, s0) = perp[upper - 1];
        let (t1, s1) = perp[upper];
        let f = ((t - t0) / (t1 - t0).max(1e-15)).clamp(0.0, 1.0);
        s0 + (s1 - s0) * f
    };
    (0..TANGENCY_STATIONS)
        .map(|i| {
            let t = i as f64 / (TANGENCY_STATIONS - 1) as f64;
            let p = eval_sextic_row(row0, t);
            let (u, v) = uv_at((half as f64 + t) * 0.5);
            let normal = if surface.point_at(u, v).distance_to(&p) < 1.0e-9 {
                surface.normal_at(u, v).normalized()
            } else {
                project_point_to_uv(surface, &p)
                    .and_then(|(u, v)| surface.normal_at(u, v).normalized())
            };
            (t, normal.unwrap_or(Vec3::zero()), scale_at(t))
        })
        .collect()
}

#[cfg(test)]
#[path = "kp20_paper_border_ring_fit_tests.rs"]
mod tests;
