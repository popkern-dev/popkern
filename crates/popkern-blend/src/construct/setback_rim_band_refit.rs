//! Junction-continuous rim band refit (PK-584): after the G0 landing,
//! re-fit every rim row on its support with the end bands free, and
//! couple the band cells across fan junctions so the rim stays C1 (and
//! keeps its second-order relation) through each junction point.
//!
//! Why: the G0 landing (`setback_rim_reconcile`) translates the three
//! CPs at each end as one block (the cross-seam relations it must not
//! break are translation-invariant), and fits only k = 3..6. On a
//! variable-radius strip the deviation field is not constant across the
//! band. The landed row then keeps a ~3e-8 residual whose slope tilts
//! the rim's own tangent 1e-6..1.2e-5 off the support along the whole
//! row. The cap normal is T × C and the jet pass only corrects C, so
//! that tilt was the variable-radius cap-strip G1 floor (saddle
//! 1.1-1.75e-5; constant 8e-8). An unconstrained degree-9 fit of the
//! same field reaches 1e-13.
//!
//! Construction: only row 0 moves; the jet and curvature passes re-true
//! the cross jets afterwards. The end point P stays where the landing
//! pinned it.
//!
//! - Two rim sides from different pieces meeting at P (a fan junction),
//!   on one support, with anti-collinear neighbour vectors
//!   a = A₁ − P ∥ b = P − B₁, a = λb: ONE shared e₁ moves A₁ by e₁ and
//!   B₁ by −e₁/λ, so a = λb holds exactly (the spoke's G1 at the rim
//!   end). ONE shared e₂ moves A₂ by e₂ and B₂ by the displacement that
//!   preserves the measured second-order relation a₂ − λ²b₂ − μa
//!   (μ from the current geometry, a₂ = A₂ − 2A₁ + P).
//! - Every other end (a setback corner, a one-sided end, a junction
//!   whose neighbour vectors are not collinear) moves its k = 1, 2
//!   cells per side.
//!
//! Measured on the saddle (variable radius): uncoupled band freedom
//! took cap-strip to ~1e-7 but disturbed the spokes to 1e-8; the
//! coupling is what keeps them at their class.

use std::collections::BTreeMap;

use popkern_core::geometry::gregory::GregoryPatch;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::uv_projection::{project_point_to_uv, project_point_to_uv_near};
use popkern_core::geometry::{Point3, Vec3};

use super::bezier9_eval;
use super::dense_spd_solve::VecLsq;
use super::setback_corner_dual_jet::{coord_for, row_for_side};
use super::setback_rim_stations::side_support;

/// Stations per rim side in the refit.
const STATIONS: usize = 64;
/// Weight of the end-station tangent rows (end tangency held
/// near-exactly; see the row loop).
const END_TANGENT_WEIGHT: f64 = 100.0;
/// Refit passes (the projection field is mildly nonlinear).
const PASSES: usize = 3;
/// Two side ends are one junction when their end points agree to this
/// (the fan's shared net makes them bit-identical; the band allows for
/// the landing's per-side end projection).
const JUNCTION_POINT: f64 = 1.0e-9;
/// Neighbour vectors count as anti-collinear (the rim C1 through the
/// junction) below this sine.
const COLLINEAR_SINE: f64 = 1.0e-7;

/// Measured outcome of [`refit_rim_bands`].
#[derive(Debug, Clone, Copy, Default)]
pub(super) struct RimBandRefit {
    /// Rim sides refitted.
    pub(super) sides: usize,
    /// Junctions whose band cells were coupled.
    pub(super) coupled_junctions: usize,
    /// Worst support deviation of the rim rows before the refit.
    pub(super) worst_before: f64,
    /// And after it.
    pub(super) worst_after: f64,
}

/// One rim side taking part in the refit.
struct RefitSide {
    patch: usize,
    side: usize,
    support: usize,
    /// First unknown of the four interior cells k = 3..6.
    interior: usize,
}

/// How one end's band cells (the end's k = 1 and k = 2 neighbours)
/// map onto unknowns: `(unknown, coefficient)` lists for each cell.
type BandMap = [Vec<(usize, f64)>; 2];

/// Re-fit the rim rows' end bands (see module docs). `sides` is the
/// rim-side registry.
pub(super) fn refit_rim_bands(
    patches: &mut [GregoryPatch],
    candidates: &[SurfaceGeometry],
    strip_count: usize,
    sides: &[(usize, usize)],
) -> RimBandRefit {
    let mut report = RimBandRefit::default();
    let mut refit: Vec<RefitSide> = Vec::new();
    let mut next = 0usize;
    for &(patch, side) in sides {
        let Some(row0) = row_for_side(patches, patch, side, 0) else {
            continue;
        };
        let Some(support) = side_support(candidates, strip_count, &row0) else {
            continue;
        };
        refit.push(RefitSide {
            patch,
            side,
            support,
            interior: next,
        });
        next += 4;
    }
    report.sides = refit.len();
    if refit.is_empty() {
        return report;
    }
    let (bands, unknowns, coupled) = band_maps(patches, &refit, next);
    report.coupled_junctions = coupled;
    report.worst_before = worst_deviation(patches, candidates, &refit);
    for _ in 0..PASSES {
        let mut lsq = VecLsq::new(unknowns);
        for (si, rs) in refit.iter().enumerate() {
            let Some(row0) = row_for_side(patches, rs.patch, rs.side, 0) else {
                continue;
            };
            let Some(field) = deviation_field(&candidates[rs.support], &row0) else {
                continue;
            };
            for (m, (dev, normal)) in field.iter().enumerate() {
                let t = m as f64 / (STATIONS - 1) as f64;
                // Position row (the support deviation) and tangent row
                // (the rim tangent's normal component, the G1 floor this
                // pass exists for), in length units: B′/9.
                for (basis, tangent_row, target) in [
                    (bern9 as fn(usize, f64) -> f64, false, Some(*dev)),
                    (
                        bern9_speed,
                        true,
                        normal.map(|n| {
                            let tangent = super::setback_rim_stations::bezier9_derive(&row0, t);
                            n * (-n.dot(&tangent) / 9.0)
                        }),
                    ),
                ] {
                    let Some(target) = target else {
                        continue;
                    };
                    // End tangency is the rim's G1 at the corner and the
                    // neighbour rim's cross jet there: held near-exactly.
                    let w = if tangent_row && (m == 0 || m == STATIONS - 1) {
                        END_TANGENT_WEIGHT
                    } else {
                        1.0
                    };
                    let target = target * w;
                    let mut entries: Vec<(usize, f64)> = Vec::new();
                    for q in 0..4 {
                        entries.push((rs.interior + q, basis(3 + q, t)));
                    }
                    for (end, cells) in [(0usize, [1usize, 2]), (1, [8, 7])] {
                        for (c, &k) in cells.iter().enumerate() {
                            for &(u, w) in &bands[&(si, end)][c] {
                                entries.push((u, w * basis(k, t)));
                            }
                        }
                    }
                    for e in entries.iter_mut() {
                        e.1 *= w;
                    }
                    lsq.add_vec_row(&entries, target);
                }
            }
        }
        let Some(solution) = lsq.solve() else {
            break;
        };
        for (si, rs) in refit.iter().enumerate() {
            let mut delta = [Vec3::zero(); 10];
            for q in 0..4 {
                delta[3 + q] = solution[rs.interior + q];
            }
            for (end, cells) in [(0usize, [1usize, 2]), (1, [8, 7])] {
                for (c, &k) in cells.iter().enumerate() {
                    delta[k] = bands[&(si, end)][c]
                        .iter()
                        .fold(Vec3::zero(), |acc, &(u, w)| acc + solution[u] * w);
                }
            }
            if let Some(grid) = patches[rs.patch].binonic_grid_mut() {
                for (k, d) in delta.iter().enumerate() {
                    let (i, j) = coord_for(rs.side, 0, k);
                    let g = grid[i][j];
                    grid[i][j] = Point3::new(g.x + d.x, g.y + d.y, g.z + d.z);
                }
            }
        }
    }
    report.worst_after = worst_deviation(patches, candidates, &refit);
    report
}

/// Unknown maps for every side end. Returns the maps, the total
/// unknown count, and how many junctions were coupled.
fn band_maps(
    patches: &[GregoryPatch],
    refit: &[RefitSide],
    first_free: usize,
) -> (BTreeMap<(usize, usize), BandMap>, usize, usize) {
    // End records: (side index, end, P, first neighbour, second).
    let mut ends: Vec<(usize, usize, Point3, Point3, Point3)> = Vec::new();
    for (si, rs) in refit.iter().enumerate() {
        let Some(row) = row_for_side(patches, rs.patch, rs.side, 0) else {
            continue;
        };
        ends.push((si, 0, row[0], row[1], row[2]));
        ends.push((si, 1, row[9], row[8], row[7]));
    }
    let mut maps: BTreeMap<(usize, usize), BandMap> = BTreeMap::new();
    let mut next = first_free;
    let mut coupled = 0usize;
    let mut used = vec![false; ends.len()];
    for a in 0..ends.len() {
        if used[a] {
            continue;
        }
        let partners: Vec<usize> = (0..ends.len())
            .filter(|&b| b != a && ends[b].2.distance_to(&ends[a].2) <= JUNCTION_POINT)
            .collect();
        let pair = match partners.as_slice() {
            [b] if !used[*b]
                && refit[ends[*b].0].patch != refit[ends[a].0].patch
                && refit[ends[*b].0].support == refit[ends[a].0].support =>
            {
                junction_coupling(&ends[a], &ends[*b])
            }
            _ => None,
        };
        match (pair, partners.as_slice()) {
            (Some((lambda, mu)), [b]) => {
                let (e1, e2) = (next, next + 1);
                next += 2;
                maps.insert((ends[a].0, ends[a].1), [vec![(e1, 1.0)], vec![(e2, 1.0)]]);
                let l2 = lambda * lambda;
                maps.insert(
                    (ends[*b].0, ends[*b].1),
                    [
                        vec![(e1, -1.0 / lambda)],
                        vec![(e2, 1.0 / l2), (e1, -(2.0 + mu) / l2 - 2.0 / lambda)],
                    ],
                );
                used[a] = true;
                used[*b] = true;
                coupled += 1;
            }
            _ => {
                maps.insert(
                    (ends[a].0, ends[a].1),
                    [vec![(next, 1.0)], vec![(next + 1, 1.0)]],
                );
                next += 2;
                used[a] = true;
            }
        }
    }
    (maps, next, coupled)
}

/// `(λ, μ)` of a junction whose rim is C1 through P: `a = A₁ − P`,
/// `b = P − B₁` anti-collinear with `a = λb`, and μ the tangential part
/// of `a₂ − λ²b₂`. `None` when the rim is not C1 there.
fn junction_coupling(
    a: &(usize, usize, Point3, Point3, Point3),
    b: &(usize, usize, Point3, Point3, Point3),
) -> Option<(f64, f64)> {
    let p = a.2;
    let va = a.3 - p;
    let vb = p - b.3;
    let (la, lb) = (va.length(), vb.length());
    if la < 1e-15 || lb < 1e-15 || va.dot(&vb) <= 0.0 {
        return None;
    }
    if va.cross(&vb).length() > COLLINEAR_SINE * la * lb {
        return None;
    }
    let lambda = la / lb;
    let a2 = (a.4 - a.3) - (a.3 - p);
    let b2 = (b.4 - b.3) - (b.3 - p);
    let mu = (a2 - b2 * (lambda * lambda)).dot(&va) / (la * la);
    Some((lambda, mu))
}

/// Support deviation of a rim row at the refit stations, with the
/// support normal at each foot point, projected by continuation (each
/// station seeded from the previous one's uv).
fn deviation_field(
    support: &SurfaceGeometry,
    row0: &[Point3; 10],
) -> Option<Vec<(Vec3, Option<Vec3>)>> {
    let mut out = Vec::with_capacity(STATIONS);
    let mut prev: Option<(f64, f64)> = None;
    for m in 0..STATIONS {
        let p = bezier9_eval(row0, m as f64 / (STATIONS - 1) as f64);
        let (u, v) = match prev {
            Some((u, v)) => project_point_to_uv_near(support, &p, u, v)
                .or_else(|| project_point_to_uv(support, &p))?,
            None => project_point_to_uv(support, &p)?,
        };
        prev = Some((u, v));
        let h = 1e-6;
        let du = support.point_at(u + h, v) - support.point_at(u - h, v);
        let dv = support.point_at(u, v + h) - support.point_at(u, v - h);
        out.push((support.point_at(u, v) - p, du.cross(&dv).normalized()));
    }
    Some(out)
}

fn worst_deviation(
    patches: &[GregoryPatch],
    candidates: &[SurfaceGeometry],
    refit: &[RefitSide],
) -> f64 {
    refit
        .iter()
        .filter_map(|rs| {
            let row0 = row_for_side(patches, rs.patch, rs.side, 0)?;
            deviation_field(&candidates[rs.support], &row0)
        })
        .flatten()
        .map(|(d, _)| d.length())
        .fold(0.0, f64::max)
}

fn bern9(k: usize, t: f64) -> f64 {
    const C: [f64; 10] = [1.0, 9.0, 36.0, 84.0, 126.0, 126.0, 84.0, 36.0, 9.0, 1.0];
    C[k] * t.powi(k as i32) * (1.0 - t).powi((9 - k) as i32)
}

/// `B′ₖ(t) / 9` of the degree-9 Bernstein basis: the derivative in the
/// same length units as the position rows.
fn bern9_speed(k: usize, t: f64) -> f64 {
    const C: [f64; 9] = [1.0, 8.0, 28.0, 56.0, 70.0, 56.0, 28.0, 8.0, 1.0];
    let b8 = |j: usize| C[j] * t.powi(j as i32) * (1.0 - t).powi((8 - j) as i32);
    let lower = if k >= 1 { b8(k - 1) } else { 0.0 };
    let upper = if k <= 8 { b8(k) } else { 0.0 };
    lower - upper
}
