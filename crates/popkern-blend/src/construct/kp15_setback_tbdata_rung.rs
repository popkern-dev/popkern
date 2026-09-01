//! KP15 boundary-faithful `TBdata` rung of the polar backstop
//! (PK-090): assemble the KP15 bi-6 G² cap directly from the
//! FIXED-CORNER setback border — the same two-piece-per-sector border
//! the KP20 paper candidate consumes.
//!
//! Why this closes where the stage-4 record declared a trade-off:
//! the old TBdata input was the crude spoke-proxy border, whose
//! cross-sector junctions carried inconsistent jets (eq (2) ≈ 1.97).
//! The setback border's pieces join BETWEEN sectors at the arc
//! midpoints `m` — single-sweep data, C∞-consistent — so the KP15
//! ring-consistency requirement is met by construction wherever the
//! trim corners are exact. Measured (PK-090 probe, 2026-07-15):
//! pentagon5_uniform closes boundary interpolation EXACTLY with
//! eq (2) = 1.1e-14 and eq (3) = 1.6e-9; hexagon6_variable closes
//! eq (2) = 1.3e-14 with eq (3) opening only to the border's own
//! skew-corner inconsistency (`sal` = 9.1e-4 → eq (3) = 1.1e-2).
//!
//! All-or-nothing: the rung commits only a fold-free cap whose
//! eq (2)/(3) residuals close at the paper level up to a measured
//! multiple of the input inconsistency `sal` (and never above the
//! 5e-2 commercial ceiling); otherwise it declines with the measured
//! reason and the shared-net rung / F_κ ladder proceed unchanged.
//! Folds are a free-set placement property (boundary rows are
//! assigned before the free set is consulted), so a folded seed is
//! retried through deterministic seed variants and a bounded
//! coordinate descent over the paper free set — both boundary-neutral.

use popkern_core::geometry::Point3;
use popkern_core::geometry::gregory::{
    Bi6DeepInteriorFree, Bi6Reparams, GregoryPatch, RectangularGregoryBiNonic,
    SectorTensorBorderBi3, build_bi6_cap_from_tensor_border, elevate_bi6_grid_to_bi9,
};

use crate::fillet_surface::FilletSurfaceResult;

const SELF_INT_TOL: f64 = 1e-6;
/// Bounded descent budget: rebuilds are linear solves, but the fold
/// detector dominates; cap the search rather than the apex latency.
const MAX_DESCENT_BUILDS: usize = 240;

/// Try the TBdata rung. `Ok` carries the slot-oriented bi-9 patches
/// (two-stripe sub-patch frame, one per ring slot); `Err` carries the
/// measured decline reason.
pub(super) fn try_kp15_setback_tbdata_rung(
    central_vertex: Point3,
    outer_corners: &[Point3],
    fillet_surfaces: &[FilletSurfaceResult],
) -> Result<Vec<GregoryPatch>, String> {
    let n = outer_corners.len();
    if n < 5 {
        return Err(format!("TBdata rung is n >= 5; got {n}"));
    }
    if fillet_surfaces.len() != n {
        return Err(format!(
            "TBdata rung needs one fillet per sector ({} != {n})",
            fillet_surfaces.len()
        ));
    }
    let (border, _extract) =
        crate::construct::kp20_paper_border_extract::kp20_paper_border_from_fillet_setback(
            fillet_surfaces,
            central_vertex,
        )
        .map_err(|e| format!("setback border extraction: {e}"))?;
    if border.sectors.len() != n {
        return Err(format!(
            "border sector count {} != n = {n}",
            border.sectors.len()
        ));
    }

    // ---- ring alignment: border sector k ↔ ring slot s -------------
    // Sector k's outer corner is piece0's far end (rows[0][3] of
    // [m^k → c_k]). Match the corner cycle against the ring's outer
    // corners; forward rotation only — a reversed match means the
    // border winds against the ring and the rung declines (fail
    // closed) rather than committing mirrored geometry.
    let corner_of = |k: usize| border.sectors[k][0].rows[0][6];
    let spacing = (0..n)
        .map(|s| outer_corners[s].distance_to(&outer_corners[(s + 1) % n]))
        .fold(f64::INFINITY, f64::min);
    let band = 0.25 * spacing;
    let rotation = (0..n)
        .find(|&r| (0..n).all(|s| corner_of((s + r) % n).distance_to(&outer_corners[s]) < band));
    let Some(rotation) = rotation else {
        let d0 = (0..n)
            .map(|r| corner_of(r).distance_to(&outer_corners[0]))
            .fold(f64::INFINITY, f64::min);
        return Err(format!(
            "border/ring corner cycles do not align forward (nearest corner-0 match {d0:.3e}, \
             band {band:.3e}) — reversed winding declines"
        ));
    };

    // ---- TBdata edges, corner-first (probe-verified orientation) ---
    // tad(s) = piece1 of sector s (c_s → m^{s+1}), tak(s) = reversed
    // piece0 (c_s → m^s). Shared spoke data between adjacent sectors
    // then comes from one sweep's rows split at m — the consistency
    // the construction requires.
    let mut tad_src: Vec<SectorTensorBorderBi3> = Vec::with_capacity(n);
    let mut tak_src: Vec<SectorTensorBorderBi3> = Vec::with_capacity(n);
    for s in 0..n {
        let k = (s + rotation) % n;
        // Cubic reductions that keep the spoke junction C²: piece 0
        // [m → c] carries m at its start, piece 1 [c → m'] at its end.
        tak_src.push(reversed_rows(
            &border.sectors[k][0].bi3_rows_matching_end_jet(false),
        ));
        tad_src.push(border.sectors[k][1].bi3_rows_matching_end_jet(true));
    }
    let params = Bi6Reparams::for_valence(n);

    // ---- free-set seeds --------------------------------------------
    let centroid = {
        let mut acc = [0.0_f64; 3];
        let mut count = 0.0;
        for src in tad_src.iter().chain(tak_src.iter()) {
            for cp in &src[0] {
                acc[0] += cp.x;
                acc[1] += cp.y;
                acc[2] += cp.z;
                count += 1.0;
            }
        }
        Point3::new(acc[0] / count, acc[1] / count, acc[2] / count)
    };
    let star = lerp(centroid, central_vertex, 0.25);
    let scale = tad_src
        .iter()
        .chain(tak_src.iter())
        .flat_map(|src| src.iter().flatten())
        .map(|cp| cp.distance_to(&star))
        .fold(0.0_f64, f64::max)
        .max(f64::MIN_POSITIVE);

    // Seed variants: mid-cap free CPs aimed at the corner, the sector
    // arc midpoint, or their average. Cheap retries before descent.
    let corner_dir: Vec<Point3> = (0..n).map(|s| tad_src[s][0][0]).collect();
    let mid_dir: Vec<Point3> = (0..n)
        .map(|s| {
            let a = tad_src[s][0][3];
            let b = tak_src[s][0][3];
            Point3::new(0.5 * (a.x + b.x), 0.5 * (a.y + b.y), 0.5 * (a.z + b.z))
        })
        .collect();
    let avg_dir: Vec<Point3> = (0..n)
        .map(|s| lerp(corner_dir[s], mid_dir[s], 0.5))
        .collect();

    let mut builds = 0usize;
    let mut best: Option<ScoredGrids> = None; // (grids, sal, folds)
    let mut free_best: Option<Bi6DeepInteriorFree> = None;
    for targets in [&mid_dir, &avg_dir, &corner_dir] {
        let free = seed_free(star, targets, n);
        let Ok((grids, sal, folds)) =
            build_and_score(&params, &tad_src, &tak_src, &free, &mut builds)
        else {
            continue;
        };
        let better = best.as_ref().map(|(_, _, bf)| folds < *bf).unwrap_or(true);
        if better {
            best = Some((grids, sal, folds));
            free_best = Some(free);
        }
        if best.as_ref().map(|(_, _, f)| *f == 0).unwrap_or(false) {
            break;
        }
    }
    let (Some((mut grids, mut sal, mut folds)), Some(mut free)) = (best, free_best) else {
        return Err("TBdata assembly failed for every seed variant".to_string());
    };

    // ---- boundary-neutral fold descent over the free set -----------
    if folds > 0 {
        let steps = [0.12, -0.12, 0.05, -0.05];
        'outer: for _pass in 0..2 {
            for slot in 0..free_slot_count(n) {
                for &step in &steps {
                    if builds >= MAX_DESCENT_BUILDS {
                        break 'outer;
                    }
                    let mut trial = free.clone();
                    nudge_free_slot(&mut trial, slot, step * scale, star);
                    if let Ok((g, s2, f2)) =
                        build_and_score(&params, &tad_src, &tak_src, &trial, &mut builds)
                        && f2 < folds
                    {
                        grids = g;
                        sal = s2;
                        folds = f2;
                        free = trial;
                        if folds == 0 {
                            break 'outer;
                        }
                    }
                }
            }
        }
    }
    if folds > 0 {
        return Err(format!(
            "fold descent exhausted at {folds}/{n} folded sectors ({builds} builds)"
        ));
    }

    // ---- all-or-nothing gates ---------------------------------------
    // eq (2)/(3) close structurally up to the input inconsistency:
    // measured multipliers ~10x (eq 2) / ~30x (eq 3) of `sal` across
    // pentagon5 (sal 1e-10 → eq3 1.6e-9) and hexagon6_variable
    // (sal 9.1e-4 → eq3 1.1e-2). The 5e-2 commercial ceiling is
    // absolute — an input bad enough to breach it declines.
    let (eq2, eq3) = worst_eq2_eq3(&grids, &params);
    let gate_eq2 = (1.0e-6 * scale).max(10.0 * sal);
    let gate_eq3 = (1.0e-6 * scale).max(30.0 * sal);
    if eq2 > gate_eq2.min(5.0e-2) || eq3 > gate_eq3.min(5.0e-2) {
        return Err(format!(
            "eq gates open: eq2={eq2:.3e} (gate {gate_eq2:.3e}) eq3={eq3:.3e} \
             (gate {gate_eq3:.3e}, ceiling 5e-2, sal={sal:.3e})"
        ));
    }
    let mut worst_corner = 0.0_f64;
    for s in 0..n {
        worst_corner = worst_corner
            .max(grids[s][0][0].distance_to(&tad_src[s][0][0]))
            .max(grids[s][6][0].distance_to(&tad_src[s][0][3]))
            .max(grids[s][0][6].distance_to(&tak_src[s][0][3]));
    }
    if worst_corner > 1.0e-9 * scale.max(1.0) {
        return Err(format!(
            "boundary interpolation broke: worst corner {worst_corner:.3e}"
        ));
    }

    // ---- slot frame ---------------------------------------------------
    // Script grid → standard ([6-i][6-j]) puts the apex at [0][0];
    // elevate to bi-9; then the two-stripe slot u-flip
    // (grid'[i][j] = grid[9-i][j]) — same derivation as the shared-net
    // commit path.
    Ok(grids
        .iter()
        .map(|script| {
            let standard: [[Point3; 7]; 7] =
                core::array::from_fn(|i| core::array::from_fn(|j| script[6 - i][6 - j]));
            let bi9 = elevate_bi6_grid_to_bi9(standard);
            let flipped: [[Point3; 10]; 10] = core::array::from_fn(|i| bi9[9 - i]);
            GregoryPatch::RectangularBiNonic(Box::new(RectangularGregoryBiNonic { grid: flipped }))
        })
        .collect())
}

fn lerp(a: Point3, b: Point3, t: f64) -> Point3 {
    Point3::new(
        a.x + (b.x - a.x) * t,
        a.y + (b.y - a.y) * t,
        a.z + (b.z - a.z) * t,
    )
}

fn reversed_rows(rows: &[[Point3; 4]; 3]) -> [[Point3; 4]; 3] {
    core::array::from_fn(|d| core::array::from_fn(|i| rows[d][3 - i]))
}

fn seed_free(star: Point3, targets: &[Point3], n: usize) -> Bi6DeepInteriorFree {
    Bi6DeepInteriorFree {
        p00: star,
        tan: [
            lerp(star, targets[0], 0.06),
            lerp(star, targets[n - 1], 0.06),
        ],
        two: [
            lerp(star, targets[0], 0.14),
            lerp(star, targets[0], 0.10),
            lerp(star, targets[n - 1], 0.14),
        ],
        p220: lerp(star, targets[0], 0.22),
        p210: lerp(star, targets[0], 0.18),
        p23: (0..n).map(|s| lerp(star, targets[s], 0.55)).collect(),
        p33: (0..n).map(|s| lerp(star, targets[s], 0.35)).collect(),
    }
}

fn free_slot_count(n: usize) -> usize {
    // p00, tan×2, two×3, p220, p210, then per-sector p23/p33.
    8 + 2 * n
}

/// Nudge one free CP along its star-radial direction (a shape move
/// that cannot touch the boundary rows).
fn nudge_free_slot(free: &mut Bi6DeepInteriorFree, slot: usize, amount: f64, star: Point3) {
    let cp: &mut Point3 = match slot {
        0 => &mut free.p00,
        1 | 2 => &mut free.tan[slot - 1],
        3..=5 => &mut free.two[slot - 3],
        6 => &mut free.p220,
        7 => &mut free.p210,
        s => {
            let k = s - 8;
            let n = free.p23.len();
            if k < n {
                &mut free.p23[k]
            } else {
                &mut free.p33[k - n]
            }
        }
    };
    let d = *cp - star;
    let len = d.length().max(1e-12);
    *cp = Point3::new(
        cp.x + amount * d.x / len,
        cp.y + amount * d.y / len,
        cp.z + amount * d.z / len,
    );
}

type ScoredGrids = (Vec<[[Point3; 7]; 7]>, f64, usize);

fn build_and_score(
    params: &Bi6Reparams,
    tad_src: &[SectorTensorBorderBi3],
    tak_src: &[SectorTensorBorderBi3],
    free: &Bi6DeepInteriorFree,
    builds: &mut usize,
) -> Result<ScoredGrids, String> {
    *builds += 1;
    let out = build_bi6_cap_from_tensor_border(params, tad_src, tak_src, free)?;
    let sal = (out.sal.x * out.sal.x + out.sal.y * out.sal.y + out.sal.z * out.sal.z).sqrt();
    let folds = out
        .sectors
        .iter()
        .filter(|g| {
            let standard: [[Point3; 7]; 7] =
                core::array::from_fn(|i| core::array::from_fn(|j| g[6 - i][6 - j]));
            let bi9 = elevate_bi6_grid_to_bi9(standard);
            GregoryPatch::RectangularBiNonic(Box::new(RectangularGregoryBiNonic { grid: bi9 }))
                .self_intersects_with_tolerance(SELF_INT_TOL)
        })
        .count();
    Ok((out.sectors, sal, folds))
}

// ---- eq (2)/(3) residuals, script convention ------------------------
// Port of the oracle's `sector_residuals` (same formula as the core
// `kp15_bi6_polar_cap_net_tests` harness); kernel-side so the rung can
// gate on it.

fn bern_eval(coeffs: &[f64], u: f64) -> f64 {
    let mut w = coeffs.to_vec();
    let m = w.len();
    for r in 1..m {
        for i in 0..(m - r) {
            w[i] = (1.0 - u) * w[i] + u * w[i + 1];
        }
    }
    w[0]
}

fn bern_diff(coeffs: &[f64]) -> Vec<f64> {
    let d = (coeffs.len() - 1) as f64;
    (0..coeffs.len() - 1)
        .map(|i| d * (coeffs[i + 1] - coeffs[i]))
        .collect()
}

#[allow(clippy::too_many_arguments)]
fn residuals_channel(
    grids: &[[[Point3; 7]; 7]],
    s: usize,
    n: usize,
    c: f64,
    tau: f64,
    u: f64,
    chan: usize,
) -> (f64, f64) {
    let t = (s + 1) % n;
    let comp = |p: Point3| -> f64 {
        match chan {
            0 => p.x,
            1 => p.y,
            _ => p.z,
        }
    };
    let cp_f = |i: usize, j: usize| comp(grids[s][6 - j][6 - i]);
    let cp_g = |i: usize, j: usize| comp(grids[t][6 - i][6 - j]);
    let f_row0: Vec<f64> = (0..7).map(|i| cp_f(i, 0)).collect();
    let fv: Vec<f64> = (0..7).map(|i| 6.0 * (cp_f(i, 1) - cp_f(i, 0))).collect();
    let fvv: Vec<f64> = (0..7)
        .map(|i| 30.0 * (cp_f(i, 2) - 2.0 * cp_f(i, 1) + cp_f(i, 0)))
        .collect();
    let gv: Vec<f64> = (0..7).map(|i| 6.0 * (cp_g(i, 1) - cp_g(i, 0))).collect();
    let gvv: Vec<f64> = (0..7)
        .map(|i| 30.0 * (cp_g(i, 2) - 2.0 * cp_g(i, 1) + cp_g(i, 0)))
        .collect();
    let val_fu = bern_eval(&bern_diff(&f_row0), u);
    let val_fuu = bern_eval(&bern_diff(&bern_diff(&f_row0)), u);
    let val_fv = bern_eval(&fv, u);
    let val_fuv = bern_eval(&bern_diff(&fv), u);
    let val_fvv = bern_eval(&fvv, u);
    let val_gv = bern_eval(&gv, u);
    let val_gvv = bern_eval(&gvv, u);
    let om = 1.0 - u;
    let g1 = val_fv + val_gv - 2.0 * c * om * om * val_fu;
    let a = 24.0 * c * (tau - 1.0) / (1.0 - u + tau * u);
    let g2 = -val_gvv + val_fvv - 4.0 * c * om * om * val_fuv
        + 4.0 * c * c * om.powi(4) * val_fuu
        + c * om.powi(3) * (a * u - 8.0 * c - a) * val_fu
        + a * om * om * val_fv;
    (g1, g2)
}

fn worst_eq2_eq3(grids: &[[[Point3; 7]; 7]], params: &Bi6Reparams) -> (f64, f64) {
    let n = params.n;
    let (mut w2, mut w3) = (0.0_f64, 0.0_f64);
    for s in 0..n {
        for chan in 0..3 {
            for step in 1..=9 {
                let u = step as f64 / 10.0;
                let (g1, g2) = residuals_channel(grids, s, n, params.c, params.tau, u, chan);
                w2 = w2.max(g1.abs());
                w3 = w3.max(g2.abs());
            }
        }
    }
    (w2, w3)
}

#[cfg(test)]
#[path = "kp15_setback_tbdata_rung_tests.rs"]
mod tests;
