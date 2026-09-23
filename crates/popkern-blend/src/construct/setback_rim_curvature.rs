//! Rim transverse-curvature reconciliation (PK-515): after the G0
//! rim landing and the G1 cross-jet reconcile, move each rim side's
//! DEPTH-2 row so the cap's transverse normal curvature along the
//! rim follows the strip's, measured pointwise on both sides.
//!
//! Why this exists: the border derivation encodes the strips' cross
//! curvature into the bi-6 fan through ONE shared Hermite block per
//! corner, and two neighbouring fillet strips are different surfaces
//! whose cross fields agree only to second order at the corner — the
//! block cannot carry both (measured floor 3.4e-4 in row units on the
//! G2 pentagon flagship at ANY objective weight, ≈ 1 in transverse
//! curvature against the 5e-2 commercial gate; re-weighting only
//! relocated the error into jet tilt and corner twist). The committed
//! pieces are bi-9: their ten-CP depth-2 rows have the freedom the
//! bi-6 block lacks, and a depth-2-only move leaves rows 0 and 1 —
//! G0, tangency, and the rim twist — bit-identical.
//!
//! Cross-seam discipline (as the jet reconcile): the three CPs
//! nearest each end take ONE shared per-corner delta (translation-
//! invariant G rows); the interior four CPs fit the measured field.

use popkern_core::geometry::gregory::GregoryPatch;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::uv_projection::project_point_to_uv;
use popkern_core::geometry::{Point3, Vec3};

use super::bezier9_eval;
use super::setback_rim_reconcile::distance_to_polylines;

const M: usize = 33;
const RIM_NEAR: f64 = 2.0e-3;
const CORNER_CLUSTER: f64 = 1.0e-3;

fn coord_for(side: usize, d: usize, k: usize) -> (usize, usize) {
    match side {
        0 => (d, k),
        1 => (9 - d, k),
        2 => (k, d),
        _ => (k, 9 - d),
    }
}

fn row_of(patches: &[GregoryPatch], patch: usize, side: usize, d: usize) -> Option<[Point3; 10]> {
    let GregoryPatch::RectangularBiNonic(bp) = &patches[patch] else {
        return None;
    };
    Some(core::array::from_fn(|k| {
        let (i, j) = coord_for(side, d, k);
        bp.grid[i][j]
    }))
}

fn bezier9_derive(row: &[Point3; 10], t: f64) -> Vec3 {
    let d: [Point3; 9] = core::array::from_fn(|k| {
        let v = row[k + 1] - row[k];
        Point3::new(9.0 * v.x, 9.0 * v.y, 9.0 * v.z)
    });
    let mut p = d;
    for level in (1..9).rev() {
        for i in 0..level {
            p[i] = Point3::new(
                p[i].x + (p[i + 1].x - p[i].x) * t,
                p[i].y + (p[i + 1].y - p[i].y) * t,
                p[i].z + (p[i + 1].z - p[i].z) * t,
            );
        }
    }
    Vec3::new(p[0].x, p[0].y, p[0].z)
}

/// Strip-side normal + transverse normal curvature at `p` in
/// direction `t_perp`, from the nearest strip within the rim band.
fn strip_normal_curvature(
    strips: &[SurfaceGeometry],
    p: &Point3,
    t_perp: Vec3,
) -> Option<(Vec3, f64)> {
    let mut best: Option<(f64, usize, (f64, f64))> = None;
    for (idx, sfc) in strips.iter().enumerate() {
        let Some((u, v)) = project_point_to_uv(sfc, p) else {
            continue;
        };
        let d = sfc.point_at(u, v).distance_to(p);
        if d < RIM_NEAR && best.map(|(bd, _, _)| d < bd).unwrap_or(true) {
            best = Some((d, idx, (u, v)));
        }
    }
    let (_, idx, (u, v)) = best?;
    let s = &strips[idx];
    let h = 1e-6;
    let su = (s.point_at(u + h, v) - s.point_at(u - h, v)) * (0.5 / h);
    let sv = (s.point_at(u, v + h) - s.point_at(u, v - h)) * (0.5 / h);
    let n = su.cross(&sv).normalized()?;
    let (e, f_, g) = (su.dot(&su), su.dot(&sv), sv.dot(&sv));
    let (r1, r2) = (t_perp.dot(&su), t_perp.dot(&sv));
    let det = e * g - f_ * f_;
    if det.abs() < 1e-30 {
        return None;
    }
    let (a, b) = ((r1 * g - r2 * f_) / det, (r2 * e - r1 * f_) / det);
    let estimate = |step: f64| -> Option<f64> {
        let p0 = s.point_at(u, v);
        let p1 = s.point_at(u + a * step, v + b * step);
        let p2 = s.point_at(u - a * step, v - b * step);
        let second = (p1 - p0) + (p2 - p0);
        let first = (p1 - p2).length() * 0.5;
        (first >= 1e-30).then(|| second.dot(&n) / (first * first))
    };
    let k1 = estimate(1e-4)?;
    let k2 = estimate(5e-5)?;
    if (k1 - k2).abs() > 0.1 * (1.0 + k2.abs()) {
        return None;
    }
    Some((n, (4.0 * k2 - k1) / 3.0))
}

/// One station's measured correction: `delta_n_bvv` is the change of
/// `n·B_vv` (bi-9 cross second derivative) that closes the cap↔strip
/// transverse-curvature gap at this station.
struct Station {
    t: f64,
    normal: Vec3,
    delta_n_bvv: f64,
}

fn side_stations(
    patches: &[GregoryPatch],
    strips: &[SurfaceGeometry],
    pi: usize,
    side: usize,
) -> Option<Vec<Station>> {
    let row0 = row_of(patches, pi, side, 0)?;
    let row1 = row_of(patches, pi, side, 1)?;
    let row2 = row_of(patches, pi, side, 2)?;
    let mut out = Vec::with_capacity(M);
    for m in 0..M {
        let t = m as f64 / (M - 1) as f64;
        let p0 = bezier9_eval(&row0, t);
        let b_u = bezier9_derive(&row0, t);
        let p1 = bezier9_eval(&row1, t);
        let p2 = bezier9_eval(&row2, t);
        let b_v = (p1 - p0) * 9.0;
        let b_vv = ((p2 - p1) - (p1 - p0)) * 72.0;
        let b_uv = (bezier9_derive(&row1, t) - b_u) * 9.0;
        let b_uu = {
            let h = 1e-4;
            let tc = t.clamp(h, 1.0 - h);
            (bezier9_derive(&row0, tc + h) - bezier9_derive(&row0, tc - h)) * (0.5 / h)
        };
        // Shared frame: rim tangent, strip normal, and its transverse
        // (the first curvature call only fetches the normal; the
        // second measures in the true t_perp).
        let tang = b_u.normalized()?;
        let (n_strip, _) = strip_normal_curvature(strips, &p0, tang)?;
        let t_perp = n_strip.cross(&tang).normalized()?;
        let (_, kappa_strip_perp) = strip_normal_curvature(strips, &p0, t_perp)?;
        // Cap curvature in t_perp from the patch jets: t_perp = a·B_u
        // + b·B_v in the tangent plane (least squares), then
        // II(t_perp) = a²·n·B_uu + 2ab·n·B_uv + b²·n·B_vv over
        // |t_perp|² = 1.
        let (e, f_, g) = (b_u.dot(&b_u), b_u.dot(&b_v), b_v.dot(&b_v));
        let det = e * g - f_ * f_;
        if det.abs() < 1e-30 {
            return None;
        }
        let (r1, r2) = (t_perp.dot(&b_u), t_perp.dot(&b_v));
        let (a, b) = ((r1 * g - r2 * f_) / det, (r2 * e - r1 * f_) / det);
        if b.abs() < 1e-9 {
            return None;
        }
        let kappa_cap = a * a * n_strip.dot(&b_uu)
            + 2.0 * a * b * n_strip.dot(&b_uv)
            + b * b * n_strip.dot(&b_vv);
        let delta_kappa = kappa_strip_perp - kappa_cap;
        out.push(Station {
            t,
            normal: n_strip,
            delta_n_bvv: delta_kappa / (b * b),
        });
    }
    Some(out)
}

/// Register rim sides exactly like the jet reconcile: near a support
/// along the whole row AND on the setback border loop.
fn rim_sides(
    patches: &[GregoryPatch],
    candidates: &[SurfaceGeometry],
    border_loop: &[Vec<Point3>],
    rim_off_loop_allowance: f64,
) -> Vec<(usize, usize)> {
    let mut sides: Vec<(usize, usize)> = Vec::new();
    for (pi, p) in patches.iter().enumerate() {
        let GregoryPatch::RectangularBiNonic(bp) = p else {
            continue;
        };
        for side in 0..4 {
            let row: [Point3; 10] = core::array::from_fn(|k| {
                let (i, j) = coord_for(side, 0, k);
                bp.grid[i][j]
            });
            let on_support = (0..=16).all(|i| {
                let pt = bezier9_eval(&row, i as f64 / 16.0);
                candidates.iter().any(|sfc| {
                    project_point_to_uv(sfc, &pt)
                        .map(|(u, v)| sfc.point_at(u, v).distance_to(&pt) < RIM_NEAR)
                        .unwrap_or(false)
                })
            });
            if !on_support {
                continue;
            }
            let worst_off_loop = (0..=16)
                .map(|i| distance_to_polylines(bezier9_eval(&row, i as f64 / 16.0), border_loop))
                .fold(0.0_f64, f64::max);
            if worst_off_loop > rim_off_loop_allowance.max(RIM_NEAR) {
                continue;
            }
            sides.push((pi, side));
        }
    }
    sides
}

/// Run the rim transverse-curvature reconciliation over the fan.
pub(super) fn reconcile_rim_curvature(
    patches: &mut [GregoryPatch],
    candidates: &[SurfaceGeometry],
    strip_count: usize,
    border_loop: &[Vec<Point3>],
    rim_off_loop_allowance: f64,
) {
    let strips = &candidates[..strip_count.min(candidates.len())];
    let sides = rim_sides(patches, candidates, border_loop, rim_off_loop_allowance);
    if sides.is_empty() {
        return;
    }
    // Per-side station fields; a side whose measurement declines at
    // any station is skipped whole (a partial correction would kink
    // the row against its clamped ends).
    let mut fields: Vec<((usize, usize), Vec<Station>)> = Vec::new();
    for &(pi, side) in &sides {
        if let Some(stations) = side_stations(patches, strips, pi, side) {
            fields.push(((pi, side), stations));
        }
    }
    // ONE joint least-squares over the whole ring: unknowns are the
    // shared per-corner deltas (one Vec3 per corner cluster - the
    // three CPs nearest each end of every incident side take it) and
    // each side's four interior CPs. Solving jointly picks the corner
    // constant that best serves the clamp-region FIELD of both
    // incident sides - an endpoint-pinned constant left the measured
    // residual peaked at 0.68 exactly inside the clamp region.
    struct SideLayout {
        pi: usize,
        side: usize,
        start_cluster: usize,
        end_cluster: usize,
        interior_base: usize,
    }
    let mut centers: Vec<Point3> = Vec::new();
    let mut cluster_of = |p: Point3| -> usize {
        centers
            .iter()
            .position(|c| c.distance_to(&p) < CORNER_CLUSTER)
            .unwrap_or_else(|| {
                centers.push(p);
                centers.len() - 1
            })
    };
    let mut layouts: Vec<SideLayout> = Vec::new();
    for ((pi, side), _) in &fields {
        let Some(row0) = row_of(patches, *pi, *side, 0) else {
            continue;
        };
        layouts.push(SideLayout {
            pi: *pi,
            side: *side,
            start_cluster: cluster_of(bezier9_eval(&row0, 0.0)),
            end_cluster: cluster_of(bezier9_eval(&row0, 1.0)),
            interior_base: 0,
        });
    }
    let n_clusters = centers.len();
    for (i, layout) in layouts.iter_mut().enumerate() {
        layout.interior_base = n_clusters + 4 * i;
    }
    let unknowns = n_clusters + 4 * layouts.len();
    let bern9 = |k: usize, t: f64| -> f64 {
        const C: [f64; 10] = [1.0, 9.0, 36.0, 84.0, 126.0, 126.0, 84.0, 36.0, 9.0, 1.0];
        C[k] * t.powi(k as i32) * (1.0 - t).powi((9 - k) as i32)
    };
    // Corner-compatibility gate: the shared clamp constant serves
    // BOTH incident sides, so it is safe only where their clamp-band
    // fields agree (a uniform-radius apex, by symmetry). On mixed
    // radii the sides genuinely need different corrections; the
    // averaged constant pushed both off their strips (measured:
    // pyramid5_mixed cap_outer g0 6.7e-6 against its 1e-9 bar). A
    // disagreeing cluster's delta is frozen at zero — that corner
    // keeps its pre-pass rows and only the interiors correct.
    let mut cluster_sides: Vec<Vec<Vec3>> = vec![Vec::new(); n_clusters];
    for (layout, ((_, _), stations)) in layouts.iter().zip(&fields) {
        for (cluster, end) in [(layout.start_cluster, false), (layout.end_cluster, true)] {
            let mut acc = Vec3::zero();
            let mut mass = 0.0_f64;
            for station in stations {
                let clamp_mass: f64 = if end {
                    (7..=9).map(|k| bern9(k, station.t)).sum()
                } else {
                    (0..=2).map(|k| bern9(k, station.t)).sum()
                };
                acc = acc + station.normal * (station.delta_n_bvv / 72.0 * clamp_mass);
                mass += clamp_mass;
            }
            if mass > 0.0 {
                cluster_sides[cluster].push(acc * (1.0 / mass));
            }
        }
    }
    let cluster_frozen: Vec<bool> = cluster_sides
        .iter()
        .map(|sides| {
            let magnitude = sides.iter().map(|v| v.length()).fold(0.0_f64, f64::max);
            let disagreement = sides
                .iter()
                .flat_map(|a| sides.iter().map(move |b| (*a - *b).length()))
                .fold(0.0_f64, f64::max);
            sides.len() < 2 || disagreement > 0.5 * magnitude.max(f64::MIN_POSITIVE)
        })
        .collect();
    let mut ata = vec![vec![0.0_f64; unknowns]; unknowns];
    let mut atb = vec![[0.0_f64; 3]; unknowns];
    for (layout, ((_, _), stations)) in layouts.iter().zip(&fields) {
        for station in stations {
            let target = station.normal * (station.delta_n_bvv / 72.0);
            // Row: clamp mass on the two corner unknowns + interior.
            let mut entries: Vec<(usize, f64)> = Vec::with_capacity(6);
            let start_mass: f64 = (0..=2).map(|k| bern9(k, station.t)).sum();
            let end_mass: f64 = (7..=9).map(|k| bern9(k, station.t)).sum();
            if !cluster_frozen[layout.start_cluster] {
                entries.push((layout.start_cluster, start_mass));
            }
            if !cluster_frozen[layout.end_cluster] {
                entries.push((layout.end_cluster, end_mass));
            }
            for k in 3..=6_usize {
                entries.push((layout.interior_base + (k - 3), bern9(k, station.t)));
            }
            for &(i, wi) in &entries {
                for &(j, wj) in &entries {
                    ata[i][j] += wi * wj;
                }
                atb[i][0] += wi * target.x;
                atb[i][1] += wi * target.y;
                atb[i][2] += wi * target.z;
            }
        }
    }
    let Some(solution) = solve_dense_spd(&mut ata, &atb) else {
        return;
    };
    // Protected cells: any registered rim side's d = 0 (G0) and
    // d = 1 (tangency/twist) rows are untouchable. The CORNER piece
    // carries TWO rim sides meeting at `c`, and each side's d = 2 row
    // crosses the OTHER side's shallow rows there — writing those
    // cells moved a genuine rim boundary row and opened the
    // cap↔strip seam (measured: pyramid5_mixed cap_outer g0 6.1e-6
    // against its 1e-9 bar). Cells claimed by two sides' d = 2 rows
    // (the (7,7)-class corner cell) take the MEAN of the claims.
    let mut protected: std::collections::BTreeSet<(usize, usize, usize)> =
        std::collections::BTreeSet::new();
    for layout in &layouts {
        for d in 0..=1_usize {
            for k in 0..10_usize {
                let (i, j) = coord_for(layout.side, d, k);
                protected.insert((layout.pi, i, j));
            }
        }
    }
    let mut writes: std::collections::BTreeMap<(usize, usize, usize), (Vec3, usize)> =
        std::collections::BTreeMap::new();
    for layout in &layouts {
        for k in 0..10_usize {
            let delta = match k {
                0..=2 if cluster_frozen[layout.start_cluster] => Vec3::zero(),
                0..=2 => solution[layout.start_cluster],
                7..=9 if cluster_frozen[layout.end_cluster] => Vec3::zero(),
                7..=9 => solution[layout.end_cluster],
                _ => solution[layout.interior_base + (k - 3)],
            };
            let (i2, j2) = coord_for(layout.side, 2, k);
            if protected.contains(&(layout.pi, i2, j2)) {
                continue;
            }
            let entry = writes
                .entry((layout.pi, i2, j2))
                .or_insert((Vec3::zero(), 0));
            entry.0 = entry.0 + delta;
            entry.1 += 1;
        }
    }
    for ((pi, i2, j2), (sum, count)) in writes {
        let GregoryPatch::RectangularBiNonic(bp) = &mut patches[pi] else {
            continue;
        };
        let delta = sum * (1.0 / count.max(1) as f64);
        let g = bp.grid[i2][j2];
        bp.grid[i2][j2] = Point3::new(g.x + delta.x, g.y + delta.y, g.z + delta.z);
    }
    if popkern_core::trace::trace_enabled("setback") {
        for layout in &layouts {
            let d_max = (0..10_usize)
                .map(|k| match k {
                    0..=2 if cluster_frozen[layout.start_cluster] => 0.0,
                    0..=2 => solution[layout.start_cluster].length(),
                    7..=9 if cluster_frozen[layout.end_cluster] => 0.0,
                    7..=9 => solution[layout.end_cluster].length(),
                    _ => solution[layout.interior_base + (k - 3)].length(),
                })
                .fold(0.0_f64, f64::max);
            eprintln!(
                "  [rim-curvature] applied patch {} side {} |delta|max {:.3e} frozen=({},{})",
                layout.pi,
                layout.side,
                d_max,
                cluster_frozen[layout.start_cluster],
                cluster_frozen[layout.end_cluster]
            );
        }
        let worst = fields
            .iter()
            .flat_map(|(_, s)| s.iter())
            .map(|s| (s.delta_n_bvv / 72.0).abs())
            .fold(0.0_f64, f64::max);
        eprintln!(
            "  [rim-curvature] {} sides, {} corner clusters: worst |drow2| target {worst:.3e} applied (joint solve)",
            layouts.len(),
            n_clusters
        );
    }
}

/// Dense symmetric positive-definite solve (Cholesky, in place),
/// three right-hand sides; `None` on a non-SPD system.
fn solve_dense_spd(ata: &mut [Vec<f64>], atb: &[[f64; 3]]) -> Option<Vec<Vec3>> {
    let n = ata.len();
    for (i, row) in ata.iter_mut().enumerate() {
        row[i] += 1e-12 * (1.0 + row[i].abs());
    }
    for j in 0..n {
        let head: f64 = ata[j].iter().take(j).map(|v| v * v).sum();
        let d = ata[j][j] - head;
        if d <= 0.0 || !d.is_finite() {
            return None;
        }
        let d = d.sqrt();
        ata[j][j] = d;
        for i in (j + 1)..n {
            let cross: f64 = ata[i]
                .iter()
                .zip(ata[j].iter())
                .take(j)
                .map(|(a, b)| a * b)
                .sum();
            ata[i][j] = (ata[i][j] - cross) / d;
        }
    }
    // All three axes at once with Vec3 arithmetic.
    let b: Vec<Vec3> = atb.iter().map(|r| Vec3::new(r[0], r[1], r[2])).collect();
    let mut y = vec![Vec3::zero(); n];
    for i in 0..n {
        let mut s = b[i];
        for (k, yk) in y.iter().enumerate().take(i) {
            s = s - *yk * ata[i][k];
        }
        y[i] = s * (1.0 / ata[i][i]);
    }
    let mut x = vec![Vec3::zero(); n];
    for i in (0..n).rev() {
        let mut s = y[i];
        for (k, xk) in x.iter().enumerate().skip(i + 1) {
            s = s - *xk * ata[k][i];
        }
        x[i] = s * (1.0 / ata[i][i]);
    }
    Some(x)
}
