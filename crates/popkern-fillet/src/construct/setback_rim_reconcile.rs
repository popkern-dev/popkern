//! Rim reconciliation for the setback cap fan: land the bi-9 cap
//! pieces' rim boundary rows on the live fillet sweeps.
//!
//! Nested inside `setback_vertex_topology` (shares its private
//! bezier9 helpers); extracted to hold the parent under its frozen
//! file-size ratchet.

use popkern_core::geometry::Point3;
use popkern_core::geometry::surface::SurfaceEval;

use super::bezier9_eval;

/// Project the rim-side boundary rows of the bi-9 cap pieces onto the
/// fillet sweeps (see call site). A grid side is a rim side when every
/// sampled point lies within `RIM_NEAR` of some fillet surface; the
/// projection deltas are applied to the boundary row and, in parallel,
/// to the two rows behind it. Deltas at shared corners are computed
/// from identical input points, so adjacent pieces stay sewn.
pub(super) fn reconcile_cap_rim_to_sweeps(
    patches: &mut [popkern_core::geometry::gregory::GregoryPatch],
    fillet_surfaces: &[crate::fillet_surface::FilletSurfaceResult],
    lateral_planes: &[popkern_core::geometry::SurfaceGeometry],
    border_loop: &[Vec<Point3>],
    rim_off_loop_allowance: f64,
) {
    use popkern_core::geometry::gregory::GregoryPatch;
    use popkern_core::geometry::uv_projection::project_point_to_uv;
    const RIM_NEAR: f64 = 2.0e-3;

    // The rim's neighbors are the sweeps along the arcs AND the
    // lateral planes along the corner segments (past the trim
    // crossings X the border leaves the fillet surface). Candidates
    // are enumerated so callers can track WHICH surface a sample
    // used: min-distance alone jumps between the two planes meeting
    // at a corner c, and interpolating the resulting discontinuous
    // field rings at the millimetre level.
    let candidates: Vec<popkern_core::geometry::SurfaceGeometry> = fillet_surfaces
        .iter()
        .map(|f| f.surface.clone())
        .chain(lateral_planes.iter().cloned())
        .collect();
    let project_on = |idx: usize, p: Point3| -> Option<(f64, Point3)> {
        let sfc = &candidates[idx];
        let (u, v) = project_point_to_uv(sfc, &p)?;
        let q = sfc.point_at(u, v);
        Some((q.distance_to(&p), q))
    };
    let project_best = |p: Point3| -> Option<(usize, f64, Point3)> {
        let mut best: Option<(usize, f64, Point3)> = None;
        for idx in 0..candidates.len() {
            if let Some((d, q)) = project_on(idx, p)
                && best.map(|(_, bd, _)| d < bd).unwrap_or(true)
            {
                best = Some((idx, d, q));
            }
        }
        best.filter(|(_, d, _)| *d < RIM_NEAR)
    };
    // Continuation projection: keep the previous sample's surface
    // while it stays near; switch only when it drifts (the genuine
    // arc→plane handoff at X happens where both surfaces agree).
    let project_seq = |pts: &[Point3]| -> Option<(Vec<Point3>, Vec<usize>)> {
        let mut out = Vec::with_capacity(pts.len());
        let mut used = Vec::with_capacity(pts.len());
        let mut current: Option<usize> = None;
        for &p in pts {
            let pick = match current {
                Some(idx) => match project_on(idx, p) {
                    Some((d, q)) if d < RIM_NEAR => Some((idx, d, q)),
                    _ => project_best(p),
                },
                None => project_best(p),
            };
            let (idx, _, q) = pick?;
            current = Some(idx);
            used.push(idx);
            out.push(q);
        }
        Some((out, used))
    };

    let probe = std::env::var_os("PK_SETBACK_PROBE").is_some();
    for (patch_index, patch) in patches.iter_mut().enumerate() {
        let GregoryPatch::RectangularBiNonic(bp) = patch else {
            continue;
        };
        // Four sides: (fixed index, along index, is_row, fixed at 0 or 9).
        // Represent as closures over (k) returning grid coordinates for
        // depth d in 0..3.
        for side in 0..4 {
            // Gather the boundary row for this side.
            let coord = |d: usize, k: usize| -> (usize, usize) {
                match side {
                    0 => (d, k),     // row i=0 boundary, depth into +i
                    1 => (9 - d, k), // row i=9 boundary, depth into -i
                    2 => (k, d),     // col j=0 boundary, depth into +j
                    _ => (k, 9 - d), // col j=9 boundary, depth into -j
                }
            };
            // Rim test on CURVE SAMPLES: control points legitimately
            // sit 1-2e-3 off their own curve, so a CP-based gate
            // wrongly rejects true rim sides.
            let side_row: [Point3; 10] = core::array::from_fn(|k| {
                let (i, j) = coord(0, k);
                bp.grid[i][j]
            });
            let gate_pts: Vec<Point3> = (0..=16)
                .map(|i| bezier9_eval(&side_row, i as f64 / 16.0))
                .collect();
            if project_seq(&gate_pts).is_none() {
                continue;
            }
            // Semantic gate: every sample of a rim side lies on the
            // setback border loop (see call site). Spoke seams near
            // the sweeps fail here at half the loop diameter.
            let worst_off_loop = gate_pts
                .iter()
                .map(|&p| distance_to_polylines(p, border_loop))
                .fold(0.0_f64, f64::max);
            if worst_off_loop > rim_off_loop_allowance.max(RIM_NEAR) {
                if probe {
                    eprintln!(
                        "  [rim-gate] patch {patch_index} side {side}: near sweeps but \
                         {worst_off_loop:.3e} off the border loop — not a rim side, skipped"
                    );
                }
                continue;
            }
            // Land the boundary CURVE on the sweep: projecting control
            // points directly leaves the curve ~3e-4 off (CPs are not
            // curve points). Fit the CP deltas to the sampled
            // deviation field by DAMPED least squares with pinned
            // endpoints, and iterate (the projection field is mildly
            // nonlinear). Exact 10-node degree-9 interpolation is
            // forbidden here: the equispaced Bernstein Vandermonde is
            // ill-conditioned and the deviation field kinks at the
            // arc→plane handoff, so measured 1.8e-3 node deviations
            // slung interior CPs by 2.8e-1 on hexagon6_variable —
            // tearing shared corner CPs 3.4e-1 apart and folding rim
            // normals (cap_outer g1 7.8e-1). The damping escalates
            // until the repair is at most DELTA_AMPLIFICATION_CAP ×
            // the deviation it repairs, keeping the move a bounded
            // local correction by construction.
            for _pass in 0..3 {
                let row: [Point3; 10] = core::array::from_fn(|k| {
                    let (i, j) = coord(0, k);
                    bp.grid[i][j]
                });
                const M: usize = 64;
                let ts: Vec<f64> = (0..M).map(|m| m as f64 / (M - 1) as f64).collect();
                let node_pts: Vec<Point3> = ts.iter().map(|&t| bezier9_eval(&row, t)).collect();
                let Some((projected, _)) = project_seq(&node_pts) else {
                    break;
                };
                let dev: Vec<[f64; 3]> = node_pts
                    .iter()
                    .zip(&projected)
                    .map(|(p, q)| [q.x - p.x, q.y - p.y, q.z - p.z])
                    .collect();
                let worst_dev = dev
                    .iter()
                    .map(|r| (r[0] * r[0] + r[1] * r[1] + r[2] * r[2]).sqrt())
                    .fold(0.0_f64, f64::max);
                // Endpoint CPs are curve points (B_0^9(0)=B_9^9(1)=1):
                // pin them to the exact endpoint deviations so shared
                // corners land on the sweep and adjacent solves agree.
                let d_end = [dev[0], dev[M - 1]];
                // DIAGNOSTIC (PK-208/209, probe-gated): the residual
                // FLOOR of this row — what an unconstrained fit of the
                // post-repair row against the projected field would
                // leave. Tells whether the seam-jet pinning (6 of 10
                // CPs) is what holds the rim above the tessellator's
                // 1e-3 off-surface gate, or the degree-9 row simply
                // cannot represent the field.
                if std::env::var_os("PK_SETBACK_PROBE").is_some() && _pass == 0 {
                    // Residual of the CURRENT row against the field is
                    // worst_dev. The floor: interpolate the projected
                    // points' deviation field with ALL 10 CPs free via
                    // damped LSQ (no pinning), then measure what the
                    // repaired row would miss.
                    let basis: Vec<[f64; 10]> = ts.iter().map(|&t| bernstein9_all(t)).collect();
                    let mut gram10 = [[0.0_f64; 10]; 10];
                    let mut mom10 = [[0.0_f64; 3]; 10];
                    for (b, r) in basis.iter().zip(&dev) {
                        for i in 0..10 {
                            for j in 0..10 {
                                gram10[i][j] += b[i] * b[j];
                            }
                            for c in 0..3 {
                                mom10[i][c] += b[i] * r[c];
                            }
                        }
                    }
                    for (i, row_g) in gram10.iter_mut().enumerate() {
                        row_g[i] += 1e-10;
                    }
                    if let Some(free) = solve_10(&gram10, &mom10) {
                        let mut floor = 0.0_f64;
                        for (b, r) in basis.iter().zip(&dev) {
                            let mut miss = 0.0_f64;
                            for c in 0..3 {
                                let mut fit = 0.0;
                                for i in 0..10 {
                                    fit += b[i] * free[i][c];
                                }
                                miss += (r[c] - fit) * (r[c] - fit);
                            }
                            floor = floor.max(miss.sqrt());
                        }
                        eprintln!(
                            "  [rim-floor] patch {patch_index} side {side}: field worst                              {worst_dev:.3e}, unconstrained-fit floor {floor:.3e}"
                        );
                    }
                }
                // The THREE CPs nearest each end (k=0..2, k=7..9) are
                // pinned to the SAME endpoint delta. These columns are
                // the cross-seam jets of the spoke/median seams meeting
                // the rim at the shared corners — k=1/k=8 the G1 jets,
                // k=2/k=7 the G2 rows: the paper cap's cross-seam
                // relations are affine (translation-invariant) per
                // transversal station, so translating boundary, jet,
                // AND second-jet CPs by one shared corner delta
                // preserves eq (8)/(9) and the G2 rows exactly no
                // matter what the neighbouring patch's own rim solve
                // does. Independently fitted k=1/k=8 deltas were the
                // measured author of the wide-fan seam creases (PK-094:
                // pentagon r=0.16 fan seam g1 1.478e-1); independently
                // fitted k=2/k=7 deltas were the measured author of the
                // fan-wide TRANSVERSE curvature jump (PK-187: g2t
                // 5.3e-1 normalized at the rim end of every spoke,
                // decaying to ~0 at the pole — the exact footprint of a
                // free G2 row over the three repaired depth rows, while
                // seam G1 stayed at 6e-5). The rim row keeps 4 free
                // interior CPs (k=3..6); the corner-region residual
                // this leaves against the sweeps is budgeted by the
                // cap-outer G0 gate (5e-3), two orders looser than the
                // seam gates this preserves.
                // Least squares for the 4 interior CP deltas against
                // the pinned-corrected deviation samples.
                let basis: Vec<[f64; 10]> = ts.iter().map(|&t| bernstein9_all(t)).collect();
                let mut gram = [[0.0_f64; 4]; 4];
                let mut moment = [[0.0_f64; 3]; 4];
                for (b, r) in basis.iter().zip(&dev) {
                    let rc: [f64; 3] = core::array::from_fn(|c| {
                        r[c] - (b[0] + b[1] + b[2]) * d_end[0][c]
                            - (b[7] + b[8] + b[9]) * d_end[1][c]
                    });
                    for i in 0..4 {
                        for j in 0..4 {
                            gram[i][j] += b[i + 3] * b[j + 3];
                        }
                        for c in 0..3 {
                            moment[i][c] += b[i + 3] * rc[c];
                        }
                    }
                }
                let mean_diag = (0..4).map(|i| gram[i][i]).sum::<f64>() / 4.0;
                const DELTA_AMPLIFICATION_CAP: f64 = 20.0;
                let cap = DELTA_AMPLIFICATION_CAP * worst_dev.max(1e-15);
                let mut interior: Option<[[f64; 3]; 4]> = None;
                let mut damping = 0.0_f64;
                loop {
                    let mut damped = gram;
                    for (i, row) in damped.iter_mut().enumerate() {
                        row[i] += damping * mean_diag;
                    }
                    if let Some(solution) = solve_4(&damped, &moment) {
                        let worst_delta = solution
                            .iter()
                            .map(|d| (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt())
                            .fold(0.0_f64, f64::max);
                        if worst_delta <= cap {
                            if probe && damping > 0.0 {
                                eprintln!(
                                    "  [rim-solve] patch {patch_index} side {side} pass {_pass}: \
                                     damping {damping:.1e} caps repair at {worst_delta:.3e} \
                                     (deviation {worst_dev:.3e})"
                                );
                            }
                            interior = Some(solution);
                            break;
                        }
                    }
                    damping = if damping == 0.0 {
                        1e-10
                    } else {
                        damping * 10.0
                    };
                    if damping > 1e3 {
                        break;
                    }
                }
                let Some(interior) = interior else {
                    break;
                };
                let mut deltas = [[0.0_f64; 3]; 10];
                deltas[0] = d_end[0];
                deltas[1] = d_end[0];
                deltas[2] = d_end[0];
                deltas[7] = d_end[1];
                deltas[8] = d_end[1];
                deltas[9] = d_end[1];
                deltas[3..7].copy_from_slice(&interior);
                for (k, delta) in deltas.iter().enumerate() {
                    for d in 0..3 {
                        let (i, j) = coord(d, k);
                        let g = bp.grid[i][j];
                        bp.grid[i][j] = Point3::new(g.x + delta[0], g.y + delta[1], g.z + delta[2]);
                    }
                }
            }
        }
    }
}

fn bernstein9_all(t: f64) -> [f64; 10] {
    const C: [f64; 10] = [1.0, 9.0, 36.0, 84.0, 126.0, 126.0, 84.0, 36.0, 9.0, 1.0];
    let mut out = [0.0_f64; 10];
    let u = 1.0 - t;
    let mut tp = 1.0;
    for (i, slot) in out.iter_mut().enumerate() {
        *slot = C[i] * tp * u.powi(9 - i as i32);
        tp *= t;
    }
    out
}

/// Solve the 6×6 linear system `mat · x = rhs` (3 right-hand sides)
/// by Gaussian elimination with partial pivoting.
// Column loops index two rows of the same matrix (pivot and target); the
// iterator form would need split borrows and obscure the elimination math.
#[allow(clippy::needless_range_loop)]
fn solve_10(mat: &[[f64; 10]; 10], rhs: &[[f64; 3]; 10]) -> Option<[[f64; 3]; 10]> {
    let mut a = *mat;
    let mut b = *rhs;
    for col in 0..10 {
        let mut piv = col;
        for r in (col + 1)..10 {
            if a[r][col].abs() > a[piv][col].abs() {
                piv = r;
            }
        }
        if a[piv][col].abs() < 1e-14 {
            return None;
        }
        a.swap(col, piv);
        b.swap(col, piv);
        let pivot = a[col];
        let pivot_b = b[col];
        for r in (col + 1)..10 {
            let f = a[r][col] / pivot[col];
            for (c, pv) in pivot.iter().enumerate().skip(col) {
                a[r][c] -= f * pv;
            }
            for (c, pv) in pivot_b.iter().enumerate() {
                b[r][c] -= f * pv;
            }
        }
    }
    let mut x = [[0.0_f64; 3]; 10];
    for row in (0..10).rev() {
        for c in 0..3 {
            let mut acc = b[row][c];
            for col in (row + 1)..10 {
                acc -= a[row][col] * x[col][c];
            }
            x[row][c] = acc / a[row][row];
        }
    }
    Some(x)
}

fn solve_4(mat: &[[f64; 4]; 4], rhs: &[[f64; 3]; 4]) -> Option<[[f64; 3]; 4]> {
    let mut a = *mat;
    let mut b = *rhs;
    for col in 0..4 {
        let mut piv = col;
        for r in (col + 1)..4 {
            if a[r][col].abs() > a[piv][col].abs() {
                piv = r;
            }
        }
        if a[piv][col].abs() < 1e-14 {
            return None;
        }
        a.swap(col, piv);
        b.swap(col, piv);
        let pivot = a[col];
        let pivot_b = b[col];
        for r in (col + 1)..4 {
            let f = a[r][col] / pivot[col];
            for (c, pv) in pivot.iter().enumerate().skip(col) {
                a[r][c] -= f * pv;
            }
            for (c, pv) in pivot_b.iter().enumerate() {
                b[r][c] -= f * pv;
            }
        }
    }
    let mut x = [[0.0_f64; 3]; 4];
    for row in (0..4).rev() {
        for c in 0..3 {
            let mut acc = b[row][c];
            for col in (row + 1)..4 {
                acc -= a[row][col] * x[col][c];
            }
            x[row][c] = acc / a[row][row];
        }
    }
    Some(x)
}

pub(super) fn eval_cubic(cps: &[Point3; 4], t: f64) -> Point3 {
    let u = 1.0 - t;
    let (b0, b1, b2, b3) = (u * u * u, 3.0 * u * u * t, 3.0 * u * t * t, t * t * t);
    Point3::new(
        b0 * cps[0].x + b1 * cps[1].x + b2 * cps[2].x + b3 * cps[3].x,
        b0 * cps[0].y + b1 * cps[1].y + b2 * cps[2].y + b3 * cps[3].y,
        b0 * cps[0].z + b1 * cps[1].z + b2 * cps[2].z + b3 * cps[3].z,
    )
}

/// Min distance from `p` to any segment of the sampled polylines.
fn distance_to_polylines(p: Point3, polylines: &[Vec<Point3>]) -> f64 {
    let mut best = f64::INFINITY;
    for polyline in polylines {
        for pair in polyline.windows(2) {
            let (a, b) = (pair[0], pair[1]);
            let ab = (b.x - a.x, b.y - a.y, b.z - a.z);
            let ap = (p.x - a.x, p.y - a.y, p.z - a.z);
            let len2 = ab.0 * ab.0 + ab.1 * ab.1 + ab.2 * ab.2;
            let s = if len2 > 0.0 {
                ((ap.0 * ab.0 + ap.1 * ab.1 + ap.2 * ab.2) / len2).clamp(0.0, 1.0)
            } else {
                0.0
            };
            let q = Point3::new(a.x + s * ab.0, a.y + s * ab.1, a.z + s * ab.2);
            best = best.min(p.distance_to(&q));
        }
    }
    best
}
