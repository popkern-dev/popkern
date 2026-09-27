//! Rim-side support measurement shared by the rim reconciles (split
//! from `setback_rim_curvature` for the 500-line cap): along a rim
//! side's rows, the cap's transverse normal curvature against the
//! support's, as the depth-2 correction that closes the gap; which
//! support a rim side lies on (whole-row attribution); and the
//! end-tangent alignment that puts each rim's end tangents in its
//! support's tangent plane (PK-576).

use popkern_core::geometry::gregory::GregoryPatch;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::uv_projection::project_point_to_uv;
use popkern_core::geometry::{Point3, Vec3};

use super::bezier9_eval;
use super::setback_corner_dual_jet::{RIM_NEAR, row_for_side};

/// Stations per rim side.
const M: usize = 33;

fn row_of(patches: &[GregoryPatch], patch: usize, side: usize, d: usize) -> Option<[Point3; 10]> {
    row_for_side(patches, patch, side, d)
}

pub(super) fn bezier9_derive(row: &[Point3; 10], t: f64) -> Vec3 {
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
pub(super) struct Station {
    pub(super) t: f64,
    pub(super) normal: Vec3,
    pub(super) delta_n_bvv: f64,
}

pub(super) fn side_stations(
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

/// A row within this distance of a strip along its whole length lies
/// ON that strip: the landed rims measure 1e-10 (constant radius) and
/// 3.5e-8 (variable), while the saddle bridge leaves every strip.
const ON_STRIP: f64 = 1.0e-6;

/// The support a rim side lies on, attributed over the WHOLE row. The
/// strips (the first `strip_count` candidates) come first: a row that
/// lies on a strip belongs to it, even where a plane touches the strip
/// along the same row (a spring line — the step's rims tie a support
/// plane to 1e-15 there, and targeting the plane's zero
/// cross-curvature declined its corner). Otherwise, the candidate with
/// the smallest worst distance over the row, within the rim band.
/// Per-point nearest-strip lookup mis-targets a rim side that lies on
/// a support PLANE (the saddle's reflex-face bridge): a strip is
/// within the band at every sample, so the side was tilted toward a
/// surface it does not lie on (PK-576).
pub(super) fn side_support(
    candidates: &[SurfaceGeometry],
    strip_count: usize,
    row0: &[Point3; 10],
) -> Option<usize> {
    let worst = |sfc: &SurfaceGeometry| -> f64 {
        (0..=16)
            .map(|i| {
                let p = bezier9_eval(row0, i as f64 / 16.0);
                project_point_to_uv(sfc, &p)
                    .map(|(u, v)| sfc.point_at(u, v).distance_to(&p))
                    .unwrap_or(f64::INFINITY)
            })
            .fold(0.0, f64::max)
    };
    let fits: Vec<(f64, usize)> = candidates
        .iter()
        .enumerate()
        .map(|(i, s)| (worst(s), i))
        .collect();
    let best = |on_strip: bool| {
        fits.iter()
            .filter(|(_, i)| !on_strip || *i < strip_count)
            .copied()
            .min_by(|a, b| a.0.total_cmp(&b.0))
    };
    if let Some((d, i)) = best(true)
        && d < ON_STRIP
    {
        return Some(i);
    }
    best(false).filter(|(d, _)| *d < RIM_NEAR).map(|(_, i)| i)
}

/// Rotate every rim side's END TANGENTS into its support's tangent
/// plane (PK-576). The G0 landing pins each end's three CPs to one
/// translation, so a rim's end tangent `row0[1] − row0[0]` stays the
/// border's. On a constant-radius strip that is tangent to 1e-7; on
/// a variable-radius strip it sat 1e-6..5e-5 off the plane. That is
/// the cap-strip G1 floor at every rim end, and at a setback corner
/// the neighbouring rim's cross-jet cell is exactly this CP, so no
/// later reconcile could reach it.
///
/// Seam safety: the move is the linear projection `c ↦ c − n(n·c)`
/// with ONE normal per landed endpoint. The fan pieces meeting at a
/// spoke junction share the endpoint and the support, so their
/// opposite cross vectors stay opposite (parametric C1) and collinear
/// vectors stay collinear (G1). Returns the worst tilt sine removed.
pub(super) fn align_rim_end_tangents(
    patches: &mut [GregoryPatch],
    candidates: &[SurfaceGeometry],
    strip_count: usize,
    sides: &[(usize, usize)],
) -> f64 {
    use super::setback_corner_dual_jet::coord_for;
    let mut worst = 0.0f64;
    for &(pi, side) in sides {
        let Some(row0) = row_of(patches, pi, side, 0) else {
            continue;
        };
        let Some(support) = side_support(candidates, strip_count, &row0) else {
            continue;
        };
        let sfc = &candidates[support];
        for (end, next) in [(0usize, 1usize), (9, 8)] {
            let Some(n) = support_normal_at(sfc, &row0[end]) else {
                continue;
            };
            let c = row0[next] - row0[end];
            let tilt = c.dot(&n);
            worst = worst.max(tilt.abs() / c.length().max(1e-300));
            let (i, j) = coord_for(side, 0, next);
            if let Some(grid) = patches[pi].binonic_grid_mut() {
                let g = grid[i][j];
                grid[i][j] = Point3::new(g.x - n.x * tilt, g.y - n.y * tilt, g.z - n.z * tilt);
            }
        }
    }
    worst
}

/// Unit normal of `sfc` at the projection of `p` (central
/// differences; `None` at a degenerate chart point).
fn support_normal_at(sfc: &SurfaceGeometry, p: &Point3) -> Option<Vec3> {
    let (u, v) = project_point_to_uv(sfc, p)?;
    let h = 1e-6;
    let du = sfc.point_at(u + h, v) - sfc.point_at(u - h, v);
    let dv = sfc.point_at(u, v + h) - sfc.point_at(u, v - h);
    let (lu, lv) = (du.length(), dv.length());
    if lu.min(lv) < 1e-3 * lu.max(lv) {
        return None;
    }
    du.cross(&dv).normalized()
}
