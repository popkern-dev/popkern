//! Implicit-distance rescue seeding for analytic × NURBS pairs.
//!
//! The patch-midpoint Newton seeder accepts converged samples
//! wherever they land, and on a NURBS whose extended basis
//! self-overlaps or re-crosses the analytic operand, EVERY midpoint
//! solve can fall into an extrapolation-phantom basin: a thread root
//! band's linear profile extension crosses a chamfer cone 0.26 off
//! the real band, all 1428 candidate solves converged there (or
//! failed), zero seeds survived, and the pair emitted nothing — the
//! merge then stranded every headed fastener (PK-460; whether it
//! happened flipped on a 5e-17 cone-axis tilt, i.e. basin luck).
//!
//! This seeder cannot escape a basin because it does not iterate: it
//! grids the NURBS at knot-span resolution INSIDE its query window,
//! evaluates the analytic operand's EXACT implicit signed distance at
//! each grid point, and mints a seed wherever the sign changes
//! between neighbours — a sign change on a real surface sample chain
//! is a crossing certificate, not a guess. Newton then only POLISHES
//! from an on-surface point next to the true curve.
//!
//! Scoped as a RESCUE: the caller invokes it only when the standard
//! seeder produced no usable transversal seed, so tuned pairs keep
//! their seed sets byte-identical.

use crate::geometry::Point3;
use crate::geometry::analytic_surfaces::analytic_signed_distance;
use crate::geometry::intersection::marching::{MarchSample, march_project};
use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};
use crate::geometry::uv_projection::project_point_to_uv_multisheet;
use crate::tolerance::Tolerance;

/// Cap on grid cells per axis (a 198-span thread wall grids at 96).
const MAX_CELLS: usize = 96;
/// Cap on rescued seeds (a long intersection needs only a few — the
/// march traces the rest).
const MAX_SEEDS: usize = 64;

/// Push rescue seeds for an analytic × NURBS pair. No-op when neither
/// orientation matches or no sign change exists in the window.
#[allow(clippy::too_many_arguments)]
pub(super) fn push_implicit_grid_seeds(
    seeds: &mut Vec<(MarchSample, bool)>,
    surface1: &SurfaceGeometry,
    u1_range: (f64, f64),
    v1_range: (f64, f64),
    surface2: &SurfaceGeometry,
    u2_range: (f64, f64),
    v2_range: (f64, f64),
    tolerance: &Tolerance,
) {
    let s1_analytic = analytic_signed_distance(surface1, &surface1.point_at(0.0, 0.0)).is_some();
    let s2_analytic = analytic_signed_distance(surface2, &surface2.point_at(0.0, 0.0)).is_some();
    let s1_nurbs = matches!(surface1, SurfaceGeometry::Nurbs(_));
    let s2_nurbs = matches!(surface2, SurfaceGeometry::Nurbs(_));
    let (analytic_is_1, analytic, nurbs, n_ur, n_vr) = if s1_analytic && s2_nurbs {
        (true, surface1, surface2, u2_range, v2_range)
    } else if s2_analytic && s1_nurbs {
        (false, surface2, surface1, u1_range, v1_range)
    } else {
        return;
    };

    // Knot-span-resolving grid over the NURBS window.
    let spans = |knots: &[f64], lo: f64, hi: f64| -> usize {
        knots
            .windows(2)
            .filter(|w| w[1] > w[0] + 1e-12 && w[1] > lo && w[0] < hi)
            .count()
            .max(1)
    };
    let SurfaceGeometry::Nurbs(ns) = nurbs else {
        return;
    };
    let nu = (2 * spans(&ns.knots_u, n_ur.0, n_ur.1)).clamp(8, MAX_CELLS);
    let nv = (2 * spans(&ns.knots_v, n_vr.0, n_vr.1)).clamp(8, MAX_CELLS);

    // Cheap prefilter: most seedless pairs are AABB-overlapping but
    // GENUINELY disjoint, and the full grid costs ~9k evaluations. A
    // 9×9 scan bounds how close the NURBS comes to the analytic
    // surface; when even the closest coarse sample is further than a
    // coarse cell can hide, no sign change exists and the rescue
    // skips (the full grid runs only for genuinely near pairs).
    {
        let mut min_abs = f64::MAX;
        let mut coarse_pts: Vec<Point3> = Vec::with_capacity(81);
        for i in 0..=8 {
            for j in 0..=8 {
                let u = n_ur.0 + (n_ur.1 - n_ur.0) * i as f64 / 8.0;
                let v = n_vr.0 + (n_vr.1 - n_vr.0) * j as f64 / 8.0;
                let pnt = nurbs.point_at(u, v);
                if let Some(d) = analytic_signed_distance(analytic, &pnt) {
                    min_abs = min_abs.min(d.abs());
                }
                coarse_pts.push(pnt);
            }
        }
        // Longest coarse-cell diagonal bounds the surface excursion a
        // coarse sample can miss.
        let mut max_cell = 0.0_f64;
        for i in 0..8 {
            for j in 0..8 {
                let a = coarse_pts[i * 9 + j];
                let b = coarse_pts[(i + 1) * 9 + j + 1];
                max_cell = max_cell.max(a.distance_to(&b));
            }
        }
        if min_abs > max_cell {
            return;
        }
    }

    let d_at = |i: usize, j: usize| -> (f64, f64, Option<f64>) {
        let u = n_ur.0 + (n_ur.1 - n_ur.0) * i as f64 / nu as f64;
        let v = n_vr.0 + (n_vr.1 - n_vr.0) * j as f64 / nv as f64;
        let p = nurbs.point_at(u, v);
        (u, v, analytic_signed_distance(analytic, &p))
    };
    let mut grid = vec![(0.0_f64, 0.0_f64, None::<f64>); (nu + 1) * (nv + 1)];
    for i in 0..=nu {
        for j in 0..=nv {
            grid[i * (nv + 1) + j] = d_at(i, j);
        }
    }

    let mut pushed: Vec<Point3> = Vec::new();
    let mut try_seed = |u_lo: (f64, f64, f64), u_hi: (f64, f64, f64)| {
        if pushed.len() >= MAX_SEEDS {
            return;
        }
        let (ua, va, da) = u_lo;
        let (ub, vb, db) = u_hi;
        // Linear zero of the implicit distance along the grid edge.
        let t = da / (da - db);
        let (gu, gv) = (ua + (ub - ua) * t, va + (vb - va) * t);
        let p = nurbs.point_at(gu, gv);
        let Some((au, av)) = project_point_to_uv_multisheet(analytic, &p) else {
            return;
        };
        // Polish with the standard equality-constrained Newton,
        // anchored AT the certificate point (no foreign basin).
        let n1;
        let n2;
        let (u1, v1, u2, v2) = if analytic_is_1 {
            n1 = analytic.normal_at(au, av);
            n2 = nurbs.normal_at(gu, gv);
            (au, av, gu, gv)
        } else {
            n1 = nurbs.normal_at(gu, gv);
            n2 = analytic.normal_at(au, av);
            (gu, gv, au, av)
        };
        let dir = {
            let c = n1.cross(&n2);
            match c.normalized() {
                Some(d) => d,
                None => return, // tangent contact — not this seeder's class
            }
        };
        let Some(sample) = march_project(
            if analytic_is_1 { analytic } else { nurbs },
            u1,
            v1,
            if analytic_is_1 { nurbs } else { analytic },
            u2,
            v2,
            &p,
            &dir,
            tolerance,
        ) else {
            return;
        };
        // Only accept polishes that stayed near the certificate: a
        // polish that flew off is the phantom class again.
        if sample.point.distance_to(&p) > 100.0 * tolerance.linear() {
            return;
        }
        let dedup = tolerance.linear() * 100.0;
        if pushed.iter().any(|q| q.distance_to(&sample.point) < dedup) {
            return;
        }
        pushed.push(sample.point);
        seeds.push((sample, false));
    };

    for i in 0..=nu {
        for j in 0..nv {
            let a = grid[i * (nv + 1) + j];
            let b = grid[i * (nv + 1) + j + 1];
            if let (Some(da), Some(db)) = (a.2, b.2)
                && da.signum() != db.signum()
                && da != 0.0
            {
                try_seed((a.0, a.1, da), (b.0, b.1, db));
            }
        }
    }
    for j in 0..=nv {
        for i in 0..nu {
            let a = grid[i * (nv + 1) + j];
            let b = grid[(i + 1) * (nv + 1) + j];
            if let (Some(da), Some(db)) = (a.2, b.2)
                && da.signum() != db.signum()
                && da != 0.0
            {
                try_seed((a.0, a.1, da), (b.0, b.1, db));
            }
        }
    }
    if !pushed.is_empty() && crate::trace::trace_enabled("seed") {
        eprintln!(
            "  [seed] implicit-grid rescue: {} seed(s) from {}x{} grid",
            pushed.len(),
            nu,
            nv
        );
    }
}
