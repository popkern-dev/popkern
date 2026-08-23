//! Shared rim-side machinery for G1 hole caps (PK-224, PK-248).
//!
//! Both the single-patch cap and the polar-wedge fan meet the
//! surviving walls the same way: each rim side becomes a quintic
//! Bezier with ZERO end legs (`[S, S, A, B, E, E]` — the corner
//! parameterization cusps so the adjacent side's corner leg cannot
//! pollute this side's cross field), and the cap's cross-derivative
//! along the side is carried by ribbon directions `n × t` oriented
//! in-plane away from the neighbour face. This module holds that
//! machinery once so the two rungs cannot drift apart.

use crate::geometry::surface::SurfaceGeometry;
use crate::geometry::{FillBoundaryEdge, Point3, Vec3};
use crate::topology::BRep;

use super::hole_fill::DetectedHole;

/// Why a rim side could not supply its quintic row or ribbon
/// directions. Callers map these into their own decline enums.
#[derive(Debug, thiserror::Error)]
pub(super) enum RimRibbonError {
    /// The side's quintic fit missed the sampled rim curve.
    #[error("side {side} quintic fit residual {residual:.3e} exceeds {budget:.3e}")]
    FitResidual {
        /// Side index in ring order.
        side: usize,
        /// Worst sample distance to the fitted Bezier.
        residual: f64,
        /// Linear tolerance budget scaled for fitting.
        budget: f64,
    },
    /// The side's ribbon direction could not be oriented.
    #[error("side {side} has no usable ribbon direction")]
    RibbonDegenerate {
        /// Side index in ring order.
        side: usize,
    },
}

/// Quintic Bezier `[S, S, A, B, E, E]` for a rim side. Zero end legs
/// by construction (see module doc); `A`, `B` interpolate the rim at
/// the 1/3 and 2/3 stations, and a GEOMETRIC residual gate (sampled
/// point-to-curve distance — the cusped parameterization cannot and
/// need not match the rim's own parameter) verifies containment.
/// `orient_loop` flips a side's endpoints without flipping its dense
/// samples, so the station lookup re-orients them by endpoint
/// distance first.
pub(super) fn fit_quintic_zero_leg(
    frame: &FillBoundaryEdge,
    side: usize,
    tolerance: f64,
) -> Result<[Point3; 6], RimRibbonError> {
    let bern = |i: usize, t: f64| -> f64 {
        let s = 1.0 - t;
        let c = [1.0, 5.0, 10.0, 10.0, 5.0, 1.0][i];
        c * t.powi(i as i32) * s.powi(5 - i as i32)
    };

    // Dense rim samples in THIS side's orientation.
    let mut samples: Vec<(f64, Point3)> = frame.extra_samples.clone();
    let reversed = samples
        .first()
        .map(|&(_, p)| p.distance_to(&frame.end) < p.distance_to(&frame.start))
        .unwrap_or(false);
    if reversed {
        for (f, _) in &mut samples {
            *f = 1.0 - *f;
        }
        samples.reverse();
    }
    samples.push((0.5, frame.mid));

    // Station targets at 1/3 and 2/3 (nearest dense samples).
    let station_point = |f0: f64| -> Point3 {
        samples
            .iter()
            .min_by(|a, b| (a.0 - f0).abs().total_cmp(&(b.0 - f0).abs()))
            .map(|&(_, p)| p)
            .unwrap_or(frame.mid)
    };
    let (t1, t2) = (1.0 / 3.0, 2.0 / 3.0);
    let p1 = station_point(t1);
    let p2 = station_point(t2);

    // Solve A, B so the curve interpolates p1 at t1 and p2 at t2.
    let a11 = bern(2, t1);
    let a12 = bern(3, t1);
    let a21 = bern(2, t2);
    let a22 = bern(3, t2);
    let det = a11 * a22 - a12 * a21;
    if det.abs() < 1e-18 {
        return Err(RimRibbonError::FitResidual {
            side,
            residual: f64::INFINITY,
            budget: tolerance,
        });
    }
    let fixed = |t: f64, sc: f64, ec: f64| -> f64 {
        (bern(0, t) + bern(1, t)) * sc + (bern(4, t) + bern(5, t)) * ec
    };
    let solve_coord = |ps: f64, pe: f64, v1: f64, v2: f64| -> (f64, f64) {
        let r1 = v1 - fixed(t1, ps, pe);
        let r2 = v2 - fixed(t2, ps, pe);
        ((a22 * r1 - a12 * r2) / det, (a11 * r2 - a21 * r1) / det)
    };
    let (ax, bx) = solve_coord(frame.start.x, frame.end.x, p1.x, p2.x);
    let (ay, by) = solve_coord(frame.start.y, frame.end.y, p1.y, p2.y);
    let (az, bz) = solve_coord(frame.start.z, frame.end.z, p1.z, p2.z);
    let cps = [
        frame.start,
        frame.start,
        Point3::new(ax, ay, az),
        Point3::new(bx, by, bz),
        frame.end,
        frame.end,
    ];

    // Geometric residual: every dense sample must lie on the curve
    // (min distance over a dense evaluation — parameterization-free).
    let eval = |t: f64| -> Point3 {
        let mut p = Point3::new(0.0, 0.0, 0.0);
        for (i, cp) in cps.iter().enumerate() {
            let w = bern(i, t);
            p = Point3::new(p.x + w * cp.x, p.y + w * cp.y, p.z + w * cp.z);
        }
        p
    };
    let budget = (tolerance * 10.0).max(1e-9);
    let mut residual = 0.0_f64;
    for &(_, p) in &samples {
        // Coarse scan + windowed refinement (parameterization-free
        // point-to-curve distance; a coarse polyline minimum alone
        // reads the sampling spacing, not the geometry).
        let mut t_min = 0.0_f64;
        let mut d_min = f64::INFINITY;
        for k in 0..=64 {
            let t = k as f64 / 64.0;
            let d = eval(t).distance_to(&p);
            if d < d_min {
                d_min = d;
                t_min = t;
            }
        }
        let mut window = 1.0 / 64.0;
        for _ in 0..12 {
            let lo = (t_min - window).max(0.0);
            let hi = (t_min + window).min(1.0);
            for k in 0..=8 {
                let t = lo + (hi - lo) * k as f64 / 8.0;
                let d = eval(t).distance_to(&p);
                if d < d_min {
                    d_min = d;
                    t_min = t;
                }
            }
            window *= 0.4;
        }
        residual = residual.max(d_min);
    }
    if residual > budget {
        return Err(RimRibbonError::FitResidual {
            side,
            residual,
            budget,
        });
    }
    Ok(cps)
}

/// Ribbon directions for a rim side at its 1/3 and 2/3 stations:
/// `n × t` unit vectors in the neighbour wall's tangent plane,
/// oriented across the rim away from the neighbour face by the
/// centroid witness. These carry the cap's cross-derivative along the
/// side, making the cap tangent to the wall mid-side by construction.
pub(super) fn rim_ribbon_directions(
    brep: &BRep,
    hole: &DetectedHole,
    frame: &FillBoundaryEdge,
    side: usize,
) -> Result<[Vec3; 2], RimRibbonError> {
    let witness = neighbour_witness(brep, hole, side);
    let mut directions = [Vec3::new(0.0, 0.0, 0.0); 2];
    for (station, f) in [(0usize, 1.0 / 3.0), (1usize, 2.0 / 3.0)] {
        let n = interpolate_frame_vector(frame.normal_start, frame.normal_mid, frame.normal_end, f);
        let t =
            interpolate_frame_vector(frame.tangent_start, frame.tangent_mid, frame.tangent_end, f);
        let mut w = n.cross(&t);
        let w_len = w.length();
        if w_len < 1e-12 {
            return Err(RimRibbonError::RibbonDegenerate { side });
        }
        w = w * (1.0 / w_len);
        if let Some(wit) = witness
            && w.dot(&wit) < 0.0
        {
            w = w * -1.0;
        }
        directions[station] = w;
    }
    Ok(directions)
}

/// The neighbour face's surface beside ring side `s`.
pub(super) fn neighbour_surface_geometry<'a>(
    brep: &'a BRep,
    hole: &DetectedHole,
    side: usize,
) -> Option<&'a SurfaceGeometry> {
    let edge = brep.edges.get(*hole.boundary_edges.get(side)?)?;
    let he = brep.half_edges.get(edge.half_edge)?;
    let face = brep.loops.get(he.loop_)?.face;
    Some(&brep.faces.get(face)?.surface.as_ref()?.geometry)
}

/// Witness direction pointing across the rim away from the neighbour
/// face: from the neighbour's outer-loop centroid toward the rim
/// side's midpoint.
fn neighbour_witness(brep: &BRep, hole: &DetectedHole, side: usize) -> Option<Vec3> {
    let edge = brep.edges.get(*hole.boundary_edges.get(side)?)?;
    let he = brep.half_edges.get(edge.half_edge)?;
    let face = brep.loops.get(he.loop_)?.face;
    let positions = brep.face_outer_loop_positions(face);
    if positions.is_empty() {
        return None;
    }
    let n = positions.len() as f64;
    let centroid = Point3::new(
        positions.iter().map(|p| p.x).sum::<f64>() / n,
        positions.iter().map(|p| p.y).sum::<f64>() / n,
        positions.iter().map(|p| p.z).sum::<f64>() / n,
    );
    let curve = edge.curve.as_ref()?;
    let mid = {
        use crate::geometry::curve::CurveEval;
        (&curve.geometry as &dyn CurveEval).point_at(0.5 * (curve.t_start + curve.t_end))
    };
    (mid - centroid).normalized()
}

/// Linear blend of the three sampled frame vectors at fraction `f`.
fn interpolate_frame_vector(start: Vec3, mid: Vec3, end: Vec3, f: f64) -> Vec3 {
    let v = if f <= 0.5 {
        let t = f * 2.0;
        start * (1.0 - t) + mid * t
    } else {
        let t = (f - 0.5) * 2.0;
        mid * (1.0 - t) + end * t
    };
    v.normalized().unwrap_or(mid)
}
