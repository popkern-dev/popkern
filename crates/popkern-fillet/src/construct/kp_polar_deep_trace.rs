//! Deep-trace diagnostic helper for the polar bi-6 path.
//!
//! Gated on `POPKERN_KP_POLAR_DEEP_TRACE=1`. For every sub-patch
//! in the committed (or rejected) cap it prints:
//! - `self_intersects_with_tolerance` result after bi-9 elevation
//! - minimum `|∂u × ∂v|` sampled on a 10×10 interior grid (0 = fold)
//! - CP coordinates at key positions: `grid[0][0]` (star),
//!   `grid[3][3]` (centre), `grid[6][6]` (outer corner).
//!
//! **Key finding surfaced by this helper (2026-04-21, refined
//! 2026-05-31):** the K-P 2015 cap assembly only populates rows
//! 3..7 at columns 3..7 via caller-supplied `interior_filler`,
//! and the F_5 solver only writes `(2, 3)` and `(3, 3)`. The
//! remaining deep-interior CPs now default to the cap star rather
//! than the origin, but N=6 still folds because the open-interface
//! row and deep corner are not paper-authored construction data.

#![allow(dead_code)]

use popkern_core::geometry::Point3;
use popkern_core::geometry::gregory::{
    GregoryPatch, RectangularGregoryBiNonic, elevate_bi6_grid_to_bi9,
};

const SELF_INT_TOL: f64 = 1e-6;

pub(super) fn deep_trace_bi6_output(patches_bi6: &[[[Point3; 7]; 7]], n: usize) {
    eprintln!("[kp_polar_bi6 deep-trace] N={n} per-sub-patch diagnostics:");
    for (k, &bi6) in patches_bi6.iter().enumerate().take(n) {
        let bi9 = elevate_bi6_grid_to_bi9(bi6);
        let gp =
            GregoryPatch::RectangularBiNonic(Box::new(RectangularGregoryBiNonic { grid: bi9 }));
        let folded = gp.self_intersects_with_tolerance(SELF_INT_TOL);
        let mut min_mag = f64::INFINITY;
        let mut min_at = (0.0, 0.0);
        for ui in 0..10 {
            for vi in 0..10 {
                let u = 0.05 + 0.90 * (ui as f64) / 9.0;
                let v = 0.05 + 0.90 * (vi as f64) / 9.0;
                let du = gp.partial_u(u, v);
                let dv = gp.partial_v(u, v);
                let c = du.cross(&dv);
                let m = (c.x * c.x + c.y * c.y + c.z * c.z).sqrt();
                if m < min_mag {
                    min_mag = m;
                    min_at = (u, v);
                }
            }
        }
        let star = bi6[0][0];
        let mid = bi6[3][3];
        let outer = bi6[6][6];
        eprintln!(
            "  sub-patch {k}: folded={folded}, min|du×dv|={min_mag:.3e} at ({:.2},{:.2}) \
             star=({:.3},{:.3},{:.3}) mid=({:.3},{:.3},{:.3}) outer=({:.3},{:.3},{:.3})",
            min_at.0,
            min_at.1,
            star.x,
            star.y,
            star.z,
            mid.x,
            mid.y,
            mid.z,
            outer.x,
            outer.y,
            outer.z,
        );
    }
}
