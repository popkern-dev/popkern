//! Constructors for [`RectangularGregoryBiquinticFull`].
//!
//! Extracted from `build.rs` to keep the latter under the 500-
//! line structural cap. Mirrors the deg5 constructor semantics
//! with two additions: the 12 free interior CPs are seeded from
//! bilinear interpolation (to be refined by the Phase 6 spoke
//! G1 solver), and `from_deg5` promotes an existing deg5 patch
//! in place.

use super::build::GregoryBuildError;
use super::rectangular_biquintic_full::RectangularGregoryBiquinticFull;
use super::rectangular_deg5::RectangularGregoryDeg5;
use crate::geometry::Point3;

impl RectangularGregoryBiquinticFull {
    /// Build a biquintic-full rectangular Gregory patch from 4
    /// quintic-Bezier boundary curves and 4 inner CP pairs (same
    /// corner convention as [`RectangularGregoryDeg5`]).
    ///
    /// The 12 free interior CPs at the non-corner-adjacent
    /// interior Bezier indices are seeded by **bilinear
    /// interpolation** of the 4 corner positions at the Bezier
    /// index `(i/5, j/5)`. This is a no-twist baseline suitable
    /// for the pre-solver path; the spoke-G1 joint solver
    /// replaces these with values that enforce cross-tangent
    /// continuity along shared spokes.
    ///
    /// The 4 corner-adjacent positions `(1,1), (1,4), (4,1),
    /// (4,4)` remain placeholders — evaluation picks them up
    /// from the inner-CP blend.
    pub fn from_boundary_beziers(
        boundaries: [[Point3; 6]; 4],
        corner_inner_cps: [(Point3, Point3); 4],
        tol_linear: f64,
    ) -> Result<RectangularGregoryBiquinticFull, GregoryBuildError> {
        let deg5 = RectangularGregoryDeg5::from_boundary_beziers(
            boundaries,
            corner_inner_cps,
            tol_linear,
        )?;
        Ok(Self::from_deg5(&deg5))
    }

    /// Promote an existing [`RectangularGregoryDeg5`] to the
    /// biquintic-full form. The 4 corner inner-CP pairs are
    /// preserved verbatim; the 12 free interior CPs are seeded
    /// from the bilinear interpolation of the 4 corner positions
    /// at Bezier index `(i/5, j/5)`. This keeps the corner
    /// tangent-plane behaviour identical to the source deg5
    /// patch while introducing the 12 extra DOF needed by the
    /// spoke-G1 solver.
    pub fn from_deg5(deg5: &RectangularGregoryDeg5) -> RectangularGregoryBiquinticFull {
        let mut grid = deg5.boundary_cps;
        let c00 = grid[0][0];
        let c05 = grid[0][5];
        let c55 = grid[5][5];
        let c50 = grid[5][0];
        #[allow(clippy::needless_range_loop)]
        for i in 1..=4 {
            for j in 1..=4 {
                if matches!((i, j), (1, 1) | (1, 4) | (4, 1) | (4, 4)) {
                    continue;
                }
                let s = (i as f64) / 5.0;
                let t = (j as f64) / 5.0;
                grid[i][j] = Point3::new(
                    (1.0 - s) * (1.0 - t) * c00.x
                        + s * (1.0 - t) * c50.x
                        + s * t * c55.x
                        + (1.0 - s) * t * c05.x,
                    (1.0 - s) * (1.0 - t) * c00.y
                        + s * (1.0 - t) * c50.y
                        + s * t * c55.y
                        + (1.0 - s) * t * c05.y,
                    (1.0 - s) * (1.0 - t) * c00.z
                        + s * (1.0 - t) * c50.z
                        + s * t * c55.z
                        + (1.0 - s) * t * c05.z,
                );
            }
        }
        RectangularGregoryBiquinticFull {
            grid,
            corner_inner_cps: deg5.inner_cps,
        }
    }
}
