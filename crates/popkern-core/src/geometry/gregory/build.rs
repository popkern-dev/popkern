//! Constructors for rational Gregory patches from boundary + twist data.
//!
//! The Gregory evaluator in [`super::rectangular`] / [`super::triangular`]
//! is a pure function of control-point inputs. This module provides the
//! construction helpers that translate higher-level geometric data —
//! four (or three) cubic-Bezier boundary curves plus optional per-corner
//! twist-compatibility adjustments — into the flat CP arrays the
//! evaluators consume.
//!
//! # Input convention
//!
//! For the rectangular (N=4) form, the 4 boundary Bezier curves are
//! supplied in counter-clockwise order starting from the `(u=0, v=0)`
//! corner:
//!
//! ```text
//!                          corner 1
//!   boundaries[0]   ┌────────────┐ boundaries[1]
//!   (u = 0 edge)    │            │ (v = 1 edge)
//!                   │            │
//!                   │            │
//!                   │            │
//!   corner 0        └────────────┘ corner 2
//!                  boundaries[3]    corner 3
//!                  (v = 0 edge)
//!                     boundaries[2] (u = 1 edge)
//! ```
//!
//! Each `boundaries[k]` is a 4-tuple of `Point3` giving the cubic
//! Bezier control polygon of that edge, **oriented in the
//! counter-clockwise traversal direction** — so:
//!
//! - `boundaries[0]` runs from corner 0 → corner 1.
//! - `boundaries[1]` runs from corner 1 → corner 2.
//! - `boundaries[2]` runs from corner 2 → corner 3.
//! - `boundaries[3]` runs from corner 3 → corner 0.
//!
//! Shared corners must agree: `boundaries[0][3] == boundaries[1][0]`
//! etc. — the constructor enforces this within `tolerance.linear()`.

use super::rectangular::RectangularGregory;
use super::triangular::TriangularGregory;
use crate::geometry::Point3;

/// Error returned by a Gregory-patch constructor when its inputs are
/// inconsistent.
#[derive(Debug, Clone, PartialEq)]
pub enum GregoryBuildError {
    /// Two boundary curves disagreed on their shared corner position
    /// by more than the supplied tolerance.
    CornerMismatch {
        /// Which corner index (0..N) the mismatch occurred at.
        corner: usize,
        /// Observed distance between the two candidate positions.
        distance: f64,
    },
}

impl std::fmt::Display for GregoryBuildError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            GregoryBuildError::CornerMismatch { corner, distance } => write!(
                f,
                "Gregory patch corner {corner} mismatch: boundary endpoints \
                 differ by {distance:.6e}",
            ),
        }
    }
}

impl std::error::Error for GregoryBuildError {}

impl RectangularGregory {
    /// Build a rectangular Gregory patch from 4 cubic-Bezier boundary
    /// curves plus 4 inner control-point pairs.
    ///
    /// Each `boundaries[k]` is a cubic Bezier control polygon
    /// `[P_0, P_1, P_2, P_3]` along edge `k`, oriented CCW around the
    /// patch (edge 0 → 1 → 2 → 3 → 0).
    ///
    /// Shared-corner coincidence is enforced to within `tol_linear`:
    /// `boundaries[k][3] == boundaries[(k + 1) % 4][0]` for every
    /// `k`. Violation returns [`GregoryBuildError::CornerMismatch`].
    ///
    /// # Corner ordering
    ///
    /// `inner_cps[k]` corresponds to the patch corner at Bezier index
    /// `(0, 0), (0, 3), (3, 3), (3, 0)` for `k = 0, 1, 2, 3`
    /// respectively — matching the `(u = 0, v = 0) → (u = 0, v = 1)
    /// → (u = 1, v = 1) → (u = 1, v = 0)` CCW traversal.
    pub fn from_boundary_beziers(
        boundaries: [[Point3; 4]; 4],
        inner_cps: [(Point3, Point3); 4],
        tol_linear: f64,
    ) -> Result<RectangularGregory, GregoryBuildError> {
        // Corner consistency check: each corner is shared by two
        // consecutive boundaries. We check all 4 shared-corner pairs
        // before assembling the grid.
        for k in 0..4 {
            let end = boundaries[k][3];
            let start = boundaries[(k + 1) % 4][0];
            let dx = end.x - start.x;
            let dy = end.y - start.y;
            let dz = end.z - start.z;
            let dist = (dx * dx + dy * dy + dz * dz).sqrt();
            if dist > tol_linear {
                // The "shared corner" lands between boundaries `k`
                // and `(k + 1) % 4`. The corner index in our CCW
                // ordering is `(k + 1) % 4` (the starting corner of
                // the *next* boundary).
                return Err(GregoryBuildError::CornerMismatch {
                    corner: (k + 1) % 4,
                    distance: dist,
                });
            }
        }

        // Assemble the 4×4 boundary-CP grid. Recall the convention:
        // boundary_cps[i][j] is at Bezier index (i, j) where i is
        // the u-index and j is the v-index.
        //
        // Corner-index-to-Bezier-index map (CCW):
        //   corner 0 → (0, 0)
        //   corner 1 → (0, 3)
        //   corner 2 → (3, 3)
        //   corner 3 → (3, 0)
        //
        // Edges:
        //   boundary 0 (corner 0 → corner 1) walks v from 0 to 3 at u=0
        //       → boundary_cps[0][0..4] = boundaries[0]
        //   boundary 1 (corner 1 → corner 2) walks u from 0 to 3 at v=3
        //       → boundary_cps[0..4][3] = boundaries[1]
        //   boundary 2 (corner 2 → corner 3) walks v from 3 to 0 at u=3
        //       → boundary_cps[3][0..4] in reverse = boundaries[2]
        //       → boundary_cps[3][3 - j] = boundaries[2][j]
        //   boundary 3 (corner 3 → corner 0) walks u from 3 to 0 at v=0
        //       → boundary_cps[0..4][0] in reverse = boundaries[3]
        //       → boundary_cps[3 - i][0] = boundaries[3][i]
        let mut grid = [[Point3::new(0.0, 0.0, 0.0); 4]; 4];
        for j in 0..4 {
            grid[0][j] = boundaries[0][j];
        }
        for i in 0..4 {
            grid[i][3] = boundaries[1][i];
        }
        for j in 0..4 {
            grid[3][3 - j] = boundaries[2][j];
        }
        for i in 0..4 {
            grid[3 - i][0] = boundaries[3][i];
        }

        Ok(RectangularGregory {
            boundary_cps: grid,
            inner_cps,
        })
    }

    /// Build a patch from the four boundary Beziers using a plain
    /// bicubic-Bezier inner layer — both candidates in each pair are
    /// set to the **same** interior point. The interior point at
    /// Bezier index `(i, j)` is chosen to match a flat bilinear
    /// interpolation of the four corner points; this reproduces an
    /// ordinary bicubic Bezier and is the correct choice when the
    /// four boundaries carry no corner twist incompatibility.
    ///
    /// Useful as a baseline: any deviation from a plain bicubic
    /// (e.g., from incident G1-continuous neighbours with mismatched
    /// twist planes) must be injected via explicit `inner_cps` with
    /// unequal candidates.
    pub fn from_boundary_beziers_no_twist(
        boundaries: [[Point3; 4]; 4],
        tol_linear: f64,
    ) -> Result<RectangularGregory, GregoryBuildError> {
        // Start with the boundary grid assembled via the main
        // constructor, then compute the 4 inner-CP pairs from the
        // four corners by bilinear interpolation — a standard
        // bicubic-Bezier interior that has no twist adjustment.
        let tmp = RectangularGregory::from_boundary_beziers(
            boundaries,
            [(Point3::new(0.0, 0.0, 0.0), Point3::new(0.0, 0.0, 0.0)); 4],
            tol_linear,
        )?;
        let c00 = tmp.boundary_cps[0][0];
        let c03 = tmp.boundary_cps[0][3];
        let c33 = tmp.boundary_cps[3][3];
        let c30 = tmp.boundary_cps[3][0];
        // Bilinear interpolation for interior Bezier CPs at
        // (i, j) = (1, 1), (1, 2), (2, 1), (2, 2), with parameters
        // (i / 3, j / 3) in [0, 1]:
        //   P(s, t) = (1-s)(1-t)·c00 + s(1-t)·c30 + st·c33 + (1-s)t·c03
        let lerp = |s: f64, t: f64| -> Point3 {
            Point3::new(
                (1.0 - s) * (1.0 - t) * c00.x
                    + s * (1.0 - t) * c30.x
                    + s * t * c33.x
                    + (1.0 - s) * t * c03.x,
                (1.0 - s) * (1.0 - t) * c00.y
                    + s * (1.0 - t) * c30.y
                    + s * t * c33.y
                    + (1.0 - s) * t * c03.y,
                (1.0 - s) * (1.0 - t) * c00.z
                    + s * (1.0 - t) * c30.z
                    + s * t * c33.z
                    + (1.0 - s) * t * c03.z,
            )
        };
        let p11 = lerp(1.0 / 3.0, 1.0 / 3.0);
        let p12 = lerp(1.0 / 3.0, 2.0 / 3.0);
        let p22 = lerp(2.0 / 3.0, 2.0 / 3.0);
        let p21 = lerp(2.0 / 3.0, 1.0 / 3.0);
        // inner_cps[k] corresponds to the inner CP reachable from
        // corner k — (1, 1) near corner 0, (1, 2) near corner 1,
        // (2, 2) near corner 2, (2, 1) near corner 3 — with both
        // candidates in the pair equal to the plain bicubic value.
        Ok(RectangularGregory {
            boundary_cps: tmp.boundary_cps,
            inner_cps: [(p11, p11), (p12, p12), (p22, p22), (p21, p21)],
        })
    }
}

impl TriangularGregory {
    /// Build a triangular Gregory patch from 3 cubic-Bezier boundary
    /// curves plus 3 inner control-point candidates.
    ///
    /// Each `boundaries[k]` is a cubic Bezier control polygon of the
    /// edge opposite corner `k`: boundary 0 runs from corner 1 to
    /// corner 2, boundary 1 from corner 2 to corner 0, boundary 2
    /// from corner 0 to corner 1. Shared endpoints must agree to
    /// within `tol_linear`.
    ///
    /// `inner_cps[k]` is the Longhi interior-CP candidate associated
    /// with corner `k`; it drives the cross-derivative along the two
    /// edges incident to that corner.
    pub fn from_boundary_beziers(
        boundaries: [[Point3; 4]; 3],
        inner_cps: [Point3; 3],
        tol_linear: f64,
    ) -> Result<TriangularGregory, GregoryBuildError> {
        // Corner consistency:
        //   corner 0 = boundaries[2][0] = boundaries[1][3]  (edge 2 starts + edge 1 ends at corner 0)
        //   corner 1 = boundaries[0][0] = boundaries[2][3]
        //   corner 2 = boundaries[1][0] = boundaries[0][3]
        let pairs = [
            (0usize, boundaries[2][0], boundaries[1][3]),
            (1usize, boundaries[0][0], boundaries[2][3]),
            (2usize, boundaries[1][0], boundaries[0][3]),
        ];
        for (corner, a, b) in pairs {
            let dx = a.x - b.x;
            let dy = a.y - b.y;
            let dz = a.z - b.z;
            let dist = (dx * dx + dy * dy + dz * dz).sqrt();
            if dist > tol_linear {
                return Err(GregoryBuildError::CornerMismatch {
                    corner,
                    distance: dist,
                });
            }
        }

        // Storage slots in TriangularGregory::boundary_cps:
        //   [P_300, P_030, P_003, P_210, P_120, P_021, P_012, P_102, P_201]
        // Bezier index (i, j, k) with i+j+k=3, barycentric coords (u, v, w).
        //
        // Corner indices in our CCW ordering map to (i, j, k):
        //   corner 0 → (3, 0, 0)
        //   corner 1 → (0, 3, 0)
        //   corner 2 → (0, 0, 3)
        //
        // Edge 0 (opposite corner 0) connects corners 1 and 2 with
        // interior CPs P_021 (j=2, k=1) and P_012 (j=1, k=2).
        //   boundaries[0] = [corner 1, P_021, P_012, corner 2]
        //                 = [P_030, P_021, P_012, P_003]
        //
        // Edge 1 (opposite corner 1) connects corners 2 and 0 with
        // interior CPs P_102 (i=1, k=2) and P_201 (i=2, k=1).
        //   boundaries[1] = [corner 2, P_102, P_201, corner 0]
        //                 = [P_003, P_102, P_201, P_300]
        //
        // Edge 2 (opposite corner 2) connects corners 0 and 1 with
        // interior CPs P_210 (i=2, j=1) and P_120 (i=1, j=2).
        //   boundaries[2] = [corner 0, P_210, P_120, corner 1]
        //                 = [P_300, P_210, P_120, P_030]
        let boundary_cps = [
            boundaries[2][0], // P_300 (corner 0)
            boundaries[0][0], // P_030 (corner 1)
            boundaries[1][0], // P_003 (corner 2)
            boundaries[2][1], // P_210
            boundaries[2][2], // P_120
            boundaries[0][1], // P_021
            boundaries[0][2], // P_012
            boundaries[1][1], // P_102
            boundaries[1][2], // P_201
        ];

        Ok(TriangularGregory {
            boundary_cps,
            inner_cps,
        })
    }
}
