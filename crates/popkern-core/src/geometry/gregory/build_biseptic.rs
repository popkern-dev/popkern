//! Constructors for [`RectangularGregoryBiseptic`].
//!
//! The primary entry point is `from_biquintic_full`, which
//! promotes an existing [`RectangularGregoryBiquinticFull`] by
//! degree-elevating its boundary Bezier curves from 5 to 7 and
//! bilinear-seeding the 32 free interior CPs. The 4 corner
//! inner-CP pairs carry over verbatim — the rational blend is
//! degree-independent.
//!
//! Degree-elevation for a degree-`n` Bernstein control polygon
//! `P_0, …, P_n` to degree `n + 1` is
//! `Q_i = (i/(n+1)) · P_{i-1} + ((n+1-i)/(n+1)) · P_i`
//! with `Q_0 = P_0`, `Q_{n+1} = P_n`. Applied twice produces the
//! degree-5 → degree-7 elevation. For circular-arc boundaries
//! (exact in degree-5 to f64) the elevation is lossless.

use super::rectangular_biquintic_full::RectangularGregoryBiquinticFull;
use super::rectangular_biseptic::RectangularGregoryBiseptic;
use crate::geometry::Point3;

impl RectangularGregoryBiseptic {
    /// Promote a [`RectangularGregoryBiquinticFull`] patch to
    /// biseptic (degree 7) by degree-elevating the four boundary
    /// Bezier curves and bilinear-seeding the 32 free interior
    /// CPs. Preserves corner positions, boundary curves
    /// (lossless degree elevation), and the 4 corner inner-CP
    /// pairs that drive the corner twist-compatibility.
    ///
    /// The bilinear interior seed is refined by the spoke-G1
    /// solver in the calling path; the seed only needs to be
    /// a reasonable starting point.
    pub fn from_biquintic_full(src: &RectangularGregoryBiquinticFull) -> Self {
        // Extract the 4 boundary Bezier curves from the biquintic
        // grid. Each is a degree-5 curve with 6 CPs.
        let boundary_u0: [Point3; 6] = src.grid[0];
        let boundary_u1: [Point3; 6] = src.grid[5];
        let boundary_v0: [Point3; 6] = std::array::from_fn(|i| src.grid[i][0]);
        let boundary_v1: [Point3; 6] = std::array::from_fn(|i| src.grid[i][5]);

        let boundary_u0_7 = elevate_5_to_7(&boundary_u0);
        let boundary_u1_7 = elevate_5_to_7(&boundary_u1);
        let boundary_v0_7 = elevate_5_to_7(&boundary_v0);
        let boundary_v1_7 = elevate_5_to_7(&boundary_v1);

        let mut grid = [[Point3::new(0.0, 0.0, 0.0); 8]; 8];
        grid[0].copy_from_slice(&boundary_u0_7);
        grid[7].copy_from_slice(&boundary_u1_7);
        for i in 0..8 {
            grid[i][0] = boundary_v0_7[i];
            grid[i][7] = boundary_v1_7[i];
        }

        let c00 = grid[0][0];
        let c07 = grid[0][7];
        let c70 = grid[7][0];
        let c77 = grid[7][7];
        #[allow(clippy::needless_range_loop)]
        for i in 1..=6 {
            for j in 1..=6 {
                if matches!((i, j), (1, 1) | (1, 6) | (6, 1) | (6, 6)) {
                    continue;
                }
                let s = (i as f64) / 7.0;
                let t = (j as f64) / 7.0;
                grid[i][j] = Point3::new(
                    (1.0 - s) * (1.0 - t) * c00.x
                        + s * (1.0 - t) * c70.x
                        + s * t * c77.x
                        + (1.0 - s) * t * c07.x,
                    (1.0 - s) * (1.0 - t) * c00.y
                        + s * (1.0 - t) * c70.y
                        + s * t * c77.y
                        + (1.0 - s) * t * c07.y,
                    (1.0 - s) * (1.0 - t) * c00.z
                        + s * (1.0 - t) * c70.z
                        + s * t * c77.z
                        + (1.0 - s) * t * c07.z,
                );
            }
        }

        RectangularGregoryBiseptic {
            grid,
            corner_inner_cps: src.corner_inner_cps,
        }
    }
}

/// Degree-elevate a degree-5 Bernstein control polygon to
/// degree-7 via two successive `n → n + 1` elevations.
fn elevate_5_to_7(p5: &[Point3; 6]) -> [Point3; 8] {
    let p6 = elevate_one(p5.as_slice(), 5);
    let p7 = elevate_one(&p6, 6);
    [p7[0], p7[1], p7[2], p7[3], p7[4], p7[5], p7[6], p7[7]]
}

fn elevate_one(src: &[Point3], n: usize) -> Vec<Point3> {
    let n1 = n + 1;
    let mut out: Vec<Point3> = Vec::with_capacity(n1 + 1);
    out.push(src[0]);
    for i in 1..n1 {
        let a = (i as f64) / (n1 as f64);
        let b = ((n1 - i) as f64) / (n1 as f64);
        out.push(Point3::new(
            a * src[i - 1].x + b * src[i].x,
            a * src[i - 1].y + b * src[i].y,
            a * src[i - 1].z + b * src[i].z,
        ));
    }
    out.push(src[n]);
    out
}
