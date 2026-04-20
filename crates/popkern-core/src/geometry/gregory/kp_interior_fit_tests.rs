//! Tests for [`super::fit_interior_cps_to_guide`].
//!
//! Verification goals:
//!
//! 1. **Zero-guide + zero-frame → zero-interior**: degenerate
//!    sanity check.
//! 2. **Full degree-elevated case**: if the pinned frame CPs
//!    equal the degree-elevated guide's CPs, the least-squares
//!    interior should also equal the degree-elevated interior
//!    (recover the guide exactly).
//! 3. **Flat plane**: flat guide + flat frame produces flat
//!    interior.
//! 4. **Idempotency**: calling the fit twice produces the same
//!    result (no drift).

#![allow(clippy::needless_range_loop)]

use super::*;
use crate::geometry::gregory::kp_guide_build::{GuideCornerData, build_guide_from_hermite};
use crate::geometry::{Point3, Vec3};

const EPS: f64 = 1e-8;

fn flat_guide() -> KpGuidePatch {
    let corners = [
        GuideCornerData {
            position: Point3::new(0.0, 0.0, 0.0),
            tangent_u: Vec3::new(1.0, 0.0, 0.0),
            tangent_v: Vec3::new(0.0, 1.0, 0.0),
            twist: Vec3::zero(),
        },
        GuideCornerData {
            position: Point3::new(1.0, 0.0, 0.0),
            tangent_u: Vec3::new(1.0, 0.0, 0.0),
            tangent_v: Vec3::new(0.0, 1.0, 0.0),
            twist: Vec3::zero(),
        },
        GuideCornerData {
            position: Point3::new(1.0, 1.0, 0.0),
            tangent_u: Vec3::new(1.0, 0.0, 0.0),
            tangent_v: Vec3::new(0.0, 1.0, 0.0),
            twist: Vec3::zero(),
        },
        GuideCornerData {
            position: Point3::new(0.0, 1.0, 0.0),
            tangent_u: Vec3::new(1.0, 0.0, 0.0),
            tangent_v: Vec3::new(0.0, 1.0, 0.0),
            twist: Vec3::zero(),
        },
    ];
    build_guide_from_hermite(corners)
}

/// Degree-elevate a bicubic guide to bidegree 9 directly — for
/// comparison with the least-squares fit.
fn elevate_guide_to_bidegree9(guide: &KpGuidePatch) -> RectangularGregoryBiNonic {
    use crate::geometry::gregory::kp_degree_elevate::elevate_cubic_point_to_nine;
    // Elevate each u-row of the guide (4 CPs) to degree 9.
    let mut intermediate = [[Point3::new(0.0, 0.0, 0.0); 10]; 4];
    for k in 0..4 {
        let row: [Point3; 4] = guide.grid[k];
        intermediate[k] = elevate_cubic_point_to_nine(row);
    }
    // Transpose, elevate each column, transpose back.
    let mut grid = [[Point3::new(0.0, 0.0, 0.0); 10]; 10];
    for j in 0..10 {
        let col: [Point3; 4] = [
            intermediate[0][j],
            intermediate[1][j],
            intermediate[2][j],
            intermediate[3][j],
        ];
        let elevated_col = elevate_cubic_point_to_nine(col);
        for i in 0..10 {
            grid[i][j] = elevated_col[i];
        }
    }
    RectangularGregoryBiNonic { grid }
}

fn point_max(a: Point3, b: Point3) -> f64 {
    ((a.x - b.x).abs())
        .max((a.y - b.y).abs())
        .max((a.z - b.z).abs())
}

#[test]
fn zero_guide_zero_frame_produces_zero_interior() {
    let zero_guide = KpGuidePatch::zeros();
    let mut patch = RectangularGregoryBiNonic::zeros();
    fit_interior_cps_to_guide(&mut patch, &zero_guide);
    for i in 3..=6 {
        for j in 3..=6 {
            assert!(
                patch.grid[i][j].x.abs() < EPS
                    && patch.grid[i][j].y.abs() < EPS
                    && patch.grid[i][j].z.abs() < EPS,
                "interior[{i}][{j}] = {:?}",
                patch.grid[i][j]
            );
        }
    }
}

#[test]
fn degree_elevated_frame_recovers_degree_elevated_interior() {
    // If the pinned frame CPs exactly equal the degree-elevated
    // guide's CPs, the least-squares fit for the interior should
    // also match the degree-elevated CPs. This is because the
    // full degree-elevated surface IS the guide, so the L² fit
    // achieves zero residual by setting interior = degree-elevated
    // interior.
    let corners = [
        GuideCornerData {
            position: Point3::new(0.0, 0.0, 0.0),
            tangent_u: Vec3::new(1.0, 0.2, 0.3),
            tangent_v: Vec3::new(0.1, 1.1, -0.2),
            twist: Vec3::new(0.2, -0.1, 0.3),
        },
        GuideCornerData {
            position: Point3::new(1.0, 0.1, 0.2),
            tangent_u: Vec3::new(1.1, -0.2, 0.25),
            tangent_v: Vec3::new(0.15, 1.0, 0.1),
            twist: Vec3::new(-0.15, 0.2, -0.1),
        },
        GuideCornerData {
            position: Point3::new(1.2, 1.1, 0.1),
            tangent_u: Vec3::new(0.9, 0.15, -0.1),
            tangent_v: Vec3::new(-0.1, 0.95, 0.2),
            twist: Vec3::new(0.1, 0.15, -0.2),
        },
        GuideCornerData {
            position: Point3::new(0.1, 1.0, -0.05),
            tangent_u: Vec3::new(1.05, -0.1, 0.15),
            tangent_v: Vec3::new(0.0, 1.05, -0.1),
            twist: Vec3::new(-0.05, -0.2, 0.15),
        },
    ];
    let guide = build_guide_from_hermite(corners);
    let elevated = elevate_guide_to_bidegree9(&guide);
    // Start with a patch whose 64 frame CPs equal the degree-
    // elevated CPs; interior at zero.
    let mut patch = elevated.clone();
    for i in 3..=6 {
        for j in 3..=6 {
            patch.grid[i][j] = Point3::new(0.0, 0.0, 0.0);
        }
    }
    fit_interior_cps_to_guide(&mut patch, &guide);
    // The fitted interior should match the degree-elevated
    // interior to machine precision.
    for i in 3..=6 {
        for j in 3..=6 {
            let diff = point_max(patch.grid[i][j], elevated.grid[i][j]);
            assert!(
                diff < 1e-7,
                "interior[{i}][{j}]: fitted {:?} vs elevated {:?}, diff {diff}",
                patch.grid[i][j],
                elevated.grid[i][j]
            );
        }
    }
}

#[test]
fn flat_guide_and_flat_frame_produce_flat_interior() {
    // Flat guide on z = 0 + flat frame on z = 0 → fit gives
    // interior on z = 0 too.
    let guide = flat_guide();
    // Build a flat frame: set grid[i][j].z = 0 and x, y = linear
    // in u = i/9, v = j/9.
    let mut patch = RectangularGregoryBiNonic::zeros();
    for i in 0..10 {
        for j in 0..10 {
            // Skip interior — it's set by the fit.
            if (3..=6).contains(&i) && (3..=6).contains(&j) {
                continue;
            }
            patch.grid[i][j] = Point3::new(i as f64 / 9.0, j as f64 / 9.0, 0.0);
        }
    }
    fit_interior_cps_to_guide(&mut patch, &guide);
    for i in 3..=6 {
        for j in 3..=6 {
            assert!(
                patch.grid[i][j].z.abs() < EPS,
                "interior[{i}][{j}].z = {}",
                patch.grid[i][j].z
            );
        }
    }
}

#[test]
fn interior_fit_is_idempotent() {
    // Calling the fit twice should give the same result. No
    // drift.
    let corners = [
        GuideCornerData {
            position: Point3::new(0.0, 0.0, 0.0),
            tangent_u: Vec3::new(1.0, 0.0, 0.0),
            tangent_v: Vec3::new(0.0, 1.0, 0.3),
            twist: Vec3::zero(),
        },
        GuideCornerData {
            position: Point3::new(1.0, 0.0, 0.2),
            tangent_u: Vec3::new(1.0, 0.0, 0.0),
            tangent_v: Vec3::new(0.0, 1.0, 0.0),
            twist: Vec3::zero(),
        },
        GuideCornerData {
            position: Point3::new(1.0, 1.0, 0.0),
            tangent_u: Vec3::new(1.0, 0.0, 0.0),
            tangent_v: Vec3::new(0.0, 1.0, 0.0),
            twist: Vec3::zero(),
        },
        GuideCornerData {
            position: Point3::new(0.0, 1.0, 0.15),
            tangent_u: Vec3::new(1.0, 0.0, 0.0),
            tangent_v: Vec3::new(0.0, 1.0, 0.0),
            twist: Vec3::zero(),
        },
    ];
    let guide = build_guide_from_hermite(corners);
    let elevated = elevate_guide_to_bidegree9(&guide);
    let mut patch = elevated.clone();
    // Wipe interior.
    for i in 3..=6 {
        for j in 3..=6 {
            patch.grid[i][j] = Point3::new(0.0, 0.0, 0.0);
        }
    }
    fit_interior_cps_to_guide(&mut patch, &guide);
    let snapshot: [[Point3; 10]; 10] = patch.grid;
    // Wipe interior again and re-fit.
    let mut patch2 = patch.clone();
    for i in 3..=6 {
        for j in 3..=6 {
            patch2.grid[i][j] = Point3::new(0.0, 0.0, 0.0);
        }
    }
    fit_interior_cps_to_guide(&mut patch2, &guide);
    for i in 3..=6 {
        for j in 3..=6 {
            let diff = point_max(patch2.grid[i][j], snapshot[i][j]);
            assert!(diff < 1e-10, "idempotency: [{i}][{j}] drifted by {diff}");
        }
    }
}
