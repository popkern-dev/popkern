//! [`SpokeG1Patch`] trait implementations for the three Gregory
//! patch types the solver supports (biquintic-full, biseptic,
//! Peters two-stripe). Extracted from `spoke_g1_solver.rs` to
//! keep that file under the 500-line structural cap.

use super::rectangular_biquintic_full::RectangularGregoryBiquinticFull;
use super::rectangular_biseptic::RectangularGregoryBiseptic;
use super::rectangular_gregory_two_stripe::RectangularGregoryTwoStripe;
use super::spoke_g1_solver::{Side, SpokeG1Patch};
use crate::geometry::{Point3, Vec3};

const BIQUINTIC_UZERO: &[(usize, usize)] = &[(1, 2), (1, 3)];
const BIQUINTIC_UONE: &[(usize, usize)] = &[(4, 2), (4, 3)];
const BIQUINTIC_VZERO: &[(usize, usize)] = &[(2, 1), (3, 1)];
const BIQUINTIC_VONE: &[(usize, usize)] = &[(2, 4), (3, 4)];

impl SpokeG1Patch for RectangularGregoryBiquinticFull {
    fn free_indices_for_side(side: Side) -> &'static [(usize, usize)] {
        match side {
            Side::UZero => BIQUINTIC_UZERO,
            Side::UOne => BIQUINTIC_UONE,
            Side::VZero => BIQUINTIC_VZERO,
            Side::VOne => BIQUINTIC_VONE,
        }
    }
    fn cross_tangent_coeff(side: Side, free_i: usize, free_j: usize, t: f64) -> f64 {
        let (inner_index, sign) = match side {
            Side::UZero => (free_j, 1.0),
            Side::UOne => (free_j, -1.0),
            Side::VZero => (free_i, 1.0),
            Side::VOne => (free_i, -1.0),
        };
        sign * 5.0 * bernstein_quintic(inner_index, t)
    }
    fn grid_get(&self, i: usize, j: usize) -> Point3 {
        self.grid[i][j]
    }
    fn grid_set(&mut self, i: usize, j: usize, p: Point3) {
        self.grid[i][j] = p;
    }
    fn partial_u(&self, u: f64, v: f64) -> Vec3 {
        self.partial_u(u, v)
    }
    fn partial_v(&self, u: f64, v: f64) -> Vec3 {
        self.partial_v(u, v)
    }
}

const BISEPTIC_UZERO: &[(usize, usize)] = &[(1, 2), (1, 3), (1, 4), (1, 5)];
const BISEPTIC_UONE: &[(usize, usize)] = &[(6, 2), (6, 3), (6, 4), (6, 5)];
const BISEPTIC_VZERO: &[(usize, usize)] = &[(2, 1), (3, 1), (4, 1), (5, 1)];
const BISEPTIC_VONE: &[(usize, usize)] = &[(2, 6), (3, 6), (4, 6), (5, 6)];

impl SpokeG1Patch for RectangularGregoryBiseptic {
    fn free_indices_for_side(side: Side) -> &'static [(usize, usize)] {
        match side {
            Side::UZero => BISEPTIC_UZERO,
            Side::UOne => BISEPTIC_UONE,
            Side::VZero => BISEPTIC_VZERO,
            Side::VOne => BISEPTIC_VONE,
        }
    }
    fn cross_tangent_coeff(side: Side, free_i: usize, free_j: usize, t: f64) -> f64 {
        let (inner_index, sign) = match side {
            Side::UZero => (free_j, 1.0),
            Side::UOne => (free_j, -1.0),
            Side::VZero => (free_i, 1.0),
            Side::VOne => (free_i, -1.0),
        };
        sign * 7.0 * super::rectangular_biseptic::bernstein_septic(inner_index, t)
    }
    fn grid_get(&self, i: usize, j: usize) -> Point3 {
        self.grid[i][j]
    }
    fn grid_set(&mut self, i: usize, j: usize, p: Point3) {
        self.grid[i][j] = p;
    }
    fn partial_u(&self, u: f64, v: f64) -> Vec3 {
        self.partial_u(u, v)
    }
    fn partial_v(&self, u: f64, v: f64) -> Vec3 {
        self.partial_v(u, v)
    }
}

// Two-stripe: 6 free CPs per side (row 1 j=1..=6 etc.) — full
// Bernstein-basis coverage `{B_1^7..B_6^7}`.
const TWO_STRIPE_UZERO: &[(usize, usize)] = &[(1, 1), (1, 2), (1, 3), (1, 4), (1, 5), (1, 6)];
const TWO_STRIPE_UONE: &[(usize, usize)] = &[(6, 1), (6, 2), (6, 3), (6, 4), (6, 5), (6, 6)];
const TWO_STRIPE_VZERO: &[(usize, usize)] = &[(1, 1), (2, 1), (3, 1), (4, 1), (5, 1), (6, 1)];
const TWO_STRIPE_VONE: &[(usize, usize)] = &[(1, 6), (2, 6), (3, 6), (4, 6), (5, 6), (6, 6)];

impl SpokeG1Patch for RectangularGregoryTwoStripe {
    fn free_indices_for_side(side: Side) -> &'static [(usize, usize)] {
        match side {
            Side::UZero => TWO_STRIPE_UZERO,
            Side::UOne => TWO_STRIPE_UONE,
            Side::VZero => TWO_STRIPE_VZERO,
            Side::VOne => TWO_STRIPE_VONE,
        }
    }
    fn cross_tangent_coeff(side: Side, free_i: usize, free_j: usize, t: f64) -> f64 {
        let (inner_index, sign) = match side {
            Side::UZero => (free_j, 1.0),
            Side::UOne => (free_j, -1.0),
            Side::VZero => (free_i, 1.0),
            Side::VOne => (free_i, -1.0),
        };
        sign * 7.0 * super::rectangular_biseptic::bernstein_septic(inner_index, t)
    }
    fn grid_get(&self, i: usize, j: usize) -> Point3 {
        self.grid[i][j]
    }
    fn grid_set(&mut self, i: usize, j: usize, p: Point3) {
        self.grid[i][j] = p;
    }
    fn partial_u(&self, u: f64, v: f64) -> Vec3 {
        self.partial_u(u, v)
    }
    fn partial_v(&self, u: f64, v: f64) -> Vec3 {
        self.partial_v(u, v)
    }
}

#[inline]
fn bernstein_quintic(i: usize, t: f64) -> f64 {
    let s = 1.0 - t;
    match i {
        0 => s * s * s * s * s,
        1 => 5.0 * t * s * s * s * s,
        2 => 10.0 * t * t * s * s * s,
        3 => 10.0 * t * t * t * s * s,
        4 => 5.0 * t * t * t * t * s,
        5 => t * t * t * t * t,
        _ => 0.0,
    }
}
