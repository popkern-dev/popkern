//! The tie ladder of radial ordering.
//!
//! Fins whose faces leave an edge in the same first-order direction
//! (within the angular tolerance) are tied: tangent faces, a cylinder
//! resting on a plane, a torus lying on a plane. The ladder separates
//! them in order of rigour:
//!
//! 1. **Curvature.** The plane normal to the edge tangent `t` holds the
//!    fin's inward direction `d` and the face normal (the edge lies in
//!    the face's tangent plane), so the face's cross-section there is
//!    its normal section along `d`, of curvature `κ_n(d) = II(d,d) /
//!    I(d,d)`. Signed toward `m = t × d` (counter-clockwise about `t`),
//!    it is the rate at which the face turns away from `d`: plane 0,
//!    a cylinder of radius r resting on it ±1/r. Exact from the
//!    surface's second derivatives.
//! 2. **Offset probe.** Fins of equal curvature (faces osculating to
//!    second order) are followed a short way into their faces and the
//!    probe points compared in the normal plane: the higher-order
//!    separation.
//! 3. **Otherwise** the fins are coincident faces nothing merged, and
//!    ordering declines with the keys it measured. It never guesses.
//!
//! Each rung can be switched off (`TieLadder`) so that negative
//! controls can show the rung is what orders a fixture.

use crate::geometry::Vec3;
use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};
use crate::healing::curvature::fundamental_forms::{
    first_fundamental_form, second_fundamental_form,
};

/// The rungs of the tie ladder in use. Both are on in production.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct TieLadder {
    /// Separate first-order ties by cross-section curvature.
    pub curvature: bool,
    /// Separate remaining ties by the offset probe.
    pub offset_probe: bool,
}

impl Default for TieLadder {
    fn default() -> Self {
        Self {
            curvature: true,
            offset_probe: true,
        }
    }
}

/// What the ladder measured for one tied fin.
#[derive(Debug, Clone, PartialEq)]
pub struct TieEvidence {
    /// The fin, as an index into the source's fins.
    pub fin: usize,
    /// Its first-order angle about the tangent from the station's
    /// reference fin, in [0, 2π) (or just below zero for a tie
    /// wrapping past 2π).
    pub angle: f64,
    /// Its signed cross-section curvature, where the rung ran.
    pub curvature: Option<f64>,
    /// The signed angle of its offset probe point from `angle`, where
    /// the rung ran.
    pub probe: Option<f64>,
}

/// Signed curvature of the face's cross-section by the plane normal to
/// `t` at `(u, v)`, positive when it bends toward `t × d`. `None` where
/// the chart is singular or `d` is not tangent to the face.
pub fn cross_section_curvature(
    surface: &SurfaceGeometry,
    (u, v): (f64, f64),
    d: Vec3,
    t: Vec3,
) -> Option<f64> {
    let ff1 = first_fundamental_form(surface, u, v);
    let ff2 = second_fundamental_form(surface, u, v)?;
    let det = ff1.determinant();
    if det.is_nan() || det.abs() <= 1e-30 {
        return None;
    }
    // d = a·∂u + b·∂v (least squares onto the tangent plane).
    let (p, q) = (
        d.dot(&surface.partial_u(u, v)),
        d.dot(&surface.partial_v(u, v)),
    );
    let a = (ff1.g * p - ff1.f * q) / det;
    let b = (ff1.e * q - ff1.f * p) / det;
    let first = ff1.e * a * a + 2.0 * ff1.f * a * b + ff1.g * b * b;
    if first.is_nan() || first <= 0.0 {
        return None;
    }
    let kappa = (ff2.l * a * a + 2.0 * ff2.m * a * b + ff2.n_val * b * b) / first;
    let n = surface.normal_at(u, v).normalized()?;
    Some(kappa * n.dot(&t.cross(&d)))
}

/// Split `members` (sorted by `key`) into runs whose consecutive keys
/// differ by at most `tie`.
pub fn runs<T>(members: Vec<T>, key: impl Fn(&T) -> f64, tie: f64) -> Vec<Vec<T>> {
    let mut out: Vec<Vec<T>> = Vec::new();
    for m in members {
        match out.last_mut() {
            Some(run) if run.last().is_some_and(|l| key(&m) - key(l) <= tie) => run.push(m),
            _ => out.push(vec![m]),
        }
    }
    out
}

#[cfg(test)]
#[path = "ladder_tests.rs"]
mod tests;
