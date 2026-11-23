//! Tangent crossings of a plane and a torus within the tolerance (W3 V
//! cells, tangency design §4.8, PK-1292).
//!
//! Where the plane is tangent to the torus at a hyperbolic point the
//! section crosses itself there; a plane opened from that tangency by
//! less than the tolerance cuts a waist (two branches of one arc passing
//! near) or a gap (two arcs ending near), never a crossing, and the
//! material between the branches is thinner than the tolerance: within
//! it the faces touch, and the result pinches as at the crossing. The
//! crossing is read off the exact gap `g(u) = |A| − K` between the plane
//! and the tube circle at axis angle `u`: a saddle of the faces' gap is
//! a local maximum of `g`, met where the plane normal's radial component
//! `m(u)` is extreme (`u = ψ`, `ψ + π`: a single tangent point lies on
//! the symmetry meridian) or where `F`'s quadratic in `m` has its vertex
//! (a symmetric pair). A maximum with `|g|` within the linear tolerance
//! is spelled as the crossing: the level is shifted there by `g`, over a
//! window flat to second order (so the shifted gap keeps the maximum,
//! now zero: a double root, the branches meeting at the point), and the
//! samples stand halfway between the shifted level and the plane, within
//! `|g| / 2` of both surfaces.

use super::SpiricFrame;
use crate::tolerance::Tolerance;
use popkern_math::PortableMath;
use std::f64::consts::{PI, TAU};

/// One level shift: `shift` at axis angle `u`, falling to zero over
/// `width` either side with zero slope and curvature at `u`.
#[derive(Debug, Clone, Copy)]
pub(super) struct LevelBump {
    u: f64,
    shift: f64,
    width: f64,
}

impl LevelBump {
    /// The shift at `u`: `1 − 6t⁴ + 8t⁶ − 3t⁸` of it at `t = Δu / width`
    /// (one at the crossing, flat to second order there, zero with zero
    /// slope at the window's ends).
    fn at(&self, u: f64) -> f64 {
        let d = u - self.u;
        let t = (d - TAU * (d / TAU).round()) / self.width;
        if t.abs() >= 1.0 {
            return 0.0;
        }
        let t4 = t * t * t * t;
        self.shift * (1.0 - 6.0 * t4 + 8.0 * t4 * t * t - 3.0 * t4 * t4)
    }
}

/// Step of the gap's finite differences in the axis angle.
const STEP: f64 = 1.0e-4;

impl SpiricFrame {
    /// The level shift at `u`.
    pub(super) fn shift_at(&self, u: f64) -> f64 {
        self.bumps.iter().map(|b| b.at(u)).sum()
    }

    /// The largest level shift: the samples lie within half of it of
    /// both surfaces.
    pub(in crate::geometry::intersection) fn largest_shift(&self) -> f64 {
        self.bumps.iter().map(|b| b.shift.abs()).fold(0.0, f64::max)
    }

    /// The exact section's gap at `u`: the plane's distance from the
    /// tube circle there, negative where it cuts the circle.
    fn gap_at(&self, u: f64) -> f64 {
        let m = self.m_at(u);
        let a = self.d0 - self.big_r * m;
        a.abs() - self.r * (m * m + self.na * self.na).sqrt()
    }

    fn gap_slope(&self, u: f64) -> f64 {
        (self.gap_at(u + STEP) - self.gap_at(u - STEP)) / (2.0 * STEP)
    }

    fn gap_curvature(&self, u: f64) -> f64 {
        (self.gap_at(u + STEP) - 2.0 * self.gap_at(u) + self.gap_at(u - STEP)) / (STEP * STEP)
    }

    /// The gap's stationary point near `u0`, by Newton on its slope (a
    /// flat extremum is located by its derivative's root, not its
    /// values); `None` when it wanders off.
    fn gap_stationary(&self, u0: f64) -> Option<f64> {
        let mut u = u0;
        for _ in 0..8 {
            let (s, c) = (self.gap_slope(u), self.gap_curvature(u));
            if c == 0.0 || !c.is_finite() || !s.is_finite() {
                break;
            }
            u -= s / c;
        }
        ((u - u0).abs() < 1.0e-2 && u.is_finite()).then_some(u)
    }

    /// Spell every tangent crossing within the tolerance: shift the level
    /// there and return the section's roots with the crossing's root in
    /// place of the roots its window held (two for a gap, none for a
    /// waist). `None` when there is none to spell.
    pub(in crate::geometry::intersection) fn spell_tangent_crossings(
        &mut self,
        roots: &[f64],
        tolerance: &Tolerance,
    ) -> Option<Vec<f64>> {
        let big_m = self.n1.phypot(self.n2);
        let psi = self.n2.patan2(self.n1);
        let mut candidates = vec![psi, psi + PI];
        let a2 = self.big_r * self.big_r - self.r * self.r;
        if a2.abs() > f64::EPSILON * self.big_r * self.big_r && big_m > 0.0 {
            let vertex = self.d0 * self.big_r / a2;
            if vertex.abs() < big_m {
                let half = (vertex / big_m).pacos();
                candidates.extend(
                    [psi + half, psi - half]
                        .into_iter()
                        .filter_map(|u| self.gap_stationary(u)),
                );
            }
        }
        let mut kept = roots.to_vec();
        let mut spelled = Vec::new();
        for u in candidates {
            let u = crate::geometry::chart::ops::wrap_into(u, 0.0, TAU);
            let near = |x: f64| {
                let d = x - u;
                (d - TAU * (d / TAU).round()).abs()
            };
            if spelled.iter().any(|&b: &f64| near(b) < 1.0e-6) {
                continue;
            }
            let gap = self.gap_at(u);
            let curvature = self.gap_curvature(u);
            // A saddle of the faces' gap: the tube circle's nearest point
            // (a minimum round it) at a maximum along the axis angle.
            if !tolerance.within_linear(gap.abs()) || curvature.is_nan() || curvature >= 0.0 {
                continue;
            }
            let width = (16.0 * gap.abs() / -curvature).sqrt().max(1.0e-12);
            if width > 0.25 {
                continue;
            }
            // A gap's window holds the two roots ending its arcs, a
            // waist's none; anything else is another structure nearby.
            let inside = kept.iter().filter(|&&x| near(x) < width).count();
            if (gap < 0.0 && inside != 0) || inside > 2 {
                continue;
            }
            let m = self.m_at(u);
            let a = self.d0 - self.big_r * m;
            let k = self.r * (m * m + self.na * self.na).sqrt();
            self.bumps.push(LevelBump {
                u,
                shift: a - a.signum() * k,
                width,
            });
            kept.retain(|&x| near(x) >= width);
            kept.push(u);
            spelled.push(u);
        }
        if spelled.is_empty() {
            return None;
        }
        kept.sort_by(f64::total_cmp);
        Some(kept)
    }
}

#[cfg(test)]
#[path = "plane_torus_spiric_crossing_tests.rs"]
mod tests;
