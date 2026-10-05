//! Per-sample seam continuity profile (PK-575): where along a seam
//! the G1/G2 residual lives, and in which tangent direction.
//!
//! The evidence instrument ([`crate::blend_continuity_evidence`])
//! keeps one worst number per seam. Driving a residual down needs its
//! location and shape: the fraction along the seam, and the full jump
//! of the second fundamental form across it. The instrument samples
//! the curvature jump along, across and on ONE diagonal. Those three
//! signed readings determine the 2×2 jump form `ΔII` on the shared
//! tangent plane (basis: seam direction `e₁`, transverse `e₂`):
//!
//! ```text
//! ΔII = [[a, b], [b, c]],  a = Δκ(e₁), c = Δκ(e₂),
//! b = Δκ((e₁+e₂)/√2) − (a + c)/2
//! ```
//!
//! Its spectral norm `|a+c|/2 + √(((a−c)/2)² + b²)` is the largest
//! normal-curvature jump over EVERY tangent direction. The three
//! sampled directions alone can under-read it, e.g. a pure twist jump
//! peaks on the other diagonal.

use popkern_core::geometry::Point3;
use popkern_core::topology::{BRep, EdgeHandle};
use popkern_math::PortableMath;

use crate::cap_seam_measure::{edge_seam_faces, resolve_edge_seam, sample_edge_seam};

/// One sample of a seam's continuity profile.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SeamProfileSample {
    /// Fraction along the seam, in edge-curve direction.
    pub fraction: f64,
    /// The seam point (first side's image).
    pub point: Point3,
    /// Worse side's distance to the edge curve.
    pub g0: f64,
    /// `|n_a × n_b|`, when both sides give a normal.
    pub g1: Option<f64>,
    /// Signed normal-curvature jump along the seam.
    pub along: Option<f64>,
    /// Signed normal-curvature jump transverse to the seam.
    pub transverse: Option<f64>,
    /// Signed normal-curvature jump on the `along + transverse`
    /// diagonal.
    pub diagonal: Option<f64>,
    /// Spectral norm of the second-fundamental-form jump: the largest
    /// normal-curvature jump over every tangent direction (needs all
    /// three directional readings).
    pub form_jump: Option<f64>,
    /// Angle (degrees, in `(-90, 90]`) from the seam direction to the
    /// tangent direction where the form jump peaks.
    pub peak_direction_deg: Option<f64>,
}

impl SeamProfileSample {
    /// The largest of the three SAMPLED directional jumps (what the
    /// evidence instrument's transverse term reads, plus along).
    pub fn sampled_jump(&self) -> Option<f64> {
        [self.along, self.transverse, self.diagonal]
            .into_iter()
            .flatten()
            .map(f64::abs)
            .reduce(f64::max)
    }
}

/// The 2×2 jump form from signed along, transverse and diagonal
/// readings: `(spectral norm, peak direction in degrees from along)`.
pub fn second_form_jump(along: f64, transverse: f64, diagonal: f64) -> (f64, f64) {
    let (a, c) = (along, transverse);
    let b = diagonal - 0.5 * (a + c);
    let mean = 0.5 * (a + c);
    let radius = (0.25 * (a - c) * (a - c) + b * b).sqrt();
    let (l1, l2) = (mean + radius, mean - radius);
    let (norm, lambda) = if l1.abs() >= l2.abs() {
        (l1.abs(), l1)
    } else {
        (l2.abs(), l2)
    };
    // Eigenvector of [[a, b], [b, c]] for `lambda`.
    let (x, y) = if b.abs() > 1e-300 {
        (b, lambda - a)
    } else if (lambda - a).abs() <= (lambda - c).abs() {
        (1.0, 0.0)
    } else {
        (0.0, 1.0)
    };
    let mut deg = y.patan2(x).to_degrees();
    if deg <= -90.0 {
        deg += 180.0;
    } else if deg > 90.0 {
        deg -= 180.0;
    }
    (norm, deg)
}

/// Sample `edge`'s continuity profile at `samples` fractions evenly
/// spaced over `[end_margin, 1 − end_margin]`, with the evidence
/// instrument's readings. The instrument itself insets by 5%; a
/// profile at `end_margin = 0` reads the seam's corner bands too,
/// where vertex-blend corner residuals live. `None` when the edge is
/// not a two-face seam with pcurves and surfaces on both sides.
pub fn seam_profile(
    brep: &BRep,
    edge: EdgeHandle,
    samples: usize,
    end_margin: f64,
) -> Option<Vec<SeamProfileSample>> {
    let jets = crate::cap_seam_side::SurfaceJetCache::default();
    let edge_data = brep.edges.get(edge)?;
    let adjacency = edge_seam_faces(brep, edge_data)?;
    let seam = resolve_edge_seam(brep, edge_data, &adjacency, &jets)?;
    Some(
        sample_edge_seam(&seam, samples, end_margin.clamp(0.0, 0.49))
            .into_iter()
            .map(|r| {
                let form = match (r.jump_along, r.jump_perp, r.jump_diag) {
                    (Some(a), Some(c), Some(d)) => Some(second_form_jump(a, c, d)),
                    _ => None,
                };
                SeamProfileSample {
                    fraction: r.fraction,
                    point: r.point,
                    g0: r.g0,
                    g1: r.g1,
                    along: r.jump_along,
                    transverse: r.jump_perp,
                    diagonal: r.jump_diag,
                    form_jump: form.map(|f| f.0),
                    peak_direction_deg: form.map(|f| f.1),
                }
            })
            .collect(),
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn form_jump_recovers_axis_aligned_and_twist_jumps() {
        // Pure transverse jump 1/r (a rolling-ball rail): peak across.
        let (n, deg) = second_form_jump(0.0, 2.0, 1.0);
        assert!((n - 2.0).abs() < 1e-12 && (deg.abs() - 90.0).abs() < 1e-9);
        // Pure along jump.
        let (n, deg) = second_form_jump(3.0, 0.0, 1.5);
        assert!((n - 3.0).abs() < 1e-12 && deg.abs() < 1e-9);
        // Pure twist b = 1: along/transverse read 0, the sampled
        // diagonal reads +1, the other diagonal −1; the norm is 1 at
        // ±45°.
        let (n, deg) = second_form_jump(0.0, 0.0, 1.0);
        assert!((n - 1.0).abs() < 1e-12 && (deg - 45.0).abs() < 1e-9);
        // A twist the SAMPLED directions under-read: a = 1, c = −1,
        // diagonal 0 → b = 0, norm 1 (sampled max 1: agreement); but
        // a = c = 0.5, diagonal = −0.5 → b = −1, eigenvalues 1.5 and
        // −0.5: the norm 1.5 exceeds every sampled |jump| (0.5).
        let (n, _) = second_form_jump(0.5, 0.5, -0.5);
        assert!((n - 1.5).abs() < 1e-12);
    }
}
