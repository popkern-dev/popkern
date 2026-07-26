//! Measured cap-quality postcondition for the exact-boundary route.
//!
//! # Why this exists
//!
//! Acceptance used to read
//!
//! ```text
//! let accepted = committed && diagnostics.is_empty() && pcurve_pairs > 0;
//! ```
//!
//! — i.e. a candidate was accepted when nothing had *complained*, not
//! when anything had been *measured*. Two things follow from that, and
//! both were observed:
//!
//! - Recording an informational diagnostic vetoed the route, so
//!   observability was penalised by construction.
//! - A cap could ship folded as long as it stayed silent. Measured on
//!   `n_gonal_pyramid` mixed-radius apexes, the committed Coons fan
//!   ran to 89.9° of normal misalignment across its own internal
//!   seams while the fallback route managed 9.8° — the worse cap won
//!   because it complained less (PK-197, PK-201).
//!
//! This module adds the missing measurement: after commit, check that
//! the cap is actually continuous across its own seams, and decline if
//! it is not so a better route can serve.
//!
//! # Where it samples, and why that matters
//!
//! Only the **middle** of each seam. Adjacent cap sectors meet at a
//! pole where one partial derivative vanishes; a finite-difference
//! normal there is direction noise, and gating on it would reject
//! good caps for a chart artifact. The distinction is not cosmetic —
//! it is what separates the two routes here. The Coons fan is bad
//! *mid-seam* (46.7° at n=3, 60.2° at n=4), which is a real crease a
//! user sees; the Gregory fallback is ~0° mid-seam and only degrades
//! at the very ends.

use popkern_core::topology::{BRep, FaceHandle};

use crate::cap_seam_continuity::seam_normal_profile;

/// Fractions along each seam that count as "the middle".
const MID_SEAM_FRACTIONS: [f64; 3] = [0.4, 0.5, 0.6];

/// Worst mid-seam normal misalignment, as `|n_a x n_b|`, over every
/// seam shared by two of `cap_faces`.
///
/// Returns `None` when the cap has no internal seam to measure — a
/// one-patch cap is not "perfectly continuous", it is simply not
/// described by this metric, and callers must not read a missing
/// measurement as a pass.
pub(crate) fn worst_mid_seam_misalignment(brep: &BRep, cap_faces: &[FaceHandle]) -> Option<f64> {
    let mut worst: Option<f64> = None;
    for (i, &a) in cap_faces.iter().enumerate() {
        for &b in &cap_faces[i + 1..] {
            let Some(profile) = seam_normal_profile(brep, a, b, &MID_SEAM_FRACTIONS) else {
                continue;
            };
            for (_f, angle_deg) in profile {
                let sin = angle_deg.to_radians().sin().abs();
                worst = Some(worst.map_or(sin, |w: f64| w.max(sin)));
            }
        }
    }
    worst
}

#[cfg(test)]
#[path = "cap_quality_postcondition_tests.rs"]
mod tests;
