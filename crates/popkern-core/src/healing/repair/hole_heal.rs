//! End-to-end hole healing: detect → fill → commit (PK-199).
//!
//! The one entry point that actually closes a body. Per detected
//! hole the ladder runs single-patch G1 (4 sides, PK-224), then the
//! polar-wedge G1 fan (5+ sides, PK-248), then the G¹ PnS cap, the
//! Coons fan, and finally the always-constructible ruled fan; each
//! rung commits through its transactional operator. Every attempt —
//! committed or declined — is recorded with its reason in the
//! returned report; there is no silent fallback.

use crate::topology::BRep;

use super::hole_commit::commit_hole_cap;
use super::hole_commit_geometry::HoleCapGeometry;
use super::hole_fill::detect_holes;
use super::hole_polar_fan::polar_fan_g1_cap;
use super::hole_ruled_fan::ruled_fan_cap;
use super::hole_single_patch::commit_single_patch_hole_cap;
use super::n_sided_fill::fill_n_sided_hole;
use super::n_sided_fill_g1::fill_n_sided_hole_g1;

/// Which filler produced a committed cap.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HoleFillerKind {
    /// The single-patch degree-5 Gregory cap (4-sided holes only —
    /// no internal seams, ribbon-G1 to the neighbours, PK-224).
    SinglePatchG1,
    /// The polar-wedge degree-5 Gregory G1 fan (5+ sides, measured
    /// wall and cross-spoke G1 gates, PK-248).
    PolarFanG1,
    /// The G¹ PnS bi-4 cap.
    PnsG1,
    /// The Coons fan.
    CoonsFan,
    /// The ruled pole↔rim fan (G⁰ across spokes, exact containment).
    RuledFan,
}

/// Outcome of one hole's heal attempt.
#[derive(Debug, Clone)]
pub struct HoleHealAttempt {
    /// Sides of the detected hole.
    pub num_sides: usize,
    /// The filler whose cap committed, if any.
    pub committed_by: Option<HoleFillerKind>,
    /// Why the single-patch path did not commit (PK-224; `None` also
    /// when the hole is not 4-sided).
    pub single_patch_decline: Option<String>,
    /// Why the polar-fan path did not commit (PK-248; build or
    /// commit decline).
    pub polar_fan_decline: Option<String>,
    /// Why the G¹ PnS path did not commit (fill or commit decline).
    pub pns_decline: Option<String>,
    /// Why the Coons path did not commit (fill or commit decline).
    pub coons_decline: Option<String>,
    /// Why the ruled fan did not commit (build or commit decline).
    pub ruled_decline: Option<String>,
    /// Worst L2 disagreement of the committed cap's coherence gate.
    pub worst_l2_gap: Option<f64>,
}

/// Structured report of a [`heal_holes`] pass.
#[derive(Debug, Clone)]
pub struct HoleHealReport {
    /// Holes found by detection.
    pub holes_detected: usize,
    /// Holes whose cap committed.
    pub holes_committed: usize,
    /// Per-hole outcome, in detection order (smallest first).
    pub attempts: Vec<HoleHealAttempt>,
}

/// Detect every hole and commit a cap for each one that a filler can
/// serve, G¹ PnS first, Coons fan as the recorded fallback.
///
/// `tolerance` is the linear tolerance threaded through detection
/// matching, both fillers, and the commit gates — explicit, never
/// global state. Holes detected in one pass have disjoint edge sets,
/// so commits do not invalidate each other's rings; each commit is
/// individually transactional, and a hole whose fills both decline
/// leaves its boundary edges untouched.
pub fn heal_holes(brep: &mut BRep, tolerance: f64) -> HoleHealReport {
    let holes = detect_holes(brep);
    let mut attempts = Vec::with_capacity(holes.len());
    let mut committed = 0usize;

    for hole in &holes {
        let mut attempt = HoleHealAttempt {
            num_sides: hole.num_sides,
            committed_by: None,
            single_patch_decline: None,
            polar_fan_decline: None,
            pns_decline: None,
            coons_decline: None,
            ruled_decline: None,
            worst_l2_gap: None,
        };

        // Single-patch rung first (PK-224): its commit is fully
        // self-contained (own topology op, gates, transaction), so it
        // does not flow through the fan-shaped HoleCapGeometry path.
        match commit_single_patch_hole_cap(brep, hole, tolerance) {
            Ok(commit) => {
                attempt.committed_by = Some(HoleFillerKind::SinglePatchG1);
                attempt.worst_l2_gap = Some(commit.worst_l2_gap);
                committed += 1;
                attempts.push(attempt);
                continue;
            }
            Err(error) => attempt.single_patch_decline = Some(error.to_string()),
        }

        type CapBuilder<'a> = Box<dyn Fn(&BRep) -> Result<HoleCapGeometry, String> + 'a>;
        let stages: [(HoleFillerKind, CapBuilder); 4] = [
            (
                HoleFillerKind::PolarFanG1,
                Box::new(|brep: &BRep| {
                    polar_fan_g1_cap(brep, hole, tolerance).map_err(|error| error.to_string())
                }),
            ),
            (
                HoleFillerKind::PnsG1,
                Box::new(|brep: &BRep| {
                    let fill = fill_n_sided_hole_g1(brep, hole, tolerance)
                        .map_err(|error| error.to_string())?;
                    HoleCapGeometry::from_pns_fill(&fill, tolerance)
                        .map_err(|error| error.to_string())
                }),
            ),
            (
                HoleFillerKind::CoonsFan,
                Box::new(|brep: &BRep| {
                    let patches = fill_n_sided_hole(brep, hole, tolerance)
                        .ok_or_else(|| "Coons fan produced no patch set".to_string())?;
                    HoleCapGeometry::from_coons_patches(&patches, hole.num_sides, tolerance)
                        .map_err(|error| error.to_string())
                }),
            ),
            (
                HoleFillerKind::RuledFan,
                Box::new(|brep: &BRep| {
                    ruled_fan_cap(brep, hole, tolerance).map_err(|error| error.to_string())
                }),
            ),
        ];

        for (kind, build) in &stages {
            let outcome = build(brep).and_then(|geometry| {
                commit_hole_cap(brep, hole, &geometry, tolerance).map_err(|error| error.to_string())
            });
            match outcome {
                Ok(commit) => {
                    attempt.committed_by = Some(*kind);
                    attempt.worst_l2_gap = Some(commit.worst_l2_gap);
                    break;
                }
                Err(reason) => match kind {
                    // The single-patch rung never reaches this loop —
                    // it committed or recorded its decline above.
                    HoleFillerKind::SinglePatchG1 => {}
                    HoleFillerKind::PolarFanG1 => attempt.polar_fan_decline = Some(reason),
                    HoleFillerKind::PnsG1 => attempt.pns_decline = Some(reason),
                    HoleFillerKind::CoonsFan => attempt.coons_decline = Some(reason),
                    HoleFillerKind::RuledFan => attempt.ruled_decline = Some(reason),
                },
            }
        }

        if attempt.committed_by.is_some() {
            committed += 1;
        }
        attempts.push(attempt);
    }

    HoleHealReport {
        holes_detected: holes.len(),
        holes_committed: committed,
        attempts,
    }
}

#[cfg(test)]
#[path = "hole_heal_tests.rs"]
mod tests;
