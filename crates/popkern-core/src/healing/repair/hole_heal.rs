//! End-to-end hole healing: detect → fill → commit (PK-199).
//!
//! The one entry point that actually closes a body. Per detected
//! hole it tries the G¹ PnS cap first (tangent-continuous, valence
//! 3–12), falls back to the Coons fan, and commits whichever fill
//! survives through [`super::hole_commit::commit_hole_cap`]. Every
//! attempt — committed or declined — is recorded with its reason in
//! the returned report; there is no silent fallback.

use crate::topology::BRep;

use super::hole_commit::commit_hole_cap;
use super::hole_commit_geometry::HoleCapGeometry;
use super::hole_fill::detect_holes;
use super::hole_ruled_fan::ruled_fan_cap;
use super::n_sided_fill::fill_n_sided_hole;
use super::n_sided_fill_g1::fill_n_sided_hole_g1;

/// Which filler produced a committed cap.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HoleFillerKind {
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
            pns_decline: None,
            coons_decline: None,
            ruled_decline: None,
            worst_l2_gap: None,
        };

        type CapBuilder<'a> = Box<dyn Fn(&BRep) -> Result<HoleCapGeometry, String> + 'a>;
        let stages: [(HoleFillerKind, CapBuilder); 3] = [
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
