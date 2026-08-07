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

//!
//! # Commit-time vs serving-time (PK-220)
//!
//! The commit-time number is measured MID-PIPELINE, through pcurves
//! that are not yet final (the fillet pipeline recomputes every
//! pcurve at its tail). Measured on the welded pyramid fixtures, the
//! commit-time seam resolution mis-parameterizes the pole half of
//! each spoke: the profile reads clean for fractions < 0.5, exactly
//! zero at 0.5, and climbs to tens of degrees above 0.5 — while the
//! served body, with identical control points (checksum-equal), is
//! clean along the whole seam. The recorded 2.8e-1 (pyramid3),
//! 2.9e-1 (pyramid4) and 5.9e-1 (oblique N=3) were all exactly the
//! 0.6-fraction artifact, not geometry. The serving-boundary refresh
//! below re-measures on the final body and appends
//! `served_mid_seam_normal` to the diagnostic; gates pin the served
//! value.

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

/// Serving-boundary mid-seam refresh (PK-220): re-measure every
/// committed multi-sector cap on the FINAL body and append the
/// served value to the commit-time diagnostic (both the route report
/// row and the top-level apex diagnostic copy). Runs at the fillet
/// pipeline tail, after every pcurve author, because the commit-time
/// number is measured through not-yet-final pcurves and over-reads
/// the pole half of each spoke (see module doc). Caps whose faces
/// did not survive to serving (absorbed or re-partitioned) keep the
/// commit-time value — a missing re-measurement must stay visible,
/// not silently pass. Returns the number of diagnostics refreshed.
pub(crate) fn refresh_served_mid_seam_diagnostics(
    brep: &BRep,
    repair_report: &mut crate::types::BlendRepairReport,
) -> usize {
    const CODE: &str = "exact_boundary_cap.mid_seam_normal";
    let mut refreshed = 0;
    for route in &mut repair_report.apex_routes {
        let faces = &route.exact_boundary_report.committed_faces;
        if faces.len() < 2 {
            continue;
        }
        let Some(served) = worst_mid_seam_misalignment(brep, faces) else {
            continue;
        };
        let served_field = format!("served_mid_seam_normal={served:.3e}");
        // The route-report row is the authoritative surface; capture
        // its commit-time detail so the top-level apex-diagnostic
        // copy can be matched EXACTLY (several routes can commit the
        // same face count).
        let mut commit_details: Vec<String> = Vec::new();
        for diagnostic in &mut route.exact_boundary_report.diagnostics {
            if diagnostic.code != CODE || diagnostic.stage != CODE {
                continue;
            }
            match diagnostic.detail_code.as_mut() {
                Some(detail) if !detail.contains("served_mid_seam_normal=") => {
                    commit_details.push(detail.clone());
                    detail.push(';');
                    detail.push_str(&served_field);
                }
                Some(_) => {}
                None => diagnostic.detail_code = Some(served_field.clone()),
            }
            diagnostic.residual = Some(served);
            refreshed += 1;
        }
        for diagnostic in &mut repair_report.apex_diagnostics {
            let crate::types::ApexConstructionDiagnostic::ExactBoundaryCap(d) = diagnostic else {
                continue;
            };
            if d.code != CODE {
                continue;
            }
            let matched = d
                .detail_code
                .as_deref()
                .is_some_and(|detail| commit_details.iter().any(|c| c == detail));
            if matched && let Some(detail) = d.detail_code.as_mut() {
                detail.push(';');
                detail.push_str(&served_field);
                d.residual = Some(served);
            }
        }
    }
    refreshed
}

#[cfg(test)]
#[path = "cap_quality_postcondition_tests.rs"]
mod tests;

/// PK-250 serving-boundary check: a direct-apex corner's region — the
/// built cap faces plus every edge-adjacent neighbour — must be free
/// of face-face self-intersection on the SERVED body. No construction
/// gate can see one (measured on the mixed [0.15, 0.1, 0.1] cube
/// corner: the candidate committed cleanly and served a real
/// strip×cap-sector contact that only `validate_all` caught), and a
/// mid-pipeline check is premature by the PK-165/PK-220 lessons —
/// wired as an acceptance gate inside the candidate facade it flipped
/// 20 serving fixtures whose final bodies are clean. The fillet
/// pipeline therefore measures this AFTER its tail passes and retries
/// with the vertex excluded from promotion when it fires. Returns the
/// offending face pair as a display string so the retry record can
/// name it.
pub(crate) fn validate_candidate_region_self_intersection(
    brep: &BRep,
    committed_faces: &[FaceHandle],
) -> Result<(), String> {
    let mut region: Vec<FaceHandle> = committed_faces.to_vec();
    for &face in committed_faces {
        let Some(loop_start) = brep
            .faces
            .get(face)
            .and_then(|f| brep.loops.get(f.outer_loop))
            .map(|l| l.half_edge)
        else {
            continue;
        };
        let mut cursor = loop_start;
        while let Some(half_edge) = brep.half_edges.get(cursor) {
            if let Some(neighbour) = brep
                .half_edges
                .get(half_edge.twin)
                .and_then(|twin| brep.loops.get(twin.loop_))
                .map(|l| l.face)
                && !region.contains(&neighbour)
            {
                region.push(neighbour);
            }
            cursor = half_edge.next;
            if cursor == loop_start {
                break;
            }
        }
    }
    brep.validate_self_intersection_among(&region)
        .map_err(|error| format!("{error}"))
}
