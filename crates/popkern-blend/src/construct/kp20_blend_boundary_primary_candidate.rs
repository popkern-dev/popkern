//! Blend-boundary KP20 paper cap as a primary apex candidate
//! (PK-090 step 2).
//!
//! The paper cap natively accepts two-piece-per-sector borders. The
//! collapsed n-arc adapters (midpoint-split tensor border, setback
//! extractor) both discard the honest sector corner: the apex blend
//! boundary's 2n alternating segments — one fillet-end arc and one
//! face-edge bridge per sector — ARE the two genuine surface pieces
//! meeting at `c` with their true corner tangents. This candidate
//! feeds that boundary to the literal extractor
//! (`kp20_paper_border_from_blend_boundary`, whose macro-edge refit
//! preserves the true corner tangents at `c`) and assembles with the
//! same hard-G² conditioning as the setback candidate. Acceptance is
//! the shared paper-candidate gate, including the
//! `input_corner_mismatch` corner gate.

use super::super::kp_polar_bi6_promote::SELF_INT_TOL;
use super::super::kp20_paper_border_extract::kp20_paper_border_from_blend_boundary;
use super::super::kp20_polar_fold_summary::{
    gregory_patches_from_macro_patches, summarize_candidate_intersections,
};
use super::paper_primary_candidate::Kp20PaperPrimaryCandidate;
use crate::blend::BlendBoundary;
use crate::fillet_surface::FilletSurfaceResult;
use popkern_core::geometry::Point3;
use popkern_core::geometry::gregory::{Kp20StripEnforcement, assemble_kp20_paper_cap_conditioned};
use popkern_core::topology::BRep;

/// Fraction of the way from the boundary-loop centroid toward the
/// original apex vertex where the guide star point is placed. Same
/// insensitive seed as the setback candidate (the assembly's internal
/// least squares makes the gates indifferent across 0.25..0.75).
const STAR_LIFT_FRACTION: f64 = 0.25;

/// Build the paper cap from the honest 2n-segment blend boundary and
/// package it as a primary apex candidate.
pub(super) fn try_kp20_blend_boundary_primary_candidate(
    boundary: &BlendBoundary,
    fillet_surfaces: &[FilletSurfaceResult],
    brep: &BRep,
) -> Result<Kp20PaperPrimaryCandidate, String> {
    let (border, extract_report) =
        kp20_paper_border_from_blend_boundary(boundary, fillet_surfaces, brep)?;

    // Star point: centroid of the boundary rim, lifted toward the
    // original vertex (mirrors the setback candidate's seed).
    let centroid = {
        let mut acc = [0.0_f64; 3];
        let mut count = 0.0;
        for sector in &border.sectors {
            for piece in sector {
                for cp in &piece.rows[0] {
                    acc[0] += cp.x;
                    acc[1] += cp.y;
                    acc[2] += cp.z;
                    count += 1.0;
                }
            }
        }
        Point3::new(acc[0] / count, acc[1] / count, acc[2] / count)
    };
    let vertex = boundary.vertex_position;
    let star = Point3::new(
        centroid.x + (vertex.x - centroid.x) * STAR_LIFT_FRACTION,
        centroid.y + (vertex.y - centroid.y) * STAR_LIFT_FRACTION,
        centroid.z + (vertex.z - centroid.z) * STAR_LIFT_FRACTION,
    );

    // Paper-faithful assembly (no guide fairing) with hard eq (8)/(9)
    // — the same conditioning the setback candidate uses.
    let report =
        assemble_kp20_paper_cap_conditioned(&border, star, Kp20StripEnforcement::HardG2, 0.0)
            .map_err(|e| e.to_string())?;

    let macro_patches = report.sectors.clone();
    let patches = gregory_patches_from_macro_patches(&macro_patches);
    let intersections = patches
        .iter()
        .enumerate()
        .filter_map(|(idx, gp)| {
            let self_int = gp.self_intersection_report(SELF_INT_TOL);
            self_int.intersects().then_some((idx, self_int))
        })
        .collect::<Vec<_>>();
    let summary = summarize_candidate_intersections(&patches, &intersections);
    let scale = border
        .sectors
        .iter()
        .flat_map(|sector| sector.iter())
        .flat_map(|piece| piece.rows.iter().flatten())
        .map(|cp| cp.distance_to(&star))
        .fold(0.0_f64, f64::max)
        .max(extract_report.max_loop_gap);

    Ok(Kp20PaperPrimaryCandidate {
        patches,
        summary,
        report,
        scale,
    })
}
