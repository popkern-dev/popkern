//! Setback-border KP20 cap as the primary apex candidate.
//!
//! The legacy paper candidate derives its border from the two-stripe
//! sub-patch spokes (`TensorBorderBi3`). This entry derives it from
//! the FILLET SURFACES themselves via the classical corner-setback
//! boundary — n arcs on the G² sweeps meeting at trim-intersection
//! corners — which is the geometry the cap must actually bound in
//! the blended BRep. With the exact C²-at-m junction conditioning
//! (PK-093), the paper assembly closes every gate at machine zero on
//! this input class.

use super::super::kp_polar_bi6_promote::SELF_INT_TOL;
use super::super::kp20_paper_border_extract::kp20_paper_border_from_fillet_setback;
use super::super::kp20_polar_fold_summary::{
    gregory_patches_from_macro_patches, summarize_candidate_intersections,
};
use super::paper_primary_candidate::Kp20PaperPrimaryCandidate;
use crate::fillet_surface::FilletSurfaceResult;
use popkern_core::geometry::Point3;
use popkern_core::geometry::gregory::{Kp20StripEnforcement, assemble_kp20_paper_cap_conditioned};

/// Fraction of the way from the setback-loop centroid toward the
/// original apex vertex where the guide star point is placed. The
/// assembly's internal least squares makes the gates insensitive to
/// this choice (measured across 0.25..0.75 on the pyramid apex); it
/// only seeds the guide fit.
const STAR_LIFT_FRACTION: f64 = 0.25;

/// Build the paper cap from the fillet surfaces' setback boundary
/// and package it as a primary apex candidate.
pub(super) fn try_kp20_setback_primary_candidate(
    fillet_surfaces: &[FilletSurfaceResult],
    vertex_position: Point3,
) -> Result<Kp20PaperPrimaryCandidate, String> {
    let (border, extract_report) =
        kp20_paper_border_from_fillet_setback(fillet_surfaces, vertex_position)?;

    // Star point: centroid of the setback loop, lifted toward the
    // original vertex.
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
    let star = Point3::new(
        centroid.x + (vertex_position.x - centroid.x) * STAR_LIFT_FRACTION,
        centroid.y + (vertex_position.y - centroid.y) * STAR_LIFT_FRACTION,
        centroid.z + (vertex_position.z - centroid.z) * STAR_LIFT_FRACTION,
    );

    // Paper-faithful assembly (no guide fairing) with hard eq (8)/(9).
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

#[cfg(test)]
mod tests {
    use super::*;
    use popkern_core::tolerance::Tolerance;
    use popkern_core::topology::BRep;

    /// End-to-end production gate on real geometry: G² fillets on the
    /// pentagonal pyramid apex must yield an ACCEPTED setback paper
    /// candidate (all paper gates at relative 1e-6, fold-free).
    #[test]
    fn setback_candidate_accepted_on_pentagonal_pyramid_g2() {
        let tol = Tolerance::default();
        let pyramid = BRep::n_gonal_pyramid(5, 1.0, 1.2, tol.clone()).expect("pyramid");
        let apex = pyramid
            .vertices
            .iter()
            .max_by(|a, b| {
                a.1.position
                    .z
                    .partial_cmp(&b.1.position.z)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .map(|(h, v)| (h, v.position))
            .expect("apex");
        let apex_edges: Vec<_> = pyramid
            .edges
            .iter()
            .filter_map(|(eh, _)| {
                let he = pyramid.edges[eh].half_edge;
                let origin = pyramid.half_edges[he].origin;
                let twin_origin = pyramid.half_edges[pyramid.half_edges[he].twin].origin;
                (origin == apex.0 || twin_origin == apex.0).then_some(eh)
            })
            .collect();
        let specs: Vec<crate::FilletSpec> = apex_edges
            .iter()
            .map(|&edge| crate::FilletSpec {
                edge,
                radius: crate::FilletRadius::Constant(0.16),
            })
            .collect();
        let preview = crate::preview_fillet_with_continuity(
            &pyramid,
            &specs,
            tol,
            crate::ContinuityLevel::G2,
        )
        .expect("G2 preview");
        let fillet_surfaces: Vec<FilletSurfaceResult> = preview
            .per_edge
            .iter()
            .map(|eg| FilletSurfaceResult {
                surface: eg.blend_surface.clone(),
                trim_curve_1: eg.trim_curve_1.clone(),
                trim_curve_2: eg.trim_curve_2.clone(),
                continuity_downgrade: None,
            })
            .collect();
        let candidate = try_kp20_setback_primary_candidate(&fillet_surfaces, apex.1)
            .expect("setback candidate");
        eprintln!("setback candidate: {}", candidate.detail());
        assert!(
            candidate.is_acceptable(),
            "setback candidate rejected: {}",
            candidate.detail(),
        );
        assert_eq!(candidate.patches.len(), 20);
    }
}
