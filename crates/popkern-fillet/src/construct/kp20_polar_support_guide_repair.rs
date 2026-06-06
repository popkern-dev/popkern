//! Bounded KP20 guide-construction repair candidates.

use popkern_core::geometry::gregory::{
    Kp20Bi6MacroCapInput, Kp20Bi6ValidationReport, Kp20GuideRing, Kp20Tri6GuideSector,
    TRI6_CP_COUNT, TensorBorderBi3, assemble_kp20_bi6_macro_cap,
    kp20_positive_nonoverlap_c2_sector_piece_maps, sample_kp20_macro_patch_from_piece_maps,
    solve_kp20_guide_lsq_with_prior, tri6_flat_index,
};
use popkern_core::geometry::{Point3, Vec3};

use super::super::super::kp20_polar_fold_summary::{
    Kp20FoldSite, Kp20SamplerFoldSummary, candidate_summary_from_cap, validation_is_clean,
};
use super::super::super::kp20_polar_support_reparameterization::{
    Kp20GuideQualityReport, Kp20GuideRepairReport,
};
use super::guide_quality::measure_kp20_guide_quality;

pub(super) fn transactional_kp20_guide_repair(
    guide: &Kp20GuideRing,
    tensor_border: &TensorBorderBi3,
    star_point: Point3,
    baseline: &(Kp20SamplerFoldSummary, Kp20Bi6ValidationReport),
) -> Result<(Kp20GuideRing, Kp20GuideQualityReport, Kp20GuideRepairReport), String> {
    let start = measure_kp20_guide_quality(guide);
    let mut report = Kp20GuideRepairReport {
        evaluated: 0,
        accepted: 0,
        rejected: 0,
        start_reversals: start.reversals,
        end_reversals: start.reversals,
        blend: 0.0,
        decline_reason: "guide_repair_no_candidate",
    };
    let radial = radial_guide_ring(tensor_border, star_point);
    let mut best_guide = guide.clone();
    let mut best_quality = start.clone();
    let mut downstream_rejected = false;
    for weight in [1.0e-6, 1.0e-5, 1.0e-4, 1.0e-3, 1.0e-2, 1.0e-1] {
        if evaluate_guide_prior_candidate(
            guide,
            tensor_border,
            star_point,
            &radial,
            weight,
            true,
            &mut best_guide,
            &mut best_quality,
            &mut report,
            &mut downstream_rejected,
        )? {
            break;
        }
    }
    if report.accepted == 0 {
        for scale in [-0.10, 0.10, -0.25, 0.25] {
            let fold_prior = fold_site_guide_ring(
                guide,
                tensor_border,
                star_point,
                &baseline.0.fold_sites,
                scale,
            );
            if evaluate_guide_prior_candidate(
                guide,
                tensor_border,
                star_point,
                &fold_prior,
                1.0e-3,
                false,
                &mut best_guide,
                &mut best_quality,
                &mut report,
                &mut downstream_rejected,
            )? {
                break;
            }
        }
    }
    report.rejected = report.evaluated - report.accepted;
    report.end_reversals = best_quality.reversals;
    if report.accepted == 0 && downstream_rejected {
        report.decline_reason = "guide_repair_downstream_fold_remaining";
    }
    Ok((best_guide, best_quality, report))
}

fn evaluate_guide_prior_candidate(
    guide: &Kp20GuideRing,
    tensor_border: &TensorBorderBi3,
    star_point: Point3,
    prior: &Kp20GuideRing,
    weight: f64,
    require_orientation_improvement: bool,
    best_guide: &mut Kp20GuideRing,
    best_quality: &mut Kp20GuideQualityReport,
    report: &mut Kp20GuideRepairReport,
    downstream_rejected: &mut bool,
) -> Result<bool, String> {
    report.evaluated += 1;
    let candidate =
        solve_kp20_guide_lsq_with_prior(tensor_border, star_point, &prior.sectors, weight)
            .map(|sectors| Kp20GuideRing {
                n: guide.n,
                sectors,
            })
            .map_err(|e| e.to_string())?;
    let quality = measure_kp20_guide_quality(&candidate);
    if require_orientation_improvement && !is_better_guide_quality(&quality, best_quality) {
        report.rejected += 1;
        report.decline_reason = "guide_repair_no_orientation_improvement";
        return Ok(false);
    }
    if !require_orientation_improvement && guide_quality_regressed(&quality, best_quality) {
        report.rejected += 1;
        report.decline_reason = "guide_repair_orientation_regressed";
        return Ok(false);
    }
    let summary = candidate_summary_for_guide(tensor_border, &candidate)?;
    if !validation_is_clean(summary.1) {
        report.rejected += 1;
        report.decline_reason = "guide_repair_downstream_validation_not_clean";
        return Ok(false);
    }
    if summary.0.fold_count > 0 {
        report.rejected += 1;
        *downstream_rejected = true;
        report.decline_reason = "guide_repair_downstream_fold_remaining";
        return Ok(false);
    }
    *best_guide = candidate;
    *best_quality = quality;
    report.accepted = 1;
    report.blend = weight;
    report.decline_reason = "none";
    Ok(true)
}

fn is_better_guide_quality(
    candidate: &Kp20GuideQualityReport,
    best: &Kp20GuideQualityReport,
) -> bool {
    candidate.reversals < best.reversals
        || (candidate.reversals == best.reversals && candidate.min_normal_cos > best.min_normal_cos)
}

fn guide_quality_regressed(
    candidate: &Kp20GuideQualityReport,
    best: &Kp20GuideQualityReport,
) -> bool {
    candidate.reversals > best.reversals
        || (candidate.reversals == best.reversals && candidate.min_normal_cos < best.min_normal_cos)
}

fn candidate_summary_for_guide(
    tensor_border: &TensorBorderBi3,
    guide: &Kp20GuideRing,
) -> Result<(Kp20SamplerFoldSummary, Kp20Bi6ValidationReport), String> {
    let sampled = guide
        .sectors
        .iter()
        .map(|sector| {
            sample_kp20_macro_patch_from_piece_maps(
                sector,
                kp20_positive_nonoverlap_c2_sector_piece_maps(),
            )
        })
        .collect::<Vec<_>>();
    let cap = assemble_kp20_bi6_macro_cap(&Kp20Bi6MacroCapInput {
        tensor_border: tensor_border.clone(),
        sampled_macro_patches: sampled,
    })
    .map_err(|e| e.to_string())?;
    candidate_summary_from_cap(tensor_border, cap.macro_patches)
}

fn radial_guide_ring(tensor_border: &TensorBorderBi3, star_point: Point3) -> Kp20GuideRing {
    Kp20GuideRing {
        n: tensor_border.sectors.len(),
        sectors: tensor_border
            .sectors
            .iter()
            .map(|sector| radial_sector(star_point, sector[0]))
            .collect(),
    }
}

fn fold_site_guide_ring(
    guide: &Kp20GuideRing,
    tensor_border: &TensorBorderBi3,
    star_point: Point3,
    fold_sites: &[Kp20FoldSite],
    scale: f64,
) -> Kp20GuideRing {
    let mut out = guide.clone();
    for site in fold_sites {
        let Some([u, v]) = guide_uv_for_fold_site(*site) else {
            continue;
        };
        if u < 0.0 || v < 0.0 || u + v > 1.0 {
            continue;
        }
        let direction = sector_normal(tensor_border, star_point, site.sector);
        let span = sector_span(tensor_border, star_point, site.sector);
        let severity = (1.0 - site.min_jacobian_normal_cos).clamp(0.0, 2.0);
        apply_fold_site_guide_lift(
            &mut out.sectors[site.sector],
            [u, v],
            direction * (scale * span * severity),
        );
    }
    out
}

fn guide_uv_for_fold_site(site: Kp20FoldSite) -> Option<[f64; 2]> {
    let maps = kp20_positive_nonoverlap_c2_sector_piece_maps();
    let (_, corner_maps) = maps.into_iter().find(|(piece, _)| *piece == site.piece)?;
    let map = corner_maps[0].1;
    Some([
        map.origin[0] + site.u * map.du[0] + site.v * map.dv[0],
        map.origin[1] + site.u * map.du[1] + site.v * map.dv[1],
    ])
}

fn sector_normal(tensor_border: &TensorBorderBi3, star_point: Point3, sector: usize) -> Vec3 {
    let row = tensor_border.sectors[sector][0];
    let spoke_a = row[0] - star_point;
    let spoke_b = row[3] - star_point;
    spoke_a
        .cross(&spoke_b)
        .normalized()
        .unwrap_or(Vec3::new(0.0, 0.0, 1.0))
}

fn sector_span(tensor_border: &TensorBorderBi3, star_point: Point3, sector: usize) -> f64 {
    let row = tensor_border.sectors[sector][0];
    let spoke_a = row[0] - star_point;
    let spoke_b = row[3] - star_point;
    0.5 * (spoke_a.length() + spoke_b.length())
}

fn apply_fold_site_guide_lift(sector: &mut Kp20Tri6GuideSector, uv: [f64; 2], lift: Vec3) {
    for i in 0..=6 {
        for j in 0..=(6 - i) {
            let level = i + j;
            if !(2..6).contains(&level) {
                continue;
            }
            let cp_u = i as f64 / 6.0;
            let cp_v = j as f64 / 6.0;
            let du = cp_u - uv[0];
            let dv = cp_v - uv[1];
            let radius2 = 0.22_f64 * 0.22_f64;
            let locality = (1.0 - (du * du + dv * dv) / radius2).max(0.0).powi(2);
            if locality == 0.0 {
                continue;
            }
            let idx = tri6_flat_index(6, i, j);
            sector.cps[idx] = sector.cps[idx] + lift * locality;
        }
    }
}

fn radial_sector(star_point: Point3, boundary_cubic: [Point3; 4]) -> Kp20Tri6GuideSector {
    let boundary = elevate_cubic_to_degree6(boundary_cubic);
    let mut cps = [Point3::new(0.0, 0.0, 0.0); TRI6_CP_COUNT];
    for i in 0..=6 {
        for j in 0..=(6 - i) {
            let level = i + j;
            let idx = tri6_flat_index(6, i, j);
            cps[idx] = if level == 0 {
                star_point
            } else {
                let radial = level as f64 / 6.0;
                let boundary_idx = ((6.0 * j as f64 / level as f64).round() as usize).min(6);
                lerp(star_point, boundary[boundary_idx], radial)
            };
        }
    }
    Kp20Tri6GuideSector { cps }
}

fn elevate_cubic_to_degree6(row: [Point3; 4]) -> [Point3; 7] {
    const C3: [f64; 4] = [1.0, 3.0, 3.0, 1.0];
    const C6: [f64; 7] = [1.0, 6.0, 15.0, 20.0, 15.0, 6.0, 1.0];
    let mut out = [Point3::new(0.0, 0.0, 0.0); 7];
    for i in 0_usize..7 {
        let j_min = i.saturating_sub(3);
        let j_max = i.min(3);
        let mut sum = Point3::new(0.0, 0.0, 0.0);
        for j in j_min..=j_max {
            let weight = C3[j] * binom3(i - j) / C6[i];
            sum.x += weight * row[j].x;
            sum.y += weight * row[j].y;
            sum.z += weight * row[j].z;
        }
        out[i] = sum;
    }
    out
}

fn lerp(a: Point3, b: Point3, t: f64) -> Point3 {
    Point3::new(
        (1.0 - t) * a.x + t * b.x,
        (1.0 - t) * a.y + t * b.y,
        (1.0 - t) * a.z + t * b.z,
    )
}

const fn binom3(i: usize) -> f64 {
    match i {
        0 | 3 => 1.0,
        1 | 2 => 3.0,
        _ => 0.0,
    }
}
