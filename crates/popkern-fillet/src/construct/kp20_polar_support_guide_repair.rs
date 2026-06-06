//! Bounded KP20 guide-construction repair candidates.

use popkern_core::geometry::Point3;
use popkern_core::geometry::gregory::{
    Kp20Bi6MacroCapInput, Kp20Bi6ValidationReport, Kp20GuideRing, Kp20Tri6GuideSector,
    TRI6_CP_COUNT, TensorBorderBi3, assemble_kp20_bi6_macro_cap,
    kp20_nonoverlap_coupled_sector_piece_maps, sample_kp20_macro_patch_from_piece_maps,
    tri6_flat_index,
};

use super::super::super::kp20_polar_fold_summary::{
    Kp20SamplerFoldSummary, candidate_summary_from_cap, is_better_summary, validation_is_clean,
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
    let (candidate_guide, candidate_quality, mut report) =
        repair_kp20_guide_orientation(guide, tensor_border, star_point, &start);
    if report.accepted == 0 {
        return Ok((guide.clone(), start, report));
    }
    let candidate_sampled = candidate_guide
        .sectors
        .iter()
        .map(|sector| {
            sample_kp20_macro_patch_from_piece_maps(
                sector,
                kp20_nonoverlap_coupled_sector_piece_maps(),
            )
        })
        .collect::<Vec<_>>();
    let candidate_cap = assemble_kp20_bi6_macro_cap(&Kp20Bi6MacroCapInput {
        tensor_border: tensor_border.clone(),
        sampled_macro_patches: candidate_sampled,
    })
    .map_err(|e| e.to_string())?;
    let candidate_summary = candidate_summary_from_cap(tensor_border, candidate_cap.macro_patches)?;
    if validation_is_clean(candidate_summary.1)
        && is_better_summary(&candidate_summary, Some(baseline))
    {
        return Ok((candidate_guide, candidate_quality, report));
    }
    report.accepted = 0;
    report.rejected = report.evaluated;
    report.end_reversals = report.start_reversals;
    report.blend = 0.0;
    report.decline_reason = "guide_repair_downstream_no_fold_reduction";
    Ok((guide.clone(), start, report))
}

fn repair_kp20_guide_orientation(
    guide: &Kp20GuideRing,
    tensor_border: &TensorBorderBi3,
    star_point: Point3,
    start: &Kp20GuideQualityReport,
) -> (Kp20GuideRing, Kp20GuideQualityReport, Kp20GuideRepairReport) {
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
    for blend in [0.25, 0.50, 0.75, 1.00] {
        report.evaluated += 1;
        let candidate = blend_guides(guide, &radial, blend);
        let quality = measure_kp20_guide_quality(&candidate);
        if !is_better_guide_quality(&quality, &best_quality) {
            report.rejected += 1;
            report.decline_reason = "guide_repair_no_orientation_improvement";
            continue;
        }
        best_guide = candidate;
        best_quality = quality;
        report.accepted = 1;
        report.blend = blend;
        report.decline_reason = "none";
    }
    report.rejected = report.evaluated - report.accepted;
    report.end_reversals = best_quality.reversals;
    (best_guide, best_quality, report)
}

fn is_better_guide_quality(
    candidate: &Kp20GuideQualityReport,
    best: &Kp20GuideQualityReport,
) -> bool {
    candidate.reversals < best.reversals
        || (candidate.reversals == best.reversals && candidate.min_normal_cos > best.min_normal_cos)
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

fn blend_guides(base: &Kp20GuideRing, target: &Kp20GuideRing, blend: f64) -> Kp20GuideRing {
    Kp20GuideRing {
        n: base.n,
        sectors: base
            .sectors
            .iter()
            .zip(&target.sectors)
            .map(|(base_sector, target_sector)| Kp20Tri6GuideSector {
                cps: core::array::from_fn(|idx| {
                    lerp(base_sector.cps[idx], target_sector.cps[idx], blend)
                }),
            })
            .collect(),
    }
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
