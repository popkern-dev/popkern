//! Structured guide-quality diagnostics for KP20 support construction.

use popkern_core::geometry::Vec3;
use popkern_core::geometry::gregory::Kp20GuideRing;

use super::super::super::kp20_polar_support_reparameterization::Kp20GuideQualityReport;

pub(super) fn measure_kp20_guide_quality(guide: &Kp20GuideRing) -> Kp20GuideQualityReport {
    let mut report = Kp20GuideQualityReport {
        samples: 0,
        reversals: 0,
        min_normal_cos: f64::INFINITY,
        worst_sector: 0,
        worst_u: 0.0,
        worst_v: 0.0,
        classification: "guide_orientation_clean",
    };
    for (sector_index, sector) in guide.sectors.iter().enumerate() {
        let reference = guide_normal(sector, 1.0 / 3.0, 1.0 / 3.0).unwrap_or(Vec3::zero());
        for [u, v] in GUIDE_SAMPLES {
            let Some(normal) = guide_normal(sector, u, v) else {
                report.reversals += 1;
                report.samples += 1;
                report.classification = "guide_degenerate_normal";
                continue;
            };
            let cos = normal.dot(&reference).clamp(-1.0, 1.0);
            report.samples += 1;
            if cos < 0.0 {
                report.reversals += 1;
            }
            if cos < report.min_normal_cos {
                report.min_normal_cos = cos;
                report.worst_sector = sector_index;
                report.worst_u = u;
                report.worst_v = v;
            }
        }
    }
    if report.reversals > 0 && report.classification == "guide_orientation_clean" {
        report.classification = "guide_orientation_reversal";
    }
    report
}

const GUIDE_SAMPLES: [[f64; 2]; 15] = [
    [0.125, 0.125],
    [0.250, 0.125],
    [0.375, 0.125],
    [0.500, 0.125],
    [0.625, 0.125],
    [0.125, 0.250],
    [0.250, 0.250],
    [0.375, 0.250],
    [0.500, 0.250],
    [0.125, 0.375],
    [0.250, 0.375],
    [0.375, 0.375],
    [0.125, 0.500],
    [0.250, 0.500],
    [0.125, 0.625],
];

fn guide_normal(
    sector: &popkern_core::geometry::gregory::Kp20Tri6GuideSector,
    u: f64,
    v: f64,
) -> Option<Vec3> {
    let jet = sector.jet3(u, v);
    let du = Vec3::new(jet.h[1][0].x, jet.h[1][0].y, jet.h[1][0].z);
    let dv = Vec3::new(jet.h[0][1].x, jet.h[0][1].y, jet.h[0][1].z);
    du.cross(&dv).normalized()
}
