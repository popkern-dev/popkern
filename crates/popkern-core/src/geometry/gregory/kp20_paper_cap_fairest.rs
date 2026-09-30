//! The fairest KP20 cap over candidate guide fairings (PK-613).
//!
//! The KP18 guide is a total-degree-6 least-squares fit to the border.
//! Where the border is not exactly representable (low, steep and
//! mixed-radius apexes), the plain fit spreads its misfit as interior
//! oscillation, and both the bi-6 sampling and the bi-9 resampling copy
//! it into the cap's star pieces. The recorded shape-conditioning
//! fairing ([`super::kp20_paper_guide::solve_kp20_paper_guide_with_fairing`])
//! trades guide misfit for a fair interior; the border bands are
//! overwritten from the input border afterwards, so the rims are
//! unaffected.
//!
//! No single weight serves every corner. Measured over the corner
//! census (159 KP20 fans, interior max principal curvature times the
//! smallest radius): weight 0.1 improved 56 fans (hexagon h1.2 variable
//! 4.5e5 → 3.96, octagon h0.4 7.5 → 2.3, pentagon h3 variable 226 → 17)
//! and worsened 20, mostly at valence 3 and 4, where the guide system is
//! singular and the null-space selector already steers it. So the cap is
//! assembled at each candidate weight and the fairest served interior is
//! kept; ties keep the paper-faithful weight 0.

use super::kp_degree_elevate::elevate_bi6_grid_to_bi9;
use super::kp20_bi9_fair::sample_basis;
use super::kp20_bi9_seams::Net;
use super::kp20_paper_cap::{
    Kp20PaperCapError, Kp20PaperCapReport, assemble_kp20_paper_cap_conditioned,
};
use super::kp20_paper_cap_strip::Kp20StripEnforcement;
use super::kp20_paper_guide::Kp20PaperGuideBorder;
use crate::geometry::Point3;

/// The guide-fairing weights tried, paper-faithful first. Measured on
/// the census: {0, 0.1} made 56 fans fairer and none worse, and lets 16
/// more corners serve a fan. Adding 0.03, 0.3 and 1 improved 15 more by
/// up to 3x but cost 30 % more time and let candidates whose committed
/// fan diverges from the assembly's reading serve (octagon h3 at up to
/// 2.4e3 x 1/r_min; the PK-612 ladder hands those to the fit).
pub const KP20_GUIDE_FAIRING_CANDIDATES: [f64; 2] = [0.0, 0.1];

/// The served cap's bi-9 nets: the resampling when it served, else the
/// bi-6 pieces elevated.
fn served_nets(report: &Kp20PaperCapReport) -> Vec<[Net; 4]> {
    match &report.bi9 {
        Ok(resample) => resample
            .pieces
            .chunks(4)
            .filter_map(|chunk| <[Net; 4]>::try_from(chunk).ok())
            .collect(),
        Err(_) => report
            .sectors
            .iter()
            .map(|m| core::array::from_fn(|k| elevate_bi6_grid_to_bi9(m.pieces[k].cps)))
            .collect(),
    }
}

/// Smallest per-piece min/mean `|S_u × S_v|` a candidate may reach
/// before it ranks as folded: the variational fit's fold floor (E6).
const RANK_MIN_REGULARITY: f64 = 2.0e-2;

/// The served cap's interior bend for ranking: the largest principal
/// curvature over interior samples of every piece, or infinity when a
/// piece nearly folds (per-piece regularity under
/// [`RANK_MIN_REGULARITY`]) or has a degenerate frame. The fold reading
/// matters: a pinched candidate reads small curvature beside its pinch
/// and zero at a degenerate sample, so ranked by curvature alone it can
/// outrank a fair cap.
pub fn kp20_cap_interior_curvature(report: &Kp20PaperCapReport) -> f64 {
    let table = sample_basis();
    let mut worst = 0.0_f64;
    for sector in served_nets(report) {
        for net in &sector {
            let (mut jac_min, mut jac_sum, mut count) = (f64::INFINITY, 0.0_f64, 0.0_f64);
            for bu in table {
                for bv in table {
                    let Some((kappa, jac)) = point_reading(net, bu, bv) else {
                        return f64::INFINITY;
                    };
                    worst = worst.max(kappa);
                    jac_min = jac_min.min(jac);
                    jac_sum += jac;
                    count += 1.0;
                }
            }
            if count > 0.0 && jac_min < RANK_MIN_REGULARITY * (jac_sum / count) {
                return f64::INFINITY;
            }
        }
    }
    worst
}

/// The largest principal curvature and `|S_u × S_v|` of a bi-9 net at
/// the parameter whose basis is `bu × bv`; `None` at a degenerate frame.
fn point_reading(net: &Net, bu: &[[f64; 3]; 10], bv: &[[f64; 3]; 10]) -> Option<(f64, f64)> {
    use crate::geometry::Vec3;
    let d = |p: usize, q: usize| {
        let mut acc = Vec3::new(0.0, 0.0, 0.0);
        for i in 0..10 {
            for j in 0..10 {
                let x = net[i][j];
                acc = acc + Vec3::new(x.x, x.y, x.z) * (bu[i][p] * bv[j][q]);
            }
        }
        acc
    };
    let (su, sv, suu, suv, svv) = (d(1, 0), d(0, 1), d(2, 0), d(1, 1), d(0, 2));
    let cross = su.cross(&sv);
    let jac = cross.length();
    let n = cross.normalized()?;
    let (e, f, g) = (su.dot(&su), su.dot(&sv), sv.dot(&sv));
    let (l, m, o) = (n.dot(&suu), n.dot(&suv), n.dot(&svv));
    let det = e * g - f * f;
    if det.abs() < 1e-30 {
        return None;
    }
    let (kg, h) = (
        (l * o - m * m) / det,
        (e * o + g * l - 2.0 * f * m) / (2.0 * det),
    );
    Some((h.abs() + (h * h - kg).max(0.0).sqrt(), jac))
}

/// Assemble the cap at every weight in `candidates`, fairest served
/// interior first; equal readings keep the candidates' order (the
/// paper-faithful weight first). A candidate that fails to assemble is
/// skipped; when every candidate fails, the first error is returned. An
/// empty list assembles the paper-faithful weight 0.
///
/// The caller commits the head and falls back along the list: the
/// fairest interior is not always the cap whose later commit
/// postconditions pass (pyramid-3 h3 wide at G1: the faired fan's seams
/// read g1 1.5e-1 where the paper fan's pass).
pub fn assemble_kp20_paper_cap_ranked(
    border: &Kp20PaperGuideBorder,
    star_point: Point3,
    enforcement: Kp20StripEnforcement,
    candidates: &[f64],
) -> Result<Vec<Kp20PaperCapReport>, Kp20PaperCapError> {
    if candidates.is_empty() {
        return assemble_kp20_paper_cap_conditioned(border, star_point, enforcement, 0.0)
            .map(|report| vec![report]);
    }
    let mut ranked: Vec<(f64, Kp20PaperCapReport)> = Vec::new();
    let mut first_error = None;
    for &weight in candidates {
        match assemble_kp20_paper_cap_conditioned(border, star_point, enforcement, weight) {
            Ok(report) => ranked.push((kp20_cap_interior_curvature(&report), report)),
            Err(e) => {
                first_error.get_or_insert(e);
            }
        }
    }
    // Stable: equal readings keep the candidates' order.
    ranked.sort_by(|a, b| a.0.total_cmp(&b.0));
    match (ranked.is_empty(), first_error) {
        (true, Some(e)) => Err(e),
        _ => Ok(ranked.into_iter().map(|(_, report)| report).collect()),
    }
}
