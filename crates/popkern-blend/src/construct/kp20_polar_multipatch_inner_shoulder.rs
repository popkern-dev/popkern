//! Exact-worst origin-piece shaping for KP20 multipatch candidates.

// Index loops follow the coefficient-grid index algebra (PK-795).
#![allow(clippy::needless_range_loop)]

use popkern_core::geometry::Vec3;
use popkern_core::geometry::gregory::Kp20Bi6Patch;

use super::kp20_polar_fold_summary::Kp20FoldSite;
use super::kp20_polar_support_surface_split::evaluate_bi6_normal;

#[derive(Debug, Clone, Copy)]
pub(super) enum Kp20OriginLiftMode {
    CorePositive,
    CoreNegative,
    ShoulderPositive,
    ShoulderNegative,
}

pub(super) fn exact_worst_origin_surface(
    accepted: &Kp20Bi6Patch,
    mut surface: Kp20Bi6Patch,
    site: Kp20FoldSite,
    mode: Kp20OriginLiftMode,
    interior_scale: f64,
) -> Kp20Bi6Patch {
    let span = average_inner_span(&surface);
    let normal = evaluate_bi6_normal(accepted, site.u, site.v).unwrap_or(Vec3::new(0.0, 0.0, 1.0));
    let sign = match mode {
        Kp20OriginLiftMode::CoreNegative | Kp20OriginLiftMode::ShoulderNegative => -1.0,
        _ => 1.0,
    };
    match mode {
        Kp20OriginLiftMode::CorePositive | Kp20OriginLiftMode::CoreNegative => {
            move_core(&mut surface, site, normal, span * interior_scale * sign);
        }
        Kp20OriginLiftMode::ShoulderPositive | Kp20OriginLiftMode::ShoulderNegative => {
            move_corner_shoulder(&mut surface, site, normal, span * interior_scale * sign);
        }
    }
    surface
}

fn move_core(surface: &mut Kp20Bi6Patch, site: Kp20FoldSite, normal: Vec3, scale: f64) {
    for row in 2..=4 {
        for column in 2..=4 {
            let weight = local_weight(site.u, site.v, row as f64 / 6.0, column as f64 / 6.0, 0.45);
            surface.cps[row][column] = surface.cps[row][column] + normal * (scale * weight);
        }
    }
}

fn move_corner_shoulder(surface: &mut Kp20Bi6Patch, site: Kp20FoldSite, normal: Vec3, scale: f64) {
    for (row, column) in [(4, 4), (4, 5)] {
        let weight = local_weight(site.u, site.v, row as f64 / 6.0, column as f64 / 6.0, 0.55);
        surface.cps[row][column] = surface.cps[row][column] + normal * (scale * weight);
    }
}

fn local_weight(site_u: f64, site_v: f64, control_u: f64, control_v: f64, radius: f64) -> f64 {
    let du = control_u - site_u;
    let dv = control_v - site_v;
    let r2 = du * du + dv * dv;
    (1.0 - r2 / (radius * radius)).max(0.0).powi(2)
}

fn average_inner_span(surface: &Kp20Bi6Patch) -> f64 {
    let u_span = surface.cps[3][1].distance_to(&surface.cps[3][5]);
    let v_span = surface.cps[1][3].distance_to(&surface.cps[5][3]);
    ((u_span + v_span) * 0.5).max(1.0e-9)
}

#[cfg(test)]
mod tests {
    use popkern_core::geometry::Point3;
    use popkern_core::geometry::gregory::Kp20PieceId;

    use super::*;

    #[test]
    fn exact_worst_origin_lift_preserves_boundary_and_adjacent_rows() {
        let mut accepted = Kp20Bi6Patch::zero();
        let mut surface = Kp20Bi6Patch::zero();
        for row in 0..7 {
            for column in 0..7 {
                let point = Point3::new(row as f64, column as f64, 0.0);
                accepted.cps[row][column] = point;
                surface.cps[row][column] = point;
            }
        }
        let lifted = exact_worst_origin_surface(
            &accepted,
            surface,
            fold_site(),
            Kp20OriginLiftMode::CorePositive,
            1.0,
        );

        assert_eq!(lifted.cps[0][3], surface.cps[0][3]);
        assert_eq!(lifted.cps[1][3], surface.cps[1][3]);
        assert_eq!(lifted.cps[5][3], surface.cps[5][3]);
        assert_eq!(lifted.cps[6][3], surface.cps[6][3]);
        assert!(lifted.cps[3][3].z > surface.cps[3][3].z);
    }

    #[test]
    fn shoulder_lift_moves_corner_shoulder_without_adjacent_rows() {
        let mut accepted = Kp20Bi6Patch::zero();
        let mut surface = Kp20Bi6Patch::zero();
        for row in 0..7 {
            for column in 0..7 {
                let point = Point3::new(row as f64, column as f64, 0.0);
                accepted.cps[row][column] = point;
                surface.cps[row][column] = point;
            }
        }
        let lifted = exact_worst_origin_surface(
            &accepted,
            surface,
            fold_site(),
            Kp20OriginLiftMode::ShoulderPositive,
            1.0,
        );

        assert_eq!(lifted.cps[5][4], surface.cps[5][4]);
        assert!(lifted.cps[4][5].z > surface.cps[4][5].z);
    }

    fn fold_site() -> Kp20FoldSite {
        Kp20FoldSite {
            sector: 0,
            piece: Kp20PieceId::NearOriginLeft,
            u: 0.5,
            v: 0.5,
            min_jacobian_ratio: -0.1,
            min_jacobian_normal_cos: -0.1,
        }
    }
}
