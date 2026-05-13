//! KP20 degree-6 guide construction from available polar input.
//!
//! The full KP18 paper path starts from a complete cap-net. In the
//! current popkern polar pipeline the durable input available at this
//! layer is the bi-3 tensor-border plus a star point. This module
//! builds a total-degree-6 triangular guide sector per tensor-border
//! sector by degree-raising the outer bi-3 boundary and radially
//! blending it to the star point in triangular BB form. The builder is
//! deterministic, affine invariant, and side-by-side with KP15.

use super::bernstein_tri6::{TRI6_CP_COUNT, flat_index};
use super::kp_polar_boundary_input::TensorBorderBi3;
use super::kp_polar_shared::weighted_sum;
use super::kp20_bi6_guide_lsq::solve_kp20_guide_lsq;
use super::kp20_bi6_guide_sample::Kp20Tri6GuideSector;
use crate::geometry::Point3;

/// Input for the KP20 tensor-border guide builder.
#[derive(Debug, Clone)]
pub struct Kp20GuideBuildInput {
    /// Star point used at the extraordinary vertex.
    pub star_point: Point3,
    /// Bi-3 tensor-border around the polar hole.
    pub tensor_border: TensorBorderBi3,
}

/// Built KP20 guide ring.
#[derive(Debug, Clone, PartialEq)]
pub struct Kp20GuideRing {
    /// Valence.
    pub n: usize,
    /// One total-degree-6 triangular guide sector per polar sector.
    pub sectors: Vec<Kp20Tri6GuideSector>,
}

impl Kp20GuideRing {
    /// Access sector `k` modulo valence.
    pub fn sector(&self, k: usize) -> &Kp20Tri6GuideSector {
        &self.sectors[k % self.n]
    }
}

/// Build a KP20 guide ring from the tensor-border position rows.
pub fn build_kp20_guide_ring(input: &Kp20GuideBuildInput) -> Result<Kp20GuideRing, String> {
    let n = input.tensor_border.sectors.len();
    if n < 3 {
        return Err(format!("KP20 guide needs at least 3 sectors, got {n}"));
    }
    if let Ok(sectors) = solve_kp20_guide_lsq(&input.tensor_border, input.star_point) {
        return Ok(Kp20GuideRing { n, sectors });
    }
    let mut sectors = Vec::with_capacity(n);
    for (sector_index, sector) in input.tensor_border.sectors.iter().enumerate() {
        validate_position_row(sector_index, &sector[0])?;
        sectors.push(build_sector(input.star_point, sector[0]));
    }
    Ok(Kp20GuideRing { n, sectors })
}

fn build_sector(star: Point3, boundary_cubic: [Point3; 4]) -> Kp20Tri6GuideSector {
    let boundary = elevate_cubic_to_degree6(boundary_cubic);
    let mut cps = [Point3::new(0.0, 0.0, 0.0); TRI6_CP_COUNT];
    for i in 0..=6 {
        for j in 0..=(6 - i) {
            let level = i + j;
            let idx = flat_index(6, i, j);
            cps[idx] = if level == 0 {
                star
            } else {
                let radial = level as f64 / 6.0;
                let boundary_idx = ((6.0 * j as f64 / level as f64).round() as usize).min(6);
                weighted_sum(&[(1.0 - radial, star), (radial, boundary[boundary_idx])])
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
        let mut terms = Vec::with_capacity(j_max - j_min + 1);
        for j in j_min..=j_max {
            terms.push((C3[j] * binom3(i - j) / C6[i], row[j]));
        }
        out[i] = weighted_sum(&terms);
    }
    out
}

const fn binom3(i: usize) -> f64 {
    match i {
        0 | 3 => 1.0,
        1 | 2 => 3.0,
        _ => 0.0,
    }
}

fn validate_position_row(sector_index: usize, row: &[Point3; 4]) -> Result<(), String> {
    for (i, p) in row.iter().enumerate() {
        if !p.x.is_finite() || !p.y.is_finite() || !p.z.is_finite() {
            return Err(format!(
                "KP20 guide sector {sector_index} boundary CP {i} is non-finite"
            ));
        }
    }
    Ok(())
}

#[cfg(test)]
#[path = "kp20_bi6_guide_build_tests.rs"]
mod tests;
