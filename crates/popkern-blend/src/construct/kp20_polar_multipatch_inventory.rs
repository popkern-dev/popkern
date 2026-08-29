//! Candidate-only KP20 support-strip/inner-pole ownership inventory.

use popkern_core::geometry::gregory::Kp20PieceId;

use super::kp20_polar_fold_summary::Kp20SamplerFoldSummary;
use super::kp20_polar_support_surface_split::support_owner_piece;

#[derive(Debug, Clone)]
pub(super) struct Kp20MultiPatchInventoryReport {
    pub(super) attempted_sector_count: usize,
    pub(super) attempted_strip_patch_count: usize,
    pub(super) attempted_inner_patch_count: usize,
    pub(super) attempted_interface_count: usize,
    pub(super) strip_control_point_count: usize,
    pub(super) inner_control_point_count: usize,
    pub(super) inventory_mismatch: Option<&'static str>,
    pub(super) sector_reports: Vec<Kp20MultiPatchSectorInventory>,
}

#[derive(Debug, Clone)]
pub(super) struct Kp20MultiPatchSectorInventory {
    pub(super) sector_index: usize,
    pub(super) sector_name: String,
    pub(super) strip_sides: Vec<Kp20PieceId>,
    pub(super) strip_names: Vec<String>,
    pub(super) inner_patch_name: String,
    pub(super) interface_names: Vec<String>,
}

pub(super) fn multipatch_inventory_from_fold_summary(
    summary: &Kp20SamplerFoldSummary,
) -> Kp20MultiPatchInventoryReport {
    let mut sector_sides = Vec::<(usize, Kp20PieceId)>::new();
    for site in &summary.fold_sites {
        let side = support_owner_piece(site.piece);
        if !sector_sides.contains(&(site.sector, side)) {
            sector_sides.push((site.sector, side));
        }
    }
    sector_sides.sort_by_key(|(sector, side)| (*sector, side.index()));

    let mut sector_reports = Vec::new();
    for sector in unique_sectors(&sector_sides) {
        let sides = sector_sides
            .iter()
            .filter_map(|(side_sector, side)| (*side_sector == sector).then_some(*side))
            .collect::<Vec<_>>();
        let sector_name = format!("kp20_polar_macro.sector_{sector}");
        let inner_patch_name = format!("{sector_name}.inner_pole_patch");
        let strip_names = sides
            .iter()
            .map(|side| format!("{sector_name}.{}_support_strip", side_label(*side)))
            .collect::<Vec<_>>();
        let interface_names = sides
            .iter()
            .map(|side| format!("{sector_name}.{}_support_interface", side_label(*side)))
            .collect::<Vec<_>>();
        sector_reports.push(Kp20MultiPatchSectorInventory {
            sector_index: sector,
            sector_name,
            strip_sides: sides,
            strip_names,
            inner_patch_name,
            interface_names,
        });
    }

    let attempted_strip_patch_count = sector_sides.len();
    let attempted_inner_patch_count = sector_reports.len();
    let attempted_interface_count = attempted_strip_patch_count;
    Kp20MultiPatchInventoryReport {
        attempted_sector_count: attempted_inner_patch_count,
        attempted_strip_patch_count,
        attempted_inner_patch_count,
        attempted_interface_count,
        strip_control_point_count: attempted_strip_patch_count * BI6_CONTROL_POINTS,
        inner_control_point_count: attempted_inner_patch_count * BI6_CONTROL_POINTS,
        inventory_mismatch: inventory_mismatch(&sector_reports),
        sector_reports,
    }
}

fn unique_sectors(sector_sides: &[(usize, Kp20PieceId)]) -> Vec<usize> {
    let mut sectors = sector_sides
        .iter()
        .map(|(sector, _)| *sector)
        .collect::<Vec<_>>();
    sectors.sort_unstable();
    sectors.dedup();
    sectors
}

fn inventory_mismatch(sector_reports: &[Kp20MultiPatchSectorInventory]) -> Option<&'static str> {
    for sector in sector_reports {
        if sector.strip_names.is_empty() {
            return Some("missing_strip");
        }
        if sector.interface_names.len() != sector.strip_names.len() {
            return Some("missing_interface");
        }
        if has_duplicate_name(&sector.strip_names) || has_duplicate_name(&sector.interface_names) {
            return Some("duplicate_name");
        }
    }
    None
}

fn has_duplicate_name(names: &[String]) -> bool {
    for (idx, name) in names.iter().enumerate() {
        if names.iter().skip(idx + 1).any(|other| other == name) {
            return true;
        }
    }
    false
}

fn side_label(piece: Kp20PieceId) -> &'static str {
    match piece {
        Kp20PieceId::NearMidpointLeft => "left",
        Kp20PieceId::NearMidpointRight => "right",
        Kp20PieceId::NearOriginLeft | Kp20PieceId::NearOriginRight => "inner",
    }
}

const BI6_CONTROL_POINTS: usize = 49;

#[cfg(test)]
mod tests {
    use popkern_core::geometry::gregory::Kp20PieceId;

    use super::*;
    use crate::construct::kp20_polar_fold_summary::{Kp20FoldSite, Kp20SamplerFoldSummary};

    #[test]
    fn inventory_uses_semantic_support_strip_and_inner_names() {
        let summary = Kp20SamplerFoldSummary {
            fold_count: 2,
            min_jacobian_normal_cos: -0.1,
            piece_fold_counts: [1, 0, 1, 0],
            fold_details: Vec::new(),
            fold_sites: vec![
                fold_site(4, Kp20PieceId::NearOriginLeft),
                fold_site(4, Kp20PieceId::NearMidpointRight),
            ],
        };

        let inventory = multipatch_inventory_from_fold_summary(&summary);

        assert_eq!(inventory.attempted_sector_count, 1);
        assert_eq!(inventory.attempted_strip_patch_count, 2);
        assert_eq!(inventory.attempted_inner_patch_count, 1);
        assert_eq!(inventory.attempted_interface_count, 2);
        assert_eq!(inventory.strip_control_point_count, 98);
        assert_eq!(inventory.inner_control_point_count, 49);
        assert_eq!(inventory.inventory_mismatch, None);
        let sector = &inventory.sector_reports[0];
        assert_eq!(sector.sector_index, 4);
        assert_eq!(sector.sector_name, "kp20_polar_macro.sector_4");
        assert_eq!(
            sector.strip_sides,
            vec![
                Kp20PieceId::NearMidpointRight,
                Kp20PieceId::NearMidpointLeft
            ]
        );
        assert!(
            sector
                .strip_names
                .iter()
                .any(|name| name.ends_with("_support_strip"))
        );
        assert_eq!(
            sector.inner_patch_name,
            "kp20_polar_macro.sector_4.inner_pole_patch"
        );
        assert!(
            sector
                .interface_names
                .iter()
                .any(|name| name.ends_with("_support_interface"))
        );
    }

    fn fold_site(sector: usize, piece: Kp20PieceId) -> Kp20FoldSite {
        Kp20FoldSite {
            sector,
            piece,
            u: 0.5,
            v: 0.5,
            min_jacobian_ratio: -0.1,
            min_jacobian_normal_cos: -0.1,
        }
    }
}
