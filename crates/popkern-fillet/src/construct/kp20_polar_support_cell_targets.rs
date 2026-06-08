use popkern_core::geometry::gregory::Kp20PieceId;

use super::super::super::kp20_polar_fold_summary::Kp20FoldSite;
use super::super::super::kp20_polar_support_surface_split::support_owner_piece;

pub(super) fn support_targets_by_sector(
    fold_sites: &[Kp20FoldSite],
    target_piece: Kp20PieceId,
) -> Vec<Kp20FoldSite> {
    let mut targets = Vec::new();
    for site in fold_sites
        .iter()
        .copied()
        .filter(|site| support_owner_piece(site.piece) == target_piece)
    {
        if let Some(existing) = targets
            .iter_mut()
            .find(|existing: &&mut Kp20FoldSite| existing.sector == site.sector)
        {
            if site.min_jacobian_normal_cos < existing.min_jacobian_normal_cos {
                *existing = site;
            }
        } else {
            targets.push(site);
        }
    }
    targets
}
