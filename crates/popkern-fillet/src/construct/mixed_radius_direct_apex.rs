//! Mixed-radius direct-apex dispatch policy.

use super::EdgeGeometry;

pub(super) fn should_promote_mixed_radius_gregory_apex(
    edge_indices: &[usize],
    edge_geometries: &[EdgeGeometry],
    use_rational_gregory: bool,
    use_kp_binonic_gregory: bool,
) -> bool {
    if !use_rational_gregory || edge_indices.len() != 3 {
        return false;
    }
    let radii = edge_indices
        .iter()
        .filter_map(|&index| {
            edge_geometries
                .get(index)
                .map(|geometry| geometry.effective_radius)
        })
        .collect::<Vec<_>>();
    let Some((&first, rest)) = radii.split_first() else {
        return false;
    };
    if rest.iter().all(|radius| (*radius - first).abs() <= 1.0e-9) {
        return false;
    }
    if use_kp_binonic_gregory {
        let min_radius = radii.iter().copied().fold(f64::INFINITY, f64::min);
        let max_radius = radii.iter().copied().fold(0.0_f64, f64::max);
        if min_radius <= 0.0 {
            return false;
        }
        // Exact-boundary direct-apex remains the preferred mixed-radius
        // construction when the radius band is modest. Explicit K-P/binonic
        // requests own high-asymmetry trihedra so the legacy K-P ladder stays
        // measurable instead of being masked by automatic direct-apex routing.
        return max_radius / min_radius <= 2.0;
    }
    true
}
