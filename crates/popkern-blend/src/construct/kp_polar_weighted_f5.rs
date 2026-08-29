//! Weighted F5 fallback profiles for polar K-P bi-6 promotion.

use popkern_core::geometry::gregory::{
    Bi6CapFreeParams, Bi6CapInput, GregoryPatch, RectangularGregoryBiNonic,
    f_kappa_solve_with_weights_optimal_bi6_cap,
};

use super::kp_polar_bi6_diagnostics::trace_polar;
use super::kp_polar_bi6_promote::SELF_INT_TOL;

/// Option E — per-valence weight-tuning profiles. Try several
/// empirical re-weightings of F5's defaults and commit the first
/// non-folding cap.
pub(super) fn try_weighted_f5_profiles(
    cap_input: &Bi6CapInput,
    n: usize,
    last_free_params: &mut Option<Bi6CapFreeParams>,
) -> Option<Vec<GregoryPatch>> {
    let profiles: [(&str, [f64; 6]); 3] = [
        ("relax-pure", [0.5, 1.0, 1.0, 1.0, 1.0, 0.5]),
        ("amplify-mixed", [1.0, 1.0, 2.0, 2.0, 1.0, 1.0]),
        ("relax-mixed", [1.0, 0.5, 0.5, 0.5, 0.5, 1.0]),
    ];
    for (name, weights) in profiles.iter() {
        let (free_params, cap_out) =
            match f_kappa_solve_with_weights_optimal_bi6_cap(cap_input, 5, weights) {
                Ok(r) => r,
                Err(e) => {
                    trace_polar(&format!("F_5[{name}] solve error"), n, Some(e));
                    continue;
                }
            };
        *last_free_params = Some(free_params);
        let bi9_grids = cap_out.patches_bi9;
        if bi9_grids.len() != n {
            continue;
        }
        let candidate: Vec<GregoryPatch> = bi9_grids
            .into_iter()
            .map(|grid| {
                GregoryPatch::RectangularBiNonic(Box::new(RectangularGregoryBiNonic { grid }))
            })
            .collect();
        if candidate
            .iter()
            .any(|gp| gp.self_intersects_with_tolerance(SELF_INT_TOL))
        {
            trace_polar(&format!("F_5[{name}] self-intersection"), n, None);
            continue;
        }
        trace_polar(
            &format!("F_5[{name}] accepted"),
            n,
            Some("option E saved the ladder".into()),
        );
        return Some(candidate);
    }
    None
}
