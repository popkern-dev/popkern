//! Public KP20 re-export surface for the Gregory module.

pub use super::kp20_bi6_assemble::{
    Kp20Bi6AssembleError, Kp20Bi6CNetCapInput, Kp20Bi6CapInput, Kp20Bi6MacroCapInput,
    Kp20Bi6MacroCapOutput, assemble_kp20_bi6_cap, assemble_kp20_bi6_cap_from_cnet,
    assemble_kp20_bi6_macro_cap,
};
pub use super::kp20_bi6_cap::{
    Kp20Bi6Patch, Kp20BoundaryHalf, Kp20MacroPatch, Kp20PieceId, Kp20SplitTensorBorder,
    initialise_kp20_macro_patches_from_tensor_border, kp20_apply_g2_strip,
    kp20_propagate_c2_upper_from_lower, split_kp20_tensor_border,
};
pub use super::kp20_bi6_cnet::{
    CUBIC_BSPLINE_TO_BEZIER, Kp20CapNet, Kp20RegularCNet, cubic_bspline_grid_to_bezier,
    kp20_star_point_from_cap_net, kp20_tensor_border_from_regular_cnet,
};
pub use super::kp20_bi6_guide_build::{Kp20GuideBuildInput, Kp20GuideRing, build_kp20_guide_ring};
pub use super::kp20_bi6_guide_continuity::{
    kp18_apply_c0_join_degree6, kp18_apply_c1_join_degree6, kp18_apply_c2_join_degree6,
    kp18_apply_c3_join_degree6, kp18_apply_cubic_join_degree6, kp18_guide_join_weights,
};
pub use super::kp20_bi6_guide_lsq::{
    Kp20GuideLsqError, Kp20GuideLsqOutput, Kp20GuideLsqReport, solve_kp20_guide_lsq,
    solve_kp20_guide_lsq_with_prior, solve_kp20_guide_lsq_with_report,
};
pub use super::kp20_bi6_guide_sample::{
    Kp20AffineGuideMap, Kp20GuideJet3, Kp20PatchCorner, Kp20Tri6GuideSector,
    kp20_corner_block_from_jet, kp20_sample_corner_block, kp20_sample_macro_patch_from_guide,
    kp20_sample_patch_from_guide,
};
pub use super::kp20_bi6_math::{
    HERMITE_TO_BI5_W, HERMITE_TO_BI6_W, Kp20CharacteristicScalars, Kp20TauParams,
    hermite_to_bi5_corner_block, hermite_to_bi6_corner_block, kappa_for_valence, kp18_c2_join_left,
    kp18_c2_join_right, kp18_guide_boundary_free_count, kp18_guide_interior_free_count,
    kp18_guide_lsq_free_count, kp18_guide_total_free_count, kp20_b_coeffs_near_m,
    kp20_b_coeffs_near_o, kp20_cc_subdominant_eigenvalue, kp20_g2_strip_weights,
    kp20_limit_weights, kp20_refined_s0_free_count, kp20_refinement_added_free_count,
    kp20_s0_free_count, kp20_stam_characteristic_eigenvector,
};
pub use super::kp20_bi6_piece_maps::{
    Kp20PieceCornerMaps, Kp20PieceMapDomainReport, Kp20SectorPieceMaps,
    kp20_c2_stitched_sector_piece_maps, kp20_coupled_stitched_sector_piece_maps,
    kp20_default_piece_map_domain_report, kp20_default_sector_piece_maps,
    kp20_nonoverlap_coupled_sector_piece_maps, kp20_orientation_valid_sector_piece_maps,
    kp20_paper_oriented_sector_piece_maps, kp20_piece_map_domain_report,
    kp20_piece_map_stitch_report, kp20_positive_nonoverlap_c2_sector_piece_maps,
    sample_kp20_macro_patch_from_default_maps, sample_kp20_macro_patch_from_piece_maps,
    sample_kp20_macro_patch_from_tau_maps,
};
pub use super::kp20_bi6_s0_solve::{
    Kp20S0Projection, Kp20S0Sample, Kp20S0SolveError, build_kp20_g2_strip_samples,
    project_kp20_g2_strip_pair, solve_kp20_s0_projection,
};
pub use super::kp20_bi6_tau::{
    KP20_TAU_CP_COUNT, Kp20TauJet3, Kp20TauMap, kp20_compose_guide_tau_jet,
    kp20_sample_corner_block_with_tau, kp20_sample_corner_block_with_tau_affine,
};
pub use super::kp20_bi6_validate::{
    Kp20Bi6ValidationError, Kp20Bi6ValidationReport, measure_kp20_bi6_macro_cap,
};
