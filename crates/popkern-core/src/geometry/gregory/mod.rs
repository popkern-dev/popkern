//! Rational Gregory patches — N=4 rectangular and N=3 triangular.
//!
//! Gregory patches (Gregory 1974, Chiyokura-Kimura 1983) resolve
//! the **twist compatibility problem** that arises at corners of
//! ordinary bicubic Bezier patches: at a shared corner, the two
//! incident boundaries each impose a different cross-tangent
//! direction, and an ordinary tensor-product Bezier cannot satisfy
//! both simultaneously. A rational Gregory patch stores *two*
//! candidate inner control points at each corner and blends them
//! rationally by a weight that vanishes on one boundary and grows
//! on the other — giving each boundary its correct cross-tangent
//! plane and thereby G1 continuity with arbitrary neighbours.
//!
//! # Representation
//!
//! * [`RectangularGregory`] is a bicubic-Bezier-domain patch on
//!   `[0, 1] × [0, 1]` with 4×4 = 16 boundary control points and 4
//!   **pairs** of inner control points (one pair per corner — the
//!   Chiyokura-Kimura form). Used for N=4 corner resolution and
//!   for three-quad-split N=3.
//! * [`TriangularGregory`] is a cubic-triangular-Bezier-domain
//!   patch on the barycentric simplex `{(u, v, w) : u + v + w = 1,
//!   u, v, w ≥ 0}` with 9 boundary control points and 3 inner
//!   **candidates** (one per corner — the Longhi form). Used for
//!   native N=3 corner resolution.
//!
//! In both forms the inner control point at evaluation is a
//! rational blend — two-way for rectangular (per corner), three-
//! way for triangular (global, weighted by barycentric products).
//! The blending weights are C∞-smooth on the open interior and
//! have removable singularities at the corners where the
//! Bernstein basis of the interior control point also vanishes,
//! so the inner CP does not contribute and the limit is chosen
//! as a tie-break without affecting the patch value.
//!
//! # References
//!
//! * Gregory, J. A. (1974). "Smooth Interpolation Without Twist
//!   Constraints." *Computer Aided Geometric Design*.
//! * Chiyokura, H., & Kimura, F. (1983). "Design of Solids with
//!   Free-form Surfaces." *Computer Graphics*, 17(3).
//! * Farin, G. (1993). "Curves and Surfaces for CAGD", Ch. 18.
//! * Piegl, L., & Tiller, W. (1997). "The NURBS Book", §9.
use crate::geometry::{Point3, Vec3};
use serde::{Deserialize, Serialize};
mod bernstein_bi3;
mod bernstein_bi5;
mod bernstein_bi6;
mod bernstein_binonic;
mod bernstein_tri5;
mod bernstein_tri6;
mod build;
mod build_biquintic_full;
mod build_biseptic;
mod build_two_stripe;
mod cc_characteristic_ring;
mod cc_characteristic_tensor_border;
mod direct_blend_coeffs;
mod direct_build;
mod direct_corner_c;
mod direct_corner_outer;
mod direct_errors;
mod direct_input;
mod direct_linalg;
mod direct_options;
mod direct_row1_g1;
mod direct_row1_ring;
mod direct_self_check;
mod direct_types;
mod kp13_bi5_assemble;
mod kp13_bi5_bi3_bi5_ops;
mod kp13_bi5_boundary;
mod kp13_bi5_cap;
mod kp13_bi5_cap_connect;
mod kp13_bi5_cap_f1;
mod kp13_bi5_cap_final;
mod kp13_bi5_cap_inner;
mod kp13_bi5_cap_outer;
mod kp13_bi5_chi;
mod kp13_bi5_circulant;
mod kp13_bi5_connect;
mod kp13_bi5_eq6;
mod kp13_bi5_free_params;
mod kp13_bi5_gtri3;
mod kp13_bi5_gtri5;
mod kp13_bi5_guide;
#[cfg(test)]
#[path = "kp13_bi5_integration_tests.rs"]
mod kp13_bi5_integration_tests;
mod kp13_bi5_jet;
mod kp13_bi5_low_valence;
mod kp13_bi5_lsq;
mod kp13_bi5_lsq_cc;
mod kp13_bi5_reparam;
mod kp13_bi5_sigma;
mod kp13_bi5_sigma_solve;
mod kp13_bi5_stam_chi;
mod kp13_bi5_star_quadratic;
mod kp15_bi6_cap;
mod kp15_bi6_cap_f5_solve;
mod kp15_bi6_circulant;
mod kp15_bi6_deep_interior;
mod kp15_bi6_f5;
mod kp15_bi6_local_constraint;
mod kp15_bi6_mu;
mod kp15_bi6_polar_bridge;
mod kp15_bi6_reparam;
mod kp15_bi6_rule_s3;
mod kp15_bi6_rules;
mod kp15_bi6_tensor_border_extract;
mod kp20_bi6_assemble;
mod kp20_bi6_cap;
mod kp20_bi6_cc_ring;
mod kp20_bi6_cnet;
mod kp20_bi6_guide_build;
mod kp20_bi6_guide_continuity;
mod kp20_bi6_guide_lsq;
mod kp20_bi6_guide_sample;
mod kp20_bi6_local_constraint;
mod kp20_bi6_math;
mod kp20_bi6_piece_maps;
mod kp20_bi6_s0_solve;
#[cfg(test)]
mod kp20_bi6_stam_eval;
mod kp20_bi6_tau;
mod kp20_bi6_tau_f3_solve;
mod kp20_bi6_tau_solve_assemble;
mod kp20_bi6_validate;
mod kp20_border_degree;
mod kp20_paper_border_piece;
mod kp20_paper_cap;
mod kp20_paper_cap_measure;
mod kp20_paper_cap_sample;
mod kp20_paper_cap_strip;
mod kp20_paper_guide;
mod kp20_paper_guide_fairing;
mod kp20_paper_tau;
mod kp20_paper_tau_measure;
mod kp_boundary_consistency;
mod kp_degree_elevate;
mod kp_guide_build;
mod kp_guide_central_hermite;
mod kp_guide_central_twist;
mod kp_guide_hex_hermite;
mod kp_guide_outer_hermite;
mod kp_guide_outer_twist;
mod kp_guide_patch;
mod kp_guide_patch_binonic;
mod kp_guide_patch_biquintic;
mod kp_guide_ring;
mod kp_guide_trait;
mod kp_interior_fit;
mod kp_polar_boundary_input;
mod kp_polar_shared;
mod kp_ring_build;
mod kp_row1;
mod kp_row2;
mod normal_curvature;
mod patch_dispatch;
mod pns_beta;
mod pns_bi3_twists;
mod pns_bi4_to_nurbs;
mod pns_boundary;
mod pns_cap_assemble;
mod pns_center;
mod pns_central_quad;
mod pns_f_kappa;
mod pns_g1_constraints;
mod pns_gram;
mod pns_sigma;
mod rectangular;
mod rectangular_binonic;
mod rectangular_binonic_dual_jet;
mod rectangular_binonic_dual_jet_to_nurbs;
mod rectangular_binonic_dual_twist;
mod rectangular_binonic_to_nurbs;
mod rectangular_biquintic_full;
mod rectangular_biquintic_full_to_nurbs;
mod rectangular_biseptic;
mod rectangular_biseptic_to_nurbs;
mod rectangular_deg5;
mod rectangular_deg5_to_nurbs;
mod rectangular_gregory_two_stripe;
mod rectangular_gregory_two_stripe_to_nurbs;
mod rectangular_to_nurbs;
mod ring_pair_crossing;
mod spoke_g1_solver;
mod spoke_g1_solver_impls;
mod spoke_poly_g1_linalg;
mod spoke_poly_g1_ring;
mod spoke_poly_g1_solver;
mod triangular;
mod triangular_to_nurbs;
pub mod vertex_enclosure;
pub use bernstein_bi3::{
    bernstein_bi3, bernstein_bi3_all, bernstein_bi3_deriv, bernstein_bi3_deriv_all,
    bernstein_bi3_second_deriv, bernstein_bi3_second_deriv_all, bernstein_bi3_third_deriv,
};
pub use bernstein_bi5::{
    bernstein_bi2_all, bernstein_bi4_all, bernstein_bi5, bernstein_bi5_all, bernstein_bi5_deriv,
    bernstein_bi5_deriv_all, bernstein_bi5_second_deriv, bernstein_bi5_second_deriv_all,
    bernstein_bi5_third_deriv, bernstein_bi5_third_deriv_all,
};
pub use bernstein_bi6::{
    bernstein_bi6, bernstein_bi6_all, bernstein_bi6_deriv, bernstein_bi6_deriv_all,
    bernstein_bi6_fifth_deriv_all, bernstein_bi6_fourth_deriv_all, bernstein_bi6_second_deriv,
    bernstein_bi6_second_deriv_all, bernstein_bi6_third_deriv, bernstein_bi6_third_deriv_all,
};
pub use bernstein_binonic::{
    bernstein_binonic, bernstein_binonic_all, bernstein_binonic_deriv, bernstein_binonic_deriv_all,
    bernstein_binonic_second_deriv, bernstein_binonic_second_deriv_all,
};
pub use bernstein_tri5::{
    TRI3_CP_COUNT, TRI4_CP_COUNT, TRI5_CP_COUNT, bernstein_tri3_all, bernstein_tri4_all,
    bernstein_tri5_all, bernstein_tri5_partial_u, bernstein_tri5_partial_v,
    bernstein_tri5_second_partial_uu, bernstein_tri5_second_partial_uv,
    bernstein_tri5_second_partial_vv, flat_index, tri_cp_count,
};
pub use bernstein_tri6::{
    TRI6_CP_COUNT, bernstein_tri6_all, bernstein_tri6_all_extended, bernstein_tri6_partial_u,
    bernstein_tri6_partial_v, bernstein_tri6_second_partial_uu, bernstein_tri6_second_partial_uv,
    bernstein_tri6_second_partial_vv, bernstein_tri6_third_partial_uuu,
    bernstein_tri6_third_partial_uuv, bernstein_tri6_third_partial_uvv,
    bernstein_tri6_third_partial_vvv, flat_index as tri6_flat_index,
};
pub use build::GregoryBuildError;
pub use cc_characteristic_ring::{
    CcCharacteristicRing, CcComplex, cc_characteristic_ring, cc_extended_net_len,
    cc_extended_subdivision_rows, cc_level0_subdivision_rows, cc_picking_vectors,
    cc_refined_net_len, cc_subdominant_lambda,
};
pub use cc_characteristic_tensor_border::{
    CcBorderPiece, CcCharacteristicTensorBorder, CcSectorBorder, cc_characteristic_tensor_border,
    cc_lifted_ring_border_pieces,
};
pub use direct_blend_coeffs::{SpokeBlendCoeffs, compute_ring_blend_coeffs};
pub use direct_build::{build, can_build};
pub use direct_corner_c::{
    CentralCornerInnerCps, central_plane_from_input, solve_central_vertex_corner,
};
pub use direct_corner_outer::{OuterCornerOutput, OuterCornerTwistCps, solve_outer_corners};
pub use direct_errors::DirectConstructionError;
pub use direct_input::{
    DirectBoundaryInput, HexBoundaryData, OuterCornerData, SpokeBoundaryData, TangentPlane,
};
pub use direct_options::DirectConstructionOptions;
pub use direct_row1_g1::{SingleSpokeRow1Input, SingleSpokeRow1Output, solve_single_spoke_row1_g1};
pub use direct_row1_ring::{RingRow1Result, solve_ring_row1_g1};
pub use direct_self_check::{
    PerSpokeQuality, RingQualityMeasurement, SelfCheckOptions, measure_ring_quality,
    measure_ring_quality_per_spoke,
};
pub use direct_types::{
    CornerDiagnostic, CornerKind, DirectConstructionResult, FeasibilityReport, FeasibilityWarning,
    FormulaFamily, SpokeDiagnostic, ValidationLevel,
};
pub use kp_boundary_consistency::{
    BoundaryConsistency, CornerConsistency, check_boundary_consistency,
};
pub use kp_degree_elevate::{
    elevate_bi5_grid_to_bi9, elevate_bi6_grid_to_bi9, elevate_cubic_point_to_nine,
    elevate_cubic_vec_to_nine, elevate_point_curve_to_nine, elevate_quintic_point_to_nine,
    elevate_quintic_vec_to_nine, elevate_sextic_point_to_nine, elevate_sextic_vec_to_nine,
};
pub use kp_guide_build::{GuideCornerData, build_guide_from_hermite};
pub use kp_guide_central_hermite::extract_central_vertex_hermite;
pub use kp_guide_central_twist::extract_central_vertex_hermite_with_twist;
pub use kp_guide_hex_hermite::{SubPatchHexCorners, extract_hex_corners};
pub use kp_guide_outer_hermite::{OuterVertexHermite, extract_outer_corner_hermite};
pub use kp_guide_outer_twist::extract_hex_corners_with_twist;
pub use kp_guide_patch::KpGuidePatch;
pub use kp_guide_patch_binonic::KpGuidePatchBiNonic;
pub use kp_guide_patch_biquintic::KpGuidePatchBiQuintic;
pub use kp_guide_ring::{
    build_guide_ring, build_guide_ring_from_centrals, build_guide_ring_from_hermite,
    build_guide_ring_with_twist,
};
pub use kp_guide_trait::KpGuide;
pub use kp_interior_fit::{
    fit_interior_cps_to_binonic_guide, fit_interior_cps_to_biquintic_guide,
    fit_interior_cps_to_guide,
};
pub use kp_polar_boundary_input::{
    PolarBoundaryInput, TensorBorder, TensorBorderBi3, TensorBorderBi5, TensorBorderKind,
};
pub use kp_polar_shared::catmull_clark_star_point;
pub use kp_ring_build::{
    build_binonic_ring, build_binonic_ring_from_binonic_guides,
    build_binonic_ring_from_biquintic_guides, build_binonic_ring_from_guides,
};
pub use kp_row1::{Row1Frame, derive_row1_frame};
pub use kp_row2::{Row2Frame, derive_row2_frame};
pub use kp13_bi5_assemble::{Bi5AssembleError, Bi5CapInput, Bi5CapOutput, assemble_bi5_macro_cap};
pub use kp13_bi5_boundary::{
    Bi3TensorBorder as Bi5CapBi3Patch, TensorBorder as Bi5CapTensorBorder,
};
pub use kp15_bi6_cap::{
    Bi6CapInput, Bi6CapOutput, assemble_bi6_polar_cap, assemble_bi6_polar_cap_via_paper_p_table,
    assemble_bi6_polar_cap_via_table,
};
pub use kp15_bi6_cap_f5_solve::{
    Bi6CapFreeParams, Bi6CapPaperFreeParams, f_kappa_solve_optimal_bi6_cap,
    f_kappa_solve_paper_p_bi6_cap, f_kappa_solve_with_weights_optimal_bi6_cap,
    f5_solve_optimal_bi6_cap,
};
pub use kp15_bi6_circulant::{
    CirculantSectorRhs, S2Eq7Result, solve_bi6_circulant_s2_eq7, solve_bi6_circulant_s2_eq8,
};
pub use kp15_bi6_deep_interior::{
    Bi6CapScriptGrids, Bi6DeepInteriorFree, SectorTensorBorderBi3, assemble_bi6_polar_cap_from_net,
    build_bi6_cap_from_net, build_bi6_cap_from_tensor_border,
};
pub use kp15_bi6_f5::{
    BI6_CP_COUNT, f_kappa_quadratic_form, f_kappa_quadratic_form_weighted, f5_evaluate,
    f5_quadratic_form,
};
pub use kp15_bi6_local_constraint::{Kp15Bi6GeneratingTable, apply_beta_reparameterization};
pub use kp15_bi6_mu::{mu_bar, mu_ddot, mu_dot, mu_tdot, mu_tilde, w_scalars};
pub use kp15_bi6_polar_bridge::{convert_polar_to_bi6_cap_fitted, convert_polar_to_bi6_cap_input};
pub use kp15_bi6_reparam::{Bi6Reparams, gamma as bi6_gamma, tau_for_valence};
pub use kp15_bi6_rule_s3::{S3Inputs, S3Outputs, rule_s3};
pub use kp15_bi6_rules::{
    S1Outputs, S4Inputs, S4Outputs, rule_5_reduce_edge, rule_s1_star_propagation, rule_s2_p12,
    rule_s2_p21, rule_s2_p22, rule_s4,
};
pub use kp15_bi6_tensor_border_extract::extract_bi3_tensor_border_from_polar_input;
// KP20 macro-cap public surface. The full item list lives in
// `kp20_bi6_exports.rs` — a dedicated file whose only content is the
// explicit `pub use` surface (kept out of this file for the 500-line
// module budget; rustdoc inlines it identically).
mod kp20_bi6_exports;
pub use kp20_bi6_exports::*;
pub use kp20_border_degree::{
    SexticRow, bernstein6, elevate_bi3_rows_to_sextic, elevate_cubic_row_to_sextic,
    eval_sextic_row, reduce_sextic_row_to_cubic, reduce_sextic_rows_to_bi3, reverse_sextic_rows,
    split_sextic_row, split_sextic_rows,
};
pub use kp20_paper_cap::{
    Kp20PaperCapError, Kp20PaperCapReport, assemble_kp20_paper_cap,
    assemble_kp20_paper_cap_conditioned, assemble_kp20_paper_cap_with_enforcement,
    assemble_kp20_paper_macro_cap,
};
pub use kp20_paper_cap_strip::Kp20StripEnforcement;
pub use kp20_paper_guide::{
    Kp20PaperGuideBorder, Kp20PaperGuideBorderPiece, Kp20PaperGuideOutput, solve_kp20_paper_guide,
};
pub use kp20_paper_tau::{
    Kp20PaperTauReport, kp20_paper_border_corner_mismatch, solve_kp20_paper_tau,
};
pub use normal_curvature::{NormalCurvatureComponents, NormalCurvatureError};
pub use normal_curvature::{normal_curvature_along, normal_curvature_components_along};
pub use patch_dispatch::GregoryPatch;
pub use pns_beta::BetaReparam;
pub use pns_bi4_to_nurbs::pns_bi4_to_nurbs;
pub use pns_boundary::{Bi2BorderStrip, Bi2TensorBorderRing};
pub use pns_cap_assemble::{
    CapAssembleError, CapAssembleInput, CapAssembleResult, assemble_pns_bi4_cap, pole_from_boundary,
};
pub use pns_central_quad::QuadMapJet;
pub use pns_sigma::{SigmaCenterJet, SigmaJetError, sigma_center_jet, sigma_center_jet_default};
pub use rectangular::RectangularGregory;
pub use rectangular_binonic::RectangularGregoryBiNonic;
pub use rectangular_binonic_dual_jet::{
    DUAL_JET_CORNER_POSITIONS, DualJetCell, DualJetCorner, RectangularGregoryBiNonicDualJet,
    dual_jet_cell_of, dual_jet_cell_position,
};
pub use rectangular_binonic_dual_twist::RectangularGregoryBiNonicDualTwist;
pub use rectangular_biquintic_full::RectangularGregoryBiquinticFull;
pub use rectangular_biseptic::RectangularGregoryBiseptic;
pub use rectangular_deg5::RectangularGregoryDeg5;
pub use rectangular_gregory_two_stripe::RectangularGregoryTwoStripe;
pub use ring_pair_crossing::{
    BinonicRingCrossing, RingPairRelation, binonic_ring_crossing,
    binonic_ring_feasible_start_crossing,
};
pub use self_intersection::GregorySelfIntersectionReport;
pub use spoke_g1_solver::{
    Side, Spoke, SpokeEnd, SpokeG1Patch, spoke_g1_solve, spoke_g1_solve_with_corner_samples,
    spoke_g1_solve_with_corner_samples_damped,
};
pub use spoke_poly_g1_ring::{RingPolyG1Result, solve_ring_poly_g1};
pub use triangular::TriangularGregory;

#[cfg(test)]
#[path = "build_tests.rs"]
mod build_tests;
#[cfg(test)]
#[path = "rectangular_binonic_dual_twist_tests.rs"]
mod rectangular_binonic_dual_twist_tests;
#[cfg(test)]
#[path = "rectangular_biquintic_full_tests.rs"]
mod rectangular_biquintic_full_tests;
#[cfg(test)]
#[path = "rectangular_biseptic_tests.rs"]
mod rectangular_biseptic_tests;
#[cfg(test)]
#[path = "rectangular_deg5_tests.rs"]
mod rectangular_deg5_tests;
#[cfg(test)]
#[path = "rectangular_gregory_two_stripe_tests.rs"]
mod rectangular_gregory_two_stripe_tests;
#[cfg(test)]
#[path = "rectangular_to_nurbs_tests.rs"]
mod rectangular_to_nurbs_tests;
#[cfg(test)]
#[path = "ring_pair_crossing_tests.rs"]
mod ring_pair_crossing_tests;
mod self_intersection;
#[cfg(test)]
#[path = "spoke_g1_solver_higher_degree_tests.rs"]
mod spoke_g1_solver_higher_degree_tests;
#[cfg(test)]
#[path = "spoke_g1_solver_tests.rs"]
mod spoke_g1_solver_tests;
#[cfg(test)]
#[path = "../gregory_tests.rs"]
mod tests;
#[cfg(test)]
#[path = "triangular_to_nurbs_tests.rs"]
mod triangular_to_nurbs_tests;
