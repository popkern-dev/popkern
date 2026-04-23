//! Rational Gregory patches — N=4 rectangular and N=3 triangular.
//!
//! Gregory patches (Gregory 1974, Chiyokura-Kimura 1983) resolve
//! the **twist compatibility problem** that arises at corners of
//! ordinary bicubic Bezier patches: at a shared corner, the two
//! incident boundaries each impose a different cross-tangent
//! direction, and an ordinary tensor-product Bezier cannot
//! satisfy both simultaneously. A rational Gregory patch stores
//! *two* candidate inner control points at each corner and blends
//! them rationally by a weight that vanishes on one boundary and
//! grows on the other — giving each boundary its correct cross-
//! tangent plane and thereby G1 continuity with arbitrary
//! neighbours, regardless of valence.
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
mod build;
mod build_biquintic_full;
mod build_biseptic;
mod build_two_stripe;
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
mod kp15_bi6_f5;
mod kp15_bi6_mu;
mod kp15_bi6_polar_bridge;
mod kp15_bi6_reparam;
mod kp15_bi6_rule_s3;
mod kp15_bi6_rules;
mod kp15_bi6_tensor_border_extract;
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
mod pns_beta;
mod pns_center;
mod pns_g1_constraints;
mod rectangular;
mod rectangular_binonic;
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
mod spoke_g1_solver;
mod spoke_g1_solver_impls;
mod spoke_poly_g1_linalg;
mod spoke_poly_g1_ring;
mod spoke_poly_g1_solver;
mod triangular;
mod triangular_to_nurbs;

pub use build::GregoryBuildError;
pub use rectangular::RectangularGregory;
pub use rectangular_binonic::RectangularGregoryBiNonic;
pub use rectangular_biquintic_full::RectangularGregoryBiquinticFull;
pub use rectangular_biseptic::RectangularGregoryBiseptic;
pub use rectangular_deg5::RectangularGregoryDeg5;
pub use rectangular_gregory_two_stripe::RectangularGregoryTwoStripe;
pub use spoke_g1_solver::{
    Side, Spoke, SpokeEnd, SpokeG1Patch, spoke_g1_solve, spoke_g1_solve_with_corner_samples,
    spoke_g1_solve_with_corner_samples_damped,
};
pub use spoke_poly_g1_ring::{RingPolyG1Result, solve_ring_poly_g1};

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
pub use direct_blend_coeffs::{SpokeBlendCoeffs, compute_ring_blend_coeffs};
pub use direct_build::{build, can_build};
pub use direct_corner_c::{
    CentralCornerInnerCps, central_plane_from_input, solve_central_vertex_corner,
};
pub use direct_corner_outer::{OuterCornerOutput, OuterCornerTwistCps, solve_outer_corners};
pub use direct_errors::DirectConstructionError;
pub use direct_input::{DirectBoundaryInput, OuterCornerData, SpokeBoundaryData, TangentPlane};
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
pub use kp15_bi6_cap::{Bi6CapInput, Bi6CapOutput, assemble_bi6_polar_cap};
pub use kp15_bi6_cap_f5_solve::{
    Bi6CapFreeParams, f_kappa_solve_optimal_bi6_cap, f_kappa_solve_with_weights_optimal_bi6_cap,
    f5_solve_optimal_bi6_cap,
};
pub use kp15_bi6_circulant::{
    CirculantSectorRhs, S2Eq7Result, solve_bi6_circulant_s2_eq7, solve_bi6_circulant_s2_eq8,
};
pub use kp15_bi6_f5::{
    BI6_CP_COUNT, f_kappa_quadratic_form, f_kappa_quadratic_form_weighted, f5_evaluate,
    f5_quadratic_form,
};
pub use kp15_bi6_mu::{mu_bar, mu_ddot, mu_dot, mu_tdot, mu_tilde, w_scalars};
pub use kp15_bi6_polar_bridge::convert_polar_to_bi6_cap_input;
pub use kp15_bi6_reparam::{Bi6Reparams, gamma as bi6_gamma, tau_for_valence};
pub use kp15_bi6_rule_s3::{S3Inputs, S3Outputs, rule_s3};
pub use kp15_bi6_rules::{
    S1Outputs, S4Inputs, S4Outputs, rule_5_reduce_edge, rule_s1_star_propagation, rule_s2_p12,
    rule_s2_p21, rule_s2_p22, rule_s4,
};
pub use kp15_bi6_tensor_border_extract::extract_bi3_tensor_border_from_polar_input;
pub use triangular::TriangularGregory;

/// A Gregory patch of either rectangular (N=4) or triangular
/// (N=3) topology.
///
/// Boxed because the fixed-size control-point arrays on both
/// variants are large (~400–600 bytes) and including them inline
/// would inflate every `SurfaceGeometry` enum instance to match;
/// Gregory patches are rare (only at corners) so the one extra
/// pointer indirection per evaluation is free compared to the
/// per-NURBS-face cache cost of a bloated enum.
///
/// Use [`SurfaceGeometry::Gregory`](crate::geometry::SurfaceGeometry::Gregory)
/// to attach to a face.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum GregoryPatch {
    /// N=4 rectangular cubic patch (bidegree 3×3).
    Rectangular(Box<RectangularGregory>),
    /// N=4 rectangular quintic patch (bidegree 5×5). Same
    /// Chiyokura-Kimura 4-inner-CP-pair structure as
    /// [`Rectangular`](GregoryPatch::Rectangular); higher
    /// boundary degree closes the 2% cubic-Bezier arc
    /// approximation error inherent at the lower degree.
    RectangularDeg5(Box<RectangularGregoryDeg5>),
    /// Karciauskas-Peters biquintic patch with 12 free-curvature
    /// interior CPs in addition to the 4 corner inner-CP pairs.
    /// Extends
    /// [`RectangularDeg5`](GregoryPatch::RectangularDeg5) with
    /// the degrees of freedom required for **geometric G1 along
    /// the entire shared spoke** between adjacent sub-patches
    /// (not just at the central vertex).
    RectangularBiquinticFull(Box<RectangularGregoryBiquinticFull>),
    /// Biseptic (degree 7) Karciauskas-Peters patch with 32
    /// free-curvature interior CPs. 4 free DOF per shared
    /// spoke (up from 2 at biquintic-full) gives the
    /// spoke-G1 solver full Bernstein-basis coverage for
    /// row `i = 1`, eliminating the structural out-of-span
    /// mismatch visible in biquintic-full's scorecard.
    /// Phase 10 extends this variant's solver with G2
    /// (curvature) continuity by leveraging row `i = 2`'s
    /// 6 free interior CPs per spoke.
    RectangularBiseptic(Box<RectangularGregoryBiseptic>),
    /// Peters two-stripe bidegree-7 Bezier patch. Drops the
    /// rational blend at the 4 corner-adjacent positions
    /// `(1,1), (1,6), (6,1), (6,6)` — they become plain-
    /// Bezier free CPs. The spoke-G1 solver now has 6 DOF
    /// per spoke side (row `i = 1` j ∈ {1..=6}), giving full
    /// Bernstein-basis coverage `{B_1^7..B_6^7}`. Corner
    /// tangent-plane compatibility that was previously
    /// enforced by the rational blend is handled by explicit
    /// solver constraints. This is the structural fix for
    /// the biseptic scorecard's 84° worst-case.
    RectangularTwoStripe(Box<RectangularGregoryTwoStripe>),
    /// Phase H — Karciauskas-Peters bidegree-(9, 9)
    /// guided-surface construction. 100 plain-Bezier CPs
    /// (no rational blend) with row-1 / row-2 frames derived
    /// from a bicubic guide that enforces G1 and G2 across
    /// every shared spoke by construction. Interior CPs
    /// (rows 3..=6, cols 3..=6) are least-squares-fit to
    /// the guide. See
    /// `docs/decisions/2026-04-19-phase-h-kp-binonic-formula-choice.md`.
    RectangularBiNonic(Box<RectangularGregoryBiNonic>),
    /// N=3 triangular patch.
    Triangular(Box<TriangularGregory>),
}

impl GregoryPatch {
    /// Evaluate at the given parameters. `(u, v)` is the Cartesian
    /// parameter for [`Rectangular`](GregoryPatch::Rectangular) and
    /// [`RectangularDeg5`](GregoryPatch::RectangularDeg5); for
    /// [`Triangular`](GregoryPatch::Triangular), `(u, v)` map into
    /// barycentric coords via `(u, v, 1 − u − v)` and must satisfy
    /// `u ≥ 0`, `v ≥ 0`, `u + v ≤ 1`.
    pub fn point_at(&self, u: f64, v: f64) -> Point3 {
        match self {
            GregoryPatch::Rectangular(r) => r.point_at(u, v),
            GregoryPatch::RectangularDeg5(r) => r.point_at(u, v),
            GregoryPatch::RectangularBiquinticFull(r) => r.point_at(u, v),
            GregoryPatch::RectangularBiseptic(r) => r.point_at(u, v),
            GregoryPatch::RectangularTwoStripe(r) => r.point_at(u, v),
            GregoryPatch::RectangularBiNonic(r) => r.point_at(u, v),
            GregoryPatch::Triangular(t) => t.point_at_uv(u, v),
        }
    }

    /// Partial derivative with respect to `u`. See
    /// [`point_at`](Self::point_at) for the parameter convention.
    pub fn partial_u(&self, u: f64, v: f64) -> Vec3 {
        match self {
            GregoryPatch::Rectangular(r) => r.partial_u(u, v),
            GregoryPatch::RectangularDeg5(r) => r.partial_u(u, v),
            GregoryPatch::RectangularBiquinticFull(r) => r.partial_u(u, v),
            GregoryPatch::RectangularBiseptic(r) => r.partial_u(u, v),
            GregoryPatch::RectangularTwoStripe(r) => r.partial_u(u, v),
            GregoryPatch::RectangularBiNonic(r) => r.partial_u(u, v),
            GregoryPatch::Triangular(t) => t.partial_u_uv(u, v),
        }
    }

    /// Partial derivative with respect to `v`. See
    /// [`point_at`](Self::point_at) for the parameter convention.
    pub fn partial_v(&self, u: f64, v: f64) -> Vec3 {
        match self {
            GregoryPatch::Rectangular(r) => r.partial_v(u, v),
            GregoryPatch::RectangularDeg5(r) => r.partial_v(u, v),
            GregoryPatch::RectangularBiquinticFull(r) => r.partial_v(u, v),
            GregoryPatch::RectangularBiseptic(r) => r.partial_v(u, v),
            GregoryPatch::RectangularTwoStripe(r) => r.partial_v(u, v),
            GregoryPatch::RectangularBiNonic(r) => r.partial_v(u, v),
            GregoryPatch::Triangular(t) => t.partial_v_uv(u, v),
        }
    }

    /// §5.3 item #10 — detect self-intersection of the patch to
    /// within `tolerance`. A patch self-intersects when (a) its
    /// Jacobian `|∂u × ∂v|` changes sign on the domain (the
    /// surface folds over itself), or (b) two non-adjacent
    /// sub-quadrant AABBs overlap deeper than `tolerance` after
    /// binary subdivision (the surface curls without folding).
    ///
    /// Cost: O(1) Jacobian grid + O(K²) subdivision-AABB check
    /// where K is the sub-quadrant count (16 at depth 2). Cheap
    /// enough for per-face runtime validation; see the
    /// [`validate_all`](crate::topology::BRep::validate_all)
    /// extension for the BRep-wide check.
    ///
    /// Returns `true` when a self-intersection is detected.
    pub fn self_intersects_with_tolerance(&self, tolerance: f64) -> bool {
        self_intersection::gregory_self_intersects(self, tolerance)
    }
}

mod self_intersection;

#[cfg(test)]
#[path = "../gregory_tests.rs"]
mod tests;

#[cfg(test)]
#[path = "build_tests.rs"]
mod build_tests;

#[cfg(test)]
#[path = "rectangular_deg5_tests.rs"]
mod rectangular_deg5_tests;

#[cfg(test)]
#[path = "rectangular_biquintic_full_tests.rs"]
mod rectangular_biquintic_full_tests;

#[cfg(test)]
#[path = "rectangular_biseptic_tests.rs"]
mod rectangular_biseptic_tests;

#[cfg(test)]
#[path = "rectangular_gregory_two_stripe_tests.rs"]
mod rectangular_gregory_two_stripe_tests;

#[cfg(test)]
#[path = "spoke_g1_solver_tests.rs"]
mod spoke_g1_solver_tests;

#[cfg(test)]
#[path = "spoke_g1_solver_higher_degree_tests.rs"]
mod spoke_g1_solver_higher_degree_tests;

#[cfg(test)]
#[path = "rectangular_to_nurbs_tests.rs"]
mod rectangular_to_nurbs_tests;

#[cfg(test)]
#[path = "triangular_to_nurbs_tests.rs"]
mod triangular_to_nurbs_tests;
