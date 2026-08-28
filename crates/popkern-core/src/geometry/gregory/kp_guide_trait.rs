//! the bi-9 hardening pass (Option C step 2) — abstract-guide trait for the K-P
//! row-1 / row-2 frame construction.
//!
//! The bidegree-9 row-1 / row-2 formulas don't care about the
//! native degree of the guide surface — they consume the guide's
//! point + partial-derivative samples and its cross-derivative
//! CPs **expressed at degree 9**. This trait captures exactly
//! that interface, letting the same derivation code work for a
//! bicubic guide (`KpGuidePatch`, 4×4 CPs) and a bi-quintic guide
//! (`KpGuidePatchBiQuintic`, 6×6 CPs).
//!
//! Each impl elevates its native cross-derivative representation
//! to the shared degree-9 target — the bicubic goes 3 → 9 (elevate
//! by 6), the bi-quintic goes 5 → 9 (elevate by 4). Callers see a
//! uniform 10-CP Bezier and don't need to know the guide's native
//! degree.

use super::kp_degree_elevate::{elevate_cubic_vec_to_nine, elevate_quintic_vec_to_nine};
use super::kp_guide_patch::KpGuidePatch;
use super::kp_guide_patch_binonic::KpGuidePatchBiNonic;
use super::kp_guide_patch_biquintic::KpGuidePatchBiQuintic;
use crate::geometry::{Point3, Vec3};

/// Abstract guide surface for the K-P row-1 / row-2 construction.
///
/// Implementors provide direct point / partial evaluation at any
/// `(u, v)` plus the boundary cross-derivative CPs pre-elevated
/// to degree 9. The K-P formulas index these CPs as `[F(v)]_j^9`
/// where `j ∈ 0..=9` — the degree-9 Bernstein representation.
pub trait KpGuide {
    /// `G(u, v)`.
    fn point_at(&self, u: f64, v: f64) -> Point3;

    /// `∂G/∂u(u, v)`.
    fn partial_u(&self, u: f64, v: f64) -> Vec3;

    /// `∂G/∂v(u, v)`.
    fn partial_v(&self, u: f64, v: f64) -> Vec3;

    /// `∂²G/∂u²(u, v)`.
    fn second_partial_uu(&self, u: f64, v: f64) -> Vec3;

    /// `∂²G/∂u∂v(u, v)`.
    fn second_partial_uv(&self, u: f64, v: f64) -> Vec3;

    /// `∂²G/∂v²(u, v)`.
    fn second_partial_vv(&self, u: f64, v: f64) -> Vec3;

    /// 10-CP degree-9 Bernstein representation of `∂G/∂u(0, v)`.
    fn cross_deriv_cps_u0_degree9(&self) -> [Vec3; 10];
    /// 10-CP degree-9 Bernstein representation of `∂G/∂u(1, v)`.
    fn cross_deriv_cps_u1_degree9(&self) -> [Vec3; 10];
    /// 10-CP degree-9 Bernstein representation of `∂G/∂v(u, 0)`.
    fn cross_deriv_cps_v0_degree9(&self) -> [Vec3; 10];
    /// 10-CP degree-9 Bernstein representation of `∂G/∂v(u, 1)`.
    fn cross_deriv_cps_v1_degree9(&self) -> [Vec3; 10];

    /// 10-CP degree-9 Bernstein representation of `∂²G/∂u²(0, v)`.
    fn cross_second_deriv_cps_uu_u0_degree9(&self) -> [Vec3; 10];
    /// 10-CP degree-9 Bernstein representation of `∂²G/∂u²(1, v)`.
    fn cross_second_deriv_cps_uu_u1_degree9(&self) -> [Vec3; 10];
    /// 10-CP degree-9 Bernstein representation of `∂²G/∂v²(u, 0)`.
    fn cross_second_deriv_cps_vv_v0_degree9(&self) -> [Vec3; 10];
    /// 10-CP degree-9 Bernstein representation of `∂²G/∂v²(u, 1)`.
    fn cross_second_deriv_cps_vv_v1_degree9(&self) -> [Vec3; 10];
}

impl KpGuide for KpGuidePatch {
    fn point_at(&self, u: f64, v: f64) -> Point3 {
        KpGuidePatch::point_at(self, u, v)
    }
    fn partial_u(&self, u: f64, v: f64) -> Vec3 {
        KpGuidePatch::partial_u(self, u, v)
    }
    fn partial_v(&self, u: f64, v: f64) -> Vec3 {
        KpGuidePatch::partial_v(self, u, v)
    }
    fn second_partial_uu(&self, u: f64, v: f64) -> Vec3 {
        KpGuidePatch::second_partial_uu(self, u, v)
    }
    fn second_partial_uv(&self, u: f64, v: f64) -> Vec3 {
        KpGuidePatch::second_partial_uv(self, u, v)
    }
    fn second_partial_vv(&self, u: f64, v: f64) -> Vec3 {
        KpGuidePatch::second_partial_vv(self, u, v)
    }
    fn cross_deriv_cps_u0_degree9(&self) -> [Vec3; 10] {
        elevate_cubic_vec_to_nine(self.cross_deriv_cps_u0())
    }
    fn cross_deriv_cps_u1_degree9(&self) -> [Vec3; 10] {
        elevate_cubic_vec_to_nine(self.cross_deriv_cps_u1())
    }
    fn cross_deriv_cps_v0_degree9(&self) -> [Vec3; 10] {
        elevate_cubic_vec_to_nine(self.cross_deriv_cps_v0())
    }
    fn cross_deriv_cps_v1_degree9(&self) -> [Vec3; 10] {
        elevate_cubic_vec_to_nine(self.cross_deriv_cps_v1())
    }
    fn cross_second_deriv_cps_uu_u0_degree9(&self) -> [Vec3; 10] {
        elevate_cubic_vec_to_nine(self.cross_second_deriv_cps_uu_u0())
    }
    fn cross_second_deriv_cps_uu_u1_degree9(&self) -> [Vec3; 10] {
        elevate_cubic_vec_to_nine(self.cross_second_deriv_cps_uu_u1())
    }
    fn cross_second_deriv_cps_vv_v0_degree9(&self) -> [Vec3; 10] {
        elevate_cubic_vec_to_nine(self.cross_second_deriv_cps_vv_v0())
    }
    fn cross_second_deriv_cps_vv_v1_degree9(&self) -> [Vec3; 10] {
        elevate_cubic_vec_to_nine(self.cross_second_deriv_cps_vv_v1())
    }
}

impl KpGuide for KpGuidePatchBiNonic {
    fn point_at(&self, u: f64, v: f64) -> Point3 {
        KpGuidePatchBiNonic::point_at(self, u, v)
    }
    fn partial_u(&self, u: f64, v: f64) -> Vec3 {
        KpGuidePatchBiNonic::partial_u(self, u, v)
    }
    fn partial_v(&self, u: f64, v: f64) -> Vec3 {
        KpGuidePatchBiNonic::partial_v(self, u, v)
    }
    fn second_partial_uu(&self, u: f64, v: f64) -> Vec3 {
        KpGuidePatchBiNonic::second_partial_uu(self, u, v)
    }
    fn second_partial_uv(&self, u: f64, v: f64) -> Vec3 {
        KpGuidePatchBiNonic::second_partial_uv(self, u, v)
    }
    fn second_partial_vv(&self, u: f64, v: f64) -> Vec3 {
        KpGuidePatchBiNonic::second_partial_vv(self, u, v)
    }
    // Native degree 9 — no elevation needed. The cross-derivative
    // CPs returned by the patch are already in the target degree-9
    // representation. Directly usable by row-1 / row-2 formulas.
    fn cross_deriv_cps_u0_degree9(&self) -> [Vec3; 10] {
        self.cross_deriv_cps_u0()
    }
    fn cross_deriv_cps_u1_degree9(&self) -> [Vec3; 10] {
        self.cross_deriv_cps_u1()
    }
    fn cross_deriv_cps_v0_degree9(&self) -> [Vec3; 10] {
        self.cross_deriv_cps_v0()
    }
    fn cross_deriv_cps_v1_degree9(&self) -> [Vec3; 10] {
        self.cross_deriv_cps_v1()
    }
    fn cross_second_deriv_cps_uu_u0_degree9(&self) -> [Vec3; 10] {
        self.cross_second_deriv_cps_uu_u0()
    }
    fn cross_second_deriv_cps_uu_u1_degree9(&self) -> [Vec3; 10] {
        self.cross_second_deriv_cps_uu_u1()
    }
    fn cross_second_deriv_cps_vv_v0_degree9(&self) -> [Vec3; 10] {
        self.cross_second_deriv_cps_vv_v0()
    }
    fn cross_second_deriv_cps_vv_v1_degree9(&self) -> [Vec3; 10] {
        self.cross_second_deriv_cps_vv_v1()
    }
}

impl KpGuide for KpGuidePatchBiQuintic {
    fn point_at(&self, u: f64, v: f64) -> Point3 {
        KpGuidePatchBiQuintic::point_at(self, u, v)
    }
    fn partial_u(&self, u: f64, v: f64) -> Vec3 {
        KpGuidePatchBiQuintic::partial_u(self, u, v)
    }
    fn partial_v(&self, u: f64, v: f64) -> Vec3 {
        KpGuidePatchBiQuintic::partial_v(self, u, v)
    }
    fn second_partial_uu(&self, u: f64, v: f64) -> Vec3 {
        KpGuidePatchBiQuintic::second_partial_uu(self, u, v)
    }
    fn second_partial_uv(&self, u: f64, v: f64) -> Vec3 {
        KpGuidePatchBiQuintic::second_partial_uv(self, u, v)
    }
    fn second_partial_vv(&self, u: f64, v: f64) -> Vec3 {
        KpGuidePatchBiQuintic::second_partial_vv(self, u, v)
    }
    fn cross_deriv_cps_u0_degree9(&self) -> [Vec3; 10] {
        elevate_quintic_vec_to_nine(self.cross_deriv_cps_u0())
    }
    fn cross_deriv_cps_u1_degree9(&self) -> [Vec3; 10] {
        elevate_quintic_vec_to_nine(self.cross_deriv_cps_u1())
    }
    fn cross_deriv_cps_v0_degree9(&self) -> [Vec3; 10] {
        elevate_quintic_vec_to_nine(self.cross_deriv_cps_v0())
    }
    fn cross_deriv_cps_v1_degree9(&self) -> [Vec3; 10] {
        elevate_quintic_vec_to_nine(self.cross_deriv_cps_v1())
    }
    fn cross_second_deriv_cps_uu_u0_degree9(&self) -> [Vec3; 10] {
        elevate_quintic_vec_to_nine(self.cross_second_deriv_cps_uu_u0())
    }
    fn cross_second_deriv_cps_uu_u1_degree9(&self) -> [Vec3; 10] {
        elevate_quintic_vec_to_nine(self.cross_second_deriv_cps_uu_u1())
    }
    fn cross_second_deriv_cps_vv_v0_degree9(&self) -> [Vec3; 10] {
        elevate_quintic_vec_to_nine(self.cross_second_deriv_cps_vv_v0())
    }
    fn cross_second_deriv_cps_vv_v1_degree9(&self) -> [Vec3; 10] {
        elevate_quintic_vec_to_nine(self.cross_second_deriv_cps_vv_v1())
    }
}

#[cfg(test)]
#[path = "kp_guide_trait_tests.rs"]
mod tests;
