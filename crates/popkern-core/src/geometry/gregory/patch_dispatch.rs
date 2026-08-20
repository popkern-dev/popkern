//! [`GregoryPatch`] — the variant dispatcher over every Gregory patch
//! form the kernel builds (rectangular degrees 3/5/7/9, two-stripe,
//! dual-twist, triangular), with unified evaluation, derivatives, and
//! self-intersection entry points.

use super::*;

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
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
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
    /// Bi-nonic with Chiyokura-Kimura dual-twist corners (see
    /// `docs/decisions/2026-06-11-kp20-binonic-dual-twist-corners.md`).
    RectangularBiNonicDualTwist(Box<RectangularGregoryBiNonicDualTwist>),
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
            GregoryPatch::RectangularBiNonicDualTwist(r) => r.point_at(u, v),
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
            GregoryPatch::RectangularBiNonicDualTwist(r) => r.partial_u(u, v),
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
            GregoryPatch::RectangularBiNonicDualTwist(r) => r.partial_v(u, v),
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

    /// Return the measured terms behind
    /// [`self_intersects_with_tolerance`](Self::self_intersects_with_tolerance).
    pub fn self_intersection_report(&self, tolerance: f64) -> GregorySelfIntersectionReport {
        self_intersection::gregory_self_intersection_report(self, tolerance)
    }

    /// Return only the fast sampled-Jacobian self-intersection terms.
    ///
    /// This is a conservative screening primitive for optimizers that need
    /// cheap fold ordering before running full subdivision validation on the
    /// candidate they intend to accept.
    pub fn jacobian_self_intersection_report(&self) -> GregorySelfIntersectionReport {
        self_intersection::gregory_jacobian_self_intersection_report(self)
    }
}
