//! §5.3 Round 6/7 rational-Gregory construction path.
//!
//! `try_rational_gregory_blend` builds a true rational Gregory patch
//! (Longhi triangular for N=3, Chiyokura-Kimura rectangular for N=4
//! — follow-up) from a [`BlendBoundary`] plus the incident fillet
//! strip data.
//!
//! **Round 7 (this commit)** landed [`extract_triangular_inner_cps`],
//! which replaces the earlier `[v_pos; 3]` placeholder by computing
//! each Longhi inner candidate from the adjacent fillet strip's
//! surface normal at the boundary midpoint — a geometrically-
//! grounded placement that actually exercises the rational-inner-CP
//! degree of freedom. The dispatch in
//! [`super::surface::build_smooth_blend`] is still gated pending
//! curvature-verified regression tests (Round 7 C4b) — once the
//! builder's output is confirmed to match or beat the polar-NURBS
//! fallback on both uniform-radius and mixed-radius cube corners,
//! the dispatcher will call it ahead of the polar NURBS.
//!
//! Kept in its own file so the surface dispatcher stays under the
//! 500-line structural cap and so the dispatch-activation follow-up
//! has a clean module boundary to plug into.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::gregory::{GregoryPatch, TriangularGregory};
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::uv_projection::project_point_to_uv;
use popkern_core::geometry::{NurbsCurve, Point3, Vec3};

use super::boundary::{BlendBoundary, BoundarySegment};
use super::surface::BlendSurfaceResult;
use crate::fillet_surface::FilletSurfaceResult;
use crate::types::{ContinuityLevel, CornerMethod};

/// §5.3 Round 6 — build a true rational Gregory patch for a
/// vertex-blend corner with `N ∈ {3, 4}` analytic-strip fillets of
/// uniform radius.
///
/// Currently handles **N=3 triangular Longhi Gregory** only; N=4 is
/// a follow-up commit once a 4-fillet test fixture exists. For N=3:
///
/// 1. Verify all incident fillet surfaces are `Torus` or `Cylinder`
///    and share one minor radius.
/// 2. Extract the 3 boundary cubic Beziers by fitting each
///    `BoundarySegment::FilletEnd` curve (assumed nearly-cubic on
///    the fillet-end arc scale) to a cubic Bezier control polygon.
/// 3. Place the 3 inner candidates at the vertex position — the
///    rational blend then passes through the vertex at the patch
///    interior while the 3 boundary curves stitch through the
///    incident fillet strips via the cross-derivative encoding in
///    the second Bezier control point.
/// 4. Return the patch as [`SurfaceGeometry::Gregory`].
///
/// **Activation caveat:** the dispatcher does not yet call this
/// function — see the module-level docs for the twist-vector
/// prerequisite.
#[allow(dead_code)]
pub(crate) fn try_rational_gregory_blend(
    boundary: &BlendBoundary,
    fillet_surfaces: &[FilletSurfaceResult],
    boundary_curves: &[NurbsCurve],
) -> Option<BlendSurfaceResult> {
    let n_fillets = boundary
        .segments
        .iter()
        .filter(|s| matches!(s, BoundarySegment::FilletEnd { .. }))
        .count();

    if n_fillets != 3 {
        return None;
    }

    // Verify uniform radii across analytic strips.
    let radii: Vec<f64> = fillet_surfaces
        .iter()
        .filter_map(|fr| match &fr.surface {
            SurfaceGeometry::Torus { minor_radius, .. } => Some(*minor_radius),
            SurfaceGeometry::Cylinder { radius, .. } => Some(*radius),
            _ => None,
        })
        .collect();
    if radii.len() != fillet_surfaces.len() || radii.is_empty() {
        return None;
    }
    let (lo, hi) = radii.iter().fold((f64::INFINITY, 0.0_f64), |(lo, hi), &r| {
        (lo.min(r), hi.max(r))
    });
    if hi - lo > 1e-6 {
        return None;
    }

    let beziers: Vec<[Point3; 4]> = boundary_curves
        .iter()
        .map(fit_cubic_bezier_through)
        .collect();
    if beziers.len() != 3 {
        return None;
    }

    // §5.3 Round 7 — twist-vector extraction. Each Longhi inner
    // candidate `P_k` controls the interior cross-derivative along
    // the boundary opposite corner `k` (at the midpoint of that
    // boundary the rational blend reduces to `P_k` alone). Place
    // each `P_k` by starting at the midpoint of its associated
    // boundary and displacing **inward** along the incident fillet
    // strip's outward surface normal at that midpoint. Step size
    // is one average radius — the leading-order curvature match
    // for a cubic Bezier in the radial direction.
    let avg_radius = radii.iter().sum::<f64>() / radii.len() as f64;
    let inner_cps = extract_triangular_inner_cps(
        boundary_curves,
        fillet_surfaces,
        avg_radius,
        boundary.vertex_position,
    )?;
    let bez_array: [[Point3; 4]; 3] = [beziers[0], beziers[1], beziers[2]];

    let patch = TriangularGregory::from_boundary_beziers(bez_array, inner_cps, 1e-3).ok()?;

    Some(BlendSurfaceResult {
        surface: SurfaceGeometry::Gregory(GregoryPatch::Triangular(Box::new(patch))),
        boundary_curves: boundary_curves.to_vec(),
        method: CornerMethod::Gregory,
        continuity: ContinuityLevel::G2,
    })
}

/// Fit a cubic Bezier control polygon to a `NurbsCurve` by
/// interpolating at `t = 0, 1/3, 2/3, 1` and solving for the inner
/// 2 control points so the cubic Bezier passes exactly through the
/// interior sample points.
///
/// The system is `[12 6; 6 12] [P_1; P_2] = [27 Q_1 - 8 P_0 - P_3;
///                                           27 Q_2 - P_0 - 8 P_3]`
/// with determinant `108`, yielding `P_1 = (2A - B) / 18` and
/// `P_2 = (2B - A) / 18` where `A` and `B` are the right-hand sides.
#[allow(dead_code)]
pub(crate) fn fit_cubic_bezier_through(curve: &NurbsCurve) -> [Point3; 4] {
    let (t0, t1) = curve.parameter_range();
    let span = t1 - t0;
    let p0 = curve.point_at(t0);
    let p3 = curve.point_at(t1);
    let q1 = curve.point_at(t0 + span / 3.0);
    let q2 = curve.point_at(t0 + 2.0 * span / 3.0);

    let axis = |q1c: f64, q2c: f64, p0c: f64, p3c: f64| -> (f64, f64) {
        let a = 27.0 * q1c - 8.0 * p0c - p3c;
        let b = 27.0 * q2c - p0c - 8.0 * p3c;
        ((2.0 * a - b) / 18.0, (2.0 * b - a) / 18.0)
    };
    let (p1x, p2x) = axis(q1.x, q2.x, p0.x, p3.x);
    let (p1y, p2y) = axis(q1.y, q2.y, p0.y, p3.y);
    let (p1z, p2z) = axis(q1.z, q2.z, p0.z, p3.z);

    [
        p0,
        Point3::new(p1x, p1y, p1z),
        Point3::new(p2x, p2y, p2z),
        p3,
    ]
}

/// Compute the 3 Longhi inner-CP candidates for a triangular
/// Gregory vertex blend from the 3 boundary Bezier arcs and the
/// incident fillet-strip data.
///
/// The `k`-th candidate corresponds to the boundary opposite
/// corner `k`; at the midpoint of that boundary the Longhi
/// rational blend reduces to `P_k` alone, so positioning `P_k`
/// directly controls the patch's interior behaviour across that
/// edge.
///
/// **Placement rule.** Let `M_k` be the 3D midpoint of boundary
/// `k` (evaluated at the mid-parameter) and `n̂_k` the outward
/// surface normal of the incident fillet strip at that midpoint.
/// The candidate is placed one step inward:
///
/// `P_k = M_k − step_k · n̂_k`
///
/// where `step_k = avg_radius` (one local radius of curvature) is
/// the leading-order match for a cubic Bezier in the radial
/// direction. The sign of the displacement is disambiguated by
/// checking that the resulting `P_k` moves **toward** the vertex
/// position (the patch interior) rather than away from it; if the
/// projected normal points the wrong way (e.g., the strip is
/// parameterised with reversed orientation) we flip the step.
///
/// Returns `None` if UV projection fails on any incident strip —
/// the caller falls through to the polar-NURBS approximation.
pub(crate) fn extract_triangular_inner_cps(
    boundary_curves: &[NurbsCurve],
    fillet_surfaces: &[FilletSurfaceResult],
    avg_radius: f64,
    vertex_position: Point3,
) -> Option<[Point3; 3]> {
    let mut inner = [Point3::new(0.0, 0.0, 0.0); 3];
    for (k, slot) in inner.iter_mut().enumerate() {
        let curve = boundary_curves.get(k)?;
        let strip = fillet_surfaces.get(k)?;
        let (t0, t1) = curve.parameter_range();
        let tmid = 0.5 * (t0 + t1);
        let mid = curve.point_at(tmid);

        // Strip surface normal at midpoint.
        let (u, v) = project_point_to_uv(&strip.surface, &mid)?;
        let n_raw = strip.surface.normal_at(u, v);
        let n_hat = n_raw.normalized()?;

        // Disambiguate sign — we want "inward from the strip" =
        // toward the vertex position (a reliable proxy for the
        // patch interior on the common convex-corner case).
        //
        // Desired: `P_k = M_k + step · d̂` where `d̂` points from
        // `M_k` toward the patch interior. Convention: if `n̂`
        // already points toward the vertex (dot with `to_vertex`
        // positive), use `d̂ = n̂`; otherwise flip to `−n̂`.
        let to_vertex = Vec3::new(
            vertex_position.x - mid.x,
            vertex_position.y - mid.y,
            vertex_position.z - mid.z,
        );
        let sign = if n_hat.dot(&to_vertex) >= 0.0 {
            1.0
        } else {
            -1.0
        };
        let step = sign * avg_radius;

        *slot = Point3::new(
            mid.x + step * n_hat.x,
            mid.y + step * n_hat.y,
            mid.z + step * n_hat.z,
        );
    }
    Some(inner)
}
