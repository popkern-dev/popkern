//! Wire the Levin pipeline (`quadric → pencil → degenerate-fiber →
//! plane decomposition`) into the EdgeTriple dispatch.
//!
//! Strategy:
//! 1. Build quadric forms `A`, `B` via `surface_to_quadric`.
//! 2. Find a degenerate fibre `D` in the pencil `λA + μB`.
//! 3. Decompose `D` into one or two planes (rank ≤ 2 case).
//! 4. Pick the plane that contains both anchors.
//! 5. Reuse the existing `plane × {cyl, sphere, cone}` analytic
//!    dispatchers — once for `surf_a`, once for `surf_b` — and stitch
//!    the resulting `curve_3d` + per-side pcurves into one
//!    `EdgeTriple` whose `surf_a` / `surf_b` remain the user's
//!    original quadrics.
//!
//! Closes the commercial-kernel quadric-pair gap for **rank ≤ 2** fibres,
//! which empirically covers cylinder×cylinder skew, cone×cone
//! non-coaxial, and cylinder×cone skew. Rank-3 (proper cone fibre)
//! is deferred to stage 2b.35.

use super::decompose::{DegenerateShape, classify_degenerate, plane_to_surface};
use super::pencil::find_degenerate_fiber;
use super::quadric::surface_to_quadric;
use crate::geometry::Point3;
use crate::geometry::curve::CurveEval;
use crate::geometry::edge_triple::analytic_trig_cone::try_plane_cone_perpendicular;
use crate::geometry::edge_triple::analytic_trig_plane_cone_hyperbola::try_plane_cone_oblique_hyperbola;
use crate::geometry::edge_triple::analytic_trig_plane_cone_oblique::try_plane_cone_oblique_ellipse;
use crate::geometry::edge_triple::analytic_trig_plane_cone_parabola::try_plane_cone_oblique_parabola;
use crate::geometry::edge_triple::analytic_trig_plane_cyl::try_plane_cylinder_perpendicular;
use crate::geometry::edge_triple::analytic_trig_plane_cyl_oblique::try_plane_cylinder_oblique;
use crate::geometry::edge_triple::analytic_trig_plane_sphere::try_plane_sphere_polar;
use crate::geometry::edge_triple::analytic_trig_plane_sphere_general::try_plane_sphere_general;
use crate::geometry::edge_triple::error::EdgeTripleError;
use crate::geometry::edge_triple::inputs::EdgeAnchors;
use crate::geometry::edge_triple::types::{
    AnalyticPair, BranchSelector, EdgeTriple, EdgeTripleProvenance, EdgeTripleQuality, ParameterMap,
};
use crate::geometry::surface::SurfaceGeometry;
use crate::tolerance::Tolerance;

/// Try to compute the intersection of two analytic quadric surfaces
/// via Levin's algorithm. Returns `None` if Levin doesn't apply
/// (non-quadric surfaces, no real degenerate fibre, rank ≥ 3 fibre,
/// no plane through both anchors, etc.).
pub(crate) fn try_levin_pair(
    surf_a: &SurfaceGeometry,
    surf_b: &SurfaceGeometry,
    anchors: &EdgeAnchors,
    tol: &Tolerance,
) -> Option<Result<EdgeTriple, EdgeTripleError>> {
    let a_q = surface_to_quadric(surf_a)?;
    let b_q = surface_to_quadric(surf_b)?;
    let fiber = find_degenerate_fiber(&a_q, &b_q)?;
    let scale = matrix_max_abs(&a_q.0).max(matrix_max_abs(&b_q.0)).max(1.0);
    let shape = classify_degenerate(&fiber.quadric, scale);
    let planes: Vec<[f64; 4]> = match shape {
        DegenerateShape::SinglePlane(p) => vec![p],
        DegenerateShape::TwoPlanes(p1, p2) => vec![p1, p2],
        _ => return None,
    };

    let (start, end) = match anchors {
        EdgeAnchors::Open { start, end } => (*start, *end),
        _ => return None,
    };
    let plane_distance =
        |p: &[f64; 4], q: Point3| -> f64 { (p[0] * q.x + p[1] * q.y + p[2] * q.z + p[3]).abs() };
    let anchor_tol = 100.0 * tol.linear();
    let plane = *planes
        .iter()
        .find(|p| plane_distance(p, start) <= anchor_tol && plane_distance(p, end) <= anchor_tol)?;
    let plane_surf = plane_to_surface(plane);

    // Reuse plane × {cyl, sphere, cone} dispatchers for both surfaces.
    let triple_a = intersect_plane_with_target(&plane_surf, surf_a, anchors, tol)?.ok()?;
    let triple_b = intersect_plane_with_target(&plane_surf, surf_b, anchors, tol)?.ok()?;

    // The composed triple uses triple_a's curve_3d as the canonical
    // parameterisation. triple_a's pcurve_a is parameterised
    // consistently with that. triple_b's pcurve_a is parameterised
    // by triple_b's *own* curve_3d (which traces the same physical
    // ellipse but with a different basis), so we need a parameter
    // map t_a → t_b. Both ellipses lie in the shared plane and
    // trace the same point set; the map is `α_b = a + b · α_a`
    // with `b ∈ {+1, −1}`. We recover (a, b) from the start anchor
    // and one mid-curve sample.
    let (a_b_coef, b_b_coef) = derive_linear_affine_b(&triple_a, &triple_b)?;

    let triple = EdgeTriple {
        curve_3d: triple_a.curve_3d.clone(),
        pcurve_a: triple_a.pcurve_a.clone(),
        pcurve_b: triple_b.pcurve_a.clone(),
        t_start: triple_a.t_start,
        t_end: triple_a.t_end,
        topology: triple_a.topology,
        parameter_map: ParameterMap::LinearAffine {
            a_a: 0.0,
            b_a: 1.0,
            a_b: a_b_coef,
            b_b: b_b_coef,
        },
        provenance: EdgeTripleProvenance::AnalyticIntersection {
            kind: provenance_kind(surf_a, surf_b),
            branch: BranchSelector::OnlyBranch,
        },
        quality: EdgeTripleQuality::default(),
    };
    Some(super::super::analytic::finalize_quality_pub(
        triple, surf_a, surf_b, tol,
    ))
}

fn intersect_plane_with_target(
    plane_surf: &SurfaceGeometry,
    target: &SurfaceGeometry,
    anchors: &EdgeAnchors,
    tol: &Tolerance,
) -> Option<Result<EdgeTriple, EdgeTripleError>> {
    let (p_o, p_n, p_u, p_v) = match plane_surf {
        SurfaceGeometry::Plane {
            origin,
            normal,
            u_axis,
            v_axis,
        } => (*origin, *normal, *u_axis, *v_axis),
        _ => return None,
    };
    match target {
        SurfaceGeometry::Cylinder {
            origin,
            axis,
            radius,
        } => try_plane_cylinder_perpendicular(
            p_o, p_n, p_u, p_v, true, *origin, *axis, *radius, anchors, tol,
        )
        .or_else(|| {
            try_plane_cylinder_oblique(
                p_o, p_n, p_u, p_v, true, *origin, *axis, *radius, anchors, tol,
            )
        }),
        SurfaceGeometry::Sphere { center, radius } => {
            try_plane_sphere_polar(p_o, p_n, p_u, p_v, true, *center, *radius, anchors, tol)
                .or_else(|| {
                    try_plane_sphere_general(
                        p_o, p_n, p_u, p_v, true, *center, *radius, anchors, tol,
                    )
                })
        }
        SurfaceGeometry::Cone {
            apex,
            axis,
            half_angle,
        } => try_plane_cone_perpendicular(
            p_o,
            p_n,
            p_u,
            p_v,
            true,
            *apex,
            *axis,
            *half_angle,
            anchors,
            tol,
        )
        .or_else(|| {
            try_plane_cone_oblique_ellipse(
                p_o,
                p_n,
                p_u,
                p_v,
                true,
                *apex,
                *axis,
                *half_angle,
                anchors,
                tol,
            )
        })
        .or_else(|| {
            try_plane_cone_oblique_parabola(
                p_o,
                p_n,
                p_u,
                p_v,
                true,
                *apex,
                *axis,
                *half_angle,
                anchors,
                tol,
            )
        })
        .or_else(|| {
            try_plane_cone_oblique_hyperbola(
                p_o,
                p_n,
                p_u,
                p_v,
                true,
                *apex,
                *axis,
                *half_angle,
                anchors,
                tol,
            )
        }),
        _ => None,
    }
}

/// Find `(a_b, b_b)` such that `t_b = a_b + b_b · t_a` makes
/// `triple_b.pcurve_a(t_b)` evaluate to the same physical point as
/// `triple_a.curve_3d(t_a)` for all `t_a` in `[triple_a.t_start,
/// triple_a.t_end]`.
///
/// Both `triple_a.curve_3d` and `triple_b.curve_3d` are ellipses in
/// the same plane tracing the same point set. The map `α_b = α_b(α_a)`
/// is linear with slope `±1`: rotation when `b = +1`, reflection when
/// `b = -1`. We recover `(a, b)` from the start anchor (which both
/// ellipses pass through at their respective `t_start`) plus a
/// mid-curve sample to disambiguate the slope sign.
fn derive_linear_affine_b(
    triple_a: &EdgeTriple,
    triple_b: &EdgeTriple,
) -> Option<(f64, f64)> {
    use crate::geometry::curve::CurveGeometry;
    let CurveGeometry::Ellipse {
        center: cb,
        u_basis: ub,
        v_basis: vb,
        ..
    } = &triple_b.curve_3d
    else {
        return None;
    };
    // Sample triple_a's curve at start and midpoint.
    let alpha_a_start = triple_a.t_start;
    let alpha_a_mid = (triple_a.t_start + triple_a.t_end) * 0.5;
    let p_start = triple_a.curve_3d.point_at(alpha_a_start);
    let p_mid = triple_a.curve_3d.point_at(alpha_a_mid);

    // Invert each onto the curve_3d_b parameterisation by solving
    // `p − cb = cos β · u_b + sin β · v_b` for `β`. This is a 2×2
    // linear system in (cos β, sin β).
    let beta_start = invert_ellipse_param(*cb, *ub, *vb, p_start)?;
    let beta_mid = invert_ellipse_param(*cb, *ub, *vb, p_mid)?;

    // Wrap delta_b into (-π, π].
    let mut delta_b = beta_mid - beta_start;
    while delta_b > std::f64::consts::PI {
        delta_b -= std::f64::consts::TAU;
    }
    while delta_b < -std::f64::consts::PI {
        delta_b += std::f64::consts::TAU;
    }
    let delta_a = alpha_a_mid - alpha_a_start;
    if delta_a.abs() < 1e-15 {
        return None;
    }
    let b_b = delta_b / delta_a;
    // Snap to ±1 (the only true values for parameterisation-preserving
    // change of basis on a planar ellipse).
    let b_b_snap = if b_b > 0.0 { 1.0 } else { -1.0 };
    let a_b = beta_start - b_b_snap * alpha_a_start;
    Some((a_b, b_b_snap))
}

/// Invert `p − center = cos β · u + sin β · v` for `β`, where `u`
/// and `v` are the ellipse's (possibly non-orthogonal) basis vectors
/// spanning the plane in which the ellipse lies. Returns `None` if
/// `(u, v)` are degenerate.
fn invert_ellipse_param(
    center: crate::geometry::Point3,
    u: crate::geometry::Vec3,
    v: crate::geometry::Vec3,
    p: crate::geometry::Point3,
) -> Option<f64> {
    // Solve [[u·u, u·v], [u·v, v·v]] · (cos β, sin β) = (u·d, v·d).
    let d = p - center;
    let uu = u.dot(&u);
    let uv = u.dot(&v);
    let vv = v.dot(&v);
    let det = uu * vv - uv * uv;
    if det.abs() < 1e-15 {
        return None;
    }
    let ud = u.dot(&d);
    let vd = v.dot(&d);
    let cos_b = (vv * ud - uv * vd) / det;
    let sin_b = (uu * vd - uv * ud) / det;
    Some(sin_b.atan2(cos_b))
}

fn matrix_max_abs(m: &[[f64; 4]; 4]) -> f64 {
    let mut max = 0.0_f64;
    for row in m.iter() {
        for &v in row {
            max = max.max(v.abs());
        }
    }
    max
}

fn provenance_kind(a: &SurfaceGeometry, b: &SurfaceGeometry) -> AnalyticPair {
    use SurfaceGeometry::*;
    match (a, b) {
        (Cylinder { .. }, Cylinder { .. }) => AnalyticPair::CylinderCylinder,
        (Cone { .. }, Cone { .. }) => AnalyticPair::ConeCone,
        (Cylinder { .. }, Cone { .. }) | (Cone { .. }, Cylinder { .. }) => {
            AnalyticPair::CylinderCone
        }
        (Cylinder { .. }, Sphere { .. }) | (Sphere { .. }, Cylinder { .. }) => {
            AnalyticPair::CylinderSphere
        }
        (Cone { .. }, Sphere { .. }) | (Sphere { .. }, Cone { .. }) => AnalyticPair::SphereCone,
        (Sphere { .. }, Sphere { .. }) => AnalyticPair::SphereSphere,
        (Plane { .. }, Cylinder { .. }) | (Cylinder { .. }, Plane { .. }) => {
            AnalyticPair::PlaneCylinder
        }
        (Plane { .. }, Sphere { .. }) | (Sphere { .. }, Plane { .. }) => AnalyticPair::PlaneSphere,
        (Plane { .. }, Cone { .. }) | (Cone { .. }, Plane { .. }) => AnalyticPair::PlaneCone,
        (Plane { .. }, Plane { .. }) => AnalyticPair::PlanePlane,
        _ => AnalyticPair::CylinderCylinder, // shouldn't reach for non-quadric inputs
    }
}
