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
use popkern_math::PortableMath;

pub(crate) use super::dispatch_nested::try_levin_pair_nested;

/// Symbolic-first tangent-intersection path. When the pencil
/// `λA + B` contains a rank-1 fibre `D = c · n · nᵀ`, the surfaces
/// touch tangentially along the plane `nᵀp = 0`; the intersection
/// is the conic `surf_a ∩ plane = surf_b ∩ plane`. We recover `n`
/// exactly in `Q` (no FP fudge) via `tangent::surfaces_are_tangent`,
/// build a `Plane` from it, and reuse the plane×surf dispatchers
/// to produce one `EdgeTriple` whose recorded `surf_a`/`surf_b`
/// remain the original quadrics.
///
/// Stage 2b.41: closes the cyl-with-inscribed-sphere class
/// (canonical CAD tangent case) which the f64 Levin path didn't
/// exercise because the pair has no `try_cyl_sphere_*` analytic
/// fallback. Returns `None` for non-rational input, no rational
/// rank-1 fibre, or anchors not on the recovered plane.
pub(crate) fn try_tangent_pair(
    surf_a: &SurfaceGeometry,
    surf_b: &SurfaceGeometry,
    anchors: &EdgeAnchors,
    tol: &Tolerance,
) -> Option<Result<EdgeTriple, EdgeTripleError>> {
    use num_traits::ToPrimitive;
    let plane_rat = super::tangent::surfaces_are_tangent(surf_a, surf_b)?;
    let mut plane: [f64; 4] = [
        plane_rat[0].to_f64()?,
        plane_rat[1].to_f64()?,
        plane_rat[2].to_f64()?,
        plane_rat[3].to_f64()?,
    ];
    let n_len2 = plane[0] * plane[0] + plane[1] * plane[1] + plane[2] * plane[2];
    if n_len2 < 1e-30 {
        return None;
    }
    let n_len = n_len2.sqrt();
    for v in plane.iter_mut() {
        *v /= n_len;
    }
    let plane_surf = super::decompose::plane_to_surface(plane);

    let triple_a = intersect_plane_with_target(&plane_surf, surf_a, anchors, tol)?.ok()?;
    let triple_b = intersect_plane_with_target(&plane_surf, surf_b, anchors, tol)?.ok()?;

    let cardinality = Some(crate::geometry::edge_triple::types::BranchCardinality {
        total: 1,
        kinds: vec![crate::geometry::edge_triple::types::BranchKind::from_curve(
            &triple_a.curve_3d,
        )],
        // Stage 2b.42: this path is reached only via the symbolic
        // rank-1 fibre detection in `tangent::find_tangent_plane`,
        // so the intersection is tangent by construction.
        tangent: true,
    });
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
            cardinality,
        },
        quality: EdgeTripleQuality::default(),
    };
    Some(super::super::analytic::finalize_quality_pub(
        triple, surf_a, surf_b, tol,
    ))
}

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

    // Cardinality: total branches = #planes the pencil decomposed
    // into. We classify each branch's conic type by sampling the
    // branch curve via plane × surf_a (just for type identification —
    // a single sample is enough).
    let total_branches = planes.len() as u8;
    let mut branch_kinds = Vec::with_capacity(planes.len());
    for p in &planes {
        if let Some(kind) = classify_branch_kind(p, surf_a, anchors, tol) {
            branch_kinds.push(kind);
        }
    }
    let cardinality = Some(crate::geometry::edge_triple::types::BranchCardinality {
        total: total_branches,
        kinds: branch_kinds,
        // Levin's rank-2 plane decomposition handles transverse
        // intersections (two distinct real planes); rank-1 tangent
        // fibres are routed to `try_tangent_pair` upstream of this
        // function, so reaching here implies a transverse pair.
        tangent: false,
    });

    // The composed triple uses triple_a's curve_3d as the canonical
    // parameterisation. triple_a's pcurve_a is parameterised
    // consistently with that; triple_b's is parameterised by triple_b's
    // own curve_3d. Both trace the same physical curve in the shared
    // plane, so the parameter relation is `t_b = a + b · t_a` with
    // `b ∈ {+1, −1}`. derive_linear_affine_b recovers (a, b).
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
            cardinality,
        },
        quality: EdgeTripleQuality::default(),
    };
    Some(super::super::analytic::finalize_quality_pub(
        triple, surf_a, surf_b, tol,
    ))
}

pub(super) fn intersect_plane_with_target(
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
/// Both `triple_a.curve_3d` and `triple_b.curve_3d` trace the same
/// physical curve in the shared plane (Arc / Ellipse / Parabola /
/// Hyperbola depending on the surface pair). When the two curves use
/// the *same* parameterisation type, the relation is `t_b = a + b·t_a`
/// with `b ∈ {+1, −1}`. We recover `(a, b)` from the start anchor
/// (where both `t_start`s map to the same physical point) plus a
/// mid-curve sample to disambiguate the slope sign.
pub(super) fn derive_linear_affine_b(
    triple_a: &EdgeTriple,
    triple_b: &EdgeTriple,
) -> Option<(f64, f64)> {
    let alpha_a_start = triple_a.t_start;
    let alpha_a_mid = (triple_a.t_start + triple_a.t_end) * 0.5;
    let p_start = triple_a.curve_3d.point_at(alpha_a_start);
    let p_mid = triple_a.curve_3d.point_at(alpha_a_mid);

    let beta_start = invert_curve_param(&triple_b.curve_3d, p_start)?;
    let beta_mid = invert_curve_param(&triple_b.curve_3d, p_mid)?;

    let delta_a = alpha_a_mid - alpha_a_start;
    if delta_a.abs() < 1e-15 {
        return None;
    }
    // For Ellipse / Arc, the parameter is angular: wrap delta_b into
    // (-π, π]. Hyperbola / Parabola use signed real parameters and
    // don't need wrapping.
    let needs_angle_wrap = matches!(
        triple_b.curve_3d,
        crate::geometry::curve::CurveGeometry::Ellipse { .. }
            | crate::geometry::curve::CurveGeometry::Arc { .. }
    );
    let mut delta_b = beta_mid - beta_start;
    if needs_angle_wrap {
        while delta_b > std::f64::consts::PI {
            delta_b -= std::f64::consts::TAU;
        }
        while delta_b < -std::f64::consts::PI {
            delta_b += std::f64::consts::TAU;
        }
    }
    let b_b = delta_b / delta_a;
    // Snap to ±1 (the only true values for parameterisation-preserving
    // change of basis on the same shared curve).
    let b_b_snap = if b_b > 0.0 { 1.0 } else { -1.0 };
    let a_b = beta_start - b_b_snap * alpha_a_start;
    Some((a_b, b_b_snap))
}

/// Invert a 3-D curve parameterisation: given a point `p` on the
/// curve, find the parameter value `t` such that `curve.point_at(t) =
/// p`. Supports Arc (circle), Ellipse, Hyperbola, and Parabola.
fn invert_curve_param(
    curve: &crate::geometry::curve::CurveGeometry,
    p: crate::geometry::Point3,
) -> Option<f64> {
    use crate::geometry::curve::CurveGeometry;
    match curve {
        CurveGeometry::Arc { center, normal, .. } => {
            let d = p - *center;
            let (u_dir, v_dir) = crate::geometry::analytic_curves::circle_basis(*normal);
            Some(d.dot(&v_dir).patan2(d.dot(&u_dir)))
        }
        CurveGeometry::Ellipse {
            center,
            u_basis,
            v_basis,
            ..
        } => {
            // Solve d = α·u + β·v with α=cos t, β=sin t.
            let (alpha, beta) = solve_2basis(*center, *u_basis, *v_basis, p)?;
            Some(beta.patan2(alpha))
        }
        CurveGeometry::Hyperbola {
            center,
            semi_major,
            semi_minor,
            ..
        } => {
            // Solve d = α·u + β·v with α=cosh t, β=sinh t.
            let (alpha, beta) = solve_2basis(*center, *semi_major, *semi_minor, p)?;
            // t = asinh(β) (sign carries the branch); cosh(asinh β) = √(1+β²)
            // matches α up to sign, indicating which nappe.
            Some(beta.pasinh().copysign(alpha))
        }
        CurveGeometry::Parabola {
            vertex,
            linear,
            quadratic,
            ..
        } => {
            // Solve d = t·linear + t²·quadratic. The 2-basis solve
            // gives (a, b) = (t, t²). We use a as t (the linear part
            // is more numerically stable when t is small).
            let (a, _b) = solve_2basis(*vertex, *linear, *quadratic, p)?;
            Some(a)
        }
        _ => None,
    }
}

/// Decompose `p − origin` against the (possibly non-orthogonal)
/// basis `(u, v)` spanning the curve's plane. Returns `(α, β)` such
/// that `p − origin ≈ α·u + β·v`, assuming `p` lies in the plane.
fn solve_2basis(
    origin: crate::geometry::Point3,
    u: crate::geometry::Vec3,
    v: crate::geometry::Vec3,
    p: crate::geometry::Point3,
) -> Option<(f64, f64)> {
    let d = p - origin;
    let uu = u.dot(&u);
    let uv = u.dot(&v);
    let vv = v.dot(&v);
    let det = uu * vv - uv * uv;
    if det.abs() < 1e-15 {
        return None;
    }
    let ud = u.dot(&d);
    let vd = v.dot(&d);
    let alpha = (vv * ud - uv * vd) / det;
    let beta = (uu * vd - uv * ud) / det;
    Some((alpha, beta))
}

/// Identify the conic kind of one Levin branch (plane × surf_a)
/// without relying on anchors lying on this specific branch — we
/// just need the curve type. We probe each plane×{cyl,sphere,cone}
/// dispatcher with a tolerance-permissive call and read the
/// resulting `curve_3d` variant. Returns `None` if no dispatcher
/// applies to this plane × surface combination (e.g., the plane
/// doesn't actually intersect the surface for this branch — likely
/// the "plane at infinity" stub from rank-2 mixed-sign decomposition).
pub(super) fn classify_branch_kind(
    plane: &[f64; 4],
    target: &SurfaceGeometry,
    anchors: &EdgeAnchors,
    tol: &Tolerance,
) -> Option<crate::geometry::edge_triple::types::BranchKind> {
    let plane_surf = plane_to_surface(*plane);
    // Skip the "plane at infinity" stub (zero spatial normal).
    let SurfaceGeometry::Plane { normal, .. } = &plane_surf else {
        return None;
    };
    if normal.length() < 1e-10 {
        return None;
    }
    let triple = intersect_plane_with_target(&plane_surf, target, anchors, tol)?.ok()?;
    Some(crate::geometry::edge_triple::types::BranchKind::from_curve(
        &triple.curve_3d,
    ))
}

pub(super) fn matrix_max_abs(m: &[[f64; 4]; 4]) -> f64 {
    let mut max = 0.0_f64;
    for row in m.iter() {
        for &v in row {
            max = max.max(v.abs());
        }
    }
    max
}

pub(super) fn provenance_kind(a: &SurfaceGeometry, b: &SurfaceGeometry) -> AnalyticPair {
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
