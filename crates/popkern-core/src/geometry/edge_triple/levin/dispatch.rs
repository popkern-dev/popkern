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
    // Note: these go through `decompose_2x2` for the elliptical pcurve
    // construction, which has a known bug when the projection matrix
    // has m22 > m11 or negative determinant. End-to-end Levin tests
    // are gated on a `decompose_2x2` fix (see stage 2b.34a/35 plan).
    let triple_a = intersect_plane_with_target(&plane_surf, surf_a, anchors, tol)?.ok()?;
    let triple_b = intersect_plane_with_target(&plane_surf, surf_b, anchors, tol)?.ok()?;

    let triple = EdgeTriple {
        curve_3d: triple_a.curve_3d.clone(),
        pcurve_a: triple_a.pcurve_a.clone(),
        pcurve_b: triple_b.pcurve_a.clone(),
        t_start: triple_a.t_start,
        t_end: triple_a.t_end,
        topology: triple_a.topology,
        parameter_map: ParameterMap::Identity,
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
