//! Trig-pair analytic dispatch (plane × cylinder, plane × sphere,
//! plane × cone). Produces analytic primitives for `curve_3d` and
//! both pcurves with `ParameterMap::Identity` — matches commercial-kernel
//! and OpenCASCADE practice for the canonical analytic-pair
//! cases. See ADR §13.

use crate::geometry::edge_triple::analytic_trig_cone::try_plane_cone_perpendicular;
use crate::geometry::edge_triple::analytic_trig_cone_cone::try_cone_cone_coaxial;
use crate::geometry::edge_triple::analytic_trig_cone_torus_axial::try_cone_torus_axial;
use crate::geometry::edge_triple::analytic_trig_cyl_cone::try_cyl_cone_coaxial;
use crate::geometry::edge_triple::analytic_trig_cyl_cone_general::try_cyl_cone_general;
use crate::geometry::edge_triple::analytic_trig_cyl_cyl::try_cyl_cyl_parallel;
use crate::geometry::edge_triple::analytic_trig_cyl_cyl_perp::try_cyl_cyl_perpendicular;
use crate::geometry::edge_triple::analytic_trig_cyl_sphere::try_cyl_sphere_axial;
use crate::geometry::edge_triple::analytic_trig_plane_cone_hyperbola::try_plane_cone_oblique_hyperbola;
use crate::geometry::edge_triple::analytic_trig_plane_cone_oblique::try_plane_cone_oblique_ellipse;
use crate::geometry::edge_triple::analytic_trig_plane_cone_parabola::try_plane_cone_oblique_parabola;
use crate::geometry::edge_triple::analytic_trig_plane_cyl::try_plane_cylinder_perpendicular;
use crate::geometry::edge_triple::analytic_trig_plane_cyl_oblique::try_plane_cylinder_oblique;
use crate::geometry::edge_triple::analytic_trig_plane_sphere::try_plane_sphere_polar;
use crate::geometry::edge_triple::analytic_trig_plane_sphere_general::try_plane_sphere_general;
use crate::geometry::edge_triple::analytic_trig_plane_torus::try_plane_torus_perpendicular;
use crate::geometry::edge_triple::analytic_trig_plane_torus_oblique::try_plane_torus_oblique;
use crate::geometry::edge_triple::analytic_trig_plane_torus_villarceau::try_plane_torus_villarceau;
use crate::geometry::edge_triple::analytic_trig_sphere_cone::try_sphere_cone_axial;
use crate::geometry::edge_triple::analytic_trig_sphere_cone_general::try_sphere_cone_general;
use crate::geometry::edge_triple::analytic_trig_sphere_sphere::try_sphere_sphere_polar;
use crate::geometry::edge_triple::analytic_trig_sphere_sphere_general::try_sphere_sphere_general;
use crate::geometry::edge_triple::analytic_trig_sphere_torus_axial::try_sphere_torus_axial;
use crate::geometry::edge_triple::error::EdgeTripleError;
use crate::geometry::edge_triple::inputs::{EdgeAnchors, EdgeHint};
use crate::geometry::edge_triple::types::EdgeTriple;
use crate::geometry::surface::SurfaceGeometry;
use crate::tolerance::Tolerance;

/// Try the analytic trig-pair dispatch.
///
/// Returns `Some(Ok(triple))` on success, `Some(Err(...))` on
/// real failure (anchors off circle, etc.), or `None` if the pair
/// isn't a trig pair we handle here (caller should fall through to
/// the NURBS-via-`intersect_surfaces` path).
pub(super) fn try_dispatch_trig(
    surf_a: &SurfaceGeometry,
    surf_b: &SurfaceGeometry,
    anchors: &EdgeAnchors,
    hint: Option<&EdgeHint>,
    tol: &Tolerance,
) -> Option<Result<EdgeTriple, EdgeTripleError>> {
    use SurfaceGeometry::*;
    match (surf_a, surf_b) {
        (
            Plane {
                origin: po,
                normal: pn,
                u_axis: pu,
                v_axis: pv,
            },
            Cylinder {
                origin: co,
                axis: ca,
                radius: cr,
            },
        ) => {
            try_plane_cylinder_perpendicular(*po, *pn, *pu, *pv, false, *co, *ca, *cr, anchors, tol)
                .or_else(|| {
                    try_plane_cylinder_oblique(
                        *po, *pn, *pu, *pv, false, *co, *ca, *cr, anchors, tol,
                    )
                })
        }
        (
            Cylinder {
                origin: co,
                axis: ca,
                radius: cr,
            },
            Plane {
                origin: po,
                normal: pn,
                u_axis: pu,
                v_axis: pv,
            },
        ) => {
            try_plane_cylinder_perpendicular(*po, *pn, *pu, *pv, true, *co, *ca, *cr, anchors, tol)
                .or_else(|| {
                    try_plane_cylinder_oblique(
                        *po, *pn, *pu, *pv, true, *co, *ca, *cr, anchors, tol,
                    )
                })
        }
        (
            Plane {
                origin: po,
                normal: pn,
                u_axis: pu,
                v_axis: pv,
            },
            Sphere { center, radius },
        ) => try_plane_sphere_polar(*po, *pn, *pu, *pv, false, *center, *radius, anchors, tol)
            .or_else(|| {
                try_plane_sphere_general(*po, *pn, *pu, *pv, false, *center, *radius, anchors, tol)
            }),
        (
            Sphere { center, radius },
            Plane {
                origin: po,
                normal: pn,
                u_axis: pu,
                v_axis: pv,
            },
        ) => try_plane_sphere_polar(*po, *pn, *pu, *pv, true, *center, *radius, anchors, tol)
            .or_else(|| {
                try_plane_sphere_general(*po, *pn, *pu, *pv, true, *center, *radius, anchors, tol)
            }),
        (
            Plane {
                origin: po,
                normal: pn,
                u_axis: pu,
                v_axis: pv,
            },
            Cone {
                apex,
                axis,
                half_angle,
            },
        ) => try_plane_cone_perpendicular(
            *po,
            *pn,
            *pu,
            *pv,
            false,
            *apex,
            *axis,
            *half_angle,
            anchors,
            tol,
        )
        .or_else(|| {
            try_plane_cone_oblique_ellipse(
                *po,
                *pn,
                *pu,
                *pv,
                false,
                *apex,
                *axis,
                *half_angle,
                anchors,
                tol,
            )
        })
        .or_else(|| {
            try_plane_cone_oblique_parabola(
                *po,
                *pn,
                *pu,
                *pv,
                false,
                *apex,
                *axis,
                *half_angle,
                anchors,
                tol,
            )
        })
        .or_else(|| {
            try_plane_cone_oblique_hyperbola(
                *po,
                *pn,
                *pu,
                *pv,
                false,
                *apex,
                *axis,
                *half_angle,
                anchors,
                tol,
            )
        }),
        (
            Cone {
                apex,
                axis,
                half_angle,
            },
            Plane {
                origin: po,
                normal: pn,
                u_axis: pu,
                v_axis: pv,
            },
        ) => try_plane_cone_perpendicular(
            *po,
            *pn,
            *pu,
            *pv,
            true,
            *apex,
            *axis,
            *half_angle,
            anchors,
            tol,
        )
        .or_else(|| {
            try_plane_cone_oblique_ellipse(
                *po,
                *pn,
                *pu,
                *pv,
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
                *po,
                *pn,
                *pu,
                *pv,
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
                *po,
                *pn,
                *pu,
                *pv,
                true,
                *apex,
                *axis,
                *half_angle,
                anchors,
                tol,
            )
        }),
        (
            Sphere {
                center: ca,
                radius: ra,
            },
            Sphere {
                center: cb,
                radius: rb,
            },
        ) => try_sphere_sphere_polar(*ca, *ra, *cb, *rb, anchors, tol)
            .or_else(|| try_sphere_sphere_general(*ca, *ra, *cb, *rb, anchors, tol)),
        (
            Cylinder {
                origin: co,
                axis: ca,
                radius: cr,
            },
            Sphere { center, radius },
        ) => try_cyl_sphere_axial(*co, *ca, *cr, false, *center, *radius, anchors, hint, tol),
        (
            Sphere { center, radius },
            Cylinder {
                origin: co,
                axis: ca,
                radius: cr,
            },
        ) => try_cyl_sphere_axial(*co, *ca, *cr, true, *center, *radius, anchors, hint, tol),
        (
            Cylinder {
                origin: oa,
                axis: aa,
                radius: ra,
            },
            Cylinder {
                origin: ob,
                axis: ab,
                radius: rb,
            },
        ) => try_cyl_cyl_parallel(*oa, *aa, *ra, *ob, *ab, *rb, anchors, hint, tol)
            .or_else(|| try_cyl_cyl_perpendicular(*oa, *aa, *ra, *ob, *ab, *rb, anchors, tol)),
        (
            Plane {
                origin: po,
                normal: pn,
                u_axis: pu,
                v_axis: pv,
            },
            Torus {
                center: tc,
                axis: ta,
                major_radius: t_big,
                minor_radius: t_small,
            },
        ) => try_plane_torus_perpendicular(
            *po, *pn, *pu, *pv, false, *tc, *ta, *t_big, *t_small, anchors, tol,
        )
        .or_else(|| {
            try_plane_torus_villarceau(
                *po, *pn, *pu, *pv, false, *tc, *ta, *t_big, *t_small, anchors, tol,
            )
        })
        .or_else(|| {
            try_plane_torus_oblique(
                *po, *pn, *pu, *pv, false, *tc, *ta, *t_big, *t_small, anchors, tol,
            )
        }),
        (
            Torus {
                center: tc,
                axis: ta,
                major_radius: t_big,
                minor_radius: t_small,
            },
            Plane {
                origin: po,
                normal: pn,
                u_axis: pu,
                v_axis: pv,
            },
        ) => try_plane_torus_perpendicular(
            *po, *pn, *pu, *pv, true, *tc, *ta, *t_big, *t_small, anchors, tol,
        )
        .or_else(|| {
            try_plane_torus_villarceau(
                *po, *pn, *pu, *pv, true, *tc, *ta, *t_big, *t_small, anchors, tol,
            )
        })
        .or_else(|| {
            try_plane_torus_oblique(
                *po, *pn, *pu, *pv, true, *tc, *ta, *t_big, *t_small, anchors, tol,
            )
        }),
        (
            Cylinder {
                origin: co,
                axis: ca,
                radius: cr,
            },
            Cone {
                apex,
                axis,
                half_angle,
            },
        ) => try_cyl_cone_coaxial(
            *co,
            *ca,
            *cr,
            false,
            *apex,
            *axis,
            *half_angle,
            anchors,
            tol,
        )
        .or_else(|| {
            try_cyl_cone_general(
                *co,
                *ca,
                *cr,
                false,
                *apex,
                *axis,
                *half_angle,
                anchors,
                tol,
            )
        }),
        (
            Cone {
                apex,
                axis,
                half_angle,
            },
            Cylinder {
                origin: co,
                axis: ca,
                radius: cr,
            },
        ) => try_cyl_cone_coaxial(*co, *ca, *cr, true, *apex, *axis, *half_angle, anchors, tol)
            .or_else(|| {
                try_cyl_cone_general(*co, *ca, *cr, true, *apex, *axis, *half_angle, anchors, tol)
            }),
        (
            Sphere { center, radius },
            Cone {
                apex,
                axis,
                half_angle,
            },
        ) => try_sphere_cone_axial(
            *center,
            *radius,
            false,
            *apex,
            *axis,
            *half_angle,
            anchors,
            tol,
        )
        .or_else(|| {
            try_sphere_cone_general(
                *center,
                *radius,
                false,
                *apex,
                *axis,
                *half_angle,
                anchors,
                tol,
            )
        }),
        (
            Cone {
                apex,
                axis,
                half_angle,
            },
            Sphere { center, radius },
        ) => try_sphere_cone_axial(
            *center,
            *radius,
            true,
            *apex,
            *axis,
            *half_angle,
            anchors,
            tol,
        )
        .or_else(|| {
            try_sphere_cone_general(
                *center,
                *radius,
                true,
                *apex,
                *axis,
                *half_angle,
                anchors,
                tol,
            )
        }),
        (
            Cone {
                apex: aa,
                axis: axa,
                half_angle: ha,
            },
            Cone {
                apex: ab,
                axis: axb,
                half_angle: hb,
            },
        ) => try_cone_cone_coaxial(*aa, *axa, *ha, *ab, *axb, *hb, anchors, tol),
        (
            Sphere { center, radius },
            Torus {
                center: tc,
                axis: ta,
                major_radius: tor_big,
                minor_radius: tor_small,
            },
        ) => try_sphere_torus_axial(
            *center, *radius, false, *tc, *ta, *tor_big, *tor_small, anchors, tol,
        ),
        (
            Torus {
                center: tc,
                axis: ta,
                major_radius: tor_big,
                minor_radius: tor_small,
            },
            Sphere { center, radius },
        ) => try_sphere_torus_axial(
            *center, *radius, true, *tc, *ta, *tor_big, *tor_small, anchors, tol,
        ),
        (
            Cone {
                apex,
                axis,
                half_angle,
            },
            Torus {
                center: tc,
                axis: ta,
                major_radius: tor_big,
                minor_radius: tor_small,
            },
        ) => try_cone_torus_axial(
            *apex,
            *axis,
            *half_angle,
            false,
            *tc,
            *ta,
            *tor_big,
            *tor_small,
            anchors,
            tol,
        ),
        (
            Torus {
                center: tc,
                axis: ta,
                major_radius: tor_big,
                minor_radius: tor_small,
            },
            Cone {
                apex,
                axis,
                half_angle,
            },
        ) => try_cone_torus_axial(
            *apex,
            *axis,
            *half_angle,
            true,
            *tc,
            *ta,
            *tor_big,
            *tor_small,
            anchors,
            tol,
        ),
        _ => None,
    }
}
