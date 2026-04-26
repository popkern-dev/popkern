//! Rigid-transform helpers for BRep geometry.
//!
//! Extracted from `brep.rs` to stay under the 500-line structural
//! cap. Kept `pub(super)` because the callers live in the same
//! topology module.

use crate::geometry::Transform3;
use crate::geometry::curve::CurveGeometry;
use crate::geometry::gregory::{
    GregoryPatch, RectangularGregory, RectangularGregoryBiNonic, RectangularGregoryBiquinticFull,
    RectangularGregoryBiseptic, RectangularGregoryDeg5, RectangularGregoryTwoStripe,
    TriangularGregory,
};
use crate::geometry::surface::SurfaceGeometry;

/// Transform a `SurfaceGeometry` by a rigid transform.
pub(super) fn transform_surface(surface: &SurfaceGeometry, t: &Transform3) -> SurfaceGeometry {
    match surface {
        SurfaceGeometry::Plane {
            origin,
            normal,
            u_axis,
            v_axis,
        } => SurfaceGeometry::Plane {
            origin: t.apply_point(origin),
            normal: t.apply_direction(normal),
            u_axis: t.apply_direction(u_axis),
            v_axis: t.apply_direction(v_axis),
        },
        SurfaceGeometry::Cylinder {
            origin,
            axis,
            radius,
        } => SurfaceGeometry::Cylinder {
            origin: t.apply_point(origin),
            axis: t.apply_direction(axis),
            radius: *radius,
        },
        SurfaceGeometry::Cone {
            apex,
            axis,
            half_angle,
        } => SurfaceGeometry::Cone {
            apex: t.apply_point(apex),
            axis: t.apply_direction(axis),
            half_angle: *half_angle,
        },
        SurfaceGeometry::Sphere { center, radius } => SurfaceGeometry::Sphere {
            center: t.apply_point(center),
            radius: *radius,
        },
        SurfaceGeometry::Torus {
            center,
            axis,
            major_radius,
            minor_radius,
        } => SurfaceGeometry::Torus {
            center: t.apply_point(center),
            axis: t.apply_direction(axis),
            major_radius: *major_radius,
            minor_radius: *minor_radius,
        },
        SurfaceGeometry::Nurbs(nurbs) => {
            let control_points = nurbs
                .control_points
                .iter()
                .map(|p| t.apply_point(p))
                .collect();
            SurfaceGeometry::Nurbs(crate::geometry::NurbsSurface {
                degree_u: nurbs.degree_u,
                degree_v: nurbs.degree_v,
                knots_u: nurbs.knots_u.clone(),
                knots_v: nurbs.knots_v.clone(),
                control_points,
                num_u: nurbs.num_u,
                num_v: nurbs.num_v,
                weights: nurbs.weights.clone(),
            })
        }
        SurfaceGeometry::Gregory(g) => SurfaceGeometry::Gregory(transform_gregory(g, t)),
    }
}

/// Apply a rigid transform to a Gregory patch by transforming every
/// control point (boundary + inner candidates).
fn transform_gregory(g: &GregoryPatch, t: &Transform3) -> GregoryPatch {
    match g {
        GregoryPatch::Rectangular(r) => {
            let mut boundary = r.boundary_cps;
            for row in boundary.iter_mut() {
                for cell in row.iter_mut() {
                    *cell = t.apply_point(cell);
                }
            }
            let inner_cps = r
                .inner_cps
                .map(|(p_u, p_v)| (t.apply_point(&p_u), t.apply_point(&p_v)));
            GregoryPatch::Rectangular(Box::new(RectangularGregory {
                boundary_cps: boundary,
                inner_cps,
            }))
        }
        GregoryPatch::RectangularDeg5(r) => {
            let mut boundary = r.boundary_cps;
            for row in boundary.iter_mut() {
                for cell in row.iter_mut() {
                    *cell = t.apply_point(cell);
                }
            }
            let inner_cps = r
                .inner_cps
                .map(|(p_u, p_v)| (t.apply_point(&p_u), t.apply_point(&p_v)));
            GregoryPatch::RectangularDeg5(Box::new(RectangularGregoryDeg5 {
                boundary_cps: boundary,
                inner_cps,
            }))
        }
        GregoryPatch::RectangularBiquinticFull(r) => {
            let mut grid = r.grid;
            for row in grid.iter_mut() {
                for cell in row.iter_mut() {
                    *cell = t.apply_point(cell);
                }
            }
            let corner_inner_cps = r
                .corner_inner_cps
                .map(|(p_u, p_v)| (t.apply_point(&p_u), t.apply_point(&p_v)));
            GregoryPatch::RectangularBiquinticFull(Box::new(RectangularGregoryBiquinticFull {
                grid,
                corner_inner_cps,
            }))
        }
        GregoryPatch::RectangularBiseptic(r) => {
            let mut grid = r.grid;
            for row in grid.iter_mut() {
                for cell in row.iter_mut() {
                    *cell = t.apply_point(cell);
                }
            }
            let corner_inner_cps = r
                .corner_inner_cps
                .map(|(p_u, p_v)| (t.apply_point(&p_u), t.apply_point(&p_v)));
            GregoryPatch::RectangularBiseptic(Box::new(RectangularGregoryBiseptic {
                grid,
                corner_inner_cps,
            }))
        }
        GregoryPatch::RectangularTwoStripe(r) => {
            let mut grid = r.grid;
            for row in grid.iter_mut() {
                for cell in row.iter_mut() {
                    *cell = t.apply_point(cell);
                }
            }
            GregoryPatch::RectangularTwoStripe(Box::new(RectangularGregoryTwoStripe { grid }))
        }
        GregoryPatch::RectangularBiNonic(r) => {
            let mut grid = r.grid;
            for row in grid.iter_mut() {
                for cell in row.iter_mut() {
                    *cell = t.apply_point(cell);
                }
            }
            GregoryPatch::RectangularBiNonic(Box::new(RectangularGregoryBiNonic { grid }))
        }
        GregoryPatch::Triangular(tr) => {
            let boundary_cps = tr.boundary_cps.map(|p| t.apply_point(&p));
            let inner_cps = tr.inner_cps.map(|p| t.apply_point(&p));
            GregoryPatch::Triangular(Box::new(TriangularGregory {
                boundary_cps,
                inner_cps,
            }))
        }
    }
}

/// Transform a `CurveGeometry` by a rigid transform.
pub(super) fn transform_curve(curve: &CurveGeometry, t: &Transform3) -> CurveGeometry {
    match curve {
        CurveGeometry::Line { origin, direction } => CurveGeometry::Line {
            origin: t.apply_point(origin),
            direction: t.apply_direction(direction),
        },
        CurveGeometry::Circle {
            center,
            normal,
            radius,
        } => CurveGeometry::Circle {
            center: t.apply_point(center),
            normal: t.apply_direction(normal),
            radius: *radius,
        },
        CurveGeometry::Arc {
            center,
            normal,
            radius,
            start_angle,
            end_angle,
        } => CurveGeometry::Arc {
            center: t.apply_point(center),
            normal: t.apply_direction(normal),
            radius: *radius,
            start_angle: *start_angle,
            end_angle: *end_angle,
        },
        CurveGeometry::Ellipse {
            center,
            u_basis,
            v_basis,
            start_angle,
            end_angle,
        } => CurveGeometry::Ellipse {
            center: t.apply_point(center),
            u_basis: t.apply_direction(u_basis),
            v_basis: t.apply_direction(v_basis),
            start_angle: *start_angle,
            end_angle: *end_angle,
        },
        CurveGeometry::Nurbs(nurbs) => {
            let control_points = nurbs
                .control_points
                .iter()
                .map(|p| t.apply_point(p))
                .collect();
            CurveGeometry::Nurbs(crate::geometry::NurbsCurve {
                degree: nurbs.degree,
                knots: nurbs.knots.clone(),
                control_points,
                weights: nurbs.weights.clone(),
            })
        }
    }
}
