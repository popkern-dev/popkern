//! Rigid-transform helpers for BRep geometry.
//!
//! Extracted from `brep.rs` to stay under the 500-line structural
//! cap. Kept `pub(super)` because the callers live in the same
//! topology module.

use crate::geometry::Transform3;
use crate::geometry::curve::CurveGeometry;
use crate::geometry::gregory::{
    GregoryPatch, RectangularGregory, RectangularGregoryBiNonic,
    RectangularGregoryBiNonicDualTwist, RectangularGregoryBiquinticFull,
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
        GregoryPatch::RectangularBiNonicDualTwist(r) => {
            let mut grid = r.grid;
            for row in grid.iter_mut() {
                for cell in row.iter_mut() {
                    *cell = t.apply_point(cell);
                }
            }
            let inner_cps = r
                .inner_cps
                .map(|(p_u, p_v)| (t.apply_point(&p_u), t.apply_point(&p_v)));
            GregoryPatch::RectangularBiNonicDualTwist(Box::new(
                RectangularGregoryBiNonicDualTwist { grid, inner_cps },
            ))
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

/// Transform an edge's curve COHERENTLY with its trim parameters.
///
/// `Circle` and `Arc` parameterize angle against a basis DERIVED
/// from the normal (`circle_basis`), which is not covariant: rotating
/// the normal silently changes where every angle points, so a
/// transformed rim edge's stored t-range (and an arc's stored
/// angles) referenced the OLD basis and the evaluated endpoints
/// drifted around the rim by the basis shift (measured: rotate_y of
/// a plain cylinder Part failed loop closure with gap = r·√2 — the
/// chord of the 90° basis rotation; the corpus cross-drilled boss
/// declined at the rotate). The angular kinds get their angles and
/// trim range remapped into the NEW derived basis; reflections flip
/// the sweep faithfully.
pub(super) fn transform_curve_on_edge(
    curve: &mut crate::geometry::curve::CurveOnEdge,
    t: &Transform3,
) {
    use crate::geometry::analytic_curves::circle_basis;
    let old_normal = match &curve.geometry {
        CurveGeometry::Circle { normal, .. } | CurveGeometry::Arc { normal, .. } => Some(*normal),
        _ => None,
    };
    curve.geometry = transform_curve(&curve.geometry, t);
    let Some(old_normal) = old_normal else {
        return;
    };
    let new_normal = match &curve.geometry {
        CurveGeometry::Circle { normal, .. } | CurveGeometry::Arc { normal, .. } => *normal,
        _ => return,
    };
    let (u0, v0) = circle_basis(old_normal);
    let (u1, v1) = circle_basis(new_normal);
    let tu0 = t.apply_direction(&u0);
    let tv0 = t.apply_direction(&v0);
    let (c, s) = (tu0.dot(&u1), tu0.dot(&v1));
    let delta = s.atan2(c);
    // Proper rotations map v0 to (−sinδ·u1 + cosδ·v1); reflections
    // land on the opposite handedness and the sweep direction flips.
    let flipped = tv0.dot(&v1) * c - tv0.dot(&u1) * s < 0.0;
    let map = |theta: f64| {
        if flipped {
            delta - theta
        } else {
            theta + delta
        }
    };
    curve.t_start = map(curve.t_start);
    curve.t_end = map(curve.t_end);
    if let CurveGeometry::Arc {
        start_angle,
        end_angle,
        ..
    } = &mut curve.geometry
    {
        *start_angle = map(*start_angle);
        *end_angle = map(*end_angle);
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
        CurveGeometry::Parabola {
            vertex,
            linear,
            quadratic,
            t_start,
            t_end,
        } => CurveGeometry::Parabola {
            vertex: t.apply_point(vertex),
            linear: t.apply_direction(linear),
            quadratic: t.apply_direction(quadratic),
            t_start: *t_start,
            t_end: *t_end,
        },
        CurveGeometry::Hyperbola {
            center,
            semi_major,
            semi_minor,
            t_start,
            t_end,
        } => CurveGeometry::Hyperbola {
            center: t.apply_point(center),
            semi_major: t.apply_direction(semi_major),
            semi_minor: t.apply_direction(semi_minor),
            t_start: *t_start,
            t_end: *t_end,
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

impl crate::topology::BRep {
    /// Apply a rigid transformation to the entire BRep, returning a new copy.
    ///
    /// Transforms vertex positions, surface geometry on every face, edge
    /// curve geometry, and NURBS control points on pcurves (for non-planar
    /// faces; pcurves on planar faces live in UV space and are
    /// parameterization-invariant for rigid transforms that keep the plane
    /// canonical form, so they are regenerated via `attach_planar_pcurves`).
    pub fn transform(
        &self,
        t: &crate::geometry::transform::Transform3,
    ) -> Result<crate::topology::BRep, crate::error::CoreError> {
        let mut result = self.clone();
        for (_, vertex) in &mut result.vertices {
            vertex.position = t.apply_point(&vertex.position);
        }
        // Rotational charts (cylinder/cone/sphere/torus) anchor U to
        // orthonormal_basis(axis) — a convention of the axis
        // direction, not of the material. A rigid rotation moves the
        // anchor, so every stored pcurve U on such a face is offset
        // by the constant angle between the ROTATED old basis and
        // the new axis's own basis (measured: a rotated arc-wall
        // gusset's pcurves sat 2π-and-a-fraction off its chart and
        // the UV clip kept a phantom intersection line while
        // rejecting nothing real). Compute the shift per face before
        // overwriting the surface, then re-anchor its pcurves.
        let mut face_u_shift: Vec<(crate::topology::FaceHandle, f64)> = Vec::new();
        for (fh, face) in &result.faces {
            let Some(sof) = face.surface.as_ref() else {
                continue;
            };
            let axis_old = match &sof.geometry {
                SurfaceGeometry::Cylinder { axis, .. }
                | SurfaceGeometry::Cone { axis, .. }
                | SurfaceGeometry::Torus { axis, .. } => *axis,
                _ => continue,
            };
            let axis_new = t.apply_direction(&axis_old);
            let (e1_old, _) = crate::geometry::analytic_surfaces::orthonormal_basis(axis_old);
            let (f1, f2) = crate::geometry::analytic_surfaces::orthonormal_basis(axis_new);
            let e1_img = t.apply_direction(&e1_old);
            let delta = e1_img.dot(&f2).atan2(e1_img.dot(&f1));
            if delta.abs() > 1.0e-14 {
                face_u_shift.push((fh, delta));
            }
        }
        for (fh, delta) in face_u_shift {
            let face = &result.faces[fh];
            let mut loops = vec![face.outer_loop];
            loops.extend(face.inner_loops.iter().copied());
            let mut hes = Vec::new();
            for lh in loops {
                let start = result.loops[lh].half_edge;
                if !result.half_edges.contains_key(start) {
                    continue;
                }
                let mut he = start;
                loop {
                    hes.push(he);
                    he = result.half_edges[he].next;
                    if he == start {
                        break;
                    }
                }
            }
            for he in hes {
                if let Some(pc) = result.half_edges[he].pcurve.as_mut() {
                    pc.shift_u(delta);
                }
            }
        }
        for (_, face) in &mut result.faces {
            if let Some(ref mut sof) = face.surface {
                sof.geometry = transform_surface(&sof.geometry, t);
            }
        }
        for (_, edge) in &mut result.edges {
            if let Some(ref mut curve_on_edge) = edge.curve {
                transform_curve_on_edge(curve_on_edge, t);
            }
        }
        result.attach_planar_pcurves()?;
        Ok(result)
    }
}
