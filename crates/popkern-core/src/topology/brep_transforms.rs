//! Rigid-transform helpers for BRep geometry.
//!
//! Extracted from `brep.rs` to stay under the 500-line structural
//! cap. Kept `pub(super)` because the callers live in the same
//! topology module.

use crate::geometry::Transform3;
use crate::geometry::curve::CurveGeometry;
use crate::geometry::surface::SurfaceGeometry;
use popkern_math::PortableMath;

#[path = "brep_transform_gregory.rs"]
mod brep_transform_gregory;
use brep_transform_gregory::transform_gregory;

/// Normalize a transformed direction field. The caller has already
/// verified `t` is a similarity (non-singular linear part), so a
/// transformed non-degenerate direction cannot vanish; the scaled
/// image is kept as a defensive fallback rather than panicking.
///
/// Rigid transforms keep the image length within ~1e-16 of unit; the
/// division is SKIPPED there so translate/rotate stays bit-identical
/// to the pre-similarity-contract behavior (an ULP of direction
/// noise on every transform is enough to flip marginal welds — the
/// Bishop-frame seed lesson).
fn unit_direction(t: &Transform3, v: &crate::geometry::Vec3) -> crate::geometry::Vec3 {
    let image = t.apply_direction(v);
    let len = image.length();
    if (len - 1.0).abs() <= 1.0e-9 {
        return image;
    }
    image.normalized().unwrap_or(image)
}

/// Transform a `SurfaceGeometry` by a similarity transform with
/// uniform scale factor `scale` (as returned by
/// [`Transform3::similarity_scale`]).
///
/// Direction fields (plane normal and UV axes, cylinder/cone/torus
/// axes) stay UNIT LENGTH — downstream offset math multiplies them
/// by radii and distances, and a scaled normal displaces every such
/// offset by the scale factor (PK-231: the notch cube's 0.6-scaled
/// floor normal put a fillet cylinder 0.048 off its tangent plane
/// instead of 0.08). Metric fields (radii) carry the scale instead.
/// Chart consequences of axis normalization (cylinder/cone V is arc
/// length) are handled by the caller via [`Pcurve::scale_v`].
pub(crate) fn transform_surface(
    surface: &SurfaceGeometry,
    t: &Transform3,
    scale: f64,
) -> SurfaceGeometry {
    match surface {
        SurfaceGeometry::Plane {
            origin,
            normal,
            u_axis,
            v_axis,
        } => SurfaceGeometry::Plane {
            origin: t.apply_point(origin),
            normal: unit_direction(t, normal),
            u_axis: unit_direction(t, u_axis),
            v_axis: unit_direction(t, v_axis),
        },
        SurfaceGeometry::Cylinder {
            origin,
            axis,
            radius,
        } => SurfaceGeometry::Cylinder {
            origin: t.apply_point(origin),
            axis: unit_direction(t, axis),
            radius: *radius * scale,
        },
        SurfaceGeometry::Cone {
            apex,
            axis,
            half_angle,
        } => SurfaceGeometry::Cone {
            apex: t.apply_point(apex),
            axis: unit_direction(t, axis),
            half_angle: *half_angle,
        },
        SurfaceGeometry::Sphere {
            center,
            radius,
            axis,
        } => SurfaceGeometry::Sphere {
            center: t.apply_point(center),
            radius: *radius * scale,
            axis: unit_direction(t, axis),
        },
        SurfaceGeometry::Torus {
            center,
            axis,
            major_radius,
            minor_radius,
        } => SurfaceGeometry::Torus {
            center: t.apply_point(center),
            axis: unit_direction(t, axis),
            major_radius: *major_radius * scale,
            minor_radius: *minor_radius * scale,
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
        SurfaceGeometry::PolyhedralSpline(p) => {
            // Transform every control point and keep the provenance:
            // a moved cap sector is still that cap's sector.
            let mut net = p.net;
            for row in net.iter_mut() {
                for cp in row.iter_mut() {
                    *cp = t.apply_point(cp);
                }
            }
            SurfaceGeometry::PolyhedralSpline(Box::new(
                crate::geometry::polyhedral_spline_patch::PolyhedralSplinePatch::new(
                    net,
                    p.provenance,
                ),
            ))
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
pub(crate) fn transform_curve_on_edge(
    curve: &mut crate::geometry::curve::CurveOnEdge,
    t: &Transform3,
    scale: f64,
) {
    use crate::geometry::analytic_curves::circle_basis;
    let old_normal = match &curve.geometry {
        CurveGeometry::Circle { normal, .. } | CurveGeometry::Arc { normal, .. } => Some(*normal),
        _ => None,
    };
    curve.geometry = transform_curve(&curve.geometry, t, scale);
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
    let delta = s.patan2(c);
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

/// Transform a `CurveGeometry` by a similarity transform with
/// uniform scale factor `scale`.
///
/// Circle/Arc radii carry the scale and their normals stay unit
/// (angle parameters are scale-invariant). Parameterization-covariant
/// fields (line direction, ellipse bases, parabola/hyperbola
/// coefficient vectors, NURBS control points) take the scaled image
/// directly — their parameter ranges remain valid because evaluation
/// multiplies the stored vectors by the parameter.
pub(super) fn transform_curve(curve: &CurveGeometry, t: &Transform3, scale: f64) -> CurveGeometry {
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
            normal: unit_direction(t, normal),
            radius: *radius * scale,
        },
        CurveGeometry::Arc {
            center,
            normal,
            radius,
            start_angle,
            end_angle,
        } => CurveGeometry::Arc {
            center: t.apply_point(center),
            normal: unit_direction(t, normal),
            radius: *radius * scale,
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
    /// Apply a similarity transformation (rigid motion, optional
    /// mirror, uniform scale) to the entire BRep, returning a new copy.
    ///
    /// Transforms vertex positions, surface geometry on every face, edge
    /// curve geometry, and NURBS control points on pcurves (for non-planar
    /// faces; pcurves on planar faces live in UV space and are
    /// parameterization-invariant for rigid transforms that keep the plane
    /// canonical form, so they are regenerated via `attach_planar_pcurves`).
    ///
    /// Under a uniform scale `s ≠ 1`, direction fields stay unit,
    /// radii scale by `s`, cylinder/cone pcurve V (arc length along
    /// the unit axis) scales by `s`, and planar inner-loop pcurves
    /// scale isotropically. Vertex/model tolerances are NOT scaled —
    /// callers scaling far from 1 should re-derive their tolerance
    /// from the new model size.
    ///
    /// # Errors
    ///
    /// Returns [`crate::error::CoreError::PreconditionFailed`] when
    /// `t` is not a similarity (non-uniform scale or shear): analytic
    /// surfaces do not stay analytic under such maps (a cylinder
    /// becomes an elliptic cylinder `SurfaceGeometry` cannot
    /// represent), and the previous behavior silently corrupted
    /// direction lengths and radii instead (PK-231).
    pub fn transform(
        &self,
        t: &crate::geometry::transform::Transform3,
    ) -> Result<crate::topology::BRep, crate::error::CoreError> {
        let scale = t.similarity_scale().ok_or_else(|| {
            crate::error::CoreError::PreconditionFailed(
                "BRep::transform requires a similarity transform (rigid motion, optional \
                 mirror, uniform scale); non-uniform scale and shear do not preserve \
                 analytic surface forms"
                    .to_string(),
            )
        })?;
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
                | SurfaceGeometry::Sphere { axis, .. }
                | SurfaceGeometry::Torus { axis, .. } => *axis,
                _ => continue,
            };
            let axis_new = t.apply_direction(&axis_old);
            let (e1_old, _) = crate::geometry::analytic_surfaces::orthonormal_basis(axis_old);
            let (f1, f2) = crate::geometry::analytic_surfaces::orthonormal_basis(axis_new);
            let e1_img = t.apply_direction(&e1_old);
            let delta = e1_img.dot(&f2).patan2(e1_img.dot(&f1));
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
                sof.geometry = transform_surface(&sof.geometry, t, scale);
            }
        }
        for (_, edge) in &mut result.edges {
            if let Some(ref mut curve_on_edge) = edge.curve {
                transform_curve_on_edge(curve_on_edge, t, scale);
            }
        }
        // Chart metric follow-up for uniform scale: cylinder/cone V is
        // arc length along the (unit) axis and plane UV is length
        // along unit axes, so both stretch by `scale`; sphere/torus
        // charts are all-angular and NURBS/Gregory charts live in
        // knot space, all scale-invariant.
        if (scale - 1.0).abs() > 1.0e-12 {
            let face_handles: Vec<crate::topology::FaceHandle> = result.faces.keys().collect();
            for fh in face_handles {
                let Some(sof) = result.faces[fh].surface.as_ref() else {
                    continue;
                };
                enum ChartScale {
                    V,
                    Uv,
                }
                let chart = match &sof.geometry {
                    SurfaceGeometry::Cylinder { .. } | SurfaceGeometry::Cone { .. } => {
                        ChartScale::V
                    }
                    SurfaceGeometry::Plane { .. } => ChartScale::Uv,
                    _ => continue,
                };
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
                        let exact = match chart {
                            ChartScale::V => pc.scale_v(scale),
                            ChartScale::Uv => pc.scale_uv(scale),
                        };
                        if !exact {
                            return Err(crate::error::CoreError::PreconditionFailed(format!(
                                "BRep::transform: pcurve on face {fh:?} cannot be scaled \
                                 exactly in its stored form (uniform scale {scale}); \
                                 re-author the pcurve before scaling"
                            )));
                        }
                    }
                }
            }
        }
        result.attach_planar_pcurves()?;
        Ok(result)
    }
}
