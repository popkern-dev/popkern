//! Exact rational NURBS of cylinder, cone and sphere patches (PK-1142).
//!
//! `to_nurbs_in_range` sampled these into a 9 × 9 grid used as the
//! CONTROL net (degree 3), a surface that leaves the analytic one (a
//! full-turn cylinder of radius 5 read 4.394 from its axis at a quarter
//! turn): the SSI seeder culls patch pairs by those hulls, and a NURBS
//! twin built that way is not the surface. Each is a surface of
//! revolution in the chart's own frame: the angle `u` in rational
//! quadratic arcs of at most 90° (the torus's pattern, knots at the
//! chart angles), the profile in `v` (a line for the cylinder and the
//! cone, a meridian arc in colatitude for the sphere). At every knot the
//! parameters are the chart's; between them the rational
//! parameterisation departs from the angle, but every point lies on the
//! surface to rounding.

use super::surface_torus_nurbs::{AxisArcPattern, axis_arc_pattern};
use crate::geometry::analytic_surfaces::orthonormal_basis;
use crate::geometry::nurbs_surface::NurbsSurface;
use crate::geometry::{Point3, Vec3};

/// The exact patch of a cylinder, cone or sphere over the window (an
/// angular span within one period; a sphere's colatitude too); `None`
/// for other surfaces and windows.
pub(super) fn exact_quadric(
    s: &crate::geometry::surface::SurfaceGeometry,
    u_range: (f64, f64),
    v_range: (f64, f64),
) -> Option<NurbsSurface> {
    use crate::geometry::surface::SurfaceGeometry;
    use popkern_math::PortableMath;
    let within = |r: (f64, f64)| r.1 - r.0 > 1.0e-9 && r.1 - r.0 < std::f64::consts::TAU + 1.0e-9;
    if !within(u_range) || v_range.1 <= v_range.0 {
        return None;
    }
    match s {
        SurfaceGeometry::Cylinder {
            origin,
            axis,
            radius,
        } => Some(exact_cylinder_nurbs(
            *origin, *axis, *radius, u_range, v_range,
        )),
        SurfaceGeometry::Cone {
            apex,
            axis,
            half_angle,
        } => Some(exact_cone_nurbs(
            *apex,
            *axis,
            half_angle.ptan(),
            u_range,
            v_range,
        )),
        SurfaceGeometry::Sphere {
            center,
            radius,
            axis,
        } if within(v_range) => Some(exact_sphere_nurbs(
            *center, *radius, *axis, u_range, v_range,
        )),
        _ => None,
    }
}

/// The tensor product of the angular arcs `pu` (in the frame `e1, e2`)
/// with a profile given as control points `(radial, axial, weight)`
/// over `knots_v`: `base + (e1 cos + e2 sin) radial + axis axial`.
fn revolve(
    base: Point3,
    (e1, e2, axis): (Vec3, Vec3, Vec3),
    pu: &AxisArcPattern,
    profile: &[(f64, f64, f64)],
    (degree_v, knots_v): (usize, Vec<f64>),
) -> NurbsSurface {
    let mut control_points = Vec::with_capacity(pu.plane.len() * profile.len());
    let mut weights = Vec::with_capacity(pu.plane.len() * profile.len());
    for (cu, wu) in pu.plane.iter().zip(pu.weights.iter()) {
        for &(radial, axial, wv) in profile {
            control_points.push(base + (e1 * cu.0 + e2 * cu.1) * radial + axis * axial);
            weights.push(wu * wv);
        }
    }
    NurbsSurface {
        degree_u: 2,
        degree_v,
        knots_u: pu.knots.clone(),
        knots_v,
        num_u: pu.plane.len(),
        num_v: profile.len(),
        control_points,
        weights,
    }
}

/// A cylinder's patch: `origin + axis v + r (cos u p1 + sin u p2)`.
fn exact_cylinder_nurbs(
    origin: Point3,
    axis: Vec3,
    radius: f64,
    u_range: (f64, f64),
    v_range: (f64, f64),
) -> NurbsSurface {
    let (p1, p2) = orthonormal_basis(axis);
    let profile = [(radius, v_range.0, 1.0), (radius, v_range.1, 1.0)];
    let knots = vec![v_range.0, v_range.0, v_range.1, v_range.1];
    revolve(
        origin,
        (p1, p2, axis),
        &axis_arc_pattern(u_range),
        &profile,
        (1, knots),
    )
}

/// A cone's patch: `apex + axis v + v tan α (cos u p1 + sin u p2)`.
fn exact_cone_nurbs(
    apex: Point3,
    axis: Vec3,
    tan: f64,
    u_range: (f64, f64),
    v_range: (f64, f64),
) -> NurbsSurface {
    let (p1, p2) = orthonormal_basis(axis);
    let profile = [
        (v_range.0 * tan, v_range.0, 1.0),
        (v_range.1 * tan, v_range.1, 1.0),
    ];
    let knots = vec![v_range.0, v_range.0, v_range.1, v_range.1];
    revolve(
        apex,
        (p1, p2, axis),
        &axis_arc_pattern(u_range),
        &profile,
        (1, knots),
    )
}

/// A sphere's patch in its chart frame (the pole `axis` at colatitude
/// `v = 0`): the meridian arc `(R sin v, R cos v)` revolved.
fn exact_sphere_nurbs(
    center: Point3,
    radius: f64,
    axis: Vec3,
    u_range: (f64, f64),
    v_range: (f64, f64),
) -> NurbsSurface {
    // The chart's frame (`sphere_to_world`): the world axes for a +z
    // pole, else the orthonormal basis of the pole.
    let (e1, e2) = if axis == Vec3::Z {
        (Vec3::new(1.0, 0.0, 0.0), Vec3::new(0.0, 1.0, 0.0))
    } else {
        orthonormal_basis(axis)
    };
    let pv = axis_arc_pattern(v_range);
    // (cos v, sin v) control coordinates give (axial, radial).
    let profile: Vec<(f64, f64, f64)> = pv
        .plane
        .iter()
        .zip(pv.weights.iter())
        .map(|(c, w)| (radius * c.1, radius * c.0, *w))
        .collect();
    revolve(
        center,
        (e1, e2, axis),
        &axis_arc_pattern(u_range),
        &profile,
        (2, pv.knots),
    )
}

#[cfg(test)]
#[path = "surface_quadric_nurbs_tests.rs"]
mod tests;
