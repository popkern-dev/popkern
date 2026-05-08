use popkern_core::geometry::nurbs_curve::NurbsCurve;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{CurveEval, Point3, Vec3};
use popkern_core::topology::{BRep, FaceHandle, HalfEdgeHandle, VertexHandle};
use std::f64::consts::TAU;

/// Find a cylinder surface in the BRep (across all faces)
/// whose surface is a `Cylinder` and whose surface contains
/// both `v_a` and `v_b`. Returns the surface (cloned) on
/// success.
pub(super) fn find_shared_cylinder(
    brep: &BRep,
    v_a: VertexHandle,
    v_b: VertexHandle,
) -> Option<SurfaceGeometry> {
    let p_a = brep.vertices.get(v_a)?.position;
    let p_b = brep.vertices.get(v_b)?.position;
    for (_, face) in brep.faces.iter() {
        let Some(sof) = &face.surface else { continue };
        let SurfaceGeometry::Cylinder {
            origin,
            axis,
            radius,
        } = &sof.geometry
        else {
            continue;
        };
        if point_on_cylinder(*origin, *axis, *radius, p_a)
            && point_on_cylinder(*origin, *axis, *radius, p_b)
        {
            return Some(sof.geometry.clone());
        }
    }
    None
}

fn point_on_cylinder(origin: Point3, axis: Vec3, radius: f64, p: Point3) -> bool {
    let axis_unit = match axis.normalized() {
        Some(a) => a,
        None => return false,
    };
    let v = Vec3::new(p.x - origin.x, p.y - origin.y, p.z - origin.z);
    let along = v.dot(&axis_unit);
    let projection = Point3::new(
        origin.x + axis_unit.x * along,
        origin.y + axis_unit.y * along,
        origin.z + axis_unit.z * along,
    );
    let radial = p.distance_to(&projection);
    (radial - radius).abs() < (radius * 0.01).max(1e-6)
}

/// Walk `face`'s outer loop and find a pair `(he_a, he_b)`
/// of half-edges where `he_a.origin == v_a` and
/// `he_b.origin == v_b`. Both half-edges live in the same
/// loop (a precondition for `mef`). Returns `None` if either
/// vertex isn't on the loop.
pub(super) fn find_loop_he_pair(
    brep: &BRep,
    face: FaceHandle,
    v_a: VertexHandle,
    v_b: VertexHandle,
) -> Option<(HalfEdgeHandle, HalfEdgeHandle)> {
    let face_data = brep.faces.get(face)?;
    let outer_loop = face_data.outer_loop;
    let loop_data = brep.loops.get(outer_loop)?;
    let start = loop_data.half_edge;
    let mut he_a = None;
    let mut he_b = None;
    let mut cursor = start;
    for _ in 0..32 {
        let he = brep.half_edges.get(cursor)?;
        if he.origin == v_a {
            he_a = Some(cursor);
        }
        if he.origin == v_b {
            he_b = Some(cursor);
        }
        cursor = he.next;
        if cursor == start {
            break;
        }
    }
    Some((he_a?, he_b?))
}

/// The vertex one step ahead in the loop (the destination of
/// `he_a`). Used to identify the footprint sandwiched between
/// the kept-vertex pair.
pub(super) fn he_after_origin_in_loop(brep: &BRep, he: HalfEdgeHandle) -> Option<VertexHandle> {
    let next = brep.half_edges.get(he)?.next;
    Some(brep.half_edges.get(next)?.origin)
}

/// Sample the 3-D arc, project each sample to UV on `surface`,
/// and fit a NURBS pcurve through the samples. When `reversed`
/// is true, samples are taken from arc(t_hi) to arc(t_lo) for
/// the twin half-edge direction.
pub(super) fn fit_pcurve_along_arc(
    surface: &SurfaceGeometry,
    arc: &NurbsCurve,
    reversed: bool,
) -> Option<NurbsCurve> {
    use popkern_core::geometry::fitting::fit_curve;
    use popkern_core::geometry::surface::SurfaceEval;
    use popkern_core::geometry::uv_projection::project_point_to_uv;

    let (t0, t1) = arc.parameter_range();
    // 48 samples + degree 9: half-circles need higher fidelity
    // because projection to sphere/cylinder UV charts is
    // non-polynomial and may cross the wrapping seam.
    let n_samples = 48;
    let mut points: Vec<Point3> = Vec::with_capacity(n_samples);
    let mut params: Vec<f64> = Vec::with_capacity(n_samples);
    let surf_has_wrapping_u = matches!(
        surface,
        SurfaceGeometry::Sphere { .. } | SurfaceGeometry::Cylinder { .. }
    );
    let mut prev_u: Option<f64> = None;
    for i in 0..n_samples {
        let s = i as f64 / (n_samples - 1) as f64;
        let t = if reversed {
            t0 + (1.0 - s) * (t1 - t0)
        } else {
            t0 + s * (t1 - t0)
        };
        let p = arc.point_at(t);
        let (mut u, v) = project_point_to_uv(surface, &p)?;
        if surf_has_wrapping_u && let Some(u_prev) = prev_u {
            let mut delta = u - u_prev;
            while delta > std::f64::consts::PI {
                u -= TAU;
                delta -= TAU;
            }
            while delta < -std::f64::consts::PI {
                u += TAU;
                delta += TAU;
            }
        }
        prev_u = Some(u);
        points.push(Point3::new(u, v, 0.0));
        params.push(s);
    }
    let degree = 9usize;
    let n_cps = (degree + 2).min(points.len());
    let mut fitted = fit_curve(&points, &params, degree, n_cps).ok()?;
    let ((_u_lo, _u_hi), (v_lo, v_hi)) = surface.parameter_range();
    for cp in &mut fitted.control_points {
        if v_lo.is_finite() && v_hi.is_finite() {
            cp.y = cp.y.clamp(v_lo, v_hi);
        }
    }
    Some(fitted)
}

/// Build a degree-2 rational NURBS half-circle on the cap
/// sphere/cylinder intersection, going from `v_a` to `v_b`.
pub(super) fn build_cap_tangent_half_circle(
    brep: &BRep,
    _cap_centre: Point3,
    _cap_radius: f64,
    cyl_surf: &SurfaceGeometry,
    v_a: VertexHandle,
    v_b: VertexHandle,
) -> NurbsCurve {
    let (cyl_origin, cyl_axis, cyl_radius) = match cyl_surf {
        SurfaceGeometry::Cylinder {
            origin,
            axis,
            radius,
        } => (*origin, *axis, *radius),
        _ => unreachable!("caller validates cylinder surface"),
    };
    let Some(cyl_axis_unit) = cyl_axis.normalized() else {
        return line_curve_between(brep, v_a, v_b);
    };

    // The cap-tangent circle is the per-fillet sphere/cylinder
    // intersection, in a plane perpendicular to the cylinder axis.
    let p_a = brep.vertices[v_a].position;
    let along = Vec3::new(
        p_a.x - cyl_origin.x,
        p_a.y - cyl_origin.y,
        p_a.z - cyl_origin.z,
    )
    .dot(&cyl_axis_unit);
    let circle_centre = Point3::new(
        cyl_origin.x + cyl_axis_unit.x * along,
        cyl_origin.y + cyl_axis_unit.y * along,
        cyl_origin.z + cyl_axis_unit.z * along,
    );
    let helper = if cyl_axis_unit.x.abs() < 0.9 {
        Vec3::new(1.0, 0.0, 0.0)
    } else {
        Vec3::new(0.0, 1.0, 0.0)
    };
    let Some(u_axis) = cyl_axis_unit.cross(&helper).normalized() else {
        return line_curve_between(brep, v_a, v_b);
    };
    let v_axis = cyl_axis_unit.cross(&u_axis);

    let p_b = brep.vertices[v_b].position;
    let theta = |p: Point3| -> f64 {
        let d = Vec3::new(
            p.x - circle_centre.x,
            p.y - circle_centre.y,
            p.z - circle_centre.z,
        );
        d.dot(&v_axis).atan2(d.dot(&u_axis))
    };
    let theta_a = theta(p_a);
    let theta_b = theta(p_b);

    let mut delta = theta_b - theta_a;
    while delta > std::f64::consts::PI {
        delta -= TAU;
    }
    while delta < -std::f64::consts::PI {
        delta += TAU;
    }
    arc_as_rational_quadratic(
        circle_centre,
        cyl_radius,
        u_axis,
        v_axis,
        theta_a,
        theta_a + delta,
    )
}

fn line_curve_between(brep: &BRep, v_a: VertexHandle, v_b: VertexHandle) -> NurbsCurve {
    NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![brep.vertices[v_a].position, brep.vertices[v_b].position],
        weights: vec![1.0, 1.0],
    }
}

/// Single rational-quadratic arc spanning theta_a..theta_b.
/// CPs are `[p0, shoulder, p2]` with weights
/// `[1, cos(alpha), 1]`, where `alpha` is the half-angle.
fn arc_as_rational_quadratic(
    centre: Point3,
    radius: f64,
    u_axis: Vec3,
    v_axis: Vec3,
    theta_a: f64,
    theta_b: f64,
) -> NurbsCurve {
    let half_span = (theta_b - theta_a) * 0.5;
    let cos_half = half_span.cos();
    let theta_mid = theta_a + half_span;
    let circle_eval = |theta: f64| -> Point3 {
        Point3::new(
            centre.x + radius * (theta.cos() * u_axis.x + theta.sin() * v_axis.x),
            centre.y + radius * (theta.cos() * u_axis.y + theta.sin() * v_axis.y),
            centre.z + radius * (theta.cos() * u_axis.z + theta.sin() * v_axis.z),
        )
    };
    let p0 = circle_eval(theta_a);
    let p_mid = circle_eval(theta_mid);
    let p2 = circle_eval(theta_b);
    let p1 = if cos_half.abs() > 1e-12 {
        Point3::new(
            centre.x + (p_mid.x - centre.x) / cos_half,
            centre.y + (p_mid.y - centre.y) / cos_half,
            centre.z + (p_mid.z - centre.z) / cos_half,
        )
    } else {
        p_mid
    };

    NurbsCurve {
        degree: 2,
        control_points: vec![p0, p1, p2],
        weights: vec![1.0, cos_half.abs(), 1.0],
        knots: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
    }
}
