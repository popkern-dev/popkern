//! Curve and pcurve geometry helpers for topology repair.

use super::*;

/// Re-project every half-edge pcurve on `face`'s outer loop into
/// the face's surface UV space. Used after retargeting a face's
/// surface (cylinder → plane, etc.) so the pcurves stay consistent.
pub(crate) fn reproject_loop_pcurves_into_plane(brep: &mut BRep, face: FaceHandle) {
    let Some(face_data) = brep.faces.get(face) else {
        return;
    };
    let outer_loop = face_data.outer_loop;
    let plane_geom = match face_data.surface.as_ref().map(|s| &s.geometry) {
        Some(SurfaceGeometry::Plane {
            origin,
            u_axis,
            v_axis,
            ..
        }) => (*origin, *u_axis, *v_axis),
        _ => return,
    };
    let (plane_origin, plane_u, plane_v) = plane_geom;

    let Some(loop_data) = brep.loops.get(outer_loop) else {
        return;
    };
    let start_he = loop_data.half_edge;
    let mut hes: Vec<HalfEdgeHandle> = Vec::new();
    let mut cursor = start_he;
    let mut safety = 64;
    loop {
        if safety == 0 {
            break;
        }
        safety -= 1;
        hes.push(cursor);
        let Some(he) = brep.half_edges.get(cursor) else {
            break;
        };
        cursor = he.next;
        if cursor == start_he {
            break;
        }
    }

    for he in hes {
        // Project the underlying edge's 3D curve into plane UV.
        let edge = brep.half_edges[he].edge;
        let (curve_3d, window) = match brep.edges.get(edge).and_then(|e| e.curve.as_ref()) {
            Some(c) => (c.geometry.clone(), (c.t_start, c.t_end)),
            None => continue,
        };
        // A control-point projection images the WHOLE 3D curve; the
        // edge may use only a window of it (a rail split at a junction
        // keeps its curve and narrows the window). Planar projection
        // preserves the parameter, so the pcurve takes that window —
        // exactly. Taking the curve's full range instead traced the
        // unsplit rail past both vertices (PK-583: 0.237 on the step's
        // saddle cap planes at a G2 request, where rails are fitted
        // NURBS; G1 rails are lines and take the sampled path).
        let (new_pcurve, projected_window) =
            match project_curve_into_plane_uv(&curve_3d, plane_origin, plane_u, plane_v) {
                Some(p) => (p, Some(window)),
                None => match sample_project_curve_into_plane_uv(
                    &curve_3d,
                    window,
                    plane_origin,
                    plane_u,
                    plane_v,
                ) {
                    Some(p) => (p, None),
                    None => continue,
                },
            };
        // Half-edge orientation: if this half-edge runs the same
        // way as the underlying edge, use the projected pcurve as-is;
        // otherwise reverse (which mirrors the window through the
        // curve's range).
        let edge_he = brep.edges[edge].half_edge;
        let forward = he == edge_he;
        let (t_min, t_max) = new_pcurve.parameter_range();
        let nurbs = if forward {
            new_pcurve
        } else {
            reverse_nurbs(&new_pcurve)
        };
        let mut pcurve: popkern_core::geometry::pcurve::Pcurve = nurbs.into();
        if let Some((w0, w1)) = projected_window
            && w0.is_finite()
            && w1.is_finite()
        {
            (pcurve.t_start, pcurve.t_end) = if forward {
                (w0, w1)
            } else {
                (t_min + t_max - w1, t_min + t_max - w0)
            };
        }
        brep.half_edges[he].pcurve = Some(pcurve);
    }
}

/// Sampling fallback for non-NURBS curves (arcs, lines with finite
/// windows): project 33 samples of the 3D curve into plane UV and
/// return a degree-1 interpolating pcurve (a UV polyline — exact to
/// chord error ~(arc/32)^2, well under render/export tolerance).
/// Straight-chord Linear pcurves for genuinely curved edges were the
/// last structural source of tessellation cracks (PK-095: 0.04-0.1
/// gaps at the base-rim fillet-end arcs).
pub(super) fn sample_project_curve_into_plane_uv(
    curve: &CurveGeometry,
    window: (f64, f64),
    plane_origin: Point3,
    plane_u: popkern_core::geometry::Vec3,
    plane_v: popkern_core::geometry::Vec3,
) -> Option<NurbsCurve> {
    use popkern_core::geometry::curve::CurveEval;
    let (t0, t1) = window;
    if !(t0.is_finite() && t1.is_finite()) || t0 == t1 {
        return None;
    }
    const N: usize = 33;
    let mut cps = Vec::with_capacity(N);
    for i in 0..N {
        let t = t0 + (t1 - t0) * i as f64 / (N - 1) as f64;
        let p = curve.point_at(t);
        let rel = popkern_core::geometry::Vec3::new(
            p.x - plane_origin.x,
            p.y - plane_origin.y,
            p.z - plane_origin.z,
        );
        let u_param = rel.dot(&plane_u) / plane_u.dot(&plane_u);
        let v_param = rel.dot(&plane_v) / plane_v.dot(&plane_v);
        cps.push(Point3::new(u_param, v_param, 0.0));
    }
    // Degree-1 interpolating NURBS: knots [0,0,1/(N-1)..,1,1].
    let mut knots = Vec::with_capacity(N + 2);
    knots.push(0.0);
    for i in 0..N {
        knots.push(i as f64 / (N - 1) as f64);
    }
    knots.push(1.0);
    Some(NurbsCurve {
        degree: 1,
        knots,
        control_points: cps,
        weights: vec![1.0; N],
    })
}

/// Project a 3D curve into a plane's (u, v) parameter space.
/// Plane parameterised as `origin + u·u_axis + v·v_axis`. Returns
/// the pcurve as a NURBS where each control point is `(u, v, 0)`.
pub(super) fn project_curve_into_plane_uv(
    curve: &CurveGeometry,
    plane_origin: Point3,
    plane_u: popkern_core::geometry::Vec3,
    plane_v: popkern_core::geometry::Vec3,
) -> Option<NurbsCurve> {
    let nurbs = match curve {
        CurveGeometry::Nurbs(n) => n.clone(),
        // For non-NURBS curves, sample and fit. Skip for the first
        // pass — analytic projection is what we need.
        _ => return None,
    };
    let projected_cps: Vec<Point3> = nurbs
        .control_points
        .iter()
        .map(|p| {
            let rel = popkern_core::geometry::Vec3::new(
                p.x - plane_origin.x,
                p.y - plane_origin.y,
                p.z - plane_origin.z,
            );
            let u = rel.dot(&plane_u) / plane_u.length().max(1e-12).powi(2) * plane_u.length();
            let v = rel.dot(&plane_v) / plane_v.length().max(1e-12).powi(2) * plane_v.length();
            // The plane is parameterised by point = origin + u·u_axis + v·v_axis,
            // so (u, v) = (rel · û, rel · v̂) only when u_axis, v_axis are unit.
            // For the general case, dot products give projection lengths and we
            // need to divide by axis-length squared to get the parameter.
            let _ = (u, v);
            let u_param = rel.dot(&plane_u) / plane_u.dot(&plane_u);
            let v_param = rel.dot(&plane_v) / plane_v.dot(&plane_v);
            Point3::new(u_param, v_param, 0.0)
        })
        .collect();
    Some(NurbsCurve {
        degree: nurbs.degree,
        knots: nurbs.knots,
        control_points: projected_cps,
        weights: nurbs.weights,
    })
}

pub(super) fn reverse_nurbs(curve: &NurbsCurve) -> NurbsCurve {
    let (t_min, t_max) = curve.parameter_range();
    let mut reversed_knots: Vec<f64> = curve
        .knots
        .iter()
        .rev()
        .map(|&k| t_min + t_max - k)
        .collect();
    reversed_knots.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let mut reversed_cps = curve.control_points.clone();
    reversed_cps.reverse();
    let mut reversed_weights = curve.weights.clone();
    reversed_weights.reverse();
    NurbsCurve {
        degree: curve.degree,
        knots: reversed_knots,
        control_points: reversed_cps,
        weights: reversed_weights,
    }
}

pub(super) fn face_loop_contains_vertex(
    brep: &BRep,
    face: FaceHandle,
    vertex: VertexHandle,
) -> bool {
    let Some(f) = brep.faces.get(face) else {
        return false;
    };
    let Some(loop_data) = brep.loops.get(f.outer_loop) else {
        return false;
    };
    let start_he = loop_data.half_edge;
    let mut cursor = start_he;
    let mut safety = 64;
    loop {
        if safety == 0 {
            return false;
        }
        safety -= 1;
        let Some(he) = brep.half_edges.get(cursor) else {
            return false;
        };
        if he.origin == vertex {
            return true;
        }
        cursor = he.next;
        if cursor == start_he {
            return false;
        }
    }
}

pub(super) fn slice_curve(curve: &NurbsCurve, lo: f64, hi: f64) -> Option<NurbsCurve> {
    let (full_min, full_max) = curve.parameter_range();
    if !(full_min..=full_max).contains(&lo) || !(full_min..=full_max).contains(&hi) {
        return None;
    }
    if hi <= lo + 1e-12 {
        return None;
    }
    let (_, right) = curve.split_at(lo).ok()?;
    let (rem_min, _) = right.parameter_range();
    let cut = rem_min + (hi - lo);
    let (left, _) = right.split_at(cut).ok()?;
    Some(left)
}

#[allow(dead_code)]
fn _force_use(_e: EdgeHandle) {}

pub(super) fn off_cyl_distance(
    p: &Point3,
    origin: Point3,
    axis: popkern_core::geometry::Vec3,
    radius: f64,
) -> f64 {
    let len = axis.length().max(1e-12);
    let axis_unit = popkern_core::geometry::Vec3::new(axis.x / len, axis.y / len, axis.z / len);
    let rel = popkern_core::geometry::Vec3::new(p.x - origin.x, p.y - origin.y, p.z - origin.z);
    let axial = rel.dot(&axis_unit);
    let radial = popkern_core::geometry::Vec3::new(
        rel.x - axial * axis_unit.x,
        rel.y - axial * axis_unit.y,
        rel.z - axial * axis_unit.z,
    );
    radial.length() - radius
}

pub(super) fn anchor_off_cyl_distance(p: &Point3, surf: &SurfaceGeometry) -> f64 {
    if let SurfaceGeometry::Cylinder {
        origin,
        axis,
        radius,
    } = surf
    {
        off_cyl_distance(p, *origin, *axis, *radius).abs()
    } else {
        0.0
    }
}

/// Linearly rescale a NURBS curve's knot vector so its parameter
/// range becomes [0, 1]. The curve geometry is unchanged — only
/// the parameterisation. Required because `split_face` hardcodes
/// `CurveOnEdge { t_start: 0.0, t_end: 1.0, .. }` for the new
/// arc edge, regardless of the input curve's actual range.
pub(super) fn remap_to_unit_range(curve: &NurbsCurve) -> NurbsCurve {
    let (t_min, t_max) = curve.parameter_range();
    if (t_max - t_min).abs() < 1e-15 {
        return curve.clone();
    }
    let new_knots: Vec<f64> = curve
        .knots
        .iter()
        .map(|&k| (k - t_min) / (t_max - t_min))
        .collect();
    NurbsCurve {
        degree: curve.degree,
        knots: new_knots,
        control_points: curve.control_points.clone(),
        weights: curve.weights.clone(),
    }
}

/// Sample `curve_3d` at N parameter values, project each sample
/// onto the cylinder's (u, v) space using `closest_point`, then
/// fit a degree-3 NURBS interpolant. This produces a pcurve that —
/// unlike `cylinder_uv_curve` from `plane_cylinder.rs` — preserves
/// the cylinder mapping at every interior point of the arc, not
/// just at the control points. Required for ShellGap-clean output.
pub(super) fn fit_cylinder_uv_pcurve(
    curve_3d: &NurbsCurve,
    cyl: &SurfaceGeometry,
) -> Option<NurbsCurve> {
    use popkern_core::geometry::fitting::fit_curve;
    use popkern_core::geometry::surface::SurfaceEval;

    let SurfaceGeometry::Cylinder { .. } = cyl else {
        return None;
    };

    let n_samples = 17;
    let (t_min, t_max) = curve_3d.parameter_range();
    let mut uv_points: Vec<Point3> = Vec::with_capacity(n_samples);
    let mut params: Vec<f64> = Vec::with_capacity(n_samples);

    let mut prev_u: Option<f64> = None;
    for i in 0..n_samples {
        let s = i as f64 / (n_samples - 1) as f64;
        let t = t_min + s * (t_max - t_min);
        let p = curve_3d.point_at(t);
        let (mut u, v) = cyl.closest_point(&p)?;
        // Unwrap u so the sample sequence is continuous (no TAU jumps).
        if let Some(prev) = prev_u {
            let tau = std::f64::consts::TAU;
            while u - prev > tau * 0.5 {
                u -= tau;
            }
            while u - prev < -tau * 0.5 {
                u += tau;
            }
        }
        prev_u = Some(u);
        uv_points.push(Point3::new(u, v, 0.0));
        params.push(s);
    }

    let degree = 3;
    let num_cps = degree + 2;
    fit_curve(&uv_points, &params, degree, num_cps).ok()
}
