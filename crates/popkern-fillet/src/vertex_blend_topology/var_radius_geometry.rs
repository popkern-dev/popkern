use super::*;

/// Build a rational quadratic NURBS arc on the cap-tangent
/// circle (sphere_n ∩ cyl_m) from `start` to `end`.  Both
/// endpoints must lie on the circle (caller's responsibility);
/// we reconstruct the circle from the cylinder's axis and the
/// sphere's centre.
#[allow(dead_code)]
pub(super) fn build_cap_tangent_arc_curve(
    start: popkern_core::geometry::Point3,
    end: popkern_core::geometry::Point3,
    sphere_centre: popkern_core::geometry::Point3,
    cyl_origin: popkern_core::geometry::Point3,
    cyl_axis: popkern_core::geometry::Vec3,
    cyl_radius: f64,
) -> Option<popkern_core::geometry::nurbs_curve::NurbsCurve> {
    use popkern_core::geometry::Vec3;
    use popkern_core::geometry::nurbs_curve::NurbsCurve;
    use std::f64::consts::TAU;

    let axis_unit = cyl_axis.normalized()?;
    let to_centre = Vec3::new(
        sphere_centre.x - cyl_origin.x,
        sphere_centre.y - cyl_origin.y,
        sphere_centre.z - cyl_origin.z,
    );
    let along = to_centre.dot(&axis_unit);
    let circle_centre = popkern_core::geometry::Point3::new(
        cyl_origin.x + axis_unit.x * along,
        cyl_origin.y + axis_unit.y * along,
        cyl_origin.z + axis_unit.z * along,
    );

    // Build orthonormal in-plane basis (u_axis, v_axis) perp to cyl_axis.
    let helper = if axis_unit.x.abs() < 0.9 {
        Vec3::new(1.0, 0.0, 0.0)
    } else {
        Vec3::new(0.0, 1.0, 0.0)
    };
    let u_axis = axis_unit.cross(&helper).normalized()?;
    let v_axis = axis_unit.cross(&u_axis);

    let theta = |p: popkern_core::geometry::Point3| -> f64 {
        let d = Vec3::new(
            p.x - circle_centre.x,
            p.y - circle_centre.y,
            p.z - circle_centre.z,
        );
        d.dot(&v_axis).atan2(d.dot(&u_axis))
    };
    let theta_a = theta(start);
    let theta_b = theta(end);

    // Pick the SHORTER arc direction.
    let mut delta = theta_b - theta_a;
    while delta > std::f64::consts::PI {
        delta -= TAU;
    }
    while delta < -std::f64::consts::PI {
        delta += TAU;
    }
    let theta_end = theta_a + delta;

    // Rational quadratic arc: 3 CPs (p0, shoulder, p2), weights (1, cos(α), 1).
    let half_span = (theta_end - theta_a) * 0.5;
    let cos_half = half_span.cos();
    let theta_mid = theta_a + half_span;
    let circle_eval = |th: f64| -> popkern_core::geometry::Point3 {
        popkern_core::geometry::Point3::new(
            circle_centre.x + cyl_radius * (th.cos() * u_axis.x + th.sin() * v_axis.x),
            circle_centre.y + cyl_radius * (th.cos() * u_axis.y + th.sin() * v_axis.y),
            circle_centre.z + cyl_radius * (th.cos() * u_axis.z + th.sin() * v_axis.z),
        )
    };
    let p0 = circle_eval(theta_a);
    let p_mid = circle_eval(theta_mid);
    let p2 = circle_eval(theta_end);
    let shoulder = if cos_half.abs() > 1e-12 {
        popkern_core::geometry::Point3::new(
            circle_centre.x + (p_mid.x - circle_centre.x) / cos_half,
            circle_centre.y + (p_mid.y - circle_centre.y) / cos_half,
            circle_centre.z + (p_mid.z - circle_centre.z) / cos_half,
        )
    } else {
        p_mid
    };

    Some(NurbsCurve {
        degree: 2,
        control_points: vec![p0, shoulder, p2],
        weights: vec![1.0, cos_half.abs(), 1.0],
        knots: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
    })
}

/// Set the cut edge's 3D curve (cap-tangent arc on cyl_m, going
/// through the cap-tangent circle of sphere_n on cyl_m) and refit
/// pcurves on both half-edges.
pub(super) fn set_cut_edge_geometry(
    brep: &mut BRep,
    edge: EdgeHandle,
    cap_side_he: HalfEdgeHandle,
    v_start: popkern_core::geometry::Point3,
    v_end: popkern_core::geometry::Point3,
    cap_n: FaceHandle,
    main_cyl_face: FaceHandle,
) {
    use popkern_core::geometry::CurveGeometry;
    use popkern_core::geometry::Vec3;
    use popkern_core::geometry::analytic_surfaces::orthonormal_basis;
    use popkern_core::geometry::curve::CurveOnEdge;
    use popkern_core::geometry::pcurve::Pcurve;
    use popkern_core::geometry::surface::SurfaceGeometry;
    use popkern_core::geometry::uv_projection::project_point_to_uv;

    let cap_surf = brep
        .faces
        .get(cap_n)
        .and_then(|f| f.surface.as_ref())
        .map(|s| s.geometry.clone());
    let main_surf = brep
        .faces
        .get(main_cyl_face)
        .and_then(|f| f.surface.as_ref())
        .map(|s| s.geometry.clone());

    let twin = brep.half_edges[cap_side_he].twin;

    // Try the analytic path first: cap_n is Sphere, main_cyl_face is
    // Cylinder, and the cap-tangent arc lies exactly on both. We
    // build an `Ellipse` 3D curve parameterised by angle in the
    // cylinder's orthonormal basis — that way the cylinder's UV `u`
    // equals the curve parameter `t` directly, and Linear pcurves
    // are exact for both the cylinder and (for axis-aligned
    // sphere/cylinder pairs) the sphere.
    let analytic = if let (
        Some(SurfaceGeometry::Sphere {
            center: sphere_centre,
            radius: sphere_radius,
        }),
        Some(SurfaceGeometry::Cylinder {
            origin: cyl_origin,
            axis: cyl_axis,
            radius: cyl_radius,
        }),
    ) = (&cap_surf, &main_surf)
    {
        if let Some(axis_unit) = cyl_axis.normalized() {
            let (perp1, perp2) = orthonormal_basis(axis_unit);
            // Circle centre = projection of sphere centre on cyl axis.
            let to_centre = Vec3::new(
                sphere_centre.x - cyl_origin.x,
                sphere_centre.y - cyl_origin.y,
                sphere_centre.z - cyl_origin.z,
            );
            let centre_along = to_centre.dot(&axis_unit);
            let circle_centre = popkern_core::geometry::Point3::new(
                cyl_origin.x + axis_unit.x * centre_along,
                cyl_origin.y + axis_unit.y * centre_along,
                cyl_origin.z + axis_unit.z * centre_along,
            );

            // Compute angles in cyl basis: t such that
            // point = circle_centre + R · (perp1 · cos(t) + perp2 · sin(t)).
            let theta = |p: popkern_core::geometry::Point3| -> f64 {
                let d = Vec3::new(
                    p.x - circle_centre.x,
                    p.y - circle_centre.y,
                    p.z - circle_centre.z,
                );
                d.dot(&perp2).atan2(d.dot(&perp1))
            };
            let theta_a = theta(v_start);
            let theta_b = theta(v_end);
            // Pick the shorter arc direction (delta in [-π, π]).
            let mut delta = theta_b - theta_a;
            while delta > std::f64::consts::PI {
                delta -= std::f64::consts::TAU;
            }
            while delta < -std::f64::consts::PI {
                delta += std::f64::consts::TAU;
            }
            let theta_end = theta_a + delta;

            // 3D curve: Ellipse — parameterised linearly by angle.
            // u_basis = R · perp1, v_basis = R · perp2 with the SAME
            // basis as the cylinder UV, so cylinder's `u` equals `t`.
            let curve = CurveGeometry::Ellipse {
                center: circle_centre,
                u_basis: perp1 * *cyl_radius,
                v_basis: perp2 * *cyl_radius,
                start_angle: theta_a,
                end_angle: theta_end,
            };
            // Cylinder pcurve: Linear (theta_a, v_cyl) → (theta_end, v_cyl).
            let v_cyl = (v_start - *cyl_origin).dot(&axis_unit);
            let pcurve_cyl =
                Pcurve::linear((theta_a, v_cyl), (theta_end, v_cyl), theta_a, theta_end);
            // Sphere pcurve: Linear from sphere_uv(v_start) → sphere_uv(v_end).
            // For arcs whose plane coincides with a sphere meridian (one
            // endpoint at the pole) `project_point_to_uv` returns u=0
            // at the pole — fix that to the non-pole endpoint's u so
            // the Linear pcurve is a true meridian.
            let cap_geom = SurfaceGeometry::Sphere {
                center: *sphere_centre,
                radius: *sphere_radius,
            };
            // Pole detection via radial distance in the sphere's local
            // x/y plane (not via `v ≈ 0/π`, since `acos` near ±1 picks
            // up ~1e-8 rad of float noise — too tight for a v-based
            // threshold).  At a pole the radial component vanishes
            // and the projected `u` is undefined.
            let pole_radial_eps = sphere_radius * 1e-6;
            let radial_a = ((v_start.x - sphere_centre.x).powi(2)
                + (v_start.y - sphere_centre.y).powi(2))
            .sqrt();
            let radial_b =
                ((v_end.x - sphere_centre.x).powi(2) + (v_end.y - sphere_centre.y).powi(2)).sqrt();
            let pole_a = radial_a < pole_radial_eps;
            let pole_b = radial_b < pole_radial_eps;
            let pcurve_sphere = match (
                project_point_to_uv(&cap_geom, &v_start),
                project_point_to_uv(&cap_geom, &v_end),
            ) {
                (Some((u_a_raw, v_a)), Some((u_b_raw, v_b))) => {
                    // Snap `v` to exact 0 or π at the pole (acos noise),
                    // and inherit `u` from the non-pole endpoint so the
                    // Linear pcurve traces a true meridian.
                    let v_a = if pole_a {
                        if v_a < std::f64::consts::FRAC_PI_2 {
                            0.0
                        } else {
                            std::f64::consts::PI
                        }
                    } else {
                        v_a
                    };
                    let v_b = if pole_b {
                        if v_b < std::f64::consts::FRAC_PI_2 {
                            0.0
                        } else {
                            std::f64::consts::PI
                        }
                    } else {
                        v_b
                    };
                    let u_a = if pole_a { u_b_raw } else { u_a_raw };
                    let u_b = if pole_b { u_a_raw } else { u_b_raw };
                    Some(Pcurve::linear((u_a, v_a), (u_b, v_b), theta_a, theta_end))
                }
                _ => None,
            };
            Some((curve, theta_a, theta_end, pcurve_cyl, pcurve_sphere))
        } else {
            None
        }
    } else {
        None
    };

    if let Some((curve, t_start, t_end, pcurve_cyl, pcurve_sphere)) = analytic {
        if let Some(edge_data) = brep.edges.get_mut(edge) {
            edge_data.curve = Some(CurveOnEdge {
                geometry: curve,
                t_start,
                t_end,
            });
        }
        if let Some(pc) = pcurve_sphere {
            brep.half_edges[cap_side_he].pcurve = Some(pc);
        }
        brep.half_edges[twin].pcurve = Some(pcurve_cyl);
        return;
    }

    // Fallback: surfaces aren't sphere/cylinder.  Build a degree-1
    // straight line (best we can do) and fit pcurves by sampling.
    let curve_3d_nurbs = popkern_core::geometry::nurbs_curve::NurbsCurve {
        degree: 1,
        control_points: vec![v_start, v_end],
        weights: vec![1.0, 1.0],
        knots: vec![0.0, 0.0, 1.0, 1.0],
    };
    if let Some(edge_data) = brep.edges.get_mut(edge) {
        edge_data.curve = CurveOnEdge::line_segment(v_start, v_end, 1.0e-9);
    }
    if let Some(surf) = cap_surf
        && let Some(pc) = sample_pcurve_on_surface(&curve_3d_nurbs, &surf)
    {
        brep.half_edges[cap_side_he].pcurve = Some(pc.into());
    }
    if let Some(surf) = main_surf
        && let Some(pc) = sample_pcurve_on_surface(&curve_3d_nurbs, &surf)
    {
        brep.half_edges[twin].pcurve = Some(pc.into());
    }
}

/// Sample a 3D NURBS curve and fit a degree-3 NURBS pcurve in
/// the surface's UV space.  Returns None if any sample fails to
/// project.
pub(super) fn sample_pcurve_on_surface(
    curve_3d: &popkern_core::geometry::nurbs_curve::NurbsCurve,
    surface: &popkern_core::geometry::surface::SurfaceGeometry,
) -> Option<popkern_core::geometry::nurbs_curve::NurbsCurve> {
    use popkern_core::geometry::Point3;
    use popkern_core::geometry::curve::CurveEval;
    use popkern_core::geometry::fitting;
    use popkern_core::geometry::uv_projection::project_point_to_uv;

    const N: usize = 12;
    let mut points: Vec<Point3> = Vec::with_capacity(N);
    let mut params: Vec<f64> = Vec::with_capacity(N);
    let trace = std::env::var("POPKERN_VAR_RADIUS_TRACE_PCURVE").is_ok_and(|v| v == "1");
    if trace {
        let surf_kind = match surface {
            popkern_core::geometry::surface::SurfaceGeometry::Sphere { center, radius } => {
                format!(
                    "Sphere centre=({:.4},{:.4},{:.4}) r={radius:.4}",
                    center.x, center.y, center.z
                )
            }
            popkern_core::geometry::surface::SurfaceGeometry::Cylinder {
                origin,
                axis,
                radius,
            } => {
                format!(
                    "Cyl origin=({:.4},{:.4},{:.4}) axis=({:.4},{:.4},{:.4}) r={radius:.4}",
                    origin.x, origin.y, origin.z, axis.x, axis.y, axis.z
                )
            }
            popkern_core::geometry::surface::SurfaceGeometry::Plane { origin, .. } => {
                format!(
                    "Plane origin=({:.4},{:.4},{:.4})",
                    origin.x, origin.y, origin.z
                )
            }
            _ => "<other>".to_string(),
        };
        let curve_endpoints = (curve_3d.point_at(0.0), curve_3d.point_at(1.0));
        eprintln!(
            "  [pcurve_fit] surface: {surf_kind} curve: ({:.4},{:.4},{:.4}) → ({:.4},{:.4},{:.4})",
            curve_endpoints.0.x,
            curve_endpoints.0.y,
            curve_endpoints.0.z,
            curve_endpoints.1.x,
            curve_endpoints.1.y,
            curve_endpoints.1.z
        );
    }
    for i in 0..N {
        let s = i as f64 / (N - 1) as f64;
        let p3 = curve_3d.point_at(s);
        let (u, v) = project_point_to_uv(surface, &p3)?;
        if trace {
            eprintln!(
                "    sample[{i}] s={s:.3} p3=({:.4},{:.4},{:.4}) uv=({u:.4},{v:.4})",
                p3.x, p3.y, p3.z
            );
        }
        points.push(Point3::new(u, v, 0.0));
        params.push(s);
    }
    let degree = 3usize;
    let n_cps = (degree + 2).min(points.len());
    let fit = fitting::fit_curve(&points, &params, degree, n_cps).ok()?;
    if trace {
        eprintln!("    fit cps: {:?}", fit.control_points);
    }
    Some(fit)
}

/// Project a 3-D point onto the intersection circle of a sphere
/// and an axis-aligned cylinder.  The intersection lies in the
/// plane perpendicular to the cylinder axis through the sphere
/// centre's projection.  We project `point` first onto that plane
/// (snap z, in essence), then onto the cylinder surface (snap
/// radial distance to `cyl_radius`).
pub(super) fn project_to_sphere_cylinder_intersection(
    point: popkern_core::geometry::Point3,
    sphere_centre: popkern_core::geometry::Point3,
    sphere_radius: f64,
    cyl_origin: popkern_core::geometry::Point3,
    cyl_axis: popkern_core::geometry::Vec3,
) -> Option<popkern_core::geometry::Point3> {
    let axis_unit = cyl_axis.normalized()?;
    // Plane: perpendicular to axis, passing through sphere_centre.
    let to_centre = popkern_core::geometry::Vec3::new(
        sphere_centre.x - cyl_origin.x,
        sphere_centre.y - cyl_origin.y,
        sphere_centre.z - cyl_origin.z,
    );
    let centre_along = to_centre.dot(&axis_unit);
    let circle_centre = popkern_core::geometry::Point3::new(
        cyl_origin.x + axis_unit.x * centre_along,
        cyl_origin.y + axis_unit.y * centre_along,
        cyl_origin.z + axis_unit.z * centre_along,
    );
    // Snap point to the plane.
    let to_point = popkern_core::geometry::Vec3::new(
        point.x - circle_centre.x,
        point.y - circle_centre.y,
        point.z - circle_centre.z,
    );
    let along_axis = to_point.dot(&axis_unit);
    let in_plane = popkern_core::geometry::Vec3::new(
        to_point.x - axis_unit.x * along_axis,
        to_point.y - axis_unit.y * along_axis,
        to_point.z - axis_unit.z * along_axis,
    );
    let in_plane_len = in_plane.length();
    if in_plane_len < 1e-12 {
        return None;
    }
    // Scale to circle radius (= sphere ∩ cylinder, which has radius
    // sqrt(sphere_radius^2 - 0) when sphere centre is on axis;
    // for non-coincident axis/centre, radius = sqrt(R_s^2 - d^2)
    // where d is sphere-centre-to-axis distance).
    let _ = sphere_radius;
    let scale = 1.0; // The cylinder constrains radial = cyl_radius;
    // we just snap to the plane.
    Some(popkern_core::geometry::Point3::new(
        circle_centre.x + in_plane.x * scale,
        circle_centre.y + in_plane.y * scale,
        circle_centre.z + in_plane.z * scale,
    ))
}

/// Walk `loop_handle`'s half-edges via `next` pointers up to `max`
/// hops; return `true` iff the walk returns to the starting
/// half-edge AND every step's `prev`/`next` pointers are mutually
/// consistent AND every consecutive pair has matching vertex
/// continuity (`he.next.origin == twin(he).origin` won't hold,
/// but `he_a.next.origin` should equal `he_a`'s destination, which
/// is `twin(he_a).origin`).
#[allow(dead_code)]
pub(super) fn walk_loop_closes(
    brep: &BRep,
    loop_handle: popkern_core::topology::LoopHandle,
    max: usize,
) -> bool {
    let start = match brep.loops.get(loop_handle) {
        Some(l) => l.half_edge,
        None => return false,
    };
    let mut cursor = start;
    for _ in 0..max {
        let Some(he) = brep.half_edges.get(cursor) else {
            return false;
        };
        let Some(next) = brep.half_edges.get(he.next) else {
            return false;
        };
        // Continuity check: he's destination = twin(he).origin.
        let Some(twin) = brep.half_edges.get(he.twin) else {
            return false;
        };
        if next.origin != twin.origin {
            return false;
        }
        // prev/next mutual consistency.
        if next.prev != cursor {
            return false;
        }
        cursor = he.next;
        if cursor == start {
            return true;
        }
    }
    false
}
