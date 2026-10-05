//! Exact carrier geometry predicates and closing-curve authoring for the
//! side-support support-cell solve (split out of `side_support_cell_solve.rs`
//! for the 500-line budget).

use popkern_core::geometry::curve::CurveOnEdge;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{CurveGeometry, Point3};
use popkern_math::PortableMath;

pub(super) fn same_carrier(a: &SurfaceGeometry, b: &SurfaceGeometry) -> bool {
    match (a, b) {
        (
            SurfaceGeometry::Plane {
                origin: oa,
                normal: na,
                ..
            },
            SurfaceGeometry::Plane {
                origin: ob,
                normal: nb,
                ..
            },
        ) => {
            na.cross(nb).length() < 1.0e-9
                && (*ob - *oa).dot(na).abs() < 1.0e-9 * (1.0 + na.length())
        }
        (
            SurfaceGeometry::Cylinder {
                origin: oa,
                axis: aa,
                radius: ra,
            },
            SurfaceGeometry::Cylinder {
                origin: ob,
                axis: ab,
                radius: rb,
            },
        ) => {
            let (Some(ua), Some(ub)) = (aa.normalized(), ab.normalized()) else {
                return false;
            };
            if (ra - rb).abs() > 1.0e-9 || ua.cross(&ub).length() > 1.0e-9 {
                return false;
            }
            let d = *ob - *oa;
            let perp = d - ua * d.dot(&ua);
            perp.length() < 1.0e-9
        }
        _ => false,
    }
}

pub(super) fn point_on(surface: &SurfaceGeometry, p: &Point3, gate: f64) -> bool {
    match surface {
        SurfaceGeometry::Plane { origin, normal, .. } => {
            let Some(n) = normal.normalized() else {
                return false;
            };
            (*p - *origin).dot(&n).abs() <= gate
        }
        SurfaceGeometry::Cylinder {
            origin,
            axis,
            radius,
        } => {
            let Some(a) = axis.normalized() else {
                return false;
            };
            let d = *p - *origin;
            let perp = d - a * d.dot(&a);
            (perp.length() - radius).abs() <= gate
        }
        _ => false,
    }
}

pub(super) fn closing_curve(
    surface: &SurfaceGeometry,
    start: Point3,
    end: Point3,
    gate: f64,
    linear_tol: f64,
) -> Result<CurveOnEdge, String> {
    match surface {
        SurfaceGeometry::Plane { .. } => CurveOnEdge::line_segment(start, end, linear_tol)
            .ok_or_else(|| "degenerate planar closing segment".to_string()),
        SurfaceGeometry::Cylinder {
            origin,
            axis,
            radius,
        } => {
            let a = axis
                .normalized()
                .ok_or_else(|| "cylinder carrier has a degenerate axis".to_string())?;
            let s_axial = (start - *origin).dot(&a);
            let e_axial = (end - *origin).dot(&a);
            if (s_axial - e_axial).abs() > gate {
                // Not an end-circle arc; a straight chord is only exact if
                // both points share a generator.
                let chord_mid = Point3::new(
                    0.5 * (start.x + end.x),
                    0.5 * (start.y + end.y),
                    0.5 * (start.z + end.z),
                );
                if point_on(surface, &chord_mid, gate) {
                    return CurveOnEdge::line_segment(start, end, linear_tol)
                        .ok_or_else(|| "degenerate generator closing segment".to_string());
                }
                return Err(format!(
                    "cylinder closing needs a geodesic (axial span {:.3e})",
                    (s_axial - e_axial).abs()
                ));
            }
            let center = *origin + a * s_axial;
            let u = start - center;
            if (u.length() - radius).abs() > gate {
                return Err("closing endpoint is off the cylinder end circle".to_string());
            }
            let mut v = a.cross(&u);
            let w = end - center;
            let cos = u.dot(&w) / (u.length() * w.length()).max(1.0e-300);
            let sin = v.dot(&w) / (v.length() * w.length()).max(1.0e-300);
            let mut sweep = sin.patan2(cos);
            if sweep < 0.0 {
                sweep += std::f64::consts::TAU;
            }
            // Both windings lie on the carrier, so the on-carrier
            // postcondition cannot distinguish them — take the SHORT
            // arc (mirror the v basis to measure the angle the other
            // way round) so the closing stays inside the pocket
            // instead of wrapping the far side of the cylinder.
            if sweep > std::f64::consts::PI {
                v = -v;
                sweep = std::f64::consts::TAU - sweep;
            }
            Ok(CurveOnEdge {
                geometry: CurveGeometry::Ellipse {
                    center,
                    u_basis: u,
                    v_basis: v,
                    start_angle: 0.0,
                    end_angle: sweep,
                },
                t_start: 0.0,
                t_end: sweep,
            })
        }
        _ => Err("unsupported carrier kind for closing curve".to_string()),
    }
}

/// Verify a half-edge's pcurve agrees with its edge curve on the given
/// surface at start/mid/end (surface-eval'd pcurve vs curve points).
pub(super) fn pcurve_agrees_with_curve(
    brep: &popkern_core::topology::BRep,
    he: popkern_core::topology::HalfEdgeHandle,
    gate: f64,
) -> Option<bool> {
    use popkern_core::geometry::curve::CurveEval;
    use popkern_core::geometry::surface::SurfaceEval;
    let data = brep.half_edges.get(he)?;
    let pcurve = data.pcurve.as_ref()?;
    let surface = brep
        .loops
        .get(data.loop_)
        .and_then(|l| brep.faces.get(l.face))
        .and_then(|f| f.surface.as_ref())?;
    let curve = brep.edges.get(data.edge).and_then(|e| e.curve.as_ref())?;
    let (plo, phi) = pcurve.parameter_range();
    // Direction-strict: the pcurve must follow the half-edge's own
    // traversal (origin -> dest), which the shell validators assume.
    let origin_pos = brep.vertices.get(data.origin)?.position;
    let forward = origin_pos.distance_to(&curve.geometry.point_at(curve.t_start))
        <= origin_pos.distance_to(&curve.geometry.point_at(curve.t_end));
    let ok = (0..3).all(|i| {
        let s = i as f64 / 2.0;
        let curve_s = if forward { s } else { 1.0 - s };
        let cp = curve
            .geometry
            .point_at(curve.t_start + curve_s * (curve.t_end - curve.t_start));
        let uv = pcurve.point_at(plo + s * (phi - plo));
        surface.geometry.point_at(uv.x, uv.y).distance_to(&cp) <= gate
    });
    Some(ok)
}

/// Author an exact pcurve for `he` by projecting its edge curve onto the
/// owning face's surface (line-on-plane gets an analytic Linear pcurve;
/// everything else a sampled interpolation with nearest-to-previous
/// continuity). Returns false when projection fails.
pub(super) fn author_exact_pcurve(
    brep: &mut popkern_core::topology::BRep,
    he: popkern_core::topology::HalfEdgeHandle,
    gate: f64,
) -> bool {
    use popkern_core::geometry::curve::CurveEval;
    use popkern_core::geometry::pcurve::Pcurve;
    use popkern_core::geometry::surface::SurfaceEval;
    use popkern_core::geometry::uv_projection::{project_point_to_uv, project_point_to_uv_near};
    let Some(data) = brep.half_edges.get(he) else {
        return false;
    };
    let Some(surface) = brep
        .loops
        .get(data.loop_)
        .and_then(|l| brep.faces.get(l.face))
        .and_then(|f| f.surface.as_ref())
        .map(|s| s.geometry.clone())
    else {
        return false;
    };
    let Some(curve) = brep
        .edges
        .get(data.edge)
        .and_then(|e| e.curve.as_ref())
        .cloned()
    else {
        return false;
    };
    // A pcurve follows its OWN half-edge's traversal; when the
    // half-edge runs against the edge curve, author in reverse.
    let origin_pos = brep
        .vertices
        .get(data.origin)
        .map(|v| v.position)
        .unwrap_or_else(|| curve.geometry.point_at(curve.t_start));
    let forward = origin_pos.distance_to(&curve.geometry.point_at(curve.t_start))
        <= origin_pos.distance_to(&curve.geometry.point_at(curve.t_end));
    let line_on_plane = matches!(curve.geometry, CurveGeometry::Line { .. })
        && matches!(surface, SurfaceGeometry::Plane { .. });
    let n = if line_on_plane { 2 } else { 17 };
    let mut uv_points = Vec::with_capacity(n);
    let mut params = Vec::with_capacity(n);
    let mut prev: Option<(f64, f64)> = None;
    for i in 0..n {
        let s = i as f64 / (n - 1) as f64;
        let curve_s = if forward { s } else { 1.0 - s };
        let p = curve
            .geometry
            .point_at(curve.t_start + curve_s * (curve.t_end - curve.t_start));
        let uv = match prev {
            Some((su, sv)) => project_point_to_uv_near(&surface, &p, su, sv)
                .or_else(|| project_point_to_uv(&surface, &p)),
            None => project_point_to_uv(&surface, &p),
        };
        let Some((u, v)) = uv else {
            return false;
        };
        if surface.point_at(u, v).distance_to(&p) > gate {
            return false;
        }
        uv_points.push(popkern_core::geometry::Point3::new(u, v, 0.0));
        params.push(s);
        prev = Some((u, v));
    }
    let pcurve = if line_on_plane {
        Pcurve::linear(
            (uv_points[0].x, uv_points[0].y),
            (uv_points[1].x, uv_points[1].y),
            0.0,
            1.0,
        )
    } else {
        match popkern_core::geometry::fitting::interpolate_curve(&uv_points, &params, 3) {
            Ok(fit) => Pcurve::from_nurbs(fit),
            Err(_) => return false,
        }
    };
    brep.half_edges[he].pcurve = Some(pcurve);
    true
}
