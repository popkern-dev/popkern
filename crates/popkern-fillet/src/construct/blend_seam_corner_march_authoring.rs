//! PK-152 stage 2 authoring: chart polylines/pcurves for the marched
//! seam, the exact cylinder×plane stub ellipse, and the rail edge
//! shortened along its own line (see `blend_seam_corner_march` for
//! the repair orchestration).

use popkern_core::geometry::analytic_surfaces::orthonormal_basis;
use popkern_core::geometry::fitting::interpolate_curve;
use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::geometry::{CurveEval, CurveGeometry, CurveOnEdge, Point3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle, HalfEdgeHandle, VertexHandle};

use crate::types::BlendSeamCornerDeclineReason as Reason;

use super::super::blend_seam_corner_repair::commit_edge_authoring;
use super::super::terminal_arc_trim::StripCylinder;
use super::AnchorEdge;
use super::solver::SupportPlane;

/// Chart polyline of the marched points in one strip's chart,
/// sequentially period-unwrapped from an optional anchor.
pub(super) fn chart_polyline(
    cyl: &StripCylinder,
    pts: &[Point3],
    anchor_u: Option<f64>,
) -> Vec<(f64, f64)> {
    let (p1, p2) = orthonormal_basis(cyl.axis);
    let mut out = Vec::with_capacity(pts.len());
    let mut prev_u: Option<f64> = anchor_u;
    for p in pts {
        let d = *p - cyl.origin;
        let radial = d - cyl.axis * d.dot(&cyl.axis);
        let mut u = radial.dot(&p2).atan2(radial.dot(&p1));
        if let Some(pu) = prev_u {
            while u - pu > std::f64::consts::PI {
                u -= std::f64::consts::TAU;
            }
            while pu - u > std::f64::consts::PI {
                u += std::f64::consts::TAU;
            }
        }
        prev_u = Some(u);
        out.push((u, d.dot(&cyl.axis)));
    }
    out
}

/// Degree-3 interpolating pcurve through chart samples with the
/// given normalized parameters.
pub(in crate::construct) fn chart_pcurve(
    uv: &[(f64, f64)],
    params: &[f64],
) -> Result<Pcurve, Reason> {
    let pts: Vec<Point3> = uv.iter().map(|(u, v)| Point3::new(*u, *v, 0.0)).collect();
    let curve = interpolate_curve(&pts, params, 3).map_err(|_| Reason::PcurveClosureFailed)?;
    Ok(Pcurve::from_nurbs(curve))
}

/// Chord-length parameters normalized to [0, 1].
pub(in crate::construct) fn chord_params(pts: &[Point3]) -> Vec<f64> {
    let mut acc = vec![0.0_f64];
    for w in pts.windows(2) {
        let d = (w[1] - w[0]).length();
        acc.push(acc.last().copied().unwrap_or(0.0) + d.max(1e-15));
    }
    let total = *acc.last().unwrap_or(&1.0);
    acc.iter().map(|a| a / total).collect()
}

/// Exact cylinder×plane ellipse arc between two points on both
/// surfaces, using the conjugate-diameter form (chart-angle
/// parameterized, short branch).
pub(super) fn cylinder_plane_arc(
    cyl: &StripCylinder,
    plane: &SupportPlane,
    from: Point3,
    to: Point3,
    tolerance: &Tolerance,
) -> Result<(CurveOnEdge, f64, f64), Reason> {
    let (p1, p2) = orthonormal_basis(cyl.axis);
    let denom = plane.normal.dot(&cyl.axis);
    if denom.abs() < 1e-9 {
        return Err(Reason::StubUnresolved);
    }
    let t_c = plane.normal.dot(&(plane.origin - cyl.origin)) / denom;
    let center = cyl.origin + cyl.axis * t_c;
    let u_basis = (p1 - cyl.axis * (plane.normal.dot(&p1) / denom)) * cyl.radius;
    let v_basis = (p2 - cyl.axis * (plane.normal.dot(&p2) / denom)) * cyl.radius;
    let angle_of = |p: Point3| -> f64 {
        let d = p - cyl.origin;
        let radial = d - cyl.axis * d.dot(&cyl.axis);
        radial.dot(&p2).atan2(radial.dot(&p1))
    };
    let th_from = angle_of(from);
    let mut th_to = angle_of(to);
    while th_to - th_from > std::f64::consts::PI {
        th_to -= std::f64::consts::TAU;
    }
    while th_from - th_to > std::f64::consts::PI {
        th_to += std::f64::consts::TAU;
    }
    if (th_to - th_from).abs() < 1e-9 {
        return Err(Reason::StubUnresolved);
    }
    let arc = CurveOnEdge {
        geometry: CurveGeometry::Ellipse {
            center,
            u_basis,
            v_basis,
            start_angle: th_from,
            end_angle: th_to,
        },
        t_start: th_from,
        t_end: th_to,
    };
    let pos_tol = (tolerance.linear() * 10.0).max(1e-9);
    if (arc.geometry.point_at(th_from) - from).length() > pos_tol
        || (arc.geometry.point_at(th_to) - to).length() > pos_tol
    {
        return Err(Reason::StubUnresolved);
    }
    Ok((arc, th_from, th_to))
}

/// Sinusoidal pcurve of a plane in a cylinder chart over
/// `[u_from, u_to]` (see stage 1: canonical plane×cylinder form).
pub(super) fn plane_sinusoid(
    cyl: &StripCylinder,
    plane: &SupportPlane,
    u_from: f64,
    u_to: f64,
) -> Result<Pcurve, Reason> {
    let (p1, p2) = orthonormal_basis(cyl.axis);
    let denom = plane.normal.dot(&cyl.axis);
    if denom.abs() < 1e-9 {
        return Err(Reason::StubUnresolved);
    }
    let n = plane.normal;
    let v_const = -n.dot(&(cyl.origin - plane.origin)) / denom;
    let v_cos = -cyl.radius * n.dot(&p1) / denom;
    let v_sin = -cyl.radius * n.dot(&p2) / denom;
    Ok(Pcurve::sinusoidal(0.0, v_const, v_cos, v_sin, u_from, u_to))
}

pub(super) fn params_reversed(params: &[f64]) -> Vec<f64> {
    let mut rev: Vec<f64> = params.iter().map(|t| 1.0 - t).collect();
    rev.reverse();
    rev
}

/// Author the stub edge's pcurves on the larger strip (sinusoid) and
/// the support plane (projected endpoints; elliptical when the chart
/// image is axis-orthogonal, sampled otherwise).
#[allow(clippy::too_many_arguments)]
pub(super) fn author_stub_pcurves(
    brep: &BRep,
    stub: &AnchorEdge,
    arc: &CurveOnEdge,
    strip_off: FaceHandle,
    cyl_off: &StripCylinder,
    plane: &SupportPlane,
    tolerance: &Tolerance,
) -> Result<Vec<(FaceHandle, HalfEdgeHandle, Pcurve)>, Reason> {
    let e = brep.edges.get(stub.edge).ok_or(Reason::StubUnresolved)?;
    let he = brep
        .half_edges
        .get(e.half_edge)
        .ok_or(Reason::StubUnresolved)?;
    let tw = brep.half_edges.get(he.twin).ok_or(Reason::StubUnresolved)?;
    let curve = arc;
    let mut staged = Vec::with_capacity(2);
    for (heh, hh) in [(e.half_edge, he), (he.twin, tw)] {
        let face = brep
            .loops
            .get(hh.loop_)
            .map(|l| l.face)
            .ok_or(Reason::StubUnresolved)?;
        // Sample the committed arc curve along this half-edge's
        // direction and project into the face chart.
        let (t0, t1) = (curve.t_start, curve.t_end);
        let forward = e.half_edge == heh;
        let n = 17;
        let mut samples3d = Vec::with_capacity(n);
        for k in 0..n {
            let f = k as f64 / (n - 1) as f64;
            let f = if forward { f } else { 1.0 - f };
            samples3d.push(curve.geometry.point_at(t0 + (t1 - t0) * f));
        }
        let pcurve = if face == strip_off {
            let uv = chart_polyline(cyl_off, &samples3d, None);
            plane_sinusoid(
                cyl_off,
                plane,
                uv.first().ok_or(Reason::StubUnresolved)?.0,
                uv.last().ok_or(Reason::StubUnresolved)?.0,
            )?
        } else {
            // Support-plane chart: project the samples exactly.
            let uv: Vec<Point3> = samples3d
                .iter()
                .map(|p| {
                    let d = *p - plane.origin;
                    Point3::new(d.dot(&plane.u_axis), d.dot(&plane.v_axis), 0.0)
                })
                .collect();
            let params: Vec<f64> = (0..n).map(|k| k as f64 / (n - 1) as f64).collect();
            let fit = interpolate_curve(&uv, &params, 3).map_err(|_| Reason::StubUnresolved)?;
            Pcurve::from_nurbs(fit)
        };
        // Sample gate: the pcurve trace must reproduce the arc.
        let pos_tol = (tolerance.linear() * 10.0).max(1e-9);
        let surface = brep
            .faces
            .get(face)
            .and_then(|f| f.surface.as_ref())
            .map(|s| s.geometry.clone())
            .ok_or(Reason::StubUnresolved)?;
        use popkern_core::geometry::SurfaceEval;
        let (pt0, pt1) = pcurve.parameter_range();
        for k in 0..=8 {
            let t = pt0 + (pt1 - pt0) * k as f64 / 8.0;
            let uvp = pcurve.point_at(t);
            let p3 = surface.point_at(uvp.x, uvp.y);
            let f = k as f64 / 8.0;
            let expect = samples3d[((f * (n - 1) as f64).round() as usize).min(n - 1)];
            if (p3 - expect).length() > pos_tol.max(1e-7) {
                return Err(Reason::StubUnresolved);
            }
        }
        staged.push((face, heh, pcurve));
    }
    Ok(staged)
}

/// Re-author the rail edge after the near anchor moved along its own
/// line: same line curve, Linear pcurves refreshed on both flanks.
#[allow(clippy::too_many_arguments)]
pub(super) fn author_rail(
    brep: &mut BRep,
    rail: &AnchorEdge,
    anchor: VertexHandle,
    far_pos: Point3,
    w: Point3,
    strip_on: FaceHandle,
    cyl_on: &StripCylinder,
    plane: &SupportPlane,
    tolerance: &Tolerance,
) -> Result<(), Reason> {
    let e = brep.edges.get(rail.edge).ok_or(Reason::StubUnresolved)?;
    let primary = e.half_edge;
    let he = brep.half_edges.get(primary).ok_or(Reason::StubUnresolved)?;
    let tw = brep.half_edges.get(he.twin).ok_or(Reason::StubUnresolved)?;
    let primary_from_anchor = he.origin == anchor;
    let (p_start, p_end) = if primary_from_anchor {
        (w, far_pos)
    } else {
        (far_pos, w)
    };
    let line = CurveOnEdge::line_segment(p_start, p_end, tolerance.linear())
        .ok_or(Reason::StubUnresolved)?;
    let mut staged: Vec<(FaceHandle, HalfEdgeHandle, Pcurve)> = Vec::with_capacity(2);
    for (heh, hh) in [(primary, he), (he.twin, tw)] {
        let face = brep
            .loops
            .get(hh.loop_)
            .map(|l| l.face)
            .ok_or(Reason::StubUnresolved)?;
        let (from3, to3) = if hh.origin == anchor {
            (w, far_pos)
        } else {
            (far_pos, w)
        };
        let uv_of = |p: Point3| -> (f64, f64) {
            if face == strip_on {
                let d = p - cyl_on.origin;
                let radial = d - cyl_on.axis * d.dot(&cyl_on.axis);
                (
                    radial
                        .dot(&orthonormal_basis(cyl_on.axis).1)
                        .atan2(radial.dot(&orthonormal_basis(cyl_on.axis).0)),
                    d.dot(&cyl_on.axis),
                )
            } else {
                let d = p - plane.origin;
                (d.dot(&plane.u_axis), d.dot(&plane.v_axis))
            }
        };
        let mut from_uv = uv_of(from3);
        let to_uv = uv_of(to3);
        if face == strip_on {
            // Rail is an axial generator: both endpoints share one
            // chart angle; unwrap the moved end onto the far end's
            // branch and pin it exactly (iso-u by construction).
            let mut u = from_uv.0;
            while u - to_uv.0 > std::f64::consts::PI {
                u -= std::f64::consts::TAU;
            }
            while to_uv.0 - u > std::f64::consts::PI {
                u += std::f64::consts::TAU;
            }
            if (u - to_uv.0).abs() * cyl_on.radius > (tolerance.linear() * 100.0).max(1e-6) {
                return Err(Reason::StubUnresolved);
            }
            from_uv = (to_uv.0, from_uv.1);
        }
        staged.push((face, heh, Pcurve::linear(from_uv, to_uv, 0.0, 1.0)));
    }
    commit_edge_authoring(brep, rail.edge, line, staged);
    Ok(())
}
