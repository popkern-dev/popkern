//! Terminal stub settle for CONE strips (circular chamfers on
//! interrupted rims).
//!
//! `repair_corner` has arms for Cylinder (terminal arc), Plane
//! (terminal chord), and Nurbs (swept fit); Cone strips declined as
//! `NonCylinderStrip`, so an interrupted circular rim's stale corner
//! vertex and its two stub edges survived untrimmed (measured: the
//! notched flange's rim chamfer shipped a top-plane stub bounding
//! Plane/Cone with a 1.0 curve-vs-cone deviation — check-9 catch at
//! serving). The terminal boundary of a cone strip against a planar
//! cap is the cone∩plane conic; both stubs are settled onto it,
//! split at the conic point nearest the stale corner.

use popkern_core::geometry::pcurve::{Pcurve, PcurveGeometry};
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::uv_projection::project_point_to_uv;
use popkern_core::geometry::{CurveEval, CurveGeometry, CurveOnEdge, Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle, HalfEdgeHandle};
use popkern_math::PortableMath;

use crate::types::TerminalArcTrimDeclineReason;

use super::terminal_arc_author::StagedStubArc;
use super::terminal_arc_trim::{CornerCandidate, incident_edge_count};

/// Cone strip parameters (already validated by the dispatcher).
pub(super) struct StripCone {
    pub(super) apex: Point3,
    pub(super) axis: Vec3,
    pub(super) half_angle: f64,
}

/// Samples per stub for the conic interpolation. The samples are
/// closed-form exact on both the cone and the cap plane; a cubic
/// interpolation through 17 of them carries error far below the
/// serving tolerance for terminal-scale spans.
const STUB_SAMPLES: usize = 17;

pub(super) fn repair_corner_conic(
    brep: &mut BRep,
    c: &CornerCandidate,
    cone: &StripCone,
    tolerance: &Tolerance,
) -> Result<(), TerminalArcTrimDeclineReason> {
    let axis = cone
        .axis
        .normalized()
        .ok_or(TerminalArcTrimDeclineReason::NonCylinderStrip)?;
    let tan_h = cone.half_angle.ptan();
    if tan_h.abs() < 1e-12 {
        return Err(TerminalArcTrimDeclineReason::NonCylinderStrip);
    }

    // Both twins must bound one shared planar cap face.
    let twin_face = |he: HalfEdgeHandle| -> Option<FaceHandle> {
        let tw = brep.half_edges.get(brep.half_edges.get(he)?.twin)?;
        let l = brep.loops.get(tw.loop_)?;
        Some(l.face)
    };
    let cap_face = twin_face(c.he_in).ok_or(TerminalArcTrimDeclineReason::NoCommonPlanarCap)?;
    if twin_face(c.he_out) != Some(cap_face) {
        return Err(TerminalArcTrimDeclineReason::NoCommonPlanarCap);
    }
    let (cap_origin, cap_normal, cap_u, cap_v) = match brep
        .faces
        .get(cap_face)
        .and_then(|f| f.surface.as_ref())
        .map(|s| &s.geometry)
    {
        Some(SurfaceGeometry::Plane {
            origin,
            normal,
            u_axis,
            v_axis,
        }) => (*origin, *normal, *u_axis, *v_axis),
        _ => return Err(TerminalArcTrimDeclineReason::NoCommonPlanarCap),
    };
    let cap_n = cap_normal
        .normalized()
        .ok_or(TerminalArcTrimDeclineReason::NoCommonPlanarCap)?;

    // The corner vertex must be used by exactly the two stub edges.
    if incident_edge_count(brep, c.corner) != 2 {
        return Err(TerminalArcTrimDeclineReason::VertexNotIsolated);
    }

    // Trusted geometric inputs: the two anchor vertex positions.
    let anchor_of = |he: HalfEdgeHandle, at_origin: bool| -> Option<Point3> {
        let data = brep.half_edges.get(he)?;
        let vertex = if at_origin {
            data.origin
        } else {
            brep.half_edges.get(data.next)?.origin
        };
        Some(brep.vertices.get(vertex)?.position)
    };
    let anchor_in =
        anchor_of(c.he_in, true).ok_or(TerminalArcTrimDeclineReason::ArcBranchUnresolved)?;
    let anchor_out =
        anchor_of(c.he_out, false).ok_or(TerminalArcTrimDeclineReason::ArcBranchUnresolved)?;

    // Anchors are trusted TOPOLOGY but their positions carry
    // terminal-construction noise (measured: the notched flange's
    // triple points sat 2.1e-2 off the cone). They are projected
    // onto the conic below; only gross misplacement declines.
    let pos_tol = (tolerance.linear() * 10.0).max(1e-9);
    let trace = popkern_core::trace::trace_enabled("terminal-arc");
    let span = anchor_in.distance_to(&anchor_out).max(pos_tol);
    let anchor_band = span * 0.2;

    // Conic sampling between the anchors: at each axial station the
    // cone is a circle; intersect it with the cap plane and keep the
    // branch continuous with the chord. Exact per sample.
    let s_in = (anchor_in - cone.apex).dot(&axis);
    let s_out = (anchor_out - cone.apex).dot(&axis);
    let conic_at = |f: f64| -> Option<Point3> {
        let station = s_in + (s_out - s_in) * f;
        let center = cone.apex + axis * station;
        let radius = station.abs() * tan_h.abs();
        // Circle plane ∩ cap plane: line through the circle's plane.
        // Solve directly in the circle's basis.
        let (bu, bv) = {
            let any = if axis.x.abs() < 0.9 {
                Vec3::new(1.0, 0.0, 0.0)
            } else {
                Vec3::new(0.0, 1.0, 0.0)
            };
            let bu = axis.cross(&any).normalized()?;
            (bu, axis.cross(&bu))
        };
        // Points: center + r(cosθ·bu + sinθ·bv), constrained to the
        // cap plane: (P − cap_origin)·n = 0 →
        // a·cosθ + b·sinθ = c with a = r·bu·n, b = r·bv·n.
        let a = radius * bu.dot(&cap_n);
        let b = radius * bv.dot(&cap_n);
        let rhs = (cap_origin - center).dot(&cap_n);
        let amp = (a * a + b * b).sqrt();
        if amp < 1e-300 || rhs.abs() > amp {
            return None;
        }
        let base = b.patan2(a);
        let delta = (rhs / amp).pacos();
        let chord_ref = anchor_in + (anchor_out - anchor_in) * f;
        let mut best: Option<(f64, Point3)> = None;
        for theta in [base + delta, base - delta] {
            let p = center + bu * (radius * theta.pcos()) + bv * (radius * theta.psin());
            let d = p.distance_to(&chord_ref);
            if best.map(|(bd, _)| d < bd).unwrap_or(true) {
                best = Some((d, p));
            }
        }
        best.map(|(_, p)| p)
    };

    // Project a point onto the conic parameterization (search over
    // an extended window so slightly-off endpoint anchors still find
    // their true triple points).
    let project_f = |target: Point3, lo: f64, hi: f64| -> Option<(f64, Point3, f64)> {
        let mut best: Option<(f64, Point3, f64)> = None;
        for k in 0..=96 {
            let f = lo + (hi - lo) * k as f64 / 96.0;
            let Some(p) = conic_at(f) else { continue };
            let d = p.distance_to(&target);
            if best.map(|(_, _, bd)| d < bd).unwrap_or(true) {
                best = Some((f, p, d));
            }
        }
        best
    };
    // The true anchor positions are the conic points AT the raw
    // anchors' own axial stations (f = 0 and f = 1 by the
    // parameterization) — the cone ∩ cap ∩ adjacent-support triple
    // points. Projecting by 3D distance instead slides the anchor
    // down the cone off its support plane (measured: 8.3e-3 off the
    // top plane, which dragged the notch edge into a 2.1e-3 lie).
    let (f_in, f_out) = (0.0, 1.0);
    let anchor_in_on = conic_at(f_in).ok_or(TerminalArcTrimDeclineReason::ArcBranchUnresolved)?;
    let anchor_out_on = conic_at(f_out).ok_or(TerminalArcTrimDeclineReason::ArcBranchUnresolved)?;
    let (d_in, d_out) = (
        anchor_in_on.distance_to(&anchor_in),
        anchor_out_on.distance_to(&anchor_out),
    );
    if d_in > anchor_band || d_out > anchor_band {
        if trace {
            eprintln!(
                "[tat]   conic anchors off by {d_in:.3e}/{d_out:.3e} (band {anchor_band:.3e})"
            );
        }
        return Err(TerminalArcTrimDeclineReason::AnchorsOffSurface);
    }
    // Settled corner: conic point nearest the stale corner.
    let old_corner = brep
        .vertices
        .get(c.corner)
        .ok_or(TerminalArcTrimDeclineReason::ArcBranchUnresolved)?
        .position;
    let (best_f, new_corner, _) = project_f(old_corner, f_in, f_out)
        .ok_or(TerminalArcTrimDeclineReason::ArcBranchUnresolved)?;
    // Sanity: the settle stays local.
    if (new_corner - old_corner).length() > anchor_in.distance_to(&anchor_out).max(pos_tol) {
        return Err(TerminalArcTrimDeclineReason::ArcBranchUnresolved);
    }

    // Author both stubs as interpolated conic segments.
    let cone_surface = SurfaceGeometry::Cone {
        apex: cone.apex,
        axis: cone.axis,
        half_angle: cone.half_angle,
    };
    let mut staged: Vec<StagedStubArc> = Vec::with_capacity(2);
    for (he_strip, f_a, f_b) in [(c.he_in, f_in, best_f), (c.he_out, best_f, f_out)] {
        staged.push(author_conic_stub(
            brep,
            he_strip,
            &conic_at,
            f_a,
            f_b,
            &cone_surface,
            (cap_origin, cap_u, cap_v),
            tolerance,
        )?);
    }
    if let Some(v) = brep.vertices.get_mut(c.corner) {
        v.position = new_corner;
    }
    // Snap the anchor vertices onto their true triple points too —
    // the terminal construction placed them off-cone by the same
    // noise the settle exists to repair.
    let anchor_vertex = |he: HalfEdgeHandle, at_origin: bool| -> Option<_> {
        let data = brep.half_edges.get(he)?;
        Some(if at_origin {
            data.origin
        } else {
            brep.half_edges.get(data.next)?.origin
        })
    };
    for (vh, target) in [
        (anchor_vertex(c.he_in, true), anchor_in_on),
        (anchor_vertex(c.he_out, false), anchor_out_on),
    ] {
        let Some(vh) = vh else { continue };
        if let Some(v) = brep.vertices.get_mut(vh) {
            v.position = target;
        }
        // Every incident edge curve must follow the moved vertex or
        // the adjacent loops open by the snap distance.
        super::topology_repair::snap_incident_edge_curve_endpoints_robust(brep, vh, target);
    }
    for stub in staged {
        if let Some(e) = brep.edges.get_mut(stub.edge) {
            e.curve = Some(stub.curve);
        }
        if let Some(he) = brep.half_edges.get_mut(stub.strip_half_edge) {
            he.pcurve = Some(stub.strip_pcurve.clone());
        }
        brep.set_face_use_pcurve(c.strip_face, stub.strip_half_edge, stub.strip_pcurve);
        if let Some(tw) = brep.half_edges.get_mut(stub.cap_half_edge) {
            tw.pcurve = Some(stub.cap_pcurve.clone());
        }
        brep.set_face_use_pcurve(cap_face, stub.cap_half_edge, stub.cap_pcurve);
    }
    Ok(())
}

/// Interpolate one stub's conic segment plus strip/cap NURBS pcurves.
#[allow(clippy::too_many_arguments)]
fn author_conic_stub(
    brep: &BRep,
    he_strip: HalfEdgeHandle,
    conic_at: &dyn Fn(f64) -> Option<Point3>,
    f_a: f64,
    f_b: f64,
    cone_surface: &SurfaceGeometry,
    cap_frame: (Point3, Vec3, Vec3),
    tolerance: &Tolerance,
) -> Result<StagedStubArc, TerminalArcTrimDeclineReason> {
    let (he_twin, edge, primary) = {
        let he = brep
            .half_edges
            .get(he_strip)
            .ok_or(TerminalArcTrimDeclineReason::ArcBranchUnresolved)?;
        let edge = he.edge;
        let primary = brep
            .edges
            .get(edge)
            .ok_or(TerminalArcTrimDeclineReason::ArcBranchUnresolved)?
            .half_edge;
        (he.twin, edge, primary)
    };
    // The edge's stored curve is oriented along its primary half-edge.
    let (g_a, g_b) = if primary == he_strip {
        (f_a, f_b)
    } else {
        (f_b, f_a)
    };

    let mut pts_3d = Vec::with_capacity(STUB_SAMPLES);
    let mut pts_strip_uv = Vec::with_capacity(STUB_SAMPLES);
    let mut pts_cap_uv = Vec::with_capacity(STUB_SAMPLES);
    let mut params = Vec::with_capacity(STUB_SAMPLES);
    let mut prev_u: Option<f64> = None;
    for i in 0..STUB_SAMPLES {
        let t = i as f64 / (STUB_SAMPLES - 1) as f64;
        let f = g_a + (g_b - g_a) * t;
        let p = conic_at(f).ok_or(TerminalArcTrimDeclineReason::ArcBranchUnresolved)?;
        pts_3d.push(p);
        params.push(t);
        let (mut u, v) =
            project_point_to_uv(cone_surface, &p).unwrap_or_else(|| (prev_u.unwrap_or(0.0), 0.0));
        if let Some(pu) = prev_u {
            while u - pu > std::f64::consts::PI {
                u -= std::f64::consts::TAU;
            }
            while pu - u > std::f64::consts::PI {
                u += std::f64::consts::TAU;
            }
        }
        prev_u = Some(u);
        pts_strip_uv.push(Point3::new(u, v, 0.0));
        let (co, cu, cv) = cap_frame;
        let d = p - co;
        pts_cap_uv.push(Point3::new(d.dot(&cu), d.dot(&cv), 0.0));
    }

    let fit3 = popkern_core::geometry::fitting::interpolate_curve(&pts_3d, &params, 3)
        .map_err(|_| TerminalArcTrimDeclineReason::ArcBranchUnresolved)?;
    let fit_strip = popkern_core::geometry::fitting::interpolate_curve(&pts_strip_uv, &params, 3)
        .map_err(|_| TerminalArcTrimDeclineReason::ArcBranchUnresolved)?;
    let fit_cap = popkern_core::geometry::fitting::interpolate_curve(&pts_cap_uv, &params, 3)
        .map_err(|_| TerminalArcTrimDeclineReason::ArcBranchUnresolved)?;

    let curve = CurveOnEdge {
        geometry: CurveGeometry::Nurbs(fit3),
        t_start: 0.0,
        t_end: 1.0,
    };
    // Postcondition: the authored curve must land on both surfaces
    // along its whole span (this is exactly what serving-boundary
    // check 9 verifies later — fail HERE, transactionally, instead).
    let pos_tol = (tolerance.linear() * 10.0).max(1e-9);
    for i in 0..=8 {
        let t = i as f64 / 8.0;
        let p = curve.geometry.point_at(t);
        let on_cone = project_point_to_uv(cone_surface, &p)
            .map(|(u, v)| cone_surface.point_at(u, v).distance_to(&p))
            .unwrap_or(f64::INFINITY);
        if on_cone > pos_tol * 10.0 {
            return Err(TerminalArcTrimDeclineReason::ArcBranchUnresolved);
        }
    }

    let (strip_t0, strip_t1, cap_t0, cap_t1) = (0.0, 1.0, 0.0, 1.0);
    Ok(StagedStubArc {
        edge,
        curve,
        strip_half_edge: he_strip,
        strip_pcurve: Pcurve {
            geometry: PcurveGeometry::Nurbs(fit_strip),
            t_start: strip_t0,
            t_end: strip_t1,
        },
        cap_half_edge: he_twin,
        cap_pcurve: Pcurve {
            geometry: PcurveGeometry::Nurbs(fit_cap),
            t_start: cap_t0,
            t_end: cap_t1,
        },
    })
}
