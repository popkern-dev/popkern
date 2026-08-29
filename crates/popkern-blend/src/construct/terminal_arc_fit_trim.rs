//! Fitted terminal-arc trim for SWEPT (NURBS) blend strips — the
//! variable-radius counterpart of the PK-148 terminal arc trim.
//!
//! A variable-radius sweep's chain end keeps the support-edge stubs
//! running through the original corner, `(√2−1)·r_end` off the strip
//! (measured 4.2e-2 on the r=0.10 end of a 0.02→0.10 cube fillet —
//! the dominant welded-open source for variable blends). Unlike the
//! constant case there is no analytic strip cylinder; here the strip
//! SURFACE and its stub pcurves are the construction authority (the
//! sweep interpolates its circular sections exactly), so:
//!
//! 1. sample the strip surface along both stub pcurve spans,
//! 2. fit the terminal circle in the planar cap face's chart
//!    (algebraic fit, residual-gated — a non-circular section
//!    declines honestly),
//! 3. settle the corner onto the circle, re-author both stub 3-D
//!    curves as exact arcs and the cap-side pcurves as ellipticals,
//!    and KEEP the strip-side pcurves untouched (they are already
//!    the correct trim in the sweep chart).

use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{CurveEval, CurveGeometry, CurveOnEdge, Point3, SurfaceEval, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle, HalfEdgeHandle};

use crate::types::TerminalArcTrimDeclineReason;

use super::terminal_arc_author::StagedStubArc;
use super::terminal_arc_trim::{CornerCandidate, incident_edge_count};
use super::terminal_chord_trim::StripPlane;
use super::terminal_corner_recognition::strip_pcurve_span;

pub(super) fn repair_corner_swept(
    brep: &mut BRep,
    c: &CornerCandidate,
    tolerance: &Tolerance,
) -> Result<(), TerminalArcTrimDeclineReason> {
    let trace = popkern_core::trace::trace_enabled("terminal-arc");
    let strip_surface = match brep
        .faces
        .get(c.strip_face)
        .and_then(|f| f.surface.as_ref())
        .map(|s| s.geometry.clone())
    {
        Some(s @ SurfaceGeometry::Nurbs(_)) => s,
        _ => return Err(TerminalArcTrimDeclineReason::NonCylinderStrip),
    };

    // Both twins must bound one shared planar cap face.
    let twin_face = |he: HalfEdgeHandle| -> Option<FaceHandle> {
        let tw = brep.half_edges.get(brep.half_edges.get(he)?.twin)?;
        Some(brep.loops.get(tw.loop_)?.face)
    };
    let cap_face = twin_face(c.he_in).ok_or(TerminalArcTrimDeclineReason::NoCommonPlanarCap)?;
    if twin_face(c.he_out) != Some(cap_face) {
        if trace {
            eprintln!(
                "[tat]   swept stub twins bound different faces: in->{cap_face:?} out->{:?}",
                twin_face(c.he_out)
            );
        }
        return Err(TerminalArcTrimDeclineReason::NoCommonPlanarCap);
    }
    let cap = match brep
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
        }) => StripPlane {
            origin: *origin,
            normal: *normal,
            u_axis: *u_axis,
            v_axis: *v_axis,
        },
        other => {
            if trace {
                eprintln!(
                    "[tat]   swept cap face {cap_face:?} is not a Plane: {:?}",
                    other.map(std::mem::discriminant)
                );
            }
            return Err(TerminalArcTrimDeclineReason::NoCommonPlanarCap);
        }
    };
    if incident_edge_count(brep, c.corner) != 2 {
        return Err(TerminalArcTrimDeclineReason::VertexNotIsolated);
    }

    // Anchors (trusted vertex positions) — extracted first because
    // they also pin the terminal iso below.
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

    // Sample the strip surface along the TRUE terminal iso — the
    // sweep chart is the authority for the terminal section, but the
    // stale corner's chart shadow bends the stub pcurves slightly off
    // the iso (measured 1.4e-3 in u), so the pcurve spans only seed
    // the iso choice. The iso COORDINATE itself is refined from the
    // trusted anchors' converged chart projections: after an in-place
    // surface swap the spans are re-authored by generic projection
    // and carry enough jitter to pull the sampled section off the
    // terminal circle (measured residual 6.4e-4 on an exact-section
    // sweep, PK-179).
    let (a_start, a_end) = strip_pcurve_span(brep, c.strip_face, c.he_in)
        .ok_or(TerminalArcTrimDeclineReason::ArcBranchUnresolved)?;
    let (_, b_end) = strip_pcurve_span(brep, c.strip_face, c.he_out)
        .ok_or(TerminalArcTrimDeclineReason::ArcBranchUnresolved)?;
    let iso_is_u = (a_start.0 - b_end.0).abs() <= (a_start.1 - b_end.1).abs();
    let mut fixed = if iso_is_u {
        0.5 * (a_start.0 + b_end.0)
    } else {
        0.5 * (a_start.1 + b_end.1)
    };
    let span_extent = if iso_is_u {
        (a_start.1 - b_end.1).abs()
    } else {
        (a_start.0 - b_end.0).abs()
    };
    if let (Some(pi), Some(po)) = (
        strip_surface.closest_point(&anchor_in),
        strip_surface.closest_point(&anchor_out),
    ) {
        let (fi, fo) = if iso_is_u { (pi.0, po.0) } else { (pi.1, po.1) };
        // Accept the refinement only when both projections agree on
        // one iso near the seeded value — a wrong-branch Gauss-Newton
        // landing must not hijack the section.
        if (fi - fo).abs() <= 0.05 * span_extent.max(1e-12)
            && (0.5 * (fi + fo) - fixed).abs() <= 0.25 * span_extent.max(1e-12)
        {
            fixed = 0.5 * (fi + fo);
        }
    }
    let free_of = |uv: (f64, f64)| if iso_is_u { uv.1 } else { uv.0 };
    let at_iso = |free: f64| -> Point3 {
        if iso_is_u {
            strip_surface.point_at(fixed, free)
        } else {
            strip_surface.point_at(free, fixed)
        }
    };
    let mut samples: Vec<Point3> = Vec::new();
    for (lo, hi) in [
        (free_of(a_start), free_of(a_end)),
        (free_of(a_end), free_of(b_end)),
    ] {
        for k in 0..=6 {
            let t = k as f64 / 6.0;
            samples.push(at_iso(lo + (hi - lo) * t));
        }
    }
    let pos_tol = (tolerance.linear() * 10.0).max(1e-9);
    // All samples must lie in the cap plane (the terminal section is
    // the cap trim).
    let uv_in_cap = |p: Point3| -> (f64, f64) {
        let d = p - cap.origin;
        (d.dot(&cap.u_axis), d.dot(&cap.v_axis))
    };
    for p in &samples {
        if (*p - cap.origin).dot(&cap.normal).abs() > pos_tol * 10.0 {
            if trace {
                eprintln!(
                    "[tat]   swept sample {p:?} off cap plane by {:.3e}",
                    (*p - cap.origin).dot(&cap.normal).abs()
                );
            }
            return Err(TerminalArcTrimDeclineReason::NoCommonPlanarCap);
        }
    }

    // Anchors (trusted vertex positions).
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

    // Terminal circle EXACTLY through the two anchors and the strip
    // surface's point at the stubs' shared chart corner (the deepest
    // section point) — anchored construction, so the authored arcs
    // hit the anchors by definition. The sampled points then serve
    // only to verify the section really is circular.
    let deep = at_iso(free_of(a_end));
    let (cu, cv, radius) =
        circumcircle_2d(uv_in_cap(anchor_in), uv_in_cap(deep), uv_in_cap(anchor_out)).ok_or_else(
            || {
                if trace {
                    eprintln!("[tat]   swept anchors/deep collinear — no terminal circle");
                }
                TerminalArcTrimDeclineReason::ArcBranchUnresolved
            },
        )?;
    let mut residual = 0.0_f64;
    for p in &samples {
        let (u, v) = uv_in_cap(*p);
        residual = residual.max((((u - cu).powi(2) + (v - cv).powi(2)).sqrt() - radius).abs());
    }
    if residual > (pos_tol * 10.0).max(radius * 1.0e-3) {
        // PK-189: a non-circular terminal section is the conic-profile
        // strip class. Fit the terminal conic (the stale corner IS the
        // tangent-triangle apex on plane supports), settle the corner
        // onto its shoulder, and author the stubs as exact rational
        // halves. A failed fit keeps the honest decline.
        if let Some(conic) = super::terminal_conic_fit_trim::fit_terminal_conic(
            brep, c.corner, anchor_in, anchor_out, &samples, tolerance,
        ) {
            let staged = [
                super::terminal_conic_fit_trim::compute_swept_stub_conic(
                    brep, c.he_in, &conic, true, &cap, tolerance,
                )?,
                super::terminal_conic_fit_trim::compute_swept_stub_conic(
                    brep, c.he_out, &conic, false, &cap, tolerance,
                )?,
            ];
            if trace {
                eprintln!(
                    "[tat]   swept terminal section is a conic (circular residual {residual:.3e}); settling on the fitted shoulder"
                );
            }
            let shoulder = conic.shoulder;
            if let Some(v) = brep.vertices.get_mut(c.corner) {
                v.position = shoulder;
            }
            for stub in staged {
                if let Some(e) = brep.edges.get_mut(stub.edge) {
                    e.curve = Some(stub.curve);
                }
                if let Some(tw) = brep.half_edges.get_mut(stub.cap_half_edge) {
                    tw.pcurve = Some(stub.cap_pcurve.clone());
                }
                brep.set_face_use_pcurve(cap_face, stub.cap_half_edge, stub.cap_pcurve);
            }
            return Ok(());
        }
        if trace {
            eprintln!("[tat]   swept terminal section not circular: residual {residual:.3e}");
        }
        return Err(TerminalArcTrimDeclineReason::ArcBranchUnresolved);
    }
    let center = cap.origin + cap.u_axis * cu + cap.v_axis * cv;

    // Settled corner: the circle point nearest the stale corner.
    let old_corner = brep
        .vertices
        .get(c.corner)
        .ok_or(TerminalArcTrimDeclineReason::ArcBranchUnresolved)?
        .position;
    let radial = old_corner - center;
    let radial = radial - cap.normal * radial.dot(&cap.normal);
    let new_corner = center
        + radial
            .normalized()
            .ok_or(TerminalArcTrimDeclineReason::ArcBranchUnresolved)?
            * radius;
    if (new_corner - old_corner).length() > radius {
        if trace {
            eprintln!(
                "[tat]   swept corner move {:.3e} exceeds radius {radius:.3e}",
                (new_corner - old_corner).length()
            );
        }
        return Err(TerminalArcTrimDeclineReason::ArcBranchUnresolved);
    }

    // Author both stub arcs; strip-side pcurves stay as-is (they are
    // the sweep-chart authority that located the terminal section).
    let mut staged: Vec<StagedStubArc> = Vec::with_capacity(2);
    for (he_strip, start_pos, end_pos) in [
        (c.he_in, anchor_in, new_corner),
        (c.he_out, new_corner, anchor_out),
    ] {
        staged.push(compute_swept_stub_arc(
            brep, he_strip, center, radius, &cap, start_pos, end_pos, tolerance,
        )?);
    }
    if let Some(v) = brep.vertices.get_mut(c.corner) {
        v.position = new_corner;
    }
    for stub in staged {
        if let Some(e) = brep.edges.get_mut(stub.edge) {
            e.curve = Some(stub.curve);
        }
        if let Some(tw) = brep.half_edges.get_mut(stub.cap_half_edge) {
            tw.pcurve = Some(stub.cap_pcurve.clone());
        }
        brep.set_face_use_pcurve(cap_face, stub.cap_half_edge, stub.cap_pcurve);
    }
    Ok(())
}

/// Circumcircle through three 2-D points; returns (cu, cv, r).
fn circumcircle_2d(a: (f64, f64), b: (f64, f64), c: (f64, f64)) -> Option<(f64, f64, f64)> {
    let d = 2.0 * (a.0 * (b.1 - c.1) + b.0 * (c.1 - a.1) + c.0 * (a.1 - b.1));
    if d.abs() < 1e-30 {
        return None;
    }
    let a2 = a.0 * a.0 + a.1 * a.1;
    let b2 = b.0 * b.0 + b.1 * b.1;
    let c2 = c.0 * c.0 + c.1 * c.1;
    let cu = (a2 * (b.1 - c.1) + b2 * (c.1 - a.1) + c2 * (a.1 - b.1)) / d;
    let cv = (a2 * (c.0 - b.0) + b2 * (a.0 - c.0) + c2 * (b.0 - a.0)) / d;
    let r = ((a.0 - cu).powi(2) + (a.1 - cv).powi(2)).sqrt();
    Some((cu, cv, r))
}

/// Author one swept-strip stub as an exact arc of the fitted terminal
/// circle plus the cap-side elliptical pcurve. The strip-side pcurve
/// is not produced — the existing sweep-chart pcurve stays.
#[allow(clippy::too_many_arguments)]
fn compute_swept_stub_arc(
    brep: &BRep,
    he_strip: HalfEdgeHandle,
    center: Point3,
    radius: f64,
    cap: &StripPlane,
    start_pos: Point3,
    end_pos: Point3,
    tolerance: &Tolerance,
) -> Result<StagedStubArc, TerminalArcTrimDeclineReason> {
    let trace = popkern_core::trace::trace_enabled("terminal-arc");
    let (he_twin, edge, primary, strip_pcurve) = {
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
        let strip_face = brep
            .loops
            .get(he.loop_)
            .map(|l| l.face)
            .ok_or(TerminalArcTrimDeclineReason::ArcBranchUnresolved)?;
        let pc = brep
            .face_use_pcurve(strip_face, he_strip)
            .or(he.pcurve.as_ref())
            .cloned()
            .ok_or(TerminalArcTrimDeclineReason::ArcBranchUnresolved)?;
        (he.twin, edge, primary, pc)
    };
    let (p_start, p_end) = if primary == he_strip {
        (start_pos, end_pos)
    } else {
        (end_pos, start_pos)
    };

    let axis = cap.normal;
    let u_basis = p_start - center;
    let mut v_basis = axis.cross(&u_basis);
    let d = p_end - center;
    let mut raw = d.dot(&v_basis).atan2(d.dot(&u_basis));
    if raw < 0.0 {
        v_basis = -v_basis;
        raw = -raw;
    }
    if !(1e-9..=std::f64::consts::PI - 1e-9).contains(&raw) {
        if trace {
            eprintln!("[tat]   swept {he_strip:?}: sweep {raw:.4} outside minor-arc range");
        }
        return Err(TerminalArcTrimDeclineReason::ArcBranchUnresolved);
    }
    let arc = CurveOnEdge {
        geometry: CurveGeometry::Ellipse {
            center,
            u_basis,
            v_basis,
            start_angle: 0.0,
            end_angle: raw,
        },
        t_start: 0.0,
        t_end: raw,
    };
    let pos_tol = (tolerance.linear() * 10.0).max(1e-9);
    let hit_start = (arc.geometry.point_at(arc.t_start) - p_start).length();
    let hit_end = (arc.geometry.point_at(arc.t_end) - p_end).length();
    if hit_start > pos_tol || hit_end > pos_tol {
        if trace {
            eprintln!(
                "[tat]   swept {he_strip:?}: arc endpoint miss {hit_start:.3e}/{hit_end:.3e}"
            );
        }
        return Err(TerminalArcTrimDeclineReason::ArcBranchUnresolved);
    }

    // Cap-side pcurve: elliptical (circular) about the fitted centre.
    let uv_of = |p: Point3| -> (f64, f64) {
        let d = p - cap.origin;
        (d.dot(&cap.u_axis), d.dot(&cap.v_axis))
    };
    let center_uv = uv_of(center);
    let (cap_from, cap_to) = (end_pos, start_pos);
    let angle_of = |p: Point3| -> f64 {
        let (u, v) = uv_of(p);
        (v - center_uv.1).atan2(u - center_uv.0)
    };
    let cap_start = angle_of(cap_from);
    let mut cap_end = angle_of(cap_to);
    while cap_end - cap_start > std::f64::consts::PI {
        cap_end -= std::f64::consts::TAU;
    }
    while cap_start - cap_end > std::f64::consts::PI {
        cap_end += std::f64::consts::TAU;
    }
    let cap_pcurve =
        Pcurve::elliptical(center_uv, radius, radius, 0.0, cap_start, cap_end, 0.0, 1.0);
    let plane_point =
        |uv: (f64, f64)| -> Point3 { cap.origin + cap.u_axis * uv.0 + cap.v_axis * uv.1 };
    let cap_p0 = plane_point({
        let p = cap_pcurve.point_at(0.0);
        (p.x, p.y)
    });
    let cap_p1 = plane_point({
        let p = cap_pcurve.point_at(1.0);
        (p.x, p.y)
    });
    if (cap_p0 - cap_from).length() > pos_tol || (cap_p1 - cap_to).length() > pos_tol {
        if trace {
            eprintln!(
                "[tat]   swept {he_strip:?}: cap pcurve endpoint miss {:.3e}/{:.3e}",
                (cap_p0 - cap_from).length(),
                (cap_p1 - cap_to).length()
            );
        }
        return Err(TerminalArcTrimDeclineReason::ArcBranchUnresolved);
    }

    Ok(StagedStubArc {
        edge,
        curve: arc,
        strip_half_edge: he_strip,
        strip_pcurve,
        cap_half_edge: he_twin,
        cap_pcurve,
    })
}

const _: fn() = || {
    let _ = Vec3::new;
};
