//! PK-148 terminal-arc stub authoring (see `terminal_arc_trim`).
//!
//! Pure computation of one stub edge's exact minor-arc 3-D curve,
//! strip-side `Linear` pcurve, and cap-side `Elliptical` pcurve from
//! trusted geometric inputs. Nothing here mutates the BRep — the
//! caller commits staged results transactionally.

use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::uv_projection::project_point_to_uv;
use popkern_core::geometry::{CurveEval, CurveGeometry, CurveOnEdge, Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, HalfEdgeHandle};

use crate::types::TerminalArcTrimDeclineReason;

use super::terminal_arc_trim::StripCylinder;
use super::terminal_chord_trim::StripPlane;

/// One computed-but-uncommitted stub re-authoring.
pub(super) struct StagedStubArc {
    pub(super) edge: EdgeHandle,
    pub(super) curve: CurveOnEdge,
    pub(super) strip_half_edge: HalfEdgeHandle,
    pub(super) strip_pcurve: Pcurve,
    pub(super) cap_half_edge: HalfEdgeHandle,
    pub(super) cap_pcurve: Pcurve,
}

/// Author one stub edge as an exact minor arc from `start_pos` to
/// `end_pos` (positions given in the half-edge's own direction),
/// plus the strip-side `Linear` pcurve and the cap-side `Elliptical`
/// pcurve. Pure computation — nothing is committed here.
#[allow(clippy::too_many_arguments)]
pub(super) fn compute_stub_arc(
    brep: &BRep,
    he_strip: HalfEdgeHandle,
    cyl: &StripCylinder,
    center: Point3,
    start_pos: Point3,
    end_pos: Point3,
    plane: (Point3, Vec3, Vec3),
    tolerance: &Tolerance,
) -> Result<StagedStubArc, TerminalArcTrimDeclineReason> {
    let trace = std::env::var_os("PK_TERMINAL_ARC_TRACE").is_some();
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
    let (p_start, p_end) = if primary == he_strip {
        (start_pos, end_pos)
    } else {
        (end_pos, start_pos)
    };

    // Exact minor arc with explicit basis: the terminal arc at a
    // chain-end corner always subtends less than π, so the branch
    // is the minor one by construction.
    let u_basis = p_start - center;
    let mut v_basis = cyl.axis.cross(&u_basis);
    let d = p_end - center;
    let mut raw = d.dot(&v_basis).atan2(d.dot(&u_basis));
    if raw < 0.0 {
        v_basis = -v_basis;
        raw = -raw;
    }
    if !(1e-9..=std::f64::consts::PI - 1e-9).contains(&raw) {
        if trace {
            eprintln!("[tat]   {he_strip:?}: sweep {raw:.4} outside minor-arc range");
        }
        return Err(TerminalArcTrimDeclineReason::ArcBranchUnresolved);
    }
    let end_angle = raw;
    let arc = CurveOnEdge {
        geometry: CurveGeometry::Ellipse {
            center,
            u_basis,
            v_basis,
            start_angle: 0.0,
            end_angle,
        },
        t_start: 0.0,
        t_end: end_angle,
    };
    // Postcondition before committing: the authored curve must hit
    // both endpoint positions.
    let hit_start = (arc.geometry.point_at(arc.t_start) - p_start).length();
    let hit_end = (arc.geometry.point_at(arc.t_end) - p_end).length();
    let pos_tol = (tolerance.linear() * 10.0).max(1e-9);
    if hit_start > pos_tol || hit_end > pos_tol {
        if trace {
            eprintln!(
                "[tat]   {he_strip:?}: postcondition miss start {hit_start:.3e} end {hit_end:.3e}"
            );
        }
        return Err(TerminalArcTrimDeclineReason::ArcBranchUnresolved);
    }

    // Strip-side pcurve: exact iso-v Linear span in the CURRENT
    // chart, projected from the trusted endpoint positions (the old
    // pcurve may live in a stale chart). Unwrap the end u to the
    // period branch nearest the start.
    let (su, sv) = project_point_to_uv(
        &SurfaceGeometry::Cylinder {
            origin: cyl.origin,
            axis: cyl.axis,
            radius: cyl.radius,
        },
        &start_pos,
    )
    .ok_or(TerminalArcTrimDeclineReason::ArcBranchUnresolved)?;
    let (eu_raw, ev) = project_point_to_uv(
        &SurfaceGeometry::Cylinder {
            origin: cyl.origin,
            axis: cyl.axis,
            radius: cyl.radius,
        },
        &end_pos,
    )
    .ok_or(TerminalArcTrimDeclineReason::ArcBranchUnresolved)?;
    let mut eu = eu_raw;
    while eu - su > std::f64::consts::PI {
        eu -= std::f64::consts::TAU;
    }
    while su - eu > std::f64::consts::PI {
        eu += std::f64::consts::TAU;
    }
    let v_chart = 0.5 * (sv + ev);
    // Strip pcurve is stored per half-edge (origin → destination);
    // start_pos/end_pos already follow he_strip, so no reversal.
    let strip_pcurve = Pcurve::linear((su, v_chart), (eu, v_chart), 0.0, 1.0);

    // Cap-side pcurve: the arc lies in the cap plane, so its UV
    // image is an axis-aligned circle of radius r (unit plane axes).
    let (plane_origin, plane_u, plane_v) = plane;
    let uv_of = |p: Point3| -> (f64, f64) {
        let d = p - plane_origin;
        (d.dot(&plane_u), d.dot(&plane_v))
    };
    let center_uv = uv_of(center);
    // Twin runs opposite to the strip half-edge.
    let (cap_from, cap_to) = (end_pos, start_pos);
    let angle_of = |p: Point3| -> f64 {
        let (u, v) = uv_of(p);
        (v - center_uv.1).atan2(u - center_uv.0)
    };
    let cap_start = angle_of(cap_from);
    let mut cap_end = angle_of(cap_to);
    // Minor-arc branch: fold the sweep into (−π, π).
    while cap_end - cap_start > std::f64::consts::PI {
        cap_end -= std::f64::consts::TAU;
    }
    while cap_start - cap_end > std::f64::consts::PI {
        cap_end += std::f64::consts::TAU;
    }
    let cap_pcurve = Pcurve::elliptical(
        center_uv, cyl.radius, cyl.radius, 0.0, cap_start, cap_end, 0.0, 1.0,
    );
    // Postcondition: the cap pcurve's UV endpoints map back onto the
    // 3-D endpoints through the plane.
    let plane_point = |uv: (f64, f64)| -> Point3 { plane_origin + plane_u * uv.0 + plane_v * uv.1 };
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
                "[tat]   {he_strip:?}: cap pcurve endpoint miss {:.3e}/{:.3e}",
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

/// Author one stub edge as an exact line segment from `start_pos` to
/// `end_pos` (positions given in the half-edge's own direction), plus
/// `Linear` pcurves on both the planar strip side and the planar cap
/// side — the chamfer counterpart of [`compute_stub_arc`]. Pure
/// computation; nothing is committed here.
pub(super) fn compute_stub_chord(
    brep: &BRep,
    he_strip: HalfEdgeHandle,
    strip_plane: &StripPlane,
    cap_plane: &StripPlane,
    start_pos: Point3,
    end_pos: Point3,
    tolerance: &Tolerance,
) -> Result<StagedStubArc, TerminalArcTrimDeclineReason> {
    let trace = std::env::var_os("PK_TERMINAL_ARC_TRACE").is_some();
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
    let (p_start, p_end) = if primary == he_strip {
        (start_pos, end_pos)
    } else {
        (end_pos, start_pos)
    };
    let line = CurveOnEdge::line_segment(p_start, p_end, tolerance.linear())
        .ok_or(TerminalArcTrimDeclineReason::ArcBranchUnresolved)?;
    let pos_tol = (tolerance.linear() * 10.0).max(1e-9);
    let hit_start = (line.geometry.point_at(line.t_start) - p_start).length();
    let hit_end = (line.geometry.point_at(line.t_end) - p_end).length();
    if hit_start > pos_tol || hit_end > pos_tol {
        if trace {
            eprintln!(
                "[tat]   {he_strip:?}: chord postcondition miss start {hit_start:.3e} end \
                 {hit_end:.3e}"
            );
        }
        return Err(TerminalArcTrimDeclineReason::ArcBranchUnresolved);
    }

    // Both pcurves are exact linear spans in unit plane frames; each
    // follows its own half-edge's direction (twin runs opposite).
    let uv_in = |plane: &StripPlane, p: Point3| -> (f64, f64) {
        let d = p - plane.origin;
        (d.dot(&plane.u_axis), d.dot(&plane.v_axis))
    };
    let strip_pcurve = Pcurve::linear(
        uv_in(strip_plane, start_pos),
        uv_in(strip_plane, end_pos),
        0.0,
        1.0,
    );
    let cap_pcurve = Pcurve::linear(
        uv_in(cap_plane, end_pos),
        uv_in(cap_plane, start_pos),
        0.0,
        1.0,
    );
    // Postcondition: both pcurves map back onto the 3-D endpoints
    // through their planes.
    let check = |plane: &StripPlane, pc: &Pcurve, from: Point3, to: Point3| -> f64 {
        let at = |t: f64| -> Point3 {
            let uv = pc.point_at(t);
            plane.origin + plane.u_axis * uv.x + plane.v_axis * uv.y
        };
        (at(0.0) - from).length().max((at(1.0) - to).length())
    };
    let strip_miss = check(strip_plane, &strip_pcurve, start_pos, end_pos);
    let cap_miss = check(cap_plane, &cap_pcurve, end_pos, start_pos);
    if strip_miss > pos_tol || cap_miss > pos_tol {
        if trace {
            eprintln!(
                "[tat]   {he_strip:?}: chord pcurve endpoint miss strip {strip_miss:.3e} cap \
                 {cap_miss:.3e}"
            );
        }
        return Err(TerminalArcTrimDeclineReason::ArcBranchUnresolved);
    }

    Ok(StagedStubArc {
        edge,
        curve: line,
        strip_half_edge: he_strip,
        strip_pcurve,
        cap_half_edge: he_twin,
        cap_pcurve,
    })
}
