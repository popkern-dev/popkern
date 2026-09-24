//! Exact boundaries for caps inserted at free (open-sheet) fillet
//! ends (PK-532).
//!
//! `insert_cap_face` splits the strip's end off with `mef` between
//! the two contacts. The new face gets the cap surface and is bounded
//! by the `mef` edge (strip | cap) and the rim L: two stubs through
//! the original corner. The `mef` edge had no curve, which check 9
//! silently skipped. The stubs kept pcurves in the strip's chart,
//! although they now bound the cap.
//!
//! **Strip | cap edge.** It is the strip's section by the plane
//! through the contacts, normal to the outward spine tangent:
//! - cylinder strips get the exact arc (`terminal_arc_oblique`);
//! - fitted strips get a marched strip ∩ plane curve
//!   (`surface_pair_march`).
//!
//! A flat cap is that plane, and a rounded cap's sphere (centre on the
//! axis at the station, radius r) meets the cylinder in the same
//! circle.
//!
//! **Rim.** For each cap kind:
//! - **Flat:** the L must lie in the cap plane (true when the rim is
//!   perpendicular to the spine). The stubs stay lines with plane-chart
//!   pcurves.
//! - **Rounded:** the quarter ball's rim is two meridians, from each
//!   contact to the outward pole, so the corner moves to the pole.
//! - **TangentMatch:** its bilinear patch has a straight contact-to-
//!   contact boundary and cannot contain the section arc, so it
//!   declines.
//!
//! **Fallback.** A decline hands the reason back. The caller restores
//! the pre-`mef` body and reports `CapInsertionDeclined`, and the end
//! is then settled open on its rim section. **Tolerance:** the
//! caller's blend tolerance.

use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{CurveEval, CurveGeometry, CurveOnEdge, Point3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle, HalfEdgeHandle, VertexHandle};

use crate::cap::CapGeometry;
use crate::types::CapKind;

use super::blend_seam_corner_march::authoring::chart_pcurve;
use super::blend_seam_corner_swept::project_seeded;
use super::surface_pair_march::{intersection_curve, intersection_pcurve, march_intersection};
use super::terminal_arc_oblique::{
    ObliqueSection, compute_oblique_stub_arc, section_arc, unwrap_near,
};
use super::terminal_arc_trim::StripCylinder;

/// Samples for a projected pcurve.
const PCURVE_SAMPLES: usize = 17;

fn face_of(brep: &BRep, he: HalfEdgeHandle) -> Option<FaceHandle> {
    Some(brep.loops.get(brep.half_edges.get(he)?.loop_)?.face)
}

fn surface_of(brep: &BRep, face: FaceHandle) -> Option<SurfaceGeometry> {
    brep.faces
        .get(face)?
        .surface
        .as_ref()
        .map(|s| s.geometry.clone())
}

fn position(brep: &BRep, v: VertexHandle) -> Option<Point3> {
    brep.vertices.get(v).map(|x| x.position)
}

/// Author the inserted cap's boundary. `cap_half_edge` is the `mef`
/// half-edge in the cap face's loop; `corner` is the original corner.
pub(super) fn author_cap_boundary(
    brep: &mut BRep,
    cap_half_edge: HalfEdgeHandle,
    corner: VertexHandle,
    cap: &CapGeometry,
    tolerance: &Tolerance,
) -> Result<(), String> {
    let pos_tol = (tolerance.linear() * 10.0).max(1e-9);
    let he_cap = brep
        .half_edges
        .get(cap_half_edge)
        .ok_or("cap half-edge missing")?;
    let (he_strip, cap_face) = (
        he_cap.twin,
        face_of(brep, cap_half_edge).ok_or("cap face missing")?,
    );
    let strip_face = face_of(brep, he_strip).ok_or("strip face missing")?;
    let strip = surface_of(brep, strip_face).ok_or("strip surface missing")?;
    let cap_surface = surface_of(brep, cap_face).ok_or("cap surface missing")?;
    // Cap loop: he_cap (to → from), then from → corner, corner → to.
    let s1 = he_cap.next;
    let s2 = brep.half_edges.get(s1).ok_or("cap loop broken")?.next;
    if brep.half_edges.get(s2).map(|h| h.next) != Some(cap_half_edge)
        || brep.half_edges.get(s2).map(|h| h.origin) != Some(corner)
    {
        return Err("the cap face is not bounded by the strip end and a two-stub rim".into());
    }
    let from = position(brep, brep.half_edges[he_strip].origin).ok_or("contact missing")?;
    let to = position(brep, he_cap.origin).ok_or("contact missing")?;
    let stale = position(brep, corner).ok_or("corner missing")?;
    let normal = cap
        .spine_tangent_outward
        .normalized()
        .ok_or("degenerate spine tangent")?;
    let mid = from + (to - from) * 0.5;

    if cap.kind == CapKind::Flat && (stale - mid).dot(&normal).abs() > pos_tol {
        return Err("the rim does not lie in the flat cap's plane".into());
    }

    // The strip | cap edge: the strip's end section.
    let section_plane = SurfaceGeometry::Plane {
        origin: mid,
        normal,
        u_axis: (to - from).normalized().ok_or("contacts coincide")?,
        v_axis: normal.cross(&(to - from).normalized().ok_or("contacts coincide")?),
    };
    let (edge, curve, strip_pcurve, cap_pcurve) = match (&strip, &cap_surface) {
        (
            SurfaceGeometry::Cylinder {
                origin,
                axis,
                radius,
            },
            _,
        ) => {
            let cyl = StripCylinder {
                origin: *origin,
                axis: *axis,
                radius: *radius,
            };
            let section = ObliqueSection::new(&cyl, mid, normal).ok_or("section degenerate")?;
            if section.miss(from) > pos_tol || section.miss(to) > pos_tol {
                return Err("contacts are off the strip's end section".into());
            }
            let (chart_u, _) =
                popkern_core::geometry::uv_projection::project_point_to_uv(&strip, &from)
                    .ok_or("contact does not project onto the strip")?;
            let u_from = unwrap_near(section.angle_of(from), chart_u);
            let u_to = unwrap_near(section.angle_of(to), u_from);
            match &cap_surface {
                SurfaceGeometry::Plane {
                    origin: po,
                    u_axis,
                    v_axis,
                    ..
                } => {
                    let staged = compute_oblique_stub_arc(
                        brep,
                        he_strip,
                        &section,
                        u_from,
                        u_to,
                        (*po, *u_axis, *v_axis),
                        tolerance,
                    )
                    .map_err(|e| format!("end section arc: {e:?}"))?;
                    (
                        staged.edge,
                        staged.curve,
                        staged.strip_pcurve,
                        staged.cap_pcurve,
                    )
                }
                _ => {
                    let arc = section_arc(brep, he_strip, &section, u_from, u_to, tolerance)
                        .map_err(|e| format!("end section arc: {e:?}"))?;
                    // The curve follows the edge's primary half-edge;
                    // the cap pcurve follows the cap half-edge.
                    let primary = brep.edges.get(arc.edge).ok_or("edge missing")?.half_edge;
                    let cap_pcurve =
                        projected_pcurve(&cap_surface, &arc.curve, primary != cap_half_edge)?;
                    (arc.edge, arc.curve, arc.strip_pcurve, cap_pcurve)
                }
            }
        }
        _ => {
            let marched = march_intersection(&strip, &section_plane, from, to, pos_tol, pos_tol)
                .map_err(|e| format!("end section march: {e:?}"))?;
            let edge = brep.half_edges[he_strip].edge;
            let primary = brep.edges.get(edge).ok_or("edge missing")?.half_edge;
            let curve = intersection_curve(&marched, primary != he_strip)
                .map_err(|e| format!("end section curve: {e:?}"))?;
            let strip_pcurve = intersection_pcurve(&strip, &marched, false, None)
                .map_err(|e| format!("end section strip pcurve: {e:?}"))?;
            let cap_pcurve = intersection_pcurve(&cap_surface, &marched, true, None)
                .map_err(|e| format!("end section cap pcurve: {e:?}"))?;
            (edge, curve, strip_pcurve, cap_pcurve)
        }
    };
    commit(
        brep,
        edge,
        curve,
        &[
            (strip_face, he_strip, strip_pcurve),
            (cap_face, cap_half_edge, cap_pcurve),
        ],
    );

    // The rim.
    match (&cap.kind, &cap_surface) {
        (CapKind::Flat, SurfaceGeometry::Plane { .. }) => {
            for (stub, a, b) in [(s1, from, stale), (s2, stale, to)] {
                rim_line(brep, stub, cap_face, &cap_surface, a, b, tolerance)?;
            }
        }
        (CapKind::Rounded, SurfaceGeometry::Sphere { center, radius }) => {
            let pole = *center + normal * *radius;
            if let Some(v) = brep.vertices.get_mut(corner) {
                v.position = pole;
            }
            for (stub, contact) in [(s1, from), (s2, to)] {
                meridian(brep, stub, cap_face, &cap_surface, *center, contact, pole)?;
            }
        }
        (kind, _) => return Err(format!("{kind:?} caps have no exact rim construction")),
    }
    for he in [cap_half_edge, s1, s2] {
        let e = brep.half_edges[he].edge;
        if brep.measure_edge_deviation(e) > pos_tol {
            return Err(format!(
                "cap boundary edge {e:?} leaves its faces by {:.3e}",
                brep.measure_edge_deviation(e)
            ));
        }
    }
    Ok(())
}

fn commit(
    brep: &mut BRep,
    edge: popkern_core::topology::EdgeHandle,
    curve: CurveOnEdge,
    pcurves: &[(FaceHandle, HalfEdgeHandle, Pcurve)],
) {
    if let Some(e) = brep.edges.get_mut(edge) {
        e.curve = Some(curve);
    }
    for (face, he, pcurve) in pcurves {
        if let Some(h) = brep.half_edges.get_mut(*he) {
            h.pcurve = Some(pcurve.clone());
        }
        brep.set_face_use_pcurve(*face, *he, pcurve.clone());
    }
}

/// A flat cap's rim stub: the line from `a` to `b` (the half-edge's
/// direction) with its plane-chart pcurve.
fn rim_line(
    brep: &mut BRep,
    stub: HalfEdgeHandle,
    cap_face: FaceHandle,
    plane: &SurfaceGeometry,
    a: Point3,
    b: Point3,
    tolerance: &Tolerance,
) -> Result<(), String> {
    let SurfaceGeometry::Plane {
        origin,
        u_axis,
        v_axis,
        ..
    } = plane
    else {
        return Err("flat cap surface is not a plane".into());
    };
    let edge = brep.half_edges[stub].edge;
    let primary = brep.edges.get(edge).ok_or("edge missing")?.half_edge;
    let (p0, p1) = if primary == stub { (a, b) } else { (b, a) };
    let line =
        CurveOnEdge::line_segment(p0, p1, tolerance.linear()).ok_or("degenerate rim stub")?;
    let uv = |p: Point3| ((p - *origin).dot(u_axis), (p - *origin).dot(v_axis));
    commit(
        brep,
        edge,
        line,
        &[(cap_face, stub, Pcurve::linear(uv(a), uv(b), 0.0, 1.0))],
    );
    Ok(())
}

/// A rounded cap's rim stub: the quarter meridian from `contact` to
/// `pole` on the sphere about `center`, with its projected pcurve.
fn meridian(
    brep: &mut BRep,
    stub: HalfEdgeHandle,
    cap_face: FaceHandle,
    sphere: &SurfaceGeometry,
    center: Point3,
    contact: Point3,
    pole: Point3,
) -> Result<(), String> {
    let edge = brep.half_edges[stub].edge;
    let primary = brep.edges.get(edge).ok_or("edge missing")?.half_edge;
    // The stub runs contact → pole along `stub`'s direction iff it
    // leaves the contact; the curve follows the primary half-edge.
    let stub_from_contact = brep
        .vertices
        .get(brep.half_edges[stub].origin)
        .is_some_and(|v| (v.position - contact).length() < 1e-12);
    let forward = (primary == stub) == stub_from_contact;
    let (u_basis, v_basis) = if forward {
        (contact - center, pole - center)
    } else {
        (pole - center, contact - center)
    };
    let curve = CurveOnEdge {
        geometry: CurveGeometry::Ellipse {
            center,
            u_basis,
            v_basis,
            start_angle: 0.0,
            end_angle: std::f64::consts::FRAC_PI_2,
        },
        t_start: 0.0,
        t_end: std::f64::consts::FRAC_PI_2,
    };
    let pcurve = projected_pcurve(sphere, &curve, primary != stub)?;
    commit(brep, edge, curve, &[(cap_face, stub, pcurve)]);
    Ok(())
}

/// The pcurve of `curve` in `surface`'s chart, sampled at uniform
/// parameter fractions, unwrapped onto each predecessor's period
/// branch, and interpolated. `reversed` builds it against the
/// curve's direction (for the half-edge opposite the primary one).
pub(super) fn projected_pcurve(
    surface: &SurfaceGeometry,
    curve: &CurveOnEdge,
    reversed: bool,
) -> Result<Pcurve, String> {
    projected_pcurve_sampled(surface, curve, reversed, PCURVE_SAMPLES)
}

/// [`projected_pcurve`] with an explicit sample count, for callers
/// that densify until their deviation gate holds.
pub(super) fn projected_pcurve_sampled(
    surface: &SurfaceGeometry,
    curve: &CurveOnEdge,
    reversed: bool,
    samples: usize,
) -> Result<Pcurve, String> {
    use popkern_core::geometry::SurfaceEval;
    let samples = samples.max(2);
    let (period_u, period_v) = surface.parameter_periodicity();
    let unwrap = |value: f64, near: f64, period: Option<f64>| match period {
        Some(t) if t > 0.0 => value - ((value - near) / t).round() * t,
        _ => value,
    };
    let mut uv: Vec<(f64, f64)> = Vec::with_capacity(samples);
    let mut params: Vec<f64> = Vec::with_capacity(samples);
    let mut hint: Option<(f64, f64)> = None;
    for k in 0..samples {
        let f = k as f64 / (samples - 1) as f64;
        let f_curve = if reversed { 1.0 - f } else { f };
        let p = curve
            .geometry
            .point_at(curve.t_start + (curve.t_end - curve.t_start) * f_curve);
        let (mut u, mut v) = project_seeded(surface, &p, hint).ok_or("sample does not project")?;
        if let Some((hu, hv)) = hint {
            u = unwrap(u, hu, period_u);
            v = unwrap(v, hv, period_v);
        }
        hint = Some((u, v));
        uv.push((u, v));
        params.push(f);
    }
    // A sample at a chart pole (a sphere's, where ∂/∂u vanishes) has
    // no meaningful u: a meridian keeps its neighbour's u there. The
    // outward pole of a rounded cap is such a point whenever the spine
    // is along the sphere chart's axis (measured 9.4e-4 pcurve miss
    // before this).
    let degenerate = |(u, v): (f64, f64)| surface.partial_u(u, v).length() < 1e-9;
    for k in 0..uv.len() {
        if degenerate(uv[k]) {
            let neighbour = if k > 0 {
                uv[k - 1].0
            } else {
                uv.get(1).map_or(uv[k].0, |n| n.0)
            };
            uv[k].0 = neighbour;
        }
    }
    chart_pcurve(&uv, &params).map_err(|e| format!("pcurve fit: {e:?}"))
}
