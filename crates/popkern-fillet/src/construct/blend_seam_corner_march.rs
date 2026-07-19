//! PK-152 stage 2: mismatched-radius pair-corner seam via the
//! marching intersector.
//!
//! With unequal radii the strip×strip intersection leaves the
//! bisector-plane family (a space quartic), and the corner topology
//! differs structurally from the equal-radius case: the near anchor
//! sits on ONE strip only, joined to the larger strip's terminal
//! circle by a Δr mismatch stub on the shared support edge. The true
//! closure:
//!
//! - the seam runs from the far anchor (rail crossing, on both
//!   strips) to the point W where the intersection curve crosses the
//!   stub's support plane — W lies on the smaller strip's rail line
//!   by construction (rail = strip ∩ support plane);
//! - the near anchor MOVES to W (its rail edge shortens along its
//!   own line, no topology change);
//! - the mismatch stub becomes the larger strip's exact
//!   cylinder×plane ellipse arc from W to its retained endpoint.
//!
//! The march is exact per sample: parameterized by the near strip's
//! chart angle, the far strip's constraint is a quadratic in the
//! axial coordinate, so there is no integrator drift; the exit is
//! bisected onto the support plane. The seam 3-D curve and chart
//! pcurves interpolate the marched samples (degree 3, endpoints
//! pinned); stub and rail re-authoring is analytic. Same
//! transactional/decline contract as stage 1.
//!
//! Submodules: `solver` (the exact march), `authoring` (chart
//! pcurves, stub ellipse, rail re-author).

use popkern_core::geometry::CurveEval;
use popkern_core::geometry::fitting::interpolate_curve;
use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::geometry::{CurveGeometry, CurveOnEdge, Point3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle, VertexHandle};

use crate::types::BlendSeamCornerDeclineReason as Reason;

use super::blend_seam_corner::{PairCornerCandidate, off_cylinder, strip_cylinder};
use super::blend_seam_corner_repair::{
    commit_edge_authoring, neighbour_chart_point, perform_flap_surgery,
};

#[path = "blend_seam_corner_march_authoring.rs"]
mod authoring;
#[path = "blend_seam_corner_march_solver.rs"]
mod solver;

use authoring::{
    author_rail, author_stub_pcurves, chart_pcurve, chart_polyline, chord_params,
    cylinder_plane_arc, params_reversed,
};
use solver::{face_plane, march_seam};

/// One edge incident to the near anchor, classified by its flanks.
struct AnchorEdge {
    edge: EdgeHandle,
    /// The flank that is NOT a strip of this corner.
    support: FaceHandle,
    /// The endpoint that is not the near anchor.
    far_vertex: VertexHandle,
}

/// Mismatched-radius pair-corner repair (stage 2).
pub(super) fn repair_pair_corner_mixed(
    brep: &mut BRep,
    c: &PairCornerCandidate,
    tolerance: &Tolerance,
) -> Result<(), Reason> {
    let cyl_n = strip_cylinder(brep, c.strip_n).ok_or(Reason::EulerRejected)?;
    let cyl_p = strip_cylinder(brep, c.strip_p).ok_or(Reason::EulerRejected)?;
    let pos_tol = (tolerance.linear() * 10.0).max(1e-9);
    let pos_of =
        |b: &BRep, v: VertexHandle| -> Option<Point3> { Some(b.vertices.get(v)?.position) };
    let a_pos = pos_of(brep, c.anchor_near).ok_or(Reason::EulerRejected)?;
    let b_pos = pos_of(brep, c.anchor_far).ok_or(Reason::EulerRejected)?;

    // The near anchor must sit on exactly one strip.
    let on_n = off_cylinder(a_pos, &cyl_n) <= pos_tol;
    let on_p = off_cylinder(a_pos, &cyl_p) <= pos_tol;
    let (strip_on, cyl_on, strip_off, cyl_off) = match (on_n, on_p) {
        (true, false) => (c.strip_n, cyl_n.clone(), c.strip_p, cyl_p.clone()),
        (false, true) => (c.strip_p, cyl_p.clone(), c.strip_n, cyl_n.clone()),
        _ => return Err(Reason::StubUnresolved),
    };

    // Classify the near anchor's other incident edges: the rail
    // (strip_on × support) and the Δr mismatch stub (strip_off ×
    // support), sharing one support face.
    let mut rail: Option<AnchorEdge> = None;
    let mut stub: Option<AnchorEdge> = None;
    let mut extra = 0usize;
    for (eh, e) in brep.edges.iter() {
        if eh == c.seam_edge {
            continue;
        }
        let Some(he) = brep.half_edges.get(e.half_edge) else {
            continue;
        };
        let Some(tw) = brep.half_edges.get(he.twin) else {
            continue;
        };
        if he.origin != c.anchor_near && tw.origin != c.anchor_near {
            continue;
        }
        let far_vertex = if he.origin == c.anchor_near {
            tw.origin
        } else {
            he.origin
        };
        let (Some(lf), Some(lt)) = (brep.loops.get(he.loop_), brep.loops.get(tw.loop_)) else {
            continue;
        };
        let flanks = [lf.face, lt.face];
        let entry = |support: FaceHandle| AnchorEdge {
            edge: eh,
            support,
            far_vertex,
        };
        if flanks.contains(&strip_on) && !flanks.contains(&strip_off) {
            let support = if flanks[0] == strip_on {
                flanks[1]
            } else {
                flanks[0]
            };
            rail = Some(entry(support));
        } else if flanks.contains(&strip_off) && !flanks.contains(&strip_on) {
            let support = if flanks[0] == strip_off {
                flanks[1]
            } else {
                flanks[0]
            };
            stub = Some(entry(support));
        } else {
            extra += 1;
        }
    }
    let rail = rail.ok_or(Reason::StubUnresolved)?;
    let stub = stub.ok_or(Reason::StubUnresolved)?;
    if extra != 0 || rail.support != stub.support {
        return Err(Reason::StubUnresolved);
    }
    let plane = face_plane(brep, stub.support).ok_or(Reason::StubUnresolved)?;
    let stub_far_pos = pos_of(brep, stub.far_vertex).ok_or(Reason::StubUnresolved)?;
    if off_cylinder(stub_far_pos, &cyl_off) > pos_tol
        || (stub_far_pos - plane.origin).dot(&plane.normal).abs() > pos_tol
    {
        return Err(Reason::StubUnresolved);
    }
    let rail_far_pos = pos_of(brep, rail.far_vertex).ok_or(Reason::StubUnresolved)?;

    // March the seam from the far anchor to the support-plane exit.
    let corner_pos = pos_of(brep, c.corner).unwrap_or(a_pos);
    let pts = march_seam(&cyl_on, &cyl_off, b_pos, corner_pos, &plane, tolerance)?;
    let w = *pts.last().ok_or(Reason::MarchUnresolved)?;
    // W replaces the near anchor; sanity: it stays within the corner
    // neighbourhood.
    if (w - a_pos).length() > cyl_off.radius + cyl_on.radius {
        return Err(Reason::MarchUnresolved);
    }

    // Author everything BEFORE surgery where possible.
    let params = chord_params(&pts);
    let seam_nurbs = interpolate_curve(&pts, &params, 3).map_err(|_| Reason::MarchUnresolved)?;
    let primary_seam = brep
        .edges
        .get(c.seam_edge)
        .ok_or(Reason::EulerRejected)?
        .half_edge;
    // After surgery he_n runs anchor_near(=W) → anchor_far; the
    // marched points run B → W.
    let seam_curve = {
        let (t0, t1) = seam_nurbs.parameter_range();
        if primary_seam == c.he_n {
            // Curve must run W → B: reverse the interpolation input.
            let mut rev = pts.clone();
            rev.reverse();
            let curve = interpolate_curve(&rev, &params_reversed(&params), 3)
                .map_err(|_| Reason::MarchUnresolved)?;
            let (r0, r1) = curve.parameter_range();
            CurveOnEdge {
                geometry: CurveGeometry::Nurbs(curve),
                t_start: r0,
                t_end: r1,
            }
        } else {
            CurveOnEdge {
                geometry: CurveGeometry::Nurbs(seam_nurbs.clone()),
                t_start: t0,
                t_end: t1,
            }
        }
    };

    // Stub arc: strip_off × support plane, W → retained endpoint.
    let (stub_arc_from, stub_arc_to) = (w, stub_far_pos);
    let primary_stub = brep
        .edges
        .get(stub.edge)
        .ok_or(Reason::StubUnresolved)?
        .half_edge;
    let stub_primary_he = brep
        .half_edges
        .get(primary_stub)
        .ok_or(Reason::StubUnresolved)?;
    let stub_primary_from_anchor = stub_primary_he.origin == c.anchor_near;
    let (sp_from, sp_to) = if stub_primary_from_anchor {
        (stub_arc_from, stub_arc_to)
    } else {
        (stub_arc_to, stub_arc_from)
    };
    let (stub_curve, _, _) = cylinder_plane_arc(&cyl_off, &plane, sp_from, sp_to, tolerance)?;

    perform_flap_surgery(brep, c)?;
    if let Some(v) = brep.vertices.get_mut(c.anchor_near) {
        v.position = w;
    }

    // Seam pcurves: chart polylines anchored at the far-anchor end.
    let seam_staged = {
        let mut staged: Vec<(FaceHandle, HalfEdgeHandle, Pcurve)> = Vec::with_capacity(2);
        for (face, cyl, he) in [(c.strip_n, &cyl_n, c.he_n), (c.strip_p, &cyl_p, c.he_p)] {
            let h = brep.half_edges.get(he).ok_or(Reason::EulerRejected)?;
            // The far anchor is he_n's destination and he_p's origin.
            let b_end_anchor = if he == c.he_n {
                neighbour_chart_point(brep, face, h.next, false)
            } else {
                neighbour_chart_point(brep, face, h.prev, true)
            };
            let uv = chart_polyline(cyl, &pts, b_end_anchor.map(|(u, _)| u));
            if let Some((au, av)) = b_end_anchor {
                let (u0c, v0c) = uv[0];
                let chart_tol = (tolerance.linear() * 100.0).max(1e-6);
                if ((u0c - au).abs() * cyl.radius).max((v0c - av).abs()) > chart_tol {
                    return Err(Reason::PcurveClosureFailed);
                }
            }
            // Orient along the half-edge: he_n runs W→B (reverse of
            // march order), he_p runs B→W (march order).
            let pcurve = if he == c.he_n {
                let mut rev = uv.clone();
                rev.reverse();
                chart_pcurve(&rev, &params_reversed(&params))?
            } else {
                chart_pcurve(&uv, &params)?
            };
            staged.push((face, he, pcurve));
        }
        staged
    };
    commit_edge_authoring(brep, c.seam_edge, seam_curve, seam_staged);

    // Stub: strip_off-side sinusoid + support-side pcurve.
    let stub_staged = author_stub_pcurves(
        brep,
        &stub,
        &stub_curve,
        strip_off,
        &cyl_off,
        &plane,
        tolerance,
    )?;
    commit_edge_authoring(brep, stub.edge, stub_curve, stub_staged);

    // Rail: same line, shortened to W.
    author_rail(
        brep,
        &rail,
        c.anchor_near,
        rail_far_pos,
        w,
        strip_on,
        &cyl_on,
        &plane,
        tolerance,
    )?;

    // Postconditions.
    for eh in [c.seam_edge, stub.edge, rail.edge] {
        if brep.measure_edge_deviation(eh) > (tolerance.linear() * 100.0).max(1e-7) {
            return Err(Reason::ValidationFailed);
        }
    }
    brep.validate().map_err(|_| Reason::ValidationFailed)?;
    Ok(())
}
