//! Different-adjacent-plane repair path for off-cylinder chord pairs.

use super::topology_repair_geometry::{
    fit_cylinder_uv_pcurve, project_curve_into_plane_uv, remap_to_unit_range, reverse_nurbs,
    slice_curve,
};
use super::topology_repair_snap::snap_incident_edge_curve_endpoints_robust;
use super::*;

/// Surgery for the different-adjacent-planes case: the off-cyl
/// vertex `cand.off_cyl_vertex` sits at the intersection of cyl,
/// plane_A, and plane_B (the three surfaces meeting at this corner).
/// Currently drifted off-cyl by a few % of radius; snap it onto the
/// analytic intersection of `cyl ∩ plane_A` and `cyl ∩ plane_B`,
/// then rewrite each chord edge's 3D curve to the proper analytic
/// arc on its respective plane.
///
/// No topology change — only vertex position + edge curve rewrite.
pub(super) fn try_repair_chord_pair_different_planes(
    brep: &mut BRep,
    cyl_face: FaceHandle,
    cand: ChordPairCandidate,
    plane_face_a: FaceHandle,
    plane_face_b: FaceHandle,
    trace: bool,
) -> bool {
    let Some(cyl_geom) = brep
        .faces
        .get(cyl_face)
        .and_then(|f| f.surface.as_ref())
        .map(|s| s.geometry.clone())
    else {
        if trace {
            eprintln!("[topo_repair]   reject (diff): cylinder face has no surface");
        }
        return false;
    };
    let plane_a_geom = match brep
        .faces
        .get(plane_face_a)
        .and_then(|f| f.surface.as_ref())
        .map(|s| s.geometry.clone())
    {
        Some(g) if matches!(g, SurfaceGeometry::Plane { .. }) => g,
        _ => {
            if trace {
                eprintln!("[topo_repair]   reject (diff): plane_face_a not a plane");
            }
            return false;
        }
    };
    let plane_b_geom = match brep
        .faces
        .get(plane_face_b)
        .and_then(|f| f.surface.as_ref())
        .map(|s| s.geometry.clone())
    {
        Some(g) if matches!(g, SurfaceGeometry::Plane { .. }) => g,
        _ => {
            if trace {
                eprintln!("[topo_repair]   reject (diff): plane_face_b not a plane");
            }
            return false;
        }
    };

    let v_off = cand.off_cyl_vertex;
    let p_off = brep.vertices[v_off].position;
    let p_anchor_in = brep.vertices[cand.anchor_in].position;
    let p_anchor_out = brep.vertices[cand.anchor_out].position;
    let tolerance = brep.tolerance;

    // #136 precondition: both anchors must lie ON the cylinder.
    // In pyramid apex-only fillets, cylinder face loops can have
    // CHAINS of off-cyl vertices that use each other as anchors;
    // surgery against off-surface anchors rewrites chord edges
    // between two off-surface points, and cumulative repairs
    // across the chain compound to ~1 m tess-nm divergence on
    // N=7/8.  See `tests/topo_repair_diff_planes_repro.rs::t_136_b3`.
    let (cyl_origin, cyl_axis, cyl_radius) = match &cyl_geom {
        SurfaceGeometry::Cylinder {
            origin,
            axis,
            radius,
        } => (*origin, *axis, *radius),
        _ => return false,
    };
    if !vertex_on_cylinder(&p_anchor_in, cyl_origin, cyl_axis, cyl_radius)
        || !vertex_on_cylinder(&p_anchor_out, cyl_origin, cyl_axis, cyl_radius)
    {
        if trace {
            eprintln!(
                "[topo_repair]   reject (diff): anchor off-cyl (chord chain) — anchor_in_on={}, anchor_out_on={}",
                vertex_on_cylinder(&p_anchor_in, cyl_origin, cyl_axis, cyl_radius),
                vertex_on_cylinder(&p_anchor_out, cyl_origin, cyl_axis, cyl_radius),
            );
        }
        return false;
    }

    // Compute cyl ∩ plane_A and cyl ∩ plane_B as analytic curves.
    let arc_a = match analytic_cyl_plane_intersection(&cyl_geom, &plane_a_geom, &tolerance) {
        Some(c) => c,
        None => {
            if trace {
                eprintln!("[topo_repair]   reject (diff): cyl ∩ plane_a failed");
            }
            return false;
        }
    };
    let arc_b = match analytic_cyl_plane_intersection(&cyl_geom, &plane_b_geom, &tolerance) {
        Some(c) => c,
        None => {
            if trace {
                eprintln!("[topo_repair]   reject (diff): cyl ∩ plane_b failed");
            }
            return false;
        }
    };

    // Find the point on arc_a closest to p_off (call it q_a) and
    // similarly q_b on arc_b. The "true" off-cyl vertex position is
    // arc_a ∩ arc_b; with both arcs near p_off (within drift), the
    // midpoint or the closer of the two projections is good enough
    // — the analytic curves themselves will hit p_off ± drift.
    let (t_a_off, q_a) = match closest_point_on_curve(&arc_a, &p_off) {
        Some((t, p)) => (t, p),
        None => return false,
    };
    let (t_b_off, q_b) = match closest_point_on_curve(&arc_b, &p_off) {
        Some((t, p)) => (t, p),
        None => return false,
    };

    // Sanity check — q_a and q_b should be within a few mm of each
    // other (both should approximate the true vertex, drift apart
    // only). If too far, this isn't the right intersection branch.
    let q_dist = q_a.distance_to(&q_b);
    let chord_scale = p_off
        .distance_to(&p_anchor_in)
        .max(p_off.distance_to(&p_anchor_out));
    if q_dist > chord_scale * 0.5 {
        if trace {
            eprintln!("[topo_repair]   reject (diff): q_a/q_b too far apart ({q_dist:.5})");
        }
        return false;
    }
    // Snap v_off to whichever arc projection is closer to the
    // current v_off — that side's chord becomes exact, and the other
    // side accepts a small residual mismatch (≤ q_dist).
    let d_a = q_a.distance_to(&p_off);
    let d_b = q_b.distance_to(&p_off);
    let v_new = if d_a <= d_b { q_a } else { q_b };
    let residual = q_dist;

    // Find the chord edges and their twin-side faces.
    let edge_in = brep.half_edges[cand.he_in].edge;
    let edge_out = brep.half_edges[cand.he_out].edge;
    // edge_in: anchor_in → v_off (on cyl), v_off → anchor_in (on plane_face_a)
    // edge_out: v_off → anchor_out (on cyl), anchor_out → v_off (on plane_face_b)

    // Move the off-cyl vertex onto the new position. Also update
    // every incident edge's 3D curve endpoint that was pointing at
    // the old position — otherwise those edges' endpoints diverge
    // from the new vertex by O(drift), producing ShellGap errors.
    if let Some(v) = brep.vertices.get_mut(v_off) {
        v.position = v_new;
    }
    snap_incident_edge_curve_endpoints_robust(brep, v_off, v_new);

    // Rewrite edge_in's 3D curve as cyl ∩ plane_A trimmed from
    // anchor_in to v_new (on plane_a side). Similarly edge_out.
    // Endpoint tolerance: residual drift between the two arcs +
    // anchor drift, with a small floor.
    let endpoint_tol = (residual + 1e-3).max(1e-3);
    if !rewrite_edge_to_arc(
        brep,
        edge_in,
        &arc_a,
        &cyl_geom,
        &plane_a_geom,
        p_anchor_in,
        v_new,
        endpoint_tol,
        trace,
        "in",
    ) {
        return false;
    }
    if !rewrite_edge_to_arc(
        brep,
        edge_out,
        &arc_b,
        &cyl_geom,
        &plane_b_geom,
        v_new,
        p_anchor_out,
        endpoint_tol,
        trace,
        "out",
    ) {
        return false;
    }

    let _ = (t_a_off, t_b_off);
    true
}

/// Walk the curve, find parameter t closest to `target`. Returns
/// (parameter, point_on_curve).
fn closest_point_on_curve(curve: &NurbsCurve, target: &Point3) -> Option<(f64, Point3)> {
    let t = curve.closest_point(target)?;
    let p = curve.point_at(t);
    Some((t, p))
}

/// Rewrite `edge`'s 3D curve to a slice of `arc` (cyl × plane
/// intersection) trimmed between `p_start` and `p_end`. Also update
/// the cyl-side and plane-side pcurves.
#[allow(clippy::too_many_arguments)]
fn rewrite_edge_to_arc(
    brep: &mut BRep,
    edge: EdgeHandle,
    arc: &NurbsCurve,
    cyl: &SurfaceGeometry,
    plane: &SurfaceGeometry,
    p_start: Point3,
    p_end: Point3,
    endpoint_tol: f64,
    trace: bool,
    label: &str,
) -> bool {
    let t_start = match arc.closest_point(&p_start) {
        Some(t) => t,
        None => {
            if trace {
                eprintln!("[topo_repair]   reject (diff/{label}): arc start projection failed");
            }
            return false;
        }
    };
    let t_end = match arc.closest_point(&p_end) {
        Some(t) => t,
        None => {
            if trace {
                eprintln!("[topo_repair]   reject (diff/{label}): arc end projection failed");
            }
            return false;
        }
    };
    let (lo, hi, reversed) = if t_start <= t_end {
        (t_start, t_end, false)
    } else {
        (t_end, t_start, true)
    };
    let sliced = match slice_curve(arc, lo, hi) {
        Some(c) => {
            if reversed {
                reverse_nurbs(&c)
            } else {
                c
            }
        }
        None => {
            if trace {
                eprintln!("[topo_repair]   reject (diff/{label}): slice_curve failed");
            }
            return false;
        }
    };
    let sliced = remap_to_unit_range(&sliced);
    // Verify the sliced arc actually starts at p_start and ends at p_end
    // (within the BRep's linear tolerance — the analytic intersection
    // should be exact, but the off-cyl-vertex move might leave residual
    // error of order anchor drift).
    let (s_min, s_max) = sliced.parameter_range();
    let arc_start = sliced.point_at(s_min);
    let arc_end = sliced.point_at(s_max);
    let endpoint_err = arc_start
        .distance_to(&p_start)
        .max(arc_end.distance_to(&p_end));
    if endpoint_err > endpoint_tol {
        if trace {
            eprintln!(
                "[topo_repair]   reject (diff/{label}): arc endpoint err {endpoint_err:.5} > tol {endpoint_tol:.5}",
            );
        }
        return false;
    }
    // Snap arc endpoints to the existing vertex positions so the
    // edge's curve_on_edge.point_at(t_start/t_end) matches vertex
    // positions exactly (validate's ShellGap check requires this).
    let mut sliced = sliced;
    if let Some(first) = sliced.control_points.first_mut() {
        *first = p_start;
    }
    if let Some(last) = sliced.control_points.last_mut() {
        *last = p_end;
    }

    // Re-project the sliced 3D arc to cyl UV (degree-3 fit) and
    // plane UV (analytic projection).
    let cyl_uv = match fit_cylinder_uv_pcurve(&sliced, cyl) {
        Some(c) => c,
        None => return false,
    };
    let (plane_origin, plane_u, plane_v) = match plane {
        SurfaceGeometry::Plane {
            origin,
            u_axis,
            v_axis,
            ..
        } => (*origin, *u_axis, *v_axis),
        _ => return false,
    };
    let plane_uv = match project_curve_into_plane_uv(
        &CurveGeometry::Nurbs(sliced.clone()),
        plane_origin,
        plane_u,
        plane_v,
    ) {
        Some(c) => c,
        None => return false,
    };

    // Apply: new 3D curve + pcurves on each half-edge. The cyl-side
    // half-edge of the edge goes "p_start → p_end" (after our
    // remapping the arc runs in this direction).
    let he = brep.edges[edge].half_edge;
    let twin = brep.half_edges[he].twin;
    let he_face = brep.loops[brep.half_edges[he].loop_].face;
    let cyl_face_handle = match brep.faces.get(he_face).and_then(|f| f.surface.as_ref()) {
        Some(sof) if matches!(sof.geometry, SurfaceGeometry::Cylinder { .. }) => he_face,
        _ => brep.loops[brep.half_edges[twin].loop_].face,
    };
    let (cyl_he, plane_he) = if he_face == cyl_face_handle {
        (he, twin)
    } else {
        (twin, he)
    };

    // Orient the sliced curve to match the edge's primary half-edge
    // direction. `sliced` currently runs p_start → p_end. The edge's
    // primary half-edge runs from edge.half_edge.origin to its
    // twin.origin. If origin matches p_start's vertex (close to
    // p_start in 3D), keep direction; otherwise reverse.
    let primary_he = brep.edges[edge].half_edge;
    let primary_origin_pos = brep.vertices[brep.half_edges[primary_he].origin].position;
    let oriented_sliced =
        if primary_origin_pos.distance_to(&p_start) <= primary_origin_pos.distance_to(&p_end) {
            sliced
        } else {
            reverse_nurbs(&sliced)
        };
    if let Some(edge_mut) = brep.edges.get_mut(edge) {
        let geometry = CurveGeometry::Nurbs(oriented_sliced);
        let geometry = geometry
            .canonical_line_in_range((0.0, 1.0), brep.tolerance.linear())
            .unwrap_or(geometry);
        edge_mut.curve = Some(CurveOnEdge {
            geometry,
            t_start: 0.0,
            t_end: 1.0,
        });
    } else {
        return false;
    }

    // Pcurve direction: the edge's "primary" half-edge runs in the
    // 3D curve direction. Our arc runs p_start → p_end. The twin
    // runs p_end → p_start; its pcurve must be reversed.
    let edge_he = brep.edges[edge].half_edge;
    let cyl_pcurve = if cyl_he == edge_he {
        cyl_uv.clone()
    } else {
        reverse_nurbs(&cyl_uv)
    };
    let plane_pcurve = if plane_he == edge_he {
        plane_uv.clone()
    } else {
        reverse_nurbs(&plane_uv)
    };
    if let Some(he_mut) = brep.half_edges.get_mut(cyl_he) {
        he_mut.pcurve = Some(cyl_pcurve.into());
    }
    if let Some(twin_mut) = brep.half_edges.get_mut(plane_he) {
        twin_mut.pcurve = Some(plane_pcurve.into());
    }
    true
}
