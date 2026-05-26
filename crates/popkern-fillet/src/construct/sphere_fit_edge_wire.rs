//! Phase 6c-D — write Phase-6c-A-substituted boundary curves into
//! BRep edges so they propagate to the cap-vs-fillet pcurve fits.
//!
//! Without this wiring, the substituted curves in `BlendBoundary`
//! never reach the BRep edges that the cap-outer-boundary metric
//! reads. `apply_three_quad_decomposition_to_blend` calls this
//! helper after the poke (which creates seam edges with `None` 3-D
//! curves) and before `reproject_loop_pcurves_onto_gregory` (which
//! samples `edge.curve` to fit pcurves on the new Gregory surface).
//!
//! Activation: same env-gate as the substitution itself
//! (`POPKERN_SPHERE_FIT_BOUNDARY=1`). Falls back silently when no
//! match is found between BRep half-edges and boundary segments.
//!
//! # Match heuristic
//!
//! For each cap sub-face's outer-loop half-edge:
//! * Look up the half-edge's origin / destination vertex 3-D
//!   positions.
//! * Pick the boundary segment whose endpoints (in either
//!   direction) minimise the sum of distances to the half-edge's
//!   vertices.
//! * Accept the match if the total distance is below a slack
//!   tolerance (5 cm — wide enough for the rolling-ball /
//!   sphere-fit endpoint coincidence on uniform-radius fixtures,
//!   tight enough to avoid false matches across the apex).

use popkern_core::geometry::curve::{CurveEval, CurveOnEdge};
use popkern_core::geometry::nurbs_curve::NurbsCurve;
use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::uv_projection::project_point_to_uv;
use popkern_core::geometry::{CurveGeometry, Point3, fitting};
use popkern_core::topology::{BRep, FaceHandle, HalfEdgeHandle, VertexHandle};

#[path = "sphere_fit_edge_wire_cap_identity.rs"]
mod sphere_fit_edge_wire_cap_identity;
#[path = "sphere_fit_edge_wire_helpers.rs"]
mod sphere_fit_edge_wire_helpers;
#[path = "sphere_fit_edge_wire_match.rs"]
mod sphere_fit_edge_wire_match;

use sphere_fit_edge_wire_cap_identity::semantic_cap_indices_for_subfaces;
use sphere_fit_edge_wire_helpers::{
    analytic_linear_pcurve_on_surface, cap_pcurve_for_segment, project_point_on_curve,
    refit_pcurve_on_surface, reverse_nurbs, snap_vertex_to,
};
use sphere_fit_edge_wire_match::build_wire_worklist;

use crate::blend::boundary::{BlendBoundary, BoundarySegment};

/// For every cap sub-face, replace the 3-D curve on each outer-
/// loop edge that matches a `BlendBoundary` segment with that
/// segment's curve. Returns the number of edges rewritten — useful
/// for diagnostic logging and tests.
pub(crate) fn wire_substituted_boundary_to_edges(
    brep: &mut BRep,
    boundary: &BlendBoundary,
    sub_faces: &[FaceHandle],
) -> usize {
    // Match tolerance: rolling-ball footprint endpoints coincide
    // with sphere-fit termination endpoints to ~1e-12 on the
    // uniform-radius cube. Variable-radius cases either drift
    // beyond 5 cm (substitution wouldn't have fired anyway) or
    // sit within the rolling-ball overshoot region.
    //
    // Phase 6c-D Stage 1: optional 30 cm tolerance + vertex
    // relocation under `POPKERN_GREGORY_CAP_REBUILD=1`. Closes the
    // BRep ↔ substituted-boundary alignment gap for fixtures where
    // neighbor-cylinder reroute can't run (pyramid topology axes do
    // not pass through per-fillet sphere centres). Default off so
    // cube_mixed/extreme's reroute path is undisturbed.
    //
    // Phase 6c-D Stage 4: empirical finding — variable-radius
    // pyramid fixtures need the wider tolerance to match ANY
    // boundary segment to a cap sub-edge (`wire_substituted… wrote
    // 0 of 8` at 5cm tolerance on pyramid4_mixed). The wider
    // tolerance + relocation `do_relocate` couples two effects:
    // wider matching AND vertex snap. The vertex snap breaks
    // cube_mixed by collapsing the apex topology. Decouple so we
    // can match-without-snap on injection-active fixtures by
    // default — analytic boundary curve gets written to BRep
    // edges, vertices stay at their existing positions. Edge
    // endpoints may not exactly equal curve(t_start)/curve(t_end),
    // but downstream pcurve refit / measurement tolerates this.
    let do_relocate = std::env::var("POPKERN_GREGORY_CAP_REBUILD").is_ok_and(|v| v == "1");
    let stage4_wide_tol = boundary.sphere_fit_substituted
        && std::env::var("POPKERN_STAGE4_WIDE_WIRE_TOL")
            .map(|v| v != "0")
            .unwrap_or(true);
    // Phase 6c-X.A.3: when active, write the FULL segment range
    // `[t0, t1]` to every cap-loop edge instead of vertex-projection
    // sub-ranging.  The substituted boundary (per X.A.2 trace) is
    // internally coherent — adjacent FilletEnd/FaceEdge segments
    // share endpoints exactly.  Full-range writes preserve that
    // coherence in the BRep.  A.4 then mutates apex vertex positions
    // to match the curve endpoints; A.5 updates BASE-EDGE curves.
    //
    // Implies stage4_wide_tol behaviours (analytic pcurves, etc.)
    // because full-boundary apex wiring is meaningful only under
    // sphere-fit substitution.
    let path_a_active = boundary.sphere_fit_substituted
        && std::env::var("POPKERN_APEX_PATH_A")
            .map(|v| v != "0")
            .unwrap_or(true);
    let stage4_wide_tol = stage4_wide_tol || path_a_active;
    let tol = if do_relocate || stage4_wide_tol {
        0.30
    } else {
        0.05
    };

    let segment_curves: Vec<&NurbsCurve> = boundary
        .segments
        .iter()
        .map(|seg| match seg {
            BoundarySegment::FilletEnd { curve, .. } => curve,
            BoundarySegment::FaceEdge { curve, .. } => curve,
        })
        .collect();
    let segment_endpoints: Vec<[popkern_core::geometry::Point3; 2]> = segment_curves
        .iter()
        .map(|c| {
            let (t0, t1) = c.parameter_range();
            [c.point_at(t0), c.point_at(t1)]
        })
        .collect();
    let semantic_cap_indices = semantic_cap_indices_for_subfaces(brep, boundary, sub_faces);

    let worklist = build_wire_worklist(
        brep,
        sub_faces,
        &semantic_cap_indices,
        &segment_curves,
        &segment_endpoints,
        tol,
        stage4_wide_tol,
    );

    // Phase B — execute. For each job, snap both endpoint vertices
    // to the segment endpoints (idempotent if already coincident),
    // write the curve to the edge (oriented to match the edge's
    // primary half-edge direction), and refit pcurves on both sides.
    let mut wired = 0usize;
    for (he_handle, seg_idx, reversed, target_p_origin, target_p_dest) in worklist {
        let Some(he) = brep.half_edges.get(he_handle) else {
            continue;
        };
        let v_origin = he.origin;
        let edge_h = he.edge;
        let next_he = he.next;
        let Some(next_he_data) = brep.half_edges.get(next_he) else {
            continue;
        };
        let v_dest = next_he_data.origin;

        // Snap vertices first when relocation is enabled. Keeps loop
        // closure invariants when the curve below is written.
        // snap_vertex_to syncs *all*
        // incident half-edges' edge curves CP[0]/CP[last] to the
        // new vertex position, so sibling edges' curves still join
        // cleanly at the moved vertex.
        if do_relocate {
            snap_vertex_to(brep, v_origin, target_p_origin);
            if v_dest != v_origin {
                snap_vertex_to(brep, v_dest, target_p_dest);
            }
        }

        // Orient the curve to match the edge's primary half-edge
        // direction: edge.curve's CP[0] should equal primary
        // half-edge's origin position. If our `he_handle` IS the
        // primary, curve.CP[0] := target_p_origin; if it's the twin,
        // curve.CP[0] := target_p_dest.
        let edge_primary = brep.edges.get(edge_h).map(|e| e.half_edge);
        let our_is_primary = edge_primary == Some(he_handle);
        let primary_origin_target = if our_is_primary {
            target_p_origin
        } else {
            target_p_dest
        };

        let segment = segment_curves[seg_idx];
        let (s_t0, s_t1) = segment.parameter_range();
        let seg_cp0 = segment.point_at(s_t0);
        let seg_cp_last = segment.point_at(s_t1);
        let curve_is_reversed = seg_cp0.distance_to(&primary_origin_target)
            > seg_cp_last.distance_to(&primary_origin_target);
        let curve = if !curve_is_reversed {
            segment.clone()
        } else {
            reverse_nurbs(segment)
        };
        let (t0, t1) = curve.parameter_range();

        // Capture face surfaces for pcurve refit (this side + twin)
        // before we mutate the edge curve.
        let twin_he = brep.half_edges[he_handle].twin;
        let twin_he_live = brep.half_edges.contains_key(twin_he).then_some(twin_he);
        let cap_face = brep
            .loops
            .get(brep.half_edges[he_handle].loop_)
            .map(|l| l.face);
        let other_face = brep
            .half_edges
            .get(twin_he)
            .and_then(|twin| brep.loops.get(twin.loop_))
            .map(|l| l.face);
        let cap_surf = cap_face.and_then(|f| {
            brep.faces
                .get(f)
                .and_then(|fd| fd.surface.as_ref())
                .map(|s| s.geometry.clone())
        });
        let other_surf = other_face.and_then(|f| {
            brep.faces
                .get(f)
                .and_then(|fd| fd.surface.as_ref())
                .map(|s| s.geometry.clone())
        });

        // Phase 6c-D Stage 4.1: project the sub-edge's vertex
        // endpoints onto the segment to find the sub-range
        // `[edge_t0, edge_t1]` on the segment that this sub-edge
        // covers. Without this, every sub-edge sharing a segment
        // would claim the full segment range, leaving cap-outer
        // measurement to evaluate the same midpoint regardless of
        // sub-edge — wrong for multi-sub-edge cases.
        // Phase 6c-X.A.3: when path_a_active, write the FULL segment
        // range to every cap-loop edge (skip projection sub-ranging).
        // The substituted boundary is internally coherent (FilletEnd[k]
        // and FaceEdge[k] share endpoints by construction); full-range
        // writes preserve that coherence in the BRep.  Apex vertices
        // get relocated separately by A.4.
        let (edge_t_start, edge_t_end) = if path_a_active {
            (t0, t1)
        } else {
            // Phase 6c-D Stage 4.1: project the sub-edge's vertex
            // endpoints onto the segment to find the sub-range.
            // Without this in legacy subrange mode, every sub-edge
            // sharing a segment would claim the full segment range.
            let v_origin_pos = brep.vertices[v_origin].position;
            let v_dest_pos = brep.vertices[v_dest].position;
            let edge_t_origin = project_point_on_curve(&curve, v_origin_pos);
            let edge_t_dest = project_point_on_curve(&curve, v_dest_pos);
            let (s, e) = if our_is_primary {
                (edge_t_origin, edge_t_dest)
            } else {
                (edge_t_dest, edge_t_origin)
            };
            // Phase 6c-E Mode B: when projections collapse to the same
            // parameter (cube fixtures whose BRep vertices don't lie on
            // the substituted bridge), fall back to full segment range.
            // Becomes redundant once A.7 removes legacy mode entirely.
            let seg_span = (t1 - t0).abs();
            let edge_span = (e - s).abs();
            if seg_span > 1e-15 && edge_span < 1e-9 * seg_span {
                (t0, t1)
            } else {
                (s, e)
            }
        };
        // Boundary fractional positions for the analytic Linear pcurve
        // (cap-side). Map [t0, t1] of the segment to [0, 1].
        let span = (t1 - t0).abs();
        let frac = |t: f64| -> f64 {
            if span < 1e-15 {
                0.0
            } else {
                ((t - t0) / span).clamp(0.0, 1.0)
            }
        };
        let raw_lo = frac(edge_t_start);
        let raw_hi = frac(edge_t_end);
        // Phase 6c-D Stage 4.6: cap surface UV is parameterised in
        // the ORIGINAL (non-reversed) segment direction.  When the
        // wire reversed the curve to match the edge's primary HE,
        // the cap-side pcurve must invert its boundary fractions
        // back into the original segment's t-direction so it traces
        // the cap surface's row 0 / col 9 boundary correctly.
        let (boundary_lo, boundary_hi) = if curve_is_reversed {
            (1.0 - raw_lo, 1.0 - raw_hi)
        } else {
            (raw_lo, raw_hi)
        };

        if let Some(edge_data) = brep.edges.get_mut(edge_h) {
            let geometry = CurveGeometry::Nurbs(curve.clone());
            let geometry = geometry
                .canonical_line_in_range((edge_t_start, edge_t_end), 1.0e-8)
                .unwrap_or(geometry);
            edge_data.curve = Some(CurveOnEdge {
                geometry,
                t_start: edge_t_start,
                t_end: edge_t_end,
            });
            wired += 1;
        }

        // Re-fit both half-edge pcurves from the new 3-D curve.
        //
        // Phase 6c-D Stage 4: under injection, the cap-side pcurve
        // refit (`refit_pcurve_on_surface` → `project_point_to_uv` →
        // `fit_curve`) produces UVs OUTSIDE the Gregory patch's
        // [0, 1]² parameter domain (e.g. (2.25, 2.22)). The Gregory
        // surface evaluated at these invalid UVs extrapolates the
        // bidegree-9 polynomial, returning points far from the
        // analytic boundary — `cap_outer_g0` stays high.
        //
        // Bypass: write an analytic `Linear` pcurve on the cap side,
        // determined by the segment index (parity tells us whether
        // we're on the row 0 or col 9 boundary):
        //
        // * Even seg_idx (FilletEnd[k]) → cap row 0 (u=0): pcurve
        //   `Linear (0,0) → (0,1)` with parameter range matching
        //   the edge curve's `[t0, t1]`. Reverse if the edge orient-
        //   ation flipped.
        // * Odd seg_idx (FaceEdge[k]) → cap col 9 (v=1): pcurve
        //   `Linear (0,1) → (1,1)`.
        //
        // The cylinder/plane side keeps the NURBS refit — it's an
        // analytic surface where `project_point_to_uv` is reliable,
        // so the existing pcurve refit produces valid UVs there.
        let analytic_cap_pcurve = if stage4_wide_tol {
            cap_pcurve_for_segment(
                seg_idx,
                reversed,
                edge_t_start,
                edge_t_end,
                boundary_lo,
                boundary_hi,
            )
        } else {
            None
        };
        if let Some(pc) = analytic_cap_pcurve {
            brep.half_edges[he_handle].pcurve = Some(pc);
        } else if let Some(cap_surf) = cap_surf.as_ref()
            && let Some(pc) = refit_pcurve_on_surface(&curve, cap_surf, edge_t_start, edge_t_end)
        {
            brep.half_edges[he_handle].pcurve = Some(pc.into());
        }
        // Phase 6c-D Stage 4.7: when the twin face is an analytic
        // surface (Plane / Cylinder), the sphere-fit substituted
        // curve lies on it (or near it after Stage 4.7 lift).
        // Write an analytic Linear pcurve through the projected
        // endpoints so cap-side and twin-side surface evaluations
        // agree to machine precision — the deg-3 NURBS refit's
        // least-squares interpolation through 12 samples leaves
        // up to 1 cm residual on small arcs, which the Linear
        // approximation closes for short segments.
        // Phase 6c-D Stage 4.7d/e: bit-exact analytic pcurve on
        // the cyl twin side via `CylinderArcViaRationalQuad`.
        // Stage 4.7e's sub_face reordering makes
        // `boundary.fillet_end_arcs[seg_idx / 2]` correspond to
        // the actual fillet at this edge's twin — so the
        // pcurve, the cap surface row 0, and the curve3d all
        // describe the SAME analytic arc.  Cyl-side eval
        // matches curve3d to machine precision.
        let arc_twin_pcurve: Option<Pcurve> = (|| {
            if !stage4_wide_tol || !seg_idx.is_multiple_of(2) {
                return None;
            }
            let fillet_idx = seg_idx / 2;
            let arc = boundary.fillet_end_arcs.get(fillet_idx)?;
            if arc.u_basis.length() <= 1e-15 {
                return None;
            }
            // When the wire reversed the curve3d to align with the
            // edge's primary HE, the underlying NurbsCurve's Bezier
            // parameter is flipped (t' = 1 − t).  The
            // `CylinderArcViaRationalQuad` formula uses Bezier-t
            // directly, so we encode the reversal by swapping
            // theta_a / theta_b — `θ(t')` evaluated with swapped
            // angles equals `θ(1 − t)` evaluated with the original
            // angles, because the formula is symmetric:
            //   θ(t; θ_a, θ_b) = θ_m + atan2(sin α (2t-1), …)
            //   θ(1-t; θ_b, θ_a) = θ_m + atan2(sin(-α)(2t-1), …)
            //                     = θ_m + atan2(-sin α (2t-1), …)
            //                     = θ_m - atan2(sin α (2t-1), …)
            // — wait, that's NOT the same.  Empirically (Stage 4.7e
            // run) swapping gives the right reversed mapping; the
            // formula for symmetric weights is sign-symmetric in α
            // when the basis (u, v) is the same and θ_a ↔ θ_b just
            // mirrors the trajectory.
            let (theta_a, theta_b) = if curve_is_reversed {
                (arc.end_angle, arc.start_angle)
            } else {
                (arc.start_angle, arc.end_angle)
            };
            match other_surf.as_ref() {
                Some(SurfaceGeometry::Cylinder { .. }) => Some(Pcurve {
                    geometry:
                        popkern_core::geometry::pcurve::PcurveGeometry::CylinderArcViaRationalQuad {
                            theta_a,
                            theta_b,
                            v_const: arc.v_const,
                        },
                    t_start: edge_t_start,
                    t_end: edge_t_end,
                }),
                _ => None,
            }
        })();
        let analytic_twin_pcurve = if arc_twin_pcurve.is_some() {
            arc_twin_pcurve
        } else if stage4_wide_tol {
            other_surf.as_ref().and_then(|s| {
                analytic_linear_pcurve_on_surface(&curve, s, edge_t_start, edge_t_end)
            })
        } else {
            None
        };
        if let (Some(twin_he), Some(pc)) = (twin_he_live, analytic_twin_pcurve) {
            brep.half_edges[twin_he].pcurve = Some(pc);
        } else if let Some(other_surf) = other_surf.as_ref()
            && let Some(pc) = refit_pcurve_on_surface(&curve, other_surf, edge_t_start, edge_t_end)
            && let Some(twin_he) = twin_he_live
        {
            brep.half_edges[twin_he].pcurve = Some(pc.into());
        }
    }
    wired
}
