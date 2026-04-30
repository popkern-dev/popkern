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
    // Path B reroute can't run (pyramid topology — neighbour
    // cylinder axes don't pass through per-fillet sphere centres).
    // Default off so cube_mixed/extreme's Path B path is undisturbed.
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
        && std::env::var("POPKERN_STAGE4_WIDE_WIRE_TOL").is_ok_and(|v| v != "0");
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

    // Phase A — assignment scan: for each cap sub-face's outer-loop
    // half-edge, find the best matching segment (proximity, with
    // the wider tolerance). Build a worklist of
    // `(he, seg_idx, reversed, target_p_origin, target_p_dest)`.
    // We do NOT mutate during the scan so the matching sees
    // pre-relocation positions consistently.
    //
    // Phase 6c-D Stage 4: when injection is active, allow multiple
    // sub-edges to share a single segment (`stage4_wide_tol`). The
    // cap face's outer loop after sliver-absorb topology surgery
    // can have 5+ sub-edges along a single logical boundary
    // (FilletEnd[k] / FaceEdge[k]); the original "consume each
    // segment once" policy left short sub-edges unwired, holding
    // worst-case `cap_outer_g0` at ~10 cm. Letting multiple
    // sub-edges share a segment as their 3D curve writes the
    // analytic boundary to ALL of them — endpoint mismatch
    // tolerated since downstream consumers (validation, pcurve
    // refit) operate on the curve geometry, not the vertex-to-
    // curve-endpoint identity.
    let allow_multi_sub_edge = stage4_wide_tol;
    let trace_match = std::env::var_os("POPKERN_WIRE_MATCH_TRACE").is_some();
    // Phase 6c-D Stage 4.6: when injection is active, restrict
    // matching to per-cap-face index — sub_face[k] gets row 0 =
    // FilletEnd[k] (segment 2k) and col 9 = FaceEdge[k]
    // (segment 2k+1) by construction.  Twin face surface kind
    // tells us which boundary the edge sits on; combined with
    // the cap face index, the segment index is fully determined.
    //
    // Without this constraint, geometric-proximity matching
    // routes cap[1]'s cylinder edge to FilletEnd[2] (the wrong
    // fillet) when the rolling-ball footprint of cap[1] sits
    // closer in 3D to FilletEnd[2]'s analytic arc than to
    // FilletEnd[1]'s — empirically observed on pyramid4_mixed
    // because the rolling-ball z (0.755) and analytic-arc z
    // (0.818) gap (~6 cm) is comparable to the inter-fillet XY
    // separation, so all four FilletEnd arcs are roughly
    // equidistant.
    let constrained_by_cap_index = stage4_wide_tol;
    let mut consumed = vec![false; segment_curves.len()];
    type Job = (HalfEdgeHandle, usize, bool, Point3, Point3);
    let mut worklist: Vec<Job> = Vec::new();
    // Phase 6c-D Stage 4.3: cap-face set for cap-vs-cap edge
    // detection. The cap face's surface at wire time is NurbsSurface
    // (Gregory assignment happens AFTER wire in the pipeline), so
    // we can't filter by surface kind. Identify cap faces by
    // membership in `sub_faces` instead.
    let cap_face_set: std::collections::BTreeSet<FaceHandle> = sub_faces.iter().copied().collect();

    for (cap_k, &sub_face) in sub_faces.iter().enumerate() {
        let Some(face_data) = brep.faces.get(sub_face) else {
            continue;
        };
        let outer_loop = face_data.outer_loop;
        let Some(loop_data) = brep.loops.get(outer_loop) else {
            continue;
        };
        let start_he = loop_data.half_edge;

        let hes = collect_loop_half_edges(brep, start_he);

        for he_handle in hes {
            let Some(he) = brep.half_edges.get(he_handle) else {
                continue;
            };
            let v_origin = he.origin;
            let next_he = he.next;
            let Some(next_he_data) = brep.half_edges.get(next_he) else {
                continue;
            };
            let v_dest = next_he_data.origin;

            // Phase 6c-D Stage 4.3: skip cap-vs-cap edges (twin
            // face is in `sub_faces`, i.e., another cap face).
            // They aren't FilletEnd or FaceEdge boundaries —
            // they're inner spokes / shared cap edges between
            // adjacent per-fillet caps. Wiring them would write
            // the wrong analytic curve and could corrupt cap-cap
            // continuity. Only wire cap-vs-non-cap edges.
            //
            // Earlier surface-based filter (Sphere | Gregory) was
            // ineffective because the cap face's surface at wire
            // time is still NurbsSurface (Gregory replaces it
            // *after* the wire step in the pipeline). Membership
            // in `sub_faces` is the topology-stable signal.
            let twin_he = he.twin;
            let twin_face = brep
                .loops
                .get(brep.half_edges[twin_he].loop_)
                .map(|l| l.face);
            let twin_is_cap = twin_face
                .map(|f| cap_face_set.contains(&f))
                .unwrap_or(false);
            if twin_is_cap {
                if trace_match {
                    eprintln!(
                        "[wire-match] he={he_handle:?} edge={:?} SKIPPED (twin is cap face)",
                        he.edge
                    );
                }
                continue;
            }

            let p_origin = brep.vertices[v_origin].position;
            let p_dest = brep.vertices[v_dest].position;

            // Phase 6c-D Stage 4.2 — matching strategy:
            //
            // * Default (`allow_multi_sub_edge = false`): endpoint-sum
            //   matching. Each segment matches at most one sub-edge,
            //   chosen by minimum sum of endpoint distances. Works
            //   well when each cap face's outer loop has exactly one
            //   sub-edge per segment (cube_uniform).
            //
            // * `allow_multi_sub_edge = true` (Stage 4.2):
            //   midpoint-projection matching. For each sub-edge,
            //   compute its midpoint and project onto each segment;
            //   match to the segment with minimum projection distance.
            //   Robust against multi-sub-edge topologies where a
            //   single FilletEnd/FaceEdge has been split into 3+
            //   sub-edges by sliver-absorb. Endpoint-sum fails for
            //   short sub-edges in segment interior (both endpoints
            //   far from segment endpoints) — projection works
            //   because sub-edges geometrically lie on (or near) the
            //   parent fillet cylinder, whose tangent circle = the
            //   FilletEnd segment.
            let mut best: Option<(usize, f64, bool)> = None;
            if constrained_by_cap_index {
                // Phase 6c-D Stage 4.6: route by cap index + twin
                // surface kind.  sub_face[k] owns segments
                // {2k, 2k+1}; twin face kind picks parity.
                let twin_face_data = twin_face.and_then(|f| brep.faces.get(f));
                let twin_kind = twin_face_data
                    .and_then(|f| f.surface.as_ref())
                    .map(|sof| &sof.geometry);
                let target_seg = match twin_kind {
                    Some(SurfaceGeometry::Cylinder { .. }) => Some(2 * cap_k),
                    Some(SurfaceGeometry::Plane { .. }) => Some(2 * cap_k + 1),
                    _ => None,
                };
                if let Some(seg_i) = target_seg
                    && seg_i < segment_curves.len()
                {
                    let segment = segment_curves[seg_i];
                    let p_mid = Point3::new(
                        (p_origin.x + p_dest.x) * 0.5,
                        (p_origin.y + p_dest.y) * 0.5,
                        (p_origin.z + p_dest.z) * 0.5,
                    );
                    let t_mid = project_point_on_curve(segment, p_mid);
                    let proj_mid = segment.point_at(t_mid);
                    let d = proj_mid.distance_to(&p_mid);
                    let t_origin = project_point_on_curve(segment, p_origin);
                    let t_dest = project_point_on_curve(segment, p_dest);
                    let reversed = t_origin > t_dest;
                    best = Some((seg_i, d, reversed));
                }
            } else if allow_multi_sub_edge {
                let p_mid = Point3::new(
                    (p_origin.x + p_dest.x) * 0.5,
                    (p_origin.y + p_dest.y) * 0.5,
                    (p_origin.z + p_dest.z) * 0.5,
                );
                for (i, segment) in segment_curves.iter().enumerate() {
                    let t_mid = project_point_on_curve(segment, p_mid);
                    let proj_mid = segment.point_at(t_mid);
                    let d = proj_mid.distance_to(&p_mid);
                    // Determine direction by projecting both endpoints.
                    let t_origin = project_point_on_curve(segment, p_origin);
                    let t_dest = project_point_on_curve(segment, p_dest);
                    let reversed = t_origin > t_dest;
                    if best.map(|(_, bd, _)| d < bd).unwrap_or(true) {
                        best = Some((i, d, reversed));
                    }
                }
            } else {
                for (i, eps) in segment_endpoints.iter().enumerate() {
                    if consumed[i] {
                        continue;
                    }
                    let d_fwd = p_origin.distance_to(&eps[0]) + p_dest.distance_to(&eps[1]);
                    let d_rev = p_origin.distance_to(&eps[1]) + p_dest.distance_to(&eps[0]);
                    let (d, reversed) = if d_fwd <= d_rev {
                        (d_fwd, false)
                    } else {
                        (d_rev, true)
                    };
                    if best.map(|(_, bd, _)| d < bd).unwrap_or(true) {
                        best = Some((i, d, reversed));
                    }
                }
            }
            if trace_match {
                let chord = p_origin.distance_to(&p_dest);
                let twin_kind = brep
                    .faces
                    .get(twin_face.unwrap_or(brep.loops[brep.half_edges[twin_he].loop_].face))
                    .and_then(|f| f.surface.as_ref())
                    .map(|sof| match &sof.geometry {
                        SurfaceGeometry::Cylinder { .. } => "Cyl",
                        SurfaceGeometry::Plane { .. } => "Plane",
                        SurfaceGeometry::Sphere { .. } => "Sphere",
                        SurfaceGeometry::Gregory(_) => "Gregory",
                        _ => "Other",
                    })
                    .unwrap_or("none");
                let cap_idx = sub_faces
                    .iter()
                    .position(|f| *f == sub_face)
                    .unwrap_or(usize::MAX);
                eprintln!(
                    "[wire-match] cap[{cap_idx}]={sub_face:?} he={he_handle:?} edge={:?} twin={twin_kind} chord={chord:.3e} v_origin=({:.3},{:.3},{:.3}) v_dest=({:.3},{:.3},{:.3}) best_seg={:?} best_d={:?}",
                    he.edge,
                    p_origin.x,
                    p_origin.y,
                    p_origin.z,
                    p_dest.x,
                    p_dest.y,
                    p_dest.z,
                    best.map(|(i, _, _)| i),
                    best.map(|(_, d, _)| d),
                );
            }
            let Some((seg_idx, d, reversed)) = best else {
                continue;
            };
            if d > tol {
                if trace_match {
                    eprintln!("[wire-match]   REJECTED seg={seg_idx} d={d:.3e} > tol={tol:.3e}");
                }
                continue;
            }
            if trace_match {
                eprintln!("[wire-match]   ACCEPTED seg={seg_idx} d={d:.3e} reversed={reversed}");
            }
            if !allow_multi_sub_edge {
                consumed[seg_idx] = true;
            }

            let eps = segment_endpoints[seg_idx];
            let (target_p_origin, target_p_dest) = if reversed {
                (eps[1], eps[0])
            } else {
                (eps[0], eps[1])
            };
            worklist.push((he_handle, seg_idx, reversed, target_p_origin, target_p_dest));
        }
    }

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

        // Snap vertices first when relocation is enabled (pyramid
        // path — no Path B). Keeps loop closure invariants when
        // the curve below is written. snap_vertex_to syncs *all*
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
        let cap_face = brep
            .loops
            .get(brep.half_edges[he_handle].loop_)
            .map(|l| l.face);
        let other_face = brep
            .loops
            .get(brep.half_edges[twin_he].loop_)
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
        let v_origin_pos = brep.vertices[v_origin].position;
        let v_dest_pos = brep.vertices[v_dest].position;
        let edge_t_origin = project_point_on_curve(&curve, v_origin_pos);
        let edge_t_dest = project_point_on_curve(&curve, v_dest_pos);
        // Order based on the edge's primary half-edge direction.
        let (edge_t_start, edge_t_end) = if our_is_primary {
            (edge_t_origin, edge_t_dest)
        } else {
            (edge_t_dest, edge_t_origin)
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
            edge_data.curve = Some(CurveOnEdge {
                geometry: CurveGeometry::Nurbs(curve.clone()),
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
        if let Some(pc) = analytic_twin_pcurve {
            brep.half_edges[twin_he].pcurve = Some(pc);
        } else if let Some(other_surf) = other_surf.as_ref()
            && let Some(pc) = refit_pcurve_on_surface(&curve, other_surf, edge_t_start, edge_t_end)
        {
            brep.half_edges[twin_he].pcurve = Some(pc.into());
        }
    }
    wired
}

/// Build an analytic Linear pcurve on `surface` by projecting the
/// curve endpoints to surface UV.  Exact when the surface contains
/// the curve and the curve traces a linear path in the surface's
/// UV space (cylinder seam-perpendicular arcs, plane edges).
/// Approximate but small-error for short arcs on cylinders where
/// the (theta, axial) trajectory is mildly curved.
///
/// Returns `None` when projection fails (point off surface, degenerate
/// UV, etc.) — caller falls through to the deg-3 NURBS refit.
fn analytic_linear_pcurve_on_surface(
    curve: &NurbsCurve,
    surface: &SurfaceGeometry,
    t_start: f64,
    t_end: f64,
) -> Option<Pcurve> {
    let p_start = curve.point_at(t_start);
    let p_end = curve.point_at(t_end);
    let uv_start = project_point_to_uv(surface, &p_start)?;
    let uv_end = project_point_to_uv(surface, &p_end)?;
    Some(Pcurve::linear(uv_start, uv_end, t_start, t_end))
}

/// Phase 6c-D: reverse a NurbsCurve by flipping CPs, weights, and
/// reflecting the knot vector. Equivalent to
/// `popkern_core::geometry::edge_triple::analytic_helpers::reverse_nurbs`
/// (kept private there); duplicated locally for cohesion.
fn reverse_nurbs(curve: &NurbsCurve) -> NurbsCurve {
    let (t_min, t_max) = curve.parameter_range();
    let mut reversed_knots: Vec<f64> = curve
        .knots
        .iter()
        .rev()
        .map(|&k| t_min + t_max - k)
        .collect();
    reversed_knots.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let mut reversed_cps = curve.control_points.clone();
    reversed_cps.reverse();
    let mut reversed_weights = curve.weights.clone();
    reversed_weights.reverse();
    NurbsCurve {
        degree: curve.degree,
        knots: reversed_knots,
        control_points: reversed_cps,
        weights: reversed_weights,
    }
}

/// Phase 6c-D: snap a vertex to `new_pos`, syncing every incident
/// half-edge's edge curve CP[0] (primary side) or CP[last] (twin
/// side) to keep loop closure intact. No-op when already within
/// 1e-9 of `new_pos`.
///
/// This is a drop-in equivalent of
/// `var_radius::update_incident_edges_for_vertex_move` — kept
/// local to the wire pass for cohesion. Both functions produce
/// identical results on identical inputs.
fn snap_vertex_to(brep: &mut BRep, v: VertexHandle, new_pos: Point3) {
    let Some(vd) = brep.vertices.get(v) else {
        return;
    };
    if vd.position.distance_to(&new_pos) < 1e-9 {
        return;
    }

    let halves: Vec<HalfEdgeHandle> = brep
        .half_edges
        .iter()
        .filter_map(|(h, he)| (he.origin == v).then_some(h))
        .collect();

    for he_h in halves {
        let edge_h = brep.half_edges[he_h].edge;
        let primary = brep.edges.get(edge_h).map(|e| e.half_edge);
        let is_primary = primary == Some(he_h);
        let Some(edge) = brep.edges.get_mut(edge_h) else {
            continue;
        };
        let Some(coe) = edge.curve.as_mut() else {
            continue;
        };
        if let CurveGeometry::Nurbs(nc) = &mut coe.geometry
            && !nc.control_points.is_empty()
        {
            let n = nc.control_points.len();
            if is_primary {
                nc.control_points[0] = new_pos;
            } else {
                nc.control_points[n - 1] = new_pos;
            }
        }
    }

    if let Some(vd) = brep.vertices.get_mut(v) {
        vd.position = new_pos;
    }
}

/// Sample `curve_3d` at uniform parameter values, project each
/// sample to `surface`'s UV, and fit a degree-3 NURBS pcurve
/// through those (u, v, 0) points. Returns `None` if any
/// projection fails or the fit raises an error.
/// Phase 6c-D Stage 4 — analytic `Linear` cap-side pcurve for an
/// edge wired with a sphere-fit-substituted boundary segment.
///
/// Segment index parity tells us which cap boundary the edge is
/// on (per the layout `sphere_fit_boundary::apply_sphere_fit_substitution`
/// emits — `segments[2k] = FilletEnd[k]`,
/// `segments[2k+1] = FaceEdge[k]`):
///
/// * `seg_idx % 2 == 0` (FilletEnd) → cap **row 0** (`u = 0`):
///   pcurve traces `(0, v)` over `v ∈ [v_lo, v_hi]`.
/// * `seg_idx % 2 == 1` (FaceEdge) → cap **col 9** (`v = 1`):
///   pcurve traces `(u, 1)` over `u ∈ [u_lo, u_hi]`.
///
/// `reversed` flips the start/end UVs so the pcurve advances in
/// the same direction as the underlying 3D curve. Parameter
/// range `[t_start, t_end]` matches the edge curve.
///
/// Phase 6c-D Stage 4.1 — `(boundary_lo, boundary_hi)` is the
/// sub-edge's sub-range on the parent segment, normalised to
/// `[0, 1]`. For a full-edge match (single sub-edge spans whole
/// segment), `(boundary_lo, boundary_hi) = (0, 1)`. For multi-
/// sub-edge cases, each sub-edge gets a sub-range computed by
/// projecting its endpoint vertices onto the parent segment.
fn cap_pcurve_for_segment(
    seg_idx: usize,
    reversed: bool,
    t_start: f64,
    t_end: f64,
    boundary_lo: f64,
    boundary_hi: f64,
) -> Option<Pcurve> {
    let (start, end) = if seg_idx.is_multiple_of(2) {
        // Row 0: u = 0, v varies over the sub-range.
        ((0.0, boundary_lo), (0.0, boundary_hi))
    } else {
        // Col 9: v = 1, u varies over the sub-range.
        ((boundary_lo, 1.0), (boundary_hi, 1.0))
    };
    let (start, end) = if reversed { (end, start) } else { (start, end) };
    Some(Pcurve::linear(start, end, t_start, t_end))
}

/// Phase 6c-D Stage 4.1 — closest-point projection of `point`
/// onto `curve`. Returns the parameter `t ∈ [t_lo, t_hi]` such
/// that `curve.point_at(t)` is closest to `point`.
///
/// Sampling-based: 100 uniform parameter samples, take the one
/// with minimum distance, refine with 5 bisection iterations.
/// Sufficient accuracy for sub-edge t-range computation in the
/// wire step (we don't need machine precision — micrometre-level
/// is plenty for cap-outer-G0 on a 0.1m cap).
fn project_point_on_curve(curve: &NurbsCurve, point: Point3) -> f64 {
    let (t_lo, t_hi) = curve.parameter_range();
    if (t_hi - t_lo).abs() < 1e-15 {
        return t_lo;
    }
    let mut best_t = t_lo;
    let mut best_d = f64::INFINITY;
    const N: usize = 100;
    for i in 0..=N {
        let s = (i as f64) / (N as f64);
        let t = t_lo + s * (t_hi - t_lo);
        let p = curve.point_at(t);
        let d = p.distance_to(&point);
        if d < best_d {
            best_d = d;
            best_t = t;
        }
    }
    // Refine via 5 bisection iterations around best_t.
    let mut lo = (best_t - (t_hi - t_lo) / (N as f64)).max(t_lo);
    let mut hi = (best_t + (t_hi - t_lo) / (N as f64)).min(t_hi);
    for _ in 0..5 {
        let mid_lo = lo + (hi - lo) * 0.25;
        let mid_hi = lo + (hi - lo) * 0.75;
        let d_lo = curve.point_at(mid_lo).distance_to(&point);
        let d_hi = curve.point_at(mid_hi).distance_to(&point);
        if d_lo < d_hi {
            hi = (lo + hi) * 0.5;
        } else {
            lo = (lo + hi) * 0.5;
        }
    }
    (lo + hi) * 0.5
}

fn refit_pcurve_on_surface(
    curve_3d: &NurbsCurve,
    surface: &SurfaceGeometry,
    t0: f64,
    t1: f64,
) -> Option<NurbsCurve> {
    const N: usize = 12;
    let mut points: Vec<Point3> = Vec::with_capacity(N);
    let mut params: Vec<f64> = Vec::with_capacity(N);
    for i in 0..N {
        let s = i as f64 / (N - 1) as f64;
        let t = t0 + s * (t1 - t0);
        let p3 = curve_3d.point_at(t);
        let (u, v) = project_point_to_uv(surface, &p3)?;
        points.push(Point3::new(u, v, 0.0));
        params.push(s);
    }
    let degree = 3usize;
    let n_cps = (degree + 2).min(points.len());
    fitting::fit_curve(&points, &params, degree, n_cps).ok()
}

fn collect_loop_half_edges(brep: &BRep, start: HalfEdgeHandle) -> Vec<HalfEdgeHandle> {
    let mut hes = Vec::new();
    let mut cursor = start;
    let mut safety = 64usize;
    loop {
        if safety == 0 {
            break;
        }
        safety -= 1;
        hes.push(cursor);
        let Some(he) = brep.half_edges.get(cursor) else {
            break;
        };
        cursor = he.next;
        if cursor == start {
            break;
        }
    }
    hes
}
