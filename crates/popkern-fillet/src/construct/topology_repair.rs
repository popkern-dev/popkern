//! Topology repair for cylinder fillet faces with off-surface vertices.
//!
//! Apex-only filleting on a pyramid (or any partial-edge filleting
//! where the cylinder ribbon ends at an unfilleted original vertex)
//! leaves the cylinder face's outer loop containing two consecutive
//! chord edges through the original (off-cylinder) vertex. The
//! chord-pair lies on the base plane of the original face, but
//! geometrically passes through the cylinder *interior* — OCCT-strict
//! slicers reject the cylinder face as non-manifold and render the
//! fillet as a sharp edge.
//!
//! ## Surgery
//!
//! For each cylinder face whose outer loop contains a chord pair
//! `axial_a → off_cyl_vertex → axial_b`, we:
//!
//! 1. Compute the cylinder × adjacent-plane intersection (an
//!    ellipse / circle / line, via `intersect_surfaces`'s analytic
//!    plane-cylinder solver — never the marching path).
//! 2. Trim the intersection between the two cylinder-surface anchors
//!    (`p_a` and `p_b`).
//! 3. Insert the new arc edge as a face split on the cylinder face,
//!    producing a sliver triangle between the new arc and the two
//!    chord edges.
//! 4. Merge the sliver with the adjacent base-plane face via `kef`
//!    on a shared chord edge.
//! 5. The remaining chord edge becomes interior to the merged face —
//!    `kev` to remove it.
//!
//! After this, the cylinder face's outer loop ends at the trim
//! anchors via the analytic intersection arc; the off-cylinder
//! vertex remains in the BRep but is no longer in the cylinder's
//! loop. The base face's loop now extends along the arc.
//!
//! ## Scope
//!
//! Only handles cylinder faces with **chord pairs to an off-cylinder
//! vertex on a planar adjacent face**. The most common case in our
//! fillet output. Other patterns (cylinder × NURBS-cap, helical end
//! caps, multi-vertex chord chains) are out of scope.

use super::cyl_plane_intersect::analytic_cyl_plane_intersection;
use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::{CurveEval, CurveGeometry, CurveOnEdge};
use popkern_core::geometry::intersection::surface_surface::{
    SurfaceSurfaceIntersection, intersect_surfaces,
};
use popkern_core::geometry::nurbs_curve::NurbsCurve;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::uv_projection::tight_bounded_parameter_range;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle, VertexHandle};

#[path = "topology_repair_candidates.rs"]
mod topology_repair_candidates;
#[path = "topology_repair_diagnostics.rs"]
mod topology_repair_diagnostics;
#[path = "topology_repair_diff_planes.rs"]
mod topology_repair_diff_planes;
#[path = "topology_repair_geometry.rs"]
mod topology_repair_geometry;
#[path = "topology_repair_snap.rs"]
mod topology_repair_snap;

use topology_repair_candidates::{
    ChordPairCandidate, find_chord_pair_candidates, is_cylinder_face, vertex_on_cylinder,
};
use topology_repair_diagnostics::dump_cylinder_face;
use topology_repair_diff_planes::try_repair_chord_pair_different_planes;
use topology_repair_geometry::{
    anchor_off_cyl_distance, face_loop_contains_vertex, fit_cylinder_uv_pcurve,
    remap_to_unit_range, reproject_loop_pcurves_into_plane, reverse_nurbs, slice_curve,
};
use topology_repair_snap::snap_incident_edge_curve_endpoints;

/// Walk every cylinder face. For each off-cylinder vertex in its
/// outer loop, try to insert a cyl × plane intersection arc that
/// bypasses the vertex, and merge the resulting sliver with the
/// adjacent planar face. Returns the number of off-cyl vertices
/// successfully repaired.
///
/// Each repair attempt runs against a *snapshot* of the BRep —
/// if the attempt produces an invalid topology (`brep.validate()`
/// fails) the snapshot is restored and we move on. This keeps a
/// single bad attempt from poisoning the rest of the build.
pub(crate) fn repair_cylinder_face_off_surface_vertices(brep: &mut BRep) -> usize {
    let trace = std::env::var("POPKERN_TOPO_REPAIR_TRACE").is_ok();
    let dump = std::env::var("POPKERN_TOPO_REPAIR_DUMP").is_ok();
    // Per-face attempt cap.  Each successful repair removes one
    // off-cyl candidate; an upper bound keeps the loop bounded
    // even if surgery somehow generates new candidates.  Override
    // via `POPKERN_TOPO_REPAIR_MAX_ATTEMPTS` for diagnostics
    // (#136 iteration-compounding hypothesis).
    let max_attempts = std::env::var("POPKERN_TOPO_REPAIR_MAX_ATTEMPTS")
        .ok()
        .and_then(|v| v.parse::<usize>().ok())
        .unwrap_or(8);
    let face_handles: Vec<FaceHandle> = brep.faces.keys().collect();
    let mut repaired = 0;
    for fh in face_handles {
        if !is_cylinder_face(brep, fh) {
            continue;
        }
        if dump {
            dump_cylinder_face(brep, fh);
        }
        for _attempt in 0..max_attempts {
            if !brep.faces.contains_key(fh) {
                break;
            }
            let candidates = find_chord_pair_candidates(brep, fh);
            let n_cand = candidates.len();
            if candidates.is_empty() {
                break;
            }
            // Try every candidate in turn — some "off" vertices are
            // unrepairable (apex-side different-adjacent-face case)
            // and we don't want to bail on the whole face just
            // because the first one failed.
            let mut any_success = false;
            for cand in candidates {
                if trace {
                    eprintln!(
                        "[topo_repair] face {fh:?}: trying repair (n_cand={n_cand}, off_cyl_v={:?})",
                        cand.off_cyl_vertex
                    );
                }
                let snapshot = brep.clone();
                let surgery_ok = try_repair_chord_pair(brep, fh, cand);
                let validate_result = brep.validate();
                if trace {
                    eprintln!(
                        "[topo_repair]   surgery_ok={surgery_ok}, validate={:?}",
                        validate_result.as_ref().err()
                    );
                }
                if surgery_ok && validate_result.is_ok() {
                    repaired += 1;
                    if trace {
                        eprintln!("[topo_repair]   SUCCESS — repair count {repaired}");
                    }
                    any_success = true;
                    break;
                } else {
                    *brep = snapshot;
                }
            }
            if !any_success {
                break;
            }
        }
    }
    if trace {
        eprintln!("[topo_repair] total repairs: {repaired}");
    }
    repaired
}

fn try_repair_chord_pair(brep: &mut BRep, cyl_face: FaceHandle, cand: ChordPairCandidate) -> bool {
    let trace = std::env::var("POPKERN_TOPO_REPAIR_TRACE").is_ok();
    // The two chord edges (he_in and he_out) sit on planar adjacent
    // faces on their twin side. There are two cases:
    //  * Same plane on both sides (e.g. base corner with chord pair
    //    fully on the base) → insert arc + sliver triangle on cyl
    //    face. Handled by `try_repair_chord_pair_same_plane`.
    //  * Different planes (e.g. chord goes from cap-cyl boundary
    //    on slant_A to base-cyl corner) → snap the off-cyl vertex
    //    to the cyl ∩ plane_A ∩ plane_B point, then rewrite each
    //    chord's 3D curve to the analytic arc on its respective
    //    plane. Handled by `try_repair_chord_pair_different_planes`.
    //
    // The different-planes surgery is currently regressing N≥8
    // pyramids (max tess-nm divergence 0.12 → 1.06) — likely the
    // vertex-snap leaves shared neighbours (other cyls, slants)
    // inconsistent. Gate it behind `POPKERN_TOPO_REPAIR_DIFF=1` so
    // the same-plane case can land safely while we work the diff
    // case separately. The same-plane case is what fixed the base-
    // corner chord pairs on N=3..7 in the first place.
    let twin_in = brep.half_edges[cand.he_in].twin;
    let twin_out = brep.half_edges[cand.he_out].twin;
    let plane_face_in = brep.loops[brep.half_edges[twin_in].loop_].face;
    let plane_face_out = brep.loops[brep.half_edges[twin_out].loop_].face;
    if plane_face_in != plane_face_out {
        if std::env::var("POPKERN_TOPO_REPAIR_DIFF").is_err() {
            if trace {
                eprintln!(
                    "[topo_repair]   skip diff-planes case ({plane_face_in:?} vs {plane_face_out:?}) — gated"
                );
            }
            return false;
        }
        return try_repair_chord_pair_different_planes(
            brep,
            cyl_face,
            cand,
            plane_face_in,
            plane_face_out,
            trace,
        );
    }
    let plane_face = plane_face_in;
    let Some(plane_sof) = brep.faces.get(plane_face).and_then(|f| f.surface.as_ref()) else {
        if trace {
            eprintln!("[topo_repair]   reject: plane face has no surface");
        }
        return false;
    };
    if !matches!(plane_sof.geometry, SurfaceGeometry::Plane { .. }) {
        if trace {
            eprintln!(
                "[topo_repair]   reject: adjacent face surface is {:?}, not a plane",
                std::mem::discriminant(&plane_sof.geometry)
            );
        }
        return false;
    }

    let Some(cyl_sof) = brep.faces.get(cyl_face).and_then(|f| f.surface.as_ref()) else {
        if trace {
            eprintln!("[topo_repair]   reject: cylinder face has no surface");
        }
        return false;
    };
    let cyl_sof = cyl_sof.clone();
    let plane_sof_clone = plane_sof.clone();
    let p_anchor_in = brep.vertices[cand.anchor_in].position;
    let p_anchor_out = brep.vertices[cand.anchor_out].position;
    let chord_len = p_anchor_in.distance_to(&p_anchor_out);
    if chord_len < 1e-9 {
        if trace {
            eprintln!("[topo_repair]   reject: chord_len < 1e-9");
        }
        return false;
    }

    // Compute cylinder × plane intersection.
    let positions_1 = brep.face_outer_loop_positions(cyl_face);
    let positions_2 = brep.face_outer_loop_positions(plane_face);
    if positions_1.is_empty() || positions_2.is_empty() {
        if trace {
            eprintln!("[topo_repair]   reject: empty positions");
        }
        return false;
    }
    let (u1_range, v1_range) = tight_bounded_parameter_range(&positions_1, &cyl_sof.geometry);
    let (u2_range, v2_range) =
        tight_bounded_parameter_range(&positions_2, &plane_sof_clone.geometry);
    let tolerance = brep.tolerance;
    let result = match intersect_surfaces(
        &cyl_sof.geometry,
        u1_range,
        v1_range,
        &plane_sof_clone.geometry,
        u2_range,
        v2_range,
        &tolerance,
    ) {
        Ok(r) => r,
        Err(e) => {
            if trace {
                eprintln!("[topo_repair]   reject: intersect_surfaces failed: {e:?}");
            }
            return false;
        }
    };

    // Pick the intersection curve passing closest to both anchors.
    let mut best: Option<(NurbsCurve, NurbsCurve, NurbsCurve, f64, f64)> = None;
    let mut best_score = f64::INFINITY;
    for inter in &result.intersections {
        let SurfaceSurfaceIntersection::Curve(curve_on_surf) = inter else {
            continue;
        };
        let Some(t1) = curve_on_surf.curve_3d.closest_point(&p_anchor_in) else {
            continue;
        };
        let Some(t2) = curve_on_surf.curve_3d.closest_point(&p_anchor_out) else {
            continue;
        };
        let d1 = curve_on_surf
            .curve_3d
            .point_at(t1)
            .distance_to(&p_anchor_in);
        let d2 = curve_on_surf
            .curve_3d
            .point_at(t2)
            .distance_to(&p_anchor_out);
        let score = d1.max(d2);
        if score < best_score {
            best_score = score;
            best = Some((
                curve_on_surf.curve_3d.clone(),
                curve_on_surf.curve_uv1.clone(),
                curve_on_surf.curve_uv2.clone(),
                t1,
                t2,
            ));
        }
    }
    let Some((curve_3d, curve_uv1, curve_uv2, t1, t2)) = best else {
        if trace {
            eprintln!("[topo_repair]   reject: no intersection curve passes near anchors");
        }
        return false;
    };
    // Relaxed score tolerance based on anchor drift.
    let off_anchor_in = anchor_off_cyl_distance(&p_anchor_in, &cyl_sof.geometry);
    let off_anchor_out = anchor_off_cyl_distance(&p_anchor_out, &cyl_sof.geometry);
    let max_anchor_drift = off_anchor_in.max(off_anchor_out);
    let score_tol = (3.0 * max_anchor_drift).max(1e-6).min(chord_len * 0.5);
    if best_score > score_tol {
        if trace {
            eprintln!(
                "[topo_repair]   reject: best_score {best_score:.5} > tol {score_tol:.5} (drifts in={:.5}, out={:.5})",
                off_anchor_in, off_anchor_out
            );
        }
        return false;
    }
    // If the chord pair traverses the intersection curve in the
    // reversed parameter direction, slice [t2, t1] and reverse so
    // the resulting arc still runs anchor_in → anchor_out.
    let (lo, hi, reversed) = if t1 <= t2 {
        (t1, t2, false)
    } else {
        (t2, t1, true)
    };
    let (sliced_3d_raw, sliced_uv1_raw, sliced_uv2_raw) = match (
        slice_curve(&curve_3d, lo, hi).map(|c| if reversed { reverse_nurbs(&c) } else { c }),
        slice_curve(&curve_uv1, lo, hi).map(|c| if reversed { reverse_nurbs(&c) } else { c }),
        slice_curve(&curve_uv2, lo, hi).map(|c| if reversed { reverse_nurbs(&c) } else { c }),
    ) {
        (Some(a), Some(b), Some(c)) => (a, b, c),
        _ => {
            if trace {
                eprintln!("[topo_repair]   reject: slice_curve failed");
            }
            return false;
        }
    };

    // `split_face` hardcodes `t_start: 0.0, t_end: 1.0` on the
    // resulting `CurveOnEdge` regardless of the input curve's actual
    // parameter range. If the curve doesn't *actually* run [0, 1],
    // validate's `point_at(0.0)` will read out-of-range and produce
    // garbage — a 0.06+ ShellGap on the new arc edge. Remap each
    // sliced curve's knots to [0, 1] to keep things consistent.
    let sliced_3d = remap_to_unit_range(&sliced_3d_raw);
    let sliced_uv1 = remap_to_unit_range(&sliced_uv1_raw);
    let sliced_uv2 = remap_to_unit_range(&sliced_uv2_raw);

    // Snap the anchor *vertex positions* to the analytic intersection
    // curve's endpoints. The anchors drift a few mm off the analytic
    // cylinder/plane intersection due to numerics in corner blend
    // insertion; correcting them onto the analytic curve makes the
    // arc lie exactly on both surfaces (no ShellGap) without
    // distorting the curve geometry. Also update any incident edge
    // curves that end at the moved vertices (chord lines being the
    // typical case) so the loops stay closed.
    let (et_start, et_end) = sliced_3d.parameter_range();
    let analytic_start = sliced_3d.point_at(et_start);
    let analytic_end = sliced_3d.point_at(et_end);
    let old_anchor_in = brep.vertices[cand.anchor_in].position;
    let old_anchor_out = brep.vertices[cand.anchor_out].position;
    if let Some(v_in) = brep.vertices.get_mut(cand.anchor_in) {
        v_in.position = analytic_start;
    }
    if let Some(v_out) = brep.vertices.get_mut(cand.anchor_out) {
        v_out.position = analytic_end;
    }
    snap_incident_edge_curve_endpoints(brep, cand.anchor_in, old_anchor_in, analytic_start);
    snap_incident_edge_curve_endpoints(brep, cand.anchor_out, old_anchor_out, analytic_end);

    // Replace `sliced_uv1` (cyl UV from intersect_surfaces) with a
    // fresh degree-3 NURBS fit through dense samples projected onto
    // the cylinder. The intersect_surfaces output reuses the 3D
    // curve's degree-2 rational structure for the cyl UV pcurve,
    // which is only an approximation — the cylinder's nonlinear
    // (cos, sin) mapping doesn't commute with NURBS evaluation, so
    // a degree-2 cyl UV with the same knots/weights as the 3D
    // ellipse misses by up to O(radius) on the oblique case.
    let fitted_cyl_uv = match fit_cylinder_uv_pcurve(&sliced_3d, &cyl_sof.geometry) {
        Some(c) => c,
        None => {
            if trace {
                eprintln!("[topo_repair]   reject: failed to fit cyl UV pcurve");
            }
            return false;
        }
    };

    // Insert the new arc edge with `split_face`.
    let split_result = match brep.split_face(cyl_face, &fitted_cyl_uv, &sliced_3d, &tolerance) {
        Ok(s) => s,
        Err(e) => {
            if trace {
                eprintln!("[topo_repair]   reject: split_face failed: {e}");
            }
            return false;
        }
    };
    let _ = sliced_uv1;

    // Identify which sub-face is the sliver (the one whose loop
    // contains `off_cyl_vertex`).
    let (sliver_face, _main_face) =
        if face_loop_contains_vertex(brep, split_result.face_a, cand.off_cyl_vertex) {
            (split_result.face_a, split_result.face_b)
        } else if face_loop_contains_vertex(brep, split_result.face_b, cand.off_cyl_vertex) {
            (split_result.face_b, split_result.face_a)
        } else {
            // Couldn't tell — bail out (the split happened but we can't
            // safely complete the surgery).
            return false;
        };

    // Set the sliver's surface to the plane (it's geometrically on
    // the plane, not the cylinder). Re-project the sliver's outer
    // loop pcurves into the plane's UV space.
    //
    // We deliberately DO NOT kef the sliver into the adjacent
    // plane face. The sliver is bounded by two chord edges (both
    // shared with `plane_face`) plus the new arc; merging would
    // require two kefs but the second one runs against an interior
    // edge (both half-edges in the merged face), which is not a
    // valid Euler-op input. Leaving the sliver as a separate co-
    // planar face is BRep-valid: multiple faces may share a planar
    // surface, and OCCT-strict slicers accept this. The cylinder
    // face's outer loop now ends at the analytic arc instead of
    // the chord pair through the off-cylinder vertex.
    if let Some(plane_geom) = brep.faces[plane_face].surface.clone() {
        brep.faces[sliver_face].surface = Some(plane_geom);
    }

    // Re-project every pcurve on the sliver's outer loop (the new
    // arc + the two chord half-edges bordering the sliver) into
    // plane UV space. The arc pcurve assigned by `split_face` was
    // in cylinder UV, which would interpret as plane parameters and
    // produce a ShellGap on the order of cylinder radius.
    reproject_loop_pcurves_into_plane(brep, sliver_face);

    let _ = curve_uv2;
    let _ = sliced_uv2;
    true
}
