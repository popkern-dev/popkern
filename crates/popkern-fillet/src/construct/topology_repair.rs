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

use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::curve::CurveGeometry;
use popkern_core::geometry::intersection::surface_surface::{
    SurfaceSurfaceIntersection, intersect_surfaces,
};
use popkern_core::geometry::nurbs_curve::NurbsCurve;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::uv_projection::tight_bounded_parameter_range;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle, VertexHandle};

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
    let face_handles: Vec<FaceHandle> = brep.faces.keys().collect();
    let mut repaired = 0;
    for fh in face_handles {
        if !is_cylinder_face(brep, fh) {
            continue;
        }
        if dump {
            dump_cylinder_face(brep, fh);
        }
        for _attempt in 0..8 {
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

fn is_cylinder_face(brep: &BRep, fh: FaceHandle) -> bool {
    brep.faces
        .get(fh)
        .and_then(|f| f.surface.as_ref())
        .is_some_and(|sof| matches!(sof.geometry, SurfaceGeometry::Cylinder { .. }))
}

#[derive(Debug, Clone, Copy)]
struct ChordPairCandidate {
    /// The off-cylinder vertex (the "elbow" of the chord pair).
    off_cyl_vertex: VertexHandle,
    /// Half-edge in the cylinder face's outer loop ending AT the
    /// off-cylinder vertex (incoming).
    he_in: HalfEdgeHandle,
    /// Half-edge in the cylinder face's outer loop starting FROM
    /// the off-cylinder vertex (outgoing).
    he_out: HalfEdgeHandle,
    /// Vertices flanking the off-cylinder vertex (anchor vertices
    /// on the cylinder surface).
    anchor_in: VertexHandle,
    anchor_out: VertexHandle,
}

fn find_chord_pair_candidates(brep: &BRep, fh: FaceHandle) -> Vec<ChordPairCandidate> {
    let Some(face) = brep.faces.get(fh) else {
        return Vec::new();
    };
    let outer_loop = face.outer_loop;
    let Some(loop_data) = brep.loops.get(outer_loop) else {
        return Vec::new();
    };
    let start_he = loop_data.half_edge;
    if !brep.half_edges.contains_key(start_he) {
        return Vec::new();
    }

    let cylinder = face.surface.as_ref().map(|sof| sof.geometry.clone());
    let Some(SurfaceGeometry::Cylinder {
        origin,
        axis,
        radius,
    }) = cylinder
    else {
        return Vec::new();
    };

    // Walk the outer loop and collect each vertex's incoming half-
    // edge along with the next half-edge (outgoing).
    let mut hes: Vec<HalfEdgeHandle> = Vec::new();
    let mut cursor = start_he;
    let mut safety = 64;
    loop {
        if safety == 0 {
            break;
        }
        safety -= 1;
        hes.push(cursor);
        let next = brep.half_edges[cursor].next;
        cursor = next;
        if cursor == start_he {
            break;
        }
    }
    if hes.len() < 3 {
        return Vec::new();
    }

    // For each vertex (the origin of each half-edge): if the vertex
    // is OFF the cylinder, it's a candidate elbow. The incoming
    // half-edge for that vertex is hes[(i + n - 1) % n] (the prev),
    // the outgoing is hes[i].
    let n = hes.len();
    let mut out = Vec::new();
    for i in 0..n {
        let he_out = hes[i];
        let he_in = hes[(i + n - 1) % n];
        let Some(out_data) = brep.half_edges.get(he_out) else {
            continue;
        };
        let Some(in_data) = brep.half_edges.get(he_in) else {
            continue;
        };
        let v = out_data.origin;
        let Some(v_data) = brep.vertices.get(v) else {
            continue;
        };
        let p = v_data.position;
        if vertex_on_cylinder(&p, origin, axis, radius) {
            continue;
        }
        let anchor_in = in_data.origin;
        let twin_out = out_data.twin;
        let Some(twin_data) = brep.half_edges.get(twin_out) else {
            continue;
        };
        let anchor_out = twin_data.origin;
        out.push(ChordPairCandidate {
            off_cyl_vertex: v,
            he_in,
            he_out,
            anchor_in,
            anchor_out,
        });
    }
    out
}

fn vertex_on_cylinder(
    p: &Point3,
    origin: Point3,
    axis: popkern_core::geometry::Vec3,
    radius: f64,
) -> bool {
    let rel = popkern_core::geometry::Vec3::new(p.x - origin.x, p.y - origin.y, p.z - origin.z);
    let axial = rel.dot(&axis) / axis.length().max(1e-12);
    let axis_unit = if axis.length() > 1e-12 {
        popkern_core::geometry::Vec3::new(
            axis.x / axis.length(),
            axis.y / axis.length(),
            axis.z / axis.length(),
        )
    } else {
        return false;
    };
    let radial = popkern_core::geometry::Vec3::new(
        rel.x - axial * axis_unit.x,
        rel.y - axial * axis_unit.y,
        rel.z - axial * axis_unit.z,
    );
    let r = radial.length();
    (r - radius).abs() <= radius * 0.05
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

    let cyl_sof = brep.faces[cyl_face].surface.as_ref().unwrap().clone();
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

/// Surgery for the different-adjacent-planes case: the off-cyl
/// vertex `cand.off_cyl_vertex` sits at the intersection of cyl,
/// plane_A, and plane_B (the three surfaces meeting at this corner).
/// Currently drifted off-cyl by a few % of radius; snap it onto the
/// analytic intersection of `cyl ∩ plane_A` and `cyl ∩ plane_B`,
/// then rewrite each chord edge's 3D curve to the proper analytic
/// arc on its respective plane.
///
/// No topology change — only vertex position + edge curve rewrite.
fn try_repair_chord_pair_different_planes(
    brep: &mut BRep,
    cyl_face: FaceHandle,
    cand: ChordPairCandidate,
    plane_face_a: FaceHandle,
    plane_face_b: FaceHandle,
    trace: bool,
) -> bool {
    let cyl_geom = brep.faces[cyl_face]
        .surface
        .as_ref()
        .unwrap()
        .geometry
        .clone();
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

/// Compute the analytic cyl×plane intersection, returning the first
/// `Curve` intersection from `intersect_surfaces`. Handles the
/// fillet-specific tangent case (cyl × slant plane = contact line)
/// that `case_parallel` rejects when `d ≈ cr`.
fn analytic_cyl_plane_intersection(
    cyl: &SurfaceGeometry,
    plane: &SurfaceGeometry,
    tolerance: &popkern_core::tolerance::Tolerance,
) -> Option<NurbsCurve> {
    // Cyl UV: u ∈ [0, 2π], v unbounded — give the solver a generous v
    // range. Plane UV: use a reasonably wide square. The analytic
    // plane×cylinder solver doesn't actually clip the result to these
    // ranges, so generous bounds are safe.
    let result = intersect_surfaces(
        cyl,
        (0.0, std::f64::consts::TAU),
        (-10.0, 10.0),
        plane,
        (-10.0, 10.0),
        (-10.0, 10.0),
        tolerance,
    )
    .ok()?;
    for inter in &result.intersections {
        if let SurfaceSurfaceIntersection::Curve(c) = inter {
            return Some(c.curve_3d.clone());
        }
    }
    // Empty result — could be a tangent case (cyl perpendicular to
    // plane normal AND axis distance equals radius). Compute the
    // contact line manually.
    cyl_plane_contact_line(cyl, plane, tolerance)
}

/// Compute the contact line where a cylinder is tangent to a plane.
/// Returns a NURBS line along the cylinder axis at the contact angle.
fn cyl_plane_contact_line(
    cyl: &SurfaceGeometry,
    plane: &SurfaceGeometry,
    tolerance: &popkern_core::tolerance::Tolerance,
) -> Option<NurbsCurve> {
    let SurfaceGeometry::Cylinder {
        origin: co,
        axis: ca,
        radius: cr,
    } = cyl
    else {
        return None;
    };
    let SurfaceGeometry::Plane {
        origin: po,
        normal: pn,
        ..
    } = plane
    else {
        return None;
    };
    // Tangent requires axis perpendicular to plane normal (within
    // angular tolerance) AND axial-distance from cyl axis to plane
    // ≈ cyl radius.
    let axis_dot_normal = ca.dot(pn);
    let ca_len = ca.length();
    let pn_len = pn.length();
    if (axis_dot_normal / (ca_len * pn_len)).abs() > tolerance.angular() {
        return None; // not parallel
    }
    let axial_dist = ((*co - *po).dot(pn)) / pn_len;
    if (axial_dist.abs() - cr).abs() > cr * 0.05 {
        return None;
    }
    // Direction from cyl axis toward plane (negative of normal·sign).
    let pn_unit = popkern_core::geometry::Vec3::new(pn.x / pn_len, pn.y / pn_len, pn.z / pn_len);
    let toward_plane_sign = -axial_dist.signum();
    let radial = popkern_core::geometry::Vec3::new(
        toward_plane_sign * pn_unit.x * *cr,
        toward_plane_sign * pn_unit.y * *cr,
        toward_plane_sign * pn_unit.z * *cr,
    );
    let line_origin = Point3::new(co.x + radial.x, co.y + radial.y, co.z + radial.z);
    let v_lo = -10.0;
    let v_hi = 10.0;
    let a = Point3::new(
        line_origin.x + ca.x * v_lo,
        line_origin.y + ca.y * v_lo,
        line_origin.z + ca.z * v_lo,
    );
    let b = Point3::new(
        line_origin.x + ca.x * v_hi,
        line_origin.y + ca.y * v_hi,
        line_origin.z + ca.z * v_hi,
    );
    Some(NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![a, b],
        weights: vec![1.0, 1.0],
    })
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
        edge_mut.curve = Some(popkern_core::geometry::curve::CurveOnEdge {
            geometry: CurveGeometry::Nurbs(oriented_sliced),
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

/// Re-project every half-edge pcurve on `face`'s outer loop into
/// the face's surface UV space. Used after retargeting a face's
/// surface (cylinder → plane, etc.) so the pcurves stay consistent.
fn reproject_loop_pcurves_into_plane(brep: &mut BRep, face: FaceHandle) {
    let Some(face_data) = brep.faces.get(face) else {
        return;
    };
    let outer_loop = face_data.outer_loop;
    let plane_geom = match face_data.surface.as_ref().map(|s| &s.geometry) {
        Some(SurfaceGeometry::Plane {
            origin,
            u_axis,
            v_axis,
            ..
        }) => (*origin, *u_axis, *v_axis),
        _ => return,
    };
    let (plane_origin, plane_u, plane_v) = plane_geom;

    let Some(loop_data) = brep.loops.get(outer_loop) else {
        return;
    };
    let start_he = loop_data.half_edge;
    let mut hes: Vec<HalfEdgeHandle> = Vec::new();
    let mut cursor = start_he;
    let mut safety = 64;
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
        if cursor == start_he {
            break;
        }
    }

    for he in hes {
        // Project the underlying edge's 3D curve into plane UV.
        let edge = brep.half_edges[he].edge;
        let curve_3d = match brep.edges.get(edge).and_then(|e| e.curve.as_ref()) {
            Some(c) => c.geometry.clone(),
            None => continue,
        };
        let new_pcurve =
            match project_curve_into_plane_uv(&curve_3d, plane_origin, plane_u, plane_v) {
                Some(p) => p,
                None => continue,
            };
        // Half-edge orientation: if this half-edge runs the same
        // way as the underlying edge, use the projected pcurve as-is;
        // otherwise reverse.
        let edge_he = brep.edges[edge].half_edge;
        let pcurve = if he == edge_he {
            new_pcurve
        } else {
            reverse_nurbs(&new_pcurve)
        };
        brep.half_edges[he].pcurve = Some(pcurve.into());
    }
}

/// Project a 3D curve into a plane's (u, v) parameter space.
/// Plane parameterised as `origin + u·u_axis + v·v_axis`. Returns
/// the pcurve as a NURBS where each control point is `(u, v, 0)`.
fn project_curve_into_plane_uv(
    curve: &CurveGeometry,
    plane_origin: Point3,
    plane_u: popkern_core::geometry::Vec3,
    plane_v: popkern_core::geometry::Vec3,
) -> Option<NurbsCurve> {
    let nurbs = match curve {
        CurveGeometry::Nurbs(n) => n.clone(),
        // For non-NURBS curves, sample and fit. Skip for the first
        // pass — analytic projection is what we need.
        _ => return None,
    };
    let projected_cps: Vec<Point3> = nurbs
        .control_points
        .iter()
        .map(|p| {
            let rel = popkern_core::geometry::Vec3::new(
                p.x - plane_origin.x,
                p.y - plane_origin.y,
                p.z - plane_origin.z,
            );
            let u = rel.dot(&plane_u) / plane_u.length().max(1e-12).powi(2) * plane_u.length();
            let v = rel.dot(&plane_v) / plane_v.length().max(1e-12).powi(2) * plane_v.length();
            // The plane is parameterised by point = origin + u·u_axis + v·v_axis,
            // so (u, v) = (rel · û, rel · v̂) only when u_axis, v_axis are unit.
            // For the general case, dot products give projection lengths and we
            // need to divide by axis-length squared to get the parameter.
            let _ = (u, v);
            let u_param = rel.dot(&plane_u) / plane_u.dot(&plane_u);
            let v_param = rel.dot(&plane_v) / plane_v.dot(&plane_v);
            Point3::new(u_param, v_param, 0.0)
        })
        .collect();
    Some(NurbsCurve {
        degree: nurbs.degree,
        knots: nurbs.knots,
        control_points: projected_cps,
        weights: nurbs.weights,
    })
}

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

fn face_loop_contains_vertex(brep: &BRep, face: FaceHandle, vertex: VertexHandle) -> bool {
    let Some(f) = brep.faces.get(face) else {
        return false;
    };
    let Some(loop_data) = brep.loops.get(f.outer_loop) else {
        return false;
    };
    let start_he = loop_data.half_edge;
    let mut cursor = start_he;
    let mut safety = 64;
    loop {
        if safety == 0 {
            return false;
        }
        safety -= 1;
        let Some(he) = brep.half_edges.get(cursor) else {
            return false;
        };
        if he.origin == vertex {
            return true;
        }
        cursor = he.next;
        if cursor == start_he {
            return false;
        }
    }
}

fn slice_curve(curve: &NurbsCurve, lo: f64, hi: f64) -> Option<NurbsCurve> {
    let (full_min, full_max) = curve.parameter_range();
    if !(full_min..=full_max).contains(&lo) || !(full_min..=full_max).contains(&hi) {
        return None;
    }
    if hi <= lo + 1e-12 {
        return None;
    }
    let (_, right) = curve.split_at(lo).ok()?;
    let (rem_min, _) = right.parameter_range();
    let cut = rem_min + (hi - lo);
    let (left, _) = right.split_at(cut).ok()?;
    Some(left)
}

#[allow(dead_code)]
fn _force_use(_e: EdgeHandle) {}

fn dump_cylinder_face(brep: &BRep, fh: FaceHandle) {
    let Some(face) = brep.faces.get(fh) else {
        return;
    };
    let cyl = face.surface.as_ref().map(|sof| &sof.geometry);
    let Some(SurfaceGeometry::Cylinder {
        origin,
        axis,
        radius,
    }) = cyl
    else {
        return;
    };
    let outer_loop = face.outer_loop;
    let Some(loop_data) = brep.loops.get(outer_loop) else {
        return;
    };
    let start_he = loop_data.half_edge;
    eprintln!(
        "[topo_dump] cyl face {fh:?} cyl(origin=({:.3},{:.3},{:.3}) axis=({:.3},{:.3},{:.3}) r={:.4}):",
        origin.x, origin.y, origin.z, axis.x, axis.y, axis.z, radius
    );
    let mut cursor = start_he;
    let mut safety = 64;
    loop {
        if safety == 0 {
            break;
        }
        safety -= 1;
        let Some(he) = brep.half_edges.get(cursor) else {
            break;
        };
        let v = he.origin;
        let p = brep.vertices[v].position;
        let twin_face = brep.loops[brep.half_edges[he.twin].loop_].face;
        let twin_surf = brep
            .faces
            .get(twin_face)
            .and_then(|f| f.surface.as_ref())
            .map(|s| match &s.geometry {
                SurfaceGeometry::Plane { .. } => "Plane",
                SurfaceGeometry::Cylinder { .. } => "Cyl",
                SurfaceGeometry::Sphere { .. } => "Sph",
                SurfaceGeometry::Cone { .. } => "Cone",
                SurfaceGeometry::Torus { .. } => "Torus",
                SurfaceGeometry::Nurbs(_) => "Nurbs",
                SurfaceGeometry::Gregory(_) => "Gregory",
            })
            .unwrap_or("?");
        let off = off_cyl_distance(&p, *origin, *axis, *radius);
        let on = off.abs() <= radius * 0.05;
        eprintln!(
            "[topo_dump]   v{:?} p=({:+.4},{:+.4},{:+.4}) off_cyl={:+.5} {} twin_face={twin_face:?}({twin_surf})",
            v,
            p.x,
            p.y,
            p.z,
            off,
            if on { "ON" } else { "OFF" }
        );
        cursor = he.next;
        if cursor == start_he {
            break;
        }
    }
}

fn off_cyl_distance(
    p: &Point3,
    origin: Point3,
    axis: popkern_core::geometry::Vec3,
    radius: f64,
) -> f64 {
    let len = axis.length().max(1e-12);
    let axis_unit = popkern_core::geometry::Vec3::new(axis.x / len, axis.y / len, axis.z / len);
    let rel = popkern_core::geometry::Vec3::new(p.x - origin.x, p.y - origin.y, p.z - origin.z);
    let axial = rel.dot(&axis_unit);
    let radial = popkern_core::geometry::Vec3::new(
        rel.x - axial * axis_unit.x,
        rel.y - axial * axis_unit.y,
        rel.z - axial * axis_unit.z,
    );
    radial.length() - radius
}

fn anchor_off_cyl_distance(p: &Point3, surf: &SurfaceGeometry) -> f64 {
    if let SurfaceGeometry::Cylinder {
        origin,
        axis,
        radius,
    } = surf
    {
        off_cyl_distance(p, *origin, *axis, *radius).abs()
    } else {
        0.0
    }
}

/// Linearly rescale a NURBS curve's knot vector so its parameter
/// range becomes [0, 1]. The curve geometry is unchanged — only
/// the parameterisation. Required because `split_face` hardcodes
/// `CurveOnEdge { t_start: 0.0, t_end: 1.0, .. }` for the new
/// arc edge, regardless of the input curve's actual range.
fn remap_to_unit_range(curve: &NurbsCurve) -> NurbsCurve {
    let (t_min, t_max) = curve.parameter_range();
    if (t_max - t_min).abs() < 1e-15 {
        return curve.clone();
    }
    let new_knots: Vec<f64> = curve
        .knots
        .iter()
        .map(|&k| (k - t_min) / (t_max - t_min))
        .collect();
    NurbsCurve {
        degree: curve.degree,
        knots: new_knots,
        control_points: curve.control_points.clone(),
        weights: curve.weights.clone(),
    }
}

/// Sample `curve_3d` at N parameter values, project each sample
/// onto the cylinder's (u, v) space using `closest_point`, then
/// fit a degree-3 NURBS interpolant. This produces a pcurve that —
/// unlike `cylinder_uv_curve` from `plane_cylinder.rs` — preserves
/// the cylinder mapping at every interior point of the arc, not
/// just at the control points. Required for ShellGap-clean output.
fn fit_cylinder_uv_pcurve(curve_3d: &NurbsCurve, cyl: &SurfaceGeometry) -> Option<NurbsCurve> {
    use popkern_core::geometry::fitting::fit_curve;
    use popkern_core::geometry::surface::SurfaceEval;

    let SurfaceGeometry::Cylinder { .. } = cyl else {
        return None;
    };

    let n_samples = 17;
    let (t_min, t_max) = curve_3d.parameter_range();
    let mut uv_points: Vec<Point3> = Vec::with_capacity(n_samples);
    let mut params: Vec<f64> = Vec::with_capacity(n_samples);

    let mut prev_u: Option<f64> = None;
    for i in 0..n_samples {
        let s = i as f64 / (n_samples - 1) as f64;
        let t = t_min + s * (t_max - t_min);
        let p = curve_3d.point_at(t);
        let (mut u, v) = cyl.closest_point(&p)?;
        // Unwrap u so the sample sequence is continuous (no TAU jumps).
        if let Some(prev) = prev_u {
            let tau = std::f64::consts::TAU;
            while u - prev > tau * 0.5 {
                u -= tau;
            }
            while u - prev < -tau * 0.5 {
                u += tau;
            }
        }
        prev_u = Some(u);
        uv_points.push(Point3::new(u, v, 0.0));
        params.push(s);
    }

    let degree = 3;
    let num_cps = degree + 2;
    fit_curve(&uv_points, &params, degree, num_cps).ok()
}

/// For every edge incident to `vh`, snap the 3D curve's endpoint
/// (the one matching v's t_start or t_end) to `new_position`.
/// More robust than the legacy version: uses half-edge orientation
/// to determine which end of the curve corresponds to `vh`, then
/// modifies the appropriate control point of the curve.
fn snap_incident_edge_curve_endpoints_robust(
    brep: &mut BRep,
    vh: VertexHandle,
    new_position: Point3,
) {
    let edge_handles: Vec<EdgeHandle> = brep.edges.keys().collect();
    for edge in edge_handles {
        let Some(edge_data) = brep.edges.get(edge) else {
            continue;
        };
        let he = edge_data.half_edge;
        let Some(he_data) = brep.half_edges.get(he) else {
            continue;
        };
        let twin = he_data.twin;
        let Some(twin_data) = brep.half_edges.get(twin) else {
            continue;
        };
        let v_origin = he_data.origin;
        let v_target = twin_data.origin;
        let edge_starts_at_vh = v_origin == vh;
        let edge_ends_at_vh = v_target == vh;
        if !edge_starts_at_vh && !edge_ends_at_vh {
            continue;
        }
        let Some(curve_on_edge) = brep.edges.get(edge).and_then(|e| e.curve.clone()) else {
            continue;
        };
        // Convert any analytic curve to NURBS first, so we have
        // control points to manipulate. Line/Arc/Circle types are
        // parameterised by analytic data (origin, direction, etc.)
        // — modifying those to relocate one endpoint is fiddly, but
        // promoting to NURBS gives us a uniform interface.
        let nurbs_full = match &curve_on_edge.geometry {
            CurveGeometry::Nurbs(n) => n.clone(),
            CurveGeometry::Line { .. }
            | CurveGeometry::Arc { .. }
            | CurveGeometry::Circle { .. }
            | CurveGeometry::Ellipse { .. }
            | CurveGeometry::Parabola { .. } => {
                let trange = (curve_on_edge.t_start, curve_on_edge.t_end);
                curve_on_edge.geometry.to_nurbs_in_range(trange)
            }
        };
        // CurveOnEdge.t_start/t_end may not match the underlying
        // NURBS's full parameter range — the curve may extend
        // beyond the edge's actual span. Slice down to [t_start,
        // t_end] so the first/last CPs of the resulting curve
        // *do* correspond to the edge endpoints.
        let nurbs_sliced =
            match slice_curve(&nurbs_full, curve_on_edge.t_start, curve_on_edge.t_end) {
                Some(s) => remap_to_unit_range(&s),
                None => nurbs_full.clone(),
            };
        // The curve direction may not match the half-edge direction
        // (edge.half_edge isn't always the curve-start half-edge).
        // So decide which CP to update based on which curve endpoint
        // is closer to the existing vertex position. Use the new
        // position, since the move can be O(drift) in either
        // direction — closer in original space stays closer.
        let mut cps = nurbs_sliced.control_points.clone();
        let snap_first = if let (Some(first), Some(last)) = (cps.first(), cps.last()) {
            first.distance_to(&new_position) <= last.distance_to(&new_position)
        } else {
            true
        };
        // Sanity: only update if exactly one end is incident to vh
        // (BRep half-edge semantics usually guarantee this).
        let _ = (edge_starts_at_vh, edge_ends_at_vh);
        if snap_first {
            if let Some(first) = cps.first_mut() {
                *first = new_position;
            }
        } else if let Some(last) = cps.last_mut() {
            *last = new_position;
        }
        let new_curve = NurbsCurve {
            degree: nurbs_sliced.degree,
            knots: nurbs_sliced.knots,
            control_points: cps,
            weights: nurbs_sliced.weights,
        };
        let (new_t_start, new_t_end) = new_curve.parameter_range();
        if let Some(edge_mut) = brep.edges.get_mut(edge) {
            edge_mut.curve = Some(popkern_core::geometry::curve::CurveOnEdge {
                geometry: CurveGeometry::Nurbs(new_curve),
                t_start: new_t_start,
                t_end: new_t_end,
            });
        }
    }
}

/// For every edge with an endpoint vertex `vh`, replace the
/// matching control point of the edge's 3D curve so the curve
/// terminates at the new vertex position. Identifies the matching
/// control point by proximity to `old_position`.
#[allow(dead_code)]
fn snap_incident_edge_curve_endpoints(
    brep: &mut BRep,
    vh: VertexHandle,
    old_position: Point3,
    new_position: Point3,
) {
    let edge_handles: Vec<EdgeHandle> = brep.edges.keys().collect();
    for edge in edge_handles {
        let Some(edge_data) = brep.edges.get(edge) else {
            continue;
        };
        let he = edge_data.half_edge;
        let Some(he_data) = brep.half_edges.get(he) else {
            continue;
        };
        let twin = he_data.twin;
        let Some(twin_data) = brep.half_edges.get(twin) else {
            continue;
        };
        let v_origin = he_data.origin;
        let v_target = twin_data.origin;
        if v_origin != vh && v_target != vh {
            continue;
        }
        let Some(curve_on_edge) = brep.edges.get(edge).and_then(|e| e.curve.clone()) else {
            continue;
        };
        let CurveGeometry::Nurbs(nurbs) = curve_on_edge.geometry else {
            continue;
        };
        let mut cps = nurbs.control_points.clone();
        let first_d = cps
            .first()
            .map(|p| p.distance_to(&old_position))
            .unwrap_or(f64::INFINITY);
        let last_d = cps
            .last()
            .map(|p| p.distance_to(&old_position))
            .unwrap_or(f64::INFINITY);
        if first_d < last_d && first_d < 1e-3 {
            if let Some(first) = cps.first_mut() {
                *first = new_position;
            }
        } else if last_d < 1e-3 {
            if let Some(last) = cps.last_mut() {
                *last = new_position;
            }
        } else {
            continue;
        }
        let new_curve = NurbsCurve {
            degree: nurbs.degree,
            knots: nurbs.knots,
            control_points: cps,
            weights: nurbs.weights,
        };
        if let Some(edge_mut) = brep.edges.get_mut(edge) {
            edge_mut.curve = Some(popkern_core::geometry::curve::CurveOnEdge {
                geometry: CurveGeometry::Nurbs(new_curve),
                t_start: curve_on_edge.t_start,
                t_end: curve_on_edge.t_end,
            });
        }
    }
}
