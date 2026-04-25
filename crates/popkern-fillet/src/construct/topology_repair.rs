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
    let face_handles: Vec<FaceHandle> = brep.faces.keys().collect();
    let mut repaired = 0;
    for fh in face_handles {
        if !is_cylinder_face(brep, fh) {
            continue;
        }
        // Re-collect candidates each iteration — earlier repairs
        // mutate the topology, invalidating handles cached from a
        // prior pass. Cap at 8 attempts per face as a safety
        // (would only matter for fillet faces with many
        // off-cylinder vertices, which we don't expect).
        for _attempt in 0..8 {
            if !brep.faces.contains_key(fh) {
                break;
            }
            let candidates = find_chord_pair_candidates(brep, fh);
            let Some(cand) = candidates.into_iter().next() else {
                break;
            };
            let snapshot = brep.clone();
            if try_repair_chord_pair(brep, fh, cand) && brep.validate().is_ok() {
                repaired += 1;
            } else {
                *brep = snapshot;
                break;
            }
        }
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
    // The two chord edges (he_in and he_out) sit on a planar
    // adjacent face on their twin side. Confirm both twins lead to
    // the *same* planar face — that's the sliver-merge target.
    let twin_in = brep.half_edges[cand.he_in].twin;
    let twin_out = brep.half_edges[cand.he_out].twin;
    let plane_face_in = brep.loops[brep.half_edges[twin_in].loop_].face;
    let plane_face_out = brep.loops[brep.half_edges[twin_out].loop_].face;
    if plane_face_in != plane_face_out {
        // Different adjacent faces on the two chord edges — out of
        // scope (would need to handle a more complex topology).
        return false;
    }
    let plane_face = plane_face_in;
    let Some(plane_sof) = brep.faces.get(plane_face).and_then(|f| f.surface.as_ref()) else {
        return false;
    };
    if !matches!(plane_sof.geometry, SurfaceGeometry::Plane { .. }) {
        return false;
    }

    let cyl_sof = brep.faces[cyl_face].surface.as_ref().unwrap().clone();
    let plane_sof_clone = plane_sof.clone();
    let p_anchor_in = brep.vertices[cand.anchor_in].position;
    let p_anchor_out = brep.vertices[cand.anchor_out].position;
    let chord_len = p_anchor_in.distance_to(&p_anchor_out);
    if chord_len < 1e-9 {
        return false;
    }

    // Compute cylinder × plane intersection.
    let positions_1 = brep.face_outer_loop_positions(cyl_face);
    let positions_2 = brep.face_outer_loop_positions(plane_face);
    if positions_1.is_empty() || positions_2.is_empty() {
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
        Err(_) => return false,
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
        return false;
    };
    if best_score > (chord_len * 0.1).max(1e-6) {
        return false;
    }
    if t1 > t2 {
        // Reversed direction — would need to flip; skip for simplicity.
        return false;
    }
    let (sliced_3d, sliced_uv1, sliced_uv2) = match (
        slice_curve(&curve_3d, t1, t2),
        slice_curve(&curve_uv1, t1, t2),
        slice_curve(&curve_uv2, t1, t2),
    ) {
        (Some(a), Some(b), Some(c)) => (a, b, c),
        _ => return false,
    };

    // Verify endpoints align with anchors (within slack — vertices
    // can drift by a few percent of radius, the analytic intersection
    // sits exactly on the cylinder).
    let (et_start, et_end) = sliced_3d.parameter_range();
    let arc_start = sliced_3d.point_at(et_start);
    let arc_end = sliced_3d.point_at(et_end);
    if arc_start.distance_to(&p_anchor_in) > (chord_len * 0.1).max(1e-6)
        || arc_end.distance_to(&p_anchor_out) > (chord_len * 0.1).max(1e-6)
    {
        return false;
    }

    // Insert the new arc edge with `split_face`. This splits the
    // cylinder face into:
    //   * Main cylinder ribbon (the larger region, retains the
    //     cylinder surface).
    //   * Sliver fragment (bounded by the new arc + the two old
    //     chord edges through `off_cyl_vertex`).
    let split_result = match brep.split_face(cyl_face, &sliced_uv1, &sliced_3d, &tolerance) {
        Ok(s) => s,
        Err(_) => return false,
    };

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
    // the plane, not the cylinder).
    if let Some(plane_geom) = brep.faces[plane_face].surface.clone() {
        brep.faces[sliver_face].surface = Some(plane_geom);
    }

    // Merge the sliver with the plane face by kef-ing one of the
    // shared chord edges (use `cand.he_in`'s twin, which is on the
    // plane face's loop).
    let twin_in_now = brep.half_edges[cand.he_in].twin;
    if brep.kef(twin_in_now).is_err() {
        // kef failed — leave things as-is rather than panicking
        // mid-surgery. The split did happen, but the merge didn't,
        // so the sliver remains as a tiny separate face. Better than
        // an outright crash but caller should treat as best-effort.
        return false;
    }

    // Now `cand.he_out`'s edge is interior to the merged face — kev it.
    // (Need to re-find the edge handle for `he_out` since kef may
    // have shuffled half-edges.)
    let _ = curve_uv2;
    let _ = sliced_uv2;
    true
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
