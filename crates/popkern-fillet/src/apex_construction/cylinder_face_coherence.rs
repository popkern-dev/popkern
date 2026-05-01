//! Phase 6c-A2.4 option 2 — cylinder face outer-loop coherence
//! detector and (eventually) repairer.
//!
//! Per ADR §25 the post-Phase-C BRep on variable-radius cube
//! fixtures (cube_mixed, cube_extreme) has cylinder faces whose
//! outer-loop vertices aren't all on the cylinder's stored
//! surface — Phase C's `merge_consumed_faces_any` "longest
//! descendant" heuristic picks the wrong consumed-face pair, so a
//! single cylinder face's loop walks across what would be other
//! fillets' natural trim segments.
//!
//! This module provides:
//!
//! * **Detection** — `find_off_cyl_loop_vertices` walks every
//!   cylinder face's outer loop and collects vertices not on the
//!   stored cylinder surface (within a 5% radius tolerance, the
//!   same threshold `topology_repair` uses).  Pure, BRep-immutable.
//! * **Reporting** — `IncoherentCylinderFace` carries the face
//!   handle, fillet radius, and the off-cyl vertex positions, so
//!   A2's dispatch can produce a precise failure-cause message
//!   rather than the generic "no shared face for mef" downstream.
//!
//! Repair (the actual `kef` + `split_edge` surgery) is a follow-up
//! slice — for now this module is detection + diagnostic only.
//! Wired into `a2_dispatch::try_a2_apex_construction` to surface
//! the precise root cause before the planner runs.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, VertexHandle};

use super::plan::ApexPlan;
use super::target_topology::TargetEdgeRole;

/// One cylinder face that has at least one outer-loop vertex off
/// its stored cylinder surface.
#[derive(Debug, Clone)]
pub(crate) struct IncoherentCylinderFace {
    pub face: FaceHandle,
    pub radius: f64,
    /// Vertices on the outer loop that are off-surface, with the
    /// distance from the vertex to the cylinder's axis (so callers
    /// can see the magnitude of the mismatch).
    pub off_cyl_vertices: Vec<(VertexHandle, Point3, f64)>,
}

/// Tolerance threshold for "on cylinder": distance from axis must
/// be within 5% of the cylinder's radius.  Matches the threshold
/// `topology_repair::vertex_on_cylinder` uses (proven adequate on
/// pyramid apex-only fillets per #136 work).
const RADIAL_TOL_FRACTION: f64 = 0.05;

/// Walk every cylinder face and collect its off-cyl outer-loop
/// vertices.  Returns the list of incoherent faces (may be empty
/// when the BRep state is already cylinder-coherent — uniform-radius
/// cube, for example).
pub(crate) fn find_off_cyl_loop_vertices(brep: &BRep) -> Vec<IncoherentCylinderFace> {
    let mut out = Vec::new();
    let face_handles: Vec<FaceHandle> = brep.faces.keys().collect();
    for fh in face_handles {
        let face = match brep.faces.get(fh) {
            Some(f) => f,
            None => continue,
        };
        let (cyl_origin, cyl_axis, cyl_radius) = match face.surface.as_ref().map(|s| &s.geometry) {
            Some(SurfaceGeometry::Cylinder {
                origin,
                axis,
                radius,
            }) => (*origin, *axis, *radius),
            _ => continue,
        };

        let outer_loop = face.outer_loop;
        let loop_data = match brep.loops.get(outer_loop) {
            Some(l) => l,
            None => continue,
        };
        let start_he = loop_data.half_edge;
        if !brep.half_edges.contains_key(start_he) {
            continue;
        }

        let mut off_cyl: Vec<(VertexHandle, Point3, f64)> = Vec::new();
        let mut cursor = start_he;
        let mut iters = 0;
        let max_iters = brep.half_edges.len() + 1;
        while let Some(he) = brep.half_edges.get(cursor) {
            let v = he.origin;
            let p = brep.vertices[v].position;
            let radial_dist = radial_distance_to_axis(&p, cyl_origin, cyl_axis);
            let tol = cyl_radius * RADIAL_TOL_FRACTION;
            if (radial_dist - cyl_radius).abs() > tol {
                off_cyl.push((v, p, radial_dist));
            }
            cursor = he.next;
            iters += 1;
            if cursor == start_he || iters >= max_iters {
                break;
            }
        }

        if !off_cyl.is_empty() {
            out.push(IncoherentCylinderFace {
                face: fh,
                radius: cyl_radius,
                off_cyl_vertices: off_cyl,
            });
        }
    }
    out
}

/// Compute the perpendicular distance from a 3D point to the
/// cylinder's axis line through `origin` along `axis`.  The point
/// is on the cylinder iff this distance equals `radius` (within
/// tolerance).
fn radial_distance_to_axis(p: &Point3, origin: Point3, axis: Vec3) -> f64 {
    let len = axis.length().max(1e-12);
    let axis_unit = Vec3::new(axis.x / len, axis.y / len, axis.z / len);
    let rel = Vec3::new(p.x - origin.x, p.y - origin.y, p.z - origin.z);
    let axial = rel.dot(&axis_unit);
    let radial = Vec3::new(
        rel.x - axial * axis_unit.x,
        rel.y - axial * axis_unit.y,
        rel.z - axial * axis_unit.z,
    );
    radial.length()
}

/// Format a list of incoherent cylinder faces into a single error
/// message suitable for `BlendWarning::A2FellBackToLegacy`'s
/// reason field.  Mentions the worst offender (largest radial
/// mismatch) and the total count.
pub(crate) fn describe_incoherence(report: &[IncoherentCylinderFace]) -> String {
    if report.is_empty() {
        return "no incoherent cylinder faces detected".to_string();
    }
    let mut total_off_cyl = 0usize;
    let mut worst_dist = 0.0f64;
    let mut worst_face: Option<FaceHandle> = None;
    let mut worst_radius = 0.0f64;
    for entry in report {
        for (_, _, d) in &entry.off_cyl_vertices {
            total_off_cyl += 1;
            let mismatch = (d - entry.radius).abs();
            if mismatch > worst_dist {
                worst_dist = mismatch;
                worst_face = Some(entry.face);
                worst_radius = entry.radius;
            }
        }
    }
    format!(
        "{} cylinder face(s) have {} off-cyl outer-loop vertex(ices); worst offender {:?} \
         (r={:.4}) has a vertex {:.4} m off-surface — Phase C `merge_consumed_faces_any` \
         heuristic produced an incoherent loop; A2.4 option 2 (rebuild cyl-face loops) \
         needed to proceed",
        report.len(),
        total_off_cyl,
        worst_face,
        worst_radius,
        worst_dist,
    )
}

// ─────────────────────────────────────────────────────────────────
// A2.4 option 2 surgery: stitch apex vertices into cylinder face
// outer loops via split-and-weld.
// ─────────────────────────────────────────────────────────────────

/// Pre-resolve apex_edge surgery: for each plan apex_edge whose
/// endpoints don't share a face's outer loop on the post-Phase-C
/// BRep, find the geometrically-natural cylinder face (whose
/// surface contains the apex_edge curve midpoint) and stitch
/// both endpoints into its outer loop via split-and-weld.
///
/// Returns the count of plan apex_edges where stitching ran (=
/// the number of cyl-side apex_edges that needed it; Gregory-side
/// edges are skipped because their curve midpoint isn't on any
/// cylinder face's surface).
///
/// Per ADR §25.4 option 2.  The split-and-weld primitive uses
/// `BRep::split_edge` to create a temporary vertex on the cyl
/// face's loop (at the parameter where the cyl loop edge is
/// closest to the apex vertex's true position), then
/// `BRep::weld_vertex_pair` to merge the temporary into the
/// existing apex vertex — which then appears as a half-edge
/// origin in the cyl face's outer loop, enabling
/// [`super::construct::find_shared_face_for_mef`] to succeed.
///
/// Important: the resulting edge curves on the loop may have
/// endpoints that don't match the apex vertex position (because
/// the curves were authored for the temporary's chord-midpoint
/// position).  Caller must rewrite those curves before
/// validation.  Today this is acceptable because A2's downstream
/// stages only need topological adjacency for `mef`.
pub(crate) fn stitch_apex_edges_into_cyl_loops(
    brep: &mut BRep,
    plan: &ApexPlan,
    apex_vertex_handles: &[VertexHandle],
) -> Result<usize, String> {
    let mut stitched = 0;
    for (i, target_edge) in plan.target_topology.apex_edges.iter().enumerate() {
        let v_origin = apex_vertex_handles
            .get(target_edge.origin_vertex_idx)
            .copied()
            .ok_or_else(|| format!("apex_edge[{i}] origin_vertex_idx out of range"))?;
        let v_dest = apex_vertex_handles
            .get(target_edge.dest_vertex_idx)
            .copied()
            .ok_or_else(|| format!("apex_edge[{i}] dest_vertex_idx out of range"))?;

        // Find the analytic face whose surface contains the curve
        // midpoint AND matches the planner's role hint (Cylinder
        // for CapVsCyl, Plane for CapVsPlane).  Without the role
        // filter, midpoints near a fillet's blend region can lie
        // within the 5% radial tolerance of two different surfaces
        // (e.g., on a cube_mixed bridge midpoint sits on both the
        // shared plane face AND the third fillet's cylinder), and
        // picking the wrong face produces a downstream `mef` on a
        // cylinder split when the planner expected a plane bridge.
        let (t_lo, t_hi) = target_edge.curve.parameter_range();
        let t_mid = (t_lo + t_hi) * 0.5;
        let curve_midpoint = target_edge.curve.point_at(t_mid);
        let role_filter = match target_edge.role {
            TargetEdgeRole::CapVsCyl { .. } => SurfaceRoleFilter::CylinderOnly,
            TargetEdgeRole::CapVsPlane { .. } => SurfaceRoleFilter::PlaneOnly,
            TargetEdgeRole::InnerSpoke => continue,
        };
        let target_face =
            match find_analytic_face_containing_position(brep, curve_midpoint, role_filter) {
                Some(f) => f,
                None => continue,
            };

        // Skip if both endpoints are already in target_face's loop —
        // any other shared face is irrelevant; the natural face for
        // mef is the one matching the planner's role.
        if vertex_is_in_face_outer_loop(brep, target_face, v_origin)
            && vertex_is_in_face_outer_loop(brep, target_face, v_dest)
        {
            continue;
        }

        // Stitch each endpoint into target_face's outer loop.
        stitch_apex_vertex_into_face_loop(brep, target_face, v_origin)
            .map_err(|e| format!("apex_edge[{i}] stitch v_origin: {e}"))?;
        stitch_apex_vertex_into_face_loop(brep, target_face, v_dest)
            .map_err(|e| format!("apex_edge[{i}] stitch v_dest: {e}"))?;

        stitched += 1;
    }
    Ok(stitched)
}

/// Surface-kind filter for `find_analytic_face_containing_position`.
#[derive(Debug, Clone, Copy)]
enum SurfaceRoleFilter {
    CylinderOnly,
    PlaneOnly,
}

/// Find any analytic face (Cylinder or Plane) whose stored surface
/// contains `position` AND whose surface kind matches `role`.
/// Cylinders are matched by radial distance (5% tolerance); planes
/// by signed-distance to the plane (1e-6 m).  Among multiple
/// candidates, prefers the face whose outer loop's vertices are
/// closest to `position` — important when several plane sub-faces
/// share the same geometric plane (Phase C's
/// `merge_consumed_faces_any` typically leaves multiple sub-faces
/// per geometric plane).
fn find_analytic_face_containing_position(
    brep: &BRep,
    position: Point3,
    role: SurfaceRoleFilter,
) -> Option<FaceHandle> {
    let plane_tol = 1e-6_f64;
    let mut best: Option<(FaceHandle, f64)> = None;
    for (fh, face) in brep.faces.iter() {
        let Some(sof) = face.surface.as_ref() else {
            continue;
        };
        let on_surface = match (&sof.geometry, role) {
            (
                SurfaceGeometry::Cylinder {
                    origin,
                    axis,
                    radius,
                },
                SurfaceRoleFilter::CylinderOnly,
            ) => {
                let radial = radial_distance_to_axis(&position, *origin, *axis);
                let tol = radius * RADIAL_TOL_FRACTION;
                (radial - radius).abs() <= tol
            }
            (SurfaceGeometry::Plane { origin, normal, .. }, SurfaceRoleFilter::PlaneOnly) => {
                let n_len = normal.length().max(1e-12);
                let nn = Vec3::new(normal.x / n_len, normal.y / n_len, normal.z / n_len);
                let d = (position.x - origin.x) * nn.x
                    + (position.y - origin.y) * nn.y
                    + (position.z - origin.z) * nn.z;
                d.abs() <= plane_tol
            }
            _ => false,
        };
        if !on_surface {
            continue;
        }
        // Score by minimum loop-vertex distance to `position`: when
        // several faces share a surface, prefer the one whose outer
        // loop is geometrically closest.
        let outer_loop = face.outer_loop;
        let Some(loop_data) = brep.loops.get(outer_loop) else {
            continue;
        };
        let start = loop_data.half_edge;
        let mut min_dist = f64::INFINITY;
        let mut cursor = start;
        let mut iters = 0;
        let max_iters = brep.half_edges.len() + 1;
        while let Some(he) = brep.half_edges.get(cursor) {
            let p = brep.vertices[he.origin].position;
            let d = position.distance_to(&p);
            if d < min_dist {
                min_dist = d;
            }
            cursor = he.next;
            iters += 1;
            if cursor == start || iters >= max_iters {
                break;
            }
        }
        match best {
            Some((_, b)) if b <= min_dist => {}
            _ => best = Some((fh, min_dist)),
        }
    }
    best.map(|(fh, _)| fh)
}

/// True iff `vertex` appears as a half-edge origin on `face`'s
/// outer loop.
fn vertex_is_in_face_outer_loop(brep: &BRep, face: FaceHandle, vertex: VertexHandle) -> bool {
    let Some(face_data) = brep.faces.get(face) else {
        return false;
    };
    let Some(loop_data) = brep.loops.get(face_data.outer_loop) else {
        return false;
    };
    let start = loop_data.half_edge;
    let mut cursor = start;
    let mut iters = 0;
    let max_iters = brep.half_edges.len() + 1;
    while let Some(he) = brep.half_edges.get(cursor) {
        if he.origin == vertex {
            return true;
        }
        cursor = he.next;
        iters += 1;
        if cursor == start || iters >= max_iters {
            return false;
        }
    }
    false
}

/// Stitch an existing apex `vertex` into `target_face`'s outer
/// loop via split-and-weld.  After this, the vertex appears as a
/// half-edge origin in the loop.
///
/// Algorithm:
/// 1. Walk the outer loop, find the loop edge whose 3D segment
///    has minimum perpendicular distance to the vertex position.
/// 2. Compute the segment-fraction projection `t ∈ (0.001, 0.999)`
///    of the vertex position onto that segment.
/// 3. `split_edge(target, t)` — creates a temporary vertex
///    `V_temp` on the loop at parameter t.
/// 4. `weld_vertex_pair(vertex, V_temp)` — half-edges that
///    originated from V_temp now originate from `vertex`;
///    V_temp is removed.
///
/// No-op when `vertex` is already in the loop.  Works for any face
/// kind (cylinder, plane, etc.) — the algorithm relies only on
/// loop-edge geometry, not the surface.
fn stitch_apex_vertex_into_face_loop(
    brep: &mut BRep,
    target_face: FaceHandle,
    vertex: VertexHandle,
) -> Result<(), String> {
    if vertex_is_in_face_outer_loop(brep, target_face, vertex) {
        return Ok(());
    }
    let v_pos = brep
        .vertices
        .get(vertex)
        .ok_or_else(|| format!("vertex {vertex:?} not found"))?
        .position;

    // Find the loop edge whose 3D line segment has minimum
    // perpendicular distance to v_pos.  Track (edge, t).
    let face = brep
        .faces
        .get(target_face)
        .ok_or_else(|| format!("target_face {target_face:?} not found"))?;
    let outer_loop = face.outer_loop;
    let loop_data = brep
        .loops
        .get(outer_loop)
        .ok_or_else(|| format!("loop {outer_loop:?} not found"))?;
    let start_he = loop_data.half_edge;

    let mut best: Option<(EdgeHandle, f64)> = None;
    let mut best_dist = f64::INFINITY;
    let mut cursor = start_he;
    let mut iters = 0;
    let max_iters = brep.half_edges.len() + 1;
    while let Some(he) = brep.half_edges.get(cursor) {
        let next_he = brep
            .half_edges
            .get(he.next)
            .ok_or_else(|| format!("next half-edge {:?} missing", he.next))?;
        let p_start = brep.vertices[he.origin].position;
        let p_end = brep.vertices[next_he.origin].position;
        let dx = p_end.x - p_start.x;
        let dy = p_end.y - p_start.y;
        let dz = p_end.z - p_start.z;
        let len_sq = dx * dx + dy * dy + dz * dz;
        if len_sq > 1e-18 {
            let raw_t = ((v_pos.x - p_start.x) * dx
                + (v_pos.y - p_start.y) * dy
                + (v_pos.z - p_start.z) * dz)
                / len_sq;
            let t = raw_t.clamp(0.001, 0.999);
            let foot = Point3::new(p_start.x + t * dx, p_start.y + t * dy, p_start.z + t * dz);
            let d = v_pos.distance_to(&foot);
            if d < best_dist {
                best_dist = d;
                best = Some((he.edge, t));
            }
        }
        cursor = he.next;
        iters += 1;
        if cursor == start_he || iters >= max_iters {
            break;
        }
    }

    let (target_edge, t) =
        best.ok_or_else(|| format!("no loop edge found near vertex {vertex:?}"))?;

    // Split the loop edge at the projected parameter, creating
    // a temporary vertex.
    let split = brep
        .split_edge(target_edge, t)
        .map_err(|e| format!("split_edge({target_edge:?}, {t}): {e:?}"))?;
    let temp_vertex = split.vertex;

    // Weld the temporary into the apex vertex.  After this, the
    // temp's outgoing half-edges all reroute to vertex; vertex
    // appears as an origin on cyl_face's outer loop.
    brep.weld_vertex_pair(vertex, temp_vertex).map_err(|e| {
        format!("weld_vertex_pair(primary={vertex:?}, secondary={temp_vertex:?}): {e:?}")
    })?;

    Ok(())
}
