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
//! Repair for the validated same-plane chord-pair class lives in
//! `construct::topology_repair` and runs in the post-blend pipeline
//! with snapshot rollback and `BRep::validate` postconditions. This
//! diagnostic remains wired into `a2_dispatch::try_a2_apex_construction`
//! because incoherence that survives that repair is a distinct
//! construction blocker, usually a multi-face ownership case rather
//! than a simple base-plane chord pair.

use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::topology::{BRep, FaceHandle, VertexHandle};

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
