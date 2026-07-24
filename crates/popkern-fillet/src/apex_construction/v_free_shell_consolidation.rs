//! V-free vestigial support-shell consolidation (PK-153).
//!
//! The V-free repeated-apex rewrite authors its final support faces as
//! face-family bridges but leaves the OLD support faces behind as
//! orphaned flaps sealed by surfaceless membranes, and the strip far
//! ends behind removed-sliver patches — a glued old shell that
//! double-covers the support planes, re-encloses removed material
//! (body volume reads ≈ the uncut solid), and ships surfaceless
//! membranes to serving. Design + measured evidence:
//! docs/architecture/v-free-support-shell-consolidation.md.
//!
//! This pass retires that shell structurally, then welds the exposed
//! coincident boundary sides:
//! 1. FLAPS — a Production face that interior-overlaps a same-plane
//!    Production partner and whose every side twins either that
//!    partner or a ConstructionResidual face is vestigial: delete it.
//! 2. SLIVERS — once flaps fired (the V-free signature), a Production
//!    plane face with ≥2 ConstructionResidual-twinned sides whose
//!    remaining twins are all non-plane Production faces is a
//!    removed-material sliver at a strip far end: delete it.
//! 3. The existing unglued-residual scaffolding retirement then
//!    removes the membranes those deletions unglue.
//! 4. WINDING — the rewritten shell is authored mirror-wound against
//!    the untouched body faces (every old↔new coincident seam pairs
//!    PARALLEL, unweldable); the reconciliation stage 2-colors faces
//!    over seam parity, reverses the mirror class by twin-exchange,
//!    and rebuilds the rings of faces contorted across both classes.
//!    See v_free_shell_orientation.rs.
//! 5. WELD — dead boundary sides with coincident endpoints are
//!    re-twinned pairwise (duplicate edge records unify; curve
//!    authority prefers the non-linear side), and coincident vertex
//!    records merge.
//!
//! Transactional: the pass validates the consolidated model and rolls
//! back entirely on failure. It is self-guarding — bodies without the
//! flap signature (every non-V-free path) see a no-op.

use std::collections::BTreeMap;

use popkern_core::geometry::Point3;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::topology::{BRep, FaceHandle, HalfEdgeHandle, TopologyRole, VertexHandle};

#[path = "v_free_shell_orientation.rs"]
mod shell_orientation;
#[path = "v_free_shell_weld.rs"]
mod shell_weld;
pub(in crate::apex_construction) use shell_weld::{
    merge_coincident_seam_vertices, weld_dead_coincident_sides,
};

/// Structured counts for the consolidation.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct VestigialShellConsolidationReport {
    /// Same-plane double-cover flap faces removed.
    pub flap_faces: usize,
    /// Removed-material sliver faces removed.
    pub sliver_faces: usize,
    /// Residual scaffolding faces removed after ungluing.
    pub scaffolding_faces: usize,
    /// Dead boundary side pairs welded.
    pub welded_side_pairs: usize,
    /// Mirror-class faces whose loops were reversed to weldable parity.
    pub reversed_faces: usize,
    /// Contorted interface faces whose rings were rebuilt.
    pub rebuilt_rings: usize,
    /// Coincident vertex records merged.
    pub merged_vertices: usize,
    /// True when the consolidated model failed validation and the
    /// input was restored.
    pub rolled_back: bool,
}

impl VestigialShellConsolidationReport {
    /// True when the pass changed the model.
    pub fn applied(&self) -> bool {
        !self.rolled_back
            && (self.flap_faces > 0
                || self.sliver_faces > 0
                || self.welded_side_pairs > 0
                || self.reversed_faces > 0)
    }
}

/// Run the consolidation. No-op unless the flap signature is present.
pub(crate) fn consolidate_vestigial_support_shell(
    brep: &mut BRep,
) -> VestigialShellConsolidationReport {
    let mut report = VestigialShellConsolidationReport::default();
    let flaps = detect_flaps(brep);
    if flaps.is_empty() {
        return report;
    }
    let rollback = brep.clone();

    for face in &flaps {
        remove_face(brep, *face);
        report.flap_faces += 1;
    }
    for face in detect_slivers(brep) {
        remove_face(brep, face);
        report.sliver_faces += 1;
    }
    let scaffolding =
        crate::residual_scaffolding::retire_unglued_construction_residual_scaffolding(brep);
    report.scaffolding_faces = scaffolding.retired_faces;

    let winding = shell_orientation::reconcile_shell_winding(brep);
    report.reversed_faces = winding.reversed_faces;
    report.rebuilt_rings = winding.rebuilt_rings;

    let (welded, mut seam_vertices) =
        weld_dead_coincident_sides(brep, &std::collections::BTreeSet::new());
    report.welded_side_pairs = welded;
    seam_vertices.extend(winding.seam_vertices);
    report.merged_vertices =
        merge_coincident_seam_vertices(brep, &seam_vertices, &std::collections::BTreeSet::new())
            .len();
    repair_vertex_outgoing(brep);

    // Heal pcurves on faces that gained welded sides, reconcile
    // periodic branch choices (a reversed or welded strip loop can
    // carry a long-way arc image whose net UV winding is ±period —
    // the full-drum area/tessellation misread), then validate.
    let faces: Vec<_> = brep.faces.keys().collect();
    for face in &faces {
        let _ = brep.recompute_face_pcurves(*face);
    }
    brep.fix_stale_pcurves();
    for face in &faces {
        brep.reconcile_periodic_pcurve_branches(*face);
    }

    if let Err(err) = brep.validate() {
        if crate::direct_apex_trace::enabled() {
            eprintln!("[direct-apex-trace] vestigial shell consolidation rolled back: {err:?}");
        }
        *brep = rollback;
        return VestigialShellConsolidationReport {
            rolled_back: true,
            ..Default::default()
        };
    }
    report
}

fn face_sides(brep: &BRep, face: FaceHandle) -> Vec<HalfEdgeHandle> {
    let Some(face_data) = brep.faces.get(face) else {
        return Vec::new();
    };
    let mut sides = Vec::new();
    let mut loops = vec![face_data.outer_loop];
    loops.extend(face_data.inner_loops.iter().copied());
    for loop_handle in loops {
        let Some(loop_data) = brep.loops.get(loop_handle) else {
            continue;
        };
        let start = loop_data.half_edge;
        let mut current = start;
        for _ in 0..brep.half_edges.len() {
            sides.push(current);
            let Some(data) = brep.half_edges.get(current) else {
                break;
            };
            current = data.next;
            if current == start {
                break;
            }
        }
    }
    sides
}

fn twin_face(brep: &BRep, side: HalfEdgeHandle) -> Option<FaceHandle> {
    let twin = brep.half_edges.get(side)?.twin;
    let loop_ = brep.half_edges.get(twin)?.loop_;
    Some(brep.loops.get(loop_)?.face)
}

fn plane_key(brep: &BRep, face: FaceHandle) -> Option<(i64, i64, i64, i64)> {
    let sof = brep.faces.get(face)?.surface.as_ref()?;
    let SurfaceGeometry::Plane { origin, normal, .. } = &sof.geometry else {
        return None;
    };
    let n = normal.normalized()?;
    let d = origin.x * n.x + origin.y * n.y + origin.z * n.z;
    // Sign-canonical quantized key.
    let (n, d) = if (n.x, n.y, n.z) < (0.0, 0.0, 0.0) {
        (-n, -d)
    } else {
        (n, d)
    };
    let q = |v: f64| (v * 1.0e7).round() as i64;
    Some((q(n.x), q(n.y), q(n.z), q(d)))
}

/// A flap: Production face whose plane hosts an interior-overlapping
/// Production partner and whose sides all twin that partner or
/// ConstructionResidual faces.
fn detect_flaps(brep: &BRep) -> Vec<FaceHandle> {
    let mut by_plane: BTreeMap<(i64, i64, i64, i64), Vec<FaceHandle>> = BTreeMap::new();
    for (face, data) in &brep.faces {
        if data.role != TopologyRole::Production {
            continue;
        }
        if let Some(key) = plane_key(brep, face) {
            by_plane.entry(key).or_default().push(face);
        }
    }
    let mut flaps = Vec::new();
    for faces in by_plane.values() {
        if faces.len() < 2 {
            continue;
        }
        for &candidate in faces {
            let sides = face_sides(brep, candidate);
            if sides.is_empty() {
                continue;
            }
            // Partner: the single Production face this candidate twins.
            let mut partner = None;
            let mut qualifies = true;
            for side in &sides {
                match twin_face(brep, *side) {
                    Some(f) if f == candidate => {}
                    Some(f)
                        if brep
                            .faces
                            .get(f)
                            .is_some_and(|d| d.role == TopologyRole::ConstructionResidual) => {}
                    Some(f) => match partner {
                        None => partner = Some(f),
                        Some(p) if p == f => {}
                        Some(_) => {
                            qualifies = false;
                            break;
                        }
                    },
                    None => {
                        qualifies = false;
                        break;
                    }
                }
            }
            let Some(partner) = partner else {
                continue;
            };
            if !qualifies
                || plane_key(brep, partner) != plane_key(brep, candidate)
                || !faces_share_plane_interior(brep, candidate, partner)
            {
                continue;
            }
            flaps.push(candidate);
        }
    }
    flaps
}

/// Interior-overlap test on the common plane (grid sampling, strict
/// interior with edge clearance) — the same instrument as the
/// duplicate-coverage serving gate.
fn faces_share_plane_interior(brep: &BRep, a: FaceHandle, b: FaceHandle) -> bool {
    let Some(sof) = brep.faces.get(a).and_then(|d| d.surface.as_ref()) else {
        return false;
    };
    let SurfaceGeometry::Plane { origin, normal, .. } = &sof.geometry else {
        return false;
    };
    let Some(n) = normal.normalized() else {
        return false;
    };
    let e1 = if n.x.abs() < 0.9 {
        popkern_core::geometry::Vec3::new(1.0, 0.0, 0.0)
    } else {
        popkern_core::geometry::Vec3::new(0.0, 1.0, 0.0)
    };
    let Some(e1) = (e1 - n * e1.dot(&n)).normalized() else {
        return false;
    };
    let e2 = n.cross(&e1);
    let ring = |face: FaceHandle| -> Vec<(f64, f64)> {
        brep.face_outer_loop_positions(face)
            .iter()
            .map(|p| {
                let r = *p - *origin;
                (r.dot(&e1), r.dot(&e2))
            })
            .collect()
    };
    let ra = ring(a);
    let rb = ring(b);
    if ra.len() < 3 || rb.len() < 3 {
        return false;
    }
    let inside = |poly: &[(f64, f64)], p: (f64, f64)| -> bool {
        let mut winding = false;
        let mut min_edge = f64::MAX;
        for i in 0..poly.len() {
            let (x0, y0) = poly[i];
            let (x1, y1) = poly[(i + 1) % poly.len()];
            if ((y0 > p.1) != (y1 > p.1)) && (p.0 < (x1 - x0) * (p.1 - y0) / (y1 - y0) + x0) {
                winding = !winding;
            }
            let (dx, dy) = (x1 - x0, y1 - y0);
            let len2 = dx * dx + dy * dy;
            let t = if len2 > 0.0 {
                (((p.0 - x0) * dx + (p.1 - y0) * dy) / len2).clamp(0.0, 1.0)
            } else {
                0.0
            };
            let (cx, cy) = (x0 + t * dx, y0 + t * dy);
            min_edge = min_edge.min(((p.0 - cx).powi(2) + (p.1 - cy).powi(2)).sqrt());
        }
        winding && min_edge > 1.0e-6
    };
    let (mut x0, mut y0, mut x1, mut y1) = (f64::MAX, f64::MAX, f64::MIN, f64::MIN);
    for &(x, y) in &ra {
        x0 = x0.min(x);
        y0 = y0.min(y);
        x1 = x1.max(x);
        y1 = y1.max(y);
    }
    for gx in 1..24 {
        for gy in 1..24 {
            let p = (
                x0 + (x1 - x0) * f64::from(gx) / 24.0,
                y0 + (y1 - y0) * f64::from(gy) / 24.0,
            );
            if inside(&ra, p) && inside(&rb, p) {
                return true;
            }
        }
    }
    false
}

/// A sliver: Production PLANE face with >=2 residual-twinned sides
/// whose remaining twins are all non-plane Production faces (the
/// strip far-end removed-material patch).
fn detect_slivers(brep: &BRep) -> Vec<FaceHandle> {
    let mut slivers = Vec::new();
    for (face, data) in &brep.faces {
        if data.role != TopologyRole::Production {
            continue;
        }
        if plane_key(brep, face).is_none() {
            continue;
        }
        let sides = face_sides(brep, face);
        if sides.is_empty() {
            continue;
        }
        let mut residual_sides = 0usize;
        let mut qualifies = true;
        for side in &sides {
            match twin_face(brep, *side) {
                Some(f)
                    if brep
                        .faces
                        .get(f)
                        .is_some_and(|d| d.role == TopologyRole::ConstructionResidual) =>
                {
                    residual_sides += 1;
                }
                Some(f)
                    if brep.faces.get(f).is_some_and(|d| {
                        d.role == TopologyRole::Production
                            && d.surface.as_ref().is_some_and(|s| {
                                !matches!(s.geometry, SurfaceGeometry::Plane { .. })
                            })
                    }) => {}
                _ => {
                    qualifies = false;
                    break;
                }
            }
        }
        if qualifies && residual_sides >= 2 && residual_sides < sides.len() {
            slivers.push(face);
        }
    }
    slivers
}

/// Remove a face, its loops, and its half-edges; surviving twin sides
/// become boundary sides (mirrors the scaffolding-retirement rules).
/// Tolerant face removal for vestigial/residual faces: half-edges
/// are dropped side-by-side, surviving twins become boundary sides,
/// fully-orphaned edges are retired. Shared with the recordless
/// corner-group unfillet (PK-174), which meets the same
/// ghost-ridden residual class this pass was built for.
pub(crate) fn remove_face(brep: &mut BRep, face: FaceHandle) {
    let sides = face_sides(brep, face);
    let Some(face_data) = brep.faces.get(face) else {
        return;
    };
    let mut loops = vec![face_data.outer_loop];
    loops.extend(face_data.inner_loops.iter().copied());
    for side in sides {
        let Some(side_data) = brep.half_edges.remove(side) else {
            continue;
        };
        if brep.half_edges.contains_key(side_data.twin) {
            brep.boundary_edges.insert(side_data.edge);
            if let Some(edge_data) = brep.edges.get_mut(side_data.edge)
                && edge_data.half_edge == side
            {
                edge_data.half_edge = side_data.twin;
            }
        } else if brep.edges.remove(side_data.edge).is_some() {
            brep.boundary_edges.remove(&side_data.edge);
        }
    }
    for loop_handle in loops {
        brep.loops.remove(loop_handle);
    }
    brep.faces.remove(face);
}

fn side_span(brep: &BRep, side: HalfEdgeHandle) -> Option<(Point3, Point3)> {
    let data = brep.half_edges.get(side)?;
    let start = brep.vertices.get(data.origin)?.position;
    let end = brep
        .vertices
        .get(brep.half_edges.get(data.next)?.origin)?
        .position;
    Some((start, end))
}

fn repair_vertex_outgoing(brep: &mut BRep) {
    let mut outgoing: BTreeMap<VertexHandle, HalfEdgeHandle> = BTreeMap::new();
    for (handle, data) in &brep.half_edges {
        outgoing.entry(data.origin).or_insert(handle);
    }
    let vertices: Vec<_> = brep.vertices.keys().collect();
    for vertex in vertices {
        match outgoing.get(&vertex) {
            Some(handle) => {
                let recorded = brep.vertices[vertex].outgoing;
                let recorded_originates_here = brep
                    .half_edges
                    .get(recorded)
                    .is_some_and(|data| data.origin == vertex);
                if !recorded_originates_here {
                    brep.vertices[vertex].outgoing = *handle;
                }
            }
            None => {
                brep.vertices.remove(vertex);
            }
        }
    }
}
