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

#[path = "v_free_shell_authority.rs"]
mod shell_authority;
#[path = "v_free_shell_flap_detection.rs"]
mod shell_flap_detection;
use shell_flap_detection::{detect_flaps, detect_slivers};
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
                || self.scaffolding_faces > 0
                || self.welded_side_pairs > 0
                || self.reversed_faces > 0
                || self.merged_vertices > 0)
    }
}

/// Run the consolidation. Self-guarding on a DEFECT signature,
/// not only the aggregate path's flap signature: the PK-167
/// in-place rewrite authors no duplicate face (`dF = 0`, flaps
/// empty by construction) yet still mints dead sides and glued
/// surfaceless membranes that only retire + weld can clear —
/// flap-gating everything made those stages unreachable on the
/// served path (PK-205/PK-206: 12 membranes, 18 dead sides,
/// half-covered support planes from dca764a4 until this guard).
/// Flap/sliver deletion itself stays flap-gated.
pub(crate) fn consolidate_vestigial_support_shell(
    brep: &mut BRep,
) -> VestigialShellConsolidationReport {
    consolidate_vestigial_support_shell_impl(brep, false)
}

/// PK-214 handover-completion entry: the same consolidation with the
/// winding reconciliation forced past its flap gate. The in-place
/// rewrite on the rational-Gregory route serves a body whose new
/// shell (caps, strips, rewritten support quads) is mirror-wound
/// against the untouched far faces — the exact parallel-dead-pair
/// class the reconciliation exists for — but authors no flap, so the
/// flap gate starves it. Callers are transactional (clone + accept
/// only on a strictly better serving-shape census), which is what
/// protects the fixtures the unconditional defect-signature firing
/// used to regress (pyramid4_uniform cap_outer 5e-10 → 7.4e-2).
pub(crate) fn consolidate_vestigial_support_shell_forcing_winding(
    brep: &mut BRep,
) -> VestigialShellConsolidationReport {
    consolidate_vestigial_support_shell_impl(brep, true)
}

fn consolidate_vestigial_support_shell_impl(
    brep: &mut BRep,
    force_winding: bool,
) -> VestigialShellConsolidationReport {
    let mut report = VestigialShellConsolidationReport::default();
    let flaps = detect_flaps(brep);
    let has_dead_sides = !brep.boundary_edges.is_empty();
    let has_glued_surfaceless_residual = brep
        .faces
        .iter()
        .any(|(_, f)| f.role == TopologyRole::ConstructionResidual && f.surface.is_none());
    if flaps.is_empty() && !has_dead_sides && !has_glued_surfaceless_residual {
        return report;
    }
    let rollback = brep.clone();

    // PK-274 dirty inventory: every mutation below contributes the
    // faces it touches, so the pcurve heal at the end can be scoped
    // to their one-ring instead of the whole model (this impl runs
    // up to 3 fixpoint rounds PER CORNER). Deleted-face handles are
    // harmless — the scoped passes ignore dead handles — but their
    // NEIGHBOURS must be noted BEFORE removal. The pipeline-tail
    // whole-model stale cross-check guards the inventory.
    let mut dirty: std::collections::BTreeSet<FaceHandle> = std::collections::BTreeSet::new();
    for face in brep.validation_region(&flaps) {
        dirty.insert(face);
    }
    for face in &flaps {
        remove_face(brep, *face);
        report.flap_faces += 1;
    }
    if !flaps.is_empty() {
        let slivers = detect_slivers(brep);
        for face in brep.validation_region(&slivers) {
            dirty.insert(face);
        }
        for face in slivers {
            remove_face(brep, face);
            report.sliver_faces += 1;
        }
    }
    let scaffolding =
        crate::residual_scaffolding::retire_unglued_construction_residual_scaffolding(brep);
    report.scaffolding_faces = scaffolding.retired_faces;

    // Winding reconciliation is specific to the aggregate path's
    // mirror-wound shell (its 2-coloring assumes the old↔new
    // parallel-seam anatomy). The in-place rewrite preserves the
    // old shell's winding by construction; running the
    // reconciliation there reverses legitimate faces (measured:
    // pyramid4_uniform cap_outer_g0 5e-10 → 7.4e-2 when it fired
    // under the widened defect guard), so it stays flap-gated.
    let mut seam_vertices_from_winding = Vec::new();
    if !flaps.is_empty() || force_winding {
        let winding = shell_orientation::reconcile_shell_winding(brep);
        report.reversed_faces = winding.reversed_faces;
        report.rebuilt_rings = winding.rebuilt_rings;
        seam_vertices_from_winding = winding.seam_vertices;
        dirty.extend(winding.touched_faces);
    }

    let (welded, mut seam_vertices) =
        weld_dead_coincident_sides(brep, &std::collections::BTreeSet::new());
    report.welded_side_pairs = welded;
    seam_vertices.extend(seam_vertices_from_winding);
    report.merged_vertices =
        merge_coincident_seam_vertices(brep, &seam_vertices, &std::collections::BTreeSet::new())
            .len();
    repair_vertex_outgoing(brep);

    // Heal pcurves on faces that gained welded sides, reconcile
    // periodic branch choices (a reversed or welded strip loop can
    // carry a long-way arc image whose net UV winding is ±period —
    // the full-drum area/tessellation misread), then validate.
    // Scoped to the dirty inventory's one-ring (PK-274): welds and
    // vertex merges touch faces incident to the seam vertices; the
    // deletions, reversals, and rebuilds contributed above. The
    // scaffolding retire is exempt by construction — it retires
    // ghost face/loop records "without touching a single half-edge".
    let seam_set: std::collections::BTreeSet<_> = seam_vertices.iter().copied().collect();
    for data in brep.half_edges.values() {
        if seam_set.contains(&data.origin)
            && let Some(loop_data) = brep.loops.get(data.loop_)
        {
            dirty.insert(loop_data.face);
        }
    }
    let dirty: Vec<FaceHandle> = dirty.into_iter().collect();
    let faces = brep.validation_region(&dirty);
    for face in &faces {
        let _ = brep.recompute_face_pcurves(*face);
    }
    brep.fix_stale_pcurves_among(&faces);
    for face in &faces {
        brep.reconcile_periodic_pcurve_branches(*face);
    }
    shell_authority::retire_stale_face_use_pcurve_authorities(brep);

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

/// PK-214 handover-completion stage 3: sever partnerless dead chord
/// sides and close the affected rings by boundary completion, then
/// weld, merge seam vertices, and heal pcurves. Returns the number
/// of rebuilt rings. Callers are transactional.
pub(crate) fn complete_partnerless_chord_rings(brep: &mut BRep) -> usize {
    let entry_faces: std::collections::BTreeSet<FaceHandle> = brep.faces.keys().collect();
    let (rebuilt, mut seam_vertices) = shell_orientation::rebuild_partnerless_chord_rings(brep);
    if rebuilt == 0 {
        return 0;
    }
    let (_welded, more) = weld_dead_coincident_sides(brep, &std::collections::BTreeSet::new());
    seam_vertices.extend(more);
    merge_coincident_seam_vertices(brep, &seam_vertices, &std::collections::BTreeSet::new());
    repair_vertex_outgoing(brep);
    // PK-274: heal pcurves over the faces this completion touched —
    // the ring-rebuild's new faces plus every face incident to a
    // reported seam vertex, dilated one ring — instead of the whole
    // model per corner (this ran once PER CORNER; measured 13/68
    // profile samples on the 5-valence matrix). The pipeline-tail
    // whole-model stale cross-check guards the derivation.
    let seam_set: std::collections::BTreeSet<_> = seam_vertices.iter().copied().collect();
    let mut dirty: std::collections::BTreeSet<FaceHandle> = brep
        .faces
        .keys()
        .filter(|f| !entry_faces.contains(f))
        .collect();
    for data in brep.half_edges.values() {
        if seam_set.contains(&data.origin)
            && let Some(loop_data) = brep.loops.get(data.loop_)
        {
            dirty.insert(loop_data.face);
        }
    }
    let dirty: Vec<FaceHandle> = dirty.into_iter().collect();
    let region = brep.validation_region(&dirty);
    for face in &region {
        let _ = brep.recompute_face_pcurves(*face);
    }
    brep.fix_stale_pcurves_among(&region);
    for face in &region {
        brep.reconcile_periodic_pcurve_branches(*face);
    }
    shell_authority::retire_stale_face_use_pcurve_authorities(brep);
    rebuilt
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
        let Some(side_data) = brep.remove_half_edge_reanchored(side) else {
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
