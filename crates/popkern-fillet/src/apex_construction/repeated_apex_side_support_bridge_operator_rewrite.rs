//! PK-167 stage 1b: in-place support-loop rewrites as the
//! constructive-winding primary for face-family bridge regions.
//!
//! Instead of authoring a NEW bridge face as the twin cycle of the
//! source fragment's dead sides (a mirror overlay whose winding the
//! consolidation must repair), each region rewrites the SOURCE
//! fragment's loop in place: its dead sides are spliced out and the
//! geometrically-coincident dead sides on the outer faces are
//! twinned into the gap (the Fix-A weld from the consolidation
//! design, executed constructively). Face identity and old-shell
//! winding are preserved; the duplicate coincident edge records are
//! unified (`dE = -removed` per region, `dF = 0`).
//!
//! Fallback boundary: any region whose dead-side chain does not
//! close on shared vertex handles — or whose rewrite declines —
//! returns `None`, and the caller runs the recorded aggregate
//! bridge path unchanged.

use popkern_core::euler::SupportLoopIntervalRewritePlan;
use popkern_core::topology::{BRep, FaceHandle, HalfEdgeHandle, LoopHandle, VertexHandle};

use super::super::repeated_apex_side_support_bridge::FaceFamilyBridgePrimitivePlan;
use super::SideSupportRegionContinuationResult;

/// Attempt the in-place rewrite for EVERY region of the primitive.
/// All-or-nothing: returns `Some(result)` only when each region's
/// chain closes and each rewrite commits; otherwise `None` and the
/// candidate must be discarded by the caller (the aggregate bridge
/// path re-clones from the source).
pub(super) fn try_in_place_support_rewrites(
    primitive: &FaceFamilyBridgePrimitivePlan,
    candidate: &mut BRep,
    plan_edges: &std::collections::BTreeSet<popkern_core::topology::EdgeHandle>,
    plan_vertices: &std::collections::BTreeSet<VertexHandle>,
) -> Option<SideSupportRegionContinuationResult> {
    let mut faces: Vec<FaceHandle> = Vec::new();
    let mut loops: Vec<LoopHandle> = Vec::new();
    let mut plans: Vec<(SupportLoopIntervalRewritePlan, FaceHandle)> = Vec::new();
    for region in primitive.regions() {
        let [cycle] = region.boundary_cycles() else {
            return None;
        };
        let support = cycle.support_interval();
        if support.len() < 2 {
            return None;
        }
        // Consumed prefix = the fragment's dead sides (the fencepost
        // convention: the last interval entry is a connection point).
        let removed: Vec<HalfEdgeHandle> = support[..support.len() - 1].to_vec();
        let first = candidate.half_edges.get(removed[0])?;
        let loop_ = first.loop_;
        let fragment_face = candidate.loops.get(loop_)?.face;
        let start_vertex = first.origin;
        let end_vertex = {
            let last = *removed.last()?;
            let next = candidate.half_edges.get(last)?.next;
            candidate.half_edges.get(next)?.origin
        };
        let replacement = closing_dead_chain(candidate, &removed, start_vertex, end_vertex)?;
        // Fix-A signature gate: the rewrite is the coincident-
        // duplicate WELD — every replacement side must span the same
        // geometry as its removed counterpart, reversed. Regions
        // whose chain closes through genuinely different boundary
        // (measured: direct-apex pyramid/cube fixtures share this
        // bridge machinery and regressed 5 gates without the gate)
        // keep the aggregate bridge path.
        if !chains_coincide(candidate, &removed, &replacement) {
            return None;
        }
        plans.push((
            SupportLoopIntervalRewritePlan {
                loop_,
                removed,
                replacement,
                pair_removed_edges: true,
            },
            fragment_face,
        ));
        faces.push(fragment_face);
        loops.push(loop_);
    }
    if plans.is_empty() {
        return None;
    }
    for (plan, _) in &plans {
        if let Err(err) = candidate.support_loop_interval_rewrite(plan) {
            if std::env::var_os("PK_SLIR_TRACE").is_some() {
                eprintln!("[slir] in-place rewrite declined: {err}");
            }
            return None;
        }
    }
    // The remaining dead sides around the corner (outside interval
    // and rails) are coincident-duplicate seams too; weld them so
    // downstream cap plans see their non-cap twins (this stage used
    // to run inside the consolidation, which the constructive path
    // legitimately short-circuits — no flap exists to detect).
    let (welded, seam_vertices) =
        super::super::v_free_shell_consolidation::weld_dead_coincident_sides(candidate, plan_edges);
    let merged = super::super::v_free_shell_consolidation::merge_coincident_seam_vertices(
        candidate,
        &seam_vertices,
        plan_vertices,
    );
    // Heal pcurves on the faces whose loops gained welded sides (the
    // weld clears the dropped side's pcurve; the consolidation stage
    // that normally recomputes them is legitimately short-circuited
    // by this constructive path — measured: the V-free acceptance
    // declined with missing_support_pcurve_edges=3 on cube_extreme).
    let healed_faces: Vec<popkern_core::topology::FaceHandle> = candidate.faces.keys().collect();
    for face in &healed_faces {
        let _ = candidate.recompute_face_pcurves(*face);
    }
    candidate.fix_stale_pcurves();
    for face in &healed_faces {
        candidate.reconcile_periodic_pcurve_branches(*face);
    }
    if std::env::var_os("PK_SLIR_TRACE").is_some() {
        eprintln!(
            "[slir] in-place rewrites committed for {} region(s); welded {} pair(s), merged {} \
             vertex(es); boundary edges left: {}",
            plans.len(),
            welded,
            merged.len(),
            candidate.boundary_edges.len()
        );
    }
    Some(
        SideSupportRegionContinuationResult::from_face_family_bridge_payload(faces, loops)
            .with_vertex_remap(merged)
            .with_current_winding_preserved(),
    )
}

/// Pairwise geometric coincidence between the removed interval and
/// the replacement chain: same count, and chain[k] spans the same
/// two positions as removed[k] (reversed direction), within the
/// model's linear tolerance.
fn chains_coincide(
    brep: &BRep,
    removed: &[HalfEdgeHandle],
    replacement: &[HalfEdgeHandle],
) -> bool {
    if removed.len() != replacement.len() {
        return false;
    }
    let tol = brep.tolerance.linear().max(1e-9) * 10.0;
    let span = |h: HalfEdgeHandle| -> Option<(
        popkern_core::geometry::Point3,
        popkern_core::geometry::Point3,
    )> {
        let d = brep.half_edges.get(h)?;
        let a = brep.vertices.get(d.origin)?.position;
        let next = brep.half_edges.get(d.next)?;
        let b = brep.vertices.get(next.origin)?.position;
        Some((a, b))
    };
    removed
        .iter()
        .zip(replacement)
        .all(|(&r, &c)| match (span(r), span(c)) {
            (Some((ra, rb)), Some((ca, cb))) => {
                ra.distance_to(&cb) <= tol && rb.distance_to(&ca) <= tol
            }
            _ => false,
        })
}

/// Walk dead (twinless) sides by vertex identity from `start` back
/// to `end`, excluding the removed interval itself. Ordered per the
/// rewrite contract: `destination(chain[0]) == start`, chaining
/// destination -> origin, ending at `end`.
fn closing_dead_chain(
    brep: &BRep,
    removed: &[HalfEdgeHandle],
    start: VertexHandle,
    end: VertexHandle,
) -> Option<Vec<HalfEdgeHandle>> {
    let dead: Vec<HalfEdgeHandle> = brep
        .half_edges
        .iter()
        .filter(|(h, data)| !brep.half_edges.contains_key(data.twin) && !removed.contains(h))
        .map(|(h, _)| h)
        .collect();
    let dest_of = |h: HalfEdgeHandle| -> Option<VertexHandle> {
        let next = brep.half_edges.get(h)?.next;
        Some(brep.half_edges.get(next)?.origin)
    };
    let mut chain: Vec<HalfEdgeHandle> = Vec::new();
    let mut want = start;
    for _ in 0..16 {
        let &h = dead
            .iter()
            .find(|&&h| dest_of(h) == Some(want) && !chain.contains(&h))?;
        chain.push(h);
        want = brep.half_edges[h].origin;
        if want == end {
            return Some(chain);
        }
    }
    None
}
