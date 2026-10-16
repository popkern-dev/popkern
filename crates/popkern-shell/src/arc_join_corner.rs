//! Gap-corner surgery. All-gap corners: peel each face run to its
//! offset corner with `split_vertex`, close detours into cap triangles
//! with `mef`, and dissolve interior connectors with `kef`, leaving one
//! k-gon sphere window (design §5 Stage 3a, PK-325). Mixed corners
//! (gap + resolvable/smooth edges at one vertex, E1): the polar cone
//! has empty interior — no sphere patch exists — so the fan peels into
//! exactly two sectors with one `split_vertex`, and the connector edge
//! is the transversal intersection curve of the two cut-wedge surfaces
//! (PK-326).

use std::collections::{BTreeMap, BTreeSet, HashMap};

use crate::arc_join::{author_corner_arc, cap_orientation, offset_corner_position};
use crate::arc_join_corner_mixed::construct_mixed_corner;
use crate::join_classification::EdgeJoinClass;
use popkern_core::euler::SplitVertexResult;
use popkern_core::geometry::Point3;
use popkern_core::geometry::surface::{SurfaceGeometry, SurfaceOnFace};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle, VertexHandle};

pub(crate) struct CornerFanEntry {
    pub(crate) out: HalfEdgeHandle,
    pub(crate) edge: EdgeHandle,
    /// The face of the wedge after this entry; `None` = sheet boundary.
    pub(crate) wedge: Option<FaceHandle>,
    pub(crate) is_seam: bool,
}

struct CornerRun {
    first: HalfEdgeHandle,
    last: HalfEdgeHandle,
    face: FaceHandle,
    has_seam: bool,
}

/// Construct one gap corner. Returns the cap face when one was created
/// (all-gap k ≥ 3), `None` for all-gap k = 2 and mixed corners (arc
/// or transversal connector only), or the blocking mechanism.
#[allow(clippy::too_many_arguments)]
pub(crate) fn construct_corner(
    skin: &mut BRep,
    v: VertexHandle,
    strip_faces: &BTreeSet<FaceHandle>,
    classes: &BTreeMap<EdgeHandle, EdgeJoinClass>,
    offset_surfaces: &HashMap<FaceHandle, SurfaceGeometry>,
    reverse_vertex: &HashMap<VertexHandle, VertexHandle>,
    reverse_face: &HashMap<FaceHandle, FaceHandle>,
    reverse_edge: &HashMap<EdgeHandle, EdgeHandle>,
    outer: &BRep,
    tolerance: &Tolerance,
    retrim: &mut BTreeSet<EdgeHandle>,
    rim_splits: &mut Vec<(VertexHandle, VertexHandle)>,
) -> Result<Option<FaceHandle>, String> {
    let p_orig = reverse_vertex
        .get(&v)
        .and_then(|ov| outer.vertices.get(*ov).map(|vd| vd.position))
        .unwrap_or_else(|| skin.vertices[v].position);

    let fan = corner_fan(skin, v).ok_or("corner fan is not walkable (boundary or broken)")?;
    // Mixed / boundary corner (E1, rim corners): classified
    // resolvable/smooth edges, boundary wedges, or retained
    // removed-face wedges (midplane openings — no offset surface)
    // survive in the fan (gap edges were replaced by strips in
    // phase 1). Strips themselves have no outer mapping and never
    // trip the unclassified test.
    let unclassified_wedge = |wedge: Option<FaceHandle>| match wedge {
        None => true,
        Some(f) if strip_faces.contains(&f) => false,
        Some(f) => reverse_face
            .get(&f)
            .is_none_or(|of| !offset_surfaces.contains_key(of)),
    };
    if fan.iter().any(|e| {
        unclassified_wedge(e.wedge)
            || matches!(
                classes.get(&e.edge),
                Some(EdgeJoinClass::Resolvable | EdgeJoinClass::Smooth)
            )
    }) {
        construct_mixed_corner(
            skin,
            v,
            &fan,
            strip_faces,
            classes,
            offset_surfaces,
            reverse_face,
            reverse_edge,
            outer,
            p_orig,
            tolerance,
            retrim,
            rim_splits,
        )?;
        return Ok(None);
    }
    let runs = fan_runs(&fan, strip_faces)?;
    if runs.len() < 2 {
        return Err("corner has fewer than two face runs between strips".into());
    }

    // Offset corner position and ball radius per run.
    let mut targets: Vec<(Point3, f64, FaceHandle)> = Vec::with_capacity(runs.len());
    for run in &runs {
        let outer_face = *reverse_face
            .get(&run.face)
            .ok_or("corner run face missing from the skin map")?;
        let (q, r) = offset_corner_position(offset_surfaces, outer_face, &p_orig)
            .ok_or("corner projection onto an offset surface failed")?;
        targets.push((q, r, outer_face));
    }
    let mean_r = targets.iter().map(|(_, r, _)| *r).sum::<f64>() / targets.len() as f64;
    let r_gate = (mean_r * 1.0e-3).max(tolerance.linear() * 100.0);
    if targets.iter().any(|(_, r, _)| (r - mean_r).abs() > r_gate) {
        return Err(format!(
            "corner ball radius varies across faces beyond {r_gate:.3e}"
        ));
    }

    // Retained run: the one holding seam edges, else the smallest outer
    // face handle (deterministic).
    let seam_count = runs.iter().filter(|r| r.has_seam).count();
    let retained_idx = match seam_count {
        0 => targets
            .iter()
            .enumerate()
            .min_by_key(|(_, (_, _, of))| *of)
            .map(|(i, _)| i)
            .unwrap_or(0),
        1 => runs.iter().position(|r| r.has_seam).unwrap_or(0),
        // CLOSED-RIM anchor (PK-344): two runs, each carrying its own
        // face's seam (the plan admits exactly the both-halves-of-one-
        // closed-edge class). Either run may keep the original vertex.
        2 if runs.len() == 2 => targets
            .iter()
            .enumerate()
            .min_by_key(|(_, (_, _, of))| *of)
            .map(|(i, _)| i)
            .unwrap_or(0),
        _ => return Err("seam edges span multiple fan runs at this corner".into()),
    };

    // Split the non-retained runs in fan order starting after the
    // retained run; the retained run keeps the original vertex.
    let k = runs.len();
    let order: Vec<usize> = (1..k).map(|i| (retained_idx + i) % k).collect();
    let mut results: Vec<SplitVertexResult> = Vec::with_capacity(k - 1);
    for &i in &order {
        let res = skin
            .split_vertex(runs[i].first, runs[i].last, targets[i].0)
            .map_err(|e| format!("split_vertex failed: {e}"))?;
        results.push(res);
    }
    skin.vertices[v].position = targets[retained_idx].0;

    // Seams ending at this corner declined their recompute while the
    // anchor carried an interim position (PK-344); now that every run
    // vertex sits at its true offset corner, re-author the stale ones
    // (endpoint-disagreement gated — consistent seams untouched).
    for entry in &fan {
        if entry.is_seam {
            crate::edge_reauthor::reauthor_stale_same_surface_edge(skin, entry.edge, tolerance)
                .map_err(|e| format!("seam re-author failed: {e}"))?;
        }
    }

    if k == 2 {
        author_corner_arc(skin, results[0].edge, p_orig, mean_r, tolerance)?;
        return Ok(None);
    }

    // Cap triangles between consecutive splits, then dissolve the
    // interior connectors.
    let mut arc_edges: Vec<EdgeHandle> = vec![results[0].edge, results[k - 2].edge];
    let mut cap_face: Option<FaceHandle> = None;
    for j in 0..k - 2 {
        let he1 = results[j].half_edge_twin;
        let he2 = skin.half_edges[results[j + 1].half_edge].next;
        let mef = skin
            .mef(he1, he2)
            .map_err(|e| format!("cap mef failed: {e}"))?;
        arc_edges.push(mef.edge);
        if cap_face.is_none() {
            cap_face = Some(mef.face);
        }
    }
    for res in results.iter().take(k - 2).skip(1) {
        skin.kef(res.half_edge)
            .map_err(|e| format!("cap connector kef failed: {e}"))?;
    }
    let cap_face = cap_face.ok_or("corner produced no cap face")?;

    for edge in arc_edges {
        author_corner_arc(skin, edge, p_orig, mean_r, tolerance)?;
    }

    let sphere = SurfaceGeometry::Sphere {
        center: p_orig,
        radius: mean_r,
        axis: popkern_core::geometry::Vec3::Z,
    };
    let orientation = cap_orientation(skin, cap_face, &sphere).unwrap_or(true);
    skin.faces[cap_face].surface = Some(SurfaceOnFace {
        geometry: sphere,
        orientation,
    });

    Ok(Some(cap_face))
}

/// Walk the closed fan at `v` via `twin(o).next`.
fn corner_fan(skin: &BRep, v: VertexHandle) -> Option<Vec<CornerFanEntry>> {
    let anchor = skin.vertices.get(v)?.outgoing;
    if skin.half_edges.get(anchor)?.origin != v {
        return None;
    }
    let mut entries = Vec::new();
    let mut cursor = anchor;
    let bound = skin.half_edges.len();
    loop {
        let he = skin.half_edges.get(cursor)?;
        let twin = skin.half_edges.get(he.twin)?;
        let my_face = skin.loops.get(he.loop_).map(|l| l.face);
        let wedge = skin.loops.get(twin.loop_).map(|l| l.face);
        entries.push(CornerFanEntry {
            out: cursor,
            edge: he.edge,
            wedge,
            is_seam: my_face.is_some() && my_face == wedge,
        });
        let next = twin.next;
        if next == anchor {
            break;
        }
        if skin.half_edges.get(next).is_none_or(|h| h.origin != v) || entries.len() > bound {
            return None;
        }
        cursor = next;
    }
    Some(entries)
}

/// Partition the cyclic fan into face runs delimited by strip wedges.
/// Run boundaries are the entries whose after-wedge is a strip face.
fn fan_runs(
    fan: &[CornerFanEntry],
    strip_faces: &BTreeSet<FaceHandle>,
) -> Result<Vec<CornerRun>, String> {
    let m = fan.len();
    let is_boundary: Vec<bool> = fan
        .iter()
        .map(|e| e.wedge.is_some_and(|f| strip_faces.contains(&f)))
        .collect();
    let Some(first_boundary) = (0..m).find(|&i| is_boundary[i]) else {
        return Err("corner fan contains no strip wedges".into());
    };
    let mut runs = Vec::new();
    let mut current: Vec<usize> = Vec::new();
    for step in 1..=m {
        let i = (first_boundary + step) % m;
        current.push(i);
        if is_boundary[i] {
            // Close the run.
            let face = if current.len() >= 2 {
                let interior: BTreeSet<FaceHandle> = current[..current.len() - 1]
                    .iter()
                    .filter_map(|&j| fan[j].wedge)
                    .collect();
                if interior.len() != 1 {
                    return Err("corner fan run spans multiple faces".into());
                }
                *interior.iter().next().ok_or("empty corner fan run")?
            } else {
                return Err("degenerate single-rail corner run".into());
            };
            let has_seam = current[..current.len() - 1].iter().any(|&j| fan[j].is_seam);
            runs.push(CornerRun {
                first: fan[current[0]].out,
                last: fan[current[current.len() - 1]].out,
                face,
                has_seam,
            });
            current = Vec::new();
        }
    }
    if !current.is_empty() {
        return Err("corner fan runs did not close cyclically".into());
    }
    Ok(runs)
}
