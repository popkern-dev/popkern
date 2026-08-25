//! Gap-corner surgery: peel each face run to its offset corner with
//! `split_vertex`, close detours into cap triangles with `mef`, and
//! dissolve interior connectors with `kef`, leaving one k-gon sphere
//! window (design §5 Stage 3a, PK-325).

use std::collections::{BTreeSet, HashMap};

use crate::arc_join::{author_corner_arc, cap_orientation, offset_corner_position};
use popkern_core::euler::SplitVertexResult;
use popkern_core::geometry::Point3;
use popkern_core::geometry::surface::{SurfaceGeometry, SurfaceOnFace};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle, VertexHandle};

struct CornerFanEntry {
    out: HalfEdgeHandle,
    wedge_face: FaceHandle,
    is_seam: bool,
}

struct CornerRun {
    first: HalfEdgeHandle,
    last: HalfEdgeHandle,
    face: FaceHandle,
    has_seam: bool,
}

/// Construct one gap corner. Returns the cap face when one was created
/// (k ≥ 3), `None` for k = 2 (arc only), or the blocking mechanism.
#[allow(clippy::too_many_arguments)]
pub(crate) fn construct_corner(
    skin: &mut BRep,
    v: VertexHandle,
    strip_faces: &BTreeSet<FaceHandle>,
    offset_surfaces: &HashMap<FaceHandle, SurfaceGeometry>,
    reverse_vertex: &HashMap<VertexHandle, VertexHandle>,
    reverse_face: &HashMap<FaceHandle, FaceHandle>,
    outer: &BRep,
    tolerance: &Tolerance,
) -> Result<Option<FaceHandle>, String> {
    let p_orig = reverse_vertex
        .get(&v)
        .and_then(|ov| outer.vertices.get(*ov).map(|vd| vd.position))
        .unwrap_or_else(|| skin.vertices[v].position);

    let fan = corner_fan(skin, v).ok_or("corner fan is not walkable (boundary or broken)")?;
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
        let my_face = skin.loops.get(he.loop_)?.face;
        let wedge_face = skin.loops.get(twin.loop_)?.face;
        entries.push(CornerFanEntry {
            out: cursor,
            wedge_face,
            is_seam: my_face == wedge_face,
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
        .map(|e| strip_faces.contains(&e.wedge_face))
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
                    .map(|&j| fan[j].wedge_face)
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
