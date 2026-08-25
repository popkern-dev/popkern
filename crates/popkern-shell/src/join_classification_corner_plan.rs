//! Corner planning for join classification (PK-325): fan-run
//! partitioning at gap-corner vertices and the retained-run rule.

use std::collections::{BTreeMap, BTreeSet, HashMap};

use super::EdgeJoinClass;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, VertexHandle};

/// Plan one gap corner: verify every incident non-seam edge is a gap
/// edge and choose the retained face run. Returns the retained OUTER
/// face, or the blocking mechanism.
pub(super) fn plan_corner(
    skin: &BRep,
    vh: VertexHandle,
    classes: &BTreeMap<EdgeHandle, EdgeJoinClass>,
    reverse_face: &HashMap<FaceHandle, FaceHandle>,
) -> Result<FaceHandle, String> {
    let fan = vertex_fan(skin, vh)
        .ok_or_else(|| "rim/boundary vertex with incident gap edges".to_string())?;

    // Runs are delimited by gap edges; each run's wedges must belong to
    // one face (seam edges of that face may sit inside the run).
    let m = fan.len();
    let gap_positions: Vec<usize> = (0..m)
        .filter(|&i| classes.get(&fan[i].edge) == Some(&EdgeJoinClass::Gap))
        .collect();
    if gap_positions.is_empty() {
        return Err("vertex reached corner planning without gap edges".into());
    }
    // Any incident non-seam edge that is NOT a gap edge makes this a
    // mixed corner (E1, PK-326).
    for entry in &fan {
        let is_gap = classes.get(&entry.edge) == Some(&EdgeJoinClass::Gap);
        if !is_gap && !entry.is_seam {
            return Err("mixed corner: gap and resolvable/smooth edges meet at this vertex".into());
        }
    }

    // Collect the face of each run (fan interval between consecutive
    // gap edges) plus whether it contains seam edges.
    let k = gap_positions.len();
    if k < 2 {
        return Err("corner with a single gap edge (dangling gap)".into());
    }
    let mut runs: Vec<(FaceHandle, bool)> = Vec::new();
    for run_idx in 0..k {
        let start = gap_positions[run_idx];
        let end = gap_positions[(run_idx + 1) % k];
        // Wedges from `start` (inclusive) to `end` (exclusive): the
        // wedge after fan[i] belongs to face fan[i].wedge_face.
        let mut faces_in_run: BTreeSet<FaceHandle> = BTreeSet::new();
        let mut has_seam = false;
        let mut i = start;
        loop {
            faces_in_run.insert(fan[i].wedge_face);
            let next = (i + 1) % m;
            if next == end {
                break;
            }
            if fan[next].is_seam {
                has_seam = true;
            } else {
                // A non-seam, non-gap edge inside a run was already
                // rejected above; a gap edge here would be a
                // bookkeeping error.
                return Err("corner fan run crosses an unexpected edge".into());
            }
            i = next;
        }
        if faces_in_run.len() != 1 {
            return Err("corner fan run spans multiple faces".into());
        }
        let face = *faces_in_run.iter().next().ok_or("empty corner fan run")?;
        let outer = reverse_face
            .get(&face)
            .copied()
            .ok_or("corner fan run face missing from the skin map")?;
        runs.push((outer, has_seam));
    }

    let seam_runs: Vec<&(FaceHandle, bool)> = runs.iter().filter(|(_, s)| *s).collect();
    match seam_runs.len() {
        0 => runs
            .iter()
            .map(|(f, _)| *f)
            .min()
            .ok_or_else(|| "corner has no fan runs".to_string()),
        1 => Ok(seam_runs[0].0),
        _ => Err("seam edges span multiple fan runs at this corner".into()),
    }
}

/// One fan entry at a vertex: the outgoing half-edge's edge, whether it
/// is a seam (same face both sides), and the face of the wedge AFTER it
/// (between it and the fan successor).
struct FanEntry {
    edge: EdgeHandle,
    is_seam: bool,
    wedge_face: FaceHandle,
}

/// Walk the closed fan at a vertex via `twin(o).next`. Returns `None`
/// when the fan crosses a boundary (missing twin) or is inconsistent.
fn vertex_fan(skin: &BRep, vh: VertexHandle) -> Option<Vec<FanEntry>> {
    let anchor = skin.vertices.get(vh)?.outgoing;
    let first = skin.half_edges.get(anchor)?;
    if first.origin != vh {
        return None;
    }
    let mut entries = Vec::new();
    let mut cursor = anchor;
    let bound = skin.half_edges.len();
    loop {
        let he = skin.half_edges.get(cursor)?;
        let twin = skin.half_edges.get(he.twin)?;
        let my_face = skin.loops.get(he.loop_)?.face;
        let twin_face = skin.loops.get(twin.loop_)?.face;
        entries.push(FanEntry {
            edge: he.edge,
            is_seam: my_face == twin_face,
            wedge_face: twin_face,
        });
        let next = twin.next;
        if next == anchor {
            break;
        }
        if skin.half_edges.get(next).is_none_or(|h| h.origin != vh) {
            return None;
        }
        if entries.len() > bound {
            return None;
        }
        cursor = next;
    }
    Some(entries)
}
