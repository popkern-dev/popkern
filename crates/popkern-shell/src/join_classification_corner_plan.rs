//! Corner planning for join classification: fan partitioning at
//! gap-corner vertices. All-gap corners take the PK-325 sphere-window
//! fan peel; mixed corners (gap + resolvable/smooth edges at one
//! vertex, E1) and free-rim boundary corners take the PK-326
//! two-sector peel. The E1 probe measured that mixed corners have an
//! empty-interior polar cone — no sphere patch exists — so the served
//! mechanism is a single `split_vertex` whose connector edge is the
//! transversal intersection curve of the two cut-wedge surfaces (or
//! the canal's end cross-section circle at a sheet boundary).

use std::collections::{BTreeMap, BTreeSet, HashMap};

use super::{CornerKind, CornerPlanEntry, CornerTarget, EdgeJoinClass};
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, VertexHandle};

/// Plan one gap corner. Returns the mechanism and retained-sector
/// target, or the blocking mechanism string (whole-skin downgrade).
pub(super) fn plan_corner(
    skin: &BRep,
    vh: VertexHandle,
    classes: &BTreeMap<EdgeHandle, EdgeJoinClass>,
    reverse_face: &HashMap<FaceHandle, FaceHandle>,
    serve_rim_corners: bool,
) -> Result<CornerPlanEntry, String> {
    let fan =
        vertex_fan(skin, vh).ok_or_else(|| "corner fan is not walkable (broken)".to_string())?;
    let has_boundary = fan.iter().any(|e| e.wedge.is_none());
    if has_boundary && !serve_rim_corners {
        return Err("rim/boundary vertex with incident gap edges".into());
    }

    let m = fan.len();
    let gap_positions: Vec<usize> = (0..m)
        .filter(|&i| classes.get(&fan[i].edge) == Some(&EdgeJoinClass::Gap))
        .collect();
    if gap_positions.is_empty() {
        return Err("vertex reached corner planning without gap edges".into());
    }
    // Any incident Irregular edge blocks the corner outright.
    if fan
        .iter()
        .any(|e| classes.get(&e.edge) == Some(&EdgeJoinClass::Irregular))
    {
        return Err("corner touches an irregular (variable-class) edge".into());
    }

    let all_gap = !has_boundary
        && fan
            .iter()
            .all(|entry| entry.is_seam || classes.get(&entry.edge) == Some(&EdgeJoinClass::Gap));
    if all_gap {
        let retained = plan_all_gap_corner(&fan, &gap_positions, reverse_face)?;
        return Ok(CornerPlanEntry {
            kind: CornerKind::AllGap,
            retained: CornerTarget::Face(retained),
        });
    }
    plan_sectored_corner(&fan, classes, reverse_face, has_boundary)
}

/// PK-325 all-gap planning: runs between consecutive gap edges must
/// each be a single face (seams may sit inside); the retained run
/// holds the seams, else the smallest outer face handle.
fn plan_all_gap_corner(
    fan: &[FanEntry],
    gap_positions: &[usize],
    reverse_face: &HashMap<FaceHandle, FaceHandle>,
) -> Result<FaceHandle, String> {
    let m = fan.len();
    let k = gap_positions.len();
    if k < 2 {
        return Err("corner with a single gap edge (dangling gap)".into());
    }
    let mut runs: Vec<(FaceHandle, bool)> = Vec::new();
    for run_idx in 0..k {
        let start = gap_positions[run_idx];
        let end = gap_positions[(run_idx + 1) % k];
        let mut faces_in_run: BTreeSet<FaceHandle> = BTreeSet::new();
        let mut has_seam = false;
        let mut i = start;
        loop {
            if let Some(f) = fan[i].wedge {
                faces_in_run.insert(f);
            }
            let next = (i + 1) % m;
            if next == end {
                break;
            }
            if fan[next].is_seam {
                has_seam = true;
            } else {
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

/// PK-326 sectored planning (mixed and boundary corners): the cyclic
/// fan is cut at gap edges, at boundary wedges (free sheet rims), and
/// at middle faces between consecutive classified interior edges. The
/// corner is served exactly when the cuts leave two sectors — one
/// `split_vertex` with a transversal (or rim end-circle) connector.
fn plan_sectored_corner(
    fan: &[FanEntry],
    classes: &BTreeMap<EdgeHandle, EdgeJoinClass>,
    reverse_face: &HashMap<FaceHandle, FaceHandle>,
    has_boundary: bool,
) -> Result<CornerPlanEntry, String> {
    let m = fan.len();
    let interior_class = |i: usize| -> Option<EdgeJoinClass> {
        if fan[i].is_seam || fan[i].is_rim {
            return None;
        }
        classes.get(&fan[i].edge).copied()
    };
    // Every non-seam, non-rim edge must carry a served class.
    for entry in fan {
        if entry.is_seam || entry.is_rim {
            continue;
        }
        match classes.get(&entry.edge) {
            Some(EdgeJoinClass::Gap | EdgeJoinClass::Resolvable | EdgeJoinClass::Smooth) => {}
            Some(EdgeJoinClass::Irregular) | None => {
                return Err("mixed corner touches an unclassified or irregular edge".into());
            }
        }
    }

    // Slot model: slot 2i = edge i, slot 2i+1 = wedge i (the wedge
    // after edge i). Cuts: gap edges, boundary wedges, middle wedges.
    let is_resolvable = |i: usize| {
        matches!(
            interior_class(i),
            Some(EdgeJoinClass::Resolvable | EdgeJoinClass::Smooth)
        )
    };
    let mut cut = vec![false; 2 * m];
    for i in 0..m {
        if interior_class(i) == Some(EdgeJoinClass::Gap) {
            cut[2 * i] = true;
        }
        if fan[i].wedge.is_none() {
            cut[2 * i + 1] = true;
        } else if is_resolvable(i) && is_resolvable((i + 1) % m) {
            cut[2 * i + 1] = true; // middle face between two resolvable edges
        }
    }
    let n_cuts = cut.iter().filter(|c| **c).count();
    if n_cuts == 0 {
        return Err("corner has no cuts (bookkeeping error)".into());
    }

    // Sectors: maximal runs of non-cut slots.
    let first_cut = (0..2 * m).find(|&s| cut[s]).ok_or("no cut slot")?;
    let mut sectors: Vec<SectorPlan> = Vec::new();
    let mut current: Option<SectorPlan> = None;
    for step in 1..=2 * m {
        let s = (first_cut + step) % (2 * m);
        if cut[s] {
            if let Some(sector) = current.take() {
                sectors.push(sector);
            }
            continue;
        }
        let sector = current.get_or_insert_with(SectorPlan::default);
        if s.is_multiple_of(2) {
            let i = s / 2;
            if is_resolvable(i) {
                let f1 = wedge_face_before(fan, i).ok_or("resolvable edge beside a boundary")?;
                let f2 = fan[i].wedge.ok_or("resolvable edge beside a boundary")?;
                sector.pairs.push((f1, f2));
            }
            if fan[i].is_seam {
                sector.has_seam = true;
            }
        } else if let Some(f) = fan[s / 2].wedge {
            sector.faces.insert(f);
        }
    }
    if let Some(sector) = current.take() {
        sectors.push(sector);
    }

    if sectors.len() != 2 {
        return Err(format!(
            "mixed corner needs {} sectors (only the 2-sector peel is \
             served; more sectors would need a polar-cone cap)",
            sectors.len()
        ));
    }

    let mut targets = Vec::new();
    for sector in &sectors {
        targets.push(sector.target(reverse_face)?);
    }

    let seam_count = sectors.iter().filter(|s| s.has_seam).count();
    let retained = match seam_count {
        1 => {
            let idx = sectors.iter().position(|s| s.has_seam).unwrap_or(0);
            targets[idx]
        }
        n if n > 1 => return Err("seam edges span multiple sectors at this corner".into()),
        _ => {
            let rank = |t: &CornerTarget| match t {
                CornerTarget::EdgePair(a, b) => (0, *a.min(b), *a.max(b)),
                CornerTarget::Face(f) => (1, *f, *f),
            };
            *targets
                .iter()
                .min_by_key(|t| rank(t))
                .ok_or("corner has no sectors")?
        }
    };
    Ok(CornerPlanEntry {
        kind: if has_boundary {
            CornerKind::BoundaryTwoSector
        } else {
            CornerKind::MixedTwoSector
        },
        retained,
    })
}

#[derive(Default)]
struct SectorPlan {
    faces: BTreeSet<FaceHandle>,
    pairs: Vec<(FaceHandle, FaceHandle)>,
    has_seam: bool,
}

impl SectorPlan {
    fn target(
        &self,
        reverse_face: &HashMap<FaceHandle, FaceHandle>,
    ) -> Result<CornerTarget, String> {
        let outer = |f: FaceHandle| -> Result<FaceHandle, String> {
            reverse_face
                .get(&f)
                .copied()
                .ok_or_else(|| "corner sector face missing from the skin map".to_string())
        };
        match self.pairs.len() {
            0 => {
                if self.faces.len() != 1 {
                    return Err("corner sector spans multiple faces without a \
                                classified edge between them"
                        .into());
                }
                let f = *self.faces.iter().next().ok_or("empty corner sector")?;
                Ok(CornerTarget::Face(outer(f)?))
            }
            1 => {
                let (f1, f2) = self.pairs[0];
                Ok(CornerTarget::EdgePair(outer(f1)?, outer(f2)?))
            }
            n => Err(format!(
                "corner sector holds {n} classified edges (expected at most 1)"
            )),
        }
    }
}

/// The face of the wedge BEFORE fan position `p` (between p−1 and p).
fn wedge_face_before(fan: &[FanEntry], p: usize) -> Option<FaceHandle> {
    let m = fan.len();
    fan[(p + m - 1) % m].wedge
}

/// One fan entry at a vertex: the outgoing half-edge's edge, whether it
/// is a seam (same face both sides) or a free-rim edge (a loopless
/// side), and the face of the wedge AFTER it (`None` = the wedge is
/// the sheet boundary).
struct FanEntry {
    edge: EdgeHandle,
    is_seam: bool,
    is_rim: bool,
    wedge: Option<FaceHandle>,
}

/// Walk the closed fan at a vertex via `twin(o).next`. Loopless
/// half-edges (materialized sheet-rim twins) participate — their
/// wedges read as boundary. Returns `None` only when the fan is
/// structurally broken (missing twins, inconsistent origins).
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
        let my_face = skin.loops.get(he.loop_).map(|l| l.face);
        let twin_face = skin.loops.get(twin.loop_).map(|l| l.face);
        entries.push(FanEntry {
            edge: he.edge,
            is_seam: my_face.is_some() && my_face == twin_face,
            is_rim: my_face.is_none() || twin_face.is_none(),
            wedge: twin_face,
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
