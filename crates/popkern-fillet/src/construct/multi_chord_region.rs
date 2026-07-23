//! PK-141 phase 2/2b — corner separators + consumed-region fallback
//! for entangled hole-loop partitions on pierced curved supports.
//! Design: `docs/architecture/curved-support-apply.md`.
//!
//! When a trim's own splits carry no descendant of its fillet edge
//! (the earlier trims' splits own the fillet edge's wall side), the
//! consumed material is the COLLECTIVE RING between the hole loop
//! and the rails — one connected face, because nothing separates
//! adjacent fillets' consumed material at the corners (measured:
//! five corner-spur crossings spanning all four hole corners).
//!
//! This pass, run after every trim's splits on a CANDIDATE clone:
//!
//! 1. **Corner separators (2b)** — for each pair of rail-end anchors
//!    from DIFFERENT trims that face each other across a hole corner
//!    (mutual-nearest within the corner gate), `mef` the face whose
//!    loop carries both, with the straight chord as the separator.
//!    The separators partition the ring per fillet (and later become
//!    the corner caps' boundary).
//! 2. **Region flood (2)** — per deferred fillet: seed at the
//!    sub-faces holding its descendants; flood without crossing any
//!    trim edge or separator; a strict-subset single-region result
//!    merges (clear-surface/merge/restore, same-surface sub-faces)
//!    into the consumed face.
//!
//! All-or-nothing: the clone commits only when every deferred trim
//! resolves; otherwise the caller's state is untouched and the
//! original failures are reported.

use std::collections::HashSet;

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::{NurbsCurve, Point3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle};

use super::multi_chord_region_separators::cut_corner_separators;
use super::phase_split::{FaceTrim, TrimSplitResult};
use crate::error::BlendFailure;

/// One trim deferred from the per-sub-chord identification.
pub(super) struct DeferredTrim {
    pub edge_idx: usize,
    pub side: u8,
    pub fillet_edge: EdgeHandle,
    /// Longest interior sub-chord's `(length, pcurve, trim_3d)`.
    pub candidate: Option<(f64, NurbsCurve, NurbsCurve)>,
    pub produced_edges: Vec<EdgeHandle>,
}
/// Resolve the deferred trims. Pushes results or failures.
#[allow(clippy::too_many_arguments)]
pub(super) fn process_deferred(
    brep: &mut BRep,
    deferred: Vec<DeferredTrim>,
    trims: &[FaceTrim],
    descendants: &std::collections::HashMap<EdgeHandle, Vec<EdgeHandle>>,
    sub_faces: &[FaceHandle],
    all_trim_edges: &[EdgeHandle],
    trim_results: &mut Vec<TrimSplitResult>,
    failures: &mut Vec<(EdgeHandle, BlendFailure)>,
    tolerance: &Tolerance,
) {
    if deferred.is_empty() {
        return;
    }
    let trace = std::env::var_os("PK_MULTICHORD_TRACE").is_some();
    let mut cand = brep.clone();
    let separators = cut_corner_separators(&mut cand, trims, sub_faces, tolerance);
    if trace {
        eprintln!("[mc-region] separators: {separators:?}");
    }
    let mut boundary: Vec<EdgeHandle> = all_trim_edges.to_vec();
    boundary.extend(separators.iter().copied());
    // Corner-territory closure: the spur edges (hole corner ↔ rail
    // anchor) let the flood tunnel through the corner triangles
    // around the separators. Block every corner-scale edge whose
    // BOTH endpoints sit within the corner gate of the rail-end
    // anchor set (rails themselves are already trim edges; hole
    // edges border non-member faces).
    {
        let mut anchor_pts: Vec<Point3> = Vec::new();
        let mut gates: Vec<f64> = Vec::new();
        for t in trims {
            let (t0, t1) = t.trim_3d.parameter_range();
            let (a, b) = (t.trim_3d.point_at(t0), t.trim_3d.point_at(t1));
            let gate = a.distance_to(&b).max(1e-12) * 0.5;
            anchor_pts.push(a);
            anchor_pts.push(b);
            gates.push(gate);
            gates.push(gate);
        }
        let near_anchor = |p: Point3| -> bool {
            anchor_pts
                .iter()
                .zip(gates.iter())
                .any(|(a, g)| a.distance_to(&p) <= *g)
        };
        for (eh, e) in cand.edges.iter() {
            if boundary.contains(&eh) {
                continue;
            }
            let Some(he) = cand.half_edges.get(e.half_edge) else {
                continue;
            };
            let Some(tw) = cand.half_edges.get(he.twin) else {
                continue;
            };
            let (Some(v0), Some(v1)) = (
                cand.vertices.get(he.origin).map(|v| v.position),
                cand.vertices.get(tw.origin).map(|v| v.position),
            ) else {
                continue;
            };
            if near_anchor(v0) && near_anchor(v1) {
                boundary.push(eh);
            }
        }
    }
    let claimed: Vec<FaceHandle> = trim_results.iter().map(|t| t.consumed).collect();
    let mut staged: Vec<TrimSplitResult> = Vec::new();
    let mut all_ok = true;
    for d in &deferred {
        let fillet_desc: &[EdgeHandle] = descendants
            .get(&d.fillet_edge)
            .map(|v| v.as_slice())
            .unwrap_or(std::slice::from_ref(&d.fillet_edge));
        // Mid-span descendant filter (PK-141 region siding): chain
        // REMNANT segments beyond the anchor cuts (corner- and
        // seam-side pieces) also carry descendant status and seed
        // the flood onto remainder faces (measured: eg0's region =
        // the whole wall including the z=0 rim). A descendant
        // qualifies as a seed only when its midpoint projects to
        // the INTERIOR of the trim span.
        // The FULL analytic trim for this (edge_idx, side) — the
        // deferred candidate's trim_3d is only the longest interior
        // sub-chord and misclassifies span fractions.
        let full_trim = trims
            .iter()
            .find(|t| t.edge_idx == d.edge_idx && t.side == d.side)
            .map(|t| &t.trim_3d);
        let mid_span: Vec<EdgeHandle> = fillet_desc
            .iter()
            .copied()
            .filter(|&e| {
                let Some((frac, dist)) = full_trim.and_then(|trim| {
                    super::multi_chord_region_seeding::descendant_trim_fraction(&cand, e, trim)
                }) else {
                    return false;
                };
                if trace {
                    eprintln!(
                        "[mc-region]   desc {e:?} for eg{} side{}: frac {frac:.3} dist {dist:.4}",
                        d.edge_idx, d.side
                    );
                }
                (0.08..=0.92).contains(&frac)
            })
            .collect();
        let seed_desc: &[EdgeHandle] = if mid_span.is_empty() {
            fillet_desc
        } else {
            &mid_span
        };
        if trace && seed_desc.len() != fillet_desc.len() {
            eprintln!(
                "[mc-region] mid-span seed filter: {} of {} descendants qualify",
                seed_desc.len(),
                fillet_desc.len()
            );
        }
        let resolved = d.candidate.as_ref().and_then(|(_, pcurve, trim_3d)| {
            consumed_region_merge(&mut cand, seed_desc, sub_faces, &boundary, &claimed).map(
                |consumed| TrimSplitResult {
                    edge_idx: d.edge_idx,
                    side: d.side,
                    consumed,
                    trim_edges: d.produced_edges.clone(),
                    pcurve: pcurve.clone(),
                    trim_3d: trim_3d.clone(),
                },
            )
        });
        match resolved {
            Some(r) => staged.push(r),
            None => {
                all_ok = false;
                break;
            }
        }
    }
    if all_ok && cand.validate().is_ok() {
        *brep = cand;
        trim_results.extend(staged);
        return;
    }
    if trace {
        eprintln!("[mc-region] batch declined (all_ok={all_ok}); original failures stand");
    }
    for d in deferred {
        failures.push((
            d.fillet_edge,
            BlendFailure::InvalidTopology {
                reason: "multi-chord trim produced no sub-face containing the fillet edge"
                    .to_string(),
                location: crate::geometry_compute::edge_midpoint_by_handle(brep, d.fillet_edge),
            },
        ));
    }
}
/// Phase 2: flood + merge the consumed region (strict subset,
/// claimed-face exclusion, single-face merge).
fn consumed_region_merge(
    brep: &mut BRep,
    fillet_desc: &[EdgeHandle],
    sub_faces: &[FaceHandle],
    boundary: &[EdgeHandle],
    claimed: &[FaceHandle],
) -> Option<FaceHandle> {
    let trace = std::env::var_os("PK_MULTICHORD_TRACE").is_some();
    // The subdivision may have gained separator pieces; membership =
    // original sub-face set PLUS any face sharing a separator/trim
    // boundary with it is approximated by seeding + flood limits.
    let member: HashSet<FaceHandle> = live_subdivision(brep, sub_faces);
    let mut seeds: Vec<FaceHandle> = Vec::new();
    for &f in &member {
        if crate::construct::merge::face_loop_contains_any_edge(brep, f, fillet_desc)
            .unwrap_or(false)
        {
            if trace {
                for &d in fillet_desc {
                    if let Some(ed) = brep.edges.get(d)
                        && let Some(he) = brep.half_edges.get(ed.half_edge)
                    {
                        let f1 = brep.loops.get(he.loop_).map(|l| l.face);
                        let tw = brep.half_edges.get(he.twin);
                        let f2 = tw.and_then(|t| brep.loops.get(t.loop_)).map(|l| l.face);
                        if f1 == Some(f) || f2 == Some(f) {
                            let a = brep.vertices.get(he.origin).map(|v| v.position);
                            eprintln!("[mc-region]   seed {f:?} matched descendant {d:?} at {a:?}");
                        }
                    }
                }
            }
            seeds.push(f);
        }
    }
    if trace {
        eprintln!("[mc-region] seeds {seeds:?}");
    }
    if seeds.is_empty() {
        return None;
    }
    // Deterministic seed order: `member` is a HashSet, and seed
    // iteration order decided the flood/merge order (measured: a
    // run-to-run defect flicker on the boss chain once the
    // transplant made multiple seeds live).
    seeds.sort();
    if trace {
        for &f in &seeds {
            let positions = brep.face_outer_loop_positions(f);
            let n = positions.len();
            let zmin = positions.iter().map(|p| p.z).fold(f64::MAX, f64::min);
            let zmax = positions.iter().map(|p| p.z).fold(f64::MIN, f64::max);
            eprintln!("[mc-region]   seed {f:?}: {n} loop vertices, z [{zmin:.3}, {zmax:.3}]");
        }
    }
    let mut seen: HashSet<FaceHandle> = seeds.iter().copied().collect();
    let mut queue = seeds;
    while let Some(f) = queue.pop() {
        let Some(fd) = brep.faces.get(f) else {
            continue;
        };
        let mut loops = vec![fd.outer_loop];
        loops.extend(fd.inner_loops.iter().copied());
        for lh in loops {
            let Some(l) = brep.loops.get(lh) else {
                continue;
            };
            let start = l.half_edge;
            let mut cur = start;
            for _ in 0..brep.half_edges.len() + 1 {
                let Some(he) = brep.half_edges.get(cur) else {
                    break;
                };
                if !boundary.contains(&he.edge)
                    && let Some(tw) = brep.half_edges.get(he.twin)
                    && let Some(nl) = brep.loops.get(tw.loop_)
                    && member.contains(&nl.face)
                    && nl.face != f
                    && seen.insert(nl.face)
                {
                    if trace {
                        eprintln!(
                            "[mc-region] cross {:?} via {:?} from {f:?}",
                            nl.face, he.edge
                        );
                    }
                    queue.push(nl.face);
                }
                cur = he.next;
                if cur == start {
                    break;
                }
            }
        }
    }
    if seen.len() >= member.len() || claimed.iter().any(|c| seen.contains(c)) {
        if trace {
            eprintln!(
                "[mc-region] DECLINE: region {}/{} (claimed hit {})",
                seen.len(),
                member.len(),
                claimed.iter().any(|c| seen.contains(c))
            );
        }
        return None;
    }
    let shared_surface = seen
        .iter()
        .find_map(|&f| brep.faces.get(f).and_then(|fd| fd.surface.clone()));
    for &f in &seen {
        if let Some(fd) = brep.faces.get_mut(f) {
            fd.surface = None;
        }
    }
    let mut region = seen;
    loop {
        let mut merge_he = None;
        'outer: for &f in &region {
            let Some(fd) = brep.faces.get(f) else {
                continue;
            };
            let mut loops = vec![fd.outer_loop];
            loops.extend(fd.inner_loops.iter().copied());
            for lh in loops {
                let Some(l) = brep.loops.get(lh) else {
                    continue;
                };
                let start = l.half_edge;
                let mut cur = start;
                for _ in 0..brep.half_edges.len() + 1 {
                    let Some(he) = brep.half_edges.get(cur) else {
                        break;
                    };
                    if !boundary.contains(&he.edge)
                        && let Some(tw) = brep.half_edges.get(he.twin)
                        && let Some(nl) = brep.loops.get(tw.loop_)
                        && nl.face != f
                        && region.contains(&nl.face)
                    {
                        merge_he = Some(cur);
                        break 'outer;
                    }
                    cur = he.next;
                    if cur == start {
                        break;
                    }
                }
            }
        }
        let Some(he) = merge_he else { break };
        let merged = match brep.merge_faces_across_edge(he) {
            Ok(m) => m,
            Err(e) => {
                if trace {
                    eprintln!("[mc-region] merge failed: {e}");
                }
                return None;
            }
        };
        region.remove(&merged.removed_face);
    }
    if region.len() != 1 {
        if trace {
            eprintln!("[mc-region] merge left {}", region.len());
        }
        return None;
    }
    let survivor = region.into_iter().next()?;
    if let Some(fd) = brep.faces.get_mut(survivor) {
        fd.surface = shared_surface;
    }
    if trace {
        eprintln!("[mc-region] consumed -> {survivor:?}");
    }
    Some(survivor)
}
/// The live subdivision: the recorded sub-faces plus separator
/// pieces the mefs minted from them (any face sharing an edge with a
/// recorded sub-face, restricted to faces carrying the SAME surface
/// family is unnecessary — the mef pieces inherit the parent's
/// surface, so one adjacency hop from the recorded set suffices).
fn live_subdivision(brep: &BRep, sub_faces: &[FaceHandle]) -> HashSet<FaceHandle> {
    let mut member: HashSet<FaceHandle> = sub_faces
        .iter()
        .copied()
        .filter(|f| brep.faces.contains_key(*f))
        .collect();
    let snapshot: Vec<FaceHandle> = member.iter().copied().collect();
    for f in snapshot {
        let Some(fd) = brep.faces.get(f) else {
            continue;
        };
        let mut loops = vec![fd.outer_loop];
        loops.extend(fd.inner_loops.iter().copied());
        for lh in loops {
            let Some(l) = brep.loops.get(lh) else {
                continue;
            };
            let start = l.half_edge;
            let mut cur = start;
            for _ in 0..brep.half_edges.len() + 1 {
                let Some(he) = brep.half_edges.get(cur) else {
                    break;
                };
                if let Some(tw) = brep.half_edges.get(he.twin)
                    && let Some(nl) = brep.loops.get(tw.loop_)
                {
                    let same_surface = brep
                        .faces
                        .get(nl.face)
                        .and_then(|a| a.surface.as_ref())
                        .zip(brep.faces.get(f).and_then(|b| b.surface.as_ref()))
                        .map(|(a, b)| {
                            std::mem::discriminant(&a.geometry)
                                == std::mem::discriminant(&b.geometry)
                        })
                        .unwrap_or(false);
                    if same_surface {
                        member.insert(nl.face);
                    }
                }
                cur = he.next;
                if cur == start {
                    break;
                }
            }
        }
    }
    member
}

/// PK-141 uniform resolution: when any trim deferred, the other
/// trims' per-sub-chord identifications came from the same entangled
/// subdivision (their phase-C kef self-pairs) — re-route every trim
/// of the face through the region path.
pub(super) fn reroute_all_when_deferred(
    deferred: &mut Vec<DeferredTrim>,
    trim_results: &mut Vec<TrimSplitResult>,
    trims: &[FaceTrim],
    fillet_edges_by_trim: &[EdgeHandle],
) {
    if deferred.is_empty() {
        return;
    }
    for r in trim_results.drain(..) {
        let fe = trims
            .iter()
            .zip(fillet_edges_by_trim.iter())
            .find(|(t, _)| t.edge_idx == r.edge_idx && t.side == r.side)
            .map(|(_, &fe)| fe);
        let Some(fillet_edge) = fe else { continue };
        deferred.push(DeferredTrim {
            edge_idx: r.edge_idx,
            side: r.side,
            fillet_edge,
            candidate: Some((0.0, r.pcurve, r.trim_3d)),
            produced_edges: r.trim_edges,
        });
    }
}
