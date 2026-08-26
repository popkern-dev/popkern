//! Stage 2: Split faces of both solids along intersection curves.
//!
//! Groups intersection records by face, then calls `BRep::split_face` on each
//! affected face to produce face fragments. Tracks provenance (which original
//! face each fragment came from) and which edges were created by splitting.

use std::collections::{HashMap, HashSet};

use popkern_core::{
    tolerance::Tolerance,
    topology::{BRep, EdgeHandle, FaceHandle},
};

use crate::error::BoolDiagnostics;
use crate::intersect::IntersectionRecord;
use popkern_core::geometry::curve::CurveEval as _;

#[path = "split_helpers.rs"]
mod helpers;
#[path = "split_helpers_uv.rs"]
mod helpers_uv;
#[path = "split_open_chain.rs"]
mod open_chain;
use helpers::{find_fragment_for_curve, split_one_face};
use open_chain::{open_chain_pre_pass, reversed_nurbs};

/// Provenance map tracking which original face each fragment came from,
/// plus which edges were created by the splitting process.
#[derive(Debug, Clone)]
pub(crate) struct FragmentMap {
    /// Maps each new fragment face handle back to the original face it was split from.
    pub origin: HashMap<FaceHandle, FaceHandle>,
    /// The set of edges created during splitting (trim curve edges).
    pub intersection_edges: HashSet<EdgeHandle>,
}

/// Result of splitting both solids along intersection curves.
#[derive(Debug)]
pub(crate) struct SplitResult {
    /// The modified solid A (with faces split).
    pub a: BRep,
    /// The modified solid B (with faces split).
    pub b: BRep,
    /// Fragment provenance for solid A.
    pub fragments_a: FragmentMap,
    /// Fragment provenance for solid B.
    pub fragments_b: FragmentMap,
}

/// Apply all intersection curves for one face of solid A or B.
///
/// Processes curves in order. The first curve targets `original_face`.
/// Subsequent curves locate the correct fragment via UV midpoint test
/// and split that fragment, with graceful fallback on failure.
#[allow(clippy::too_many_arguments)]
fn split_face_curves(
    brep: &mut BRep,
    original_face: FaceHandle,
    entries: &[(usize, bool)],
    records: &[IntersectionRecord],
    tolerance: &Tolerance,
    fragments: &mut FragmentMap,
) -> Result<(), BoolDiagnostics> {
    // Island pre-pass: the marching intersector delivers a closed
    // interior loop (e.g. a boss poking through a plate face) as a
    // CHAIN of open fragments. Each fragment alone has zero boundary
    // crossings and is not closed, so the per-curve island path
    // never fires and the hole is silently dropped (the union then
    // fails Euler–Poincaré). Chain the fragments by 3D endpoints;
    // a closed chain whose pieces all lie inside the face becomes an
    // inner loop carrying the TRUE piece curves.
    let consumed =
        island_chain_pre_pass(brep, original_face, entries, records, tolerance, fragments)?;

    // Open-chain pre-pass (the island chainer's open sibling): two
    // trim curves that DEAD-END at a shared interior point — a
    // bounded cutter's bottom and side planes meeting at its crease
    // where the crease pierces this face — can never split
    // one-at-a-time (each endpoint matches no vertex and no boundary
    // edge; the T-junction pre-split and the deferred retry both
    // need an EDGE to exist first). Chain them into one
    // boundary-to-boundary polyline and split with that.
    let chained = open_chain_pre_pass(
        brep,
        original_face,
        entries,
        records,
        &consumed,
        tolerance,
        fragments,
    )?;

    let mut deferred: Vec<(usize, bool)> = Vec::new();
    for (i, &(record_idx, use_uv1)) in entries.iter().enumerate() {
        if consumed.contains(&record_idx) || chained.contains(&record_idx) {
            continue;
        }
        let record = &records[record_idx];
        let curve_uv = if use_uv1 {
            &record.curve.curve_uv1
        } else {
            &record.curve.curve_uv2
        };

        let target_face = if i == 0 {
            original_face
        } else {
            // Multi-curve splitting: locate the fragment containing this
            // curve's UV midpoint and split that fragment.
            match find_fragment_for_curve(brep, original_face, curve_uv, fragments, tolerance) {
                Some(frag) => frag,
                None => continue,
            }
        };

        if !brep.faces.contains_key(target_face) {
            continue;
        }

        // T-junction order sensitivity: a curve that DEAD-ENDS on
        // another trim curve can only split once that other curve has
        // become an edge (the pre-split in split_one_face splits
        // boundary EDGES, not future trims). Defer failures and retry
        // after the rest of this face's curves have landed.
        if split_one_face(
            brep,
            target_face,
            &record.curve.curve_3d,
            curve_uv,
            &record.curve,
            tolerance,
            records,
            fragments,
        )
        .is_err()
        {
            deferred.push((record_idx, use_uv1));
        }
    }

    for pass in 0..2 {
        if deferred.is_empty() {
            break;
        }
        let batch = std::mem::take(&mut deferred);
        let last_pass = pass == 1;
        for (record_idx, use_uv1) in batch {
            let record = &records[record_idx];
            let curve_uv = if use_uv1 {
                &record.curve.curve_uv1
            } else {
                &record.curve.curve_uv2
            };
            let target_face = match find_fragment_for_curve(
                brep,
                original_face,
                curve_uv,
                fragments,
                tolerance,
            ) {
                Some(frag) => frag,
                None => continue,
            };
            if !brep.faces.contains_key(target_face) {
                continue;
            }
            let result = split_one_face(
                brep,
                target_face,
                &record.curve.curve_3d,
                curve_uv,
                &record.curve,
                tolerance,
                records,
                fragments,
            );
            match result {
                Ok(()) => {}
                Err(e) if last_pass => return Err(e),
                Err(_) => deferred.push((record_idx, use_uv1)),
            }
        }
    }
    Ok(())
}

/// Chain open intersection fragments into closed interior loops and
/// add each as an inner loop of `face`. Returns the record indices
/// consumed as islands.
#[allow(clippy::too_many_arguments)]
fn island_chain_pre_pass(
    brep: &mut BRep,
    face: FaceHandle,
    entries: &[(usize, bool)],
    records: &[IntersectionRecord],
    tolerance: &Tolerance,
    fragments: &mut FragmentMap,
) -> Result<HashSet<usize>, BoolDiagnostics> {
    use popkern_core::topology::trim::InnerLoopPiece;
    let mut consumed: HashSet<usize> = HashSet::new();
    if entries.len() < 2 {
        return Ok(consumed);
    }
    // Endpoints of each record's 3D curve.
    let ends: Vec<(
        popkern_core::geometry::Point3,
        popkern_core::geometry::Point3,
    )> = entries
        .iter()
        .map(|&(ri, _)| {
            let c = &records[ri].curve.curve_3d;
            let (t0, t1) = c.parameter_range();
            (c.point_at(t0), c.point_at(t1))
        })
        .collect();
    // Junction gate: marching fragments meet with up to ~1e-3 of
    // slop (measured 3.5e-4..1.2e-3 on a NURBS-cylinder × plane
    // circle), far above the linear tolerance. Scale with the
    // smallest fragment so the gate can never bridge a real gap
    // between distinct loops.
    let min_chord = entries
        .iter()
        .map(|&(ri, _)| {
            let c = &records[ri].curve.curve_3d;
            let (t0, t1) = c.parameter_range();
            c.point_at(t0).distance_to(&c.point_at(t1))
        })
        .fold(f64::INFINITY, f64::min);
    let chain_gate = (0.05 * min_chord).max(tolerance.linear() * 10.0);
    let near = |a: &popkern_core::geometry::Point3, b: &popkern_core::geometry::Point3| {
        a.distance_to(b) < chain_gate
    };
    // Deduplicate: the record list can carry the same fragment twice
    // (both march directions); a duplicate chains with its own
    // reverse into a degenerate zero-area "loop". Compare endpoint
    // pairs + midpoints.
    let mids: Vec<popkern_core::geometry::Point3> = entries
        .iter()
        .map(|&(ri, _)| {
            let c = &records[ri].curve.curve_3d;
            let (t0, t1) = c.parameter_range();
            c.point_at(0.5 * (t0 + t1))
        })
        .collect();
    let mut used = vec![false; entries.len()];
    // Accepted island geometry (sampled) — the marching stage can
    // deliver the SAME intersection loop many times, split at
    // different parameter offsets; only one copy becomes the hole.
    let mut accepted: Vec<Vec<popkern_core::geometry::Point3>> = Vec::new();
    let _ = &mids;
    for seed in 0..entries.len() {
        if used[seed] {
            continue;
        }
        // Grow a chain from the seed: entries are (list index, reversed).
        let mut chain: Vec<(usize, bool)> = vec![(seed, false)];
        let mut chain_used: Vec<bool> = used.clone();
        chain_used[seed] = true;
        let chain_start = ends[seed].0;
        let mut tip = ends[seed].1;
        let closed = loop {
            if near(&tip, &chain_start) && chain.len() >= 2 {
                break true;
            }
            let mut extended = false;
            for k in 0..entries.len() {
                if chain_used[k] {
                    continue;
                }
                if near(&ends[k].0, &tip) {
                    chain.push((k, false));
                    chain_used[k] = true;
                    tip = ends[k].1;
                    extended = true;
                    break;
                }
                if near(&ends[k].1, &tip) {
                    chain.push((k, true));
                    chain_used[k] = true;
                    tip = ends[k].0;
                    extended = true;
                    break;
                }
            }
            if !extended {
                break false;
            }
        };
        if std::env::var_os("PK_BOOL_TRACE").is_some() {
            eprintln!(
                "  [island] face {face:?}: seed {seed} chain len {} closed={closed} tip-gap {:.2e}",
                chain.len(),
                tip.distance_to(&chain_start)
            );
        }
        if !closed {
            continue;
        }
        // Every piece midpoint must lie inside the face (a chain that
        // touches the boundary belongs to the crossing path).
        let inside = chain.iter().all(|&(k, _)| {
            let ri = entries[k].0;
            let c = &records[ri].curve.curve_3d;
            let (t0, t1) = c.parameter_range();
            let mid = c.point_at(0.5 * (t0 + t1));
            let sof = brep.faces[face].surface.as_ref();
            match sof.and_then(|sof| {
                popkern_core::geometry::uv_projection::project_point_to_uv(&sof.geometry, &mid)
            }) {
                Some((u, v)) => crate::split_inner_loop::point_inside_face_uv(
                    brep,
                    face,
                    popkern_core::geometry::Point3::new(u, v, 0.0),
                    tolerance,
                ),
                None => false,
            }
        });
        if std::env::var_os("PK_BOOL_TRACE").is_some() {
            eprintln!("  [island] inside={inside}");
        }
        if !inside {
            continue;
        }
        // Geometric duplicate of an already-accepted island?
        let chain_samples: Vec<popkern_core::geometry::Point3> = chain
            .iter()
            .flat_map(|&(k, _)| {
                let ri = entries[k].0;
                let c = &records[ri].curve.curve_3d;
                let (t0, t1) = c.parameter_range();
                (0..8).map(move |q| c.point_at(t0 + (t1 - t0) * q as f64 / 8.0))
            })
            .collect();
        let duplicate = accepted.iter().any(|island| {
            chain_samples
                .iter()
                .all(|p| island.iter().any(|q| p.distance_to(q) < chain_gate * 4.0))
        });
        if duplicate {
            // Consume the records so the crossing path doesn't churn
            // on them.
            for &(k, _) in &chain {
                consumed.insert(entries[k].0);
            }
            used = chain_used;
            continue;
        }
        // Build oriented pieces.
        let mut pieces: Vec<InnerLoopPiece> = Vec::with_capacity(chain.len());
        for &(k, rev) in &chain {
            let (ri, use_uv1) = entries[k];
            let rec = &records[ri];
            let uv = if use_uv1 {
                &rec.curve.curve_uv1
            } else {
                &rec.curve.curve_uv2
            };
            let (c3, cuv) = if rev {
                (reversed_nurbs(&rec.curve.curve_3d), reversed_nurbs(uv))
            } else {
                (rec.curve.curve_3d.clone(), uv.clone())
            };
            pieces.push(InnerLoopPiece {
                curve_3d: c3,
                curve_uv: cuv,
            });
        }
        // Hole winding must oppose the outer loop in UV: compare
        // shoelace signs and reverse the whole chain if they match.
        let loop_area = |ps: &[InnerLoopPiece]| -> f64 {
            let mut pts = Vec::new();
            for p in ps {
                let (t0, t1) = p.curve_uv.parameter_range();
                for i in 0..8 {
                    let t = t0 + (t1 - t0) * i as f64 / 8.0;
                    pts.push(p.curve_uv.point_at(t));
                }
            }
            let mut a = 0.0;
            for i in 0..pts.len() {
                let j = (i + 1) % pts.len();
                a += pts[i].x * pts[j].y - pts[j].x * pts[i].y;
            }
            0.5 * a
        };
        let outer_area = outer_loop_uv_area(brep, face);
        if loop_area(&pieces).signum() == outer_area.signum() {
            pieces = pieces
                .into_iter()
                .rev()
                .map(|p| InnerLoopPiece {
                    curve_3d: reversed_nurbs(&p.curve_3d),
                    curve_uv: reversed_nurbs(&p.curve_uv),
                })
                .collect();
        }
        // Weld the piece junctions exactly: move each pair of
        // adjacent endpoint control points to their midpoint (the
        // marching slop is far below any fragment's own scale).
        let n_pieces = pieces.len();
        for i in 0..n_pieces {
            let j = (i + 1) % n_pieces;
            let end = *pieces[i]
                .curve_3d
                .control_points
                .last()
                .unwrap_or_else(|| unreachable!("piece with no control points"));
            let start = pieces[j].curve_3d.control_points[0];
            let weld = popkern_core::geometry::Point3::new(
                0.5 * (end.x + start.x),
                0.5 * (end.y + start.y),
                0.5 * (end.z + start.z),
            );
            if let Some(p) = pieces[i].curve_3d.control_points.last_mut() {
                *p = weld;
            }
            pieces[j].curve_3d.control_points[0] = weld;
            // Mirror in UV so the pcurves stay junction-consistent.
            let Some(&end_uv) = pieces[i].curve_uv.control_points.last() else {
                continue;
            };
            let start_uv = pieces[j].curve_uv.control_points[0];
            let weld_uv = popkern_core::geometry::Point3::new(
                0.5 * (end_uv.x + start_uv.x),
                0.5 * (end_uv.y + start_uv.y),
                0.5 * (end_uv.z + start_uv.z),
            );
            if let Some(p) = pieces[i].curve_uv.control_points.last_mut() {
                *p = weld_uv;
            }
            pieces[j].curve_uv.control_points[0] = weld_uv;
        }
        match brep.add_inner_loop_from_pieces(face, &pieces, tolerance) {
            Ok(result) => {
                for &edge in &result.edges {
                    fragments.intersection_edges.insert(edge);
                }
                for &(k, _) in &chain {
                    consumed.insert(entries[k].0);
                }
                used = chain_used;
                accepted.push(chain_samples);
            }
            Err(e) => {
                return Err(BoolDiagnostics::from_error(
                    crate::error::BoolError::SplitFailed {
                        face,
                        reason: format!("island inner loop: {e}"),
                    },
                    crate::error::Stage::Split,
                ));
            }
        }
    }
    Ok(consumed)
}

/// Shoelace area of a face's outer loop in its UV chart.
fn outer_loop_uv_area(brep: &BRep, face: FaceHandle) -> f64 {
    let mut pts = Vec::new();
    if let Some(boundary) = crate::clip::face_uv_boundary(brep, face) {
        for (_he, pcurve) in &boundary {
            let (t0, t1) = pcurve.parameter_range();
            for i in 0..8 {
                let t = t0 + (t1 - t0) * i as f64 / 8.0;
                pts.push(pcurve.point_at(t));
            }
        }
    }
    let mut a = 0.0;
    for i in 0..pts.len() {
        let j = (i + 1) % pts.len();
        a += pts[i].x * pts[j].y - pts[j].x * pts[i].y;
    }
    0.5 * a
}

/// Split the faces of two BRep solids along intersection curves.
///
/// For each intersection record, splits the corresponding face of A using
/// `curve_uv1` (the UV curve on surface A) and the face of B using
/// `curve_uv2` (the UV curve on surface B).
///
/// All intersection curves per face are processed in sequence. After each
/// split, subsequent curves target whichever fragment contains the curve's
/// midpoint. Faces with no intersection curves map to themselves in the
/// fragment map.
///
/// # Errors
///
/// Returns [`BoolDiagnostics`] if any `split_face` call fails.
pub(crate) fn split_along_curves(
    mut a: BRep,
    mut b: BRep,
    records: &[IntersectionRecord],
    tolerance: &Tolerance,
) -> Result<SplitResult, BoolDiagnostics> {
    let mut fragments_a = FragmentMap {
        origin: HashMap::new(),
        intersection_edges: HashSet::new(),
    };
    let mut fragments_b = FragmentMap {
        origin: HashMap::new(),
        intersection_edges: HashSet::new(),
    };

    // Group records by face_a and face_b (collect all curves per face).
    // Deterministic face processing order: the split cascade mutates
    // faces (fragments feed find_fragment_for_curve), so HashMap
    // iteration made whole boolean outcomes run-to-run
    // nondeterministic (measured: the open-chain pre-pass fired in
    // one of three identical runs).
    let mut curves_a: std::collections::BTreeMap<FaceHandle, Vec<(usize, bool)>> =
        std::collections::BTreeMap::new();
    let mut curves_b: std::collections::BTreeMap<FaceHandle, Vec<(usize, bool)>> =
        std::collections::BTreeMap::new();

    for (idx, record) in records.iter().enumerate() {
        curves_a.entry(record.face_a).or_default().push((idx, true));
        curves_b
            .entry(record.face_b)
            .or_default()
            .push((idx, false));
    }

    // Split faces in solid A.
    for (&original_face, entries) in &curves_a {
        if !a.faces.contains_key(original_face) {
            continue;
        }
        split_face_curves(
            &mut a,
            original_face,
            entries,
            records,
            tolerance,
            &mut fragments_a,
        )?;
    }

    // Split faces in solid B.
    for (&original_face, entries) in &curves_b {
        if !b.faces.contains_key(original_face) {
            continue;
        }
        split_face_curves(
            &mut b,
            original_face,
            entries,
            records,
            tolerance,
            &mut fragments_b,
        )?;
    }

    // All unsplit faces in A map to themselves.
    for face_handle in a.faces.keys() {
        fragments_a.origin.entry(face_handle).or_insert(face_handle);
    }

    // All unsplit faces in B map to themselves.
    for face_handle in b.faces.keys() {
        fragments_b.origin.entry(face_handle).or_insert(face_handle);
    }

    Ok(SplitResult {
        a,
        b,
        fragments_a,
        fragments_b,
    })
}

/// Split ONE body's faces along self-intersection curves: each record
/// names two faces of the SAME body; face `face_a` splits with the
/// record's `curve_uv1`, face `face_b` with `curve_uv2`. A face that
/// appears in both roles receives all its curves in one sided pass
/// (per-entry `use_uv1` flags), so fragment targeting stays correct.
///
/// # Errors
///
/// Returns [`BoolDiagnostics`] if any face split fails.
pub(crate) fn self_split_along_curves(
    mut body: BRep,
    records: &[IntersectionRecord],
    tolerance: &Tolerance,
) -> Result<(BRep, FragmentMap), BoolDiagnostics> {
    let mut fragments = FragmentMap {
        origin: HashMap::new(),
        intersection_edges: HashSet::new(),
    };
    let mut per_face: std::collections::BTreeMap<FaceHandle, Vec<(usize, bool)>> =
        std::collections::BTreeMap::new();
    for (idx, record) in records.iter().enumerate() {
        per_face.entry(record.face_a).or_default().push((idx, true));
        per_face
            .entry(record.face_b)
            .or_default()
            .push((idx, false));
    }
    for (&original_face, entries) in &per_face {
        if !body.faces.contains_key(original_face) {
            continue;
        }
        split_face_curves(
            &mut body,
            original_face,
            entries,
            records,
            tolerance,
            &mut fragments,
        )?;
    }
    for face_handle in body.faces.keys() {
        fragments.origin.entry(face_handle).or_insert(face_handle);
    }
    Ok((body, fragments))
}

#[cfg(test)]
#[path = "split_tests.rs"]
mod tests;

#[cfg(test)]
#[path = "split_extra_tests.rs"]
mod extra_tests;

#[cfg(test)]
#[path = "split_multicurve_tests.rs"]
mod multicurve_tests;
