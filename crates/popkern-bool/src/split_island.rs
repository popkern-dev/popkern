//! Island (closed interior chain) pre-pass for boolean face
//! splitting — split from `split.rs` (500-line cap ratchet). Chains
//! closed intersection fragments wholly interior to a face into
//! inner-loop imprints; the fresh-hole containment rule for nested
//! annular imprints lives here.

use std::collections::HashSet;

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};

use crate::error::BoolDiagnostics;
use crate::intersect::IntersectionRecord;

use super::open_chain::reversed_nurbs;
use super::{FragmentMap, outer_loop_uv_area};
use popkern_core::geometry::curve::CurveEval as _;

pub(super) fn island_chain_pre_pass(
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
        if popkern_core::trace::trace_enabled("boolean") {
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
                // REGION containment against PRE-EXISTING holes: an
                // island inside a hole the face brought into the
                // boolean imprints an invalidly imbricated loop (the
                // vent-slot class), while an island inside a hole
                // minted by THIS boolean's own earlier island is a
                // legitimately nested imprint (a ring lip's annular
                // flush contact) and must proceed.
                Some((u, v)) => {
                    crate::split_inner_loop::point_inside_face_region_excluding_fresh_holes(
                        brep,
                        face,
                        popkern_core::geometry::Point3::new(u, v, 0.0),
                        tolerance,
                        &fragments.intersection_edges,
                    )
                }
                None => false,
            }
        });
        if popkern_core::trace::trace_enabled("boolean") {
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
