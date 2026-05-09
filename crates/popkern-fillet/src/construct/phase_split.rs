//! Phase B/C data structures and helpers for single-pass fillet topology.

use popkern_core::geometry::NurbsCurve;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, VertexHandle};

/// A trim to be applied to a face during Phase B.
pub(super) struct FaceTrim {
    /// Index into the edge_geometries array (determines sort order).
    pub edge_idx: usize,
    /// Which side of the fillet (1 or 2).
    pub side: u8,
    /// Anchor vertices (start, end) created during Phase A.
    pub anchors: (VertexHandle, VertexHandle),
    /// Analytic 3D trim curve from `EdgeGeometry::trim_3d_<side>`.
    /// When the fillet is a cylinder and the adjacent face is a
    /// plane, this is the analytic cyl×plane intersection.  Plumbed
    /// through to `multi_chord` for the #136 commercial-grade fix
    /// (so the cylinder face's outer loop never contains chord-line
    /// edges through cylinder interior).  Currently populated but
    /// unused — slice 6 of the #136 plan will switch `multi_chord`
    /// to use this instead of `helpers::line_nurbs`.
    #[allow(dead_code)]
    pub trim_3d: NurbsCurve,
    /// Pcurve from `EdgeGeometry::pcurve_<side>` — the analytic
    /// trim curve projected into the planar adjacent face's UV.
    /// Will be the pcurve argument to `split_face` once slice 6
    /// of the #136 plan lands.
    #[allow(dead_code)]
    pub pcurve: NurbsCurve,
}

/// Result of Phase B face splitting for one trim.
pub(super) struct TrimSplitResult {
    /// The fillet edge index this trim belongs to.
    pub edge_idx: usize,
    /// Which side (1 or 2).
    pub side: u8,
    /// The consumed sub-face (contains the fillet edge).
    pub consumed: FaceHandle,
    /// The trim edge(s) created by split_face.
    pub trim_edges: Vec<EdgeHandle>,
    /// Pcurve on the kept face's surface, computed during Phase B before kef
    /// changes face identity.
    pub pcurve: NurbsCurve,
    /// 3D trim curve between actual anchor vertex positions.
    pub trim_3d: NurbsCurve,
}

/// Find which face in `candidates` has both vertices on its outer loop.
pub(super) fn find_face_containing_vertices(
    brep: &BRep,
    candidates: &[FaceHandle],
    v1: VertexHandle,
    v2: VertexHandle,
) -> Option<FaceHandle> {
    for &fh in candidates {
        let face = match brep.faces.get(fh) {
            Some(f) => f,
            None => continue,
        };
        let start_he = match brep.loops.get(face.outer_loop) {
            Some(l) => l.half_edge,
            None => continue,
        };
        if !brep.half_edges.contains_key(start_he) {
            continue;
        }

        let mut found_v1 = false;
        let mut found_v2 = false;
        let mut cursor = start_he;
        let limit = brep.half_edges.len() + 1;
        for _ in 0..limit {
            let origin = brep.half_edges[cursor].origin;
            if origin == v1 {
                found_v1 = true;
            }
            if origin == v2 {
                found_v2 = true;
            }
            if found_v1 && found_v2 {
                return Some(fh);
            }
            cursor = brep.half_edges[cursor].next;
            if cursor == start_he {
                break;
            }
        }
    }
    None
}

use super::{EdgeGeometry, helpers, multi_chord, pcurve_attach};
use crate::error::BlendFailure;
use popkern_core::tolerance::Tolerance;
use std::collections::HashMap;

pub(super) fn apply_face_trim_pipeline(
    brep: &BRep,
    result: &mut BRep,
    edge_geometries: &[EdgeGeometry],
    tolerance: &Tolerance,
    perp_data: &[Option<super::perpendicular_splits::PerpSplitData>],
    descendants: &mut HashMap<EdgeHandle, Vec<EdgeHandle>>,
    failures: &mut Vec<(EdgeHandle, BlendFailure)>,
) {
    // Collect trim data from Phase A results (BTreeMap for deterministic iteration).
    let mut face_trims: std::collections::BTreeMap<FaceHandle, Vec<FaceTrim>> =
        std::collections::BTreeMap::new();

    for (idx, eg) in edge_geometries.iter().enumerate() {
        let Some(ref pd) = perp_data[idx] else {
            continue;
        };

        face_trims
            .entry(eg.analysis.face_1)
            .or_default()
            .push(FaceTrim {
                edge_idx: idx,
                side: 1,
                anchors: pd.face_1_anchors,
                trim_3d: eg.trim_3d_1.clone(),
                pcurve: eg.pcurve_1.clone(),
            });
        face_trims
            .entry(eg.analysis.face_2)
            .or_default()
            .push(FaceTrim {
                edge_idx: idx,
                side: 2,
                anchors: pd.face_2_anchors,
                trim_3d: eg.trim_3d_2.clone(),
                pcurve: eg.pcurve_2.clone(),
            });
    }

    // Sort trims within each face by edge_idx for deterministic ordering.
    for trims in face_trims.values_mut() {
        trims.sort_by_key(|t| (t.edge_idx, t.side));
    }

    // ── Phase B: split faces by trim curves via `multi_chord`. ───
    let mut trim_results: Vec<TrimSplitResult> = Vec::new();
    for (original_face, trims) in &face_trims {
        let fillet_edges_by_trim: Vec<EdgeHandle> = trims
            .iter()
            .map(|t| edge_geometries[t.edge_idx].analysis.edge)
            .collect();
        let outcome = multi_chord::apply_multi_chord_for_face(
            result,
            *original_face,
            trims,
            &fillet_edges_by_trim,
            descendants,
            tolerance,
        );
        // Fold multi-chord's sub-edge splits (from inserting
        // crossing-point vertices on earlier trims' sub-chords)
        // into the global descendants map. Without this, Phase C
        // can't find the fillet edge on a face whose owning
        // sub-chord edge has been split.
        for (old_edge, edge_a, edge_b) in &outcome.extra_edge_splits {
            for desc_list in descendants.values_mut() {
                if let Some(pos) = desc_list.iter().position(|e| e == old_edge) {
                    desc_list.splice(pos..=pos, [*edge_a, *edge_b]);
                }
            }
        }
        trim_results.extend(outcome.trim_results);
        failures.extend(outcome.failures);
    }

    // ── Re-identify consumed faces from descendants (longest-sub-edge). ──
    for tr in &mut trim_results {
        let original_edge = edge_geometries[tr.edge_idx].analysis.edge;
        let fillet_descendants: &[EdgeHandle] = descendants
            .get(&original_edge)
            .map(|v| v.as_slice())
            .unwrap_or(std::slice::from_ref(&original_edge));

        // The trim's original face is face_1 or face_2 of the
        // fillet edge depending on `tr.side`.
        let eg = &edge_geometries[tr.edge_idx];
        let desired_parent = if tr.side == 1 {
            eg.analysis.face_1
        } else {
            eg.analysis.face_2
        };

        // Find the longest descendant and, of its two adjacent
        // faces, pick the one that's a descendant of
        // `desired_parent`.
        let mut longest: Option<(f64, FaceHandle, FaceHandle)> = None;
        for &d in fillet_descendants {
            let Some(edata) = result.edges.get(d) else {
                continue;
            };
            let he_h = edata.half_edge;
            let Some(he) = result.half_edges.get(he_h) else {
                continue;
            };
            let Some(twin) = result.half_edges.get(he.twin) else {
                continue;
            };
            let p0 = result.vertices[he.origin].position;
            let p1 = result.vertices[twin.origin].position;
            let len = p0.distance_to(&p1);
            let f_a = result.loops[he.loop_].face;
            let f_b = result.loops[twin.loop_].face;
            if longest.map(|(l, _, _)| len > l).unwrap_or(true) {
                longest = Some((len, f_a, f_b));
            }
        }
        if let Some((_, f_a, f_b)) = longest {
            // Prefer whichever of f_a, f_b descends from `desired_parent`.
            let pick = multi_chord::choose_side_face(result, brep, desired_parent, f_a, f_b);
            tr.consumed = pick;
        }
    }

    // ── Phase C: merge consumed face pairs via kef ────────────────
    let mut edge_trims: std::collections::BTreeMap<usize, Vec<&TrimSplitResult>> =
        std::collections::BTreeMap::new();
    for tr in &trim_results {
        edge_trims.entry(tr.edge_idx).or_default().push(tr);
    }
    let mut merged_consumed: std::collections::HashSet<FaceHandle> =
        std::collections::HashSet::new();

    for (idx, trims) in &edge_trims {
        let eg = &edge_geometries[*idx];
        let fillet_edge = eg.analysis.edge;

        let side_1 = trims.iter().find(|t| t.side == 1);
        let side_2 = trims.iter().find(|t| t.side == 2);

        let (Some(s1), Some(s2)) = (side_1, side_2) else {
            // Missing one side — consumed face(s) cleaned up below.
            continue;
        };

        let fillet_descendants: &[EdgeHandle] = descendants
            .get(&fillet_edge)
            .map(|v| v.as_slice())
            .unwrap_or(std::slice::from_ref(&fillet_edge));
        match helpers::merge_consumed_faces_any(
            result,
            s1.consumed,
            s2.consumed,
            fillet_descendants,
            &eg.surface,
        ) {
            Ok(fillet_face) => {
                merged_consumed.insert(s1.consumed);
                merged_consumed.insert(s2.consumed);

                // Use pcurves computed during Phase B (before kef).
                let params = pcurve_attach::FilletPcurveParams {
                    fillet_face,
                    fillet_surface: &eg.surface,
                    trim_edges_1: &s1.trim_edges,
                    trim_edges_2: &s2.trim_edges,
                    eg_pcurve_1: &s1.pcurve,
                    eg_pcurve_2: &s2.pcurve,
                    trim_3d_1: &s1.trim_3d,
                    trim_3d_2: &s2.trim_3d,
                };
                pcurve_attach::attach_fillet_pcurves(result, &params);
            }
            Err(e) => {
                let loc = crate::geometry_compute::edge_midpoint_by_handle(brep, fillet_edge);
                failures.push((
                    fillet_edge,
                    BlendFailure::InvalidTopology {
                        reason: e,
                        location: loc,
                    },
                ));
            }
        }
    }

    // ── Phase C cleanup: undo orphaned consumed faces ─────────────
    // If a fillet edge only had one side split, its consumed face is
    // dangling. Undo the split by kef-ing the trim edge, but only if
    // the resulting merged face stays small (≤ 8 vertices). This
    // prevents creating oversized faces that no surface can fit.
    for tr in &trim_results {
        if merged_consumed.contains(&tr.consumed) {
            continue;
        }
        if !result.faces.contains_key(tr.consumed) {
            continue;
        }
        for &trim_edge in &tr.trim_edges {
            if !result.edges.contains_key(trim_edge) {
                continue;
            }
            // Check vertex count of both faces before merging.
            let he = result.edges[trim_edge].half_edge;
            let he_face = result.loops[result.half_edges[he].loop_].face;
            let twin = result.half_edges[he].twin;
            let twin_face = result.loops[result.half_edges[twin].loop_].face;
            if he_face == twin_face {
                continue;
            }
            let verts_a = result.face_outer_loop_positions(he_face).len();
            let verts_b = result.face_outer_loop_positions(twin_face).len();
            // After kef, the merged face has approximately verts_a + verts_b - 2 vertices.
            if verts_a + verts_b > 10 {
                continue; // Would create an oversized face.
            }
            let _ = result.kef(he);
        }
    }
}
