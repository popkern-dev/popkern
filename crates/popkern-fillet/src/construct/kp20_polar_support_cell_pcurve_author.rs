//! Pcurve loop optimizer for KP20 support-cell topology.

use popkern_core::error::CoreError;
use popkern_core::geometry::gregory::GregoryPatch;
use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::topology::{BRep, FaceHandle, HalfEdgeHandle};

use crate::construct::kp20_polar_expanded_support_binding::Kp20ExpandedSupportSurfaceSlot;

#[derive(Debug, Clone)]
struct SlotUvCandidate {
    uvs: Vec<(f64, f64)>,
    vertex_residual: f64,
}

#[derive(Debug, Clone, Copy)]
struct SlotAdjacency {
    slot_a: usize,
    edge_a: usize,
    slot_b: usize,
    edge_b: usize,
}

#[derive(Debug, Clone, Copy)]
struct LoopSelectionScore {
    worst: f64,
    total: f64,
}

impl LoopSelectionScore {
    fn better_than(self, other: Self) -> bool {
        self.worst + 1.0e-12 < other.worst
            || ((self.worst - other.worst).abs() <= 1.0e-12 && self.total + 1.0e-12 < other.total)
    }
}

pub(crate) fn author_optimized_support_cell_pcurves(
    brep: &mut BRep,
    slots: &[Kp20ExpandedSupportSurfaceSlot],
    support_surfaces: &[GregoryPatch],
    source_map: &[usize],
) -> Result<(), CoreError> {
    if slots.len() != source_map.len() {
        return Err(CoreError::PreconditionFailed(format!(
            "kp20 support-cell optimized pcurves expected {} source refs, got {}",
            slots.len(),
            source_map.len()
        )));
    }
    let adjacencies = support_cell_adjacencies(brep, slots)?;
    let candidates = slot_uv_candidates(brep, slots, support_surfaces, source_map)?;
    let selected =
        select_candidate_indices(&candidates, support_surfaces, source_map, &adjacencies);
    for (slot_index, slot) in slots.iter().enumerate() {
        let half_edges = face_outer_half_edges(brep, slot.face)?;
        let uvs = &candidates[slot_index][selected[slot_index]].uvs;
        for (index, &half_edge) in half_edges.iter().enumerate() {
            let start = uvs[index];
            let end = uvs[(index + 1) % uvs.len()];
            brep.half_edges[half_edge].pcurve = Some(Pcurve::linear(start, end, 0.0, 1.0));
        }
    }
    brep.validate().map_err(|error| {
        CoreError::PreconditionFailed(format!(
            "kp20 support-cell optimized pcurve authoring produced invalid BRep: {error}"
        ))
    })
}

fn select_candidate_indices(
    candidates: &[Vec<SlotUvCandidate>],
    support_surfaces: &[GregoryPatch],
    source_map: &[usize],
    adjacencies: &[SlotAdjacency],
) -> Vec<usize> {
    let mut selected = candidates
        .iter()
        .map(|slot_candidates| {
            slot_candidates
                .iter()
                .enumerate()
                .min_by(|(_, a), (_, b)| a.vertex_residual.total_cmp(&b.vertex_residual))
                .map(|(index, _)| index)
                .unwrap_or(0)
        })
        .collect::<Vec<_>>();
    let mut best_score = selected_loop_score(
        &selected,
        candidates,
        support_surfaces,
        source_map,
        adjacencies,
    );
    for _ in 0..3 {
        let mut improved = false;
        for slot in 0..candidates.len() {
            let current = selected[slot];
            let mut slot_best = current;
            let mut slot_score = best_score;
            for candidate in 0..candidates[slot].len() {
                selected[slot] = candidate;
                let score = selected_loop_score(
                    &selected,
                    candidates,
                    support_surfaces,
                    source_map,
                    adjacencies,
                );
                if score.better_than(slot_score) {
                    slot_score = score;
                    slot_best = candidate;
                }
            }
            selected[slot] = slot_best;
            if slot_score.better_than(best_score) {
                best_score = slot_score;
                improved = true;
            }
        }
        for adjacency in adjacencies {
            let mut pair_best = (selected[adjacency.slot_a], selected[adjacency.slot_b]);
            let mut pair_score = best_score;
            for candidate_a in 0..candidates[adjacency.slot_a].len() {
                selected[adjacency.slot_a] = candidate_a;
                for candidate_b in 0..candidates[adjacency.slot_b].len() {
                    selected[adjacency.slot_b] = candidate_b;
                    let score = selected_loop_score(
                        &selected,
                        candidates,
                        support_surfaces,
                        source_map,
                        adjacencies,
                    );
                    if score.better_than(pair_score) {
                        pair_score = score;
                        pair_best = (candidate_a, candidate_b);
                    }
                }
            }
            selected[adjacency.slot_a] = pair_best.0;
            selected[adjacency.slot_b] = pair_best.1;
            if pair_score.better_than(best_score) {
                best_score = pair_score;
                improved = true;
            }
        }
        if !improved {
            break;
        }
    }
    selected
}

pub(crate) fn author_assigned_support_cell_pcurves(
    brep: &mut BRep,
    faces: &[FaceHandle],
) -> Result<(), CoreError> {
    let mut slots = Vec::with_capacity(faces.len());
    let mut surfaces = Vec::with_capacity(faces.len());
    for &face in faces {
        let Ok(half_edges) = face_outer_half_edges(brep, face) else {
            continue;
        };
        let boundary_edge_count = half_edges.len();
        if boundary_edge_count < 4 {
            continue;
        }
        let Some(SurfaceGeometry::Gregory(surface)) = brep
            .faces
            .get(face)
            .and_then(|face_data| face_data.surface.as_ref())
            .map(|surface_on_face| &surface_on_face.geometry)
        else {
            continue;
        };
        slots.push(Kp20ExpandedSupportSurfaceSlot {
            face,
            boundary_edge_count,
            assignment: if boundary_edge_count == 4 {
                crate::construct::kp20_polar_expanded_support_binding::Kp20ExpandedSupportSurfaceAssignment::FullDomainRectangularPatch
            } else {
                crate::construct::kp20_polar_expanded_support_binding::Kp20ExpandedSupportSurfaceAssignment::TrimmedMacroFace
            },
        });
        surfaces.push(surface.clone());
    }
    let source_map = (0..slots.len()).collect::<Vec<_>>();
    author_optimized_support_cell_pcurves(brep, &slots, &surfaces, &source_map)
}

fn support_cell_adjacencies(
    brep: &BRep,
    slots: &[Kp20ExpandedSupportSurfaceSlot],
) -> Result<Vec<SlotAdjacency>, CoreError> {
    let mut out = Vec::new();
    for (slot_a, slot) in slots.iter().enumerate() {
        let half_edges = face_outer_half_edges(brep, slot.face)?;
        for (edge_a, &half_edge) in half_edges.iter().enumerate() {
            let twin = brep.half_edges[half_edge].twin;
            let Some(face_b) = brep
                .loops
                .get(brep.half_edges[twin].loop_)
                .map(|loop_| loop_.face)
            else {
                continue;
            };
            let Some(slot_b) = slots.iter().position(|slot| slot.face == face_b) else {
                continue;
            };
            if slot_b <= slot_a {
                continue;
            }
            let edge_b = edge_index_for_half_edge(brep, slots[slot_b].face, twin)?;
            out.push(SlotAdjacency {
                slot_a,
                edge_a,
                slot_b,
                edge_b,
            });
        }
    }
    Ok(out)
}

fn edge_index_for_half_edge(
    brep: &BRep,
    face: popkern_core::topology::FaceHandle,
    target: HalfEdgeHandle,
) -> Result<usize, CoreError> {
    let half_edges = face_outer_half_edges(brep, face)?;
    half_edges
        .iter()
        .position(|&half_edge| half_edge == target)
        .ok_or_else(|| {
            CoreError::PreconditionFailed(
                "kp20 support-cell optimizer could not locate twin half-edge".into(),
            )
        })
}

fn slot_uv_candidates(
    brep: &BRep,
    slots: &[Kp20ExpandedSupportSurfaceSlot],
    surfaces: &[GregoryPatch],
    source_map: &[usize],
) -> Result<Vec<Vec<SlotUvCandidate>>, CoreError> {
    let mut out = Vec::with_capacity(slots.len());
    for (slot_index, slot) in slots.iter().enumerate() {
        let half_edges = face_outer_half_edges(brep, slot.face)?;
        let source = source_map[slot_index];
        let candidates = unit_domain_loop_candidates(half_edges.len())
            .into_iter()
            .map(|uvs| {
                let vertex_residual = half_edges
                    .iter()
                    .zip(uvs.iter())
                    .map(|(&half_edge, &(u, v))| {
                        let position = brep.vertices[brep.half_edges[half_edge].origin].position;
                        position.distance_to(&surfaces[source].point_at(u, v))
                    })
                    .sum();
                SlotUvCandidate {
                    uvs,
                    vertex_residual,
                }
            })
            .collect();
        out.push(candidates);
    }
    Ok(out)
}

fn selected_loop_score(
    selected: &[usize],
    candidates: &[Vec<SlotUvCandidate>],
    surfaces: &[GregoryPatch],
    source_map: &[usize],
    adjacencies: &[SlotAdjacency],
) -> LoopSelectionScore {
    let mut worst = 0.0_f64;
    let mut total = 0.0_f64;
    for adjacency in adjacencies {
        let source_a = source_map[adjacency.slot_a];
        let source_b = source_map[adjacency.slot_b];
        let uvs_a = &candidates[adjacency.slot_a][selected[adjacency.slot_a]].uvs;
        let uvs_b = &candidates[adjacency.slot_b][selected[adjacency.slot_b]].uvs;
        let gap = shared_edge_gap(
            &surfaces[source_a],
            uvs_a,
            adjacency.edge_a,
            &surfaces[source_b],
            uvs_b,
            adjacency.edge_b,
        );
        worst = worst.max(gap);
        total += gap;
    }
    LoopSelectionScore { worst, total }
}

fn shared_edge_gap(
    surface_a: &GregoryPatch,
    uvs_a: &[(f64, f64)],
    edge_a: usize,
    surface_b: &GregoryPatch,
    uvs_b: &[(f64, f64)],
    edge_b: usize,
) -> f64 {
    if uvs_a.is_empty() || uvs_b.is_empty() {
        return f64::INFINITY;
    }
    let a0 = uvs_a[edge_a % uvs_a.len()];
    let a1 = uvs_a[(edge_a + 1) % uvs_a.len()];
    let b0 = uvs_b[edge_b % uvs_b.len()];
    let b1 = uvs_b[(edge_b + 1) % uvs_b.len()];
    let mut worst = 0.0_f64;
    for index in 0..9 {
        let t = index as f64 / 8.0;
        let ua = lerp_uv(a0, a1, t);
        let ub = lerp_uv(b1, b0, t);
        let pa = surface_a.point_at(ua.0, ua.1);
        let pb = surface_b.point_at(ub.0, ub.1);
        worst = worst.max(pa.distance_to(&pb));
    }
    worst
}

fn lerp_uv(a: (f64, f64), b: (f64, f64), t: f64) -> (f64, f64) {
    (a.0 + (b.0 - a.0) * t, a.1 + (b.1 - a.1) * t)
}

fn unit_domain_loop_candidates(edge_count: usize) -> Vec<Vec<(f64, f64)>> {
    let corners = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)];
    let bases = if edge_count == 4 {
        vec![corners.to_vec()]
    } else if edge_count == 5 {
        (0..4)
            .map(|split| {
                let mut out = Vec::with_capacity(5);
                for index in 0..4 {
                    out.push(corners[index]);
                    if index == split {
                        let next = corners[(index + 1) % 4];
                        out.push((
                            (corners[index].0 + next.0) * 0.5,
                            (corners[index].1 + next.1) * 0.5,
                        ));
                    }
                }
                out
            })
            .collect()
    } else if edge_count > 5 {
        boundary_aligned_unit_square_loops(edge_count)
    } else {
        Vec::new()
    };
    oriented_loop_candidates(bases)
}

fn oriented_loop_candidates(bases: Vec<Vec<(f64, f64)>>) -> Vec<Vec<(f64, f64)>> {
    let mut candidates = Vec::new();
    for base in bases {
        for shift in 0..base.len() {
            let rotated = (0..base.len())
                .map(|index| base[(index + shift) % base.len()])
                .collect::<Vec<_>>();
            candidates.push(rotated.clone());
            candidates.push(rotated.into_iter().rev().collect());
        }
    }
    candidates
}

fn boundary_aligned_unit_square_loops(edge_count: usize) -> Vec<Vec<(f64, f64)>> {
    let extra_count = edge_count.saturating_sub(4);
    let mut loops = Vec::new();
    for bottom_extra in 0..=extra_count {
        for right_extra in 0..=extra_count - bottom_extra {
            for top_extra in 0..=extra_count - bottom_extra - right_extra {
                let left_extra = extra_count - bottom_extra - right_extra - top_extra;
                loops.push(boundary_aligned_unit_square_loop([
                    bottom_extra,
                    right_extra,
                    top_extra,
                    left_extra,
                ]));
            }
        }
    }
    loops
}

fn boundary_aligned_unit_square_loop(side_extra_counts: [usize; 4]) -> Vec<(f64, f64)> {
    let corners = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)];
    let mut out = Vec::with_capacity(4 + side_extra_counts.iter().sum::<usize>());
    for side in 0..4 {
        let start = corners[side];
        let end = corners[(side + 1) % 4];
        out.push(start);
        let extra = side_extra_counts[side];
        for split in 0..extra {
            let t = (split + 1) as f64 / (extra + 1) as f64;
            out.push(lerp_uv(start, end, t));
        }
    }
    out
}

fn face_outer_half_edges(brep: &BRep, face: FaceHandle) -> Result<Vec<HalfEdgeHandle>, CoreError> {
    let outer_loop = brep
        .faces
        .get(face)
        .ok_or_else(|| CoreError::InvalidHandle("kp20 support-cell pcurve face invalid".into()))?
        .outer_loop;
    let start = brep
        .loops
        .get(outer_loop)
        .ok_or_else(|| CoreError::InvalidHandle("kp20 support-cell pcurve loop invalid".into()))?
        .half_edge;
    let mut out = Vec::new();
    let mut cursor = start;
    for _ in 0..=brep.half_edges.len() {
        out.push(cursor);
        cursor = brep.half_edges[cursor].next;
        if cursor == start {
            return Ok(out);
        }
    }
    Err(CoreError::PreconditionFailed(
        "kp20 support-cell pcurve loop did not close".into(),
    ))
}
