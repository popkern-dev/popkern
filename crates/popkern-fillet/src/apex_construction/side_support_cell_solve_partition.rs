//! Boundary-run labeling and refit predicates for the side-support
//! support-cell solve (split out of `side_support_cell_solve.rs` for the
//! 500-line budget).

use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::topology::{BRep, FaceHandle, HalfEdgeHandle};

use super::geometry::same_carrier;
use super::loop_samples;
use super::{Carrier, HeSample};

pub(super) fn face_has_live_twin(brep: &BRep, face: FaceHandle) -> bool {
    let Ok(samples) = loop_samples(brep, face) else {
        return false;
    };
    samples.iter().any(|sample| {
        brep.half_edges
            .get(sample.he)
            .and_then(|d| brep.half_edges.get(d.twin))
            .and_then(|t| brep.loops.get(t.loop_))
            .and_then(|l| brep.faces.get(l.face))
            .is_some()
    })
}

pub(super) fn worst_boundary_misfit(brep: &BRep, face: FaceHandle) -> Option<f64> {
    use popkern_core::geometry::surface::SurfaceEval;
    let surface = brep.faces.get(face).and_then(|f| f.surface.as_ref())?;
    let samples = loop_samples(brep, face).ok()?;
    let mut worst = 0.0f64;
    for sample in &samples {
        for p in [Some(sample.origin), sample.mid].into_iter().flatten() {
            let d =
                popkern_core::geometry::uv_projection::project_point_to_uv(&surface.geometry, &p)
                    .map(|(u, v)| surface.geometry.point_at(u, v).distance_to(&p))
                    .unwrap_or(f64::INFINITY);
            worst = worst.max(d);
        }
    }
    Some(worst)
}

pub(super) fn face_surface_misfits_boundary(brep: &BRep, face: FaceHandle, gate: f64) -> bool {
    use popkern_core::geometry::surface::SurfaceEval;
    let Some(surface) = brep.faces.get(face).and_then(|f| f.surface.as_ref()) else {
        return false;
    };
    let Ok(samples) = loop_samples(brep, face) else {
        return false;
    };
    samples.iter().any(|sample| {
        [Some(sample.origin), sample.mid]
            .into_iter()
            .flatten()
            .any(|p| {
                popkern_core::geometry::uv_projection::project_point_to_uv(&surface.geometry, &p)
                    .map(|(u, v)| surface.geometry.point_at(u, v).distance_to(&p) > gate)
                    .unwrap_or(true)
            })
    })
}

pub(super) fn assign_labels(samples: &[HeSample], n_carriers: usize) -> Vec<usize> {
    let bare = n_carriers;
    let n = samples.len();
    let options: Vec<Vec<usize>> = samples
        .iter()
        .map(|s| {
            let opts: Vec<usize> = (0..n_carriers)
                .filter(|c| s.carrier_mask >> c & 1 == 1)
                .collect();
            if opts.is_empty() { vec![bare] } else { opts }
        })
        .collect();
    // DP over the circular sequence for each choice of first label.
    let mut best: Option<(usize, Vec<usize>)> = None;
    for &first in &options[0] {
        let mut cost = vec![std::collections::BTreeMap::new(); n];
        cost[0].insert(first, (0usize, first));
        for i in 1..n {
            let (prev, rest) = cost.split_at_mut(i);
            let prev = &prev[i - 1];
            for &label in &options[i] {
                let mut cheapest: Option<(usize, usize)> = None;
                for (&pl, &(pc, _)) in prev.iter() {
                    let c = pc + usize::from(pl != label);
                    if cheapest.is_none_or(|(bc, _)| c < bc) {
                        cheapest = Some((c, pl));
                    }
                }
                if let Some((c, from)) = cheapest {
                    rest[0].insert(label, (c, from));
                }
            }
        }
        // Close the circle: transition back to `first` counts.
        let Some((&last_label, &(last_cost, _))) = cost[n - 1]
            .iter()
            .min_by_key(|entry| entry.1.0 + usize::from(*entry.0 != first))
        else {
            continue;
        };
        let total = last_cost + usize::from(last_label != first);
        if best.as_ref().is_none_or(|(bc, _)| total < *bc) {
            // Reconstruct.
            let mut labels = vec![0usize; n];
            labels[n - 1] = last_label;
            for i in (1..n).rev() {
                let (_, from) = cost[i][&labels[i]];
                labels[i - 1] = from;
            }
            best = Some((total, labels));
        }
    }
    best.map(|(_, labels)| labels)
        .unwrap_or_else(|| vec![bare; n])
}

/// Choose the carrier run with the most REAL (non-authored)
/// half-edges, tie-broken by total length. Picking the FIRST run at a
/// label boundary favored 1-half-edge runs whose bigon splits are
/// size-neutral (remove one real, mint one closing) and burned the
/// split budget without converging (measured on cube_extreme's 10v5:
/// 13 splits on a fully carrier-covered loop that barely shrank).
pub(super) fn first_carrier_run(
    labels: &[usize],
    bare: usize,
    samples: &[HeSample],
    authored: &std::collections::BTreeSet<HalfEdgeHandle>,
) -> Option<(usize, usize, usize)> {
    let n = labels.len();
    // Find an index where the label changes (run boundary).
    let boundary = (0..n).find(|&i| labels[i] != labels[(i + n - 1) % n])?;
    let mut best: Option<(usize, usize, usize, usize)> = None;
    let mut i = boundary;
    let mut walked = 0usize;
    while walked < n {
        let label = labels[i];
        let mut len = 0usize;
        while len < n && labels[(i + len) % n] == label {
            len += 1;
        }
        let real = (0..len)
            .filter(|&k| !authored.contains(&samples[(i + k) % n].he))
            .count();
        if label != bare
            && real > 0
            && best.is_none_or(|(_, blen, _, breal)| real > breal || (real == breal && len > blen))
        {
            best = Some((i, len, label, real));
        }
        i = (i + len) % n;
        walked += len;
    }
    best.map(|(start, len, label, _)| (start, len, label))
}

pub(super) fn collect_carriers(brep: &BRep, bridge_faces: &[FaceHandle]) -> Vec<Carrier> {
    let mut frontier: Vec<FaceHandle> = bridge_faces.to_vec();
    let mut seen: std::collections::BTreeSet<FaceHandle> = frontier.iter().copied().collect();
    for _hop in 0..3 {
        let mut next = Vec::new();
        for &face in &frontier {
            let Some(loop_) = brep
                .faces
                .get(face)
                .and_then(|f| brep.loops.get(f.outer_loop))
            else {
                continue;
            };
            let start = loop_.half_edge;
            let mut he = start;
            for _ in 0..brep.half_edges.len() {
                let Some(data) = brep.half_edges.get(he) else {
                    break;
                };
                if let Some(twin_face) = brep
                    .half_edges
                    .get(data.twin)
                    .and_then(|t| brep.loops.get(t.loop_))
                    .map(|l| l.face)
                    && seen.insert(twin_face)
                {
                    next.push(twin_face);
                }
                he = data.next;
                if he == start {
                    break;
                }
            }
        }
        frontier = next;
    }
    let mut carriers: Vec<Carrier> = Vec::new();
    for face in seen {
        let Some(surface) = brep.faces.get(face).and_then(|f| f.surface.as_ref()) else {
            continue;
        };
        if !matches!(
            surface.geometry,
            SurfaceGeometry::Plane { .. } | SurfaceGeometry::Cylinder { .. }
        ) {
            continue;
        }
        if carriers
            .iter()
            .any(|c| same_carrier(&c.surface, &surface.geometry))
        {
            continue;
        }
        carriers.push(Carrier {
            surface: surface.geometry.clone(),
            source: face,
        });
    }
    carriers
}
