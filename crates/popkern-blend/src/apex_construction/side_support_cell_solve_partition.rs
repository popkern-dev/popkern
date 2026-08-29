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
    let surface = brep.faces.get(face).and_then(|f| f.surface.as_ref())?;
    let samples = loop_samples(brep, face).ok()?;
    let mut worst = 0.0f64;
    for sample in &samples {
        for p in [Some(sample.origin), sample.mid].into_iter().flatten() {
            worst = worst.max(domain_clamped_projection_distance(&surface.geometry, &p));
        }
    }
    Some(worst)
}

pub(super) fn face_surface_misfits_boundary(brep: &BRep, face: FaceHandle, gate: f64) -> bool {
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
            .any(|p| domain_clamped_projection_distance(&surface.geometry, &p) > gate)
    })
}

/// Distance from a point to the surface with the projected UV clamped
/// into the surface's finite parameter domain. Unclamped Gauss-Newton
/// projection on a Nurbs converges OUTSIDE the domain and `point_at`
/// extrapolates to meet the point, reporting a vacuously small misfit
/// (measured: cube_extreme cap 5v3, boundary vertex 1.67e-1 off the
/// clamped surface but ~0 off the extrapolation).
fn domain_clamped_projection_distance(
    surface: &popkern_core::geometry::SurfaceGeometry,
    p: &popkern_core::geometry::Point3,
) -> f64 {
    use popkern_core::geometry::surface::SurfaceEval;
    let Some((mut u, mut v)) =
        popkern_core::geometry::uv_projection::project_point_to_uv(surface, p)
    else {
        return f64::INFINITY;
    };
    let ((u_lo, u_hi), (v_lo, v_hi)) = surface.parameter_range();
    if u_lo.is_finite() && u_hi.is_finite() {
        u = u.clamp(u_lo, u_hi);
    }
    if v_lo.is_finite() && v_hi.is_finite() {
        v = v.clamp(v_lo, v_hi);
    }
    surface.point_at(u, v).distance_to(p)
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

/// One-shot diagnostic census (trace-gated): every edge whose two
/// half-edges are not BOTH members of their loops' next-cycles, with
/// face pair, topology roles, and origin — maps the unlive-twin
/// records left when a boundary rewrite re-homes a new cycle without
/// retiring the displaced half-edges (measured on cube_extreme:
/// 66/108 edges broken at acceptance time).
pub(super) fn trace_edge_liveness_census(brep: &BRep) {
    let cycle_live = |he: popkern_core::topology::HalfEdgeHandle| -> bool {
        let Some(data) = brep.half_edges.get(he) else {
            return false;
        };
        let Some(loop_data) = brep.loops.get(data.loop_) else {
            return false;
        };
        let head = loop_data.half_edge;
        let mut cursor = head;
        for _ in 0..brep.half_edges.len() {
            if cursor == he {
                return true;
            }
            let Some(c) = brep.half_edges.get(cursor) else {
                return false;
            };
            cursor = c.next;
            if cursor == head {
                return false;
            }
        }
        false
    };
    let mut both = 0usize;
    let mut broken = 0usize;
    for (edge_handle, edge) in brep.edges.iter() {
        let a = edge.half_edge;
        let b = brep.half_edges.get(a).map(|d| d.twin);
        let a_live = cycle_live(a);
        let b_live = b.is_some_and(cycle_live);
        if a_live && b_live {
            both += 1;
            continue;
        }
        broken += 1;
        let face_of = |he: Option<popkern_core::topology::HalfEdgeHandle>| {
            he.and_then(|h| brep.half_edges.get(h))
                .and_then(|d| brep.loops.get(d.loop_))
                .map(|l| l.face)
        };
        let origin = brep
            .half_edges
            .get(a)
            .and_then(|d| brep.vertices.get(d.origin))
            .map(|v| v.position);
        let role_of = |fh: Option<FaceHandle>| fh.and_then(|f| brep.faces.get(f)).map(|f| f.role);
        eprintln!(
            "[direct-apex-trace] census {edge_handle:?}: a_live={a_live} b_live={b_live} \
             faces {:?}/{:?} roles {:?}/{:?} edge_role={:?} at {:?}",
            face_of(Some(a)),
            face_of(b),
            role_of(face_of(Some(a))),
            role_of(face_of(b)),
            edge.role,
            origin.map(|p| format!("({:.4},{:.4},{:.4})", p.x, p.y, p.z)),
        );
    }
    eprintln!("[direct-apex-trace] census: {both} edges fully live, {broken} broken");
}
