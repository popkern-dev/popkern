//! PK-141 phase 2b separators: corner chords between rail-end
//! anchors (split from `multi_chord_region.rs`, line ratchet).

use std::collections::HashSet;

use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::CurveEval;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle};

use super::face_trim_pipeline::FaceTrim;

/// Mef separator chords between mutual-nearest rail-end
/// anchor PAIRS from different trims. Returns the separator edges.
pub(super) fn cut_corner_separators(
    brep: &mut BRep,
    trims: &[FaceTrim],
    sub_faces: &[FaceHandle],
    tolerance: &Tolerance,
) -> Vec<EdgeHandle> {
    let trace = popkern_core::trace::trace_enabled("multi-chord");
    // Rail-end anchor positions per trim.
    let mut ends: Vec<(usize, Point3)> = Vec::new();
    let mut spans: Vec<f64> = Vec::new();
    for (i, t) in trims.iter().enumerate() {
        let (t0, t1) = t.trim_3d.parameter_range();
        let (a, b) = (t.trim_3d.point_at(t0), t.trim_3d.point_at(t1));
        spans.push(a.distance_to(&b).max(1e-12));
        ends.push((i, a));
        ends.push((i, b));
    }
    let mut used: HashSet<usize> = HashSet::new();
    let mut separators = Vec::new();
    for x in 0..ends.len() {
        if used.contains(&x) {
            continue;
        }
        let (ti, pa) = ends[x];
        // Mutual nearest across trims.
        let mut best: Option<(usize, f64)> = None;
        for (y, &(tj, pb)) in ends.iter().enumerate() {
            if y == x || used.contains(&y) || tj == ti {
                continue;
            }
            let dxy = pa.distance_to(&pb);
            if best.map(|(_, bd)| dxy < bd).unwrap_or(true) {
                best = Some((y, dxy));
            }
        }
        let Some((y, d)) = best else { continue };
        let (tj, pb) = ends[y];
        // Corner gate: the pair must be corner-scale close (well
        // under either rail's span).
        let gate = spans[ti].min(spans[tj]) * 0.5;
        if d > gate || d < tolerance.linear().max(1e-9) {
            continue;
        }
        // The face whose one loop carries BOTH positions.
        let direct = loop_hes_at_positions(brep, sub_faces, pa, pb, tolerance);
        if direct.is_none() {
            // The anchors ended on different loops (splits assigned
            // them apart). Route the separator VIA the shared hole
            // corner: two chords, each within one loop.
            let mid = Point3::new(
                0.5 * (pa.x + pb.x),
                0.5 * (pa.y + pb.y),
                0.5 * (pa.z + pb.z),
            );
            let corner = nearest_vertex_position(brep, mid, d);
            let mut routed = 0;
            if let Some(pc) = corner {
                // ALL-OR-NOTHING: resolve both chords before cutting
                // either. A half-route (measured at the boss -y top
                // corner: corner->rail cut, corner->contact missing)
                // ships the lone chord as a detour edge 7e-3 off the
                // wall; when the second chord has no loop pair the
                // trims already partition the territory and no
                // separator is needed at that corner.
                let precheck = [(pa, pc), (pc, pb)].map(|(qa, qb)| {
                    loop_hes_at_positions(brep, sub_faces, qa, qb, tolerance).is_some()
                });
                if trace && precheck.iter().any(|ok| !ok) {
                    eprintln!(
                        "[mc-region] via-corner SKIPPED (chords resolvable: {precheck:?}): \
                         {pa:?} <-> {pb:?} via {pc:?}"
                    );
                }
                if precheck.iter().all(|&ok| ok) {
                    for (qa, qb) in [(pa, pc), (pc, pb)] {
                        // Re-resolve per chord: the first mef moves
                        // loop assignments.
                        if let Some((face, ha, hb)) =
                            loop_hes_at_positions(brep, sub_faces, qa, qb, tolerance)
                            && let Ok(mef) = brep.mef(ha, hb)
                        {
                            if let Some(e) = brep.edges.get_mut(mef.edge) {
                                e.curve = popkern_core::geometry::CurveOnEdge::line_segment(
                                    qa,
                                    qb,
                                    tolerance.linear(),
                                );
                            }
                            if trace {
                                eprintln!(
                                    "[mc-region] via-corner separator on {face:?}: {qa:?} -> {qb:?}"
                                );
                            }
                            separators.push(mef.edge);
                            routed += 1;
                        }
                    }
                }
            }
            if routed > 0 {
                used.insert(x);
                used.insert(y);
            } else if trace {
                eprintln!(
                    "[mc-region] separator pair SKIPPED (no common loop): {pa:?} <-> {pb:?} d={d:.4}"
                );
            }
            continue;
        }
        let Some((face, he_a, he_b)) = direct else {
            continue;
        };
        match brep.mef(he_a, he_b) {
            Ok(mef) => {
                if let Some(e) = brep.edges.get_mut(mef.edge) {
                    e.curve = popkern_core::geometry::CurveOnEdge::line_segment(
                        pa,
                        pb,
                        tolerance.linear(),
                    );
                }
                if trace {
                    eprintln!(
                        "[mc-region] separator on {face:?}: {pa:?} -> {pb:?} (edge {:?})",
                        mef.edge
                    );
                }
                separators.push(mef.edge);
                used.insert(x);
                used.insert(y);
                // Corner-seal unification (PK-141): also cut the
                // hole-corner chord (the chain-corner vertex to the
                // rail end) so the band seals from the territory
                // BEYOND the rail at this corner. The -y corners
                // carried this seal from vertex-blend routing while the +y
                // ones only had the contact separator — measured:
                // eg5's rail cut yielded band+upper as ONE piece and
                // the strip surface over-extended 0.7 to the rim.
                let mid2 = Point3::new(
                    0.5 * (pa.x + pb.x),
                    0.5 * (pa.y + pb.y),
                    0.5 * (pa.z + pb.z),
                );
                if let Some(pc) = nearest_vertex_position(brep, mid2, d)
                    && let Some((cface, ca, cb)) =
                        loop_hes_at_positions(brep, sub_faces, pc, pa, tolerance)
                    && let Ok(cmef) = brep.mef(ca, cb)
                {
                    if let Some(e) = brep.edges.get_mut(cmef.edge) {
                        e.curve = popkern_core::geometry::CurveOnEdge::line_segment(
                            pc,
                            pa,
                            tolerance.linear(),
                        );
                    }
                    if trace {
                        eprintln!("[mc-region] corner-seal on {cface:?}: {pc:?} -> {pa:?}");
                    }
                    separators.push(cmef.edge);
                }
            }
            Err(e) => {
                if trace {
                    eprintln!("[mc-region] separator mef failed on {face:?}: {e}");
                }
            }
        }
    }
    separators
}

/// Find one face (among the subdivision) with a single loop carrying
/// half-edges whose origins sit at BOTH positions.
pub(super) fn loop_hes_at_positions(
    brep: &BRep,
    sub_faces: &[FaceHandle],
    pa: Point3,
    pb: Point3,
    tolerance: &Tolerance,
) -> Option<(FaceHandle, HalfEdgeHandle, HalfEdgeHandle)> {
    // The recorded subdivision list can lag mid-surgery (re-homed
    // spurs and prior separator mefs mint pieces); fall through to
    // the full face set when the recorded list misses.
    loop_hes_at_positions_in(brep, sub_faces.iter().copied(), pa, pb, tolerance).or_else(|| {
        loop_hes_at_positions_in(
            brep,
            brep.faces.keys().collect::<Vec<_>>().into_iter(),
            pa,
            pb,
            tolerance,
        )
    })
}

fn loop_hes_at_positions_in(
    brep: &BRep,
    faces: impl Iterator<Item = FaceHandle>,
    pa: Point3,
    pb: Point3,
    tolerance: &Tolerance,
) -> Option<(FaceHandle, HalfEdgeHandle, HalfEdgeHandle)> {
    let tol = tolerance.linear().max(1e-9) * 100.0;
    for f in faces {
        let fd = brep.faces.get(f)?;
        let mut loops = vec![fd.outer_loop];
        loops.extend(fd.inner_loops.iter().copied());
        for lh in loops {
            let l = brep.loops.get(lh)?;
            let start = l.half_edge;
            let (mut ha, mut hb) = (None, None);
            let mut cur = start;
            for _ in 0..brep.half_edges.len() + 1 {
                let he = brep.half_edges.get(cur)?;
                if let Some(v) = brep.vertices.get(he.origin) {
                    if v.position.distance_to(&pa) < tol {
                        ha = Some(cur);
                    }
                    if v.position.distance_to(&pb) < tol {
                        hb = Some(cur);
                    }
                }
                cur = he.next;
                if cur == start {
                    break;
                }
            }
            if let (Some(a), Some(b)) = (ha, hb)
                && a != b
            {
                return Some((f, a, b));
            }
        }
    }
    None
}

/// Nearest existing vertex to `p` within `gate`.
fn nearest_vertex_position(brep: &BRep, p: Point3, gate: f64) -> Option<Point3> {
    let mut best: Option<(Point3, f64)> = None;
    for (_, v) in brep.vertices.iter() {
        let d = v.position.distance_to(&p);
        if d <= gate && best.map(|(_, bd)| d < bd).unwrap_or(true) {
            best = Some((v.position, d));
        }
    }
    best.map(|(q, _)| q)
}
