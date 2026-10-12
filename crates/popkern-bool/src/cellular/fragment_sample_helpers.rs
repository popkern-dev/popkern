//! Sampling helpers shared by the cellular fragment classifier
//! (PK-496 split from `fragment_classify_one.rs` under the 500-line
//! budget): edge-proximity tests, boundary chord midpoints, and the
//! tuned edge-midpoint fragment sample shortcut.

use popkern_core::topology::BRep;

use popkern_core::geometry::intersection::surface_subdivision::Aabb;
use popkern_core::tolerance::Tolerance;

use crate::merge::graph::IntersectionGraph;
use crate::merge::walk::Fragment;

pub(super) fn point_near_any_edge(
    brep: &BRep,
    p: &popkern_core::geometry::Point3,
    band: f64,
) -> bool {
    use popkern_core::geometry::curve::CurveEval as _;
    for (_, edge) in brep.edges.iter() {
        let Some(he) = brep.half_edges.get(edge.half_edge) else {
            continue;
        };
        let a = brep.vertices[he.origin].position;
        let b = brep
            .half_edges
            .get(he.twin)
            .map(|t| brep.vertices[t.origin].position)
            .unwrap_or(a);
        if let Some(c) = edge.curve.as_ref() {
            // Every chord below lies in the curve's conservative box,
            // so a point `band` outside it is near none of them (PK-743:
            // each fragment sample walked every edge of a 64-slot disc).
            if let Some(bx) = crate::coincidence::edge_box::conservative_curve_box(
                &c.geometry,
                (c.t_start, c.t_end),
            ) && (p.x < bx.min.x - band
                || p.y < bx.min.y - band
                || p.z < bx.min.z - band
                || p.x > bx.max.x + band
                || p.y > bx.max.y + band
                || p.z > bx.max.z + band)
            {
                continue;
            }
            let mut prev = c.geometry.point_at(c.t_start);
            for i in 1..=16 {
                let t = c.t_start + (c.t_end - c.t_start) * i as f64 / 16.0;
                let q = c.geometry.point_at(t);
                if seg_dist(&prev, &q, p) < band {
                    return true;
                }
                prev = q;
            }
        } else if seg_dist(&a, &b, p) < band {
            return true;
        }
    }
    false
}

pub(super) fn seg_dist(
    a: &popkern_core::geometry::Point3,
    b: &popkern_core::geometry::Point3,
    p: &popkern_core::geometry::Point3,
) -> f64 {
    let ab = *b - *a;
    let l2 = ab.dot(&ab);
    if l2 < 1e-24 {
        return p.distance_to(a);
    }
    let t = ((*p - *a).dot(&ab) / l2).clamp(0.0, 1.0);
    (*a + ab * t).distance_to(p)
}

/// Chord midpoints of the loop's non-Shared half-edges, in loop
/// order, capped at `cap`. The first entry equals the tuned
/// edge-midpoint shortcut, so a unanimous vote reproduces the uncalibrated shortcut
/// sampling exactly.
/// Midpoint of a graph half-edge, ON THE EDGE'S CURVE when it has
/// one (chord midpoint otherwise). A curved edge's chord midpoint is
/// not a point of the fragment at all — a shallow rim-arc sliver cut
/// off by a breached wall has its arc CHORD lying exactly IN the
/// cutting plane, so the chord midpoint sampled the parity coin-flip
/// point on the other body's surface while the arc's belly sits
/// decisively outside it (measured: the sliver classified Inside and
/// a subtract kept a flap beyond the body). Straight edges are
/// unchanged (chord == curve).
pub(super) fn he_curve_midpoint(
    graph: &IntersectionGraph,
    he_idx: usize,
) -> popkern_core::geometry::Point3 {
    use popkern_core::geometry::curve::CurveEval as _;
    let he = &graph.half_edges[he_idx];
    if let Some(c) = he
        .curve
        .as_ref()
        .or_else(|| graph.half_edges[he.twin].curve.as_ref())
    {
        return c.geometry.point_at(0.5 * (c.t_start + c.t_end));
    }
    let a = graph.nodes[he.origin].position;
    let b = graph.nodes[graph.half_edges[he.twin].origin].position;
    popkern_core::geometry::Point3::new((a.x + b.x) * 0.5, (a.y + b.y) * 0.5, (a.z + b.z) * 0.5)
}

pub(super) fn non_shared_chord_midpoints(
    graph: &IntersectionGraph,
    loop_hes: &[usize],
    cap: usize,
) -> Vec<popkern_core::geometry::Point3> {
    use crate::merge::graph::EdgeSource;
    let mut out = Vec::new();
    for &he_idx in loop_hes {
        if out.len() >= cap {
            break;
        }
        let he = &graph.half_edges[he_idx];
        if matches!(he.source, EdgeSource::Shared { .. }) {
            continue;
        }
        out.push(he_curve_midpoint(graph, he_idx));
    }
    out
}

pub(super) fn sample_fragment_point(
    graph: &IntersectionGraph,
    loop_hes: &[usize],
) -> popkern_core::geometry::Point3 {
    use crate::merge::graph::EdgeSource;

    // Prefer midpoint of a non-shared edge (away from intersection curve)
    for &he_idx in loop_hes {
        let he = &graph.half_edges[he_idx];
        match &he.source {
            EdgeSource::OriginalA(_)
            | EdgeSource::OriginalB(_)
            | EdgeSource::SplitA { .. }
            | EdgeSource::SplitB { .. } => {
                return he_curve_midpoint(graph, he_idx);
            }
            EdgeSource::Shared { .. } => continue,
        }
    }

    // Fallback: centroid of all loop vertices
    if loop_hes.is_empty() {
        return popkern_core::geometry::Point3::new(0.0, 0.0, 0.0);
    }
    let mut sx = 0.0;
    let mut sy = 0.0;
    let mut sz = 0.0;
    for &he_idx in loop_hes {
        let p = graph.nodes[graph.half_edges[he_idx].origin].position;
        sx += p.x;
        sy += p.y;
        sz += p.z;
    }
    let n = loop_hes.len() as f64;
    popkern_core::geometry::Point3::new(sx / n, sy / n, sz / n)
}

/// DEGENERATE loops (a lone non-closed half-edge, or a he/twin
/// bigon) enclose no area — boundary-graze residue. Keeping one
/// produced a dangling face in the assembly.
pub(super) fn is_degenerate_fragment(
    graph: &IntersectionGraph,
    fragment: &Fragment,
    tolerance: &Tolerance,
) -> bool {
    match fragment.outer_loop.len() {
        1 => {
            let he = &graph.half_edges[fragment.outer_loop[0]];
            if he.origin != graph.half_edges[he.twin].origin {
                true
            } else {
                // A lone SELF-LOOP is a real fragment boundary when
                // it is a closed intersection curve (a circle cut
                // through a face). A ZERO-EXTENT self-loop is a
                // tangency-point record (a surface grazing another
                // at one point) — keeping it assembled a point-face
                // that broke the solid's Euler count. Cull by
                // sampled curve extent against the tolerance.
                let extent = he
                    .curve
                    .as_ref()
                    .or_else(|| graph.half_edges[he.twin].curve.as_ref())
                    .map(|c| {
                        use popkern_core::geometry::curve::CurveEval;
                        let origin = graph.nodes[he.origin].position;
                        (1..8)
                            .map(|i| {
                                let t = c.t_start + (c.t_end - c.t_start) * i as f64 / 8.0;
                                c.geometry.point_at(t).distance_to(&origin)
                            })
                            .fold(0.0, f64::max)
                    })
                    .unwrap_or(0.0);
                extent < tolerance.linear() * 10.0
            }
        }
        2 => {
            let a = fragment.outer_loop[0];
            if graph.half_edges[a].twin == fragment.outer_loop[1] {
                // A he/twin bigon is graze residue only for RECORD
                // pairs. A BOUNDARY pair walked as one loop is a
                // periodic seam spine — the outer loop of a full
                // periodic chart (a full sphere's two meridian
                // copies bound the whole wall face, with the
                // imprint loops as its holes).
                matches!(
                    graph.half_edges[a].source,
                    crate::merge::graph::EdgeSource::Shared { .. }
                )
            } else {
                // An OPEN two-half-edge walk with (sampled) zero
                // area is residue of a seeded walk that swallowed
                // a tangency self-loop and dead-ended — e.g. the
                // unconsumed direction of a boundary-collinear
                // corner record plus the zero-extent point-contact
                // record at its far end. It can never assemble
                // (measured: kept as OnBoundary, it became an
                // unclosed 20 mm loop and the whole flush-arm
                // chain fell to graph-merge). Real two-arc bigons
                // (a split circle) have genuine sampled area and
                // stay.
                let start = graph.half_edges[fragment.outer_loop[0]].origin;
                let last = fragment.outer_loop[1];
                let end = graph.half_edges[graph.half_edges[last].twin].origin;
                let open = start != end;
                let area =
                    crate::merge::walk::walk_geometry::loop_newell(graph, &fragment.outer_loop)
                        .length();
                open && area < (tolerance.linear() * 10.0).powi(2)
            }
        }
        _ => false,
    }
}

/// Whether `p` lies within `band` of a face surface of `brep` (face
/// boxes pre-filter, then the exact closest point on the carrier). The
/// carrier, not the trimmed face, decides, so a hit can over-report;
/// callers use it only to switch to an interior sample, never to
/// decide a class (PK-750).
pub(super) fn point_near_any_surface(
    brep: &BRep,
    aabbs: &[(popkern_core::topology::FaceHandle, Aabb)],
    p: &popkern_core::geometry::Point3,
    band: f64,
) -> bool {
    use popkern_core::geometry::surface::SurfaceEval as _;
    let probe = Aabb { min: *p, max: *p };
    aabbs.iter().any(|(face, b)| {
        b.overlaps(&probe, band)
            && brep
                .faces
                .get(*face)
                .and_then(|f| f.surface.as_ref())
                .and_then(|s| {
                    let (u, v) = s.geometry.closest_point(p)?;
                    Some(s.geometry.point_at(u, v).distance_to(p) < band)
                })
                .unwrap_or(false)
    })
}
