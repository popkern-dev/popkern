//! Permutation walk over a face's chart subdivision, and the dual-walk
//! comparator (PK-104 stages 2–3 — see
//! docs/architecture/fragment-walk-chart-subdivision.md §3.2–3.3).
//!
//! `next()` is the classical DCEL rotational successor for
//! material-on-left traversal: among half-edges departing the branch
//! point where `h` arrives, take the one with the smallest strictly
//! positive CLOCKWISE angle from `h`'s backward direction in oriented
//! UV, twin ranked last. It is a pure function of the half-edge, so
//! loops are the orbits of a permutation — the partition cannot depend
//! on seed order.

use std::collections::BTreeSet;

use super::chart_subdivision::FaceChartSubdivision;
use crate::merge::graph::IntersectionGraph;

#[path = "chart_walk_next.rs"]
mod chart_walk_next;
pub(crate) use chart_walk_next::next_in_chart;

/// Enumerate the permutation's orbits, ascending from the lowest
/// unconsumed half-edge (deterministic, seed-order independent by
/// construction). Returns `None` when any `next` fails or an orbit
/// fails to close within the face's half-edge budget.
pub(crate) fn walk_chart_orbits(
    sub: &FaceChartSubdivision,
    graph: &IntersectionGraph,
) -> Option<Vec<Vec<usize>>> {
    let mut consumed: BTreeSet<usize> = BTreeSet::new();
    let mut consumed_by: std::collections::HashMap<usize, usize> = std::collections::HashMap::new();
    let mut orbits = Vec::new();
    let budget = sub.embeddings.len() + 1;
    for &seed in sub.embeddings.keys() {
        if consumed.contains(&seed) {
            continue;
        }
        let mut orbit = vec![seed];
        consumed.insert(seed);
        consumed_by.insert(seed, seed);
        let mut h = seed;
        loop {
            let Some(n) = next_in_chart(sub, graph, h) else {
                if popkern_core::trace::trace_enabled("fragment-walk") {
                    eprintln!("      [orbit-fail] no successor from he {h}");
                }
                return None;
            };
            if n == seed {
                break;
            }
            // A permutation orbit cannot revisit a non-seed member; if
            // it does, the substrate lied — bail to fallback.
            if consumed.contains(&n) || orbit.len() > budget {
                if popkern_core::trace::trace_enabled("fragment-walk") {
                    eprintln!(
                        "      [orbit-fail] revisit/budget: he {h} -> {n} (seed {seed}, len {}, {n} consumed by orbit {:?})",
                        orbit.len(),
                        consumed_by.get(&n)
                    );
                    if let Some(bp) = sub.end_bp.get(&h) {
                        for (hh, e) in sub.end_bp.iter() {
                            if e == bp {
                                eprintln!(
                                    "        arrival {hh}: end_dir {:?}",
                                    sub.embeddings.get(hh).and_then(|em| em.end_dir())
                                );
                            }
                        }
                        let arrive = sub.embeddings.get(&h).and_then(|e| e.end_dir());
                        eprintln!("        arrive dir {arrive:?} at bp {bp:?}");
                        if let Some(deps) = sub.departures.get(bp) {
                            for &c in deps {
                                let d = sub.embeddings.get(&c).and_then(|e| e.start_dir());
                                let uv0 = sub.embeddings.get(&c).map(|e| e.uv[0]);
                                let uvl = sub.embeddings.get(&c).and_then(|e| e.uv.last().copied());
                                eprintln!(
                                    "        dep {c}: dir {d:?} uv {uv0:?}->{uvl:?} consumed={}",
                                    consumed.contains(&c)
                                );
                            }
                        }
                    }
                }
                return None;
            }
            orbit.push(n);
            consumed.insert(n);
            consumed_by.insert(n, seed);
            h = n;
        }
        orbits.push(orbit);
    }
    Some(orbits)
}

/// A chart polyline: sampled oriented-UV points in walk order.
pub(crate) type ChartPolyline = Vec<(f64, f64)>;

/// Concatenated chart polyline of an orbit: the uv points of each
/// half-edge's embedding in walk order (consecutive embeddings share
/// their junction point — skipped after the first).
pub(crate) fn orbit_polyline(sub: &FaceChartSubdivision, orbit: &[usize]) -> ChartPolyline {
    let mut pts: Vec<(f64, f64)> = Vec::new();
    for &h in orbit {
        if let Some(emb) = sub.embeddings.get(&h) {
            let skip = usize::from(!pts.is_empty());
            pts.extend(emb.uv.iter().skip(skip).copied());
        }
    }
    pts
}

/// Signed shoelace area of a chart polyline (closed implicitly).
pub(crate) fn polyline_shoelace(pts: &[(f64, f64)]) -> f64 {
    let n = pts.len();
    let mut s = 0.0;
    for i in 0..n {
        let (x0, y0) = pts[i];
        let (x1, y1) = pts[(i + 1) % n];
        s += x0 * y1 - x1 * y0;
    }
    0.5 * s
}

/// Even-odd containment of a chart point in a chart polyline.
fn polyline_contains(pts: &[(f64, f64)], p: (f64, f64)) -> bool {
    let n = pts.len();
    let mut inside = false;
    for i in 0..n {
        let (x0, y0) = pts[i];
        let (x1, y1) = pts[(i + 1) % n];
        if (y0 > p.1) != (y1 > p.1) {
            let xi = x0 + (p.1 - y0) * (x1 - x0) / (y1 - y0);
            if xi > p.0 {
                inside = !inside;
            }
        }
    }
    inside
}

/// Role assembly for a chart-walked face's orbit partition.
///
/// Each region of the face's planar subdivision contributes one orbit
/// per boundary COMPONENT: orbits winding WITH the oriented chart
/// (positive σ-folded shoelace) are region outer components; orbits
/// winding AGAINST are hole components of the region that contains
/// them. This is by construction of the material-left walk — no 3D
/// election needed (the 3D-Newell election degenerates on faces
/// whose stored border is a seam pair: a full sphere's two
/// anti-parallel meridian copies measure |N| ≈ 0, and the seam orbit
/// lost the outer election to a record triangle on cube∪sphere).
pub(crate) struct ChartRoleAssembly {
    /// One entry per region: (outer orbit index, hole orbit indices).
    /// The FIRST entry is the face's primary fragment (largest
    /// positive region — same size rule the 3D-Newell election used).
    pub regions: Vec<(usize, Vec<usize>)>,
}

/// Assign orbit roles from σ-folded chart shoelace signs. Returns
/// `None` when no orbit winds positive (substrate defect — caller
/// falls back to the 3D-Newell election).
/// The source face's inner-loop edge set — the provenance authority
/// for the hole override in [`assemble_chart_roles`] (PK-341).
pub(crate) fn source_inner_edge_set(
    face_ref: &super::super::graph::FaceRef,
    brep_a: &popkern_core::topology::BRep,
    brep_b: &popkern_core::topology::BRep,
) -> std::collections::BTreeSet<popkern_core::topology::EdgeHandle> {
    source_loop_edge_set(face_ref, brep_a, brep_b, false)
}

/// The edge set of the source face's OUTER loop (see
/// [`source_inner_edge_set`] for the inner-loop dual).
pub(crate) fn source_outer_edge_set(
    face_ref: &super::super::graph::FaceRef,
    brep_a: &popkern_core::topology::BRep,
    brep_b: &popkern_core::topology::BRep,
) -> std::collections::BTreeSet<popkern_core::topology::EdgeHandle> {
    source_loop_edge_set(face_ref, brep_a, brep_b, true)
}

fn source_loop_edge_set(
    face_ref: &super::super::graph::FaceRef,
    brep_a: &popkern_core::topology::BRep,
    brep_b: &popkern_core::topology::BRep,
    outer: bool,
) -> std::collections::BTreeSet<popkern_core::topology::EdgeHandle> {
    let source_brep = match face_ref.operand {
        super::super::graph::Operand::A => brep_a,
        super::super::graph::Operand::B => brep_b,
    };
    let mut edges = std::collections::BTreeSet::new();
    if let Some(face) = source_brep.faces.get(face_ref.face) {
        let loops: Vec<_> = if outer {
            vec![face.outer_loop]
        } else {
            face.inner_loops.clone()
        };
        for lh in loops {
            let Some(loop_data) = source_brep.loops.get(lh) else {
                continue;
            };
            let start = loop_data.half_edge;
            let mut cursor = start;
            for _ in 0..source_brep.half_edges.len().max(1) {
                let Some(he) = source_brep.half_edges.get(cursor) else {
                    break;
                };
                edges.insert(he.edge);
                cursor = he.next;
                if cursor == start {
                    break;
                }
            }
        }
    }
    edges
}

pub(crate) fn assemble_chart_roles(
    graph: &IntersectionGraph,
    loops: &[Vec<usize>],
    polys: &[ChartPolyline],
    sigma: f64,
    u_period: Option<f64>,
    source_inner_edges: &std::collections::BTreeSet<popkern_core::topology::EdgeHandle>,
    source_outer_edges: &std::collections::BTreeSet<popkern_core::topology::EdgeHandle>,
) -> Option<ChartRoleAssembly> {
    let mut areas: Vec<f64> = polys.iter().map(|p| sigma * polyline_shoelace(p)).collect();
    // A chart-WRAPPING orbit — net u displacement of one closure
    // period across the unfolded seam — is NON-CONTRACTIBLE: it
    // cannot bound a disk, so it can never take the hole role. Its
    // implicit-closure shoelace sign is phase luck (the closing jump
    // contributes ±period·(u_wrap − mean-u of the loop)): a coil
    // wall's seam-crossing cut oval read negative at plane z=0.85,
    // filed as an inner hole, the wall never partitioned, one
    // interior sample above the plane classified the whole wall
    // Inside, and an EMPTY body SERVED (PK-452; z=0.9 served only
    // because the sign happened to land positive). Wrapping orbits
    // are region boundaries by construction; the provenance override
    // below still outranks this (sign only decides where both are
    // silent).
    // Wrap detection is the chart u-WINDING NUMBER: the embedding may
    // fold the seam crossing back into the chart (the coil oval's
    // polyline closes at u = 0 exactly), so first-vs-last displacement
    // reads 0 — but the interior samples still traverse the full
    // period with one fold jump. Minimal-image du maps fold jumps to
    // ~0 while genuine traversal accumulates ±period; a nonzero
    // rounded winding is a wrap. Embedding samples are dense pcurve
    // samples, well under half a period apart, so the minimal image
    // is faithful.
    let wrapping: Vec<bool> = match u_period {
        Some(p) if p > 0.0 => polys
            .iter()
            .map(|poly| {
                let mut acc = 0.0;
                for w in poly.windows(2) {
                    let du = w[1].0 - w[0].0;
                    acc += du - p * (du / p).round();
                }
                if let (Some(a), Some(b)) = (poly.first(), poly.last()) {
                    let du = a.0 - b.0;
                    acc += du - p * (du / p).round();
                }
                (acc / p).round() as i64 != 0
            })
            .collect(),
        _ => vec![false; polys.len()],
    };
    for (i, a) in areas.iter_mut().enumerate() {
        if wrapping[i] && *a < 0.0 {
            *a = -*a;
        }
    }
    if popkern_core::trace::trace_enabled("boolean") {
        for (i, poly) in polys.iter().enumerate() {
            eprintln!(
                "  [roles] orbit {i}: len {} area {:.5} wrap={} uv {:?}->{:?} (u_period {u_period:?})",
                loops[i].len(),
                areas[i],
                wrapping[i],
                poly.first(),
                poly.last()
            );
        }
    }
    // Provenance override: an orbit consisting ENTIRELY of the source
    // face's INNER-loop edges is a hole BY CONSTRUCTION — the chart
    // embedding's shoelace sign for unpaired closed rims was measured
    // flipping per cap orientation (PK-341: the tube's down-facing
    // cap read its r = 0.3 rim positive, minted a spurious disc
    // region over the hole, and the assembly double-covered the
    // core). Sign only decides where provenance is silent.
    for (i, orbit) in loops.iter().enumerate() {
        if areas[i] <= 0.0 || orbit.is_empty() {
            continue;
        }
        let all_inner = orbit.iter().all(|&h| {
            use super::super::graph::EdgeSource;
            match &graph.half_edges[h].source {
                EdgeSource::OriginalA(e) | EdgeSource::OriginalB(e) => {
                    source_inner_edges.contains(e)
                }
                EdgeSource::SplitA { original } | EdgeSource::SplitB { original } => {
                    source_inner_edges.contains(original)
                }
                _ => false,
            }
        });
        if all_inner {
            areas[i] = -areas[i];
        }
    }
    // Dual provenance override (PK-480): an orbit carrying ANY of the
    // source face's stored OUTER-loop pieces bounds material at the
    // face's outermost trim — it is a region OUTER by construction
    // and can never take the hole role, whatever its σ-folded sign
    // reads (a hole ring is interior to its region; the outer trim
    // cannot be). Measured on the angled tap's crest strip: the
    // stored boundary winds heterogeneously (outer pieces plus a
    // mouth window's pieces, rewired across the entry-cone circle
    // cut), the kept band's ring read −0.647 against its discard-side
    // sibling's +0.189, the band filed as the sliver's hole, both
    // discarded together, and the circle record's strip side
    // stranded (honest decline, one half-claimed side).
    for (i, orbit) in loops.iter().enumerate() {
        if areas[i] > 0.0 || orbit.is_empty() {
            continue;
        }
        let any_outer = orbit.iter().any(|&h| {
            use super::super::graph::EdgeSource;
            match &graph.half_edges[h].source {
                EdgeSource::OriginalA(e) | EdgeSource::OriginalB(e) => {
                    source_outer_edges.contains(e)
                }
                EdgeSource::SplitA { original } | EdgeSource::SplitB { original } => {
                    source_outer_edges.contains(original)
                }
                _ => false,
            }
        });
        if any_outer {
            areas[i] = -areas[i];
        }
    }
    let primary = areas
        .iter()
        .enumerate()
        .filter(|(_, a)| **a > 0.0)
        .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
        .map(|(i, _)| i)?;

    let mut order: Vec<usize> = (0..loops.len()).filter(|&i| areas[i] > 0.0).collect();
    order.sort_by_key(|&i| (i != primary, i));
    let mut holes_of: Vec<Vec<usize>> = vec![Vec::new(); loops.len()];
    for i in 0..loops.len() {
        if areas[i] > 0.0 {
            continue;
        }
        // The hole's containing region: the SMALLEST positive region
        // whose polygon contains a sample point of the hole loop. The
        // twin region (the piece inside the hole) shares the hole's
        // half-edges as twins — its polygon runs THROUGH the sample
        // point and even-odd verdicts there are coin flips, so
        // twin-sharing candidates are excluded (the true container
        // never shares an edge with a hole component of another
        // region's interior). The sample is the midpoint of the
        // hole's LONGEST segment, never a polyline vertex: a vertex
        // is a junction shared by every cell meeting there, and
        // even-odd through it is a coin flip against cells the twin
        // exclusion cannot see (measured, PK-475: the hexagram lip
        // socket's coplanar wall pairs mint residue chords that tile
        // the removed star into hexagon + tip cells; the star-outline
        // hole's first VERTEX lay on the hexagon cell's boundary, the
        // hole attached to that doomed interior cell instead of the
        // face, the kept plane fragment never received its star rim,
        // and all 12 rim record sides stranded). An edge-interior
        // point is shared only by the two sides of that edge — both
        // twin-excluded.
        let p = match polys[i]
            .windows(2)
            .max_by(|a, b| {
                let la = (a[1].0 - a[0].0).powi(2) + (a[1].1 - a[0].1).powi(2);
                let lb = (b[1].0 - b[0].0).powi(2) + (b[1].1 - b[0].1).powi(2);
                la.partial_cmp(&lb).unwrap_or(std::cmp::Ordering::Equal)
            })
            .map(|w| (0.5 * (w[0].0 + w[1].0), 0.5 * (w[0].1 + w[1].1)))
        {
            Some(p) => p,
            None => match polys[i].first() {
                Some(&p) => p,
                None => continue,
            },
        };
        let hole_twins: BTreeSet<usize> =
            loops[i].iter().map(|&h| graph.half_edges[h].twin).collect();
        let mut best: Option<(f64, usize)> = None;
        for (j, &aj) in areas.iter().enumerate() {
            if aj <= 0.0
                || loops[j].iter().any(|h| hole_twins.contains(h))
                || !polyline_contains(&polys[j], p)
            {
                continue;
            }
            if best.is_none_or(|(ba, _)| aj < ba) {
                best = Some((aj, j));
            }
        }
        holes_of[best.map(|(_, j)| j).unwrap_or(primary)].push(i);
    }

    Some(ChartRoleAssembly {
        regions: order
            .into_iter()
            .map(|i| (i, std::mem::take(&mut holes_of[i])))
            .collect(),
    })
}
