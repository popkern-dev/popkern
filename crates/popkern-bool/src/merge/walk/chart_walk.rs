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

/// The rotational successor of `h` on this face's chart.
///
/// Returns `None` only on substrate defects (missing branch point or
/// no departures) — impossible on a subdivision that passed I1–I5,
/// and treated as a face-level fallback trigger by callers.
pub(crate) fn next_in_chart(
    sub: &FaceChartSubdivision,
    graph: &IntersectionGraph,
    h: usize,
) -> Option<usize> {
    let bp = sub.end_bp.get(&h)?;
    let deps = sub.departures.get(bp)?;
    let twin = graph.half_edges[h].twin;

    // POLE branch points rotate by AZIMUTH: every incident half-edge
    // lands at the pole carrying its meridian's u, and the geometric
    // angle around the pole IS that u — local UV directions there all
    // point along ±v and cannot order the fan. The material-left
    // successor is the departure with the smallest positive
    // (u_arrive − u_depart)·dir mod period, dir folding σ with the
    // pole side (rotation sense reverses between poles); the twin
    // sits at the arrival azimuth (angle → period) and ranks last
    // naturally.
    if let Some(&dir) = sub.poles.get(&bp.node)
        && let Some(period) = sub.u_period
    {
        let emb_h = sub.embeddings.get(&h)?;
        let u_a = emb_h.uv[emb_h.uv.len() - 1].0;
        let mut best: Option<(f64, usize)> = None;
        for &c in deps {
            let Some(emb_c) = sub.embeddings.get(&c) else {
                continue;
            };
            let u_c = emb_c.uv[0].0;
            let mut a = ((u_a - u_c) * dir).rem_euclid(period);
            if a <= 1e-12 || c == twin {
                a += period;
            }
            let better = match best {
                None => true,
                Some((ba, bc)) => a < ba - 1e-15 || ((a - ba).abs() <= 1e-15 && c < bc),
            };
            if better {
                best = Some((a, c));
            }
        }
        return best.map(|(_, c)| c);
    }

    // GLOBAL WEDGE PAIRING (permutation by construction): the former
    // per-arrival min-CW rule picked each arrival's successor
    // independently, and at a slop-cluster junction two arrivals
    // could elect the SAME departure — the "orbits" then stop being a
    // permutation, the walk collides, and the face falls back to the
    // 3D-Newell walker mis-partitioned (measured: rotated torus × torus,
    // a node where the correct matching requires one arrival to
    // U-turn; the mis-partition cost 5.5% of the union volume before
    // the strand gate declined it). The DCEL-correct assignment is
    // the cyclic wedge matching over ALL slots at the branch point:
    // arrivals (at their back directions) and departures (at their
    // start directions) sort into one σ-oriented circular order and
    // match parenthesis-style sweeping clockwise — each departure
    // pairs the nearest pending arrival, injective by construction.
    let bp_key = *bp;
    let psi = |d: (f64, f64)| -> f64 {
        (sub.sigma * d.1)
            .atan2(d.0)
            .rem_euclid(std::f64::consts::TAU)
    };
    // Slots: (angle, is_departure, he). Tie order in the clockwise
    // sweep puts the DEPARTURE first so an exactly-parallel departure
    // is NOT taken by its parallel arrival (the historic strictly-
    // positive-angle rule; the twin's departure is exactly parallel
    // to its arrival's back and ranks last this way).
    let mut slots: Vec<(f64, bool, usize)> = Vec::new();
    for (hh, e) in sub.end_bp.iter() {
        if *e == bp_key
            && let Some(d) = sub.embeddings.get(hh).and_then(|em| em.end_dir())
        {
            slots.push((psi((-d.0, -d.1)), false, *hh));
        }
    }
    for &c in deps {
        if let Some(d) = sub.embeddings.get(&c).and_then(|em| em.start_dir()) {
            slots.push((psi(d), true, c));
        }
    }
    let _ = twin;
    // Clockwise = descending ψ. Deterministic ties: departures first,
    // then half-edge index.
    slots.sort_by(|a, b| {
        b.0.partial_cmp(&a.0)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| b.1.cmp(&a.1))
            .then_with(|| a.2.cmp(&b.2))
    });
    // Two clockwise passes with a pending-arrival stack: a departure
    // matches the most recent pending arrival (wedge nesting); the
    // second pass closes pairs straddling the angular origin.
    if crate::trace_enabled_wedge() {
        let names: Vec<String> = slots
            .iter()
            .map(|(a, d, hh)| format!("{}{}@{a:.3}", if *d { "dep" } else { "arr" }, hh))
            .collect();
        eprintln!("      [wedge] bp {bp_key:?} slots {}", names.join(" "));
    }
    let mut pending: Vec<usize> = Vec::new();
    let mut matched: std::collections::HashMap<usize, usize> = std::collections::HashMap::new();
    let mut used_deps: std::collections::HashSet<usize> = std::collections::HashSet::new();
    for _pass in 0..2 {
        for &(_, is_dep, hh) in &slots {
            if is_dep {
                if used_deps.contains(&hh) {
                    continue;
                }
                if let Some(a) = pending.pop() {
                    matched.insert(a, hh);
                    used_deps.insert(hh);
                }
            } else if !matched.contains_key(&hh) && !pending.contains(&hh) {
                pending.push(hh);
            }
        }
    }
    // Exact-tie degeneracy fallback: pair leftovers deterministically.
    let mut leftover_deps: Vec<usize> = deps
        .iter()
        .copied()
        .filter(|d| !used_deps.contains(d))
        .collect();
    leftover_deps.sort_unstable();
    pending.sort_unstable();
    for a in pending {
        if matched.contains_key(&a) {
            continue;
        }
        if let Some(d) = leftover_deps.pop() {
            matched.insert(a, d);
        }
    }
    // NEAR-TIE TWIN-MATCH REPAIR: when a record chain meets a face
    // boundary almost tangentially (a thread groove exiting through
    // a rim at an x.75 phase meets it 0.013 rad off the tangent),
    // the slot angles of the record pair and the boundary departure
    // tie within the pcurve-fit slop, the sweep order flips, and the
    // parenthesis matching pairs the record arrival with its OWN
    // TWIN — the walk then traverses the chain out-and-back as a
    // zero-area ring and the region beside it never gets a fragment
    // (measured: the thread cap's kept crescent vanished and the
    // whole cap classified away, stranding the mouth records —
    // PK-460). A twin pairing inside a near-tie cluster is the
    // degenerate signature; swapping with the tied pair restores
    // both rings. Far-from-tie twin pairings are LEGITIMATE U-turns
    // (the rotated torus × torus junction this matching was built
    // for) and are left alone.
    const NEAR_TIE: f64 = 0.02;
    let angle_of = |he: usize, is_dep: bool| -> Option<f64> {
        slots
            .iter()
            .find(|(_, d, hh)| *d == is_dep && *hh == he)
            .map(|(a, _, _)| *a)
    };
    let circ_close = |a: f64, b: f64| -> bool {
        let d = (a - b).abs();
        d.min(std::f64::consts::TAU - d) < NEAR_TIE
    };
    let twin_matched: Vec<(usize, usize)> = matched
        .iter()
        .filter(|(a, d)| graph.half_edges[**a].twin == **d)
        .map(|(a, d)| (*a, *d))
        .collect();
    for (a, d) in twin_matched {
        let (Some(a_ang), Some(d_ang)) = (angle_of(a, false), angle_of(d, true)) else {
            continue;
        };
        let swap = matched
            .iter()
            .filter(|(a2, _)| **a2 != a)
            .filter(|(a2, d2)| {
                let close = angle_of(**d2, true).is_some_and(|x| circ_close(x, d_ang))
                    || angle_of(**a2, false).is_some_and(|x| circ_close(x, a_ang));
                close && graph.half_edges[a].twin != **d2 && graph.half_edges[**a2].twin != d
            })
            .map(|(a2, d2)| (*a2, *d2))
            .min_by_key(|&(a2, _)| a2);
        if let Some((a2, d2)) = swap {
            matched.insert(a, d2);
            matched.insert(a2, d);
            if crate::trace_enabled_wedge() {
                eprintln!("      [wedge] near-tie twin repair at {bp_key:?}: {a}->{d2}, {a2}->{d}");
            }
        }
    }
    if crate::trace_enabled_wedge() {
        eprintln!("      [wedge] h {h} -> {:?}", matched.get(&h));
    }
    matched.get(&h).copied()
}

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
    let source_brep = match face_ref.operand {
        super::super::graph::Operand::A => brep_a,
        super::super::graph::Operand::B => brep_b,
    };
    let mut edges = std::collections::BTreeSet::new();
    if let Some(face) = source_brep.faces.get(face_ref.face) {
        for lh in face.inner_loops.iter().copied() {
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
        // whose polygon contains the hole's first point. The twin
        // region (the piece inside the hole) shares the hole's
        // half-edges as twins — its polygon runs THROUGH the sample
        // point and even-odd verdicts there are coin flips, so
        // twin-sharing candidates are excluded (the true container
        // never shares an edge with a hole component of another
        // region's interior).
        let Some(&p) = polys[i].first() else {
            continue;
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
