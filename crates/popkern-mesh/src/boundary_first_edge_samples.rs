//! Canonical per-edge sampling for the boundary-first tessellator:
//! every edge is sampled ONCE, adaptively by 3D chord height, and
//! both incident faces consume the same polyline (watertight by
//! construction). Split out of `boundary_first.rs` (500-line cap).

use std::collections::HashMap;

use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::CurveEval;
use popkern_core::topology::{BRep, EdgeHandle};

/// Canonical per-edge samples: parameter fractions (0..=1 of the
/// edge's curve window) and their 3D points. Both faces adjacent to
/// the edge consume the same set.
pub(crate) struct EdgeSamples {
    /// Curve-window fractions, ascending, including 0 and 1.
    pub fractions: Vec<f64>,
    /// 3D points at those fractions.
    pub points: Vec<Point3>,
}

/// Sample every edge with a finite curve window once, adaptively by
/// 3D chord height.
pub(crate) fn sample_all_edges(brep: &BRep, chord_height: f64) -> HashMap<EdgeHandle, EdgeSamples> {
    let mut map = HashMap::new();
    for (eh, edge) in brep.edges.iter() {
        let Some(curve) = edge.curve.as_ref() else {
            // Curveless edge: synthesize vertex-chord samples so the
            // ring builder doesn't bail (a boolean can leave
            // degenerate zero-length curveless edges at discarded
            // tool corners; a missing entry made every face touching
            // one fall back to the cracking uniform-grid path). Zero-length
            // edges get NO entry — the ring builder skips them.
            let Some((va, vb)) = edge_endpoint_positions(brep, edge.half_edge) else {
                continue;
            };
            if va.distance_to(&vb) > 1.0e-9 {
                map.insert(
                    eh,
                    EdgeSamples {
                        fractions: vec![0.0, 1.0],
                        points: vec![va, vb],
                    },
                );
            }
            continue;
        };
        if !(curve.t_start.is_finite() && curve.t_end.is_finite()) {
            continue;
        }
        let eval: &dyn CurveEval = &curve.geometry;
        let at = |f: f64| eval.point_at(curve.t_start + (curve.t_end - curve.t_start) * f);
        // Adaptive subdivision by chord height — BREADTH-FIRST under
        // the budget. The old left-to-right scan refined depth-first,
        // so hitting the cap mid-curve left everything right of the
        // cursor UNREFINED: the 6-turn coil seam capped at fraction
        // 0.5 and its ring carried one segment spanning HALF the
        // chart, whose diametral protecting disk evacuated 50% of
        // the interior lattice (measured: the ground spring meshed
        // watertight at 49% of its analytic volume). Splitting the
        // worst-deviation intervals first keeps a capped refinement
        // uniform along the curve.
        // KNOT-SEEDED initial fractions (PK-494): starting from the
        // bare endpoints, the breadth-first sagitta refinement can
        // ALIAS on periodic curves — an interval spanning a whole
        // number of coil turns has its parametric midpoint back in
        // phase with the chord, the deviation test reads ~0, and the
        // interval never splits (measured on a 41-turn crease edge:
        // one surviving segment spanned exactly 20.0 turns — 322
        // metric units of a 776-long chart — and its diametral
        // protecting disk evacuated the ENTIRE interior lattice of
        // the root-band face). Seeding a fraction at every distinct
        // interior knot bounds every initial interval by one knot
        // span (a quarter turn on the canonical helix), where the
        // sagitta test cannot alias. Non-NURBS curves keep the bare
        // endpoints (their sagitta is unimodal per interval).
        let mut fractions = vec![0.0, 1.0];
        if let popkern_core::geometry::CurveGeometry::Nurbs(nc) = &curve.geometry {
            let (t0, t1) = (curve.t_start, curve.t_end);
            let span = (t1 - t0).abs().max(1e-300);
            let mut seeds: Vec<f64> = nc
                .knots
                .iter()
                .map(|&k| (k - t0) / (t1 - t0))
                .filter(|&f| f > 1e-9 && f < 1.0 - 1e-9)
                .collect();
            seeds.dedup_by(|a, b| (*a - *b).abs() * span < 1e-12 * span.max(1.0));
            // Respect the sample budget: subsample evenly if a very
            // fine knot vector alone would exceed it.
            let budget = 4096usize;
            if seeds.len() > budget {
                let step = seeds.len() as f64 / budget as f64;
                seeds = (0..budget)
                    .map(|i| seeds[(i as f64 * step) as usize])
                    .collect();
            }
            fractions.extend(seeds);
            fractions.sort_by(f64::total_cmp);
            fractions.dedup_by(|a, b| (*a - *b).abs() < 1e-12);
            // Alias-proof coarsening: merge seeded intervals back
            // while every DISCARDED seed point stays within the chord
            // tolerance of the merged chord. Judged on the real seed
            // points — not a parametric midpoint — so a merged
            // interval can never hide interior excursions the way an
            // in-phase midpoint does, and a loose tolerance coarsens
            // back to the historical minimal sample sets (the
            // 2-triangle loose-sphere floor).
            if fractions.len() > 2 {
                let pts: Vec<Point3> = fractions.iter().map(|&f| at(f)).collect();
                let dev = |a: usize, b: usize, k: usize| -> f64 {
                    let (pa, pb, pk) = (pts[a], pts[b], pts[k]);
                    let (dx, dy, dz) = (pb.x - pa.x, pb.y - pa.y, pb.z - pa.z);
                    let l2 = dx * dx + dy * dy + dz * dz;
                    if l2 < 1e-300 {
                        return pk.distance_to(&pa);
                    }
                    let t = (((pk.x - pa.x) * dx + (pk.y - pa.y) * dy + (pk.z - pa.z) * dz) / l2)
                        .clamp(0.0, 1.0);
                    let q = Point3::new(pa.x + dx * t, pa.y + dy * t, pa.z + dz * t);
                    pk.distance_to(&q)
                };
                let mut kept: Vec<usize> = vec![0];
                let mut s_idx = 0usize;
                let mut j = 1usize;
                while j < fractions.len() {
                    let mut end = j;
                    while end + 1 < fractions.len() {
                        let candidate = end + 1;
                        let ok = (s_idx + 1..candidate)
                            .all(|k| dev(s_idx, candidate, k) <= chord_height);
                        if ok {
                            end = candidate;
                        } else {
                            break;
                        }
                    }
                    kept.push(end);
                    s_idx = end;
                    j = end + 1;
                }
                fractions = kept.into_iter().map(|i| fractions[i]).collect();
            }
        }
        // Budget and floor scale with the curve's knot-span count
        // (PK-461): a flat 512-sample cap and 1e-3 minimum span
        // under-sampled a 41-turn thread-tool crease edge to ~12
        // samples per coil turn — its sparse ring segments failed the
        // Delaunay coverage next to the per-wrap record chains, every
        // tool wall face fell back to the uniform-grid path, and the
        // served mesh tore along the whole coil (8235 open edges on a
        // body whose topology was correct). Sixteen samples per span
        // (64 per quarter-arc turn) holds sag quality; simple curves
        // converge on chord height far below either bound.
        let spans = match &curve.geometry {
            popkern_core::geometry::CurveGeometry::Nurbs(nc) => nc
                .knots
                .windows(2)
                .filter(|w| w[1] > w[0] + 1e-12)
                .count()
                .max(1),
            _ => 1,
        };
        let cap = (16usize * spans).clamp(512, 8192);
        let min_span = (0.0625 / spans as f64).min(1e-3);
        loop {
            let mut split_at: Vec<(f64, usize)> = Vec::new();
            for i in 0..fractions.len() - 1 {
                let (a, b) = (fractions[i], fractions[i + 1]);
                if (b - a) <= min_span {
                    continue;
                }
                let (pa, pm, pb) = (at(a), at(0.5 * (a + b)), at(b));
                let chord_mid = Point3::new(
                    0.5 * (pa.x + pb.x),
                    0.5 * (pa.y + pb.y),
                    0.5 * (pa.z + pb.z),
                );
                let dev = pm.distance_to(&chord_mid);
                if dev > chord_height {
                    split_at.push((dev, i));
                }
            }
            if split_at.is_empty() || fractions.len() >= cap {
                break;
            }
            let budget = cap - fractions.len();
            if split_at.len() > budget {
                split_at.sort_by(|x, y| y.0.partial_cmp(&x.0).unwrap_or(std::cmp::Ordering::Equal));
                split_at.truncate(budget);
                split_at.sort_by_key(|&(_, i)| i);
            }
            for &(_, i) in split_at.iter().rev() {
                let mid = 0.5 * (fractions[i] + fractions[i + 1]);
                fractions.insert(i + 1, mid);
            }
        }
        let mut points: Vec<Point3> = fractions.iter().map(|&f| at(f)).collect();
        // Snap the endpoint samples to the shared VERTEX positions:
        // curve-endpoint evaluations can differ from the topological
        // vertex by ~1e-8 (refit/snap noise), and two edges meeting
        // at a vertex would then produce near-miss boundary points
        // that straddle the weld rounding. The gate is the body's own
        // linear tolerance (floored at 1e-6): a curve end within the
        // tolerance of its vertex IS that vertex by the body's
        // contract (a bored sphere at the 1e-3 tolerance carries its
        // exact seam arc's end 1e-4 from the vertex the 1e-3
        // intersection placed; an absolute 1e-6 gate left the seam
        // and rim rings disagreeing there and 16 open mesh edges).
        let snap_gate = brep.tolerance.linear().max(1.0e-6);
        if let (Some(first), Some(last), Some((va, vb))) = (
            points.first().copied(),
            points.last().copied(),
            edge_endpoint_positions(brep, edge.half_edge),
        ) {
            let (pa, pb) = if first.distance_to(&va) + last.distance_to(&vb)
                <= first.distance_to(&vb) + last.distance_to(&va)
            {
                (va, vb)
            } else {
                (vb, va)
            };
            if first.distance_to(&pa) < snap_gate {
                points[0] = pa;
            }
            let n = points.len() - 1;
            if points[n].distance_to(&pb) < snap_gate {
                points[n] = pb;
            }
        }
        map.insert(eh, EdgeSamples { fractions, points });
    }
    map
}

/// Endpoint vertex positions of an edge via its canonical half-edge.
///
/// Boundary edges have a dead twin; the end vertex then comes from the
/// surviving side's `next` origin. Returns `None` when neither
/// resolution succeeds (stale canonical half-edge) — callers skip the
/// edge rather than panic on partially-retired topology.
pub(crate) fn edge_endpoint_positions(
    brep: &BRep,
    half_edge: popkern_core::topology::HalfEdgeHandle,
) -> Option<(Point3, Point3)> {
    let he = brep.half_edges.get(half_edge)?;
    let end_vertex = match brep.half_edges.get(he.twin) {
        Some(twin) => twin.origin,
        None => brep.half_edges.get(he.next)?.origin,
    };
    Some((
        brep.vertices.get(he.origin)?.position,
        brep.vertices.get(end_vertex)?.position,
    ))
}

/// Env-gated Delaunay validity audit (the delaunay-audit trace topic): count
/// circumcircle violations in a Bowyer-Watson output. Proved the
/// PK-146 seam-band fans were LEGAL (zero violations) — the defect
/// was the cocircular lattice, not the triangulation.
pub(crate) fn delaunay_audit(
    face_handle: popkern_core::topology::FaceHandle,
    bw_pts: &[(f64, f64)],
    tris: &[[usize; 3]],
) {
    if !popkern_core::trace::trace_enabled("delaunay-audit") || bw_pts.len() <= 500 {
        return;
    }
    {
        let mut violations = 0usize;
        let mut worst = 0.0_f64;
        for t in tris {
            let (a, b, c) = (bw_pts[t[0]], bw_pts[t[1]], bw_pts[t[2]]);
            let d2 = 2.0 * (a.0 * (b.1 - c.1) + b.0 * (c.1 - a.1) + c.0 * (a.1 - b.1));
            if d2.abs() < 1e-30 {
                continue;
            }
            let ux = ((a.0 * a.0 + a.1 * a.1) * (b.1 - c.1)
                + (b.0 * b.0 + b.1 * b.1) * (c.1 - a.1)
                + (c.0 * c.0 + c.1 * c.1) * (a.1 - b.1))
                / d2;
            let uy = ((a.0 * a.0 + a.1 * a.1) * (c.0 - b.0)
                + (b.0 * b.0 + b.1 * b.1) * (a.0 - c.0)
                + (c.0 * c.0 + c.1 * c.1) * (b.0 - a.0))
                / d2;
            let r2 = (a.0 - ux).powi(2) + (a.1 - uy).powi(2);
            for (k, p) in bw_pts.iter().enumerate() {
                if k == t[0] || k == t[1] || k == t[2] {
                    continue;
                }
                let dd = (p.0 - ux).powi(2) + (p.1 - uy).powi(2);
                if dd < r2 * (1.0 - 1e-9) {
                    violations += 1;
                    let depth = (r2.sqrt() - dd.sqrt()) / r2.sqrt();
                    if depth > worst {
                        worst = depth;
                    }
                    break;
                }
            }
        }
        eprintln!(
            "  [bf-delaunay-audit] face {face_handle:?}: {} tris, {violations} circumcircle violations (worst depth {worst:.3})",
            tris.len()
        );
    }
}

/// Whether a half-edge's pcurve runs ALONG the ring-walk direction —
/// the walk mirrors the pcurve window when it does not. Detected from
/// the pcurve endpoint images against the walk origin, EXCEPT closed
/// edges on APERIODIC charts (closed NURBS walls), which walk IN
/// pcurve direction by construction: the endpoint test is 0-vs-0
/// noise there (both ends image to the same 3D point) and a flip
/// lands the seam junction on the wrong geometric-closure copy,
/// which no ring unwrap repairs (PK-449). Periodic charts keep the
/// endpoint test: their ring-level unwrap absorbs copy jumps, and
/// forcing pcurve direction turned a cap-cut sphere's rim into a
/// phantom full wrap.
pub(super) fn pcurve_forward_along_walk(
    pc: &popkern_core::geometry::pcurve::Pcurve,
    surface: &popkern_core::geometry::SurfaceGeometry,
    origin: &Point3,
    closed_edge: bool,
) -> bool {
    use popkern_core::geometry::surface::SurfaceEval as _;
    if closed_edge && surface.parameter_periodicity() == (None, None) {
        return true;
    }
    let p_start = pc.point_at(pc.t_start);
    let p_end = pc.point_at(pc.t_end);
    let s_start = surface.point_at(p_start.x, p_start.y);
    let s_end = surface.point_at(p_end.x, p_end.y);
    s_start.distance_to(origin) <= s_end.distance_to(origin)
}

/// Start the ring walk at a NON-closed edge when the loop has one: a
/// closed edge's traversal branch is decided by continuity with the
/// previous exit uv (see [`ring_walk_direction`]), which needs a
/// predecessor — a loop-opening closed edge falls back to its
/// authored pcurve direction, which is phase-arbitrary for boolean
/// records (PK-454: the z = 0.45 band grind's loop began at the cut
/// oval and the ring entered it on the wrong seam branch).
pub(super) fn rotate_ring_start(
    brep: &BRep,
    start: popkern_core::topology::HalfEdgeHandle,
) -> popkern_core::topology::HalfEdgeHandle {
    let mut cand = start;
    for _ in 0..10_000 {
        // GEOMETRIC closure, not vertex-handle identity: a tube seam
        // can carry two distinct vertex handles at one position.
        if edge_endpoint_positions(brep, cand).is_some_and(|(a, b)| a.distance_to(&b) >= 1.0e-9) {
            return cand;
        }
        cand = brep.half_edges[cand].next;
        if cand == start {
            break;
        }
    }
    start
}

/// Ring-walk direction decision for one half-edge: returns
/// `(closed_edge, closed_entry_forward, forward)`.
///
/// For a CLOSED edge on an APERIODIC chart the pcurve's authored
/// direction is PHASE-ARBITRARY (a boolean record's fitted pcurve
/// runs whichever way the march seeded), and a seam-crossing edge has
/// its two pcurve ends on the two closure branches: entering on the
/// wrong one span-jumps the ring across the whole chart, the jump
/// segment is uncovered, and encroachment splits churn without
/// converging (PK-454: 1100 open mesh edges on an OCC-valid body).
/// The ENTRY end is the one branch-continuous with the previous ring
/// exit uv; the sample-walk direction is then DERIVED — never probed
/// at the closure point, where every distance is knife-edge — from
/// the pcurve-vs-curve alignment measured at 25% along:
/// `forward = entry == along`. Open edges keep the origin-distance
/// test; closed edges on periodic charts keep the early-image test
/// (their ring-level unwrap absorbs closure-copy jumps — PK-449).
pub(super) fn ring_walk_direction(
    pc: &popkern_core::geometry::pcurve::Pcurve,
    surface: &popkern_core::geometry::SurfaceGeometry,
    samples: &EdgeSamples,
    uv_last: Option<(f64, f64)>,
    origin: &Point3,
) -> Option<(bool, Option<bool>, bool)> {
    use popkern_core::geometry::surface::SurfaceEval as _;
    let closed_edge = samples
        .points
        .first()
        .zip(samples.points.last())
        .map(|(a, b)| a.distance_to(b) < 1e-9)
        .unwrap_or(false);
    let aperiodic = surface.parameter_periodicity() == (None, None);
    // `uv_last` is the previous ring exit ONLY on mixed loops (the
    // caller passes None for all-closed loops, whose pcurve family is
    // authored consistently and keeps the historical rules).
    let closed_entry_forward = if closed_edge && aperiodic {
        uv_last.map(|(lu, lv)| {
            let (pt0, pt1) = pc.parameter_range();
            let a = pc.point_at(pt0);
            let b = pc.point_at(pt1);
            let da = (a.x - lu).powi(2) + (a.y - lv).powi(2);
            let db = (b.x - lu).powi(2) + (b.y - lv).powi(2);
            da <= db
        })
    } else {
        None
    };
    let forward = if let (Some(entry_forward), true) =
        (closed_entry_forward, samples.points.len() >= 5)
    {
        let (pt0, pt1) = pc.parameter_range();
        let q_uv = pc.point_at(pt0 + (pt1 - pt0) * 0.25);
        let q = surface.point_at(q_uv.x, q_uv.y);
        let n = samples.points.len();
        let along =
            q.distance_to(&samples.points[n / 4]) <= q.distance_to(&samples.points[n - 1 - n / 4]);
        entry_forward == along
    } else if closed_edge && samples.points.len() >= 3 {
        let (pt0, pt1) = pc.parameter_range();
        let early_uv = pc.point_at(pt0 + (pt1 - pt0) * 0.05);
        let early = surface.point_at(early_uv.x, early_uv.y);
        early.distance_to(&samples.points[1])
            <= early.distance_to(&samples.points[samples.points.len() - 2])
    } else {
        samples.points.first()?.distance_to(origin) <= samples.points.last()?.distance_to(origin)
    };
    Some((closed_edge, closed_entry_forward, forward))
}

/// Refine a boundary sample's UV from its pcurve seed by point
/// inversion of the canonical 3D sample, returning
/// `(u, v, residual, gate)`. The pcurve is only a SEED — imprecise
/// pcurves (straight chords for bulging curves) still yield exact
/// boundary UVs. Chain, best image wins at every rung:
/// - analytic cylinders invert by atan2 in the surface's own basis
///   (the generic Newton wanders period branches);
/// - the seeded projection is kept only when it lands CLOSER than the
///   seed's own image — near a geometrically-closed chart's seam it
///   can stall clamped on the border (PK-454: a fit-quality seed at
///   u = 0.9839 "refined" to the u = 1.0 border 3.1e-3 off-surface);
/// - a COLD global projection retries under the same gate (a setback
///   fan's stale rim pcurves sat 1.8e-3 off and the seeded projection
///   converged near the displaced location);
/// - the knot-seeded MULTI-SHEET projector is the last resort: on a
///   many-turn closed NURBS wall the seeded projection stalls on the
///   border and the cold 5-start walks the ADJACENT COIL TURN, so
///   both read millimetres for a point ON the surface (PK-454, 1023
///   open mesh edges on an OCC-valid body).
///
/// The residual gate itself must NOT be widened by recorded entity
/// tolerance — tried twice, falsified twice: it repartitions the
/// boundary-first/uniform-grid face split and the mixed boundaries
/// crack far worse (hexagon6_blunt 0 → 482).
pub(super) fn refine_boundary_uv(
    surface: &popkern_core::geometry::SurfaceGeometry,
    seed: (f64, f64),
    target: &Point3,
    samples: &EdgeSamples,
) -> Option<(f64, f64, f64, f64)> {
    use popkern_core::geometry::surface::SurfaceEval as _;
    let (mut uu, mut vv) = seed;
    if let popkern_core::geometry::SurfaceGeometry::Cylinder { origin, axis, .. } = surface {
        let (perp1, perp2) = popkern_core::geometry::analytic_surfaces::orthonormal_basis(*axis);
        let d = *target - *origin;
        (uu, vv) = (d.dot(&perp2).atan2(d.dot(&perp1)), d.dot(axis));
    } else if surface.point_at(uu, vv).distance_to(target) > 1.0e-9
        && let Some((ru, rv)) =
            popkern_core::geometry::uv_projection::project_point_to_uv_near(surface, target, uu, vv)
    {
        let seed_res = surface.point_at(uu, vv).distance_to(target);
        let refined_res = surface.point_at(ru, rv).distance_to(target);
        if refined_res < seed_res {
            (uu, vv) = (ru, rv);
        }
    }
    let mut residual = surface.point_at(uu, vv).distance_to(target);
    let scale = samples
        .points
        .first()?
        .distance_to(samples.points.last()?)
        .max(1e-9);
    let gate = 1.0e-3 * scale.max(1.0);
    if residual > gate {
        if let Some((gu, gv)) =
            popkern_core::geometry::uv_projection::project_point_to_uv(surface, target)
        {
            let g_res = surface.point_at(gu, gv).distance_to(target);
            if g_res < residual {
                (uu, vv) = (gu, gv);
                residual = g_res;
            }
        }
        // The rescue must stay on the RING'S BRANCH: the multi-sheet
        // projector returns the best image on any closure copy, and a
        // period-shifted uv splices a foreign branch into the ring (a
        // pipe tube's rescued samples tore 127 open edges where the
        // historical gate failure fell back to the uniform grid and
        // meshed watertight). Accept only an image within half a
        // closure period of the pcurve seed per closed axis; anything
        // farther keeps the gate failure — the tuned fallback.
        if residual > gate
            && let Some((su, sv)) =
                popkern_core::geometry::uv_projection::project_point_to_uv_multisheet(
                    surface, target,
                )
        {
            let (pu, pv) = popkern_core::geometry::chart_periods::geometric_chart_periods(surface);
            let on_seed_branch = |x: f64, sx: f64, p: Option<f64>| match p {
                Some(p) if p.is_finite() && p > 0.0 => (x - sx).abs() <= 0.5 * p,
                _ => true,
            };
            let s_res = surface.point_at(su, sv).distance_to(target);
            if s_res < residual && on_seed_branch(su, seed.0, pu) && on_seed_branch(sv, seed.1, pv)
            {
                (uu, vv) = (su, sv);
                residual = s_res;
            }
        }
    }
    Some((uu, vv, residual, gate))
}
