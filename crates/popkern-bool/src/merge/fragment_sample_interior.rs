//! Interior medial sampling for all-shared fragments (extracted from
//! `fragment_sample.rs` under the 500-line file budget). See
//! [`interior_step_sample`].

use popkern_core::{geometry::Point3, topology::BRep};

use super::graph::{FaceRef, IntersectionGraph, Operand};

/// Sample the fragment interior by stepping off the longest outer
/// edge's UV midpoint perpendicular to the loop, on whichever side the
/// outer-loop polygon contains, then evaluating the face surface.
/// `None` when pcurves/surface are unavailable or neither side is
/// contained (seam-wrapped charts) — the caller falls back.
pub(crate) fn interior_step_sample(
    graph: &IntersectionGraph,
    loop_hes: &[usize],
    face: &FaceRef,
    brep_a: &BRep,
    brep_b: &BRep,
) -> Option<Point3> {
    use popkern_core::geometry::surface::SurfaceEval;
    let brep = match face.operand {
        Operand::A => brep_a,
        Operand::B => brep_b,
    };
    let surface = &brep.faces.get(face.face)?.surface.as_ref()?.geometry;

    // The pcurve in THIS face's chart, with the traversal direction of
    // half-edge `h`. Pairs store the face_a chart on the forward half
    // and the face_b chart on the twin, so `h.pcurve` alone mixes
    // charts on cross-registered loops: an all-shared wall strip's
    // longest edge carried the ARM-plane chart, the step landed at a
    // wall angle far outside the strip, and the doomed strip
    // classified Outside and served (PK-439, flush-cap prism ∪ hub
    // off-seam). `reversed` = the chart lives on the twin, so its
    // samples must run t1→t0 to follow `h`.
    let pcurve_for = |h: usize| -> Option<(&popkern_core::geometry::pcurve::Pcurve, bool)> {
        let he = &graph.half_edges[h];
        if he.face.as_ref() == Some(face) {
            return he.pcurve.as_ref().map(|p| (p, false));
        }
        let tw = &graph.half_edges[he.twin];
        if tw.face.as_ref() == Some(face) {
            return tw.pcurve.as_ref().map(|p| (p, true));
        }
        he.pcurve.as_ref().map(|p| (p, false))
    };

    let trace = popkern_core::trace::trace_enabled("loop-sampling");
    // UV run of one half-edge in THIS face's chart, following the
    // half-edge's traversal direction: the face-owned pcurve when
    // present, else the walk's chart embedding (which projects the
    // trimmed 3D curve). Split record pieces carry NO pcurves by
    // design — graph_split strips the stale full-span charts — and
    // declining on them dropped the whole interior step (measured: a
    // spring-grind section disc fell to the knife-edge centroid ON
    // its own chord, ray-classified Outside, and the subtraction
    // stranded the wall record's other side — PK-443).
    let chart_run = |h: usize, n: usize| -> Option<Vec<(f64, f64)>> {
        if let Some((pc, reversed)) = pcurve_for(h) {
            let (t0, t1) = pc.parameter_range();
            let (s0, s1) = if reversed { (t1, t0) } else { (t0, t1) };
            let mut run: Vec<(f64, f64)> = (0..n)
                .map(|i| {
                    let t = s0 + (s1 - s0) * (i as f64) / n as f64;
                    let q = pc.point_at(t);
                    (q.x, q.y)
                })
                .collect();
            // Orient by the GRAPH NODES: registration direction does
            // not guarantee a stored pcurve runs with its half-edge
            // (split shared pairs violate it — the chart embedding
            // measures orientation the same way). Trusting the
            // registration flag ran a banana fragment's third flank
            // BACKWARDS, the chart polygon bowtied, the
            // self-intersection gate declined, and the knife-edge
            // centroid fell in the banana's concave hollow — off the
            // solid — so the tip's cap disc classified Outside and
            // the subtraction stranded (PK-453, the z = 0.75 deep
            // grind whose plate top crosses the coil's cap rim).
            // Closed self-loops (origin == end) have no node
            // authority and keep the registration direction.
            let he = &graph.half_edges[h];
            let end_node = graph.half_edges[he.twin].origin;
            if he.origin != end_node && run.len() >= 2 {
                let p_start = graph.nodes[he.origin].position;
                let p_end = graph.nodes[end_node].position;
                let first = surface.point_at(run[0].0, run[0].1);
                let last = surface.point_at(run[run.len() - 1].0, run[run.len() - 1].1);
                let as_is = first.distance_to(&p_start) + last.distance_to(&p_end);
                let flipped = first.distance_to(&p_end) + last.distance_to(&p_start);
                if flipped < as_is {
                    run.reverse();
                }
            }
            return Some(run);
        }
        let emb = super::walk::chart_embedding::build_embedding(graph, face, h, surface)?;
        // The embedding's own density serves; drop its final sample so
        // concatenated runs do not duplicate junction points.
        let mut uv = emb.uv;
        if uv.len() > 1 {
            uv.pop();
        }
        Some(uv)
    };

    // UV polygon of the outer loop. Density adapts to the loop's
    // half-edge count: a whole closed intersection loop can arrive as
    // ONE self-loop half-edge, and an 8-gon of a girdle-shaped curve
    // under-covers it so badly that the 5%-of-diagonal interior step
    // lands outside the polygon on both sides (PK-440: the cylinder
    // disc inside a propeller blade fell to the knife-edge centroid
    // and classified Outside — the phantom disc served).
    let per_he = (64 / loop_hes.len().max(1)).clamp(8, 64);
    let mut runs: Vec<(usize, Vec<(f64, f64)>)> = Vec::with_capacity(loop_hes.len());
    for &h in loop_hes {
        let Some(run) = chart_run(h, per_he) else {
            if trace {
                eprintln!(
                    "  [sample] interior step decline: he {h} has no chart run in {face:?} \
                     (he.face {:?}, twin.face {:?})",
                    graph.half_edges[h].face, graph.half_edges[graph.half_edges[h].twin].face
                );
            }
            return None;
        };
        runs.push((h, run));
    }
    // Chart-continuous runs: each run shifts by whole periods to start
    // where the previous one ended. A record's stored pcurve can sit a
    // period off its neighbours (measured, PK-922: a drill hole's third
    // record ran at u = 2pi on a band of [-pi/2, pi/2]); the polygon
    // then read as wrapping, the medial search was skipped, and the
    // hole disc sampled outside the drill and was kept.
    let (run_pu, run_pv) = crate::clip::classification_chart_periods(surface);
    for i in 1..runs.len() {
        let (Some(&end), Some(&start)) = (runs[i - 1].1.last(), runs[i].1.first()) else {
            continue;
        };
        let lift = |gap: f64, period: Option<f64>| match period {
            Some(p) if p.is_finite() && p > 0.0 => (gap / p).round() * p,
            _ => 0.0,
        };
        let (du, dv) = (lift(end.0 - start.0, run_pu), lift(end.1 - start.1, run_pv));
        if du != 0.0 || dv != 0.0 {
            for q in &mut runs[i].1 {
                q.0 += du;
                q.1 += dv;
            }
        }
    }
    let polygon: Vec<(f64, f64)> = runs.iter().flat_map(|(_, r)| r.iter().copied()).collect();
    if polygon.len() < 3 {
        return None;
    }
    // A SELF-INTERSECTING chart polygon means the chart embedding is
    // not a faithful picture of the region (a sphere corner circuit's
    // azimuth jumps near the pole fold the polygon over itself), and
    // every even-odd containment answer on it is luck: the marched
    // sample landed on the far side of the sphere, the inside-cube
    // triangle read Outside, and the corner union served a WRONG body
    // (volume 0.82 vs 2.26 — PK-441). Decline to the centroid
    // fallback instead; a wrong-but-knife-edge centroid strands the
    // assembly and the boolean declines honestly.
    if let Some((i, j)) = polygon_self_intersection(&polygon) {
        if trace {
            let n = polygon.len();
            let lobe = &polygon[(i + 1).min(n - 1)..=j.min(n - 1)];
            let (mut lo, mut hi) = ((f64::MAX, f64::MAX), (f64::MIN, f64::MIN));
            for &(u, v) in lobe {
                lo = (lo.0.min(u), lo.1.min(v));
                hi = (hi.0.max(u), hi.1.max(v));
            }
            eprintln!(
                "  [sample] interior step decline: self-intersecting {n}-pt chart polygon of {face:?} (segs {i}x{j}, lobe {} pts, extent ({:.5},{:.5}))",
                lobe.len(),
                hi.0 - lo.0,
                hi.1 - lo.1
            );
        }
        return None;
    }
    let (mut lo_u, mut hi_u, mut lo_v, mut hi_v) = (f64::MAX, f64::MIN, f64::MAX, f64::MIN);
    for &(u, v) in &polygon {
        lo_u = lo_u.min(u);
        hi_u = hi_u.max(u);
        lo_v = lo_v.min(v);
        hi_v = hi_v.max(v);
    }
    let diag = ((hi_u - lo_u).powi(2) + (hi_v - lo_v).powi(2)).sqrt();
    if diag <= 0.0 {
        return None;
    }
    let delta = diag * 0.05;

    // Longest outer edge by 3D arc length of its run, whose run can
    // supply a midpoint and tangent. Endpoint-chord spans tie between
    // a curve and the chord subtending it (the double-grind cut cap's
    // wall record wraps the tube between the same two junctions its
    // disc chord connects — both chords 0.0794): the tie picked the
    // CHORD, whose midpoint lies exactly ON the region boundary, the
    // step refined down onto the knife edge, and the cap classified
    // Outside off its own disc plane (PK-441). Arc length puts the
    // step on the wrapping record's far side, robustly interior.
    // Measured in 3D, not chart units: chart metric distortion near a
    // sphere pole inflates azimuth runs that are millimetres long.
    let run_len = |r: &[(f64, f64)]| -> f64 {
        r.windows(2)
            .map(|w| {
                surface
                    .point_at(w[0].0, w[0].1)
                    .distance_to(&surface.point_at(w[1].0, w[1].1))
            })
            .sum()
    };
    let (_, run) = runs
        .iter()
        .filter(|(_, r)| r.len() >= 2)
        .max_by(|(x, rx), (y, ry)| {
            run_len(rx)
                .partial_cmp(&run_len(ry))
                .unwrap_or(std::cmp::Ordering::Equal)
                .then_with(|| y.cmp(x))
        })?;
    // Base stations along the run: a region bounded by a curve that
    // DOUBLES BACK (the double-grind crescent is fenced by one wall
    // record whose two flanks meet at a degenerate tip) has its run
    // midpoint AT the tip, where every inward step is knife-edge thin
    // (PK-441). Several stations compete below; the widest medial
    // span wins.
    let station = |idx: usize| -> Option<((f64, f64), (f64, f64))> {
        let idx = idx.min(run.len() - 1);
        let q = run[idx];
        let q_ahead = run[(idx + 1).min(run.len() - 1)];
        let q_behind = run[idx.saturating_sub(1)];
        let (tx, ty) = (q_ahead.0 - q_behind.0, q_ahead.1 - q_behind.1);
        let t_len = (tx * tx + ty * ty).sqrt();
        if t_len <= 0.0 {
            return None;
        }
        Some((q, (-ty / t_len, tx / t_len))) // normal = left of travel
    };
    let n_run = run.len();
    let mid = n_run / 2;
    let (q, (px, py)) = station(mid)?;
    let contains = |u: f64, v: f64| -> bool {
        let mut inside = false;
        let n = polygon.len();
        let mut j = n - 1;
        for i in 0..n {
            let (ui, vi) = polygon[i];
            let (uj, vj) = polygon[j];
            if ((vi > v) != (vj > v)) && (u < (uj - ui) * (v - vi) / (vj - vi) + ui) {
                inside = !inside;
            }
            j = i;
        }
        inside
    };
    // Medial sample: march inward along the step normal, find the
    // longest contiguous contained span, and sample its MIDDLE. The
    // first-contained point sits one step off the boundary — for the
    // double-grind cut crescent that shallow sample lay 1.8e-3 inside
    // the coil wall and the point-in-solid ray read it Outside, so the
    // kept cap fragment was discarded and the assembly stranded both
    // record sides (PK-441). The span middle is the deepest sample
    // this direction offers.
    // The medial search below trusts even-odd polygon containment,
    // which is meaningless when the polygon WRAPS a periodic chart
    // direction (a sphere's corner-triangle records wrap the pole
    // azimuth: the polygon "contained" a point on the OPPOSITE side
    // of the sphere, the inside-cube triangle classified Outside,
    // and the corner-sphere union served a wrong body). On wrapping
    // charts only the near-step delta ladder below runs — its small
    // steps stay in the neighbourhood where even-odd noise cannot
    // reach the far side.
    // Chart-aware periods: geometric closure of a NURBS chart counts
    // (a closed hub wall wraps exactly like an analytic cylinder).
    let (pp_u, pp_v) = crate::clip::classification_chart_periods(surface);
    let wraps = |span: f64, period: Option<f64>| {
        period.is_some_and(|p| p.is_finite() && p > 0.0 && span > 0.5 * p)
    };
    let chart_wraps = wraps(hi_u - lo_u, pp_u) || wraps(hi_v - lo_v, pp_v);

    let march = diag * 0.01;
    // (span, base, normal, sign, s_lo, s_hi)
    type MedialSpan = (f64, (f64, f64), (f64, f64), f64, f64, f64);
    let mut best_span: Option<MedialSpan> = None;
    let stations: &[usize] = if chart_wraps {
        &[]
    } else {
        &[mid, n_run / 4, (3 * n_run) / 4, n_run / 8, (7 * n_run) / 8]
    };
    for &idx in stations {
        let Some((b, (nx, ny))) = station(idx) else {
            continue;
        };
        for sign in [1.0_f64, -1.0] {
            let mut lo: Option<f64> = None;
            let mut hi = 0.0_f64;
            // The base sits ON the boundary, so a genuine inward
            // direction is contained from the very first steps. A
            // span that only begins farther out is NOT this
            // boundary's inward side — on a chart-degenerate polygon
            // (sphere records wrapping a pole) the even-odd test
            // "contains" arbitrary far regions, and accepting one
            // sampled the OPPOSITE side of the face: the corner
            // sphere ∪ cube kept its inside-cube triangle as if it
            // were the bulge and served a wrong body (measured:
            // union volume 0.82 vs 2.26).
            for k in 1..=25 {
                let s = march * k as f64;
                if contains(b.0 + nx * s * sign, b.1 + ny * s * sign) {
                    if lo.is_none() {
                        lo = Some(s);
                    }
                    hi = s;
                } else if lo.is_some() || k > 2 {
                    break;
                }
            }
            if let Some(lo) = lo
                && best_span.is_none_or(|(bs, ..)| hi - lo > bs)
            {
                best_span = Some((hi - lo, b, (nx, ny), sign, lo, hi));
            }
        }
    }
    if let Some((_, b, (nx, ny), sign, lo, hi)) = best_span {
        let s = 0.5 * (lo + hi);
        let (u, v) = (b.0 + nx * s * sign, b.1 + ny * s * sign);
        if contains(u, v) {
            if trace {
                eprintln!(
                    "  [sample] medial uv ({u:.4},{v:.4}) base ({:.4},{:.4}) span [{lo:.4},{hi:.4}] bbox u[{lo_u:.4},{hi_u:.4}] v[{lo_v:.4},{hi_v:.4}]",
                    b.0, b.1
                );
            }
            return Some(surface.point_at(u, v));
        }
    }
    // A thin fragment (a blade cross-section on a hub wall is a
    // slanted strip ~3% of the loop diagonal wide) is jumped clean
    // over by the 5% step on BOTH sides — refine the step until a
    // side is contained (PK-440). The smallest step is still orders
    // of magnitude above classification tolerance.
    for delta in [delta, delta * 0.2, delta * 0.04, delta * 0.008] {
        for sign in [1.0_f64, -1.0] {
            let (u, v) = (q.0 + px * delta * sign, q.1 + py * delta * sign);
            if contains(u, v) {
                return Some(surface.point_at(u, v));
            }
        }
    }
    if trace {
        eprintln!(
            "  [sample] interior step decline: neither side of ({:.4},{:.4})±{delta:.4} \
             contained in the {}-pt outer polygon of {face:?} bbox u[{lo_u:.4},{hi_u:.4}] \
             v[{lo_v:.4},{hi_v:.4}] first ({:.4},{:.4}) mid ({:.4},{:.4})",
            q.0,
            q.1,
            polygon.len(),
            polygon[0].0,
            polygon[0].1,
            polygon[polygon.len() / 2].0,
            polygon[polygon.len() / 2].1
        );
    }
    None
}

/// Whether any two non-adjacent polygon segments properly cross.
/// Quadratic, but only the all-shared interior-step fallback pays it,
/// on polygons of at most a few hundred points.
fn polygon_self_intersection(poly: &[(f64, f64)]) -> Option<(usize, usize)> {
    let n = poly.len();
    if n < 4 {
        return None;
    }
    let seg = |i: usize| (poly[i], poly[(i + 1) % n]);
    let cross = |o: (f64, f64), a: (f64, f64), b: (f64, f64)| {
        (a.0 - o.0) * (b.1 - o.1) - (a.1 - o.1) * (b.0 - o.0)
    };
    for i in 0..n {
        let (a, b) = seg(i);
        for j in (i + 2)..n {
            // Adjacent segments share a vertex; the closing segment
            // (n-1, 0) is adjacent to segment 0.
            if i == 0 && j == n - 1 {
                continue;
            }
            let (c, d) = seg(j);
            let d1 = cross(a, b, c);
            let d2 = cross(a, b, d);
            let d3 = cross(c, d, a);
            let d4 = cross(c, d, b);
            if ((d1 > 0.0 && d2 < 0.0) || (d1 < 0.0 && d2 > 0.0))
                && ((d3 > 0.0 && d4 < 0.0) || (d3 < 0.0 && d4 > 0.0))
            {
                return Some((i, j));
            }
        }
    }
    None
}
