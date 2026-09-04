//! Marching-record consolidation (split per file-size rule): drop
//! duplicate/covered records per face pair and weld nearly-closed
//! loops (see [`consolidate_records`]).

use super::{IntersectionRecord, endpoint_of, set_endpoint};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::FaceHandle;

/// Consolidate marching output per face pair: drop geometric
/// duplicates of the same curve piece and weld junction endpoints of
/// pieces that chain (endpoint slop up to the marching class,
/// measured 3.5e-4..1.2e-3 — far above the linear tolerance, far
/// below any fragment's own extent).
pub(super) fn consolidate_records(records: &mut Vec<IntersectionRecord>, tolerance: &Tolerance) {
    use popkern_core::geometry::curve::CurveEval as _;
    use std::collections::HashMap;
    if records.len() < 2 {
        return;
    }
    let mut groups: HashMap<(FaceHandle, FaceHandle), Vec<usize>> = HashMap::new();
    for (i, r) in records.iter().enumerate() {
        groups.entry((r.face_a, r.face_b)).or_default().push(i);
    }
    let mut drop = vec![false; records.len()];
    // Junction weld targets: (record index, is_start, new point set).
    for idxs in groups.values() {
        if idxs.len() < 2 {
            continue;
        }
        let info: Vec<(
            popkern_core::geometry::Point3,
            popkern_core::geometry::Point3,
            popkern_core::geometry::Point3,
            f64,
        )> = idxs
            .iter()
            .map(|&i| {
                let c = &records[i].curve.curve_3d;
                let (t0, t1) = c.parameter_range();
                let s = c.point_at(t0);
                let e = c.point_at(t1);
                let m = c.point_at(0.5 * (t0 + t1));
                // Extent, not chord: a CLOSED loop's start-to-end
                // distance is ~0, so a chord-based gate collapsed to
                // 10·tol — far below the ~1e-3 slop between two
                // independent marching traces of the same loop, and a
                // duplicate girdle loop survived to imprint a
                // valence-1 spur that broke the cellular pipeline
                // (PK-440: propeller blade ∪ hub). The start-to-mid
                // distance measures a loop's diameter; for open
                // records it stays in the chord's class.
                let extent = s.distance_to(&e).max(s.distance_to(&m));
                (s, e, m, extent)
            })
            .collect();
        // Median extent: a single short fragment must not collapse
        // the junction gate below the marching slop class.
        let mut extents: Vec<f64> = info.iter().map(|x| x.3).collect();
        extents.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        let median_extent = extents[extents.len() / 2].max(tolerance.linear());
        let gate = (0.05 * median_extent).max(tolerance.linear() * 10.0);
        // Duplicates: same endpoints (either direction) + same midpoint.
        for a in 0..idxs.len() {
            if drop[idxs[a]] {
                continue;
            }
            for b in (a + 1)..idxs.len() {
                if drop[idxs[b]] {
                    continue;
                }
                let fwd = info[a].0.distance_to(&info[b].0) < gate
                    && info[a].1.distance_to(&info[b].1) < gate;
                let rev = info[a].0.distance_to(&info[b].1) < gate
                    && info[a].1.distance_to(&info[b].0) < gate;
                if (fwd || rev) && info[a].2.distance_to(&info[b].2) < gate {
                    drop[idxs[b]] = true;
                }
            }
        }
        // COVERED records: junction-adjacent marching can deliver a
        // span AND overlapping micro-fragments of the same span from
        // different seeds (measured on the rotated torus × torus
        // union: pair (9v1,7v1) carried [0→0.025], [0→0.003] and
        // [0.005→0.025] of one junction arc — the walk chewed the
        // overlap into duplicate edge pairs and stranded halves).
        // A record whose start, mid and end all lie ON another
        // record's sampled polyline, with smaller-or-equal extent, is
        // duplicate authority and drops.
        for a in 0..idxs.len() {
            if drop[idxs[a]] {
                continue;
            }
            let poly_a: Vec<popkern_core::geometry::Point3> = {
                let c = &records[idxs[a]].curve.curve_3d;
                let (t0, t1) = c.parameter_range();
                (0..=32)
                    .map(|k| c.point_at(t0 + (t1 - t0) * k as f64 / 32.0))
                    .collect()
            };
            let on_poly = |p: &popkern_core::geometry::Point3| -> bool {
                poly_a.windows(2).any(|w| {
                    let d = w[1] - w[0];
                    let l2 = d.dot(&d).max(1e-300);
                    let t = ((*p - w[0]).dot(&d) / l2).clamp(0.0, 1.0);
                    let q = w[0] + d * t;
                    q.distance_to(p) < gate
                })
            };
            for b in 0..idxs.len() {
                if a == b || drop[idxs[b]] {
                    continue;
                }
                if info[b].3 <= info[a].3 + gate
                    && on_poly(&info[b].0)
                    && on_poly(&info[b].2)
                    && on_poly(&info[b].1)
                {
                    drop[idxs[b]] = true;
                }
            }
        }
        // Loop-level dedupe: the walker can deliver the SAME closed
        // intersection loop several times, split at different
        // parameter offsets (piecewise dedupe can't see those).
        // Chain the survivors by endpoints; closed chains whose
        // sampled geometry matches an earlier chain are dropped
        // whole.
        {
            let alive: Vec<usize> = idxs.iter().copied().filter(|&i| !drop[i]).collect();
            let ends_of = |i: usize| -> (
                popkern_core::geometry::Point3,
                popkern_core::geometry::Point3,
            ) {
                let c = &records[i].curve.curve_3d;
                let (t0, t1) = c.parameter_range();
                (c.point_at(t0), c.point_at(t1))
            };
            let mut used = vec![false; alive.len()];
            let mut kept_loops: Vec<Vec<popkern_core::geometry::Point3>> = Vec::new();
            for seed in 0..alive.len() {
                if used[seed] {
                    continue;
                }
                let mut chain: Vec<(usize, bool)> = vec![(seed, false)];
                let mut chain_used = used.clone();
                chain_used[seed] = true;
                let (start, mut tip) = ends_of(alive[seed]);
                // A single record can be a whole loop on its own (the
                // marcher closes the girdle in one trace): it counts
                // as closed when its ends meet AND it has real extent
                // (start to parametric midpoint beyond the gate), so a
                // second copy of the same loop is deduped even though
                // no chaining happened.
                let single_loop_extent = {
                    let c = &records[alive[seed]].curve.curve_3d;
                    let (t0, t1) = c.parameter_range();
                    c.point_at(t0).distance_to(&c.point_at(0.5 * (t0 + t1)))
                };
                let closed = loop {
                    if tip.distance_to(&start) < gate
                        && (chain.len() >= 2 || single_loop_extent > gate)
                    {
                        break true;
                    }
                    let mut ext = false;
                    for k in 0..alive.len() {
                        if chain_used[k] {
                            continue;
                        }
                        let (ks, ke) = ends_of(alive[k]);
                        if ks.distance_to(&tip) < gate {
                            chain.push((k, false));
                            chain_used[k] = true;
                            tip = ke;
                            ext = true;
                            break;
                        }
                        if ke.distance_to(&tip) < gate {
                            chain.push((k, true));
                            chain_used[k] = true;
                            tip = ks;
                            ext = true;
                            break;
                        }
                    }
                    if !ext {
                        break false;
                    }
                };
                if popkern_core::trace::trace_enabled("boolean") {
                    eprintln!(
                        "  [consolidate] chain from seed {seed}: len {} closed={closed}",
                        chain.len()
                    );
                }
                if !closed {
                    continue;
                }
                // Sample in chain-traversal order (reversed pieces
                // sample backward — otherwise the comparison
                // polyline zigzags across the loop).
                let samples: Vec<popkern_core::geometry::Point3> = chain
                    .iter()
                    .flat_map(|&(k, rev)| {
                        let c = &records[alive[k]].curve.curve_3d;
                        let (t0, t1) = c.parameter_range();
                        (0..8).map(move |q| {
                            let f = q as f64 / 8.0;
                            let f = if rev { 1.0 - f } else { f };
                            c.point_at(t0 + (t1 - t0) * f)
                        })
                    })
                    .collect();
                // Distance to the kept loop's POLYLINE (segment
                // projection) — sample-to-sample comparison misses
                // duplicates whose samples interleave.
                let dist_to_polyline =
                    |p: &popkern_core::geometry::Point3,
                     poly: &[popkern_core::geometry::Point3]| {
                        let mut best = f64::INFINITY;
                        for w in 0..poly.len() {
                            let a = poly[w];
                            let b = poly[(w + 1) % poly.len()];
                            let ab = b - a;
                            let len2 = ab.dot(&ab);
                            let t = if len2 > 1e-300 {
                                ((*p - a).dot(&ab) / len2).clamp(0.0, 1.0)
                            } else {
                                0.0
                            };
                            let q = a + ab * t;
                            best = best.min(p.distance_to(&q));
                        }
                        best
                    };
                let duplicate = kept_loops
                    .iter()
                    .any(|kl| samples.iter().all(|p| dist_to_polyline(p, kl) < gate));
                if popkern_core::trace::trace_enabled("boolean") {
                    let n = samples.len() as f64;
                    let cx = samples.iter().map(|p| p.x).sum::<f64>() / n;
                    let cy = samples.iter().map(|p| p.y).sum::<f64>() / n;
                    let cz = samples.iter().map(|p| p.z).sum::<f64>() / n;
                    let c = popkern_core::geometry::Point3::new(cx, cy, cz);
                    let rmin = samples
                        .iter()
                        .map(|p| p.distance_to(&c))
                        .fold(f64::INFINITY, f64::min);
                    let rmax = samples
                        .iter()
                        .map(|p| p.distance_to(&c))
                        .fold(0.0_f64, f64::max);
                    eprintln!(
                        "  [consolidate] duplicate={duplicate} centroid ({cx:.3},{cy:.3},{cz:.3}) r {rmin:.4}..{rmax:.4} gate {gate:.2e}"
                    );
                }
                used = chain_used;
                if duplicate {
                    for &(k, _) in &chain {
                        drop[alive[k]] = true;
                    }
                } else {
                    kept_loops.push(samples);
                }
            }
        }

        // Self-weld: a nearly-closed loop record (closure gap within
        // the marching slop) becomes exactly closed — 3D and both UV
        // charts — so the graph sees one node, not an open arc with
        // a 1e-3 gap.
        for &i in idxs.iter().filter(|&&i| !drop[i]) {
            let (s3, e3) = {
                let c = &records[i].curve.curve_3d;
                let (t0, t1) = c.parameter_range();
                (c.point_at(t0), c.point_at(t1))
            };
            let d = s3.distance_to(&e3);
            if d > tolerance.linear() && d < gate {
                let weld = popkern_core::geometry::Point3::new(
                    0.5 * (s3.x + e3.x),
                    0.5 * (s3.y + e3.y),
                    0.5 * (s3.z + e3.z),
                );
                set_endpoint(&mut records[i], false, weld);
                set_endpoint(&mut records[i], true, weld);
                for uv in [true, false] {
                    let c = if uv {
                        &mut records[i].curve.curve_uv1
                    } else {
                        &mut records[i].curve.curve_uv2
                    };
                    if let (Some(first), Some(last)) = (
                        c.control_points.first().copied(),
                        c.control_points.last().copied(),
                    ) {
                        let m = popkern_core::geometry::Point3::new(
                            0.5 * (first.x + last.x),
                            0.5 * (first.y + last.y),
                            0.5 * (first.z + last.z),
                        );
                        if let Some(cp) = c.control_points.first_mut() {
                            *cp = m;
                        }
                        if let Some(cp) = c.control_points.last_mut() {
                            *cp = m;
                        }
                    }
                }
            }
        }

        // Weld junctions among the survivors: cluster endpoints and
        // move each cluster's endpoints (3D + both UV charts stay
        // consistent because the curves are only ADJUSTED at their
        // endpoint control points, which both charts share).
        let survivors: Vec<usize> = idxs.iter().copied().filter(|&i| !drop[i]).collect();
        for a_pos in 0..survivors.len() {
            for b_pos in (a_pos + 1)..survivors.len() {
                let (ia, ib) = (survivors[a_pos], survivors[b_pos]);
                for (a_end, b_end) in [(false, true), (true, false), (false, false), (true, true)] {
                    let pa = endpoint_of(&records[ia], a_end);
                    let pb = endpoint_of(&records[ib], b_end);
                    let d = pa.distance_to(&pb);
                    if d > tolerance.linear() && d < gate {
                        let weld = popkern_core::geometry::Point3::new(
                            0.5 * (pa.x + pb.x),
                            0.5 * (pa.y + pb.y),
                            0.5 * (pa.z + pb.z),
                        );
                        set_endpoint(&mut records[ia], a_end, weld);
                        set_endpoint(&mut records[ib], b_end, weld);
                    }
                }
            }
        }
    }
    if popkern_core::trace::trace_enabled("boolean") {
        let dropped = drop.iter().filter(|&&d| d).count();
        eprintln!(
            "  [consolidate] {} records in {} groups, dropped {}",
            records.len(),
            groups.len(),
            dropped
        );
        for (k, v) in &groups {
            eprintln!("    group {k:?}: {} records", v.len());
        }
    }
    let mut keep_iter = drop.into_iter();
    records.retain(|_| !keep_iter.next().unwrap_or(false));
}
