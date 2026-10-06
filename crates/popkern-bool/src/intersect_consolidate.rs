//! Marching-record consolidation (split per file-size rule): drop
//! duplicate/covered records per face pair and weld nearly-closed
//! loops (see [`consolidate_records`]).

use super::{IntersectionRecord, endpoint_of, set_endpoint};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};

#[path = "intersect_consolidate_gap.rs"]
mod gap;

/// Consolidate marching output per face pair: drop geometric
/// duplicates of the same curve piece and weld junction endpoints of
/// pieces that chain (endpoint slop up to the marching class,
/// measured 3.5e-4..1.2e-3 — far above the linear tolerance, far
/// below any fragment's own extent).
pub(super) fn consolidate_records(
    records: &mut Vec<IntersectionRecord>,
    a: &BRep,
    b: &BRep,
    tolerance: &Tolerance,
) {
    use popkern_core::geometry::curve::CurveEval as _;
    use std::collections::HashMap;
    if records.len() < 2 {
        return;
    }
    let mut groups: HashMap<(FaceHandle, FaceHandle), Vec<usize>> = HashMap::new();
    for (i, r) in records.iter().enumerate() {
        groups.entry((r.face_a, r.face_b)).or_default().push(i);
    }
    // DETERMINISTIC group order: HashMap iteration order is
    // per-process random, and the dedupe/weld decisions below are
    // order-dependent (an earlier-processed group keeps its copy and
    // drops the later one) — the same boolean flip-flopped between
    // serve and decline across identical runs (PK-460, the flagged
    // intersect nondeterminism). Order groups by their smallest
    // record index: content-derived and stable.
    let mut ordered_groups: Vec<Vec<usize>> = groups.into_values().collect();
    ordered_groups.sort_by_key(|idxs| idxs.iter().copied().min().unwrap_or(usize::MAX));
    let n_groups = ordered_groups.len();
    let mut drop = vec![false; records.len()];
    // Junction weld targets: (record index, is_start, new point set).
    for idxs in &ordered_groups {
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
            // Arc-length position of `p`'s closest on-polyline witness
            // within the gate, or None when off the polyline.
            let arc_pos = |p: &popkern_core::geometry::Point3| -> Option<f64> {
                let mut best: Option<(f64, f64)> = None; // (dist, arc)
                let mut acc = 0.0;
                for w in poly_a.windows(2) {
                    let d = w[1] - w[0];
                    let l = d.dot(&d).max(1e-300).sqrt();
                    let t = ((*p - w[0]).dot(&d) / (l * l)).clamp(0.0, 1.0);
                    let q = w[0] + d * t;
                    let dist = q.distance_to(p);
                    if dist < gate && best.is_none_or(|(bd, _)| dist < bd) {
                        best = Some((dist, acc + t * l));
                    }
                    acc += l;
                }
                best.map(|(_, arc)| arc)
            };
            for b in 0..idxs.len() {
                if a == b || drop[idxs[b]] {
                    continue;
                }
                if info[b].3 > info[a].3 + gate {
                    continue;
                }
                let (Some(p0), Some(p2), Some(p1)) = (
                    arc_pos(&info[b].0),
                    arc_pos(&info[b].2),
                    arc_pos(&info[b].1),
                ) else {
                    continue;
                };
                // A duplicate OVERLAPS the polyline: its projections
                // spread along it comparably to its own extent. A
                // CHAINED continuation shorter than the gate projects
                // all three points onto one clamped endpoint (spread
                // ≈ 0) — dropping it severed a seam-sliver record and
                // stranded the assembly (measured: the spring-grind
                // section arc's 5.6e-4 piece across the wall seam,
                // PK-443).
                let spread = p0.max(p1).max(p2) - p0.min(p1).min(p2);
                if spread + tolerance.linear() >= 0.5 * info[b].3 {
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
        // a 1e-3 gap. Guarded by the record's OWN extent: a true loop
        // closes a gap far smaller than its diameter, while a short
        // exact CONNECTOR between two distinct junctions (rim crossing
        // to seam crossing, 5.6e-4 apart at a grazing corner —
        // PK-441) has a "gap" comparable to its whole extent and must
        // keep both endpoints.
        for &i in idxs.iter().filter(|&&i| !drop[i]) {
            let (s3, e3, ext) = {
                let c = &records[i].curve.curve_3d;
                let (t0, t1) = c.parameter_range();
                (
                    c.point_at(t0),
                    c.point_at(t1),
                    c.point_at(t0).distance_to(&c.point_at(0.5 * (t0 + t1))),
                )
            };
            let d = s3.distance_to(&e3);
            if d > tolerance.linear() && d < gate && d < 0.5 * ext {
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
        // endpoint control points, which both charts share). The weld
        // distance is capped by a quarter of the SHORTER record's own
        // extent: moving an endpoint further materially deforms the
        // record, and a short exact connector's two junctions (5.6e-4
        // apart at a grazing rim × seam corner, both Newton-refined —
        // PK-441) must not be averaged into one.
        let survivors: Vec<usize> = idxs.iter().copied().filter(|&i| !drop[i]).collect();
        let survivor_extents: Vec<f64> = survivors
            .iter()
            .map(|&i| {
                let c = &records[i].curve.curve_3d;
                let (t0, t1) = c.parameter_range();
                let s = c.point_at(t0);
                s.distance_to(&c.point_at(t1))
                    .max(s.distance_to(&c.point_at(0.5 * (t0 + t1))))
            })
            .collect();
        for a_pos in 0..survivors.len() {
            for b_pos in (a_pos + 1)..survivors.len() {
                let (ia, ib) = (survivors[a_pos], survivors[b_pos]);
                let pair_gate =
                    gate.min(0.25 * survivor_extents[a_pos].min(survivor_extents[b_pos]));
                for (a_end, b_end) in [(false, true), (true, false), (false, false), (true, true)] {
                    let pa = endpoint_of(&records[ia], a_end);
                    let pb = endpoint_of(&records[ib], b_end);
                    let d = pa.distance_to(&pb);
                    if d > tolerance.linear()
                        && d < pair_gate
                        && !gap::gap_leaves_a_face(
                            a,
                            b,
                            &records[ia],
                            a_end,
                            &records[ib],
                            b_end,
                            tolerance,
                        )
                    {
                        // KEPT-END AUTHORITY (PK-480): a kept record's
                        // ends are junction-refined at clip time (near
                        // pin, far re-trim, triple-point junctions); a
                        // trim residue's ends are deliberately RAW
                        // (refining them re-strands the m6 mouth
                        // chamfer). Averaging the two moved BOTH onto
                        // a midpoint off the exact junction — the
                        // deep-plunged tap's crest triple points
                        // landed 1.1e-3 off, past the dead-end weld
                        // band, and the cone chart stranded. A
                        // kept×residue weld adopts the kept end;
                        // same-class pairs keep the midpoint.
                        let a_res = records[ia].residue_owner_a.is_some();
                        let b_res = records[ib].residue_owner_a.is_some();
                        if a_res != b_res && popkern_core::trace::trace_enabled("boolean") {
                            eprintln!(
                                "  [weld-kept] pair ({:?},{:?}) kept ({:.6},{:.6},{:.6}) residue ({:.6},{:.6},{:.6}) d {d:.3e}",
                                (records[ia].face_a, records[ia].face_b),
                                (records[ib].face_a, records[ib].face_b),
                                if a_res { pb.x } else { pa.x },
                                if a_res { pb.y } else { pa.y },
                                if a_res { pb.z } else { pa.z },
                                if a_res { pa.x } else { pb.x },
                                if a_res { pa.y } else { pb.y },
                                if a_res { pa.z } else { pb.z },
                            );
                        }
                        let weld = if a_res == b_res {
                            popkern_core::geometry::Point3::new(
                                0.5 * (pa.x + pb.x),
                                0.5 * (pa.y + pb.y),
                                0.5 * (pa.z + pb.z),
                            )
                        } else if a_res {
                            pb
                        } else {
                            pa
                        };
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
            n_groups,
            dropped
        );
    }
    let mut keep_iter = drop.into_iter();
    records.retain(|_| !keep_iter.next().unwrap_or(false));
}
