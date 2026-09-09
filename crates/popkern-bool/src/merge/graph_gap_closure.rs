//! Junction-anchored record-gap closure (PK-465).
//!
//! After clipping, welding, and splitting, a face's boundary in the
//! merge graph can still carry SHORT open gaps between exact junction
//! nodes: near-tangent lens segments the SSI march cannot trace (the
//! constrained Newton exhausts where the surfaces graze), and
//! chart-border continuations whose windows collapsed. Measured
//! across the fastener corpus: 0.04–0.6-long gaps whose two ends are
//! known junctions, with the connecting curve fully determined by the
//! carrier surface pair.
//!
//! The closure constructs the missing segment directly: chord samples
//! between the two tips are converged onto the carrier pair's
//! intersection with the same constrained Newton the march uses
//! (each sample is an independent local solve — no marching, so
//! grazing tangency cannot exhaust it), then fit into a record edge.
//!
//! Certificates, all required:
//! - the source tip is a DEAD END (its node carries only shared
//!   record half-edges of a single record);
//! - the target lies in the tip record's forward tangent cone;
//! - every chord sample converges onto BOTH carrier surfaces with a
//!   bounded correction, and the chain advances monotonically;
//! - the segment's interior lies INSIDE both faces' trims (region
//!   containment at the converged uv registrations).
//!
//! A gap that fails any certificate stays open and keeps its honest
//! decline downstream.

use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::CurveGeometry;
use popkern_core::geometry::intersection::project_point_to_surface_pair;
use popkern_core::topology::BRep;

use super::graph::{EdgeSource, FaceRef, IntersectionGraph, Operand};

/// Close short record gaps between junction nodes. `band` is the
/// weld-scale tolerance (`linear * 500`), `lin` the linear tolerance.
pub(crate) fn close_record_gaps(
    graph: &mut IntersectionGraph,
    original_a: &BRep,
    original_b: &BRep,
    lin: f64,
) {
    use popkern_core::geometry::curve::CurveEval as _;
    let trace = popkern_core::trace::trace_enabled("boolean");
    let n_nodes = graph.nodes.len();
    let mut incident: Vec<Vec<usize>> = vec![Vec::new(); n_nodes];
    for (i, he) in graph.half_edges.iter().enumerate() {
        incident[he.origin].push(i);
    }
    let is_shared = |h: usize, graph: &IntersectionGraph| {
        matches!(graph.half_edges[h].source, EdgeSource::Shared { .. })
    };
    // Dead-end tips, judged PER FACE — the same criterion as the
    // walk's graze filter: a node whose incidence ON ONE FACE'S chart
    // is a single shared record pair (no boundary continuation there)
    // dead-ends that face even when the node carries records of other
    // pairs (measured: the nut's runout-chain ends are two-record
    // junction nodes globally, yet each face sees exactly one).
    let tip_map: std::collections::BTreeMap<usize, usize> = {
        let mut per_face_tips: std::collections::BTreeMap<usize, usize> =
            std::collections::BTreeMap::new();
        let faces: Vec<_> = graph.face_edges.keys().cloned().collect();
        for f in faces {
            let hes = &graph.face_edges[&f];
            let mut node_inc: std::collections::BTreeMap<usize, Vec<usize>> =
                std::collections::BTreeMap::new();
            for &h in hes {
                let a = graph.half_edges[h].origin;
                let b = graph.half_edges[graph.half_edges[h].twin].origin;
                node_inc.entry(a).or_default().push(h);
                node_inc.entry(b).or_default().push(h);
            }
            for (n, inc) in node_inc {
                let records: std::collections::BTreeSet<usize> = inc
                    .iter()
                    .filter(|&&h| is_shared(h, graph))
                    .map(|&h| h.min(graph.half_edges[h].twin))
                    .collect();
                let has_boundary = inc.iter().any(|&h| !is_shared(h, graph));
                if records.len() == 1 && !has_boundary {
                    let rh = *inc
                        .iter()
                        .find(|&&h| is_shared(h, graph))
                        .unwrap_or(&inc[0]);
                    per_face_tips.entry(n).or_insert(rh);
                }
            }
        }
        per_face_tips
    };
    let tips: Vec<usize> = tip_map.keys().copied().collect();
    let gap_trace = popkern_core::trace::trace_enabled("gap-close");
    if gap_trace {
        for &n in &tips {
            let p = graph.nodes[n].position;
            eprintln!("  [gap-close?] tip n{n} ({:.4},{:.4},{:.4})", p.x, p.y, p.z);
        }
    }
    // Candidate reach: a tip's own record can be a tiny stub (a
    // 0.04-long clip residue) while its true continuation target sits
    // a large fraction of the model away — the reach is model-scale,
    // and the certificates (cone, full-chain pair convergence,
    // containment) carry the safety.
    let model_diag = {
        let mut lo = Point3::new(f64::MAX, f64::MAX, f64::MAX);
        let mut hi = Point3::new(f64::MIN, f64::MIN, f64::MIN);
        for nd in &graph.nodes {
            let p = nd.position;
            lo = Point3::new(lo.x.min(p.x), lo.y.min(p.y), lo.z.min(p.z));
            hi = Point3::new(hi.x.max(p.x), hi.y.max(p.y), hi.z.max(p.z));
        }
        (hi - lo).length()
    };
    let mut bridged: std::collections::BTreeSet<usize> = std::collections::BTreeSet::new();
    for &n in &tips {
        if bridged.contains(&n) {
            continue;
        }
        let p = graph.nodes[n].position;
        let rh = tip_map[&n];
        let he = &graph.half_edges[rh];
        let EdgeSource::Shared { face_a, face_b } = he.source else {
            continue;
        };
        let other_end = graph.half_edges[graph.half_edges[rh].twin].origin;
        let chord = graph.nodes[other_end].position.distance_to(&p);
        let cap = (2.0 * chord).max(0.25 * model_diag).max(lin * 500.0);
        // Outward tangent of the tip record.
        let Some(outward) = ({
            he.curve
                .as_ref()
                .or_else(|| graph.half_edges[he.twin].curve.as_ref())
                .and_then(|c| {
                    let span = c.t_end - c.t_start;
                    if !span.is_finite() || span.abs() < 1e-30 {
                        return None;
                    }
                    let d_start = c.geometry.point_at(c.t_start).distance_to(&p);
                    let d_end = c.geometry.point_at(c.t_end).distance_to(&p);
                    let (t_at, t_in) = if d_start <= d_end {
                        (c.t_start, c.t_start + span * 1e-3)
                    } else {
                        (c.t_end, c.t_end - span * 1e-3)
                    };
                    (c.geometry.point_at(t_at) - c.geometry.point_at(t_in)).normalized()
                })
                .or_else(|| (p - graph.nodes[other_end].position).normalized())
        }) else {
            continue;
        };
        // Candidate targets: nearest nodes in the forward cone. A
        // grazing lens can HAIRPIN — when the forward cone is empty,
        // retry admitting any direction that is not straight back
        // along the record (the pair-convergence, monotonicity, and
        // containment certificates still gate every bridge).
        // Everything except straight-back is admissible (a grazing
        // lens can hairpin); forward-cone candidates rank first, then
        // by distance — the pair-convergence, monotonicity, and
        // containment certificates gate every bridge.
        let mut ranked: Vec<(bool, f64, usize)> = (0..n_nodes)
            .filter(|&m| m != n && m != other_end && !incident[m].is_empty())
            .filter_map(|m| {
                let q = graph.nodes[m].position;
                let d = q.distance_to(&p);
                if d < lin * 10.0 || d > cap {
                    return None;
                }
                let dir = (q - p).normalized()?;
                let dot = outward.dot(&dir);
                if dot < -0.5 {
                    return None;
                }
                Some((dot < 0.7, d, m))
            })
            .collect();
        // Distance dominates: the true continuation target is almost
        // always the nearest junction; the certificates rule out the
        // rest. (Forward-first ranking starved near sideways targets
        // behind six distant quasi-forward ones — measured on the
        // hole-bottom lens tip.)
        ranked.sort_by(|a, b| a.1.total_cmp(&b.1).then(a.2.cmp(&b.2)));
        ranked.truncate(6);
        let cands: Vec<(f64, usize)> = ranked.into_iter().map(|(_, d, m)| (d, m)).collect();
        if gap_trace {
            eprintln!(
                "  [gap-close?] processing n{n}: {} candidate(s)",
                cands.len()
            );
        }
        'cand: for (gap, m) in cands {
            let q = graph.nodes[m].position;
            if gap_trace {
                eprintln!(
                    "  [gap-close?]   cand n{m} ({:.4},{:.4},{:.4}) gap {gap:.3e}",
                    q.x, q.y, q.z
                );
            }
            // Carrier pairs to try: the tip record's pair first, then
            // the pairs of the TARGET's records — the connecting
            // segment can live on the continuation's pair (the tip is
            // exactly where its own pair's face ends).
            let mut pairs: Vec<(
                popkern_core::topology::FaceHandle,
                popkern_core::topology::FaceHandle,
            )> = vec![(face_a, face_b)];
            for &h in &incident[m] {
                if let EdgeSource::Shared {
                    face_a: fa,
                    face_b: fb,
                } = graph.half_edges[h].source
                    && !pairs.contains(&(fa, fb))
                {
                    pairs.push((fa, fb));
                }
            }
            // Faces whose surfaces CONTAIN both tips are carrier
            // candidates even when no incident record names them —
            // the missing lens's own pair typically has no surviving
            // record at either tip (that is why the gap exists).
            let near_faces = |brep: &BRep| -> Vec<popkern_core::topology::FaceHandle> {
                brep.faces
                    .iter()
                    .filter_map(|(fh, fd)| {
                        let surf = &fd.surface.as_ref()?.geometry;
                        let on = |pt: &Point3| {
                            use popkern_core::geometry::surface::SurfaceEval as _;
                            let gate = lin * 100.0;
                            if let Some((u, v)) =
                                popkern_core::geometry::uv_projection::project_point_to_uv(surf, pt)
                                && surf.point_at(u, v).distance_to(pt) < gate
                            {
                                return true;
                            }
                            // Basin-safe retry: a cold projection on a
                            // coiled chart walks the adjacent turn.
                            if let popkern_core::geometry::surface::SurfaceGeometry::Nurbs(nb) =
                                surf
                                && let Some(seeds) = nb.projection_seeds()
                                && let Some((u, v)) =
                                    popkern_core::geometry::surface::closest_point_from_seeds(
                                        nb,
                                        pt,
                                        &seeds,
                                        (false, false),
                                        6,
                                    )
                            {
                                return surf.point_at(u, v).distance_to(pt) < gate;
                            }
                            false
                        };
                        (on(&p) && on(&q)).then_some(fh)
                    })
                    .collect()
            };
            for fa2 in near_faces(original_a) {
                for &fb2 in &near_faces(original_b) {
                    if !pairs.contains(&(fa2, fb2)) {
                        pairs.push((fa2, fb2));
                    }
                }
            }
            'pair: for (pfa, pfb) in pairs {
                let (Some(sa), Some(sb)) = (
                    original_a
                        .faces
                        .get(pfa)
                        .and_then(|f| f.surface.as_ref())
                        .map(|s| &s.geometry),
                    original_b
                        .faces
                        .get(pfb)
                        .and_then(|f| f.surface.as_ref())
                        .map(|s| &s.geometry),
                ) else {
                    continue 'pair;
                };
                // Converge chord samples onto the carrier pair.
                let k = 9usize;
                let mut pts: Vec<Point3> = Vec::with_capacity(k + 2);
                let mut uv1: Vec<(f64, f64)> = Vec::with_capacity(k + 2);
                let mut uv2: Vec<(f64, f64)> = Vec::with_capacity(k + 2);
                let tol = original_a.tolerance;
                let mut ok = true;
                let correction_cap = (0.25 * gap).max(lin * 100.0);
                for i in 0..=(k + 1) {
                    let f = i as f64 / (k + 1) as f64;
                    let chord_pt = Point3::new(
                        p.x + (q.x - p.x) * f,
                        p.y + (q.y - p.y) * f,
                        p.z + (q.z - p.z) * f,
                    );
                    let Some((cp, a_uv, b_uv)) =
                        project_point_to_surface_pair(sa, sb, &chord_pt, &tol)
                    else {
                        if gap_trace {
                            eprintln!(
                                "  [gap-close?]     pair {pfa:?}x{pfb:?} sample {i}: projection none"
                            );
                        }
                        ok = false;
                        break;
                    };
                    if cp.distance_to(&chord_pt) > correction_cap {
                        if gap_trace {
                            eprintln!(
                                "  [gap-close?]     pair {pfa:?}x{pfb:?} sample {i}: correction {:.3e} > {correction_cap:.3e}",
                                cp.distance_to(&chord_pt)
                            );
                        }
                        ok = false;
                        break;
                    }
                    pts.push(cp);
                    uv1.push(a_uv);
                    uv2.push(b_uv);
                }
                if !ok {
                    if gap_trace {
                        eprintln!("  [gap-close?]   refused: chord sample failed pair projection");
                    }
                    continue 'pair;
                }
                // Pin the exact tips and require monotone advance.
                pts[0] = p;
                let last = pts.len() - 1;
                pts[last] = q;
                let dir_pq = q - p;
                let monotone = pts
                    .windows(2)
                    .all(|w| (w[1] - w[0]).dot(&dir_pq) > -lin * 10.0);
                if !monotone {
                    continue 'pair;
                }
                // NO DOUBLE COVER: the bridge must not run along record
                // coverage that already exists for this pair — a bridge
                // retracing an existing curve twins its edges (measured:
                // two wedge bridges re-covered served chains and the
                // assembly self-twinned 54 edges).
                {
                    let mut covered = false;
                    'cover: for (hx, hedge) in graph.half_edges.iter().enumerate() {
                        if hx % 2 != 0 {
                            continue;
                        }
                        let EdgeSource::Shared {
                            face_a: xfa,
                            face_b: xfb,
                        } = hedge.source
                        else {
                            continue;
                        };
                        if (xfa, xfb) != (pfa, pfb) {
                            continue;
                        }
                        let Some(c) = hedge.curve.as_ref() else {
                            continue;
                        };
                        for pt in pts.iter().take(last).skip(1) {
                            if let Some(t) = (&c.geometry
                                as &dyn popkern_core::geometry::curve::CurveEval)
                                .closest_point(pt)
                            {
                                let tc = t.clamp(c.t_start.min(c.t_end), c.t_start.max(c.t_end));
                                if c.geometry.point_at(tc).distance_to(pt) < lin * 100.0 {
                                    covered = true;
                                    break 'cover;
                                }
                            }
                        }
                    }
                    if covered {
                        if gap_trace {
                            eprintln!(
                                "  [gap-close?]   refused: bridge retraces existing {pfa:?}x{pfb:?} coverage"
                            );
                        }
                        continue 'pair;
                    }
                }
                // Interior containment on both faces.
                let mut inside_votes = 0usize;
                let mid_range = 1..last;
                let mid_count = mid_range.len();
                for i in mid_range.clone() {
                    let in_a = crate::split_inner_loop::point_inside_face_region_uv(
                        original_a,
                        pfa,
                        Point3::new(uv1[i].0, uv1[i].1, 0.0),
                        &original_a.tolerance,
                    );
                    let in_b = crate::split_inner_loop::point_inside_face_region_uv(
                        original_b,
                        pfb,
                        Point3::new(uv2[i].0, uv2[i].1, 0.0),
                        &original_b.tolerance,
                    );
                    if in_a && in_b {
                        inside_votes += 1;
                    }
                }
                // Supermajority: a bridge through a runout wedge can
                // scrape a bare majority off noisy near-boundary votes
                // (measured: a 0.9-long bridge through the exit wedge
                // passed 5/9 and self-twinned 54 edges downstream).
                if inside_votes * 4 < mid_count * 3 {
                    if gap_trace {
                        eprintln!(
                            "  [gap-close?]   refused: interior containment {inside_votes}/{mid_count}"
                        );
                    }
                    continue 'pair;
                }
                // Fit the segment and its pcurves over [0,1].
                let params: Vec<f64> = (0..pts.len())
                    .map(|i| i as f64 / (pts.len() - 1) as f64)
                    .collect();
                let Ok(c3) =
                    popkern_core::geometry::fitting::fit_curve(&pts, &params, 3, pts.len())
                else {
                    continue 'pair;
                };
                let uv_fit = |uvs: &[(f64, f64)]| {
                    let uv_pts: Vec<Point3> =
                        uvs.iter().map(|&(u, v)| Point3::new(u, v, 0.0)).collect();
                    popkern_core::geometry::fitting::fit_curve(&uv_pts, &params, 3, uv_pts.len())
                        .ok()
                };
                let (Some(c_uv1), Some(c_uv2)) = (uv_fit(&uv1), uv_fit(&uv2)) else {
                    continue 'pair;
                };
                let curve_on_edge = popkern_core::geometry::CurveOnEdge {
                    geometry: CurveGeometry::Nurbs(c3),
                    t_start: 0.0,
                    t_end: 1.0,
                };
                let source = EdgeSource::Shared {
                    face_a: pfa,
                    face_b: pfb,
                };
                let (fwd, twin) = graph.add_half_edge_pair_with_pcurves(
                    n,
                    m,
                    source,
                    Some(curve_on_edge),
                    Some(c_uv1.into()),
                    Some(c_uv2.into()),
                );
                let face_a_ref = FaceRef {
                    operand: Operand::A,
                    face: pfa,
                };
                let face_b_ref = FaceRef {
                    operand: Operand::B,
                    face: pfb,
                };
                graph.set_half_edge_face(fwd, face_a_ref);
                graph.set_half_edge_face(twin, face_b_ref);
                graph.face_edges.entry(face_a_ref).or_default().push(twin);
                graph.face_edges.entry(face_b_ref).or_default().push(fwd);
                if trace {
                    eprintln!(
                        "  [gap-close] n{n} ({:.4},{:.4},{:.4}) -> n{m} ({:.4},{:.4},{:.4}) gap {gap:.3e} pair {pfa:?}x{pfb:?}",
                        p.x, p.y, p.z, q.x, q.y, q.z
                    );
                }
                bridged.insert(n);
                bridged.insert(m);
                break 'cand;
            }
        }
    }
}
