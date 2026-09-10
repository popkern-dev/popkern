//! Point-loop retirement for boolean outputs (split from
//! `slit_retire.rs` under the 500-line file budget) — PK-465.

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

/// Retire ALL point-loop edge pairs (see
/// [`retire_one_point_loop_pair`]); returns how many were removed.
/// Public entry for the PRE-validation cleanup: these slivers are
/// weld residue and must retire before the serve gate measures loop
/// closure (the postprocess passes run after it).
pub fn retire_point_loop_pairs(brep: &mut BRep, tolerance: &Tolerance) -> usize {
    // Weld-slop scale, capped relative to the body's extent: at a
    // coarse model tolerance linear×500 reaches FEATURE scale and
    // retired a real hole rim (measured: the plate fixture's 0.6-r
    // hole under 1e-2 tolerance).
    let mut lo = popkern_core::geometry::Point3::new(f64::INFINITY, f64::INFINITY, f64::INFINITY);
    let mut hi = popkern_core::geometry::Point3::new(
        f64::NEG_INFINITY,
        f64::NEG_INFINITY,
        f64::NEG_INFINITY,
    );
    for (_, v) in &brep.vertices {
        lo.x = lo.x.min(v.position.x);
        lo.y = lo.y.min(v.position.y);
        lo.z = lo.z.min(v.position.z);
        hi.x = hi.x.max(v.position.x);
        hi.y = hi.y.max(v.position.y);
        hi.z = hi.z.max(v.position.z);
    }
    let diag = ((hi.x - lo.x).powi(2) + (hi.y - lo.y).powi(2) + (hi.z - lo.z).powi(2)).sqrt();
    let gate = (tolerance.linear() * 500.0).min(diag.max(0.0) * 1.0e-4);
    let mut n = 0;
    while retire_one_point_loop_pair(brep, gate) {
        n += 1;
    }
    n
}

/// Retire ONE point-loop edge pair: a twin pair whose two half-edges
/// originate at the SAME vertex and whose curve extent is below
/// `gate`. Splices both halves out of their loops (excise pattern),
/// with snapshot rollback on structural regression.
fn retire_one_point_loop_pair(brep: &mut BRep, gate: f64) -> bool {
    use popkern_core::geometry::curve::CurveEval as _;
    let candidates: Vec<_> = brep
        .edges
        .iter()
        .filter_map(|(_, d)| {
            let heh = d.half_edge;
            let he = brep.half_edges.get(heh)?;
            if he.twin == heh {
                return None; // self-twin: the stub pass owns it
            }
            let twin = brep.half_edges.get(he.twin)?;
            if he.origin != twin.origin {
                return None;
            }
            // Extent by dense sampling over BOTH the stored window
            // and the geometry's own parameter range: a closed rim
            // can arrive with a degenerate window (t_start ==
            // t_end), and a two-sample extent read it as a point —
            // retiring a real hole rim (measured: the
            // cylinder-through-plate fixture lost its top hole).
            let extent = d
                .curve
                .as_ref()
                .map(|c| {
                    let (t0, t1) = (c.t_start, c.t_end);
                    // A degenerate stored window measures nothing —
                    // treat as unmeasurable, never as a point (a
                    // closed rim can arrive with t_start == t_end).
                    if (t1 - t0).abs() < 1e-12 {
                        return f64::INFINITY;
                    }
                    let mut lo = popkern_core::geometry::Point3::new(
                        f64::INFINITY,
                        f64::INFINITY,
                        f64::INFINITY,
                    );
                    let mut hi = popkern_core::geometry::Point3::new(
                        f64::NEG_INFINITY,
                        f64::NEG_INFINITY,
                        f64::NEG_INFINITY,
                    );
                    for k in 0..=8 {
                        let t = t0 + (t1 - t0) * k as f64 / 8.0;
                        let q = c.geometry.point_at(t);
                        lo.x = lo.x.min(q.x);
                        lo.y = lo.y.min(q.y);
                        lo.z = lo.z.min(q.z);
                        hi.x = hi.x.max(q.x);
                        hi.y = hi.y.max(q.y);
                        hi.z = hi.z.max(q.z);
                    }
                    (hi.x - lo.x) + (hi.y - lo.y) + (hi.z - lo.z)
                })
                // No curve = no measurement — with one exception: a
                // masked edge's curve is deliberately stripped
                // (mask_dangling) while its true geometry can be a
                // full hole rim, so an unmeasurable pair spanning TWO
                // faces stays. A curveless point-loop on ONE face has
                // no such alibi: it is masked residue (measured: the
                // angled tap's above-plate remnants at z 9.8-10.3).
                .unwrap_or_else(|| {
                    let same_face = (|| {
                        let f1 = brep.loops.get(he.loop_)?.face;
                        let f2 = brep.loops.get(twin.loop_)?.face;
                        Some(f1 == f2)
                    })()
                    .unwrap_or(false);
                    if same_face { 0.0 } else { f64::INFINITY }
                });
            if extent > gate {
                return None;
            }
            if popkern_core::trace::trace_enabled("spur") {
                let p = brep.vertices.get(he.origin).map(|v| v.position);
                eprintln!(
                    "[spur] point-loop candidate {heh:?}: extent {extent:.3e} curve={} window=({:?}) at {:?}",
                    d.curve.is_some(),
                    d.curve.as_ref().map(|c| (c.t_start, c.t_end)),
                    p
                );
            }
            Some(heh)
        })
        .collect();
    for heh in candidates {
        if !brep.half_edges.contains_key(heh) {
            continue;
        }
        let snapshot = brep.clone();
        let twin_h = brep.half_edges[heh].twin;
        let edge_h = brep.half_edges[heh].edge;
        let vertex = brep.half_edges[heh].origin;
        let mut ok = true;
        for h in [heh, twin_h] {
            if !brep.half_edges.contains_key(h) {
                continue;
            }
            let loop_h = brep.half_edges[h].loop_;
            let next_h = brep.half_edges[h].next;
            if next_h == h {
                // The half is its entire loop: retire an inner loop.
                let Some(loop_data) = brep.loops.get(loop_h) else {
                    ok = false;
                    break;
                };
                let face_h = loop_data.face;
                let Some(face) = brep.faces.get_mut(face_h) else {
                    ok = false;
                    break;
                };
                if face.outer_loop == loop_h {
                    ok = false;
                    break;
                }
                face.inner_loops.retain(|&l| l != loop_h);
                brep.loops.remove(loop_h);
                brep.remove_half_edge_reanchored(h);
                continue;
            }
            // Splice past this half.
            let mut prev_h = next_h;
            let mut steps = 0usize;
            loop {
                let Some(prev) = brep.half_edges.get(prev_h) else {
                    ok = false;
                    break;
                };
                if prev.next == h {
                    break;
                }
                prev_h = prev.next;
                steps += 1;
                if steps > 100_000 {
                    ok = false;
                    break;
                }
            }
            if !ok {
                break;
            }
            brep.half_edges[prev_h].next = next_h;
            brep.half_edges[next_h].prev = prev_h;
            if let Some(loop_data) = brep.loops.get_mut(loop_h)
                && loop_data.half_edge == h
            {
                loop_data.half_edge = next_h;
            }
            brep.remove_half_edge_reanchored(h);
        }
        if ok {
            if !brep.half_edges.iter().any(|(_, he)| he.edge == edge_h) {
                brep.edges.remove(edge_h);
            }
            match brep
                .half_edges
                .iter()
                .find(|(_, he)| he.origin == vertex)
                .map(|(hh, _)| hh)
            {
                Some(out_he) => {
                    if let Some(v) = brep.vertices.get_mut(vertex) {
                        v.outgoing = out_he;
                    }
                }
                None => {
                    brep.vertices.remove(vertex);
                }
            }
            let before = snapshot.validate_collect_structural().len();
            let after = brep.validate_collect_structural().len();
            if after <= before {
                if popkern_core::trace::trace_enabled("spur") {
                    eprintln!("[spur] point-loop pair {heh:?} retired");
                }
                return true;
            }
        }
        *brep = snapshot;
    }
    false
}
