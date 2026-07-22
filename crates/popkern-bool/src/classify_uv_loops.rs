//! UV-loop containment for ray classification: pcurve-sampled
//! even-odd polygon tests for outer and inner face loops.
//!
//! Split from `classify.rs` to keep it within the structural line
//! ratchet.

use popkern_core::topology::{BRep, FaceHandle};

/// Sample a loop's pcurves into a UV polygon and even-odd test
/// `(u, v)` against it. Returns `None` when the loop has no usable
/// pcurve polygon (fewer than 3 samples) — the caller must treat
/// that as "no verdict", not "outside".
pub(crate) fn uv_point_in_loop(
    brep: &BRep,
    lh: popkern_core::topology::LoopHandle,
    u: f64,
    v: f64,
) -> Option<bool> {
    let loop_data = brep.loops.get(lh)?;
    let start = loop_data.half_edge;
    if !brep.half_edges.contains_key(start) {
        return None;
    }
    let mut ring: Vec<(f64, f64)> = Vec::new();
    let mut he = start;
    let cap = brep.half_edges.len() + 1;
    for _ in 0..cap {
        let d = &brep.half_edges[he];
        if let Some(pc) = d.pcurve.as_ref() {
            let (t0, t1) = pc.parameter_range();
            for i in 0..8 {
                let t = t0 + (t1 - t0) * i as f64 / 8.0;
                let p = pc.point_at(t);
                ring.push((p.x, p.y));
            }
        }
        he = d.next;
        if he == start {
            break;
        }
    }
    if ring.len() < 3 {
        return None;
    }
    let mut crossings = 0;
    let n = ring.len();
    for i in 0..n {
        let (x0, y0) = ring[i];
        let (x1, y1) = ring[(i + 1) % n];
        if (y0 > v) != (y1 > v) {
            let xc = x0 + (v - y0) / (y1 - y0) * (x1 - x0);
            if xc > u {
                crossings += 1;
            }
        }
    }
    Some(crossings % 2 == 1)
}

/// True when the surface hit at `(u, v)` lands inside the face's
/// OUTER pcurve polygon. The verdict is trusted when the outer ring
/// is COMPLETE (every half-edge carries a pcurve). Linear pcurves
/// contribute their exact endpoints; CURVED pcurves are sampled
/// densely (64 stations — sagitta ≲ 0.2% of the span), and hits
/// within a dead-band of the polygon boundary abstain (keep the
/// hit) so sampling error can never reject a rim-grazing crossing
/// (measured: an 8-station circle misjudged a 2.7 mm band at r=35
/// and broke star-socket caps). A hit clearly outside — like the
/// 0.67-off phantom on an ellipse-fronted top face that flipped a
/// tool fragment to Inside — is rejected. Partial rings
/// mid-pipeline still abstain entirely.
pub(crate) fn hit_inside_outer_loop(brep: &BRep, face: FaceHandle, u: f64, v: f64) -> bool {
    let Some(face_data) = brep.faces.get(face) else {
        return true;
    };
    let Some(loop_data) = brep.loops.get(face_data.outer_loop) else {
        return true;
    };
    let start = loop_data.half_edge;
    if !brep.half_edges.contains_key(start) {
        return true;
    }
    // Collect each half-edge's pcurve segment, then chain them
    // orientation-aware: a stored pcurve may run against the loop
    // traversal (twin-shared spans), so each segment is flipped to
    // continue from the previous tail.
    let mut has_curved = false;
    let mut segs: Vec<((f64, f64), (f64, f64))> = Vec::new();
    let mut he = start;
    let cap = brep.half_edges.len() + 1;
    for _ in 0..cap {
        let d = &brep.half_edges[he];
        let Some(pc) = d.pcurve.as_ref() else {
            return true; // incomplete ring — no verdict
        };
        let linear_content = match &pc.geometry {
            popkern_core::geometry::pcurve::PcurveGeometry::Linear { .. } => true,
            popkern_core::geometry::pcurve::PcurveGeometry::Nurbs(c) => c.degree == 1,
            _ => false,
        };
        let (t0, t1) = pc.parameter_range();
        if linear_content {
            let p0 = pc.point_at(t0);
            let p1 = pc.point_at(t1);
            segs.push(((p0.x, p0.y), (p1.x, p1.y)));
        } else {
            has_curved = true;
            // Curved pcurve: densify into chained sub-segments.
            let n = 64;
            let mut prev = pc.point_at(t0);
            for i in 1..=n {
                let t = t0 + (t1 - t0) * i as f64 / n as f64;
                let q = pc.point_at(t);
                segs.push(((prev.x, prev.y), (q.x, q.y)));
                prev = q;
            }
        }
        he = d.next;
        if he == start {
            break;
        }
    }
    if segs.len() < 3 {
        return true;
    }
    let mut ring: Vec<(f64, f64)> = Vec::with_capacity(segs.len());
    let d2 = |a: (f64, f64), b: (f64, f64)| (a.0 - b.0).powi(2) + (a.1 - b.1).powi(2);
    let (mut tail, first_tail) = {
        // Orient the first segment by whichever end the SECOND
        // segment continues from.
        let (a0, a1) = segs[0];
        let (b0, b1) = segs[1];
        let fwd = d2(a1, b0).min(d2(a1, b1));
        let rev = d2(a0, b0).min(d2(a0, b1));
        if fwd <= rev {
            ring.push(a0);
            (a1, a0)
        } else {
            ring.push(a1);
            (a0, a1)
        }
    };
    for &(p0, p1) in segs.iter().skip(1) {
        if d2(tail, p0) <= d2(tail, p1) {
            ring.push(p0);
            tail = p1;
        } else {
            ring.push(p1);
            tail = p0;
        }
    }
    let _ = first_tail;
    // Dead-band: within sampling error of the boundary, abstain
    // (keep the hit) — the polygon is approximate for curved rims.
    let n = ring.len();
    let (mut lo_u, mut hi_u, mut lo_v, mut hi_v) = (f64::MAX, f64::MIN, f64::MAX, f64::MIN);
    for &(x, y) in &ring {
        lo_u = lo_u.min(x);
        hi_u = hi_u.max(x);
        lo_v = lo_v.min(y);
        hi_v = hi_v.max(y);
    }
    // The dead-band exists for SAMPLING error only: linear rings are
    // exact and get none (a 1 mm band on an exact ring re-admitted
    // the very phantoms the gate exists to reject).
    let band = if has_curved {
        0.002 * ((hi_u - lo_u) + (hi_v - lo_v)).max(1e-9)
    } else {
        0.0
    };
    if band > 0.0 {
        let mut min_d = f64::MAX;
        for i in 0..n {
            let (x0, y0) = ring[i];
            let (x1, y1) = ring[(i + 1) % n];
            let (dx, dy) = (x1 - x0, y1 - y0);
            let l2 = (dx * dx + dy * dy).max(1e-300);
            let t = (((u - x0) * dx + (v - y0) * dy) / l2).clamp(0.0, 1.0);
            let d = ((u - x0 - dx * t).powi(2) + (v - y0 - dy * t).powi(2)).sqrt();
            min_d = min_d.min(d);
        }
        if min_d < band {
            return true;
        }
    }
    let mut crossings = 0;
    for i in 0..n {
        let (x0, y0) = ring[i];
        let (x1, y1) = ring[(i + 1) % n];
        if (y0 > v) != (y1 > v) {
            let xc = x0 + (v - y0) / (y1 - y0) * (x1 - x0);
            if xc > u {
                crossings += 1;
            }
        }
    }
    crossings % 2 == 1
}

/// True when the surface hit at `(u, v)` lands inside one of the
/// face's INNER-loop polygons (pcurve-sampled, even-odd).
pub(crate) fn hit_inside_inner_loop(brep: &BRep, face: FaceHandle, u: f64, v: f64) -> bool {
    let Some(face_data) = brep.faces.get(face) else {
        return false;
    };
    face_data
        .inner_loops
        .iter()
        .any(|&lh| uv_point_in_loop(brep, lh, u, v) == Some(true))
}

// ---------------------------------------------------------------------------
// 4. Point-in-solid with multi-ray voting
// ---------------------------------------------------------------------------
