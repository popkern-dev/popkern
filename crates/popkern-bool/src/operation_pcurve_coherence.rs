//! Post-merge pcurve coherence sweep for doubly-periodic (torus)
//! faces (extracted from `operation.rs` under the file-size ratchet).
//!
//! A record pcurve transplanted from the merge graph can run the LONG
//! way around the chart (measured: a served torus × torus face
//! carried one use spanning 5.77 rad — the loop closed only modulo a
//! full period, the exporter's branch normalization saw a 2π gap, and
//! the tessellator's trim loops failed, rendering the face as a
//! hole). A single use spanning more than three quarters of a period
//! on a doubly-periodic face is always the wrong branch; the face's
//! use pcurves are cleared (the recompute skips uses that already
//! carry one, so the wrong-branch transplant would survive untouched)
//! and rebuilt through the chart-coherent recompute (PK-286).

use popkern_core::topology::BRep;

/// Detect and rebuild wrong-branch pcurves on doubly-periodic faces.
pub(super) fn recompute_incoherent_periodic_pcurves(result: &mut BRep) {
    let mut incoherent: Vec<popkern_core::topology::FaceHandle> = Vec::new();
    for (fh, fd) in &result.faces {
        let Some(sof) = fd.surface.as_ref() else {
            continue;
        };
        use popkern_core::geometry::surface::SurfaceEval as _;
        let (pu, pv) = sof.geometry.parameter_periodicity();
        let (Some(pu), Some(pv)) = (pu, pv) else {
            continue;
        };
        let mut loops = vec![fd.outer_loop];
        loops.extend(fd.inner_loops.iter().copied());
        'face: for lh in loops {
            let Some(ld) = result.loops.get(lh) else {
                continue;
            };
            let start = ld.half_edge;
            let mut cur = start;
            for _ in 0..result.half_edges.len() + 1 {
                let Some(hd) = result.half_edges.get(cur) else {
                    break;
                };
                if let Some(pc) = hd.pcurve.as_ref() {
                    let (mut lo_u, mut hi_u, mut lo_v, mut hi_v) =
                        (f64::MAX, f64::MIN, f64::MAX, f64::MIN);
                    for k in 0..=8 {
                        let t = pc.t_start + (pc.t_end - pc.t_start) * k as f64 / 8.0;
                        let q = pc.point_at(t);
                        lo_u = lo_u.min(q.x);
                        hi_u = hi_u.max(q.x);
                        lo_v = lo_v.min(q.y);
                        hi_v = hi_v.max(q.y);
                    }
                    if hi_u - lo_u > 0.75 * pu || hi_v - lo_v > 0.75 * pv {
                        incoherent.push(fh);
                        break 'face;
                    }
                }
                cur = hd.next;
                if cur == start {
                    break;
                }
            }
        }
    }
    for fh in incoherent {
        // Clear the face's use pcurves first: the recompute skips
        // uses that already carry one, so the wrong-branch
        // transplant would survive untouched.
        let loops: Vec<popkern_core::topology::LoopHandle> = result
            .faces
            .get(fh)
            .map(|fd| {
                std::iter::once(fd.outer_loop)
                    .chain(fd.inner_loops.iter().copied())
                    .collect()
            })
            .unwrap_or_default();
        for lh in loops {
            let Some(ld) = result.loops.get(lh) else {
                continue;
            };
            let start = ld.half_edge;
            let mut cur = start;
            for _ in 0..result.half_edges.len() + 1 {
                let Some(hd) = result.half_edges.get_mut(cur) else {
                    break;
                };
                hd.pcurve = None;
                let next = hd.next;
                cur = next;
                if cur == start {
                    break;
                }
            }
        }
        let r = result.recompute_face_pcurves(fh);
        if popkern_core::trace::trace_enabled("boolean") {
            eprintln!("  [pc-coherence] recompute {fh:?}: {r:?}");
        }
    }
}

/// the cellular trace topic diagnostics: dump every edge with a non-2 use count
/// (+ position/span) from a cellular assembly that failed validation.
pub(super) fn trace_cellular_nonmanifold(r: &BRep) {
    if popkern_core::trace::trace_enabled("cellular") {
        // Dump every edge's use count + position for
        // non-manifold diagnosis.
        for (eh, edge) in r.edges.iter() {
            let mut uses = 0;
            for (_, he) in r.half_edges.iter() {
                if he.edge == eh && r.loops.contains_key(he.loop_) {
                    uses += 1;
                }
            }
            if uses != 2 {
                use popkern_core::geometry::curve::CurveEval;
                let pos = edge
                    .curve
                    .as_ref()
                    .map(|c| {
                        let p = c.geometry.point_at(0.5 * (c.t_start + c.t_end));
                        format!("({:.2},{:.2},{:.2})", p.x, p.y, p.z)
                    })
                    .unwrap_or_else(|| {
                        let p = r.vertices[r.half_edges[edge.half_edge].origin].position;
                        format!("v({:.2},{:.2},{:.2})", p.x, p.y, p.z)
                    });
                let mut endpoints = String::new();
                {
                    use std::fmt::Write as _;
                    let he0 = edge.half_edge;
                    if r.half_edges.contains_key(he0) {
                        let a = r.vertices[r.half_edges[he0].origin].position;
                        let tw = r.half_edges[he0].twin;
                        if r.half_edges.contains_key(tw) {
                            let b = r.vertices[r.half_edges[tw].origin].position;
                            let _ = write!(
                                endpoints,
                                " span ({:.1},{:.1},{:.1})->({:.1},{:.1},{:.1})",
                                a.x, a.y, a.z, b.x, b.y, b.z
                            );
                        }
                    }
                }
                eprintln!("  [nm] {eh:?} uses={uses} at {pos}{endpoints}");
            }
        }
    }
}
