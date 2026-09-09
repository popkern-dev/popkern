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

/// Detect and rebuild wrong-branch pcurves on doubly-periodic faces,
/// and IMAGE-DISHONEST pcurves on any periodic face: a use pcurve
/// whose fraction-matched image runs away from its own 3-D edge curve
/// (the nut's top mouth-chamfer ring transplanted a piece sweeping
/// azimuth 0→3.9 the MAJOR way while the curve covers 3.9→2π; the
/// loop then read chart-continuous with zero net winding, the seam
/// sides co-branched, and the annulus meshed torn — PK-464). Clearing
/// and recomputing routes the piece through the branch-honest chained
/// projection.
pub(super) fn recompute_incoherent_periodic_pcurves(result: &mut BRep) {
    use popkern_core::geometry::curve::CurveEval as _;
    let mut incoherent: Vec<popkern_core::topology::FaceHandle> = Vec::new();
    for (fh, fd) in &result.faces {
        let Some(sof) = fd.surface.as_ref() else {
            continue;
        };
        use popkern_core::geometry::surface::SurfaceEval as _;
        let (pu, pv) = sof.geometry.parameter_periodicity();
        if pu.is_none() && pv.is_none() {
            continue;
        }
        let doubly = pu.is_some() && pv.is_some();
        let image_gate = result.tolerance.linear().max(1.0e-9) * 100.0;
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
                    if doubly {
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
                        if let (Some(pu), Some(pv)) = (pu, pv)
                            && (hi_u - lo_u > 0.75 * pu || hi_v - lo_v > 0.75 * pv)
                        {
                            incoherent.push(fh);
                            break 'face;
                        }
                    }
                    // Fraction-matched image honesty against the 3-D
                    // edge curve (direction resolved by the endpoint
                    // nearer the half-edge origin).
                    if let Some(c) = result.edges.get(hd.edge).and_then(|e| e.curve.as_ref())
                        && let Some(v) = result.vertices.get(hd.origin)
                    {
                        let cs = c.geometry.point_at(c.t_start);
                        let ce = c.geometry.point_at(c.t_end);
                        let (ta, tb) = if v.position.distance_to(&cs) <= v.position.distance_to(&ce)
                        {
                            (c.t_start, c.t_end)
                        } else {
                            (c.t_end, c.t_start)
                        };
                        let mut worst = 0.0_f64;
                        for k in 1..4 {
                            let f = k as f64 / 4.0;
                            let q = pc.point_at(pc.t_start + (pc.t_end - pc.t_start) * f);
                            let s = sof.geometry.point_at(q.x, q.y);
                            let cp = c.geometry.point_at(ta + (tb - ta) * f);
                            worst = worst.max(s.distance_to(&cp));
                        }
                        if worst > image_gate {
                            if popkern_core::trace::trace_enabled("boolean") {
                                eprintln!(
                                    "  [pc-coherence] {fh:?} use pcurve image {worst:.4} off its edge curve — rebuilding face"
                                );
                            }
                            incoherent.push(fh);
                            break 'face;
                        }
                    }
                }
                cur = hd.next;
                if cur == start {
                    break;
                }
            }
        }
    }
    let mut rebuilt = false;
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
        rebuilt = true;
    }
    // The orientation flag is derived from loop winding × surface,
    // and `reconcile_orientation_flags` ran BEFORE this sweep — a
    // rebuilt face's recomputed pcurves can wind opposite the
    // assembly-installed ones the reconcile judged (measured: the
    // torus fillet strip of a filleted-cylinder union read
    // winding-dishonest on validate_all, signed area +π² against
    // flag false). Re-reconcile so the served flags match the served
    // pcurves.
    if rebuilt {
        let flips = result.reconcile_orientation_flags();
        if flips > 0 && popkern_core::trace::trace_enabled("boolean") {
            eprintln!("  [pc-coherence] re-reconciled {flips} orientation flag(s)");
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
        // Coincident-vertex census: cellular stub edges are minted
        // between DUPLICATE nodes the graph failed to unify — the
        // duplicate pair is the root, the stub only the symptom.
        let verts: Vec<_> = r.vertices.iter().collect();
        for i in 0..verts.len() {
            for j in (i + 1)..verts.len() {
                let d = verts[i].1.position.distance_to(&verts[j].1.position);
                if d < 1.0e-4 {
                    let p = verts[i].1.position;
                    eprintln!(
                        "  [nm-dup] {:?} ~ {:?} d={d:.3e} at ({:.7},{:.7},{:.7})",
                        verts[i].0, verts[j].0, p.x, p.y, p.z
                    );
                }
            }
        }
    }
}

/// Loop chart continuity on closed charts (PK-441): assembly adoption
/// and projection recompute cannot pick a seam BRANCH, and a trim
/// loop with both seam sides on one border carries chart
/// discontinuities into every consumer (the doubly-ground coil's wall
/// meshed with hundreds of open edges). Whole-period shifts only —
/// exact on geometrically closed charts.
pub(super) fn normalize_output_loop_charts(result: &mut BRep) {
    let faces: Vec<_> = result.faces.keys().collect();
    result.normalize_loop_chart_continuity(&faces);
}
