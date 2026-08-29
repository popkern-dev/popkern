//! Slit-seam dissolution (PK-295).
//!
//! A boolean that carves a periodic face down to a fragment which no
//! longer wraps the chart period can leave the original seam edge
//! INSIDE the fragment as a slit: both half-edges wired into the same
//! loop of the same face, with their pcurves co-located in the chart
//! (a REAL seam's halves live a full period apart — that is what
//! makes the wrapped chart contiguous). The slit pinches the face's
//! UV ring into duplicate points, the boundary-first tessellator can
//! never cover the pinched segments, and the face falls back to the
//! uniform grid whose boundary samples do not match the neighbours'
//! (measured: the θ=0 torus bulge of a rim-filleted cylinder minus a
//! cuboid shipped 57 open mesh edges along its spiric and rim).
//!
//! The pcurve distance between the two halves is the discriminator:
//! co-located (≪ period) means slit — dissolving it restores the
//! plain fragment outline; a period apart means load-bearing seam —
//! untouchable. Dissolution is transactional: the loop must survive
//! as a single cycle and the body must still validate, else the edge
//! is left in place.

use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle};

/// Report from a slit-seam dissolution pass.
#[derive(Debug, Clone, Copy, Default)]
pub struct SlitSeamReport {
    /// Slit seam edges dissolved.
    pub dissolved: usize,
    /// Candidates that declined (loop would split, or validation
    /// regressed) and were rolled back.
    pub declined: usize,
}

/// Dissolve seam edges whose two halves are chart-co-located slits.
pub fn dissolve_slit_seams(brep: &mut BRep, _tolerance: &Tolerance) -> SlitSeamReport {
    let mut report = SlitSeamReport::default();
    let candidates: Vec<EdgeHandle> = brep
        .edges
        .iter()
        .filter_map(|(eh, ed)| {
            let he = brep.half_edges.get(ed.half_edge)?;
            let twin = brep.half_edges.get(he.twin)?;
            // Same loop, same face — the seam/slit shape.
            if he.loop_ != twin.loop_ {
                return None;
            }
            if popkern_core::trace::trace_enabled("slit") {
                eprintln!(
                    "  [slit-scan] {eh:?} same_loop adj={}/{} pcurves={}/{}",
                    he.next == he.twin,
                    twin.next == ed.half_edge,
                    he.pcurve.is_some(),
                    twin.pcurve.is_some()
                );
            }
            // Consecutive halves (an out-and-back flap) are the slit
            // at its narrow limit and dissolve through the same
            // splice: prev-of-first bridges to next-of-second. The
            // DEGENERATE-length flavour is the spur pass's spur; a
            // macroscopic flap like the θ=0 bulge's 0.22-long seam
            // remnant belongs here.
            let face = brep.loops.get(he.loop_)?.face;
            let surface = brep.faces.get(face)?.surface.as_ref()?;
            let (pu, pv) = surface.geometry.parameter_periodicity();
            if pu.is_none() && pv.is_none() {
                return None;
            }
            match (he.pcurve.as_ref(), twin.pcurve.as_ref()) {
                (Some(pc_a), Some(pc_b)) => {
                    let mid = |pc: &popkern_core::geometry::pcurve::Pcurve| {
                        pc.point_at(0.5 * (pc.t_start + pc.t_end))
                    };
                    let (ma, mb) = (mid(pc_a), mid(pc_b));
                    // Co-located on every periodic axis ⇒ slit. A real
                    // seam reads a whole period apart on its seam axis.
                    let co_located = |d: f64, period: Option<f64>| match period {
                        Some(p) if p > 0.0 => d.abs() < p * 0.25,
                        _ => true,
                    };
                    if popkern_core::trace::trace_enabled("slit") {
                        eprintln!(
                            "  [slit] candidate {eh:?} mids ({:.4},{:.4})/({:.4},{:.4}) periods {pu:?}/{pv:?}",
                            ma.x, ma.y, mb.x, mb.y
                        );
                    }
                    if !(co_located(ma.x - mb.x, pu) && co_located(ma.y - mb.y, pv)) {
                        return None;
                    }
                }
                _ => {
                    // No pcurves to arbitrate with: admit only a
                    // geometrically ZERO-LENGTH edge — co-located by
                    // definition (junction residue: curve-less
                    // same-vertex loop edges left by graph assembly,
                    // measured as un-exportable face uses demoting
                    // the STEP shell to OPEN_SHELL).
                    let a = brep.vertices.get(he.origin)?.position;
                    let b = brep.vertices.get(twin.origin)?.position;
                    if a.distance_to(&b) > _tolerance.linear() {
                        return None;
                    }
                }
            }
            // Both endpoints must keep OTHER incident edges: a slit
            // whose tip has no other edge is a chart-necessary
            // anchor (a sphere cap's pole edge charts exactly this
            // way), not boolean residue — dissolving it strands the
            // tip vertex and breaks the chart.
            let has_other = |v: popkern_core::topology::VertexHandle| {
                brep.half_edges
                    .iter()
                    .any(|(hh, hd)| hd.origin == v && hh != ed.half_edge && hh != he.twin)
            };
            if !has_other(he.origin) || !has_other(twin.origin) {
                return None;
            }
            Some(eh)
        })
        .collect();

    for eh in candidates {
        if !brep.edges.contains_key(eh) {
            continue;
        }
        if dissolve_one(brep, eh) {
            report.dissolved += 1;
        } else {
            report.declined += 1;
        }
    }
    report
}

/// Splice one slit edge out of its loop. Transactional: the loop must
/// remain a single cycle and the body must validate, else rollback.
fn dissolve_one(brep: &mut BRep, eh: EdgeHandle) -> bool {
    let snapshot = brep.clone();
    let Some(ed) = brep.edges.get(eh) else {
        return false;
    };
    let h = ed.half_edge;
    let Some(hd) = brep.half_edges.get(h).cloned() else {
        return false;
    };
    let t = hd.twin;
    let Some(td) = brep.half_edges.get(t).cloned() else {
        return false;
    };
    let loop_h = hd.loop_;

    // A NON-consecutive slit threads the loop through the edge twice:
    // removing it necessarily splits the cycle in two — that is kemr
    // (kill edge, make ring), the island rim becoming a proper inner
    // loop. The splice below is only sound for the consecutive
    // out-and-back flap, where both halves drop from ONE surviving
    // cycle (the old cycle check declined every keyhole slit and the
    // merged torus faces kept their grid slits un-dissolved).
    let consecutive = hd.next == t || td.next == h;
    if !consecutive {
        match brep.kemr(h) {
            Ok(_) => {
                if let Err(e) = brep.validate() {
                    if popkern_core::trace::trace_enabled("slit") {
                        eprintln!("  [slit] kemr {eh:?} validate failed: {e}");
                    }
                    *brep = snapshot;
                    return false;
                }
                return true;
            }
            Err(e) => {
                if popkern_core::trace::trace_enabled("slit") {
                    eprintln!("  [slit] kemr {eh:?} declined: {e}");
                }
                *brep = snapshot;
                return false;
            }
        }
    }

    // Members of the loop before the splice (for the cycle check).
    let mut members = 0usize;
    {
        let start = brep.loops[loop_h].half_edge;
        let mut cursor = start;
        for _ in 0..brep.half_edges.len() + 1 {
            members += 1;
            cursor = brep.half_edges[cursor].next;
            if cursor == start {
                break;
            }
        }
    }

    // Splice: bypass both halves. In the consecutive (out-and-back)
    // configuration one of hn/tn IS the other removed half — route
    // around it so every surviving pointer is live.
    let (hp, hn, tp, tn) = (hd.prev, hd.next, td.prev, td.next);
    let (live_after_h, live_after_t) =
        (if tn == h { hn } else { tn }, if hn == t { tn } else { hn });
    if hp != t {
        brep.half_edges[hp].next = live_after_h;
        brep.half_edges[live_after_h].prev = hp;
    }
    if tp != h {
        brep.half_edges[tp].next = live_after_t;
        brep.half_edges[live_after_t].prev = tp;
    }
    if brep.loops[loop_h].half_edge == h || brep.loops[loop_h].half_edge == t {
        brep.loops[loop_h].half_edge = live_after_t;
    }
    brep.remove_half_edge_reanchored(h);
    brep.remove_half_edge_reanchored(t);
    brep.edges.remove(eh);
    // Vertex outgoing repair: any surviving half-edge that actually
    // ORIGINATES at the vertex (a positional guess from the splice
    // neighbourhood picks the wrong vertex in the consecutive case).
    for v in [hd.origin, td.origin] {
        let stale = brep
            .vertices
            .get(v)
            .is_some_and(|vd| !brep.half_edges.contains_key(vd.outgoing));
        if stale {
            let replacement = brep
                .half_edges
                .iter()
                .find(|(_, cand)| cand.origin == v)
                .map(|(hh, _)| hh);
            if let (Some(r), Some(vd)) = (replacement, brep.vertices.get_mut(v)) {
                vd.outgoing = r;
            }
        }
    }

    // The splice must leave ONE cycle carrying every remaining
    // member; splitting the loop in two would orphan a cycle.
    let mut walked = 0usize;
    {
        let start = brep.loops[loop_h].half_edge;
        let mut cursor = start;
        for _ in 0..brep.half_edges.len() + 1 {
            walked += 1;
            let Some(hd) = brep.half_edges.get(cursor) else {
                walked = usize::MAX;
                break;
            };
            cursor = hd.next;
            if cursor == start {
                break;
            }
        }
    }
    if walked != members - 2 {
        if popkern_core::trace::trace_enabled("slit") {
            eprintln!(
                "  [slit] declined {eh:?}: walked {walked} != members-2 {}",
                members - 2
            );
        }
        *brep = snapshot;
        return false;
    }
    if let Err(e) = brep.validate() {
        if popkern_core::trace::trace_enabled("slit") {
            eprintln!("  [slit] declined {eh:?}: validate {e}");
        }
        *brep = snapshot;
        return false;
    }
    true
}
