//! PK-141: expel foreign slit POCKETS from strip loops.
//!
//! A strip that temporarily absorbed neighbouring support territory
//! and later released the area can keep a zero-width boundary
//! excursion — a sub-path of its loop that leaves through territory
//! far off the strip surface and returns to the same vertex
//! (measured: the +y wall strip carried the wall's upper SEAM plus
//! the ENTIRE TOP RIM as such a pocket, putting a 0.7
//! curve-surface deviation on the rim edge). The pocket's edges
//! belong to the neighbour whose surface actually contains them;
//! this pass splices the pocket out of the strip's loop and into
//! that neighbour's loop at the shared anchor vertex.
//!
//! Detection: a maximal run of loop half-edges whose interior
//! vertices all sit farther than the pocket gate off the strip
//! surface, starting and ending at the same anchor vertex. The
//! receiving face is the one (other than the strip) whose loop
//! passes the anchor and whose surface fits the pocket's vertices.
//! Splice is pure pointer surgery: `dV = dE = dF = dL = 0`.

use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{Point3, SurfaceEval};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle, HalfEdgeHandle};

/// Expel pockets from every NURBS strip loop. Returns spliced runs.
pub(crate) fn expel_strip_pockets(brep: &mut BRep, tolerance: &Tolerance) -> usize {
    let trace = std::env::var_os("PK_JUNCTION_TRACE").is_some();
    let strips: Vec<FaceHandle> = brep
        .faces
        .iter()
        .filter(|(_, f)| {
            matches!(
                f.surface.as_ref().map(|s| &s.geometry),
                Some(SurfaceGeometry::Nurbs(_))
            )
        })
        .map(|(h, _)| h)
        .collect();
    let mut expelled = 0usize;
    for strip in strips {
        while let Some(n) = expel_one_pocket(brep, strip, tolerance, trace) {
            expelled += n;
        }
    }
    expelled
}

fn off_face(brep: &BRep, face: FaceHandle, p: &Point3) -> Option<f64> {
    let sof = brep.faces.get(face)?.surface.as_ref()?;
    let (u, v) = popkern_core::geometry::uv_projection::project_point_to_uv(&sof.geometry, p)?;
    Some(p.distance_to(&sof.geometry.point_at(u, v)))
}

/// Find and splice ONE pocket on `strip`; None when the loop is
/// pocket-free (loop-until-fixpoint caller).
fn expel_one_pocket(
    brep: &mut BRep,
    strip: FaceHandle,
    tolerance: &Tolerance,
    trace: bool,
) -> Option<usize> {
    let fd = brep.faces.get(strip)?;
    let lh = fd.outer_loop;
    let start = brep.loops.get(lh)?.half_edge;
    // Pocket gate: an order of magnitude beyond any honest boundary
    // staleness (the settle classes live at <= 5e-2; the measured
    // pocket vertices sit 0.7 off).
    let gate = (tolerance.linear().max(1e-9) * 1e5).min(0.3);
    // Walk the loop once, recording each half-edge and whether its
    // ORIGIN is far off the strip.
    let mut hes: Vec<(HalfEdgeHandle, bool)> = Vec::new();
    let mut cur = start;
    for _ in 0..brep.half_edges.len() + 1 {
        let he = brep.half_edges.get(cur)?;
        let pos = brep.vertices.get(he.origin)?.position;
        let far = off_face(brep, strip, &pos)
            .map(|d| d > gate)
            .unwrap_or(false);
        hes.push((cur, far));
        cur = he.next;
        if cur == start {
            break;
        }
    }
    let n = hes.len();
    if n < 4 {
        return None;
    }
    // A pocket = maximal run of half-edges [i..j] where every
    // interior origin is far, the bounding origins (i's origin and
    // j+1's origin) are near, and those bounding origins are the
    // SAME vertex.
    for i in 0..n {
        if hes[i].1 || !hes[(i + 1) % n].1 {
            continue;
        }
        // Run starts at i+1 (first far origin); extend while far.
        let mut j = (i + 1) % n;
        let mut len = 0usize;
        while hes[j].1 && len <= n {
            j = (j + 1) % n;
            len += 1;
        }
        if len == 0 || len + 2 >= n {
            continue;
        }
        // The pocket half-edges: from hes[i] (leaves the near vertex
        // into the far run) through the run, back to hes[j] whose
        // origin... splice boundary: entry = hes[i], exit = the
        // half-edge ENDING at the anchor = predecessor of hes[j].
        let entry = hes[i].0;
        let exit = hes[(j + n - 1) % n].0;
        let anchor = brep.half_edges.get(entry)?.origin;
        let after = hes[j].0;
        let ret_origin = brep.half_edges.get(after)?.origin;
        if anchor != ret_origin {
            continue;
        }
        // Receiving face: another face whose loop passes the anchor
        // and whose surface CONTAINS the pocket vertices.
        let probe_vertex = brep.half_edges.get(hes[(i + 1) % n].0)?.origin;
        let probe = brep.vertices.get(probe_vertex)?.position;
        let mut receiver: Option<(FaceHandle, HalfEdgeHandle)> = None;
        for (h, he) in brep.half_edges.iter() {
            if he.origin != anchor || he.loop_ == brep.half_edges.get(entry)?.loop_ {
                continue;
            }
            let f = brep.loops.get(he.loop_)?.face;
            if off_face(brep, f, &probe)
                .map(|d| d < gate * 0.1)
                .unwrap_or(false)
            {
                receiver = Some((f, h));
                break;
            }
        }
        let Some((rface, rhe)) = receiver else {
            continue;
        };
        // Splice: remove [entry..exit] from the strip loop; insert
        // before `rhe` in the receiver loop.
        let strip_loop = brep.half_edges.get(entry)?.loop_;
        let prev_entry = brep.half_edges.get(entry)?.prev;
        brep.half_edges[prev_entry].next = after;
        brep.half_edges[after].prev = prev_entry;
        if brep.loops.get(strip_loop).map(|l| l.half_edge) == Some(entry) || {
            let mut inside = false;
            let mut c = entry;
            loop {
                if Some(c) == brep.loops.get(strip_loop).map(|l| l.half_edge) {
                    inside = true;
                }
                if c == exit {
                    break;
                }
                c = brep.half_edges[c].next;
            }
            inside
        } {
            brep.loops[strip_loop].half_edge = after;
        }
        let rloop = brep.half_edges[rhe].loop_;
        let rprev = brep.half_edges[rhe].prev;
        brep.half_edges[rprev].next = entry;
        brep.half_edges[entry].prev = rprev;
        brep.half_edges[exit].next = rhe;
        brep.half_edges[rhe].prev = exit;
        let mut c = entry;
        loop {
            brep.half_edges[c].loop_ = rloop;
            if c == exit {
                break;
            }
            c = brep.half_edges[c].next;
        }
        if trace {
            eprintln!(
                "[pocket] expelled {len}-edge pocket from {strip:?} into {rface:?} at {anchor:?}"
            );
        }
        return Some(1);
    }
    None
}
