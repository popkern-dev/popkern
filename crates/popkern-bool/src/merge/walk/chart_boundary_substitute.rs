//! Boundary-substitute direction resolution for coincident booleans.
//!
//! The coincident weld path can REPLACE pieces of a face's outer
//! boundary with overlap intersection records: the original boundary
//! half-edge is gone from the face's registration and a Shared pair
//! spans the gap. Both record directions register on the face, but
//! only one bounds anything — the other's material side lies off-face
//! (its region belongs to the coplanar neighbour, not this face).
//!
//! Degree balance alone cannot decide which direction to keep. On the
//! flush L-union both of these balanced arithmetically:
//!   - drop one direction of the INTERIOR cut (breaks the two-region
//!     partition — a cut must serve one direction per region), and
//!   - drop the off-face direction of each boundary-lying pair (the
//!     correct resolution);
//!
//! the bounded search took whichever it enumerated first, and on the
//! fully-contained coincident faces it found NOTHING because its
//! candidate seeding stops at initially-imbalanced nodes while the
//! surplus twin ring propagates its imbalance through balanced ones.
//! A transversal tool wall flush on ALL FOUR edges (the through-slot
//! class) has no imbalance at all — both full rings are balanced —
//! and the surplus ring served as a duplicate wall fragment.
//!
//! The resolution here is decided by construction instead:
//!   1. A Shared pair is a BOUNDARY SUBSTITUTE iff its chart run lies
//!      ON the face's ORIGINAL outer boundary (sampled from the
//!      operand BRep, independent of what the weld left registered)
//!      and is not already covered by a registered boundary edge.
//!      Pairs interior to the original boundary are genuine cuts and
//!      keep both directions (one per region of the subdivision).
//!   2. The kept direction is the one whose traversal AGREES with the
//!      original outer loop's traversal along that segment: stored
//!      outer loops wind material-left by the BRep winding invariant
//!      (census E1), so direction agreement IS material-left — local,
//!      anchor-free, and exact even when every boundary edge of the
//!      face was welded away (the through-slot tool wall).
//!
//! Scope: non-periodic charts in coincident booleans; outer AND
//! inner original loops (a flush cavity union welds the pocket-hole
//! rings of annulus faces — inner-loop substitutes whose surplus
//! direction walked as a garbage hole-ring fragment sampled in the
//! void).

use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::uv_projection::project_point_to_uv;
use popkern_core::topology::BRep;

use super::chart_subdivision::{FaceChartSubdivision, FilteredRecord};
use crate::merge::graph::{EdgeSource, FaceRef, IntersectionGraph};

/// Boundary samples per original-loop edge for the UV polygon.
const POLY_SAMPLES_PER_EDGE: usize = 4;

/// Resolve boundary-substitute pair directions on `sub` in place.
///
/// Caller gates on non-empty overlap pairs (the coincident class) and
/// a non-periodic chart — the flush/coplanar weld configurations.
pub(super) fn resolve_boundary_substitutes(
    graph: &IntersectionGraph,
    face: &FaceRef,
    sub: &mut FaceChartSubdivision,
    surface: &SurfaceGeometry,
    original_brep: &BRep,
) {
    // Shared pairs with BOTH directions still registered.
    let pairs: Vec<(usize, usize)> = sub
        .embeddings
        .keys()
        .copied()
        .filter(|&h| {
            let tw = graph.half_edges[h].twin;
            h < tw
                && sub.embeddings.contains_key(&tw)
                && matches!(graph.half_edges[h].source, EdgeSource::Shared { .. })
        })
        .map(|h| (h, graph.half_edges[h].twin))
        .collect();
    if pairs.is_empty() {
        return;
    }

    // The face's ORIGINAL boundary loops as closed UV polygons, in
    // stored traversal order (material-left by the winding census —
    // inner loops wind opposite the outer, which is exactly the
    // traversal a substitute along a HOLE boundary must agree with;
    // the pocket-hole rings of a flush cavity union live here).
    let Some(original_face) = original_brep.faces.get(face.face) else {
        return;
    };
    let loops: Vec<popkern_core::topology::LoopHandle> = std::iter::once(original_face.outer_loop)
        .chain(original_face.inner_loops.iter().copied())
        .collect();
    let polys: Vec<Vec<(f64, f64)>> = loops
        .iter()
        .map(|&lh| {
            original_brep
                .loop_sample_positions(lh, POLY_SAMPLES_PER_EDGE)
                .iter()
                .filter_map(|p| project_point_to_uv(surface, p))
                .collect::<Vec<_>>()
        })
        .filter(|poly| poly.len() >= 3)
        .collect();
    if polys.is_empty() {
        return;
    }
    if std::env::var_os("PK_BSUB_TRACE").is_some() {
        for (i, poly) in polys.iter().enumerate() {
            eprintln!(
                "      [bsub] {face:?} poly {i}: {} pts {:?}…",
                poly.len(),
                &poly[..poly.len().min(6)]
            );
        }
    }
    let (mut lo, mut hi) = ((f64::MAX, f64::MAX), (f64::MIN, f64::MIN));
    for &(u, v) in polys.iter().flatten() {
        lo = (lo.0.min(u), lo.1.min(v));
        hi = (hi.0.max(u), hi.1.max(v));
    }
    let diag = ((hi.0 - lo.0).powi(2) + (hi.1 - lo.1).powi(2)).sqrt();
    let tol = 1e-7 * (1.0 + diag);
    let seg_dist2 = |p: (f64, f64), a: (f64, f64), b: (f64, f64)| -> f64 {
        let (dx, dy) = (b.0 - a.0, b.1 - a.1);
        let l2 = (dx * dx + dy * dy).max(1e-300);
        let t = (((p.0 - a.0) * dx + (p.1 - a.1) * dy) / l2).clamp(0.0, 1.0);
        (p.0 - a.0 - dx * t).powi(2) + (p.1 - a.1 - dy * t).powi(2)
    };
    // Nearest boundary segment across all loops (each closed):
    // distance² and the segment's directed endpoints.
    let nearest_seg = |p: (f64, f64)| -> (f64, (f64, f64), (f64, f64)) {
        let mut best = (f64::MAX, (0.0, 0.0), (0.0, 0.0));
        for poly in &polys {
            for i in 0..poly.len() {
                let (a, b) = (poly[i], poly[(i + 1) % poly.len()]);
                let d2 = seg_dist2(p, a, b);
                if d2 < best.0 {
                    best = (d2, a, b);
                }
            }
        }
        best
    };

    // Registered boundary runs (non-record embeddings): a record run
    // already covered by one of these is the partial-cover class the
    // I5 dedup owns, not a substitute for a missing piece.
    let registered_runs: Vec<Vec<(f64, f64)>> = sub
        .embeddings
        .iter()
        .filter(|&(&h, _)| !matches!(graph.half_edges[h].source, EdgeSource::Shared { .. }))
        .map(|(_, e)| e.uv.clone())
        .collect();
    let on_registered = |p: (f64, f64)| -> bool {
        registered_runs
            .iter()
            .any(|run| run.windows(2).any(|w| seg_dist2(p, w[0], w[1]) < tol * tol))
    };

    let mut resolved: Vec<(usize, usize)> = Vec::new(); // (kept, dropped)
    for (h, tw) in pairs {
        // Classification and direction agreement vote per interior
        // sample window (endpoints sit at junction nodes shared with
        // registered edges by construction; a record spanning a
        // polygon corner votes per-window, majority wins).
        let uv = &sub.embeddings[&h].uv;
        let trace = std::env::var_os("PK_BSUB_TRACE").is_some();
        if uv.len() < 3 {
            if trace {
                eprintln!(
                    "      [bsub] {face:?} pair ({h},{tw}): run too short ({})",
                    uv.len()
                );
            }
            continue;
        }
        let interior = &uv[1..uv.len() - 1];
        let on_poly = interior.iter().all(|&p| nearest_seg(p).0 < tol * tol);
        let covered = interior.iter().all(|&p| on_registered(p));
        if !on_poly || covered {
            if trace {
                eprintln!(
                    "      [bsub] {face:?} pair ({h},{tw}): on_poly={on_poly} covered={covered} — kept as cut"
                );
            }
            continue; // interior cut, or the I5 partial-cover class
        }
        let mut agree = 0i64;
        for w in uv.windows(2) {
            let mid = (0.5 * (w[0].0 + w[1].0), 0.5 * (w[0].1 + w[1].1));
            let (d2, a, b) = nearest_seg(mid);
            if d2 >= tol * tol {
                continue;
            }
            let dot = (w[1].0 - w[0].0) * (b.0 - a.0) + (w[1].1 - w[0].1) * (b.1 - a.1);
            agree += dot.signum() as i64;
        }
        if trace {
            eprintln!("      [bsub] {face:?} pair ({h},{tw}): substitute, agree={agree}");
        }
        match agree.cmp(&0) {
            std::cmp::Ordering::Greater => resolved.push((h, tw)),
            std::cmp::Ordering::Less => resolved.push((tw, h)),
            // Indecisive vote (degenerate or corner-straddling run):
            // leave the pair registered — the balance search and the
            // 3D-Newell walk remain its recorded fallback, and the I4
            // degree check reports the node if it stays unbalanced.
            std::cmp::Ordering::Equal => {}
        }
    }

    for (kept, dropped) in resolved {
        if sub.embeddings.remove(&dropped).is_some() {
            sub.report
                .filtered
                .push(FilteredRecord::BoundarySubstituteDirection { dropped, kept });
        }
    }
}
