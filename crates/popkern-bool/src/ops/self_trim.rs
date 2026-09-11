//! Global self-intersection trim for a single body (PK-327, design §5
//! Stage 4).
//!
//! Consumes a [`SelfIntersectionScan`]: splits both faces of every
//! detected pair along the shared curve (`self_split_along_curves` —
//! the boolean cascade's splitter driven with per-entry sides),
//! classifies every face of the split body with the caller's `keep`
//! predicate, removes the rejected region
//! ([`BRep::remove_face_region_to_boundary`]), and re-closes the rim
//! by welding the opposing boundary copies of each intersection curve
//! ([`BRep::weld_boundary_edge_pair`]).
//!
//! The split body stays CONNECTED through all of this — the E3 probe
//! outcome (recorded in the design doc): the graph-walk/assemble arm
//! of the cascade is not needed for self-trims, because splitting one
//! body along its own intersection curves never separates it; only
//! the splitter is shared with the boolean pipeline.

use std::collections::BTreeSet;

use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::CurveEval;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle};

use crate::error::BoolError;
use crate::intersect::IntersectionRecord;
use crate::ops::self_intersect::{SelfIntersectionCurve, SelfIntersectionScan};
use crate::split::self_split_along_curves;

/// Result of [`trim_self_intersections`].
#[derive(Debug)]
pub struct SelfTrimOutcome {
    /// The trimmed, re-closed body.
    pub body: BRep,
    /// Faces removed as unrealized (post-split count).
    pub removed_faces: usize,
    /// Boundary-edge pairs welded along intersection curves.
    pub welded_seams: usize,
}

/// Trim `body` along the scan's curves, keeping exactly the faces the
/// `keep` predicate accepts, and re-close the rim.
///
/// `keep` is consulted once per face of the SPLIT body (fragments
/// included) together with the pre-split face it descends from
/// (identity for unsplit faces); it owns the realization semantics
/// (for offset skins: does any sample of the face achieve its own
/// offset distance against the source body?).
///
/// # Errors
///
/// [`BoolError::SelfTrimFailed`] when the predicate rejects everything
/// (empty result), the region removal fails, or the rim does not close
/// into opposing weld pairs. [`BoolError::SplitFailed`]-class
/// diagnostics propagate from the splitter. The input body is never
/// mutated — the outcome carries the new body.
pub fn trim_self_intersections(
    body: &BRep,
    scan: &SelfIntersectionScan,
    keep: &dyn Fn(&BRep, FaceHandle, FaceHandle) -> bool,
    tolerance: &Tolerance,
) -> Result<SelfTrimOutcome, BoolError> {
    let records: Vec<IntersectionRecord> = scan
        .curves
        .iter()
        .map(|c| IntersectionRecord {
            residue_owner_a: None,
            face_a: c.face_a,
            face_b: c.face_b,
            curve: c.curve.clone(),
            is_tangent: false,
            is_coincident: false,
            start_boundary: Default::default(),
            end_boundary: Default::default(),
        })
        .collect();

    let pre_boundary: BTreeSet<EdgeHandle> = body.boundary_edges.clone();
    let (mut split_body, fragments) = self_split_along_curves(body.clone(), &records, tolerance)
        .map_err(|d| {
            BoolError::SelfTrimFailed(match d.error {
                Some(e) => format!("self-split failed: {e}"),
                None => "self-split failed".into(),
            })
        })?;
    super::self_trim_sharpen::sharpen_intersection_vertices(
        &mut split_body,
        &fragments,
        tolerance,
    )?;

    let trace = popkern_core::trace::trace_enabled("trim");
    let dead: BTreeSet<FaceHandle> = split_body
        .faces
        .keys()
        .filter(|&fh| {
            let origin = fragments.origin.get(&fh).copied().unwrap_or(fh);
            let alive = keep(&split_body, fh, origin);
            if trace {
                let n = split_body.face_outer_loop_positions(fh).len();
                let c = centroid(&split_body, fh);
                eprintln!(
                    "  [trim] face {fh:?} (origin {origin:?}, {n} verts, centroid {c:?}): {}",
                    if alive { "ALIVE" } else { "dead" }
                );
            }
            !alive
        })
        .collect();
    if dead.is_empty() {
        return Ok(SelfTrimOutcome {
            body: split_body,
            removed_faces: 0,
            welded_seams: 0,
        });
    }
    if dead.len() >= split_body.faces.len() {
        return Err(BoolError::SelfTrimFailed(
            "the whole body is unrealized — nothing would survive the trim".into(),
        ));
    }

    let removed_faces = dead.len();
    split_body
        .remove_face_region_to_boundary(&dead)
        .map_err(|e| BoolError::SelfTrimFailed(format!("region removal failed: {e}")))?;

    // Re-close: the new boundary edges are the opposing copies of the
    // intersection curves (one per split side). Pair by opposed
    // endpoints + agreeing midpoints, then weld.
    let new_boundary: Vec<EdgeHandle> = split_body
        .boundary_edges
        .iter()
        .copied()
        .filter(|eh| !pre_boundary.contains(eh))
        .collect();
    // The midpoint arm reads the surviving half-edge's PCURVE through
    // its face surface: split-authored 3D edge curves can be chord
    // polylines (their midpoint degenerates to the chord center and
    // cannot distinguish opposing halves of a circle), while the
    // pcurve is authored against the exact face chart.
    let signature = |brep: &BRep, eh: EdgeHandle| -> Option<(Point3, Point3, Point3)> {
        let he = brep.edges[eh].half_edge;
        let he_data = brep.half_edges.get(he)?;
        let start = brep.vertices[he_data.origin].position;
        let end = brep.vertices[brep.half_edges[he_data.next].origin].position;
        let pcurve_mid = he_data.pcurve.as_ref().and_then(|pc| {
            let (t0, t1) = pc.parameter_range();
            let uv = pc.point_at((t0 + t1) * 0.5);
            let loop_data = brep.loops.get(he_data.loop_)?;
            let face = brep.faces.get(loop_data.face)?;
            let sof = face.surface.as_ref()?;
            use popkern_core::geometry::surface::SurfaceEval;
            Some(sof.geometry.point_at(uv.x, uv.y))
        });
        let mid = pcurve_mid
            .or_else(|| {
                brep.edges[eh]
                    .curve
                    .as_ref()
                    .map(|c| c.geometry.point_at((c.t_start + c.t_end) * 0.5))
            })
            .unwrap_or(Point3::new(
                (start.x + end.x) * 0.5,
                (start.y + end.y) * 0.5,
                (start.z + end.z) * 0.5,
            ));
        Some((start, end, mid))
    };
    let tol = tolerance.linear() * 10.0;
    let mut welded_seams = 0usize;
    let mut unpaired: Vec<EdgeHandle> = new_boundary;
    if trace {
        for &eh in &unpaired {
            if let Some((s0, s1, sm)) = signature(&split_body, eh) {
                eprintln!("  [trim] rim {eh:?}: start {s0:?} end {s1:?} mid {sm:?}");
            }
        }
    }
    while let Some(keep_edge) = unpaired.pop() {
        if !split_body.boundary_edges.contains(&keep_edge) {
            continue; // already welded as a partner
        }
        let Some((k0, k1, km)) = signature(&split_body, keep_edge) else {
            return Err(BoolError::SelfTrimFailed(format!(
                "rim edge {keep_edge:?} has no surviving half-edge"
            )));
        };
        let partner = unpaired
            .iter()
            .copied()
            .filter(|eh| split_body.boundary_edges.contains(eh))
            .find(|&eh| {
                signature(&split_body, eh).is_some_and(|(x0, x1, xm)| {
                    (x0 - k1).length() <= tol
                        && (x1 - k0).length() <= tol
                        && (xm - km).length() <= tol
                })
            });
        let Some(kill_edge) = partner else {
            return Err(BoolError::SelfTrimFailed(format!(
                "rim edge {keep_edge:?} has no opposing weld partner — the trim region \
                 does not close (start {k0:?}, end {k1:?})"
            )));
        };
        if trace {
            eprintln!("  [trim] weld {keep_edge:?} <- {kill_edge:?}");
        }
        // Re-author the surviving seam's 3D curve from the exact scan
        // curve before welding: the splitter's provisional edge curves
        // are chord polylines (assemble re-authors them in the boolean
        // pipeline; the self-trim owns that step here). Matched by the
        // pcurve-evaluated midpoint against each scan curve.
        author_seam_curve(&mut split_body, keep_edge, (k0, k1, km), scan, tolerance);
        if let Some(sig) = signature(&split_body, kill_edge) {
            author_seam_curve(&mut split_body, kill_edge, sig, scan, tolerance);
        }
        split_body
            .weld_boundary_edge_pair(keep_edge, kill_edge, tolerance)
            .map_err(|e| BoolError::SelfTrimFailed(format!("seam weld failed: {e}")))?;
        welded_seams += 1;
    }

    // Orphan sweep: welding can strand anchor vertices whose every
    // half-edge was redirected (trace identifies producers).
    let incident: BTreeSet<_> = split_body.half_edges.values().map(|he| he.origin).collect();
    let orphaned: Vec<_> = split_body
        .vertices
        .keys()
        .filter(|vh| !incident.contains(vh))
        .collect();
    for vh in orphaned {
        if trace {
            eprintln!(
                "  [trim] orphan vertex {vh:?} at {:?}",
                split_body.vertices[vh].position
            );
        }
        split_body.vertices.remove(vh);
    }

    if trace {
        // Connected components over twin-face adjacency.
        let mut comp: std::collections::BTreeMap<FaceHandle, usize> = Default::default();
        let mut next = 0usize;
        for fh in split_body.faces.keys() {
            if comp.contains_key(&fh) {
                continue;
            }
            let mut stack = vec![fh];
            while let Some(f) = stack.pop() {
                if comp.insert(f, next).is_some() {
                    continue;
                }
                let face = &split_body.faces[f];
                for lh in std::iter::once(face.outer_loop).chain(face.inner_loops.iter().copied()) {
                    let Some(l) = split_body.loops.get(lh) else {
                        continue;
                    };
                    let start = l.half_edge;
                    let mut cur = start;
                    for _ in 0..=split_body.half_edges.len() {
                        let Some(he) = split_body.half_edges.get(cur) else {
                            break;
                        };
                        if let Some(tw) = split_body.half_edges.get(he.twin)
                            && let Some(tl) = split_body.loops.get(tw.loop_)
                            && !comp.contains_key(&tl.face)
                        {
                            stack.push(tl.face);
                        }
                        cur = he.next;
                        if cur == start {
                            break;
                        }
                    }
                }
            }
            next += 1;
        }
        eprintln!("  [trim] {} connected component(s) after welds", next);
    }

    // Batch postconditions: the whole welded body validates (welds
    // check only local consistency — intermediate slit states are
    // legal mid-batch), and the boundary set is exactly what the
    // caller started with (open-sheet rims survive; nothing leaks).
    split_body.validate().map_err(|e| {
        BoolError::SelfTrimFailed(format!(
            "welded body failed validation: {e} (V={} E={} F={} L={} HE={} shells={} boundary={})",
            split_body.vertices.len(),
            split_body.edges.len(),
            split_body.faces.len(),
            split_body.loops.len(),
            split_body.half_edges.len(),
            split_body.shells.len(),
            split_body.boundary_edges.len(),
        ))
    })?;
    if split_body.boundary_edges != pre_boundary {
        return Err(BoolError::SelfTrimFailed(format!(
            "trim changed the boundary set: {} before, {} after",
            pre_boundary.len(),
            split_body.boundary_edges.len()
        )));
    }

    Ok(SelfTrimOutcome {
        body: split_body,
        removed_faces,
        welded_seams,
    })
}

/// Face centroid of outer-loop vertex positions (trace only).
fn centroid(brep: &BRep, fh: FaceHandle) -> Point3 {
    let ps = brep.face_outer_loop_positions(fh);
    let n = ps.len().max(1) as f64;
    let mut c = Point3::new(0.0, 0.0, 0.0);
    for p in &ps {
        c.x += p.x / n;
        c.y += p.y / n;
        c.z += p.z / n;
    }
    c
}

/// Replace `eh`'s provisional curve with the exact scan curve whose
/// locus contains the edge (best midpoint match within 10× linear),
/// re-trimmed to the edge's endpoint vertices.
fn author_seam_curve(
    body: &mut BRep,
    eh: EdgeHandle,
    (start, end, mid): (Point3, Point3, Point3),
    scan: &SelfIntersectionScan,
    tolerance: &Tolerance,
) {
    let tol = tolerance.linear() * 10.0;
    // Select by IN-DOMAIN sampling: Newton closest-point can converge
    // to parameters beyond a clipped curve's declared range (the
    // clipped NURBS still evaluates the underlying full circle there),
    // which made the wrong half of a split circle claim zero distance.
    let mut best: Option<(f64, &SelfIntersectionCurve)> = None;
    for c in &scan.curves {
        let (a0, a1) = c.curve.curve_3d.parameter_range();
        let mut d = f64::INFINITY;
        for i in 0..=32 {
            let t = a0 + (a1 - a0) * (i as f64) / 32.0;
            d = d.min(c.curve.curve_3d.point_at(t).distance_to(&mid));
        }
        if best.as_ref().is_none_or(|(bd, _)| d < *bd) {
            best = Some((d, c));
        }
    }
    let trace_on = popkern_core::trace::trace_enabled("trim");
    let Some((d, chosen)) = best else {
        if trace_on {
            eprintln!("  [trim] author {eh:?}: no scan curve candidates");
        }
        return;
    };
    // Sampled-min gate: 32 samples over the span bound the chord gap;
    // the authored range is still verified pointwise below.
    let select_gate = tol.max(curve_span_estimate(chosen) * 0.01);
    if d > select_gate {
        if trace_on {
            eprintln!("  [trim] author {eh:?}: best scan match {d:.3e} > {select_gate:.3e}");
        }
        return; // no exact author found: keep the provisional curve
    }
    let geometry =
        popkern_core::geometry::curve::CurveGeometry::Nurbs(chosen.curve.curve_3d.clone());
    let (Some(t_s0), Some(t_e0)) = (geometry.closest_point(&start), geometry.closest_point(&end))
    else {
        if trace_on {
            eprintln!("  [trim] author {eh:?}: endpoint projection failed");
        }
        return;
    };
    if geometry.point_at(t_s0).distance_to(&start) > tol
        || geometry.point_at(t_e0).distance_to(&end) > tol
    {
        if trace_on {
            eprintln!(
                "  [trim] author {eh:?}: endpoint residuals {:.3e}/{:.3e}",
                geometry.point_at(t_s0).distance_to(&start),
                geometry.point_at(t_e0).distance_to(&end)
            );
        }
        return;
    }
    // Closed-curve branch disambiguation: both halves of a circle
    // share their endpoints, so the projected range can land on the
    // wrong branch, and backward traversal needs an inverted span.
    // Try every in-domain (± one period) endpoint combination and
    // keep the range whose midpoint matches the edge's
    // pcurve-evaluated midpoint. Shifted parameters outside the
    // curve's domain are skipped (NURBS evaluation clamps there).
    let (c0, c1) = chosen.curve.curve_3d.parameter_range();
    let period = c1 - c0;
    let in_domain = |t: f64| t >= c0 - 1e-12 && t <= c1 + 1e-12;
    let (mut t_s, mut t_e) = (t_s0, t_e0);
    let mut best_mid = f64::INFINITY;
    for cand_s in [t_s0, t_s0 + period, t_s0 - period] {
        if !in_domain(cand_s) || geometry.point_at(cand_s).distance_to(&start) > tol {
            continue;
        }
        for cand_e in [t_e0, t_e0 + period, t_e0 - period] {
            if !in_domain(cand_e) || geometry.point_at(cand_e).distance_to(&end) > tol {
                continue;
            }
            let m = geometry.point_at((cand_s + cand_e) * 0.5);
            let dm = m.distance_to(&mid);
            if dm < best_mid {
                best_mid = dm;
                t_s = cand_s;
                t_e = cand_e;
            }
        }
    }
    if best_mid > tol * 100.0 {
        if trace_on {
            eprintln!(
                "  [trim] author {eh:?}: no range matches pcurve mid (best {best_mid:.3e});                  t_s0={t_s0:.4} t_e0={t_e0:.4} domain=[{c0:.4},{c1:.4}] start={start:?} end={end:?}"
            );
            for cand_s in [t_s0, t_s0 + period, t_s0 - period] {
                let ds = geometry.point_at(cand_s).distance_to(&start);
                eprintln!(
                    "    cand_s {cand_s:.4}: residual {ds:.3e} in_domain={}",
                    in_domain(cand_s)
                );
            }
        }
        return; // no range matches the pcurve: keep the provisional curve
    }
    if popkern_core::trace::trace_enabled("trim") {
        let m = geometry.point_at((t_s + t_e) * 0.5);
        eprintln!(
            "  [trim] author {eh:?}: t_s={t_s:.4} t_e={t_e:.4} curve-mid {m:?} (pcurve mid {mid:?})"
        );
    }
    if let Some(edge) = body.edges.get_mut(eh) {
        edge.curve = Some(popkern_core::geometry::curve::CurveOnEdge {
            geometry,
            t_start: t_s,
            t_end: t_e,
        });
    }
}

/// Rough 3D span of a scan curve (chord sum over 8 samples).
fn curve_span_estimate(c: &SelfIntersectionCurve) -> f64 {
    let (a0, a1) = c.curve.curve_3d.parameter_range();
    let mut len = 0.0;
    let mut prev = c.curve.curve_3d.point_at(a0);
    for i in 1..=8 {
        let t = a0 + (a1 - a0) * (i as f64) / 8.0;
        let p = c.curve.curve_3d.point_at(t);
        len += p.distance_to(&prev);
        prev = p;
    }
    len
}
