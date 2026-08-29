//! Periodic UV-loop branch normalization for STEP export.
//!
//! STEP readers reconstruct a face's trim contour by chaining its
//! boundary pcurves in UV. On a periodic surface the kernel's stored
//! pcurves can sit in unrelated period branches — the two uses of a
//! seam edge carry the SAME pcurve, and opposite rims can occupy
//! windows a full period apart. The reconstructed contour then jumps
//! by ±period at every branch break, never closes, and OpenCASCADE's
//! trim of the face degenerates (measured on the flagship boss: the
//! seam-fillet torus and boss cylinder sewed to their neighbours only
//! along the chart seam line, gapping everywhere else).
//!
//! This pass walks each loop of every face on a periodic basis
//! surface and computes a per-use constant shift (k·period per
//! periodic direction) that chains each pcurve's walk-start onto the
//! previous pcurve's walk-end. The second use of a seam edge lands
//! one period away from the first by construction — producing exactly
//! the branch-separated pcurve pair STEP requires for SEAM_CURVEs.
//! Shifts are applied at pcurve emission; kernel pcurves are never
//! modified (the kernel's own consumers unwrap rings themselves).

use std::collections::HashMap;

use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::topology::{BRep, FaceHandle, HalfEdgeHandle, LoopHandle};

use crate::step::export_pcurve_synthesis::synthesize_export_half_edge_pcurve_on_surface;
use crate::step::export_topology::ExportedFaceSurface;

/// Per-face-use pcurve shift in (u, v), each a whole multiple of the
/// respective period.
pub(super) type BranchShifts = HashMap<(FaceHandle, HalfEdgeHandle), (f64, f64)>;

/// Per-face-use pcurve traversal direction chosen by the loop chain
/// (true = the walk traverses the pcurve from t_start to t_end).
pub(super) type ChainDirections = HashMap<(FaceHandle, HalfEdgeHandle), bool>;

/// One candidate loop-chain assignment, in loop-walk order: the UV
/// branch shift `(du, dv)` applied to each use and whether the use is
/// traversed forward (t_start → t_end).
type ChainAssignment = Vec<((f64, f64), bool)>;

/// Result of the periodic loop normalization pre-pass.
#[derive(Default)]
pub(crate) struct PeriodicNormalization {
    pub shifts: BranchShifts,
    pub directions: ChainDirections,
}

/// The pcurve the exporter will emit for a face use: the stored
/// face-use override or half-edge pcurve when it validates against
/// the edge, otherwise a synthesized projection. This is the SAME
/// selection `try_export_he_pcurve` performs — the shift pre-pass
/// must see the pcurves that will actually be written.
pub(super) fn export_pcurve_for_use(
    brep: &BRep,
    face: FaceHandle,
    half_edge: HalfEdgeHandle,
    surface: &SurfaceGeometry,
) -> Option<Pcurve> {
    use crate::step::export_pcurve_synthesis::{
        export_pcurve_trace_matches_edge, step_export_pcurve_association_tolerance,
    };
    use crate::step::export_topology_helpers::pcurve_endpoints_match_edge_segment;
    let he = brep.half_edges.get(half_edge)?;
    let edge = he.edge;
    let tolerance = step_export_pcurve_association_tolerance(brep);
    let matches = |pcurve: &Pcurve| {
        pcurve_endpoints_match_edge_segment(brep, edge, pcurve, surface, tolerance)
            || export_pcurve_trace_matches_edge(brep, edge, half_edge, pcurve, surface, tolerance)
    };
    let stored_raw = brep
        .face_use_pcurve(face, half_edge)
        .or(he.pcurve.as_ref())
        .cloned();
    let stored = stored_raw.clone().filter(|pcurve| matches(pcurve));
    let result = stored
        .or_else(|| synthesize_export_half_edge_pcurve_on_surface(brep, edge, half_edge, surface))
        .filter(|pcurve| matches(pcurve))
        // Final fallback: the faithful UNGATED projection of the
        // edge curve. Its residual trace distance is the edge's
        // honest serving-geometry debt (e.g. a corner chord bowing
        // off the wall mid-span) — sharing the edge with an
        // approximate pcurve keeps the shell SEWABLE and lets the
        // reader's SameParameter healing absorb the gap, where the
        // face-local fork guarantees once-used edges (measured on
        // the boss: 6 closed free-bound wires and an un-normalized
        // periodic wall loop that OCC read as self-intersecting).
        .or_else(|| {
            crate::step::export_pcurve_synthesis::synthesize_half_edge_pcurve_ungated(
                brep, edge, half_edge, surface,
            )
        });
    if result.is_none() && popkern_core::trace::trace_enabled("step-io") {
        use crate::step::export_pcurve_synthesis::sampled_edge_surface_disagreement;
        let stored_dev = stored_raw.as_ref().and_then(|pc| {
            sampled_edge_surface_disagreement(brep, edge, half_edge, surface, pc, 9)
        });
        eprintln!(
            "  [step-export] no exportable pcurve for use {half_edge:?} (edge {edge:?}): \
             stored_present={} stored_dev={stored_dev:?}, tol={tolerance:.3e}",
            stored_raw.is_some(),
        );
    }
    result
}

/// Best-effort pcurve for CHAINING only: when no candidate passes
/// the association gate (residual serving-geometry debt — e.g. a
/// corner chord bowing off its support mid-span), the branch
/// normalization still needs a UV image to chain the loop through;
/// an approximately-placed image beats abandoning the whole loop to
/// mixed branches (measured: the boss wall exported un-normalized
/// and OCC read a self-intersecting wire). Emission stays gated —
/// this pcurve is never written.
pub(super) fn chain_pcurve_for_use(
    brep: &BRep,
    face: FaceHandle,
    half_edge: HalfEdgeHandle,
    surface: &SurfaceGeometry,
) -> Option<Pcurve> {
    if let Some(pc) = export_pcurve_for_use(brep, face, half_edge, surface) {
        return Some(pc);
    }
    let he = brep.half_edges.get(half_edge)?;
    crate::step::export_pcurve_synthesis::synthesize_half_edge_pcurve_ungated(
        brep, he.edge, half_edge, surface,
    )
    .or_else(|| {
        brep.face_use_pcurve(face, half_edge)
            .or(he.pcurve.as_ref())
            .cloned()
    })
}

/// True when the pcurve's trace runs OPPOSITE to the edge's exported
/// 3D curve. STEP associated pcurves must be co-parameterized with
/// the basis curve (the EDGE_CURVE's same_sense flag relates the
/// edge to that shared parameterization); popkern stores pcurves in
/// HALF-EDGE direction, so one side of every shared edge is reversed
/// and must be flipped at emission. OpenCASCADE silently heals
/// reversed pcurves by reprojection on simple surfaces, but on
/// periodic surfaces the reprojection is branch-ambiguous — the
/// flagship seam fillet sewed only along the chart seam line.
pub(super) fn pcurve_opposes_edge_curve(
    brep: &BRep,
    half_edge: HalfEdgeHandle,
    surface: &SurfaceGeometry,
    pcurve: &Pcurve,
) -> bool {
    use popkern_core::geometry::CurveEval;
    let Some(he) = brep.half_edges.get(half_edge) else {
        return false;
    };
    let Some(curve) = brep.edges.get(he.edge).and_then(|e| e.curve.as_ref()) else {
        return false;
    };
    let ceval: &dyn CurveEval = &curve.geometry;
    let early_3d = ceval.point_at(curve.t_start + 0.05 * (curve.t_end - curve.t_start));
    let pc_at = |f: f64| {
        let uv = pcurve.point_at(pcurve.t_start + f * (pcurve.t_end - pcurve.t_start));
        surface.point_at(uv.x, uv.y)
    };
    pc_at(0.05).distance_to(&early_3d) > pc_at(0.95).distance_to(&early_3d)
}

/// Reverse a NURBS pcurve in place: control points/weights reversed,
/// knots mirrored about the parameter range.
pub(super) fn reverse_pcurve_nurbs(nurbs: &mut popkern_core::geometry::NurbsCurve) {
    nurbs.control_points.reverse();
    nurbs.weights.reverse();
    let lo = nurbs.knots.first().copied().unwrap_or(0.0);
    let hi = nurbs.knots.last().copied().unwrap_or(1.0);
    let sum = lo + hi;
    nurbs.knots.reverse();
    for k in &mut nurbs.knots {
        *k = sum - *k;
    }
}

/// Compute branch shifts + chain directions for every face on a
/// periodic basis surface.
pub(super) fn compute_periodic_branch_shifts(
    brep: &BRep,
    surface_map: &HashMap<FaceHandle, ExportedFaceSurface>,
) -> PeriodicNormalization {
    let mut norm = PeriodicNormalization::default();
    for (fh, face) in &brep.faces {
        let Some(exported) = surface_map.get(&fh) else {
            continue;
        };
        let geometry = &exported.surface.geometry;
        let (period_u, period_v) = geometry.parameter_periodicity();
        if period_u.is_none() && period_v.is_none() {
            continue;
        }
        let face_sense = exported.surface.orientation;
        let mut loops = vec![(face.outer_loop, true)];
        loops.extend(face.inner_loops.iter().map(|&lh| (lh, false)));
        for (lh, is_outer) in loops {
            // Required UV winding: outer contours run CCW when the
            // face sense is true (material left of the traversal in
            // the chart), CW when flipped; inner contours opposite.
            let want_ccw = face_sense == is_outer;
            normalize_loop(
                brep,
                fh,
                lh,
                geometry,
                (period_u, period_v),
                want_ccw,
                &mut norm.shifts,
                &mut norm.directions,
            );
        }
    }
    norm
}

/// Walk one loop, chaining each use's pcurve onto the previous one.
/// Direction choices for closed uses (both pcurve endpoints at the
/// same 3D point) are resolved by exact loop CLOSURE: greedy
/// residuals tie for closed rims (either direction chains after a
/// period shift), but only one assignment returns the cursor to the
/// walk start — the other leaves a ±k·period gap. All ambiguous
/// combinations are enumerated (loops carry at most a few closed
/// uses) and the closing one wins. The chosen per-use traversal
/// direction is recorded in `directions` — it is the loop-consistent
/// (manifold) direction `edge_same_sense` needs for closed rims,
/// where the raw pcurve orientation contract is unenforced.
#[allow(clippy::too_many_arguments)]
fn normalize_loop(
    brep: &BRep,
    face: FaceHandle,
    lh: LoopHandle,
    geometry: &SurfaceGeometry,
    periods: (Option<f64>, Option<f64>),
    want_ccw: bool,
    shifts: &mut BranchShifts,
    directions: &mut ChainDirections,
) {
    let Some(loop_) = brep.loops.get(lh) else {
        return;
    };
    let start_he = loop_.half_edge;
    if !brep.half_edges.contains_key(start_he) {
        return;
    }

    // Collect the loop's uses with their pcurve UV endpoints.
    struct Use {
        half_edge: HalfEdgeHandle,
        // Both pcurve endpoints (walk direction resolved during
        // chaining — closed rims are direction-ambiguous from
        // geometry alone).
        end_a: (f64, f64),
        end_b: (f64, f64),
        // end_a is the walk start when known from 3D (open edges).
        a_is_start: Option<bool>,
        // Sampled UV polyline in pcurve parameter order (for the
        // winding sign of candidate assignments).
        samples: Vec<(f64, f64)>,
    }
    let mut uses: Vec<Use> = Vec::new();
    let mut cursor_he = start_he;
    for _ in 0..brep.half_edges.len() + 1 {
        let Some(he) = brep.half_edges.get(cursor_he) else {
            return;
        };
        let Some(pc) = chain_pcurve_for_use(brep, face, cursor_he, geometry) else {
            // A use without a pcurve breaks the chain; leave this
            // loop un-normalized rather than guessing across the gap.
            if popkern_core::trace::trace_enabled("step-io") {
                eprintln!(
                    "  [step-export] periodic loop {lh:?} on {face:?} UN-NORMALIZED: use \
                     {cursor_he:?} has no exportable pcurve"
                );
            }
            return;
        };
        let uv0 = pc.point_at(pc.t_start);
        let uv1 = pc.point_at(pc.t_end);
        let p0 = geometry.point_at(uv0.x, uv0.y);
        let p1 = geometry.point_at(uv1.x, uv1.y);
        let origin = brep
            .vertices
            .get(he.origin)
            .map(|v| v.position)
            .unwrap_or(p0);
        let d0 = p0.distance_to(&origin);
        let d1 = p1.distance_to(&origin);
        // Open edges resolve their walk direction from 3D; closed
        // rims (both endpoints at the same 3D point) stay ambiguous
        // and are resolved against the chain cursor below.
        let a_is_start = if p0.distance_to(&p1) > 1e-9 {
            Some(d0 <= d1)
        } else {
            None
        };
        let samples: Vec<(f64, f64)> = (0..=8)
            .map(|k| {
                let t = pc.t_start + (pc.t_end - pc.t_start) * k as f64 / 8.0;
                let uv = pc.point_at(t);
                (uv.x, uv.y)
            })
            .collect();
        uses.push(Use {
            half_edge: cursor_he,
            end_a: (uv0.x, uv0.y),
            end_b: (uv1.x, uv1.y),
            a_is_start,
            samples,
        });
        cursor_he = he.next;
        if cursor_he == start_he {
            break;
        }
    }
    if uses.len() < 2 {
        // Single-use loops (a full rim bounding a hole) have nothing
        // to chain against; either direction covers the rim.
        return;
    }

    let shift_for = |from: (f64, f64), to: (f64, f64)| -> (f64, f64) {
        let du = periods
            .0
            .map(|p| ((to.0 - from.0) / p).round() * p)
            .unwrap_or(0.0);
        let dv = periods
            .1
            .map(|p| ((to.1 - from.1) / p).round() * p)
            .unwrap_or(0.0);
        (du, dv)
    };
    let dist = |a: (f64, f64), b: (f64, f64)| ((a.0 - b.0).powi(2) + (a.1 - b.1).powi(2)).sqrt();

    // Ambiguous (closed) uses get both direction candidates; open
    // uses are pinned by 3D. Enumerate assignments, keep the one
    // whose chain closes exactly.
    let ambiguous: Vec<usize> = uses
        .iter()
        .enumerate()
        .filter(|(_, u)| u.a_is_start.is_none())
        .map(|(i, _)| i)
        .collect();
    let combos = 1usize << ambiguous.len().min(12);
    // score = (chain error, flips from stored, winding mismatch) —
    // lexicographic. Closure first. Then FEWEST reversals of the
    // stored per-use pcurve directions: when the stored directions
    // themselves chain into an exactly-closed contour they ARE the
    // manifold traversal (the face-use pcurve contract), and an
    // annulus's two-fold symmetry must NOT be re-decided by chart
    // winding — a valid torus band exports with its outer contour
    // wound clockwise (orientation-true boss round-over, measured:
    // the assumed CCW-iff-face-sense convention flipped the top
    // rim's chain direction, the shared rim's EDGE_CURVE flag came
    // out .T. instead of .F., and OCC rejected the torus face with
    // SelfIntersectingWire + UnorientableShape on every corpus
    // curved-seam part). Winding stays as the LAST tiebreak for
    // loops whose stored directions do not close on their own.
    let mut best: Option<(f64, u32, bool, ChainAssignment)> = None;
    for mask in 0..combos {
        let forward_of = |idx: usize| -> bool {
            match uses[idx].a_is_start {
                Some(fwd) => fwd,
                None => {
                    let pos = ambiguous.iter().position(|&k| k == idx).unwrap_or(0);
                    mask & (1 << pos) == 0
                }
            }
        };
        let (s0, e0) = if forward_of(0) {
            (uses[0].end_a, uses[0].end_b)
        } else {
            (uses[0].end_b, uses[0].end_a)
        };
        let mut cursor = e0;
        let mut assignment: Vec<((f64, f64), bool)> = vec![((0.0, 0.0), forward_of(0))];
        let mut drift = 0.0_f64;
        for (idx, use_) in uses.iter().enumerate().skip(1) {
            let fwd = forward_of(idx);
            let (start, end) = if fwd {
                (use_.end_a, use_.end_b)
            } else {
                (use_.end_b, use_.end_a)
            };
            let (du, dv) = shift_for(start, cursor);
            drift = drift.max(dist((start.0 + du, start.1 + dv), cursor));
            assignment.push(((du, dv), fwd));
            cursor = (end.0 + du, end.1 + dv);
        }
        let closure = dist(cursor, s0);
        let score = closure.max(drift);
        // Shoelace over the traversed contour.
        let mut area2 = 0.0_f64;
        let mut prev: Option<(f64, f64)> = None;
        let mut first: Option<(f64, f64)> = None;
        for (use_, ((du, dv), fwd)) in uses.iter().zip(assignment.iter()) {
            let mut pts = use_.samples.clone();
            if !*fwd {
                pts.reverse();
            }
            for p in pts {
                let q = (p.0 + du, p.1 + dv);
                if let Some(pr) = prev {
                    area2 += pr.0 * q.1 - q.0 * pr.1;
                }
                if first.is_none() {
                    first = Some(q);
                }
                prev = Some(q);
            }
        }
        if let (Some(pr), Some(f0)) = (prev, first) {
            area2 += pr.0 * f0.1 - f0.0 * pr.1;
        }
        let winding_ok = (area2 > 0.0) == want_ccw;
        let flips = (mask as u32).count_ones();
        let better = match &best {
            None => true,
            Some((b_score, b_flips, b_wind, _)) => {
                (score, flips, !winding_ok) < (*b_score, *b_flips, !*b_wind)
            }
        };
        if better {
            best = Some((score, flips, winding_ok, assignment));
        }
        // The zero-flip assignment is unbeatable once it closes
        // exactly (flips outranks winding).
        if score < 1e-9 && mask == 0 {
            break;
        }
    }
    let Some((score, _flips, _winding_ok, assignment)) = best else {
        return;
    };
    if popkern_core::trace::trace_enabled("step-io") {
        eprintln!(
            "  [step-export] periodic loop {lh:?} on {face:?}: {} use(s), score {score:.3e}, \
             flips {_flips}, winding_ok {_winding_ok}, want_ccw {want_ccw}",
            uses.len()
        );
    }
    for (use_, ((du, dv), fwd)) in uses.iter().zip(assignment.iter()) {
        if *du != 0.0 || *dv != 0.0 {
            shifts.insert((face, use_.half_edge), (*du, *dv));
        }
        directions.insert((face, use_.half_edge), *fwd);
    }
    if score > 1e-6 && popkern_core::trace::trace_enabled("step-io") {
        eprintln!(
            "  [step-export] periodic loop on face {face:?} does not close in UV after branch normalization (gap {score:.3e})"
        );
        // Re-walk the winning assignment and print per-use spans +
        // chain gaps so the defective pcurve is identifiable.
        let mut cursor: Option<(f64, f64)> = None;
        for (use_, ((du, dv), fwd)) in uses.iter().zip(assignment.iter()) {
            let (start, end) = if *fwd {
                (use_.end_a, use_.end_b)
            } else {
                (use_.end_b, use_.end_a)
            };
            let s = (start.0 + du, start.1 + dv);
            let e = (end.0 + du, end.1 + dv);
            let gap = cursor.map(|c| ((s.0 - c.0).powi(2) + (s.1 - c.1).powi(2)).sqrt());
            eprintln!(
                "      [chain] he {:?} fwd={fwd} ({:.3},{:.3})->({:.3},{:.3}) gap {gap:?}",
                use_.half_edge, s.0, s.1, e.0, e.1
            );
            cursor = Some(e);
        }
    }
}
