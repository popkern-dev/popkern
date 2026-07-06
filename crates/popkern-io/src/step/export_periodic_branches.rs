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
    let stored = brep
        .face_use_pcurve(face, half_edge)
        .or(he.pcurve.as_ref())
        .cloned()
        .filter(|pcurve| matches(pcurve));
    stored
        .or_else(|| synthesize_export_half_edge_pcurve_on_surface(brep, edge, half_edge, surface))
        .filter(|pcurve| matches(pcurve))
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

/// Compute branch shifts for every face on a periodic basis surface.
pub(super) fn compute_periodic_branch_shifts(
    brep: &BRep,
    surface_map: &HashMap<FaceHandle, ExportedFaceSurface>,
) -> BranchShifts {
    let mut shifts = BranchShifts::new();
    for (fh, face) in &brep.faces {
        let Some(exported) = surface_map.get(&fh) else {
            continue;
        };
        let geometry = &exported.surface.geometry;
        let (period_u, period_v) = geometry.parameter_periodicity();
        if period_u.is_none() && period_v.is_none() {
            continue;
        }
        let mut loops = vec![face.outer_loop];
        loops.extend(face.inner_loops.iter().copied());
        for lh in loops {
            normalize_loop(brep, fh, lh, geometry, (period_u, period_v), &mut shifts);
        }
    }
    shifts
}

/// Walk one loop, chaining each use's pcurve onto the previous one.
fn normalize_loop(
    brep: &BRep,
    face: FaceHandle,
    lh: LoopHandle,
    geometry: &SurfaceGeometry,
    periods: (Option<f64>, Option<f64>),
    shifts: &mut BranchShifts,
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
    }
    let mut uses: Vec<Use> = Vec::new();
    let mut cursor_he = start_he;
    for _ in 0..brep.half_edges.len() + 1 {
        let Some(he) = brep.half_edges.get(cursor_he) else {
            return;
        };
        let Some(pc) = export_pcurve_for_use(brep, face, cursor_he, geometry) else {
            // A use without a pcurve breaks the chain; leave this
            // loop un-normalized rather than guessing across the gap.
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
        uses.push(Use {
            half_edge: cursor_he,
            end_a: (uv0.x, uv0.y),
            end_b: (uv1.x, uv1.y),
            a_is_start,
        });
        cursor_he = he.next;
        if cursor_he == start_he {
            break;
        }
    }
    if uses.len() < 2 {
        // Single-use loops (a full rim bounding a hole) have nothing
        // to chain against.
        return;
    }

    // Anchor on the first use; chain the rest.
    let first_start = match uses[0].a_is_start {
        Some(true) | None => uses[0].end_a,
        Some(false) => uses[0].end_b,
    };
    let mut cursor = match uses[0].a_is_start {
        Some(true) | None => uses[0].end_b,
        Some(false) => uses[0].end_a,
    };
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

    for use_ in uses.iter().skip(1) {
        // Try both endpoint assignments (respecting 3D direction
        // when known) and keep the one whose shifted start lands
        // nearest the cursor.
        let candidates: Vec<((f64, f64), (f64, f64))> = match use_.a_is_start {
            Some(true) => vec![(use_.end_a, use_.end_b)],
            Some(false) => vec![(use_.end_b, use_.end_a)],
            None => vec![(use_.end_a, use_.end_b), (use_.end_b, use_.end_a)],
        };
        let mut best: Option<((f64, f64), (f64, f64), f64)> = None;
        for (start, end) in candidates {
            let (du, dv) = shift_for(start, cursor);
            let shifted_start = (start.0 + du, start.1 + dv);
            let residual = dist(shifted_start, cursor);
            if best.map(|(_, _, r)| residual < r).unwrap_or(true) {
                best = Some(((du, dv), end, residual));
            }
        }
        let Some(((du, dv), end, _)) = best else {
            continue;
        };
        if du != 0.0 || dv != 0.0 {
            shifts.insert((face, use_.half_edge), (du, dv));
        }
        cursor = (end.0 + du, end.1 + dv);
    }

    if std::env::var_os("PK_STEP_TRACE").is_some() {
        let closure = dist(cursor, first_start);
        if closure > 1e-6 {
            eprintln!(
                "  [step-export] periodic loop on face {face:?} does not close in UV after branch normalization (gap {closure:.3e})"
            );
        }
    }
}
