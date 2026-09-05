//! Pcurve preservation through graph edge splits (PK-441).
//!
//! `split_graph_edge` used to DROP both sides' pcurves ("regenerated
//! downstream"), letting chart embeddings re-project the trimmed 3D
//! curve. For an edge lying ON the seam of a geometrically-closed
//! chart, projection cannot recover the seam BRANCH (u = 0 vs u = 1
//! image the same 3D curve): the re-projected halves landed on
//! inconsistent branches, the chart-walk substrate declined on a
//! degree imbalance at the seam junction, and the face fell back to
//! the 3D-Newell walk — whose tangent-plane angles cannot separate a
//! grazing corner's near-parallel directions (the double-grind
//! spring's 5.6e-4 rim × seam corner stranded the assembly).
//!
//! The cure is PK-391's, one layer up: cut each side's pcurve at the
//! parameter whose IMAGE is the split point, and hand each half-edge
//! its honest sub-span. On any failure (no face, no surface, image
//! never reaches the point, endpoint cut) the affected side falls
//! back to dropping its pcurves — the pre-existing behavior.

use popkern_core::geometry::Point3;
use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::geometry::pcurve_split::split_pcurve_at_image;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::topology::BRep;

use super::graph::IntersectionGraph;

/// Whether a half-edge pair's stored 3D curve runs in the FORWARD
/// half's traversal direction, resolved through the pcurve when the
/// endpoint test cannot tell (a CLOSED edge — a boolean-promoted rim
/// circle — has both curve endpoints at the same vertex, and the
/// split's head/tail windows then went to the wrong sides: the wall
/// side's traversal runs parameter-DECREASING on the ground coil's
/// rim, so the tie-broken head handed the wall a piece traversed
/// backwards, the boundary shoelace flipped sign, and the seam
/// borders mis-assigned — PK-441). Samples one step inside each: the
/// pcurve is half-edge-aligned by convention, so its image agreeing
/// with the curve's start side means co-directed. Returns `None` when
/// no pcurve/surface is available to consult.
pub(super) fn curve_co_directed_with_fwd(
    graph: &IntersectionGraph,
    he_idx: usize,
    brep_a: &BRep,
    brep_b: &BRep,
) -> Option<bool> {
    let twin_idx = graph.half_edges[he_idx].twin;
    let curve = graph.half_edges[he_idx]
        .curve
        .as_ref()
        .or_else(|| graph.half_edges[twin_idx].curve.as_ref())?;
    use popkern_core::geometry::curve::CurveEval as _;
    let (t0, t1) = (curve.t_start, curve.t_end);
    if !(t0.is_finite() && t1.is_finite()) || (t1 - t0).abs() <= 0.0 {
        return None;
    }
    let p_near_start = curve.geometry.point_at(t0 + (t1 - t0) * 0.05);
    // Resolve through the SOURCE BREP half-edge: the graph's fwd half
    // was inserted from `edge.half_edge`, so that brep half-edge's
    // traversal IS the fwd direction. Its pcurve direction is only
    // trustworthy when it is traversal-aligned — proved by LOOP
    // CONTINUITY (its chart end meets its successor's chart start;
    // edge-aligned transplants break the chain and decline the
    // vote — an untrusted vote swapped the head/tail windows of the
    // flush arm ∪ hub rims and stranded that union).
    use super::graph::EdgeSource;
    let (brep, eh) = match &graph.half_edges[he_idx].source {
        EdgeSource::OriginalA(eh) | EdgeSource::SplitA { original: eh } => (brep_a, *eh),
        EdgeSource::OriginalB(eh) | EdgeSource::SplitB { original: eh } => (brep_b, *eh),
        _ => return None,
    };
    let edge = brep.edges.get(eh)?;
    let he_a = brep.half_edges.get(edge.half_edge)?;
    let pc = he_a.pcurve.as_ref()?;
    let next = brep.half_edges.get(he_a.next)?;
    let pc_next = next.pcurve.as_ref()?;
    let e = pc.point_at(pc.t_end);
    let s2 = pc_next.point_at(pc_next.t_start);
    let chart_scale = {
        let a = pc.point_at(pc.t_start);
        ((e.x - a.x).powi(2) + (e.y - a.y).powi(2)).sqrt().max(1e-6)
    };
    if ((e.x - s2.x).powi(2) + (e.y - s2.y).powi(2)).sqrt() > chart_scale * 1e-3 {
        return None;
    }
    let face = brep.loops.get(he_a.loop_)?.face;
    let surface = &brep.faces.get(face)?.surface.as_ref()?.geometry;
    let span = pc.t_end - pc.t_start;
    if span.abs() <= 0.0 {
        return None;
    }
    let image = |t: f64| {
        let uv = pc.point_at(t);
        surface.point_at(uv.x, uv.y)
    };
    // Compare SWEEP DIRECTIONS at a common 3D point, never parameter
    // origins: a closed curve and its pcurve routinely parameterize
    // from different seam origins (the coil rim's exact-arc window
    // starts at −π/2 while its chart starts at u = 0), so
    // start-proximity reads the origin OFFSET, not the direction.
    let t_mid = 0.5 * (t0 + t1);
    let p_mid = curve.geometry.point_at(t_mid);
    let tan_curve = curve.geometry.point_at(t_mid + (t1 - t0) * 1e-3) - p_mid;
    let mut best = (f64::INFINITY, pc.t_start);
    for k in 0..=64 {
        let t = pc.t_start + span * k as f64 / 64.0;
        let d = image(t).distance_to(&p_mid);
        if d < best.0 {
            best = (d, t);
        }
    }
    // The pcurve must actually pass the reference point.
    let curve_len_scale = p_near_start.distance_to(&p_mid).max(1e-9);
    if best.0 > curve_len_scale * 0.2 {
        return None;
    }
    let dt = span * 1e-3;
    let t_lo = (best.1 - dt).max(pc.t_start.min(pc.t_end));
    let t_hi = (best.1 + dt).min(pc.t_start.max(pc.t_end));
    let tan_pc = image(t_hi) - image(t_lo);
    let dot = tan_pc.dot(&tan_curve) * span.signum();
    if !dot.is_finite() || dot.abs() <= 1e-30 {
        return None;
    }
    let verdict = Some(dot > 0.0);
    if popkern_core::trace::trace_enabled("face-split") {
        eprintln!(
            "[tie-break] he {he_idx} src {:?} verdict {verdict:?}",
            graph.half_edges[he_idx].source
        );
    }
    verdict
}

/// Split one side's pcurve at the split point's image and return
/// `(origin_side, far_side)` pieces, or `None` to drop (fallback).
///
/// Assignment prefers the ORIGIN-side 3D curve midpoint
/// (`origin_mid`) when the caller has one: the piece whose image
/// midpoint is nearest it belongs to the origin side. Midpoints are
/// unambiguous even for a CLOSED edge (a full rim circle), where the
/// outer ENDPOINTS of both pieces image to the same vertex and an
/// endpoint test ties — the tie broke wrong on the ground coil's rim,
/// the mis-assigned piece failed the embedding's curve verification,
/// and the rim fell back to a projection whose branch choice broke
/// the boundary's chart closure (σ flipped, seam borders followed).
/// `origin_pos` remains the fallback discriminator for curveless
/// splits.
fn split_one_side(
    pcurve: &Pcurve,
    surface: &SurfaceGeometry,
    split_point: &Point3,
    origin_pos: &Point3,
    origin_mid: Option<&Point3>,
    tolerance: f64,
) -> Option<(Pcurve, Pcurve)> {
    let (first, second) = split_pcurve_at_image(pcurve, surface, split_point, tolerance)?;
    let image = |pc: &Pcurve, t: f64| {
        let uv = pc.point_at(t);
        surface.point_at(uv.x, uv.y)
    };
    if let Some(mid) = origin_mid {
        let first_mid = image(&first, 0.5 * (first.t_start + first.t_end));
        let second_mid = image(&second, 0.5 * (second.t_start + second.t_end));
        return if first_mid.distance_to(mid) <= second_mid.distance_to(mid) {
            Some((first, second))
        } else {
            Some((second, first))
        };
    }
    let first_outer = image(&first, first.t_start);
    let second_outer = image(&second, second.t_end);
    if first_outer.distance_to(origin_pos) <= second_outer.distance_to(origin_pos) {
        Some((first, second))
    } else {
        Some((second, first))
    }
}

/// Distribute the split pcurve pieces of both sides of a just-split
/// graph edge. `he_idx`/`twin_idx` are the HEAD pair
/// (origin → split → origin), `new_fwd_idx`/`new_twin_idx` the TAIL
/// pair (split → end → split). `origin_pos`/`end_pos` are the node
/// positions from BEFORE the split. Sides whose pcurve cannot be
/// image-split are left dropped.
#[allow(clippy::too_many_arguments)]
pub(super) fn distribute_split_pcurves(
    graph: &mut IntersectionGraph,
    he_idx: usize,
    twin_idx: usize,
    new_fwd_idx: usize,
    new_twin_idx: usize,
    fwd_pcurve: Option<Pcurve>,
    twin_pcurve: Option<Pcurve>,
    split_point: &Point3,
    origin_pos: &Point3,
    end_pos: &Point3,
    brep_a: &BRep,
    brep_b: &BRep,
    tolerance: f64,
) {
    // 3D midpoints of the two curve pieces (the HEAD keeps the
    // origin-side window, the TAIL the far window) — the unambiguous
    // assignment references for closed edges.
    use popkern_core::geometry::curve::CurveEval as _;
    let head_mid = graph.half_edges[he_idx]
        .curve
        .as_ref()
        .map(|c| c.geometry.point_at(0.5 * (c.t_start + c.t_end)));
    let tail_mid = graph.half_edges[new_fwd_idx]
        .curve
        .as_ref()
        .map(|c| c.geometry.point_at(0.5 * (c.t_start + c.t_end)));
    // The image-split acceptance band must NOT scale with a coarse
    // body tolerance: the split point is an exact graph node and the
    // pcurve's image reaches it within the pcurve's own fidelity. At
    // linear = 1e-2 the 10× gate accepted cuts 0.1 off the node, the
    // distributed pieces were garbage, and a sphere − cube at that
    // tolerance walked every fragment into the discard set and served
    // an EMPTY body where the honest closure decline used to fire.
    let tolerance = tolerance.min(1.0e-4);
    // Scope: preserved pieces serve the CLOSED-NURBS-CHART class,
    // where projection cannot recover a seam branch. Analytic charts
    // keep the historical drop-and-reproject — their periodicity is
    // evaluable and the downstream branch machinery is tuned to
    // projections (installing pieces on a sphere − cube at coarse
    // tolerance walked every fragment into the discard set and
    // served an EMPTY body).
    let nurbs_surface_of = |idx: usize| -> Option<SurfaceGeometry> {
        let face = graph.half_edges[idx].face?;
        let (s, _) = super::walk::chart_subdivision::face_surface(&face, brep_a, brep_b)?;
        matches!(s, SurfaceGeometry::Nurbs(_)).then(|| s.clone())
    };
    let fwd_surface = nurbs_surface_of(he_idx);
    let twin_surface = nurbs_surface_of(twin_idx);
    // Forward side: origin-side piece → head fwd, far piece → tail fwd.
    if let (Some(pc), Some(surface)) = (fwd_pcurve, fwd_surface)
        && let Some((head, tail)) = split_one_side(
            &pc,
            &surface,
            split_point,
            origin_pos,
            head_mid.as_ref(),
            tolerance,
        )
    {
        graph.half_edges[he_idx].pcurve = Some(head);
        graph.half_edges[new_fwd_idx].pcurve = Some(tail);
    }
    // Twin side: its origin is the END node; origin-side piece → tail
    // twin (end → split), far piece → head twin (split → origin).
    if let (Some(pc), Some(surface)) = (twin_pcurve, twin_surface)
        && let Some((end_side, origin_side)) = split_one_side(
            &pc,
            &surface,
            split_point,
            end_pos,
            tail_mid.as_ref(),
            tolerance,
        )
    {
        graph.half_edges[new_twin_idx].pcurve = Some(end_side);
        graph.half_edges[twin_idx].pcurve = Some(origin_side);
    }
}
