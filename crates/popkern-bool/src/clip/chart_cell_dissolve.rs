//! Chart-span guard for co-surfacial face dissolution (extracted from
//! `chart_cell.rs` under the file-size ratchet).

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use super::face_chart_cell;

/// Would dissolving `edge` (merging its two incident faces) keep the
/// merged face cell-resolvable? `true` for non-doubly-periodic
/// surfaces (unchanged behaviour). For torus faces the merge
/// is admitted only while the combined chart cell stays strictly
/// inside one period per axis: dissolving every co-surfacial patch
/// produced one full-wrap mega-face whose loop closes in no UV
/// branch (measured: 31-use loop, UV gap 0.67, STEP demoted to
/// OPEN_SHELL, OCC NULL shape).
pub(crate) fn edge_dissolve_keeps_cell_resolvable(
    brep: &BRep,
    edge: popkern_core::topology::EdgeHandle,
    tolerance: &Tolerance,
) -> bool {
    let Some(edge_data) = brep.edges.get(edge) else {
        return true;
    };
    let Some(he) = brep.half_edges.get(edge_data.half_edge) else {
        return true;
    };
    let Some(twin) = brep.half_edges.get(he.twin) else {
        return true;
    };
    let (Some(loop_a), Some(loop_b)) = (brep.loops.get(he.loop_), brep.loops.get(twin.loop_))
    else {
        return true;
    };
    let face_a = loop_a.face;
    let face_b = loop_b.face;
    let Some(sof) = brep.faces.get(face_a).and_then(|f| f.surface.as_ref()) else {
        return true;
    };
    let (u_period, v_period) = sof.geometry.periodicity();
    let (Some(p_u), Some(p_v)) = (u_period, v_period) else {
        return true;
    };
    let (Some(cell_a), Some(cell_b)) = (
        face_chart_cell(brep, face_a, tolerance),
        face_chart_cell(brep, face_b, tolerance),
    ) else {
        // A winding (non-resolvable) participant must not grow.
        return false;
    };
    let align = |lo: f64, hi: f64, glo: f64, ghi: f64, period: f64| {
        let shift = ((0.5 * (glo + ghi) - 0.5 * (lo + hi)) / period).round() * period;
        (lo + shift, hi + shift)
    };
    let (blo_u, bhi_u) = align(
        cell_b.u_range.0,
        cell_b.u_range.1,
        cell_a.u_range.0,
        cell_a.u_range.1,
        p_u,
    );
    let (blo_v, bhi_v) = align(
        cell_b.v_range.0,
        cell_b.v_range.1,
        cell_a.v_range.0,
        cell_a.v_range.1,
        p_v,
    );
    let span_u = cell_a.u_range.1.max(bhi_u) - cell_a.u_range.0.min(blo_u);
    let span_v = cell_a.v_range.1.max(bhi_v) - cell_a.v_range.0.min(blo_v);
    span_u < p_u - 1.0e-9 && span_v < p_v - 1.0e-9
}
