//! Helper functions for topology manipulation during blend construction.
//!
//! Provides perpendicular edge splitting, consumed-face identification,
//! and fillet-face merging used by the main `apply_blend_geometry` orchestrator.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{NurbsCurve, Point3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, VertexHandle};

#[path = "perp_anchor_positions.rs"]
mod perp_anchor_positions;
use perp_anchor_positions::compute_perp_anchor_positions_for_face;

pub use super::merge::{identify_consumed_face_any, merge_consumed_faces_any};

/// Split perpendicular edges at the fillet offset distance.
///
/// `brep_snapshot` is an immutable clone of the BRep captured *before*
/// Phase A began. This function reads the fillet edge's topology and
/// perpendicular boundary edges from the snapshot — so even when an
/// earlier iteration of Phase A has already split one of them, the
/// anchor *positions* are computed from the original un-fragmented
/// geometry. Vertex insertion then happens in the mutable `brep` via
/// `insert_anchor_at_position`, which locates the descendant sub-edge
/// currently containing each 3D anchor point.
///
/// `radius_at_start` is used for splits at `v_start` (the fillet
/// edge's start vertex) and `radius_at_end` for splits at `v_end`.
/// For a constant-radius spec, both are equal; for a variable
/// spec, they're the per-end values from `FilletRadius::Variable`.
/// Round 5k (continued): the construct phase now honours the full
/// variable profile rather than collapsing it to a single average.
///
/// Round 5r+: `edge_splits` in the returned `SplitInfo` is now
/// populated and consumed by `apply_blend_geometry` to maintain a
/// per-fillet-edge descendants map — Phase C identifies a consumed
/// face by matching *any* descendant sub-edge of the canonical
/// fillet edge handle.
///
/// §28.H Session 2 — backwards-compatible wrapper that delegates
/// to [`split_perpendicular_edges_by_radius_with_overrides`] with
/// `None` overrides.  Retained as a public no-override entry
/// point for tests; the production caller in `apply_blend_geometry`
/// uses the `_with_overrides` variant directly so it can pass
/// planner-derived overrides when `use_direct_apex_construction` is
/// set.
#[allow(clippy::too_many_arguments)]
#[allow(dead_code)]
pub fn split_perpendicular_edges_by_radius(
    brep_snapshot: &BRep,
    brep: &mut BRep,
    face_1: FaceHandle,
    face_2: FaceHandle,
    fillet_edge: EdgeHandle,
    trim_3d_face_1: Option<&NurbsCurve>,
    trim_3d_face_2: Option<&NurbsCurve>,
    anchor_on_trim: bool,
    radius_at_start: f64,
    radius_at_end: f64,
    tolerance: &Tolerance,
) -> Result<SplitInfo, String> {
    split_perpendicular_edges_by_radius_with_overrides(
        brep_snapshot,
        brep,
        face_1,
        face_2,
        fillet_edge,
        trim_3d_face_1,
        trim_3d_face_2,
        anchor_on_trim,
        None,
        radius_at_start,
        radius_at_end,
        tolerance,
        None,
    )
}

/// §28.H Session 2 — per-anchor-position override.  When the DirectApex
/// planner pre-pass has emitted `apex_vertex` targets for a fillet
/// edge's V-end and W-end, the override map carries those targets
/// down to Phase A's anchor placement.
///
/// Each field is the planner's preferred position for the anchor
/// that perp-walk math would otherwise compute.  When `None`, the
/// perp-walk-math position is used.  When `Some` but the position
/// is interior to a face (no existing edge passes through it),
/// `insert_anchor_at_position` falls back to perp-walk math —
/// Session 3 will introduce the topology surgery (mev/mef) needed
/// to land interior anchors.
#[derive(Debug, Clone, Default)]
pub struct PerpAnchorOverride {
    /// Override for the anchor at the fillet edge's V_start side
    /// on `face_1`.
    pub v_start_face_1: Option<Point3>,
    /// Override for the anchor at V_start side on `face_2`.
    pub v_start_face_2: Option<Point3>,
    /// Override for V_end side on `face_1`.
    pub v_end_face_1: Option<Point3>,
    /// Override for V_end side on `face_2`.
    pub v_end_face_2: Option<Point3>,
}

/// §28.H Session 2 — main implementation taking optional planner
/// overrides per anchor slot.  See [`split_perpendicular_edges_by_
/// radius`] for the field-by-field doc; this variant differs only
/// in the `override_` parameter.
#[allow(clippy::too_many_arguments)]
pub fn split_perpendicular_edges_by_radius_with_overrides(
    brep_snapshot: &BRep,
    brep: &mut BRep,
    face_1: FaceHandle,
    face_2: FaceHandle,
    fillet_edge: EdgeHandle,
    trim_3d_face_1: Option<&NurbsCurve>,
    trim_3d_face_2: Option<&NurbsCurve>,
    anchor_on_trim: bool,
    bevel_plane: Option<(Point3, popkern_core::geometry::Vec3)>,
    radius_at_start: f64,
    radius_at_end: f64,
    tolerance: &Tolerance,
    override_: Option<&PerpAnchorOverride>,
) -> Result<SplitInfo, String> {
    let edge_data = brep_snapshot
        .edges
        .get(fillet_edge)
        .ok_or_else(|| format!("fillet edge {fillet_edge:?} not found in snapshot"))?;
    let he_fwd = edge_data.half_edge;
    let he_fwd_data = brep_snapshot
        .half_edges
        .get(he_fwd)
        .ok_or_else(|| "half-edge not found".to_string())?;
    let he_rev = he_fwd_data.twin;
    let v_start = he_fwd_data.origin;
    let v_end = brep_snapshot
        .half_edges
        .get(he_rev)
        .ok_or_else(|| "twin not found".to_string())?
        .origin;

    if v_start == v_end {
        return Err(format!(
            "closed fillet edges (v_start == v_end on edge {fillet_edge:?}) are not yet \
             supported by the construct phase — the cylindrical-hole rim case requires a \
             separate non-anchored construction path"
        ));
    }

    let perp_pos_1 = compute_perp_anchor_positions_for_face(
        brep_snapshot,
        face_1,
        fillet_edge,
        v_start,
        v_end,
        radius_at_start,
        radius_at_end,
        trim_3d_face_1,
        anchor_on_trim,
        bevel_plane,
        tolerance,
    )?;
    let perp_pos_2 = compute_perp_anchor_positions_for_face(
        brep_snapshot,
        face_2,
        fillet_edge,
        v_start,
        v_end,
        radius_at_start,
        radius_at_end,
        trim_3d_face_2,
        anchor_on_trim,
        bevel_plane,
        tolerance,
    )?;

    // §28.H Session 2 — apply planner overrides per anchor slot.
    // Each `Some` override replaces the perp-walk-math position;
    // `None` keeps the perp-walk-math default.  Caller passes
    // overrides keyed by V_start vs V_end and face_1 vs face_2.
    let pos_1: (Point3, Point3) = match override_ {
        Some(o) => (
            o.v_start_face_1.unwrap_or(perp_pos_1.0),
            o.v_end_face_1.unwrap_or(perp_pos_1.1),
        ),
        None => perp_pos_1,
    };
    let pos_2: (Point3, Point3) = match override_ {
        Some(o) => (
            o.v_start_face_2.unwrap_or(perp_pos_2.0),
            o.v_end_face_2.unwrap_or(perp_pos_2.1),
        ),
        None => perp_pos_2,
    };

    let mut edge_splits: Vec<(EdgeHandle, EdgeHandle, EdgeHandle)> = Vec::new();
    // §28.H Session 2 — when override-provided position lands
    // outside any existing edge (interior to a face),
    // `insert_anchor_at_position` returns Err.  Fall back to the
    // perp-walk-math position so legacy/uniform fixtures stay
    // unbroken.  Session 3 will introduce the topology surgery
    // (mev/mef) needed to land genuine interior anchors.
    let v1_start = insert_anchor_at_position(brep, pos_1.0, &mut edge_splits)
        .or_else(|_| insert_anchor_at_position(brep, perp_pos_1.0, &mut edge_splits))?;
    let v1_end = insert_anchor_at_position(brep, pos_1.1, &mut edge_splits)
        .or_else(|_| insert_anchor_at_position(brep, perp_pos_1.1, &mut edge_splits))?;
    let v2_start = insert_anchor_at_position(brep, pos_2.0, &mut edge_splits)
        .or_else(|_| insert_anchor_at_position(brep, perp_pos_2.0, &mut edge_splits))?;
    let v2_end = insert_anchor_at_position(brep, pos_2.1, &mut edge_splits)
        .or_else(|_| insert_anchor_at_position(brep, perp_pos_2.1, &mut edge_splits))?;

    let p1_a = brep.vertices[v1_start].position;
    let p1_b = brep.vertices[v1_end].position;
    let p2_a = brep.vertices[v2_start].position;
    let p2_b = brep.vertices[v2_end].position;

    let trim_3d_1 = line_nurbs(&p1_a, &p1_b);
    let trim_3d_2 = line_nurbs(&p2_a, &p2_b);

    Ok(SplitInfo {
        trim_3d_1,
        trim_3d_2,
        face_1_anchors: (v1_start, v1_end),
        face_2_anchors: (v2_start, v2_end),
        edge_splits,
    })
}

/// Scan all edges in `brep` and find the one whose 3D line segment
/// contains `point` within a small tolerance, returning the edge
/// handle and parameter along the edge in `(0.001, 0.999)`.
///
/// `pub(crate)` for §6c-DirectApex.2.b — the apex builder calls this from
/// `resolve_apex_vertices`'s create branch when no existing BRep
/// vertex matches the planner's target position; the edge passing
/// through the target is then split via [`BRep::split_edge`] to
/// insert a new vertex on it.  Production-grade callers in this
/// crate (Phase A, multi_chord) already exercise this path on
/// cube fixtures.
pub(crate) fn find_edge_and_param_at(brep: &BRep, point: Point3) -> Option<(EdgeHandle, f64)> {
    for (eh, edata) in brep.edges.iter() {
        let he_h = edata.half_edge;
        let Some(he) = brep.half_edges.get(he_h) else {
            continue;
        };
        // The second vertex is the twin's origin on manifold
        // edges, or the next half-edge's origin for boundary
        // edges (whose twin is the null sentinel).
        let v1_pos = if let Some(twin) = brep.half_edges.get(he.twin) {
            brep.vertices[twin.origin].position
        } else if let Some(next_he) = brep.half_edges.get(he.next) {
            brep.vertices[next_he.origin].position
        } else {
            continue;
        };
        let v0 = brep.vertices[he.origin].position;
        let v1 = v1_pos;
        let dx = v1.x - v0.x;
        let dy = v1.y - v0.y;
        let dz = v1.z - v0.z;
        let len_sq = dx * dx + dy * dy + dz * dz;
        if len_sq < 1e-18 {
            continue;
        }
        let t = ((point.x - v0.x) * dx + (point.y - v0.y) * dy + (point.z - v0.z) * dz) / len_sq;
        if !(0.001..=0.999).contains(&t) {
            continue;
        }
        let px = v0.x + t * dx;
        let py = v0.y + t * dy;
        let pz = v0.z + t * dz;
        let dist_sq = (point.x - px).powi(2) + (point.y - py).powi(2) + (point.z - pz).powi(2);
        if dist_sq < 1e-10 {
            return Some((eh, t));
        }
    }
    None
}

/// Insert a new vertex at `position` in the mutable BRep by finding
/// the descendant sub-edge currently containing it and splitting it.
/// Appends `(old_edge, edge_a, edge_b)` to `edge_splits` so the
/// caller can maintain a fillet-edge → descendants map.
///
/// If a vertex already exists at `position` (another Phase A
/// iteration put one there — e.g. a shared cube corner anchor), we
/// reuse it without splitting any edge.
fn insert_anchor_at_position(
    brep: &mut BRep,
    position: Point3,
    edge_splits: &mut Vec<(EdgeHandle, EdgeHandle, EdgeHandle)>,
) -> Result<VertexHandle, String> {
    for (vh, v) in brep.vertices.iter() {
        let d_sq = (v.position.x - position.x).powi(2)
            + (v.position.y - position.y).powi(2)
            + (v.position.z - position.z).powi(2);
        if d_sq < 1e-12 {
            return Ok(vh);
        }
    }
    if let Some((eh, t)) = find_edge_and_param_at(brep, position) {
        let split = brep
            .split_edge(eh, t)
            .map_err(|e| format!("split_edge for anchor failed: {e}"))?;
        edge_splits.push((eh, split.edge_a, split.edge_b));
        return Ok(split.vertex);
    }
    // Curve-aware fallback: the chord test above misses anchors on
    // CURVED adjacent edges twice over — perp-walk math offsets the
    // anchor along the tangent (off-curve by the chord sagitta), and
    // the containment test itself compares against the vertex chord
    // (measured: a boss-wall chain's anchors sat 5e-3 off the r=1
    // wall arc and Phase A failed "no edge contains anchor"). Find
    // the closest interior point on any curved edge within a
    // sagitta-scaled gate and split THERE — the realized vertex
    // lands exactly on the curve, which is more faithful than the
    // requested tangent-offset position.
    let mut best: Option<(EdgeHandle, Point3, f64)> = None;
    for (eh, edge) in brep.edges.iter() {
        let Some(curve) = edge.curve.as_ref() else {
            continue;
        };
        use popkern_core::geometry::curve::CurveEval;
        let (t0, t1) = (curve.t_start, curve.t_end);
        let mut bt = t0;
        let mut bd = f64::INFINITY;
        for k in 0..=32 {
            let t = t0 + (t1 - t0) * k as f64 / 32.0;
            let d = curve.geometry.point_at(t).distance_to(&position);
            if d < bd {
                bd = d;
                bt = t;
            }
        }
        let (mut lo, mut hi) = (
            (bt - (t1 - t0) / 32.0).max(t0),
            (bt + (t1 - t0) / 32.0).min(t1),
        );
        for _ in 0..30 {
            let m1 = lo + (hi - lo) / 3.0;
            let m2 = hi - (hi - lo) / 3.0;
            if curve.geometry.point_at(m1).distance_to(&position)
                < curve.geometry.point_at(m2).distance_to(&position)
            {
                hi = m2;
            } else {
                lo = m1;
            }
        }
        let t = 0.5 * (lo + hi);
        // Interior only — endpoint projections belong to the vertex
        // reuse path above.
        let frac = (t - t0) / (t1 - t0).max(1e-300);
        if !(0.001..=0.999).contains(&frac) {
            continue;
        }
        let proj = curve.geometry.point_at(t);
        let d = proj.distance_to(&position);
        if best.as_ref().is_none_or(|(_, _, bd2)| d < *bd2) {
            best = Some((eh, proj, d));
        }
    }
    let gate = brep.tolerance.linear().max(1.0e-9) * 1.0e4;
    if let Some((eh, proj, d)) = best
        && d < gate
    {
        let tol = brep.tolerance;
        let split = brep
            .split_edge_at_point(eh, proj, &tol)
            .map_err(|e| format!("curve-aware anchor split failed: {e}"))?;
        edge_splits.push((eh, split.edge_a, split.edge_b));
        return Ok(split.vertex);
    }
    Err(format!("no edge contains anchor position {position:?}"))
}

/// Information returned by perpendicular edge splitting.
#[derive(Debug)]
pub struct SplitInfo {
    /// Corrected 3D trim curve for face 1 (kept for compatibility).
    #[allow(dead_code)]
    pub trim_3d_1: NurbsCurve,
    /// Corrected 3D trim curve for face 2 (kept for compatibility).
    #[allow(dead_code)]
    pub trim_3d_2: NurbsCurve,
    /// Anchor vertices on face 1 (start, end matching fillet edge direction).
    pub face_1_anchors: (VertexHandle, VertexHandle),
    /// Anchor vertices on face 2 (start, end matching fillet edge direction).
    pub face_2_anchors: (VertexHandle, VertexHandle),
    /// Edge splits performed: `(old_edge, edge_a, edge_b)`.
    #[allow(dead_code)]
    pub edge_splits: Vec<(EdgeHandle, EdgeHandle, EdgeHandle)>,
}

/// Build a simple linear NURBS curve between two points.
pub(crate) fn line_nurbs(a: &Point3, b: &Point3) -> NurbsCurve {
    NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![*a, *b],
        weights: vec![1.0, 1.0],
    }
}

/// Build a pcurve (UV-space line) between two points projected onto a surface.
fn compute_pcurve_for_trim(surface: &SurfaceGeometry, pt_a: &Point3, pt_b: &Point3) -> NurbsCurve {
    let (u_a, v_a) = project_to_uv(surface, pt_a);
    let (u_b, v_b) = project_to_uv(surface, pt_b);
    NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(u_a, v_a, 0.0), Point3::new(u_b, v_b, 0.0)],
        weights: vec![1.0, 1.0],
    }
}

/// Project a 3D point onto a surface's UV parameter space.
pub(crate) fn project_to_uv(surface: &SurfaceGeometry, point: &Point3) -> (f64, f64) {
    use popkern_core::geometry::uv_projection::project_point_to_uv;
    match project_point_to_uv(surface, point) {
        Some(uv) => uv,
        None => {
            // Fallback: grid search.
            use popkern_core::geometry::surface::SurfaceEval;
            let ((u_min, u_max), (v_min, v_max)) = surface.parameter_range();
            let u_lo = if u_min.is_finite() { u_min } else { -10.0 };
            let u_hi = if u_max.is_finite() { u_max } else { 10.0 };
            let v_lo = if v_min.is_finite() { v_min } else { -10.0 };
            let v_hi = if v_max.is_finite() { v_max } else { 10.0 };
            let n = 20;
            let mut best = ((u_lo + u_hi) * 0.5, (v_lo + v_hi) * 0.5);
            let mut best_d = f64::INFINITY;
            for iu in 0..n {
                let u = u_lo + (u_hi - u_lo) * (iu as f64 / (n - 1) as f64);
                for iv in 0..n {
                    let v = v_lo + (v_hi - v_lo) * (iv as f64 / (n - 1) as f64);
                    let d = surface.point_at(u, v).distance_to(point);
                    if d < best_d {
                        best_d = d;
                        best = (u, v);
                    }
                }
            }
            best
        }
    }
}

/// Compute a pcurve for a trim curve on the given face surface.
/// Endpoint-projected straight pcurve (the trim-aware variant is
/// [`pcurve_for_face_along`]; this one remains for exact 2-point
/// trims and its unit tests).
#[cfg_attr(not(test), allow(dead_code))]
pub fn pcurve_for_face(
    brep: &BRep,
    face: FaceHandle,
    trim_start: &Point3,
    trim_end: &Point3,
) -> Result<NurbsCurve, String> {
    let face_data = brep
        .faces
        .get(face)
        .ok_or_else(|| format!("face {face:?} not found"))?;
    let surface = face_data
        .surface
        .as_ref()
        .ok_or_else(|| format!("face {face:?} has no surface"))?;
    Ok(compute_pcurve_for_trim(
        &surface.geometry,
        trim_start,
        trim_end,
    ))
}

/// Extract the fillet radius from the edge geometry analysis.
///
/// Uses the distance between a trim curve endpoint and the nearest fillet
/// edge vertex as a fallback radius estimate.
#[allow(dead_code)]
pub fn estimate_radius_from_geometry(
    brep: &BRep,
    fillet_edge: EdgeHandle,
    trim_3d: &NurbsCurve,
) -> f64 {
    let edge_data = match brep.edges.get(fillet_edge) {
        Some(e) => e,
        None => return 0.0,
    };
    let he = match brep.half_edges.get(edge_data.half_edge) {
        Some(h) => h,
        None => return 0.0,
    };
    let v_start_pos = brep.vertices[he.origin].position;

    let (t_start, _) = trim_3d.parameter_range();
    let trim_start = trim_3d.point_at(t_start);

    v_start_pos.distance_to(&trim_start)
}

#[cfg(test)]
#[path = "helpers_tests.rs"]
mod tests;

/// Pcurve for a face-splitting trim that follows a CURVED 3-D trim:
/// samples the curve and projects each sample into the face's chart
/// (degree-1 through the projections, uniform knots — matching the
/// polyline sub-trim authored by `multi_chord`). Straight trims (2
/// control points) delegate to the endpoint projection so analytic
/// splits keep their exact linear pcurves.
pub fn pcurve_for_face_along(
    brep: &BRep,
    face: FaceHandle,
    trim_3d: &NurbsCurve,
) -> Result<NurbsCurve, String> {
    use popkern_core::geometry::CurveEval;
    let face_data = brep
        .faces
        .get(face)
        .ok_or_else(|| format!("face {face:?} not found"))?;
    let surface = face_data
        .surface
        .as_ref()
        .ok_or_else(|| format!("face {face:?} has no surface"))?;
    let (t0, t1) = trim_3d.parameter_range();
    if trim_3d.control_points.len() <= 2 {
        let a = trim_3d.point_at(t0);
        let b = trim_3d.point_at(t1);
        return Ok(compute_pcurve_for_trim(&surface.geometry, &a, &b));
    }
    let n = trim_3d.control_points.len().max(3);
    let mut uv_points = Vec::with_capacity(n);
    for k in 0..n {
        let t = t0 + (t1 - t0) * k as f64 / (n - 1) as f64;
        let p = trim_3d.point_at(t);
        let (u, v) =
            popkern_core::geometry::uv_projection::project_point_to_uv(&surface.geometry, &p)
                .ok_or_else(|| format!("trim sample does not project onto face {face:?}"))?;
        uv_points.push(Point3::new(u, v, 0.0));
    }
    let degree = 1;
    let n_pts = uv_points.len();
    let mut knots = Vec::with_capacity(n_pts + 2);
    knots.push(0.0);
    for k in 0..n_pts {
        knots.push(k as f64 / (n_pts - 1) as f64);
    }
    knots.push(1.0);
    Ok(NurbsCurve {
        degree,
        knots,
        control_points: uv_points,
        weights: vec![1.0; n_pts],
    })
}
