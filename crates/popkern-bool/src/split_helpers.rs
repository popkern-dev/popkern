//! Helper functions for face splitting along intersection curves.
//!
//! Extracted from `split.rs` to keep file sizes under the 500-line limit.

use popkern_core::{
    geometry::curve::CurveEval,
    geometry::curve_on_surface::CurveOnSurface,
    geometry::nurbs_curve::NurbsCurve,
    geometry::surface::SurfaceGeometry,
    geometry::{Point3, Vec3},
    tolerance::Tolerance,
    topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle, VertexHandle},
};

use crate::clip::{UvCrossing, clip_curve_to_single_face, face_uv_boundary};
use crate::error::{BoolDiagnostics, BoolError, Stage};
use crate::intersect::IntersectionRecord;
use crate::split_inner_loop::point_inside_face_uv;

use super::FragmentMap;

/// Walk the outer loop of a face, collecting `(half_edge, origin_vertex)` pairs.
pub(super) fn walk_face_loop(
    brep: &BRep,
    face: FaceHandle,
) -> Result<Vec<(HalfEdgeHandle, VertexHandle)>, popkern_core::error::CoreError> {
    let face_data = brep
        .faces
        .get(face)
        .ok_or_else(|| popkern_core::error::CoreError::InvalidHandle("face".into()))?;
    let loop_data = brep
        .loops
        .get(face_data.outer_loop)
        .ok_or_else(|| popkern_core::error::CoreError::InvalidHandle("loop".into()))?;
    let start = loop_data.half_edge;
    let mut pairs = Vec::new();
    let mut cursor = start;
    let limit = brep.half_edges.len() + 1;
    for _ in 0..limit {
        let he = &brep.half_edges[cursor];
        pairs.push((cursor, he.origin));
        cursor = he.next;
        if cursor == start {
            break;
        }
    }
    Ok(pairs)
}

/// Compute the parameter `t` at which `point` projects onto the line segment
/// from `a` to `b`. Returns `None` if the point is not within `tol` of the
/// segment or if `t` is not in `(eps, 1-eps)` (i.e. already at an endpoint).
pub(super) fn point_on_segment(a: &Point3, b: &Point3, point: &Point3, tol: f64) -> Option<f64> {
    let ab = *b - *a;
    let len_sq = ab.length_squared();
    if len_sq < tol * tol {
        return None; // degenerate edge
    }
    let ap = *point - *a;
    let t = ap.dot(&ab) / len_sq;
    // Must be strictly interior — if close to 0 or 1, the endpoint already exists.
    let eps = 1e-8;
    if t <= eps || t >= 1.0 - eps {
        return None;
    }
    let closest = *a + Vec3::new(ab.x * t, ab.y * t, ab.z * t);
    // Adaptive predicate: determinate on-segment decision at the
    // tolerance envelope, avoids the sqrt.
    if popkern_core::adaptive_precision::adaptive_dist_sq_lt(
        (point.x, point.y, point.z),
        (closest.x, closest.y, closest.z),
        tol * tol,
    ) {
        Some(t)
    } else {
        None
    }
}

/// Ensure vertices exist on the face boundary at 3D points (planar fallback).
pub(super) fn ensure_endpoint_vertices(
    brep: &mut BRep,
    face: FaceHandle,
    endpoints: &[Point3],
    tolerance: &Tolerance,
) -> Result<bool, popkern_core::error::CoreError> {
    for endpoint in endpoints {
        // Re-walk the loop each time because split_edge mutates topology.
        let pairs = walk_face_loop(brep, face)?;

        // Check if any existing vertex already matches.
        let already_exists = pairs
            .iter()
            .any(|&(_, vh)| tolerance.points_coincident(&brep.vertices[vh].position, endpoint));
        if already_exists {
            continue;
        }

        // Find which edge the endpoint lies on.
        // Compute t relative to the EDGE's own parameterization
        // (from edge.half_edge.origin to its twin.origin), not the loop
        // half-edge direction which may be reversed.
        let mut best: Option<(EdgeHandle, f64)> = None;
        for &(he_handle, _) in &pairs {
            let he = &brep.half_edges[he_handle];
            let origin_pos = brep.vertices[he.origin].position;
            let next_he = &brep.half_edges[he.next];
            let end_pos = brep.vertices[next_he.origin].position;

            if let Some(_t_loop) =
                point_on_segment(&origin_pos, &end_pos, endpoint, tolerance.linear())
            {
                // Compute t relative to the edge's canonical direction.
                let edge_he = brep.edges[he.edge].half_edge;
                let edge_origin = brep.vertices[brep.half_edges[edge_he].origin].position;
                let edge_twin = brep.half_edges[edge_he].twin;
                let edge_end = brep.vertices[brep.half_edges[edge_twin].origin].position;
                if let Some(t_edge) =
                    point_on_segment(&edge_origin, &edge_end, endpoint, tolerance.linear())
                {
                    best = Some((he.edge, t_edge));
                }
                break;
            }
        }

        if let Some((edge, t)) = best {
            brep.split_edge(edge, t)?;
        } else {
            // Endpoint not on any boundary edge — curve extends beyond this face.
            return Ok(false);
        }
    }
    Ok(true)
}

/// Place vertices on the face boundary using UV crossing edge identity.
pub(super) fn ensure_crossing_vertices(
    brep: &mut BRep,
    face: FaceHandle,
    crossings: &[&UvCrossing],
    curve_3d: &NurbsCurve,
    tolerance: &Tolerance,
) -> Result<bool, popkern_core::error::CoreError> {
    for crossing in crossings {
        let he_handle = crossing.boundary_he;
        let crossing_3d = curve_3d.point_at(crossing.t_curve);
        // Check if vertex already exists.
        let pairs = walk_face_loop(brep, face)?;
        if pairs
            .iter()
            .any(|&(_, vh)| tolerance.points_coincident(&brep.vertices[vh].position, &crossing_3d))
        {
            continue;
        }
        // Check half-edge still exists (previous split may have changed topology).
        if !brep.half_edges.contains_key(he_handle) {
            let pairs2 = walk_face_loop(brep, face)?;
            if pairs2.iter().any(|&(_, vh)| {
                tolerance.points_coincident(&brep.vertices[vh].position, &crossing_3d)
            }) {
                continue;
            }
            return Ok(false);
        }
        let edge_handle = brep.half_edges[he_handle].edge;
        if !brep.edges.contains_key(edge_handle) {
            return Ok(false);
        }
        let canonical_he = brep.edges[edge_handle].half_edge;
        let t_edge = if canonical_he == he_handle {
            crossing.t_boundary
        } else {
            1.0 - crossing.t_boundary
        };
        let eps = 1e-8;
        if t_edge <= eps || t_edge >= 1.0 - eps {
            continue;
        }
        brep.split_edge(edge_handle, t_edge)?;
    }
    Ok(true)
}

/// Returns `true` if the face has planar surface geometry (or no geometry at all).
///
/// Non-planar faces (sphere, cylinder, etc.) are not supported by the
/// Cyrus-Beck polygon clipping used during splitting, so they must be
/// skipped and classified as whole faces instead.
pub(super) fn is_planar_face(brep: &BRep, face: FaceHandle) -> bool {
    match brep.faces.get(face) {
        Some(f) => match &f.surface {
            Some(sof) => matches!(sof.geometry, SurfaceGeometry::Plane { .. }),
            None => true, // No geometry attached — treat as planar (legacy)
        },
        None => true, // Face doesn't exist — caller will skip anyway
    }
}

/// Returns `true` if the face's outer loop half-edges all have pcurves.
pub(super) fn face_has_pcurves(brep: &BRep, face: FaceHandle) -> bool {
    face_uv_boundary(brep, face).is_some()
}

/// Build a [`BoolDiagnostics`] for a split-stage failure.
pub(super) fn split_err(
    face: FaceHandle,
    reason: String,
    records: &[IntersectionRecord],
) -> BoolDiagnostics {
    // §13.1: derive a 3D location from the first intersection curve that
    // touches this face if one is available; otherwise leave None. The
    // first record's curve midpoint is a reasonable proxy for "where on
    // this face the split hit a problem."
    let location = records
        .iter()
        .find(|r| r.face_a == face || r.face_b == face)
        .or_else(|| records.first())
        .map(|r| {
            let (t0, t1) = r.curve.curve_3d.parameter_range();
            let t_mid = (t0 + t1) * 0.5;
            r.curve.curve_3d.point_at(t_mid)
        });
    BoolDiagnostics {
        error: Some(BoolError::SplitFailed {
            face,
            reason: reason.clone(),
        }),
        intersection_curves: records.iter().map(|r| r.curve.clone()).collect(),
        face_classifications: Vec::new(),
        failed_stage: Some(Stage::Split),
        entries: vec![crate::error::DiagnosticEntry {
            severity: crate::error::Severity::Error,
            phase: Stage::Split,
            message: format!("Split failed on face {face:?}: {reason}"),
            entities: vec![crate::error::EntityRef::Face(face)],
            location,
            measured_value: None,
            remediation: Some(crate::error::Remediation {
                action: crate::error::RemediationAction::IncreaseTolerance,
                parameter: Some("linear".into()),
                value: None,
                rationale: "Face splitting failed. Try increasing linear \
                    tolerance or simplifying the input geometry."
                    .into(),
            }),
        }],
        phase_metrics: Vec::new(),
    }
}

/// Try the planar Cyrus-Beck clipping path for a face.
///
/// Returns `true` if the face was split, `false` if it should be skipped.
pub(super) fn try_planar_split(
    brep: &mut BRep,
    face_handle: FaceHandle,
    curve_uv: &NurbsCurve,
    full_curve: &CurveOnSurface,
    tolerance: &Tolerance,
    fragments: &mut FragmentMap,
    records: &[IntersectionRecord],
) -> Result<bool, BoolDiagnostics> {
    let clipped = match clip_curve_to_single_face(brep, face_handle, full_curve) {
        Some(c) => c,
        None => return Ok(false),
    };

    let (t_start, t_end) = clipped.curve_3d.parameter_range();
    let start_pt = clipped.curve_3d.point_at(t_start);
    let end_pt = clipped.curve_3d.point_at(t_end);
    match ensure_endpoint_vertices(brep, face_handle, &[start_pt, end_pt], tolerance) {
        Ok(true) => {}
        Ok(false) => return Ok(false),
        Err(e) => {
            return Err(split_err(
                face_handle,
                format!("pre-split edge failed: {e}"),
                records,
            ));
        }
    }

    // Determine which UV curve to use: for face A it's curve_uv1, for face B
    // it's curve_uv2. Compare the input curve_uv with the record's uv1 vs uv2.
    let use_uv1 = curve_uv.control_points.first().map(|p| p.x)
        == full_curve.curve_uv1.control_points.first().map(|p| p.x);
    let trim_uv = if use_uv1 {
        &clipped.curve_uv1
    } else {
        &clipped.curve_uv2
    };

    match brep.split_face(face_handle, trim_uv, &clipped.curve_3d, tolerance) {
        Ok(result) => {
            fragments.origin.insert(result.face_a, face_handle);
            fragments.origin.insert(result.face_b, face_handle);
            for &edge in &result.new_edges {
                fragments.intersection_edges.insert(edge);
            }
            Ok(true)
        }
        Err(e) => Err(split_err(face_handle, e.to_string(), records)),
    }
}

/// Split a single face along an intersection curve using UV-space or planar
/// clipping depending on the face's geometry.
///
/// Tries UV-space splitting first (for faces with pcurves), then falls back
/// to Cyrus-Beck planar clipping for planar faces without pcurves.
#[allow(clippy::too_many_arguments)]
pub(super) fn split_one_face(
    brep: &mut BRep,
    face_handle: FaceHandle,
    curve_3d: &NurbsCurve,
    curve_uv: &NurbsCurve,
    full_curve: &CurveOnSurface,
    tolerance: &Tolerance,
    records: &[IntersectionRecord],
    fragments: &mut FragmentMap,
) -> Result<(), BoolDiagnostics> {
    // DEGENERATE trims split nothing: a record can clip down to a
    // single point when a tool face meets the target face only at a
    // boundary touch (measured: a flush-sided gusset's embedded top
    // wall grazed the arm's side wall exactly at its width-step
    // vertex; the point-trim was retried until the whole boolean
    // failed with "trim curve end point does not match any vertex").
    // Consume them as no-ops.
    {
        use popkern_core::geometry::curve::CurveEval as _;
        let (t0, t1) = curve_3d.parameter_range();
        let mut len = 0.0;
        let mut prev = curve_3d.point_at(t0);
        for i in 1..=8 {
            let q = curve_3d.point_at(t0 + (t1 - t0) * i as f64 / 8.0);
            len += prev.distance_to(&q);
            prev = q;
        }
        if len < tolerance.linear() * 10.0 {
            if std::env::var_os("PK_SPLIT_TRACE").is_some() {
                eprintln!("  [split1] {face_handle:?} DEGENERATE trim (len {len:.2e}) — skipped");
            }
            return Ok(());
        }
    }
    // T-junction pre-split: a trim curve's endpoint may land MID-EDGE
    // — most commonly on an edge created by an earlier split of the
    // same face (a bounded cutter's bottom and side planes cross a
    // wall in an L; the second curve ends ON the first curve's edge).
    // The face split's entry/exit machinery matches endpoints against
    // VERTICES only, so without a vertex at the junction the whole
    // boolean failed with "trim curve end point does not match any
    // vertex". Split the boundary edge at the endpoint first.
    presplit_boundary_at_endpoints(brep, face_handle, curve_3d, tolerance);

    if std::env::var_os("PK_SPLIT_TRACE").is_some() {
        let (t0, t1) = curve_3d.parameter_range();
        let p0 = curve_3d.point_at(t0);
        let p1 = curve_3d.point_at(t1);
        eprintln!(
            "  [split1] {face_handle:?} rec_closed={} pcurves={} reflex={} trim ({:.4},{:.4},{:.4})->({:.4},{:.4},{:.4})",
            full_curve.is_closed,
            face_has_pcurves(brep, face_handle),
            crate::clip::face_polygon_is_reflex(brep, face_handle),
            p0.x,
            p0.y,
            p0.z,
            p1.x,
            p1.y,
            p1.z
        );
    }
    if face_has_pcurves(brep, face_handle) {
        // UV-space splitting for faces with pcurves (curved or planar).
        if super::helpers_uv::try_uv_split(
            brep,
            face_handle,
            curve_3d,
            curve_uv,
            tolerance,
            fragments,
            records,
        )? {
            return Ok(());
        }
        // UV split didn't produce crossings — try planar fallback if applicable.
    }

    if face_has_pcurves(brep, face_handle) && crate::clip::face_polygon_is_reflex(brep, face_handle)
    {
        // A face WITH pcurves whose UV split found no crossings has a
        // fully-interior trim (island material, handled by the island
        // pre-pass). The planar fallback below clips with Cyrus-Beck,
        // which is CONVEX-only — on a reflex outline it fabricates
        // crossings at outline edges (measured: seat-ring trims on a
        // reflex deck top split two arm edges 4 mm outside the seat,
        // and the polluted topology broke the NEXT boolean with 35
        // open edges at the seat floor).
        return Ok(());
    }
    if is_planar_face(brep, face_handle) {
        // Planar Cyrus-Beck clipping path (fallback).
        try_planar_split(
            brep,
            face_handle,
            curve_uv,
            full_curve,
            tolerance,
            fragments,
            records,
        )?;
    }
    // Non-planar face without pcurves: skip (classified as whole face).

    Ok(())
}

/// Split the face's boundary edges wherever one of the trim curve's
/// endpoints lies mid-edge (no coincident vertex). Failures are
/// non-fatal: the subsequent face split reports its own precondition
/// error if the junction stays vertex-less.
fn presplit_boundary_at_endpoints(
    brep: &mut BRep,
    face_handle: FaceHandle,
    curve_3d: &NurbsCurve,
    tolerance: &Tolerance,
) {
    use popkern_core::geometry::curve_inversion::invert_curve_at_point;
    let (t0, t1) = curve_3d.parameter_range();
    for endpoint in [curve_3d.point_at(t0), curve_3d.point_at(t1)] {
        let Some(face) = brep.faces.get(face_handle) else {
            return;
        };
        let mut loops = vec![face.outer_loop];
        loops.extend(face.inner_loops.iter().copied());
        let mut on_vertex = false;
        let mut split_target: Option<popkern_core::topology::EdgeHandle> = None;
        // PASS 1 — the COMPLETE boundary vertex scan. An existing
        // vertex within tolerance always wins over splitting an edge:
        // when a neighbouring face's trim already minted the junction
        // vertex, a single interleaved walk could reach a splittable
        // edge FIRST (loop-order dependent) and re-split at the same
        // point, minting a ~1e-13 stub edge between two coincident
        // vertices. Such a stub is a point contact, not an edge — and
        // it bridged classification flood-fill across the intersection
        // loop (measured: torus × cylinder union kept all 27
        // fragments; the y<0 junctions were saved only by loop order).
        'vertex_scan: for lh in &loops {
            let Some(loop_data) = brep.loops.get(*lh) else {
                continue;
            };
            let start = loop_data.half_edge;
            let mut he = start;
            loop {
                let he_data = &brep.half_edges[he];
                let vpos = brep.vertices[he_data.origin].position;
                if tolerance.points_coincident(&vpos, &endpoint) {
                    on_vertex = true;
                    break 'vertex_scan;
                }
                he = he_data.next;
                if he == start {
                    break;
                }
            }
        }
        // PASS 2 — only a boundary with NO coincident vertex may have
        // an edge split at the endpoint.
        'search: for lh in loops {
            if on_vertex {
                break;
            }
            let Some(loop_data) = brep.loops.get(lh) else {
                continue;
            };
            let start = loop_data.half_edge;
            let mut he = start;
            loop {
                let he_data = &brep.half_edges[he];
                let edge = he_data.edge;
                if let Some(curve) = brep.edges[edge].curve.as_ref() {
                    let inv = invert_curve_at_point(curve, endpoint, tolerance);
                    if std::env::var_os("PK_SPLIT_TRACE").is_some() {
                        eprintln!(
                            "    [presplit-inv] edge {edge:?} window ({:.3},{:.3}) -> {inv:?}",
                            curve.t_start, curve.t_end
                        );
                    }
                    if let Ok(t) = inv {
                        let lo = curve.t_start.min(curve.t_end);
                        let hi = curve.t_start.max(curve.t_end);
                        let span = (hi - lo).max(1e-12);
                        let margin = span * 1e-6;
                        if t > lo + margin
                            && t < hi - margin
                            && curve.geometry.point_at(t).distance_to(&endpoint)
                                < tolerance.linear() * 10.0
                        {
                            split_target = Some(edge);
                            break 'search;
                        }
                    }
                }
                he = he_data.next;
                if he == start {
                    break;
                }
            }
        }
        if !on_vertex && let Some(edge) = split_target {
            let r = brep.split_edge_at_point(edge, endpoint, tolerance);
            if std::env::var_os("PK_SPLIT_TRACE").is_some() {
                eprintln!(
                    "  [presplit] {face_handle:?} edge {edge:?} at ({:.1},{:.1},{:.1}) -> {r:?}",
                    endpoint.x, endpoint.y, endpoint.z
                );
            }
        } else if std::env::var_os("PK_SPLIT_TRACE").is_some() {
            eprintln!(
                "  [presplit] {face_handle:?} endpoint ({:.1},{:.1},{:.1}) on_vertex={on_vertex} no-target={}",
                endpoint.x,
                endpoint.y,
                endpoint.z,
                split_target.is_none()
            );
            if let Some(face) = brep.faces.get(face_handle)
                && let Some(loop_data) = brep.loops.get(face.outer_loop)
            {
                let start = loop_data.half_edge;
                let mut he = start;
                for _ in 0..12 {
                    let h = &brep.half_edges[he];
                    let v = brep.vertices[h.origin].position;
                    let w = brep.edges[h.edge].curve.as_ref().map(|c| {
                        use popkern_core::geometry::curve::CurveGeometry;
                        let kind = match &c.geometry {
                            CurveGeometry::Line { .. } => "Line".into(),
                            CurveGeometry::Circle { .. } => "Circle".into(),
                            CurveGeometry::Arc {
                                start_angle,
                                end_angle,
                                ..
                            } => format!("Arc[{start_angle:.3},{end_angle:.3}]"),
                            CurveGeometry::Nurbs(n) => format!(
                                "Nurbs[{:.3},{:.3}]",
                                n.knots.first().copied().unwrap_or(0.0),
                                n.knots.last().copied().unwrap_or(0.0)
                            ),
                            _ => "Other".into(),
                        };
                        (kind, c.t_start, c.t_end)
                    });
                    eprintln!(
                        "    [presplit-loop] he edge={:?} origin=({:.1},{:.1},{:.1}) window={w:?} d_end={:.3e}",
                        h.edge,
                        v.x,
                        v.y,
                        v.z,
                        v.distance_to(&endpoint)
                    );
                    he = h.next;
                    if he == start {
                        break;
                    }
                }
            }
        }
    }
}

/// Find which fragment face (from `original_face`) contains the midpoint of
/// a UV curve. Returns `None` if no fragment contains the point.
pub(super) fn find_fragment_for_curve(
    brep: &BRep,
    original_face: FaceHandle,
    curve_uv: &NurbsCurve,
    fragments: &FragmentMap,
    tolerance: &Tolerance,
) -> Option<FaceHandle> {
    let (t0, t1) = curve_uv.parameter_range();
    let mid_uv = curve_uv.point_at((t0 + t1) / 2.0);

    // BTreeSet, not `origin`'s hash order: the FIRST containing fragment
    // wins, and a midpoint near a shared boundary can pass for two (PK-300).
    let candidates: std::collections::BTreeSet<FaceHandle> = fragments
        .origin
        .iter()
        .filter(|&(_, &orig)| orig == original_face)
        .map(|(&frag, _)| frag)
        .collect();
    for candidate in &candidates {
        if !brep.faces.contains_key(*candidate) {
            continue;
        }
        if point_inside_face_uv(brep, *candidate, mid_uv, tolerance) {
            return Some(*candidate);
        }
    }

    // Fallback: if the original face still exists and no fragment matched,
    // try the original face itself.
    if brep.faces.contains_key(original_face) {
        return Some(original_face);
    }

    None
}

#[cfg(test)]
#[path = "split_helpers_tests.rs"]
mod tests;
