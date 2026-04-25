//! Post-fillet edge reconciliation.
//!
//! After topology operations and surface assignment, edges between
//! curved faces (cylinder, sphere, NURBS) and adjacent faces can
//! end up with 3-D curves that don't lie on both adjacent surfaces.
//! `multi_chord` and corner-blend `mef` calls write straight-line
//! curves that were correct *at the time* (the adjacent face was
//! still planar) but became inconsistent once the face surface was
//! later replaced by a cylindrical/curved fit.
//!
//! Symptom in OCCT-strict slicers (Bambu, FreeCAD): cylindrical
//! fillet ribbons show up with sharp boundaries because the
//! slicer detects the edge geometry doesn't lie on the cylinder
//! surface and refuses to render the face.
//!
//! This pass walks every edge, computes the analytic surface ×
//! surface intersection of its two adjacent faces, and replaces
//! the edge's 3-D curve with the intersection curve trimmed to
//! the existing vertex endpoints — when the trim succeeds and
//! the new curve hits both endpoints within tolerance.
//!
//! ## Scope
//!
//! Currently handles **cylinder × plane** only. Other combinations
//! (cylinder × cylinder, cylinder × NURBS, sphere × plane, etc.)
//! are out of scope for this first pass; those edges keep their
//! existing curves. Pcurves are also untouched here — they remain
//! whatever was assigned at construction time. A more complete
//! reconciliation pass would update both pcurves and possibly
//! topology (kev internal edges); that's tracked as future work.

use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::curve::{CurveGeometry, CurveOnEdge};
use popkern_core::geometry::intersection::surface_surface::{
    SurfaceSurfaceIntersection, intersect_surfaces,
};
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::uv_projection::tight_bounded_parameter_range;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle};

/// Walk every edge in `brep`. For each edge bordered on one side
/// by a cylindrical face and on the other by a planar face,
/// replace the edge's 3-D curve with the analytic cylinder × plane
/// intersection curve, trimmed to the existing vertex endpoints.
/// Edges where the intersection round-trip doesn't match the
/// vertex positions within tolerance are left alone.
pub(crate) fn reconcile_edge_curves(brep: &mut BRep) {
    let edge_handles: Vec<EdgeHandle> = brep.edges.keys().collect();
    let tolerance = brep.tolerance;
    for edge in edge_handles {
        try_reconcile_cyl_plane_edge(brep, edge, &tolerance);
    }
}

fn try_reconcile_cyl_plane_edge(brep: &mut BRep, edge: EdgeHandle, tolerance: &Tolerance) -> bool {
    let Some(edge_data) = brep.edges.get(edge) else {
        return false;
    };
    let he = edge_data.half_edge;
    let Some(he_data) = brep.half_edges.get(he) else {
        return false;
    };
    let twin = he_data.twin;
    let Some(twin_data) = brep.half_edges.get(twin) else {
        return false;
    };

    let face_1 = brep.loops[he_data.loop_].face;
    let face_2 = brep.loops[twin_data.loop_].face;
    if face_1 == face_2 {
        // Self-loop edge — pcurve repair only, skip 3D curve rewrite.
        return false;
    }

    let Some(sof_1) = brep.faces.get(face_1).and_then(|f| f.surface.as_ref()) else {
        return false;
    };
    let Some(sof_2) = brep.faces.get(face_2).and_then(|f| f.surface.as_ref()) else {
        return false;
    };

    // Whitelist of analytic surface-pair combos that
    // `intersect_surfaces` resolves without falling back to the
    // marching solver. Plane × plane is left alone; the line is
    // already correct. NURBS-involving combos go through a
    // separate alternating-projection fit (`approximate_curved_intersection`)
    // — `intersect_surfaces` would hit the marching path which can
    // loop indefinitely on small fillet face inputs.
    let analytic_combo = matches!(
        (&sof_1.geometry, &sof_2.geometry),
        (
            SurfaceGeometry::Cylinder { .. },
            SurfaceGeometry::Plane { .. }
        ) | (
            SurfaceGeometry::Plane { .. },
            SurfaceGeometry::Cylinder { .. }
        ) | (
            SurfaceGeometry::Sphere { .. },
            SurfaceGeometry::Plane { .. }
        ) | (
            SurfaceGeometry::Plane { .. },
            SurfaceGeometry::Sphere { .. }
        ) | (
            SurfaceGeometry::Sphere { .. },
            SurfaceGeometry::Cylinder { .. }
        ) | (
            SurfaceGeometry::Cylinder { .. },
            SurfaceGeometry::Sphere { .. }
        )
    );
    let nurbs_combo = matches!(
        (&sof_1.geometry, &sof_2.geometry),
        (SurfaceGeometry::Cylinder { .. }, SurfaceGeometry::Nurbs(_))
            | (SurfaceGeometry::Nurbs(_), SurfaceGeometry::Cylinder { .. })
            | (SurfaceGeometry::Sphere { .. }, SurfaceGeometry::Nurbs(_))
            | (SurfaceGeometry::Nurbs(_), SurfaceGeometry::Sphere { .. })
            | (SurfaceGeometry::Plane { .. }, SurfaceGeometry::Nurbs(_))
            | (SurfaceGeometry::Nurbs(_), SurfaceGeometry::Plane { .. })
    );
    if !analytic_combo && !nurbs_combo {
        return false;
    }
    if nurbs_combo {
        // Skip `intersect_surfaces` for NURBS combos and try the
        // bounded alternating-projection fit instead.
        return try_curved_fit(brep, edge, he, twin);
    }

    let v1 = he_data.origin;
    let v2 = twin_data.origin;
    let p1 = brep.vertices[v1].position;
    let p2 = brep.vertices[v2].position;

    let positions_1 = brep.face_outer_loop_positions(face_1);
    let positions_2 = brep.face_outer_loop_positions(face_2);
    if positions_1.is_empty() || positions_2.is_empty() {
        return false;
    }
    let (u1_range, v1_range) = tight_bounded_parameter_range(&positions_1, &sof_1.geometry);
    let (u2_range, v2_range) = tight_bounded_parameter_range(&positions_2, &sof_2.geometry);

    let result = match intersect_surfaces(
        &sof_1.geometry,
        u1_range,
        v1_range,
        &sof_2.geometry,
        u2_range,
        v2_range,
        tolerance,
    ) {
        Ok(r) => r,
        Err(_) => return false,
    };

    // Score each intersection curve by how close it passes to both
    // vertex positions. Pick the best.
    let mut best = None;
    let mut best_score = f64::INFINITY;
    for inter in &result.intersections {
        let SurfaceSurfaceIntersection::Curve(curve_on_surf) = inter else {
            continue;
        };
        let Some(t1) = curve_on_surf.curve_3d.closest_point(&p1) else {
            continue;
        };
        let Some(t2) = curve_on_surf.curve_3d.closest_point(&p2) else {
            continue;
        };
        let d1 = curve_on_surf.curve_3d.point_at(t1).distance_to(&p1);
        let d2 = curve_on_surf.curve_3d.point_at(t2).distance_to(&p2);
        let score = d1.max(d2);
        if score < best_score {
            best_score = score;
            best = Some((
                curve_on_surf.curve_3d.clone(),
                curve_on_surf.curve_uv1.clone(),
                curve_on_surf.curve_uv2.clone(),
                t1,
                t2,
            ));
        }
    }
    let Some((curve_3d, curve_uv1, curve_uv2, t1, t2)) = best else {
        return false;
    };

    // The vertex positions in the BRep can drift up to several
    // percent of cylinder radius from the analytic cylinder
    // surface (a known kernel issue with corner-blend vertex
    // placement). Reject the rewrite if the chord between the
    // existing line endpoints would beat the analytic intersection
    // — that means the intersection is very far from the actual
    // edge geometry and we'd produce worse output.
    let chord_len = p1.distance_to(&p2);
    if best_score > (chord_len * 0.25).max(1e-6) {
        return false;
    }

    // Trim the 3D curve and both pcurves to the parameter range
    // [t1, t2] (assumed t1 <= t2; if reversed, abort — we'd have
    // to flip the curve, which the existing helpers don't do).
    if t1 > t2 {
        return false;
    }
    let (trimmed_3d, trimmed_uv1, trimmed_uv2) = match (
        slice_curve(&curve_3d, t1, t2),
        slice_curve(&curve_uv1, t1, t2),
        slice_curve(&curve_uv2, t1, t2),
    ) {
        (Some(a), Some(b), Some(c)) => (a, b, c),
        _ => return false,
    };

    // Verify the trimmed 3D curve endpoints land at p1 and p2.
    let (et_start, et_end) = trimmed_3d.parameter_range();
    let trim_start = trimmed_3d.point_at(et_start);
    let trim_end = trimmed_3d.point_at(et_end);
    let trim_err = trim_start.distance_to(&p1).max(trim_end.distance_to(&p2));
    if trim_err > (chord_len * 0.25).max(1e-6) {
        return false;
    }

    // Apply: 3D curve on edge, pcurve_1 on the half-edge whose face
    // is face_1, pcurve_2 on the twin (face_2's half-edge).
    let new_curve = CurveOnEdge {
        geometry: CurveGeometry::Nurbs(trimmed_3d),
        t_start: et_start,
        t_end: et_end,
    };
    if let Some(edge_mut) = brep.edges.get_mut(edge) {
        edge_mut.curve = Some(new_curve);
    } else {
        return false;
    }
    if let Some(he_mut) = brep.half_edges.get_mut(he) {
        he_mut.pcurve = Some(trimmed_uv1);
    }
    if let Some(twin_mut) = brep.half_edges.get_mut(twin) {
        twin_mut.pcurve = Some(trimmed_uv2);
    }
    true
}

/// Approximate a NURBS-involving surface intersection by sampling
/// the chord between `p1` and `p2`, alternately projecting each
/// sample onto both surfaces until convergence, and fitting a
/// degree-3 NURBS through the refined points. Avoids the marching
/// solver in `intersect_surfaces` (which loops on small fillet
/// face inputs).
fn try_curved_fit(
    brep: &mut BRep,
    edge: EdgeHandle,
    he: popkern_core::topology::HalfEdgeHandle,
    twin: popkern_core::topology::HalfEdgeHandle,
) -> bool {
    use popkern_core::geometry::fitting::fit_curve;
    use popkern_core::geometry::surface::SurfaceEval;

    let Some(he_data) = brep.half_edges.get(he) else {
        return false;
    };
    let Some(twin_data) = brep.half_edges.get(twin) else {
        return false;
    };
    let face_1 = brep.loops[he_data.loop_].face;
    let face_2 = brep.loops[twin_data.loop_].face;
    let v1 = he_data.origin;
    let v2 = twin_data.origin;
    let p1 = brep.vertices[v1].position;
    let p2 = brep.vertices[v2].position;
    let chord_len = p1.distance_to(&p2);
    if chord_len < 1e-9 {
        return false;
    }

    let surf_1 = brep.faces.get(face_1).and_then(|f| f.surface.as_ref());
    let surf_2 = brep.faces.get(face_2).and_then(|f| f.surface.as_ref());
    let (Some(sof_1), Some(sof_2)) = (surf_1, surf_2) else {
        return false;
    };

    // Sample chord, project alternately. 7 interior samples gives
    // a degree-3 fit room while staying cheap.
    let n_segments = 8;
    let mut points: Vec<Point3> = Vec::with_capacity(n_segments + 1);
    let mut params: Vec<f64> = Vec::with_capacity(n_segments + 1);
    points.push(p1);
    params.push(0.0);
    for i in 1..n_segments {
        let t = i as f64 / n_segments as f64;
        let mut p = Point3::new(
            p1.x + t * (p2.x - p1.x),
            p1.y + t * (p2.y - p1.y),
            p1.z + t * (p2.z - p1.z),
        );
        let max_iters = 8;
        let mut prev = p;
        let mut converged = false;
        for _ in 0..max_iters {
            let Some((u1, v1u)) = sof_1.geometry.closest_point(&p) else {
                return false;
            };
            p = sof_1.geometry.point_at(u1, v1u);
            let Some((u2, v2u)) = sof_2.geometry.closest_point(&p) else {
                return false;
            };
            p = sof_2.geometry.point_at(u2, v2u);
            if p.distance_to(&prev) < chord_len * 1e-5 {
                converged = true;
                break;
            }
            prev = p;
        }
        if !converged {
            return false;
        }
        // Reject samples that drift far from the chord — likely
        // landed on a different intersection branch.
        let chord_pt = Point3::new(
            p1.x + t * (p2.x - p1.x),
            p1.y + t * (p2.y - p1.y),
            p1.z + t * (p2.z - p1.z),
        );
        if p.distance_to(&chord_pt) > chord_len * 0.5 {
            return false;
        }
        points.push(p);
        params.push(t);
    }
    points.push(p2);
    params.push(1.0);

    // Fit a degree-3 NURBS through the projected points.
    let degree = 3;
    let num_cps = (points.len()).clamp(degree + 1, points.len());
    let fitted = match fit_curve(&points, &params, degree, num_cps) {
        Ok(c) => c,
        Err(_) => return false,
    };

    // Verify endpoints land at p1 and p2.
    let (et_start, et_end) = fitted.parameter_range();
    let trim_start = fitted.point_at(et_start);
    let trim_end = fitted.point_at(et_end);
    let endpoint_err = trim_start.distance_to(&p1).max(trim_end.distance_to(&p2));
    if endpoint_err > (chord_len * 0.05).max(1e-6) {
        return false;
    }

    if let Some(edge_mut) = brep.edges.get_mut(edge) {
        edge_mut.curve = Some(CurveOnEdge {
            geometry: CurveGeometry::Nurbs(fitted),
            t_start: et_start,
            t_end: et_end,
        });
        return true;
    }
    false
}

/// Slice `curve` to the parameter range `[lo, hi]` using
/// `split_at` twice. Returns the middle segment, or `None` if the
/// splits fail (e.g. range is degenerate or out of bounds).
fn slice_curve(
    curve: &popkern_core::geometry::nurbs_curve::NurbsCurve,
    lo: f64,
    hi: f64,
) -> Option<popkern_core::geometry::nurbs_curve::NurbsCurve> {
    let (full_min, full_max) = curve.parameter_range();
    if !(full_min..=full_max).contains(&lo) || !(full_min..=full_max).contains(&hi) {
        return None;
    }
    if hi <= lo + 1e-12 {
        return None;
    }
    let (_, right) = curve.split_at(lo).ok()?;
    let (rem_min, _rem_max) = right.parameter_range();
    let cut = rem_min + (hi - lo);
    let (left, _) = right.split_at(cut).ok()?;
    Some(left)
}

#[allow(dead_code)]
fn pt_zero() -> Point3 {
    Point3::new(0.0, 0.0, 0.0)
}
