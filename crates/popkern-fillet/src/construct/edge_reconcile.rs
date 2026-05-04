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
use popkern_core::geometry::edge_triple::{EdgeAnchors, EdgeTripleOptions, compute_edge_geometry};
use popkern_core::geometry::surface::SurfaceGeometry;
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

    // Stage 2b.57: NURBS combos route directly to `try_curved_fit`
    // (chord-projection alternating-projection). Going through
    // `compute_edge_geometry` works (with 2b.54's timeout +
    // fallback) but adds ~4 s to the popkern-fillet test suite —
    // marching is genuinely slower than chord projection on small
    // NURBS fillet face inputs, even after 2b.56's AABB tightening.
    // The unified-API benefit doesn't justify the overhead for
    // this consumer; the marching catch-all stays available in
    // `compute_edge_geometry` for callers that prefer the
    // EdgeTriple shape over per-input optimisation.
    if nurbs_combo {
        return try_curved_fit(brep, edge, he, twin);
    }

    let v1 = he_data.origin;
    let v2 = twin_data.origin;
    let p1 = brep.vertices[v1].position;
    let p2 = brep.vertices[v2].position;

    // Analytic combos go through the EdgeTriple factory (stage
    // 2b.51). The trig dispatcher resolves cyl/sphere/plane pairs
    // in microseconds.
    let anchors = EdgeAnchors::Open { start: p1, end: p2 };
    let opts = EdgeTripleOptions::default();
    let triple = match compute_edge_geometry(
        &sof_1.geometry,
        &sof_2.geometry,
        anchors,
        None,
        tolerance,
        &opts,
    ) {
        Ok(t) => t,
        Err(_) => return false,
    };

    // Defensive G5 verification: EdgeTriple guarantees endpoint
    // accuracy by construction, but we double-check at the strict
    // BRep tolerance (anything looser opens a validate() OpenLoop).
    let strict_tol = tolerance.linear().max(1e-9);
    let trim_start = triple.curve_3d.point_at(triple.t_start);
    let trim_end = triple.curve_3d.point_at(triple.t_end);
    if trim_start.distance_to(&p1).max(trim_end.distance_to(&p2)) > strict_tol {
        return false;
    }

    // Apply: 3D curve on edge, pcurve_a on the half-edge whose face
    // is face_1, pcurve_b on the twin (face_2's half-edge).
    let geometry = triple
        .curve_3d
        .canonical_line_in_range(
            (triple.t_start, triple.t_end),
            tolerance.linear().max(1.0e-8),
        )
        .unwrap_or_else(|| triple.curve_3d.clone());
    let new_curve = CurveOnEdge {
        geometry,
        t_start: triple.t_start,
        t_end: triple.t_end,
    };
    if let Some(edge_mut) = brep.edges.get_mut(edge) {
        edge_mut.curve = Some(new_curve);
    } else {
        return false;
    }
    if let Some(he_mut) = brep.half_edges.get_mut(he) {
        he_mut.pcurve = Some(triple.pcurve_a.clone());
    }
    if let Some(twin_mut) = brep.half_edges.get_mut(twin) {
        twin_mut.pcurve = Some(triple.pcurve_b.clone());
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

#[allow(dead_code)]
fn pt_zero() -> Point3 {
    Point3::new(0.0, 0.0, 0.0)
}
