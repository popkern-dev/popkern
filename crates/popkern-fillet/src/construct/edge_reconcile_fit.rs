//! Curved-fit reconciliation and vertex-snap repair helpers.
//!
//! Split from `edge_reconcile.rs` to honour the structural line
//! ratchet: `try_curved_fit` (the NURBS-involving intersection
//! approximation) and `snap_vertex_and_adjacent_chords` (the
//! vertex-move repair used by strip end arcs and corner arc
//! boundaries) are self-contained tails of that pass.

use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::curve::{CurveGeometry, CurveOnEdge};
use popkern_core::topology::{BRep, EdgeHandle};

/// Approximate a NURBS-involving surface intersection by sampling
/// the chord between `p1` and `p2`, alternately projecting each
/// sample onto both surfaces until convergence, and fitting a
/// degree-3 NURBS through the refined points. Avoids the marching
/// solver in `intersect_surfaces` (which loops on small fillet
/// face inputs).
pub(super) fn try_curved_fit(
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
    // Deviation gate: only refit STALE curves. An unconditional
    // NURBS x NURBS refit rewrote healthy authored cap-seam curves
    // (measured: pyramid4_zoo went from 1 pinned seam to 4 seams at
    // 1.35e-2 when the strip x strip combo landed).
    if let Some(cur) = brep.edges.get(edge).and_then(|e| e.curve.as_ref()) {
        use popkern_core::geometry::curve::CurveEval;
        use popkern_core::geometry::surface::SurfaceEval;
        let mut worst = 0.0_f64;
        for i in 0..=6 {
            let t = cur.t_start + (cur.t_end - cur.t_start) * (i as f64) / 6.0;
            let pt = cur.geometry.point_at(t);
            for sof in [sof_1, sof_2] {
                if let Some((u, v)) =
                    popkern_core::geometry::uv_projection::project_point_to_uv(&sof.geometry, &pt)
                {
                    worst = worst.max(pt.distance_to(&sof.geometry.point_at(u, v)));
                }
            }
        }
        let band = brep
            .edges
            .get(edge)
            .and_then(|e| e.tolerance)
            .unwrap_or_else(|| brep.tolerance.linear())
            .max(1e-9)
            * 10.0;
        if worst <= band.max(3.3e-5) {
            // Below the refit floor the curve is already as good as
            // the surfaces' mutual agreement — but when it still
            // exceeds the edge's recorded band, record the measured
            // deviation as the edge's honest tolerance so check 9
            // gauges against it (measured: the boss corner chord at
            // 1.16e-5 — the sampled-skin arc-density limit — sat
            // above its 1e-6 record with nothing left to refit).
            if worst > band / 10.0
                && worst * 1.5 < 5.0e-5
                && worst > brep.tolerance.linear()
                && let Some(edge_mut) = brep.edges.get_mut(edge)
            {
                let recorded = edge_mut.tolerance.unwrap_or(0.0);
                if recorded < worst * 1.5 {
                    edge_mut.tolerance = Some(worst * 1.5);
                    if std::env::var_os("PK_RECONCILE_TRACE").is_some() {
                        eprintln!("[fit-tol] {edge:?} recorded floor-skip {:.3e}", worst * 1.5);
                    }
                }
            }
            return false;
        }
    }

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
        // Alternating projection converges LINEARLY near a shallow
        // (tangential) intersection — a step-size-only test with 8
        // iterations rejected exactly the strip-boundary edges the
        // pass exists for (PK-141: 2.1e-2 seam remnants on the
        // region-merged strips). Iterate deeper and accept on the
        // RESIDUAL to both surfaces, same as the swept-corner march.
        let max_iters = 200;
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
            // Tangential support-pair rescue: a fillet strip is
            // TANGENT to its own analytic support along the rail,
            // and the alternating projection converges only
            // linearly there (measured stall: d2 = 0 on the strip,
            // d1 stuck 1.09e-4 off the plane). The honest
            // intersection is the strip RAIL — clamp the STRIP's v
            // to each domain end and accept only a witness inside
            // the fit gate (same-face seam pairs and NURBS x NURBS
            // joins must NOT be dragged to rails: unconstrained,
            // this rescue shipped a 2.5e-3 seam defect).
            let analytic_nurbs_pair = matches!(
                (&sof_1.geometry, &sof_2.geometry),
                (
                    popkern_core::geometry::surface::SurfaceGeometry::Plane { .. }
                        | popkern_core::geometry::surface::SurfaceGeometry::Cylinder { .. },
                    popkern_core::geometry::surface::SurfaceGeometry::Nurbs(_)
                ) | (
                    popkern_core::geometry::surface::SurfaceGeometry::Nurbs(_),
                    popkern_core::geometry::surface::SurfaceGeometry::Plane { .. }
                        | popkern_core::geometry::surface::SurfaceGeometry::Cylinder { .. },
                )
            ) && face_1 != face_2;
            if analytic_nurbs_pair {
                let gate = (chord_len * 1e-4).max(3.3e-5);
                let strip = if matches!(
                    &sof_1.geometry,
                    popkern_core::geometry::surface::SurfaceGeometry::Nurbs(_)
                ) {
                    sof_1
                } else {
                    sof_2
                };
                if let Some((u, _)) =
                    popkern_core::geometry::uv_projection::project_point_to_uv(&strip.geometry, &p)
                {
                    let ((_, _), (v0, v1)) = strip.geometry.parameter_range();
                    for v_end in [v0, v1] {
                        if !v_end.is_finite() {
                            continue;
                        }
                        let q = strip.geometry.point_at(u, v_end);
                        let score = [sof_1, sof_2]
                            .iter()
                            .filter_map(|so| {
                                popkern_core::geometry::uv_projection::project_point_to_uv(
                                    &so.geometry,
                                    &q,
                                )
                                .map(|(uu, vv)| q.distance_to(&so.geometry.point_at(uu, vv)))
                            })
                            .fold(0.0_f64, f64::max);
                        if score < gate {
                            p = q;
                            break;
                        }
                    }
                }
            }
            // Residual acceptance: on both surfaces within the fit
            // scale even though the tangential slide has not stalled.
            let d1 = sof_1
                .geometry
                .closest_point(&p)
                .map(|(u, v)| p.distance_to(&sof_1.geometry.point_at(u, v)))
                .unwrap_or(f64::INFINITY);
            // Accept into the RECORDABLE band as well: a pair whose
            // surfaces mutually disagree by a sub-precision gap
            // (independently fitted strip halves at a seam) can
            // never satisfy the relative gate, and declining left a
            // snapped stale fit at 1.43e-5 with no recorded
            // tolerance. The achieved residual is recorded on the
            // edge below.
            if d1 > (chord_len * 1e-4).max(3.3e-5) {
                if std::env::var_os("PK_RECONCILE_TRACE").is_some() {
                    let d2 = sof_2
                        .geometry
                        .closest_point(&p)
                        .map(|(u, v)| p.distance_to(&sof_2.geometry.point_at(u, v)))
                        .unwrap_or(-1.0);
                    eprintln!(
                        "[fit-decline] {edge:?} sample residual d1={d1:.3e} d2={d2:.3e} at                          t={t:.3} p=({:.4},{:.4},{:.4})",
                        p.x, p.y, p.z
                    );
                }
                return false;
            }
        }
        // Reject samples that drift far from the chord — likely
        // landed on a different intersection branch.
        let chord_pt = Point3::new(
            p1.x + t * (p2.x - p1.x),
            p1.y + t * (p2.y - p1.y),
            p1.z + t * (p2.z - p1.z),
        );
        if p.distance_to(&chord_pt) > chord_len * 0.5 {
            if std::env::var_os("PK_RECONCILE_TRACE").is_some() {
                eprintln!(
                    "[fit-decline] {edge:?} branch drift {:.3e} at t={t:.3}",
                    p.distance_to(&chord_pt)
                );
            }
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

    // Achieved-residual tolerance: when the two surfaces themselves
    // disagree along the edge (independently fitted strip halves
    // meet with a mutual G0 gap — measured 1.43e-5 at the boss
    // cross-seam boundary, sampling-density-independent to 12
    // digits), no curve can satisfy both closer than that gap.
    // Record the achieved residual as the edge's own tolerance so
    // check 9 measures against the honest local band.
    let achieved = {
        use popkern_core::geometry::curve::CurveEval;
        let probe = CurveOnEdge {
            geometry: CurveGeometry::Nurbs(fitted.clone()),
            t_start: et_start,
            t_end: et_end,
        };
        let mut worst = 0.0_f64;
        for sof in [sof_1, sof_2] {
            use popkern_core::geometry::surface::SurfaceEval;
            // Continuation-seeded projection: a cold projection can
            // clamp to the wrong chart region on extended strip
            // charts and misreport (measured 2.105e-2 for a fit
            // check 9 gauges at 1.16e-5); each sample seeds from
            // the previous sample's uv.
            let mut prev_uv: Option<(f64, f64)> = None;
            for i in 0..=8 {
                let t = et_start + (et_end - et_start) * (i as f64) / 8.0;
                let pt = probe.geometry.point_at(t);
                let uv = match prev_uv {
                    Some((pu, pv)) => {
                        popkern_core::geometry::uv_projection::project_point_to_uv_near(
                            &sof.geometry,
                            &pt,
                            pu,
                            pv,
                        )
                        .or_else(|| {
                            popkern_core::geometry::uv_projection::project_point_to_uv(
                                &sof.geometry,
                                &pt,
                            )
                        })
                    }
                    None => popkern_core::geometry::uv_projection::project_point_to_uv(
                        &sof.geometry,
                        &pt,
                    ),
                };
                if let Some((u, v)) = uv {
                    prev_uv = Some((u, v));
                    worst = worst.max(pt.distance_to(&sof.geometry.point_at(u, v)));
                }
            }
        }
        worst
    };
    if let Some(edge_mut) = brep.edges.get_mut(edge) {
        edge_mut.curve = Some(CurveOnEdge {
            geometry: CurveGeometry::Nurbs(fitted),
            t_start: et_start,
            t_end: et_end,
        });
        // Record only SUB-PRECISION gaps (the kernel's precision
        // validator caps edge tolerances at 1e-4; a fit that leaves
        // more than that is a real defect to expose, not a band to
        // widen — measured: broad recording pushed 4.7e-2 tolerances
        // onto corner chords and tripped 10 Precision errors).
        let model_band = brep.tolerance.linear().max(1e-9);
        if std::env::var_os("PK_RECONCILE_TRACE").is_some() {
            eprintln!("[fit-tol] {edge:?} achieved {achieved:.3e}");
        }
        // Record sub-precision residuals with the safety margin
        // CAPPED at the recordable ceiling rather than rejected: an
        // achieved of 3.35e-5 measured x1.5 = 5.02e-5 and lost its
        // record by a hair while check 9 gauged the same edge at
        // 1.16e-5 (the kernel Precision validator caps edge
        // tolerances at 1e-4; 9e-5 stays under it).
        if achieved > model_band && achieved < 5.0e-5 {
            let recorded = edge_mut.tolerance.unwrap_or(0.0);
            edge_mut.tolerance = Some(recorded.max((achieved * 1.5).min(9.0e-5)));
        }
        return true;
    }
    false
}

#[allow(dead_code)]
fn pt_zero() -> Point3 {
    Point3::new(0.0, 0.0, 0.0)
}

/// Move `vertex` to `target` and rebuild adjacent edges' endpoint
/// geometry: chord Lines are rebuilt between their (moved) vertex
/// positions; NURBS/polyline curves get their nearest endpoint
/// control point snapped when it sits within the same class of
/// distance. `except` (the edge being re-trimmed) is left alone.
pub(crate) fn snap_vertex_and_adjacent_chords(
    brep: &mut BRep,
    vertex: popkern_core::topology::VertexHandle,
    target: popkern_core::geometry::Point3,
    except: EdgeHandle,
) {
    use popkern_core::geometry::curve::CurveEval;
    let Some(v) = brep.vertices.get_mut(vertex) else {
        return;
    };
    let old_pos = v.position;
    if old_pos.distance_to(&target) < 1e-12 {
        return;
    }
    v.position = target;

    // Collect edges incident to the vertex.
    let incident: Vec<EdgeHandle> = brep
        .edges
        .iter()
        .filter(|(eh, e)| {
            *eh != except && {
                let he = e.half_edge;
                brep.half_edges.get(he).is_some_and(|h| {
                    h.origin == vertex
                        || brep
                            .half_edges
                            .get(h.twin)
                            .is_some_and(|t| t.origin == vertex)
                })
            }
        })
        .map(|(eh, _)| eh)
        .collect();
    for eh in incident {
        let Some(edge) = brep.edges.get(eh) else {
            continue;
        };
        let Some(curve) = edge.curve.clone() else {
            continue;
        };
        let he = edge.half_edge;
        let (va, vb) = {
            let h = &brep.half_edges[he];
            (h.origin, brep.half_edges[h.twin].origin)
        };
        let pa = brep.vertices[va].position;
        let pb = brep.vertices[vb].position;
        match &curve.geometry {
            popkern_core::geometry::CurveGeometry::Line { .. } => {
                // Rebuild the line between the current vertex
                // positions.
                let dir = pb - pa;
                let len = dir.length();
                if len > 1e-12
                    && let Some(e) = brep.edges.get_mut(eh)
                {
                    e.curve = Some(popkern_core::geometry::CurveOnEdge {
                        geometry: popkern_core::geometry::CurveGeometry::Line {
                            origin: pa,
                            direction: dir * (1.0 / len),
                        },
                        t_start: 0.0,
                        t_end: len,
                    });
                }
            }
            popkern_core::geometry::CurveGeometry::Nurbs(n) => {
                // Snap the endpoint control point nearest the moved
                // vertex.
                let ceval: &dyn CurveEval = &curve.geometry;
                let start = ceval.point_at(curve.t_start);
                let end = ceval.point_at(curve.t_end);
                let mut n2 = n.clone();
                let snap_first = start.distance_to(&old_pos) <= end.distance_to(&old_pos);
                let moved = if snap_first {
                    n2.control_points.first_mut()
                } else {
                    n2.control_points.last_mut()
                };
                // The CP belongs to this endpoint when it sits much
                // closer to the old vertex than the move distance
                // (exact-match 1e-6 silently left junction-adjacent
                // curve starts 2e-2 behind a 5e-2 settle move —
                // measured on the boss smooth-junction settles).
                let move_len = old_pos.distance_to(&target);
                if let Some(cp) = moved
                    && cp.distance_to(&old_pos) < (move_len * 0.5).max(1e-6)
                {
                    *cp = target;
                    if let Some(e) = brep.edges.get_mut(eh) {
                        e.curve = Some(popkern_core::geometry::CurveOnEdge {
                            geometry: popkern_core::geometry::CurveGeometry::Nurbs(n2),
                            t_start: curve.t_start,
                            t_end: curve.t_end,
                        });
                    }
                }
            }
            _ => {
                // Analytic curves (arcs etc.) keep their geometry;
                // a moved vertex within the snap class stays within
                // downstream association gates.
            }
        }
    }
}
