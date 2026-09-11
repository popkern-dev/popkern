//! Junction-cluster geometry reconciliation (PK-472).
//!
//! Near a GRAZE junction cluster every crossing solve is
//! ill-conditioned (1/sin(crossing angle) amplification), and the
//! assembled topology can carry correct wiring whose edge GEOMETRY
//! disagrees at validator-gate scale: an edge curve ending a hair off
//! its vertex (the two solves of one junction differ by the graze
//! slop), or a junction-cluster micro-edge whose synthesized curve
//! sits off its adjacent surfaces by the record-fit slop (measured on
//! the m6 nut's mouth chamfer: a 1.5e-6 shell gap and a 7.9e-5
//! surface deviation on a 6e-4 crest×cone stub — both flagged by the
//! transactional validate, and the whole otherwise-correct boolean
//! declined).
//!
//! Topology is never touched. Every move is bounded and verified:
//! endpoint re-anchors are capped at the junction-slop band, and the
//! micro-edge refit only commits when it measurably reduces the
//! deviation. The transactional post-check still runs after this
//! pass — reconciliation that fails to reach the gates declines
//! exactly as before.

use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle};

/// Re-anchor edge-curve endpoints onto their vertices and refit
/// junction-cluster micro-edges onto their adjacent surfaces.
pub(crate) fn reconcile_junction_cluster_geometry(brep: &mut BRep, tolerance: &Tolerance) {
    let lin = tolerance.linear();
    let trace = popkern_core::trace::trace_enabled("boolean");

    // half-edge -> face map for adjacent-surface lookup.
    let mut he_face: std::collections::HashMap<popkern_core::topology::HalfEdgeHandle, FaceHandle> =
        std::collections::HashMap::new();
    for (fh, f) in brep.faces.iter() {
        for lh in std::iter::once(f.outer_loop).chain(f.inner_loops.iter().copied()) {
            let Some(l) = brep.loops.get(lh) else {
                continue;
            };
            let start = l.half_edge;
            let mut h = start;
            let mut n = 0usize;
            loop {
                he_face.insert(h, fh);
                let Some(he) = brep.half_edges.get(h) else {
                    break;
                };
                h = he.next;
                n += 1;
                if h == start || n > brep.half_edges.len() + 1 {
                    break;
                }
            }
        }
    }

    // VERTEX SNAP: a junction vertex solved on a FITTED boundary
    // curve inherits that curve's sag off the analytic carrier
    // (measured, PK-472: the m6 mouth pierce junction sat 7.9e-5 off
    // the crest cylinder — the fitted crest-helix edge's mid-span
    // sag — and the stub edge pinned to it failed the deviation
    // gate). Alternating projection onto every incident surface
    // pulls the vertex to their common intersection; the edge pass
    // below then re-anchors incident curves.
    {
        let mut vertex_faces: std::collections::HashMap<
            popkern_core::topology::VertexHandle,
            Vec<FaceHandle>,
        > = std::collections::HashMap::new();
        for (hh, he) in brep.half_edges.iter() {
            if let Some(f) = he_face.get(&hh) {
                let e = vertex_faces.entry(he.origin).or_default();
                if !e.contains(f) {
                    e.push(*f);
                }
            }
        }
        let snap_band = lin * 200.0;
        for (vh, faces) in &vertex_faces {
            let Some(p0) = brep.vertices.get(*vh).map(|v| v.position) else {
                continue;
            };
            let dist_to = |p: &Point3, fh: FaceHandle, brep: &BRep| -> Option<f64> {
                let sof = brep.faces.get(fh)?.surface.as_ref()?;
                let (u, v) = popkern_core::geometry::uv_projection::project_point_to_uv_multisheet(
                    &sof.geometry,
                    p,
                )?;
                Some(sof.geometry.point_at(u, v).distance_to(p))
            };
            let worst0 = faces
                .iter()
                .filter_map(|&f| dist_to(&p0, f, brep))
                .fold(0.0f64, f64::max);
            if worst0 <= lin * 2.0 || worst0 > snap_band {
                continue;
            }
            let mut p = p0;
            for _ in 0..6 {
                for &fh in faces {
                    let Some(sof) = brep.faces.get(fh).and_then(|f| f.surface.as_ref()) else {
                        continue;
                    };
                    if let Some((u, v)) =
                        popkern_core::geometry::uv_projection::project_point_to_uv_multisheet(
                            &sof.geometry,
                            &p,
                        )
                    {
                        let q = sof.geometry.point_at(u, v);
                        if q.distance_to(&p) <= snap_band {
                            p = q;
                        }
                    }
                }
            }
            let worst1 = faces
                .iter()
                .filter_map(|&f| dist_to(&p, f, brep))
                .fold(0.0f64, f64::max);
            if worst1 < worst0 && p.distance_to(&p0) <= snap_band {
                if trace {
                    eprintln!(
                        "  [junction-reconcile] vertex {vh:?} snapped {:.2e} (worst surface {worst0:.2e}->{worst1:.2e}, {} faces)",
                        p.distance_to(&p0),
                        faces.len()
                    );
                }
                if let Some(v) = brep.vertices.get_mut(*vh) {
                    v.position = p;
                }
            }
        }
    }

    let edge_handles: Vec<EdgeHandle> = brep.edges.iter().map(|(eh, _)| eh).collect();
    for eh in edge_handles {
        let Some((he_h, t_range)) = brep.edges.get(eh).and_then(|e| {
            let c = e.curve.as_ref()?;
            Some((e.half_edge, (c.t_start, c.t_end)))
        }) else {
            continue;
        };
        let (t0, t1) = t_range;
        let Some((v0, v1)) = brep.half_edges.get(he_h).and_then(|he| {
            let p0 = brep.vertices.get(he.origin)?.position;
            let twin = brep.half_edges.get(he.twin)?;
            let p1 = brep.vertices.get(twin.origin)?.position;
            Some((p0, p1))
        }) else {
            continue;
        };
        let Some(c) = brep.edges.get(eh).and_then(|e| e.curve.as_ref()) else {
            continue;
        };
        let p_start = c.geometry.point_at(t0);
        let p_end = c.geometry.point_at(t1);
        // Match curve ends to vertices as an unordered pair, exactly
        // like the watertight check.
        let natural = p_start.distance_to(&v0).max(p_end.distance_to(&v1));
        let swapped = p_start.distance_to(&v1).max(p_end.distance_to(&v0));
        let (target_start, target_end, gap) = if natural <= swapped {
            (v0, v1, natural)
        } else {
            (v1, v0, swapped)
        };

        // Adjacent surfaces (both sides), for deviation measurement.
        let surfaces: Vec<FaceHandle> = {
            let mut out = Vec::new();
            if let Some(he) = brep.half_edges.get(he_h) {
                if let Some(f) = he_face.get(&he_h) {
                    out.push(*f);
                }
                if let Some(f) = he_face.get(&he.twin)
                    && !out.contains(f)
                {
                    out.push(*f);
                }
            }
            out
        };
        let deviation_of = |pts: &[Point3], brep: &BRep| -> f64 {
            let mut worst = 0.0f64;
            for fh in &surfaces {
                let Some(sof) = brep.faces.get(*fh).and_then(|f| f.surface.as_ref()) else {
                    continue;
                };
                let mut prev_uv: Option<(f64, f64)> = None;
                for p in pts {
                    let cold =
                        popkern_core::geometry::uv_projection::project_point_to_uv_multisheet(
                            &sof.geometry,
                            p,
                        );
                    let seeded = prev_uv.and_then(|(pu, pv)| {
                        popkern_core::geometry::uv_projection::project_point_to_uv_near(
                            &sof.geometry,
                            p,
                            pu,
                            pv,
                        )
                    });
                    let best = [cold, seeded]
                        .into_iter()
                        .flatten()
                        .map(|(u, v)| ((u, v), sof.geometry.point_at(u, v).distance_to(p)))
                        .min_by(|a, b| a.1.total_cmp(&b.1));
                    if let Some(((u, v), d)) = best {
                        prev_uv = Some((u, v));
                        worst = worst.max(d);
                    }
                }
            }
            worst
        };

        // ENDPOINT RE-ANCHOR: a gap above the watertight gate but
        // inside the junction-slop band is the two-solves-of-one-
        // junction class.
        let anchor_band = lin * 200.0;
        let needs_anchor = gap > lin * 0.5 && gap <= anchor_band;

        // MICRO-EDGE CARRIER REFIT: short junction-cluster edges
        // whose curve rides off the adjacent surfaces by more than
        // the deviation gate (edge tolerance unset here, so the
        // validator gate is 10×linear).
        let chord = p_start.distance_to(&p_end);
        let n_s = 9usize;
        let samples: Vec<Point3> = (0..n_s)
            .map(|k| {
                c.geometry
                    .point_at(t0 + (t1 - t0) * k as f64 / (n_s - 1) as f64)
            })
            .collect();
        let dev = deviation_of(&samples, brep);
        let needs_refit = chord < lin * 5000.0 && dev > lin * 5.0;

        if !needs_anchor && !needs_refit {
            continue;
        }

        // Build the reconciled sample set: ends pinned to the
        // vertices; for the micro-edge class interior samples are
        // pulled onto the adjacent surfaces by alternating
        // projection.
        let mut pts = samples.clone();
        let n_pts = pts.len();
        pts[0] = target_start;
        pts[n_pts - 1] = target_end;
        if needs_refit {
            // Nearest on-surface point with a wrong-basin rescue: a
            // cold NURBS projection can converge into another turn
            // of a multi-turn chart (the PK-463 class); the
            // knot-span-seeded multi-start recovers the true basin.
            let nearest_on = |surf: &popkern_core::geometry::surface::SurfaceGeometry,
                              p: &Point3|
             -> Option<Point3> {
                let mut best: Option<Point3> = None;
                if let Some((u, v)) =
                    popkern_core::geometry::uv_projection::project_point_to_uv_multisheet(surf, p)
                {
                    best = Some(surf.point_at(u, v));
                }
                if best.as_ref().is_none_or(|q| q.distance_to(p) > anchor_band)
                    && let popkern_core::geometry::surface::SurfaceGeometry::Nurbs(nurbs) = surf
                    && let Some(seeds) = nurbs.projection_seeds()
                    && let Some((u, v)) = popkern_core::geometry::surface::closest_point_from_seeds(
                        nurbs,
                        p,
                        &seeds,
                        (false, false),
                        6,
                    )
                {
                    let q = surf.point_at(u, v);
                    if best
                        .as_ref()
                        .is_none_or(|b| q.distance_to(p) < b.distance_to(p))
                    {
                        best = Some(q);
                    }
                }
                best
            };
            for (pi, p) in pts.iter_mut().enumerate().take(n_pts - 1).skip(1) {
                for round in 0..3 {
                    for fh in &surfaces {
                        let Some(sof) = brep.faces.get(*fh).and_then(|f| f.surface.as_ref()) else {
                            continue;
                        };
                        let q = nearest_on(&sof.geometry, p);
                        if trace && round == 0 {
                            eprintln!(
                                "    [junction-reconcile] edge {eh:?} sample {pi} face {fh:?} pull d {:?}",
                                q.as_ref().map(|q| format!("{:.3e}", q.distance_to(p)))
                            );
                        }
                        if let Some(q) = q {
                            // Bounded pull: never move a sample
                            // further than the anchor band — a
                            // wrong-basin projection must not drag
                            // the curve.
                            if q.distance_to(p) <= anchor_band {
                                *p = q;
                            }
                        }
                    }
                }
            }
        }
        let params: Vec<f64> = (0..n_pts)
            .map(|k| t0 + (t1 - t0) * k as f64 / (n_pts - 1) as f64)
            .collect();
        let Ok(fitted) = popkern_core::geometry::fitting::fit_curve(&pts, &params, 3, n_pts) else {
            continue;
        };
        // Commit only when the reconciliation measurably improves
        // BOTH instruments it targets.
        let new_samples: Vec<Point3> = (0..n_s)
            .map(|k| fitted.point_at(t0 + (t1 - t0) * k as f64 / (n_s - 1) as f64))
            .collect();
        let new_dev = deviation_of(&new_samples, brep);
        let new_gap = fitted
            .point_at(t0)
            .distance_to(&target_start)
            .max(fitted.point_at(t1).distance_to(&target_end));
        let gap_ok = new_gap <= gap.max(lin * 0.5);
        let dev_ok = !needs_refit || new_dev < dev;
        if !(gap_ok && dev_ok) && trace {
            eprintln!(
                "  [junction-reconcile] edge {eh:?} REFUSED gap {gap:.2e}->{new_gap:.2e} dev {dev:.2e}->{new_dev:.2e} chord {chord:.2e} (anchor {needs_anchor} refit {needs_refit})",
            );
        }
        if gap_ok && dev_ok {
            if trace {
                eprintln!(
                    "  [junction-reconcile] edge {eh:?} gap {gap:.2e}->{new_gap:.2e} dev {dev:.2e}->{new_dev:.2e} chord {chord:.2e}",
                );
            }
            if let Some(e) = brep.edges.get_mut(eh)
                && let Some(cm) = e.curve.as_mut()
            {
                cm.geometry = popkern_core::geometry::curve::CurveGeometry::Nurbs(fitted);
                cm.t_start = t0;
                cm.t_end = t1;
            }
        }
    }
}
