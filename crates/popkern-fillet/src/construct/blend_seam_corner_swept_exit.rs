//! PK-160 swept-strip pair-corner repair, mismatched-corner-radius
//! branch (the swept analogue of the analytic stage-2 W-exit march).
//!
//! When the two variable-radius laws disagree at the shared corner,
//! the near anchor sits on ONE strip only and the true seam exits
//! through the shared support plane at a point `W` (on both strips
//! AND the plane). The repair marches the strip×strip intersection
//! from the far anchor toward the support plane (predictor step +
//! alternating-projection corrector), bisects the exit, MOVES the
//! near anchor to `W`, re-authors the Δr mismatch stub as the marched
//! `strip_off × plane` intersection and the rail as the marched
//! `strip_on × plane` intersection — all with projected chart
//! pcurves. Same transactional contract as every corner stage: any
//! unresolved step declines and the caller restores its snapshot.

use popkern_core::geometry::fitting::interpolate_curve;
use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{CurveGeometry, CurveOnEdge, Point3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle, VertexHandle};

use crate::types::BlendSeamCornerDeclineReason as Reason;

use super::blend_seam_corner::PairCornerCandidate;
use super::blend_seam_corner_march::authoring::{chart_pcurve, chord_params, params_reversed};
use super::blend_seam_corner_repair::{commit_edge_authoring, perform_flap_surgery};
use super::blend_seam_corner_swept::{project_seeded, refine_onto_both, surface_distance};

/// One edge incident to the near anchor, classified by its flanks.
struct AnchorEdge {
    edge: EdgeHandle,
    support: FaceHandle,
    far_vertex: VertexHandle,
}

/// March steps from the far anchor toward the support plane.
const MARCH_STEPS: usize = 48;

#[allow(clippy::too_many_arguments)]
pub(super) fn repair_swept_mismatch(
    brep: &mut BRep,
    c: &PairCornerCandidate,
    geom_n: &SurfaceGeometry,
    geom_p: &SurfaceGeometry,
    a_on_n: bool,
    tolerance: &Tolerance,
) -> Result<(), Reason> {
    let trace = std::env::var_os("PK_BLEND_SEAM_TRACE").is_some();
    let pos_tol = (tolerance.linear() * 10.0).max(1e-9);
    let (strip_on, geom_on, strip_off, geom_off) = if a_on_n {
        (c.strip_n, geom_n, c.strip_p, geom_p)
    } else {
        (c.strip_p, geom_p, c.strip_n, geom_n)
    };
    let pos_of = |b: &BRep, v: VertexHandle| -> Result<Point3, Reason> {
        Ok(b.vertices.get(v).ok_or(Reason::EulerRejected)?.position)
    };
    let a_pos = pos_of(brep, c.anchor_near)?;
    let b_pos = pos_of(brep, c.anchor_far)?;
    let corner_pos = pos_of(brep, c.corner).unwrap_or(a_pos);

    // Classify the near anchor's other incident edges: rail
    // (strip_on × support) and Δr mismatch stub (strip_off ×
    // support), sharing one support face.
    let (rail, stub) = classify_anchor_edges(brep, c, strip_on, strip_off)?;
    let plane_geom = brep
        .faces
        .get(stub.support)
        .and_then(|f| f.surface.as_ref())
        .map(|s| s.geometry.clone())
        .ok_or(Reason::StubUnresolved)?;
    let (p_origin, p_normal) = match &plane_geom {
        SurfaceGeometry::Plane { origin, normal, .. } => (*origin, *normal),
        _ => return Err(Reason::StubUnresolved),
    };
    let n_hat = p_normal.normalized().ok_or(Reason::StubUnresolved)?;
    let plane_d = |p: Point3| (p - p_origin).dot(&n_hat);
    let stub_far_pos = pos_of(brep, stub.far_vertex)?;
    let rail_far_pos = pos_of(brep, rail.far_vertex)?;
    if surface_distance(geom_off, &stub_far_pos, None) > pos_tol
        || plane_d(stub_far_pos).abs() > pos_tol
    {
        return Err(Reason::StubUnresolved);
    }

    // March from B toward the plane along the strip×strip
    // intersection; bisect the exit W.
    let start_side = plane_d(b_pos);
    if trace {
        eprintln!(
            "[bsc-exit] a=({:.4},{:.4},{:.4}) b=({:.4},{:.4},{:.4}) plane_d(b)={start_side:.3e} rail {:?} stub {:?}",
            a_pos.x, a_pos.y, a_pos.z, b_pos.x, b_pos.y, b_pos.z, rail.edge, stub.edge
        );
    }
    if start_side.abs() <= pos_tol {
        return Err(Reason::MarchUnresolved);
    }
    let mut pts: Vec<Point3> = vec![b_pos];
    let mut hint_on = popkern_core::geometry::uv_projection::project_point_to_uv(geom_on, &b_pos);
    let mut hint_off = popkern_core::geometry::uv_projection::project_point_to_uv(geom_off, &b_pos);
    // Predictor direction: the intersection TANGENT (cross of the
    // two surface normals) — a chord toward the corner fights the
    // corrector, which pulls straight back to the nearest
    // intersection point (measured: the march pinned at the far
    // anchor for 48 steps). Signed toward decreasing distance to the
    // support plane so the walk heads for the exit.
    use popkern_core::geometry::SurfaceEval as _;
    let tangent_at = |p: Point3,
                      hn: &Option<(f64, f64)>,
                      hp: &Option<(f64, f64)>|
     -> Option<popkern_core::geometry::Vec3> {
        let (ua, va) = (*hn)?;
        let (ub, vb) = (*hp)?;
        let na = geom_on.normal_at(ua, va);
        let nb = geom_off.normal_at(ub, vb);
        let t = na.cross(&nb);
        let _ = p;
        t.normalized()
    };
    let step = (corner_pos - b_pos).length() / (MARCH_STEPS as f64 / 3.0);
    // B is the SINGULAR point of the intersection (both strips
    // tangent to the shared far support there, normals parallel):
    // the first predictor aims at the NEAR anchor, which lies along
    // the seam and off the far support; the tangent predictor takes
    // over once the walk escapes the tangency.
    let mut dir = (a_pos - b_pos)
        .normalized()
        .ok_or(Reason::MarchUnresolved)?;
    let mut w: Option<Point3> = None;
    for _ in 0..MARCH_STEPS {
        let cur = *pts.last().ok_or(Reason::MarchUnresolved)?;
        let predicted = cur + dir * step;
        let q = refine_onto_both(
            geom_on,
            geom_off,
            predicted,
            &mut hint_on,
            &mut hint_off,
            pos_tol,
        )
        .ok_or(Reason::MarchUnresolved)?;
        // The strip may END exactly on the support plane (its rail
        // lies in it): the march then approaches the exit
        // asymptotically with no sign change — accept on proximity.
        if plane_d(q).abs() <= pos_tol {
            w = Some(q);
            break;
        }
        if plane_d(q).signum() != start_side.signum() {
            // Bisect the exit between cur and q.
            let (mut lo, mut hi) = (cur, q);
            for _ in 0..48 {
                let mid = Point3::new(
                    0.5 * (lo.x + hi.x),
                    0.5 * (lo.y + hi.y),
                    0.5 * (lo.z + hi.z),
                );
                let Some(m) =
                    refine_onto_both(geom_on, geom_off, mid, &mut hint_on, &mut hint_off, pos_tol)
                else {
                    break;
                };
                if plane_d(m).signum() == start_side.signum() {
                    lo = m;
                } else {
                    hi = m;
                }
            }
            let exit = if plane_d(hi).abs() < plane_d(lo).abs() {
                hi
            } else {
                lo
            };
            if plane_d(exit).abs() > pos_tol {
                return Err(Reason::MarchUnresolved);
            }
            w = Some(exit);
            break;
        }
        let escaped = (q - b_pos).length() > step * 0.5;
        if escaped && let Some(t) = tangent_at(q, &hint_on, &hint_off) {
            // Keep marching the same way along the curve.
            dir = if t.dot(&dir) >= 0.0 { t } else { -t };
        } else if let Some(d) = (q - cur).normalized()
            && escaped
        {
            dir = d;
        }
        pts.push(q);
    }
    if trace && w.is_none() {
        let last = pts.last().copied();
        eprintln!(
            "[bsc-exit] no exit after {} steps; last {:?}",
            pts.len(),
            last.map(|p| (p.x, p.y, p.z))
        );
    }
    let w = w.ok_or(Reason::MarchUnresolved)?;
    // Sanity: the exit stays in the corner neighbourhood.
    if (w - a_pos).length() > (b_pos - a_pos).length() * 2.0 + pos_tol {
        return Err(Reason::MarchUnresolved);
    }
    pts.push(w);
    if trace {
        eprintln!(
            "[bsc-swept] mismatch W exit ({:.4},{:.4},{:.4}) after {} marched samples",
            w.x,
            w.y,
            w.z,
            pts.len()
        );
    }

    // Author the seam (B → W in march order; stored curve follows
    // the primary half-edge, he_n runs anchor_near(=W) → anchor_far).
    let params = chord_params(&pts);
    let primary_seam = brep
        .edges
        .get(c.seam_edge)
        .ok_or(Reason::EulerRejected)?
        .half_edge;
    let seam_curve = if primary_seam == c.he_n {
        let mut rev = pts.clone();
        rev.reverse();
        fit_curve_on_edge(&rev, &params_reversed(&params))?
    } else {
        fit_curve_on_edge(&pts, &params)?
    };
    let seam_staged = stage_projected_pcurves(
        brep,
        &[(c.strip_n, c.he_n, geom_n), (c.strip_p, c.he_p, geom_p)],
        c.anchor_far,
        &pts,
        &params,
    )?;

    // Stub: marched strip_off × plane from W to the retained
    // endpoint; rail: marched strip_on × plane from W to its far
    // endpoint (handles non-straight rails from nonlinear laws).
    let stub_curve_pts = intersection_samples(geom_off, &plane_geom, w, stub_far_pos, pos_tol)?;
    let rail_curve_pts = intersection_samples(geom_on, &plane_geom, w, rail_far_pos, pos_tol)?;

    perform_flap_surgery(brep, c)?;
    if let Some(v) = brep.vertices.get_mut(c.anchor_near) {
        v.position = w;
    }
    commit_edge_authoring(brep, c.seam_edge, seam_curve, seam_staged);
    commit_marched_side_edge(
        brep,
        stub.edge,
        c.anchor_near,
        &stub_curve_pts,
        geom_off,
        &plane_geom,
        strip_off,
        stub.support,
    )?;
    commit_marched_side_edge(
        brep,
        rail.edge,
        c.anchor_near,
        &rail_curve_pts,
        geom_on,
        &plane_geom,
        strip_on,
        rail.support,
    )?;

    for eh in [c.seam_edge, stub.edge, rail.edge] {
        let dev = brep.measure_edge_deviation(eh);
        if dev > (tolerance.linear() * 100.0).max(1e-7) {
            if trace {
                eprintln!("[bsc-swept] mismatch postcondition: {eh:?} deviation {dev:.3e}");
            }
            return Err(Reason::ValidationFailed);
        }
    }
    brep.validate().map_err(|_| Reason::ValidationFailed)?;
    Ok(())
}

/// Rail/stub classification at the near anchor (surface-agnostic).
fn classify_anchor_edges(
    brep: &BRep,
    c: &PairCornerCandidate,
    strip_on: FaceHandle,
    strip_off: FaceHandle,
) -> Result<(AnchorEdge, AnchorEdge), Reason> {
    let mut rail: Option<AnchorEdge> = None;
    let mut stub: Option<AnchorEdge> = None;
    let mut extra = 0usize;
    for (eh, e) in brep.edges.iter() {
        if eh == c.seam_edge {
            continue;
        }
        let Some(he) = brep.half_edges.get(e.half_edge) else {
            continue;
        };
        let Some(tw) = brep.half_edges.get(he.twin) else {
            continue;
        };
        if he.origin != c.anchor_near && tw.origin != c.anchor_near {
            continue;
        }
        let far_vertex = if he.origin == c.anchor_near {
            tw.origin
        } else {
            he.origin
        };
        let (Some(lf), Some(lt)) = (brep.loops.get(he.loop_), brep.loops.get(tw.loop_)) else {
            continue;
        };
        let flanks = [lf.face, lt.face];
        let entry = |support: FaceHandle| AnchorEdge {
            edge: eh,
            support,
            far_vertex,
        };
        if flanks.contains(&strip_on) && !flanks.contains(&strip_off) {
            let support = if flanks[0] == strip_on {
                flanks[1]
            } else {
                flanks[0]
            };
            rail = Some(entry(support));
        } else if flanks.contains(&strip_off) && !flanks.contains(&strip_on) {
            let support = if flanks[0] == strip_off {
                flanks[1]
            } else {
                flanks[0]
            };
            stub = Some(entry(support));
        } else {
            extra += 1;
        }
    }
    let rail = rail.ok_or(Reason::StubUnresolved)?;
    let stub = stub.ok_or(Reason::StubUnresolved)?;
    if extra != 0 || rail.support != stub.support {
        return Err(Reason::StubUnresolved);
    }
    Ok((rail, stub))
}

fn fit_curve_on_edge(pts: &[Point3], params: &[f64]) -> Result<CurveOnEdge, Reason> {
    let curve = interpolate_curve(pts, params, 3).map_err(|_| Reason::MarchUnresolved)?;
    use popkern_core::geometry::CurveEval;
    let (t0, t1) = curve.parameter_range();
    Ok(CurveOnEdge {
        geometry: CurveGeometry::Nurbs(curve),
        t_start: t0,
        t_end: t1,
    })
}

/// Project marched samples into each listed face chart, oriented
/// along that side's half-edge (`he_n` runs W→B, `he_p` B→W; the
/// samples run B→W).
fn stage_projected_pcurves(
    brep: &BRep,
    sides: &[(FaceHandle, HalfEdgeHandle, &SurfaceGeometry)],
    _anchor_far: VertexHandle,
    pts: &[Point3],
    params: &[f64],
) -> Result<Vec<(FaceHandle, HalfEdgeHandle, Pcurve)>, Reason> {
    let mut staged = Vec::with_capacity(sides.len());
    for &(face, he, geom) in sides {
        // he_n's origin is the near anchor (W end).
        let runs_from_w = {
            let h = brep.half_edges.get(he).ok_or(Reason::EulerRejected)?;
            let _ = h;
            he == sides[0].1
        };
        let ordered: Vec<Point3> = if runs_from_w {
            pts.iter().rev().copied().collect()
        } else {
            pts.to_vec()
        };
        let side_params: Vec<f64> = if runs_from_w {
            params_reversed(params)
        } else {
            params.to_vec()
        };
        let mut uv = Vec::with_capacity(ordered.len());
        let mut hint: Option<(f64, f64)> = None;
        for p in &ordered {
            let (u, v) = project_seeded(geom, p, hint).ok_or(Reason::PcurveClosureFailed)?;
            hint = Some((u, v));
            uv.push((u, v));
        }
        staged.push((face, he, chart_pcurve(&uv, &side_params)?));
    }
    Ok(staged)
}

/// Marched samples of `surface × plane` from `from` to `to`.
fn intersection_samples(
    surface: &SurfaceGeometry,
    plane: &SurfaceGeometry,
    from: Point3,
    to: Point3,
    pos_tol: f64,
) -> Result<Vec<Point3>, Reason> {
    const N: usize = 17;
    let mut out = Vec::with_capacity(N);
    let mut hint_s: Option<(f64, f64)> = None;
    let mut hint_p: Option<(f64, f64)> = None;
    for i in 0..N {
        let f = i as f64 / (N - 1) as f64;
        let mut p = Point3::new(
            from.x + f * (to.x - from.x),
            from.y + f * (to.y - from.y),
            from.z + f * (to.z - from.z),
        );
        if i > 0 && i < N - 1 {
            p = refine_onto_both(surface, plane, p, &mut hint_s, &mut hint_p, pos_tol)
                .ok_or(Reason::StubUnresolved)?;
            if surface_distance(surface, &p, hint_s) > pos_tol
                || surface_distance(plane, &p, hint_p) > pos_tol
            {
                return Err(Reason::StubUnresolved);
            }
        }
        out.push(p);
    }
    Ok(out)
}

/// Re-author one anchor-incident edge (stub or rail) from marched
/// samples: fitted 3-D curve plus projected pcurves on both flanks.
#[allow(clippy::too_many_arguments)]
fn commit_marched_side_edge(
    brep: &mut BRep,
    edge: EdgeHandle,
    anchor: VertexHandle,
    pts_from_w: &[Point3],
    strip_geom: &SurfaceGeometry,
    plane_geom: &SurfaceGeometry,
    strip_face: FaceHandle,
    support_face: FaceHandle,
) -> Result<(), Reason> {
    let e = brep.edges.get(edge).ok_or(Reason::StubUnresolved)?;
    let primary = e.half_edge;
    let he = brep.half_edges.get(primary).ok_or(Reason::StubUnresolved)?;
    let tw = brep.half_edges.get(he.twin).ok_or(Reason::StubUnresolved)?;
    let primary_from_anchor = he.origin == anchor;
    let params = chord_params(pts_from_w);
    let curve = if primary_from_anchor {
        fit_curve_on_edge(pts_from_w, &params)?
    } else {
        let rev: Vec<Point3> = pts_from_w.iter().rev().copied().collect();
        fit_curve_on_edge(&rev, &params_reversed(&params))?
    };
    let mut staged: Vec<(FaceHandle, HalfEdgeHandle, Pcurve)> = Vec::with_capacity(2);
    for (heh, hh) in [(primary, he), (he.twin, tw)] {
        let face = brep
            .loops
            .get(hh.loop_)
            .map(|l| l.face)
            .ok_or(Reason::StubUnresolved)?;
        let geom = if face == strip_face {
            strip_geom
        } else if face == support_face {
            plane_geom
        } else {
            return Err(Reason::StubUnresolved);
        };
        let from_anchor = hh.origin == anchor;
        let ordered: Vec<Point3> = if from_anchor {
            pts_from_w.to_vec()
        } else {
            pts_from_w.iter().rev().copied().collect()
        };
        let side_params = if from_anchor {
            params.clone()
        } else {
            params_reversed(&params)
        };
        let mut uv = Vec::with_capacity(ordered.len());
        let mut hint: Option<(f64, f64)> = None;
        for p in &ordered {
            let (u, v) = project_seeded(geom, p, hint).ok_or(Reason::PcurveClosureFailed)?;
            hint = Some((u, v));
            uv.push((u, v));
        }
        staged.push((face, heh, chart_pcurve(&uv, &side_params)?));
    }
    commit_edge_authoring(brep, edge, curve, staged);
    Ok(())
}
