//! Fragment sample-point construction in the owning face's UV chart.
//!
//! Two samplers live here:
//! - [`loop_uv_runs`] walks a fragment loop and returns each
//!   half-edge's UV sample run in the face's chart (pcurve-first,
//!   3D-projection fallback), period-branch-aligned so consecutive
//!   runs are adjacent — the shared substrate for both the scanline
//!   interior sampler in `cellular::boolean` and the cut-offset
//!   sampler below.
//! - [`fragment_cut_offset_sample`] produces a point a small step OFF
//!   the intersection curve into the fragment's own side. It exists
//!   because a fragment cut out of a face can have its ENTIRE
//!   boundary hugging the other body (measured on the flush swoop
//!   web: cut curve on the plate underside, flush side edges on the
//!   arm walls, both ends embedded) — every boundary midpoint then
//!   misrepresents the interior, and which one the walk enumerated
//!   first decided the classification. Stepping off the cut is the
//!   one sample that is decisive by construction: the other body's
//!   surface passes through the fragment only along the cut, so an
//!   ε-offset on the fragment side is cleanly in or out.

use crate::merge::graph::{EdgeSource, FaceRef, IntersectionGraph};
use popkern_core::geometry::Point3;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};

/// One half-edge's UV sample run: (half-edge index, chart samples).
pub(crate) type UvRun = (usize, Vec<(f64, f64)>);

/// Sample each half-edge of `hes` into the UV chart of `face`,
/// returning one run per half-edge (loop order preserved, branches
/// aligned). Returns `None` when any half-edge cannot be charted.
pub(crate) fn loop_uv_runs(
    graph: &IntersectionGraph,
    face: &FaceRef,
    hes: &[usize],
    surface: &SurfaceGeometry,
) -> Option<Vec<UvRun>> {
    use popkern_core::geometry::curve::CurveEval;
    use popkern_core::geometry::uv_projection::{project_point_to_uv, project_point_to_uv_near};

    let trace = std::env::var_os("PK_SEL_TRACE").is_some();
    let u_period = surface.parameter_periodicity().0;
    let mut runs: Vec<UvRun> = Vec::new();
    let mut last_uv: Option<(f64, f64)> = None;

    for &raw in hes {
        // Pcurves live in their OWNER face's chart. A fragment loop
        // may carry either direction of an edge (the walker keeps the
        // first-walked direction), so read the side that belongs to
        // THIS face — the twin's pcurve is in the other body's chart
        // and produced garbage scanline polygons (measured: a
        // counterbore floor's star-hole fragment sampled its own hole
        // centre and was dropped).
        let h = if graph.half_edges[raw].face.as_ref() == Some(face) {
            raw
        } else {
            let twin = graph.half_edges[raw].twin;
            if graph.half_edges[twin].face.as_ref() == Some(face) {
                twin
            } else {
                raw
            }
        };
        let mut run: Vec<(f64, f64)> = Vec::new();
        let kind;
        if let Some(pc) = graph.half_edges[h].pcurve.as_ref() {
            kind = "PCURVE";
            // A split graph half-edge keeps the FULL original edge's
            // pcurve; its own 3D-curve window is the narrowed span in
            // the same parameter space. Sampling the full pcurve made
            // a 1 mm stub's polygon row span the whole 21 mm wall,
            // and the scanline midpoint landed on a different region
            // entirely (a below-body skirt fragment classified Inside
            // and was kept as a phantom face). Use the window when it
            // nests inside the pcurve's range; fall back to the full
            // range otherwise (foreign parameterization).
            let (mut s0, mut s1) = (pc.t_start, pc.t_end);
            if let Some(c) = graph.half_edges[h].curve.as_ref() {
                // The stored pcurve may cover the FULL original edge
                // while the half-edge's 3D window is a stub of it, in
                // a different parameter space (the 3D curve is
                // re-parameterized at split time). Recover the pcurve
                // sub-window by INVERTING the he's 3D endpoints onto
                // the pcurve through the surface: coarse scan +
                // ternary refinement, then verify. Keep the full
                // range when inversion fails (foreign chart).
                let p3a = c.geometry.point_at(c.t_start);
                let p3b = c.geometry.point_at(c.t_end);
                let scale3 = 1.0 + p3a.distance_to(&p3b);
                let nearest = |p3: &Point3| -> (f64, f64) {
                    let d_at = |t: f64| {
                        let uv = pc.point_at(t);
                        surface.point_at(uv.x, uv.y).distance_to(p3)
                    };
                    let n = 64;
                    let mut best = (pc.t_start, f64::INFINITY);
                    for i in 0..=n {
                        let t = pc.t_start + (pc.t_end - pc.t_start) * i as f64 / n as f64;
                        let d = d_at(t);
                        if d < best.1 {
                            best = (t, d);
                        }
                    }
                    let step = (pc.t_end - pc.t_start).abs() / n as f64;
                    let (mut lo, mut hi) = (best.0 - step, best.0 + step);
                    for _ in 0..48 {
                        let m1 = lo + (hi - lo) / 3.0;
                        let m2 = hi - (hi - lo) / 3.0;
                        if d_at(m1) <= d_at(m2) {
                            hi = m2;
                        } else {
                            lo = m1;
                        }
                    }
                    let t = 0.5 * (lo + hi);
                    (t, d_at(t))
                };
                // Scoped to NON-PERIODIC charts: on a periodic chart
                // every u-branch evaluates to the same 3D points, so
                // inversion cannot distinguish branches and the
                // recovered window can run the wrong way around the
                // circle (measured: a drill-hole rim's window emptied
                // the scanline polygon and reopened a through-hole).
                // Closed spans (rim half-edges) invert both endpoints
                // to the same parameter — skip those too, and verify
                // the midpoint.
                let u_periodic = surface.parameter_periodicity().0.is_some();
                if !u_periodic && p3a.distance_to(&p3b) > 1e-9 * scale3 {
                    let (ta, da) = nearest(&p3a);
                    let (tb, db) = nearest(&p3b);
                    let mid_ok = {
                        let t_mid_pc = 0.5 * (ta + tb);
                        let uv = pc.point_at(t_mid_pc);
                        let p3_mid = c.geometry.point_at(0.5 * (c.t_start + c.t_end));
                        surface.point_at(uv.x, uv.y).distance_to(&p3_mid) < 1e-3 * scale3
                    };
                    if da < 1e-6 * scale3 && db < 1e-6 * scale3 && (ta - tb).abs() > 1e-12 && mid_ok
                    {
                        s0 = ta;
                        s1 = tb;
                    }
                }
            }
            for i in 0..8 {
                let t = s0 + (s1 - s0) * i as f64 / 8.0;
                let p = pc.point_at(t);
                run.push((p.x, p.y));
            }
        } else if let Some(c) = graph.half_edges[h].curve.as_ref() {
            // No pcurve (split edges): project curve samples, with
            // continuation to stay on one branch.
            kind = "CURVE";
            let ceval: &dyn CurveEval = &c.geometry;
            for i in 0..8 {
                let t = c.t_start + (c.t_end - c.t_start) * i as f64 / 8.0;
                let p = ceval.point_at(t);
                let uv = match run.last().or(last_uv.as_ref()) {
                    Some(&(pu, pv)) => project_point_to_uv_near(surface, &p, pu, pv)
                        .or_else(|| project_point_to_uv(surface, &p)),
                    None => project_point_to_uv(surface, &p),
                }?;
                run.push(uv);
            }
        } else {
            // Straight graph edge without geometry: project the
            // endpoints.
            kind = "ENDPOINTS";
            let a = graph.nodes[graph.half_edges[h].origin].position;
            let b = graph.nodes[graph.half_edges[graph.half_edges[h].twin].origin].position;
            for p in [a, b] {
                run.push(project_point_to_uv(surface, &p)?);
            }
        }
        // Stored pcurves may live on ANY period branch of a periodic
        // chart (a cylinder pcurve at u−2π is the same geometry),
        // while projected samples come back on the principal branch.
        // Mixing branches in one polygon opens a phantom interval
        // spanning the branch gap. Shift each run by whole periods so
        // it lands adjacent to the previous one.
        if let (Some(p), Some((lu, _)), Some(&(fu, _))) = (u_period, last_uv, run.first()) {
            let k = ((lu - fu) / p).round();
            if k != 0.0 {
                for uv in &mut run {
                    uv.0 += k * p;
                }
            }
        }
        if trace && !run.is_empty() {
            eprintln!(
                "      [selpoly] he {raw} (used {h}) src {:?} {kind} first=({:.3},{:.3}) last=({:.3},{:.3})",
                graph.half_edges[h].source,
                run[0].0,
                run[0].1,
                run[run.len() - 1].0,
                run[run.len() - 1].1
            );
        }
        last_uv = run.last().copied();
        runs.push((raw, run));
    }
    Some(runs)
}

/// Sample a point ε off the intersection curve, on the fragment's own
/// side, in the owning face's chart. Returns `None` when the loop has
/// no Shared half-edge, cannot be charted, or no offset candidate
/// verifies inside the loop polygon.
pub(crate) fn fragment_cut_offset_sample(
    graph: &IntersectionGraph,
    fragment: &crate::merge::walk::Fragment,
    surface: &SurfaceGeometry,
) -> Option<Point3> {
    let runs = loop_uv_runs(graph, &fragment.face, &fragment.outer_loop, surface)?;
    let poly: Vec<(f64, f64)> = runs.iter().flat_map(|(_, r)| r.iter().copied()).collect();
    if poly.len() < 3 {
        return None;
    }
    // Loop winding via shoelace: interior is left of travel for CCW,
    // right for CW — self-consistent with the same polygon used for
    // verification, so no chart-orientation convention is assumed.
    let mut area2 = 0.0;
    for i in 0..poly.len() {
        let (u0, v0) = poly[i];
        let (u1, v1) = poly[(i + 1) % poly.len()];
        area2 += u0 * v1 - u1 * v0;
    }
    let (mut lo_u, mut hi_u, mut lo_v, mut hi_v) = (f64::MAX, f64::MIN, f64::MAX, f64::MIN);
    for &(u, v) in &poly {
        lo_u = lo_u.min(u);
        hi_u = hi_u.max(u);
        lo_v = lo_v.min(v);
        hi_v = hi_v.max(v);
    }
    let diag = ((hi_u - lo_u).powi(2) + (hi_v - lo_v).powi(2)).sqrt();
    if diag <= 0.0 || area2.abs() < 1e-12 * diag * diag {
        return None;
    }
    // Longest Shared run carries the most reliable local tangent.
    let run_len = |r: &[(f64, f64)]| -> f64 {
        r.windows(2)
            .map(|w| ((w[1].0 - w[0].0).powi(2) + (w[1].1 - w[0].1).powi(2)).sqrt())
            .sum()
    };
    let (_, shared_run) = runs
        .iter()
        .filter(|(h, r)| {
            matches!(graph.half_edges[*h].source, EdgeSource::Shared { .. }) && r.len() >= 2
        })
        .max_by(|a, b| run_len(&a.1).total_cmp(&run_len(&b.1)))?;
    let m = shared_run.len() / 2;
    let (pu, pv) = shared_run[m.saturating_sub(1)];
    let (qu, qv) = shared_run[(m + 1).min(shared_run.len() - 1)];
    let (mu, mv) = shared_run[m];
    let t_len = ((qu - pu).powi(2) + (qv - pv).powi(2)).sqrt();
    if t_len < 1e-15 {
        return None;
    }
    let left = (-(qv - pv) / t_len, (qu - pu) / t_len);
    let first_dir = if area2 > 0.0 { 1.0 } else { -1.0 };
    let even_odd_inside = |u: f64, v: f64| -> bool {
        let mut inside = false;
        for i in 0..poly.len() {
            let (u0, v0) = poly[i];
            let (u1, v1) = poly[(i + 1) % poly.len()];
            if (v0 > v) != (v1 > v) && u < u0 + (v - v0) / (v1 - v0) * (u1 - u0) {
                inside = !inside;
            }
        }
        inside
    };
    for eps in [1e-3, 3e-3, 1e-2] {
        for dir in [first_dir, -first_dir] {
            let u = mu + dir * left.0 * eps * diag;
            let v = mv + dir * left.1 * eps * diag;
            if even_odd_inside(u, v) {
                if std::env::var_os("PK_SEL_TRACE").is_some() {
                    eprintln!(
                        "      [cut-offset] uv=({u:.4},{v:.4}) eps={eps} dir={dir} area2={area2:.3}"
                    );
                }
                return Some(surface.point_at(u, v));
            }
        }
    }
    None
}
