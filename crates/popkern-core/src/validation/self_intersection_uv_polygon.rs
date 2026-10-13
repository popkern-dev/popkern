//! Loop UV polygon assembly for trim classification: raw pcurve
//! sampling with adaptive density, degenerate-row (pole/apex) branch
//! flags, and the full-period seam handling. Split from
//! `self_intersection_trim.rs` (500-line ratchet, PK-328).

use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};
use crate::topology::{BRep, LoopHandle};

use super::self_intersection_uv_containment::{RawLoopSample, UvPoint, assemble_periodic_polygon};

impl BRep {
    pub(in crate::validation) fn loop_uv_polygon(
        &self,
        loop_handle: LoopHandle,
        surface: &SurfaceGeometry,
        query: UvPoint,
    ) -> Option<Vec<UvPoint>> {
        let loop_ = self.loops.get(loop_handle)?;
        let start_he = loop_.half_edge;
        if !self.half_edges.contains_key(start_he) {
            return None;
        }

        // Collect RAW pcurve samples plus a flag marking loop positions
        // where the traversal crosses a parametrically DEGENERATE row
        // (sphere pole, cone apex) with a period-sized branch jump.
        // Assembly — chain-continuous unwrapping, with a raw-branch
        // fallback for full-wrap seam loops — lives in
        // `assemble_periodic_polygon` (see its doc for the measured
        // failure modes on both sides of this trade).
        //
        // Sampling density adapts to the loop's half-edge count: a
        // one-pcurve loop (a closed hole rim) sampled 6× is an
        // INSCRIBED HEXAGON whose inradius under-covers the hole by
        // 13% — a bore surface passing inside its own skin's hole rim
        // classified Interior and the phantom "crossing" was trimmed
        // into a served body (PK-328, measured on the bracket hollow
        // at thin offsets).
        const SAMPLES_PER_PCURVE: usize = 6;
        const MIN_LOOP_SAMPLES: usize = 24;
        let (periodic_u, periodic_v) = surface.is_periodic();
        let mut loop_len = 0usize;
        {
            let mut cursor = start_he;
            for _ in 0..=self.half_edges.len().max(1) {
                let he = self.half_edges.get(cursor)?;
                loop_len += 1;
                cursor = he.next;
                if cursor == start_he {
                    break;
                }
            }
        }
        let samples_per_pcurve = SAMPLES_PER_PCURVE.max(MIN_LOOP_SAMPLES / loop_len.max(1));
        // The classifier's parametric boundary band (`face_trim_classify_uv`).
        let band = self
            .tolerance
            .parametric_base()
            .max(self.tolerance.linear())
            .max(1.0e-9)
            * 10.0;
        let mut raw_samples: Vec<RawLoopSample> = Vec::new();
        let mut current = start_he;
        let walk_bound = self.half_edges.len().max(1);
        for _ in 0..=walk_bound {
            let he = self.half_edges.get(current)?;
            if let Some(pcurve) = &he.pcurve {
                let stations = query_refined_stations(
                    pcurve,
                    samples_per_pcurve,
                    query,
                    band,
                    surface.parameter_periodicity(),
                );
                for (sample, t) in stations.into_iter().enumerate() {
                    let uv = pcurve.point_at(t);
                    let raw = UvPoint { u: uv.x, v: uv.y };
                    // A chain restarts wherever a pcurve leaves a
                    // degenerate row (its start is the previous pcurve's
                    // end): across a period-sized branch jump, or on the
                    // same raw branch (a seam walked down to a pole and
                    // back up). Chains only add assembly variants for a
                    // folded polygon that tears or has no area.
                    let pole_jump = sample == 0
                        && periodic_u
                        && !raw_samples.is_empty()
                        && degenerate_row_witness(surface, raw.v, true);
                    raw_samples.push(RawLoopSample {
                        uv: raw,
                        pole_jump_before: pole_jump,
                    });
                }
            } else {
                let vertex = self.vertices.get(he.origin)?;
                let (u, v) =
                    crate::geometry::uv_projection::project_point_to_uv(surface, &vertex.position)?;
                raw_samples.push(RawLoopSample {
                    uv: UvPoint { u, v },
                    pole_jump_before: false,
                });
            }

            current = he.next;
            if current == start_he {
                break;
            }
        }

        // A period-sized closure jump is 3D-free when EITHER closure
        // endpoint sits on a degenerate row (its u is arbitrary there);
        // with both endpoints on regular rows the jump is a real tear.
        let closure_on_degenerate_rows = periodic_u
            && raw_samples.len() >= 2
            && (degenerate_row_witness(surface, raw_samples[0].uv.v, true)
                || raw_samples
                    .last()
                    .is_some_and(|s| degenerate_row_witness(surface, s.uv.v, true)));
        let polygon = assemble_periodic_polygon(
            &raw_samples,
            query,
            periodic_u,
            periodic_v,
            closure_on_degenerate_rows,
        );

        if polygon.len() >= 3 {
            Some(polygon)
        } else {
            None
        }
    }
}

/// Pcurve parameters for the loop polygon: `samples` uniform stations
/// over `[t_start, t_end]`, plus bisection of every span whose chord
/// passes within twice its own sagitta of `query`, until that sagitta
/// is under a quarter of `band`. Six chords per pcurve leave a
/// quarter arc's chord 0.86% of its radius inside the arc, and every
/// point in that sliver read Outside (PK-792: half the display
/// triangles of a rounded box's bottom face, all of them sound). The
/// refinement only touches spans near the query, so the polygon is
/// exact where the answer is decided and as cheap as before elsewhere.
/// Distances wrap by the chart's periods.
fn query_refined_stations(
    pcurve: &crate::geometry::pcurve::Pcurve,
    samples: usize,
    query: UvPoint,
    band: f64,
    periods: (Option<f64>, Option<f64>),
) -> Vec<f64> {
    const MAX_DEPTH: u32 = 24;
    let at = |t: f64| {
        let p = pcurve.point_at(t);
        (p.x, p.y)
    };
    let wrap = |d: f64, period: Option<f64>| match period {
        Some(p) if p > 0.0 => d - p * (d / p).round(),
        _ => d,
    };
    let near = |a: (f64, f64), b: (f64, f64), sagitta: f64| {
        // Query relative to a, wrapped to a's branch.
        let q = (
            wrap(query.u - a.0, periods.0),
            wrap(query.v - a.1, periods.1),
        );
        let d = (b.0 - a.0, b.1 - a.1);
        let len2 = d.0 * d.0 + d.1 * d.1;
        let t = if len2 > 0.0 {
            ((q.0 * d.0 + q.1 * d.1) / len2).clamp(0.0, 1.0)
        } else {
            0.0
        };
        let (x, y) = (q.0 - d.0 * t, q.1 - d.1 * t);
        (x * x + y * y).sqrt() <= 2.0 * sagitta + band
    };
    let mut out = Vec::new();
    let span = pcurve.t_end - pcurve.t_start;
    let mut stack: Vec<(f64, f64, u32)> = Vec::new();
    for k in (0..samples).rev() {
        let t0 = pcurve.t_start + span * k as f64 / samples as f64;
        let t1 = pcurve.t_start + span * (k + 1) as f64 / samples as f64;
        stack.push((t0, t1, 0));
    }
    while let Some((t0, t1, depth)) = stack.pop() {
        let (a, b) = (at(t0), at(t1));
        let tm = 0.5 * (t0 + t1);
        let m = at(tm);
        let sagitta =
            ((m.0 - 0.5 * (a.0 + b.0)).powi(2) + (m.1 - 0.5 * (a.1 + b.1)).powi(2)).sqrt();
        if depth < MAX_DEPTH && sagitta > 0.25 * band && near(a, b, sagitta) {
            stack.push((tm, t1, depth + 1));
            stack.push((t0, tm, depth + 1));
        } else {
            out.push(t0);
        }
    }
    // The end station too: where the next pcurve starts on another
    // branch (a meridian arriving at a pole, the next leaving it at
    // another u) the corner exists only as this end, and a chord from
    // the last interior station to the next start cut it off (PK-792:
    // a die corner patch's polar triangles read Outside). Repeats of
    // the next start are dropped in assembly.
    out.push(pcurve.t_end);
    out
}

/// Whether a whole parameter row is 3D-degenerate (sphere pole, cone
/// apex): sampling the periodic axis at thirds of the period all lands
/// on one point. A seam (u=0 ≡ u=τ) does NOT pass — interior stations
/// of the period differ there.
fn degenerate_row_witness(surface: &SurfaceGeometry, fixed: f64, vary_u: bool) -> bool {
    use std::f64::consts::TAU;
    let at = |s: f64| {
        if vary_u {
            surface.point_at(s, fixed)
        } else {
            surface.point_at(fixed, s)
        }
    };
    let p0 = at(0.0);
    let p1 = at(TAU / 3.0);
    let p2 = at(2.0 * TAU / 3.0);
    let scale = 1.0
        + p0.x
            .abs()
            .max(p0.y.abs())
            .max(p0.z.abs())
            .max(p1.x.abs().max(p1.y.abs()).max(p1.z.abs()));
    let eps = 1.0e-9 * scale;
    p0.distance_to(&p1) < eps && p0.distance_to(&p2) < eps && p1.distance_to(&p2) < eps
}
