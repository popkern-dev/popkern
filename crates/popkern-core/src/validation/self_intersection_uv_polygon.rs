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
        use std::f64::consts::TAU;
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
        let mut raw_samples: Vec<RawLoopSample> = Vec::new();
        let mut current = start_he;
        let walk_bound = self.half_edges.len().max(1);
        for _ in 0..=walk_bound {
            let he = self.half_edges.get(current)?;
            if let Some(pcurve) = &he.pcurve {
                for sample in 0..samples_per_pcurve {
                    let t = pcurve.t_start
                        + (pcurve.t_end - pcurve.t_start)
                            * (sample as f64 / samples_per_pcurve as f64);
                    let uv = pcurve.point_at(t);
                    let raw = UvPoint { u: uv.x, v: uv.y };
                    let pole_jump = sample == 0
                        && periodic_u
                        && raw_samples.last().is_some_and(|prev| {
                            (raw.u - prev.uv.u).abs() > TAU * 0.5
                                && degenerate_row_witness(surface, raw.v, true)
                        });
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
