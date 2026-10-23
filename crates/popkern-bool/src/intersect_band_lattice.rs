//! Record splits at the period lattice of multi-turn bands (PK-919).
//!
//! A thread root is a helical band on a cylinder: one face winding
//! several turns between two helices, with a seam edge only in the
//! unthreaded strip below the thread start. The thread builder splits
//! the helices at every crossing of the seam's u, so the chart walk's
//! period-lattice barriers (PK-464) cut the band into per-turn pieces,
//! each within one period. A record crossing that lattice line ABOVE
//! the seam edge gets no node there, because no edge exists to split it
//! on. It joins the pieces on both sides into one component spanning two
//! periods, and the face falls back to the 3D-Newell walk (measured: the
//! tap tool's flat-end cut, u -0.781 -> +0.781 at v 16.3, span 12.57;
//! the angled fitting's cap screw root, span 41.4).
//!
//! The pass splits such records where their own uv curve on the band
//! crosses u = seam_u (mod P) outside the seam's v-range, so the walk
//! finds a lattice node there like the helices'. The split vertex is
//! degree 2 on its result edge, and post-processing's redundant-chain
//! pass dissolves it. Only bands whose seam spans less than their v
//! extent take the pass; an ordinary wall's seam spans its height and
//! already splits its records.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::topology::{BRep, FaceHandle};

use crate::clip::clip_curve_to_crossings;
use crate::intersect::{EndpointBoundaryType, IntersectionRecord};

/// Samples per record when looking for lattice crossings.
const CROSSING_SAMPLES: usize = 64;

/// A face whose seam spans only part of its height: the seam's u,
/// the chart period, and the seam's v-range.
#[derive(Debug, Clone, Copy)]
struct PartialSeam {
    seam_u: f64,
    period: f64,
    seam_v: (f64, f64),
}

/// The face's partial seam, when it is a u-periodic analytic face
/// whose self-twinned seam covers less than its v extent.
fn partial_seam(brep: &BRep, face: FaceHandle) -> Option<PartialSeam> {
    let f = brep.faces.get(face)?;
    let surface = &f.surface.as_ref()?.geometry;
    if matches!(surface, SurfaceGeometry::Nurbs(_)) {
        return None;
    }
    let period = surface.periodicity().0?;
    let ((_, _), (v0, v1)) = brep.face_trim_uv_window(face)?;
    let mut seam: Option<(f64, f64, f64)> = None;
    for l in std::iter::once(f.outer_loop).chain(f.inner_loops.iter().copied()) {
        let start = brep.loops.get(l)?.half_edge;
        let mut h = start;
        for _ in 0..=brep.half_edges.len() {
            let he = brep.half_edges.get(h)?;
            let twin_here = brep
                .half_edges
                .get(he.twin)
                .and_then(|t| brep.loops.get(t.loop_))
                .is_some_and(|tl| tl.face == face);
            if twin_here && let Some(pc) = he.pcurve.as_ref() {
                let (t0, t1) = pc.parameter_range();
                let (a, b) = (pc.point_at(t0), pc.point_at(t1));
                if (a.x - b.x).abs() < period * 1e-6 && (a.y - b.y).abs() > 0.0 {
                    let (lo, hi) = (a.y.min(b.y), a.y.max(b.y));
                    seam = Some(match seam {
                        Some((u, l0, h0)) => (u, l0.min(lo), h0.max(hi)),
                        None => (a.x.rem_euclid(period), lo, hi),
                    });
                }
            }
            h = he.next;
            if h == start {
                break;
            }
        }
    }
    let (seam_u, lo, hi) = seam?;
    let margin = 1e-6 * (v1 - v0).abs().max(1e-12);
    (lo > v0 + margin || hi < v1 - margin).then_some(PartialSeam {
        seam_u,
        period,
        seam_v: (lo, hi),
    })
}

/// Parameters where `uv` crosses u = seam_u (mod P) at a v outside the
/// seam's range, strictly inside the record.
fn lattice_crossings(uv: &popkern_core::geometry::NurbsCurve, band: PartialSeam) -> Vec<f64> {
    let (t0, t1) = uv.parameter_range();
    let cell = |t: f64| ((uv.point_at(t).x - band.seam_u) / band.period).floor();
    let mut out = Vec::new();
    for k in 0..CROSSING_SAMPLES {
        let (a, b) = (
            t0 + (t1 - t0) * k as f64 / CROSSING_SAMPLES as f64,
            t0 + (t1 - t0) * (k + 1) as f64 / CROSSING_SAMPLES as f64,
        );
        if cell(a) == cell(b) {
            continue;
        }
        // The lattice line between the two cells.
        let line = band.seam_u + band.period * cell(a).max(cell(b));
        let (mut lo, mut hi) = (a, b);
        let side = |t: f64| uv.point_at(t).x >= line;
        for _ in 0..60 {
            let m = 0.5 * (lo + hi);
            if side(m) == side(lo) {
                lo = m;
            } else {
                hi = m;
            }
        }
        let t = 0.5 * (lo + hi);
        let v = uv.point_at(t).y;
        let inside = (t - t0).min(t1 - t) > (t1 - t0) * 1e-6;
        if inside && (v < band.seam_v.0 || v > band.seam_v.1) {
            out.push(t);
        }
    }
    out
}

/// Split every record at its lattice crossings on a multi-turn band of
/// either operand (see the module docs). Records on other faces pass
/// through unchanged.
pub(crate) fn split_records_at_band_lattice(
    records: Vec<IntersectionRecord>,
    a: &BRep,
    b: &BRep,
) -> Vec<IntersectionRecord> {
    let mut bands: std::collections::HashMap<(bool, FaceHandle), Option<PartialSeam>> =
        std::collections::HashMap::new();
    let mut out = Vec::with_capacity(records.len());
    for record in records {
        if record.curve.is_closed {
            out.push(record);
            continue;
        }
        let mut ts: Vec<f64> = Vec::new();
        for (on_a, face, uv) in [
            (true, record.face_a, &record.curve.curve_uv1),
            (false, record.face_b, &record.curve.curve_uv2),
        ] {
            let band = *bands
                .entry((on_a, face))
                .or_insert_with(|| partial_seam(if on_a { a } else { b }, face));
            if let Some(band) = band {
                ts.extend(lattice_crossings(uv, band));
            }
        }
        if ts.is_empty() {
            out.push(record);
            continue;
        }
        ts.sort_by(f64::total_cmp);
        let (t0, t1) = record.curve.curve_3d.parameter_range();
        let mut cuts = vec![t0];
        cuts.extend(ts);
        cuts.push(t1);
        let pieces: Option<Vec<IntersectionRecord>> = cuts
            .windows(2)
            .enumerate()
            .map(|(i, w)| {
                let curve = popkern_core::geometry::curve_on_surface::CurveOnSurface {
                    curve_3d: clip_curve_to_crossings(&record.curve.curve_3d, w[0], w[1])?,
                    curve_uv1: clip_curve_to_crossings(&record.curve.curve_uv1, w[0], w[1])?,
                    curve_uv2: clip_curve_to_crossings(&record.curve.curve_uv2, w[0], w[1])?,
                    is_closed: false,
                };
                let mut piece = record.clone();
                piece.curve = curve;
                if i > 0 {
                    piece.start_boundary = EndpointBoundaryType::Interior;
                }
                if i + 2 < cuts.len() {
                    piece.end_boundary = EndpointBoundaryType::Interior;
                }
                Some(piece)
            })
            .collect();
        match pieces {
            Some(pieces) => out.extend(pieces),
            None => out.push(record),
        }
    }
    out
}
