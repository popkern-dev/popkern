//! Geometric deviation measurement for BRep edges and vertices.
//!
//! Provides methods to measure actual geometric deviations between
//! a BRep's topological entities and their underlying geometry.

use super::brep::BRep;
use crate::geometry::Point3;
use crate::geometry::curve::{CurveEval, CurveOnEdge};
use crate::geometry::pcurve::Pcurve;
use crate::geometry::pcurve_trace::{TraceNearest, nearest_on_trace as pcurve_trace_nearest};
use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};
use crate::topology::brep_edge_surface_compatibility_helpers::unique_edge_half_edges;
use crate::topology::{EdgeHandle, FaceHandle, HalfEdgeHandle, VertexHandle};

/// Number of sample points along an edge curve for deviation measurement.
const DEFAULT_EDGE_SAMPLES: usize = 20;

/// Coarse bracketing samples for the refined trace-distance measurement.
const TRACE_COARSE_SAMPLES: usize = 128;

/// Per-edge report describing how well adjacent pcurves reproduce a 3D edge.
#[derive(Debug, Clone, PartialEq)]
pub struct EdgePcurveDeviationReport {
    /// Edge whose 3D curve was sampled.
    pub edge: EdgeHandle,
    /// Maximum accepted geometric disagreement over all adjacent half-edges.
    pub max_deviation: f64,
    /// Half-edge pcurve measurements contributing to this report.
    pub half_edges: Vec<HalfEdgePcurveDeviation>,
}

/// Per-half-edge pcurve agreement against the owning edge curve and face surface.
#[derive(Debug, Clone, PartialEq)]
pub struct HalfEdgePcurveDeviation {
    /// Half-edge whose pcurve was measured.
    pub half_edge: HalfEdgeHandle,
    /// Adjacent face that owns the half-edge pcurve.
    pub face: FaceHandle,
    /// Maximum accepted disagreement. When the trace direction agrees
    /// this is the minimum of same-fraction and trace-distance
    /// disagreement at each sample, maxed over samples; otherwise the
    /// trace measure is rejected and only same-fraction counts.
    pub max_deviation: f64,
    /// Maximum disagreement when edge and pcurve are evaluated at the same fraction.
    pub max_same_fraction_deviation: f64,
    /// Maximum refined distance from edge sample points to the pcurve trace.
    pub max_trace_deviation: f64,
    /// Whether the nearest-trace parameters advance forward along the
    /// oriented half-edge (seam wraps on closed traces allowed). When
    /// false the pcurve traces the right point set in the wrong
    /// direction — a reversed pcurve — and trace agreement must not
    /// excuse it (PK-287: direction honesty is a designed contract,
    /// not an accident of sampling resolution).
    pub trace_direction_agrees: bool,
    /// Fraction along the oriented half-edge at which `max_deviation` occurred.
    pub worst_fraction: f64,
    /// Oriented 3D edge parameter range used for this half-edge.
    pub edge_parameter_range: (f64, f64),
    /// Pcurve parameter range used for this half-edge.
    pub pcurve_parameter_range: (f64, f64),
}

impl BRep {
    /// Maximum deviation between an edge's 3D curve and adjacent pcurve traces.
    ///
    /// STEP and other exchange formats may carry pcurves whose parameterization
    /// differs from the 3D edge curve. This measurement accepts either
    /// same-fraction agreement or geometric trace agreement at each edge sample.
    ///
    /// Returns `0.0` if the edge has no curve or no adjacent pcurves.
    pub fn measure_edge_deviation(&self, eh: EdgeHandle) -> f64 {
        self.measure_edge_pcurve_deviation(eh)
            .map(|report| report.max_deviation)
            .unwrap_or(0.0)
    }

    /// Detailed pcurve-to-edge deviation report for a single edge.
    ///
    /// Returns `None` only when the edge handle is invalid or the edge has no
    /// 3D curve. A valid edge with no measurable adjacent pcurves returns an
    /// empty report with `max_deviation == 0.0`.
    pub fn measure_edge_pcurve_deviation(
        &self,
        eh: EdgeHandle,
    ) -> Option<EdgePcurveDeviationReport> {
        let edge = self.edges.get(eh)?;
        let curve_on_edge = edge.curve.as_ref()?;
        if !self.half_edges.contains_key(edge.half_edge) {
            return Some(EdgePcurveDeviationReport {
                edge: eh,
                max_deviation: 0.0,
                half_edges: Vec::new(),
            });
        }

        let mut max_deviation = 0.0_f64;
        let mut half_edges = Vec::new();
        for heh in unique_edge_half_edges(self, edge.half_edge) {
            let Some(report) = self.measure_half_edge_pcurve_deviation_with_curve(
                heh,
                curve_on_edge,
                edge.half_edge,
            ) else {
                continue;
            };
            max_deviation = max_deviation.max(report.max_deviation);
            half_edges.push(report);
        }

        Some(EdgePcurveDeviationReport {
            edge: eh,
            max_deviation,
            half_edges,
        })
    }

    /// Detailed pcurve-to-edge deviation for one half-edge.
    ///
    /// This is useful for diagnostics and export/import repair decisions that
    /// need face ownership and half-edge orientation instead of only the edge
    /// maximum.
    pub fn measure_half_edge_pcurve_deviation(
        &self,
        heh: HalfEdgeHandle,
    ) -> Option<HalfEdgePcurveDeviation> {
        let he = self.half_edges.get(heh)?;
        let edge = self.edges.get(he.edge)?;
        let curve_on_edge = edge.curve.as_ref()?;
        self.measure_half_edge_pcurve_deviation_with_curve(heh, curve_on_edge, edge.half_edge)
    }

    /// An upper bound on [`Self::measure_half_edge_pcurve_deviation`]'s
    /// `max_deviation` from the same samples without the trace search:
    /// the measure accepts at most each sample's same-fraction distance.
    /// A screen for callers that only act on large deviations (a
    /// 28-edge shelled cube spent 10.7 ms measuring traces it never
    /// needed — PK-826 regression).
    pub(crate) fn half_edge_pcurve_deviation_bound(&self, heh: HalfEdgeHandle) -> Option<f64> {
        let he = self.half_edges.get(heh)?;
        let edge = self.edges.get(he.edge)?;
        let curve_on_edge = edge.curve.as_ref()?;
        self.measure_half_edge_with(heh, curve_on_edge, edge.half_edge, false)
            .map(|r| r.max_same_fraction_deviation)
    }

    fn measure_half_edge_pcurve_deviation_with_curve(
        &self,
        heh: HalfEdgeHandle,
        curve_on_edge: &CurveOnEdge,
        canonical_half_edge: HalfEdgeHandle,
    ) -> Option<HalfEdgePcurveDeviation> {
        self.measure_half_edge_with(heh, curve_on_edge, canonical_half_edge, true)
    }

    fn measure_half_edge_with(
        &self,
        heh: HalfEdgeHandle,
        curve_on_edge: &CurveOnEdge,
        canonical_half_edge: HalfEdgeHandle,
        trace: bool,
    ) -> Option<HalfEdgePcurveDeviation> {
        let he = self.half_edges.get(heh)?;
        let pcurve = he.pcurve.as_ref()?;
        let loop_ = self.loops.get(he.loop_)?;
        let face = self.faces.get(loop_.face)?;
        let surface = &face.surface.as_ref()?.geometry;
        // Which way this half-edge runs the stored [t_start, t_end]:
        // on an OPEN edge its origin vertex says (producers disagree on
        // whether a descending range starts at the canonical side's
        // origin — sketch laminas store lines t [1, 0] run from t_end,
        // extruded profile arcs t [4.71, 1.57] run from t_start); a
        // closed edge has one vertex, so the canonical convention
        // decides.
        let forward = {
            use crate::geometry::curve::CurveEval;
            let a = curve_on_edge.geometry.point_at(curve_on_edge.t_start);
            let b = curve_on_edge.geometry.point_at(curve_on_edge.t_end);
            let origin = self.vertices.get(he.origin).map(|v| v.position);
            match origin {
                Some(o) if a.distance_to(&b) > self.tolerance.linear() => {
                    let (da, db) = (o.distance_to(&a), o.distance_to(&b));
                    if (da - db).abs() > 1.0e-3 * a.distance_to(&b) {
                        da < db
                    } else {
                        canonical_half_edge == heh
                    }
                }
                _ => canonical_half_edge == heh,
            }
        };
        measure_pcurve_deviation(
            heh,
            loop_.face,
            pcurve,
            surface,
            curve_on_edge,
            forward,
            DEFAULT_EDGE_SAMPLES,
            trace,
        )
    }

    /// Maximum distance from a vertex position to incident edge curves
    /// evaluated at the vertex parameter.
    ///
    /// Walks the vertex fan (outgoing half-edges) and for each incident
    /// edge that has a curve, evaluates the curve at `t_start` or `t_end`
    /// depending on which end corresponds to this vertex.
    ///
    /// Returns `0.0` if no incident edges have curves.
    pub fn measure_vertex_deviation(&self, vh: VertexHandle) -> f64 {
        let Some(vertex) = self.vertices.get(vh) else {
            return 0.0;
        };
        let pos = vertex.position;

        // Walk the vertex fan: start at outgoing, follow twin.next.
        let start = vertex.outgoing;
        if !self.half_edges.contains_key(start) {
            return 0.0;
        }

        let mut max_dev = 0.0_f64;
        let mut cursor = start;
        let bound = self.half_edges.len();
        let mut steps = 0;

        loop {
            steps += 1;
            if steps > bound {
                break; // Safety: avoid infinite loop on corrupted data.
            }

            let Some(he) = self.half_edges.get(cursor) else {
                break;
            };
            let edge_handle = he.edge;

            if let Some(edge) = self.edges.get(edge_handle)
                && let Some(ref curve_on_edge) = edge.curve
            {
                // Determine which end of the curve corresponds to this vertex.
                // The edge's canonical half-edge origin is at t_start.
                // Its twin's origin is at t_end.
                let t = if let Some(canonical_he) = self.half_edges.get(edge.half_edge) {
                    if canonical_he.origin == vh {
                        curve_on_edge.t_start
                    } else {
                        curve_on_edge.t_end
                    }
                } else {
                    cursor = self.advance_vertex_fan(cursor, start);
                    if cursor == start {
                        break;
                    }
                    continue;
                };

                let curve_pt = curve_on_edge.geometry.point_at(t);
                let dev = pos.distance_to(&curve_pt);
                max_dev = max_dev.max(dev);
            }

            cursor = self.advance_vertex_fan(cursor, start);
            if cursor == start {
                break;
            }
        }

        max_dev
    }

    /// Advance the vertex fan cursor: follow twin.next.
    ///
    /// Returns `start` if the fan is complete or a stale handle is
    /// encountered, signalling the caller to stop.
    fn advance_vertex_fan(
        &self,
        cursor: crate::topology::HalfEdgeHandle,
        start: crate::topology::HalfEdgeHandle,
    ) -> crate::topology::HalfEdgeHandle {
        let Some(he) = self.half_edges.get(cursor) else {
            return start;
        };
        let Some(twin_data) = self.half_edges.get(he.twin) else {
            return start;
        };
        twin_data.next
    }
}

fn measure_pcurve_deviation(
    half_edge: HalfEdgeHandle,
    face: FaceHandle,
    pcurve: &Pcurve,
    surface: &SurfaceGeometry,
    curve_on_edge: &CurveOnEdge,
    forward: bool,
    samples: usize,
    trace_search: bool,
) -> Option<HalfEdgePcurveDeviation> {
    let samples = samples.max(1);
    let (t_lo, t_hi) = oriented_edge_parameter_range(curve_on_edge, forward);
    let pc_range = pcurve.parameter_range();
    if !t_lo.is_finite() || !t_hi.is_finite() || !pc_range.0.is_finite() || !pc_range.1.is_finite()
    {
        return None;
    }

    let mut max_same_fraction_deviation = 0.0_f64;
    let mut max_trace_deviation = 0.0_f64;
    let mut per_sample: Vec<(f64, f64, TraceNearest)> = Vec::with_capacity(samples + 1);

    for (frac, p3d) in edge_points_for_half_edge(curve_on_edge, forward, samples) {
        let pc_t = pc_range.0 + (pc_range.1 - pc_range.0) * frac;
        let uv = pcurve.point_at(pc_t);
        let same_fraction_point = surface.point_at(uv.x, uv.y);
        let same_fraction_dev = p3d.distance_to(&same_fraction_point);
        let trace = if trace_search {
            pcurve_trace_nearest(pcurve, surface, &p3d, TRACE_COARSE_SAMPLES)
        } else {
            TraceNearest {
                distance: same_fraction_dev,
                t: pc_t,
            }
        };

        max_same_fraction_deviation = max_same_fraction_deviation.max(same_fraction_dev);
        max_trace_deviation = max_trace_deviation.max(trace.distance);
        per_sample.push((frac, same_fraction_dev, trace));
    }

    let trace_direction_agrees = trace_parameters_advance_forward(
        per_sample.iter().map(|(_, _, trace)| trace.t),
        pc_range,
        samples,
    );

    let mut max_deviation = 0.0_f64;
    let mut worst_fraction = 0.0_f64;
    for (frac, same_fraction_dev, trace) in &per_sample {
        // The pcurve's own ends must land on the edge's ends: trace
        // agreement alone accepts a pcurve whose trace merely CONTAINS
        // the edge (PK-820: a moved countersink's sub-arc edge kept its
        // parent arc's whole pcurve, 4.9 off, and its face's trim
        // polygon misread 82 of 99 sound triangles). End identity does
        // not depend on parameterization, so it is never excused.
        let at_end = *frac == 0.0 || *frac == 1.0;
        let accepted_dev = if trace_direction_agrees && !at_end {
            same_fraction_dev.min(trace.distance)
        } else {
            *same_fraction_dev
        };
        if accepted_dev > max_deviation {
            max_deviation = accepted_dev;
            worst_fraction = *frac;
        }
    }

    Some(HalfEdgePcurveDeviation {
        half_edge,
        face,
        max_deviation,
        max_same_fraction_deviation,
        max_trace_deviation,
        worst_fraction,
        trace_direction_agrees,
        edge_parameter_range: (t_lo, t_hi),
        pcurve_parameter_range: pc_range,
    })
}

/// Whether nearest-trace parameters advance forward along the oriented
/// half-edge samples. "Forward" is the pcurve's own declared direction:
/// deltas are oriented by `sign(t_end - t_start)`, so descending
/// parameter ranges (twin-side and STEP descending-domain pcurves,
/// PK-215) are first-class.
///
/// Oriented deltas are classified against the pcurve's parameter range
/// (scale-free, no global tolerance): |delta| within
/// `range / (4 · samples)` is projection noise; a backward step of at
/// least half the range is a seam wrap on a closed trace (the nearest
/// parameter legitimately jumps from one end to the other); any other
/// backward step means the pcurve is traversed against the half-edge
/// direction and trace agreement must be rejected. A fully reversed
/// pcurve steps backward by ~range/samples at every sample — squarely
/// in the rejected band at any model scale.
fn trace_parameters_advance_forward(
    trace_ts: impl Iterator<Item = f64>,
    pc_range: (f64, f64),
    samples: usize,
) -> bool {
    let span = pc_range.1 - pc_range.0;
    let range = span.abs();
    // Degenerate or NaN range: no direction to disagree with.
    if range <= 0.0 || range.is_nan() {
        return true;
    }
    let orientation = span.signum();
    let noise_slack = range / (4.0 * samples.max(1) as f64);
    let wrap_threshold = range * 0.5;
    let mut previous: Option<f64> = None;
    for t in trace_ts {
        if let Some(prev) = previous {
            let delta = (t - prev) * orientation;
            if delta < -noise_slack && delta > -wrap_threshold {
                return false;
            }
        }
        previous = Some(t);
    }
    true
}

fn oriented_edge_parameter_range(curve_on_edge: &CurveOnEdge, forward: bool) -> (f64, f64) {
    if forward {
        (curve_on_edge.t_start, curve_on_edge.t_end)
    } else {
        (curve_on_edge.t_end, curve_on_edge.t_start)
    }
}

fn edge_points_for_half_edge(
    curve_on_edge: &CurveOnEdge,
    forward: bool,
    samples: usize,
) -> Vec<(f64, Point3)> {
    let samples = samples.max(1);
    let (t_lo, t_hi) = oriented_edge_parameter_range(curve_on_edge, forward);
    (0..=samples)
        .map(|i| {
            let frac = i as f64 / samples as f64;
            let t = t_lo + (t_hi - t_lo) * frac;
            (frac, curve_on_edge.geometry.point_at(t))
        })
        .collect()
}

#[cfg(test)]
#[path = "brep_measurement_tests.rs"]
mod brep_measurement_tests;
