//! Chart embeddings: branch-resolved UV images of registered half-edges.
//!
//! PK-104 substrate, stage 1 — see
//! docs/architecture/fragment-walk-chart-subdivision.md §3.4. Every
//! half-edge registered on a face gets a UV sample run in that face's
//! chart, in the half-edge's traversal direction. Sources, in order:
//! the half-edge's own pcurve (verified against its 3D curve window),
//! the twin's pcurve reversed (cross-registered directions store their
//! pcurve in the other face's chart), and mid-out 3D projection with
//! period-branch continuation (split edges drop pcurves at graph level;
//! projection is the safety net the fragment sampler already relies
//! on post-walk).

use popkern_core::geometry::curve::{CurveEval, CurveOnEdge};
use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::uv_projection::{
    project_point_to_uv_multisheet, project_point_to_uv_near,
};

use crate::merge::graph::{FaceRef, IntersectionGraph};

/// Where a chart embedding's UV run came from (feeds the coverage
/// census in the subdivision report).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum EmbeddingSource {
    /// The registered half-edge's own pcurve, on this face's chart.
    OwnPcurve,
    /// The twin's pcurve, reversed into this half-edge's direction.
    TwinPcurve,
    /// 3D curve samples projected with period-branch continuation.
    Projection,
    /// Endpoint positions projected (no curve geometry on the pair).
    EndpointProjection,
}

/// Branch-resolved UV image of one registered half-edge on one face.
#[derive(Debug, Clone)]
pub(crate) struct ChartEmbedding {
    /// UV samples in the half-edge's traversal direction on this face.
    /// At least 2 samples; endpoints are the run's chart endpoints.
    pub uv: Vec<(f64, f64)>,
    /// Which source produced the run.
    pub source: EmbeddingSource,
    /// True when both directions of the pair are registered on this
    /// face with this face as their primary owner — a chart seam pair.
    pub seam_like: bool,
}

impl ChartEmbedding {
    /// Chart departure direction at the run start (first non-degenerate
    /// segment), or `None` when the run is a single point.
    pub fn start_dir(&self) -> Option<(f64, f64)> {
        first_dir(self.uv.iter())
    }

    /// Chart arrival direction at the run end (direction of travel of
    /// the last non-degenerate segment).
    pub fn end_dir(&self) -> Option<(f64, f64)> {
        first_dir(self.uv.iter().rev()).map(|(du, dv)| (-du, -dv))
    }

    /// Translate the whole run in u by `delta` (period-branch move).
    pub fn shift_u(&mut self, delta: f64) {
        for uv in &mut self.uv {
            uv.0 += delta;
        }
    }
}

fn first_dir<'a>(mut it: impl Iterator<Item = &'a (f64, f64)>) -> Option<(f64, f64)> {
    let &(u0, v0) = it.next()?;
    for &(u, v) in it {
        let (du, dv) = (u - u0, v - v0);
        let len = (du * du + dv * dv).sqrt();
        if len > 1e-13 {
            return Some((du / len, dv / len));
        }
    }
    None
}

/// Sample count per run (inclusive endpoints — the endpoints ARE the
/// chart incidence data the subdivision welds on).
const N_SAMPLES: usize = 16;

/// Build the chart embedding of registered half-edge `raw` on `face`.
///
/// Returns `None` when no source can produce a run — an I1 (embedding
/// totality) failure the caller records and routes to fallback.
pub(crate) fn build_embedding(
    graph: &IntersectionGraph,
    face: &FaceRef,
    raw: usize,
    surface: &SurfaceGeometry,
) -> Option<ChartEmbedding> {
    let twin = graph.half_edges[raw].twin;
    // Pcurves live in their owner face's chart: prefer the direction
    // assigned to THIS face (same resolution the fragment sampler
    // uses; see cellular/fragment_sampling.rs).
    let h = if graph.half_edges[raw].face.as_ref() == Some(face) {
        raw
    } else if graph.half_edges[twin].face.as_ref() == Some(face) {
        twin
    } else {
        raw
    };
    let seam_like = graph.half_edges[raw].face.as_ref() == Some(face)
        && graph.half_edges[graph.half_edges[raw].twin].face.as_ref() == Some(face);

    let curve = graph.half_edges[h]
        .curve
        .as_ref()
        .or_else(|| graph.half_edges[graph.half_edges[h].twin].curve.as_ref());
    let curve_is_own = graph.half_edges[h].curve.is_some();

    let mut run: Option<(Vec<(f64, f64)>, EmbeddingSource)> = None;

    if let Some(pc) = graph.half_edges[h].pcurve.as_ref()
        && let Some(samples) = pcurve_run(pc, curve, curve_is_own, surface)
    {
        let source = if h == raw {
            EmbeddingSource::OwnPcurve
        } else {
            EmbeddingSource::TwinPcurve
        };
        run = Some((samples, source));
    }
    if run.is_none()
        && let Some(c) = curve
        && let Some(samples) = projection_run(c, curve_is_own, surface)
    {
        run = Some((samples, EmbeddingSource::Projection));
    }
    if run.is_none() {
        // No curve geometry at all: project the two node positions
        // (same nearest-branch snap as `projection_run` — the `_near`
        // helper ignores hints on analytic surfaces).
        let a = graph.nodes[graph.half_edges[h].origin].position;
        let b = graph.nodes[graph.half_edges[graph.half_edges[h].twin].origin].position;
        // Sheet-correct anchor (same rationale as `projection_run`).
        let ua = project_point_to_uv_multisheet(surface, &a)?;
        let (mut ub, mut vb) = project_point_to_uv_near(surface, &b, ua.0, ua.1)
            .or_else(|| project_point_to_uv_multisheet(surface, &b))?;
        let (pp_u, pp_v) = surface.parameter_periodicity();
        if let Some(pp) = pp_u {
            ub += ((ua.0 - ub) / pp).round() * pp;
        }
        if let Some(pv) = pp_v {
            vb += ((ua.1 - vb) / pv).round() * pv;
        }
        run = Some((vec![ua, (ub, vb)], EmbeddingSource::EndpointProjection));
    }

    let (mut uv, source) = run?;
    // The resolved half-edge `h` runs opposite the registered `raw`
    // when the twin owned this face's pcurve/curve — reverse into
    // `raw`'s traversal direction (the fragment sampler's measured
    // zigzag lesson).
    if h != raw {
        uv.reverse();
    }
    if uv.len() < 2 {
        return None;
    }
    // The graph does NOT guarantee a stored curve/pcurve runs in its
    // half-edge's traversal direction — split shared pairs violate it
    // (measured: all four breached-bore mouth arcs embedded backwards,
    // which merged the seam crossing's two chart branches into one and
    // reproduced the 3D-Newell mega-loop). The 3D-Newell walker measures
    // which curve end sits at the node instead of trusting direction;
    // do the same: the graph nodes are the orientation authority.
    let p_start = graph.nodes[graph.half_edges[raw].origin].position;
    let p_end = graph.nodes[graph.half_edges[twin].origin].position;
    let first3 = surface.point_at(uv[0].0, uv[0].1);
    let last3 = surface.point_at(uv[uv.len() - 1].0, uv[uv.len() - 1].1);
    let as_is = first3.distance_to(&p_start) + last3.distance_to(&p_end);
    let flipped = first3.distance_to(&p_end) + last3.distance_to(&p_start);
    if flipped < as_is {
        uv.reverse();
    }
    if popkern_core::trace::trace_enabled("branch-pairing") {
        let (ulo, uhi) = uv.iter().fold((f64::MAX, f64::MIN), |(lo, hi), &(u, _)| {
            (lo.min(u), hi.max(u))
        });
        eprintln!(
            "  [embed] he {raw} (resolved {h}) source {source:?} n {} u[{ulo:.4},{uhi:.4}] first {:?} last {:?}",
            uv.len(),
            uv.first(),
            uv.last()
        );
    }
    Some(ChartEmbedding {
        uv,
        source,
        seam_like,
    })
}

/// Sample a pcurve into a run, verified against the half-edge's 3D
/// curve so a stale full-range pcurve on a split half-edge cannot
/// produce a wrong-span embedding. Tries the 3D curve's parameter
/// window directly on the pcurve first (split edges share the parent's
/// parameter space when consistent), then the pcurve's own full range;
/// each candidate must 3D-verify at both ends and the middle. Returns
/// `None` when neither verifies (caller falls through to projection).
fn pcurve_run(
    pc: &Pcurve,
    curve: Option<&CurveOnEdge>,
    curve_is_own: bool,
    surface: &SurfaceGeometry,
) -> Option<Vec<(f64, f64)>> {
    let verify = |s0: f64, s1: f64| -> bool {
        let Some(c) = curve else {
            // No 3D curve to verify against: accept the pcurve range.
            return true;
        };
        // The stored 3D curve runs in the OWNER's direction; when it
        // was borrowed from the twin the parameter window maps
        // reversed relative to `pc`'s owner. Both orders are checked.
        let (t0, t1) = if curve_is_own {
            (c.t_start, c.t_end)
        } else {
            (c.t_end, c.t_start)
        };
        let p3 = |t: f64| c.geometry.point_at(t);
        let span = p3(t0).distance_to(&p3(t1)).max(1e-9);
        let tol = 1e-6 * (1.0 + span);
        for (fs, ft) in [(0.0, 0.0), (0.5, 0.5), (1.0, 1.0)] {
            let uv = pc.point_at(s0 + (s1 - s0) * fs);
            let on_surface = surface.point_at(uv.x, uv.y);
            let on_curve = p3(t0 + (t1 - t0) * ft);
            if on_surface.distance_to(&on_curve) > tol {
                return false;
            }
        }
        true
    };

    let mut window: Option<(f64, f64)> = None;
    if let Some(c) = curve {
        let (w0, w1) = if curve_is_own {
            (c.t_start, c.t_end)
        } else {
            (c.t_end, c.t_start)
        };
        // A DEGENERATE 3D window (a zero-length pole edge — the
        // chart's collapsed v-border) trivially verifies at a single
        // point and would collapse the run; the pcurve's own full
        // range carries the real border span (measured: a sphere's
        // pole edge pcurve spans the whole u-line while its 3D curve
        // has zero extent).
        let span3 = c
            .geometry
            .point_at(c.t_start)
            .distance_to(&c.geometry.point_at(c.t_end));
        if span3 > 1e-9 && verify(w0, w1) {
            window = Some((w0, w1));
        }
    }
    if window.is_none() && verify(pc.t_start, pc.t_end) {
        window = Some((pc.t_start, pc.t_end));
    }
    let (s0, s1) = window?;
    let mut run = Vec::with_capacity(N_SAMPLES + 3);
    for i in 0..=N_SAMPLES {
        let t = s0 + (s1 - s0) * i as f64 / N_SAMPLES as f64;
        let p = pc.point_at(t);
        run.push((p.x, p.y));
        if i == 0 {
            // ENDPOINT-REFINED direction sample (mirrored before the
            // final point below): `start_dir`/`end_dir` are chords of
            // the first/last run segment, and a 1/16-run chord on a
            // curved record is up to ~0.06 rad off the true tangent —
            // more than the angular separation of a thread groove
            // meeting a rim near-tangentially (0.013 rad measured,
            // PK-460), so the walk's wedge matching at such junctions
            // ordered the slots wrong and fused the two adjacent
            // regions into one orbit. A δ = span/1000 neighbour makes
            // the chord direction honest at the only places it is
            // consumed.
            let td = s0 + (s1 - s0) * 1.0e-3;
            let pd = pc.point_at(td);
            run.push((pd.x, pd.y));
        }
        if i == N_SAMPLES - 1 {
            let td = s0 + (s1 - s0) * (1.0 - 1.0e-3);
            let pd = pc.point_at(td);
            run.push((pd.x, pd.y));
        }
    }
    Some(run)
}

/// Project the 3D curve into the chart, middle sample first so the
/// period branch is anchored by an unambiguous interior point (curve
/// endpoints often sit ON the chart seam, where plain projection may
/// return either border), then continued outward to both ends.
///
/// Branch continuity is enforced HERE: `project_point_to_uv_near`
/// ignores the hint for analytic surfaces (it exists to stop NURBS
/// Gauss-Newton branch jumping), so every projected sample is snapped
/// to the period branch nearest its already-projected neighbour —
/// measured on the breached-bore mouth arc, whose seam-node endpoint
/// projected to the principal border and put the run two periods wide.
fn projection_run(
    c: &CurveOnEdge,
    curve_is_own: bool,
    surface: &SurfaceGeometry,
) -> Option<Vec<(f64, f64)>> {
    let (u_period, v_period) = surface.parameter_periodicity();
    let (t0, t1) = if curve_is_own {
        (c.t_start, c.t_end)
    } else {
        (c.t_end, c.t_start)
    };
    let ts: Vec<f64> = (0..=N_SAMPLES)
        .map(|i| t0 + (t1 - t0) * i as f64 / N_SAMPLES as f64)
        .collect();
    let mid = N_SAMPLES / 2;
    let mut run = vec![(f64::NAN, f64::NAN); ts.len()];
    let p_mid = c.geometry.point_at(ts[mid]);
    // The anchor must be SHEET-CORRECT: the cold 5-start projector
    // walks the adjacent coil turn on many-turn closed NURBS walls
    // (the PK-446 class), and every subsequent sample is continued
    // NEAR its neighbour — a wrong-sheet anchor glues the whole run
    // to the wrong turn (measured, PK-452: a seam-crossing cut
    // oval's run projected onto the turn at v ≈ 0.55 instead of
    // 0.86, folded to a sub-period span, both endpoints clustered
    // into ONE seam branch, and the walk isolated the oval instead
    // of splitting the wall into bands).
    run[mid] = project_point_to_uv_multisheet(surface, &p_mid)?;
    let proj_near = |p: &popkern_core::geometry::Point3, seed: (f64, f64)| -> Option<(f64, f64)> {
        let (mut u, mut v) = project_point_to_uv_near(surface, p, seed.0, seed.1)
            .or_else(|| project_point_to_uv_multisheet(surface, p))?;
        if let Some(pp) = u_period {
            u += ((seed.0 - u) / pp).round() * pp;
        }
        // Unwrap v too: a doubly-periodic (torus) run crossing the
        // v-seam tore mid-run without it — the torn embedding flipped
        // its face's boundary shoelace and σ mis-calibrated (measured:
        // the seam-adjacent quarter of a torus operand walked σ = −1
        // and its whole partition classified as one Inside fragment).
        if let Some(pv) = v_period {
            v += ((seed.1 - v) / pv).round() * pv;
        }
        Some((u, v))
    };
    for i in (mid + 1)..ts.len() {
        let p = c.geometry.point_at(ts[i]);
        run[i] = proj_near(&p, run[i - 1])?;
    }
    for i in (0..mid).rev() {
        let p = c.geometry.point_at(ts[i]);
        run[i] = proj_near(&p, run[i + 1])?;
    }
    // Endpoint-refined direction samples (same rationale as
    // `pcurve_run`): δ = span/1000 neighbours next to each end make
    // `start_dir`/`end_dir` honest tangents instead of 1/16-run
    // chords.
    let td0 = t0 + (t1 - t0) * 1.0e-3;
    let td1 = t0 + (t1 - t0) * (1.0 - 1.0e-3);
    let d0 = proj_near(&c.geometry.point_at(td0), run[0])?;
    let d1 = proj_near(&c.geometry.point_at(td1), run[ts.len() - 1])?;
    run.insert(1, d0);
    run.insert(run.len() - 1, d1);
    Some(run)
}
