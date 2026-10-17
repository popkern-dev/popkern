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

#[path = "chart_embedding_closed.rs"]
mod chart_embedding_closed;
#[path = "chart_embedding_unfold.rs"]
mod chart_embedding_unfold;

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
    /// Endpoint-refined direction sample next to `uv[0]` (δ =
    /// span/1000 along the source curve): `start_dir`/`end_dir` are
    /// chords, and a 1/16-run chord on a curved record is up to
    /// ~0.06 rad off the true tangent — more than the angular
    /// separation of a thread groove meeting a rim near-tangentially
    /// (0.013 rad measured, PK-460), so the walk's wedge matching at
    /// such junctions ordered the slots wrong and fused two adjacent
    /// regions into one orbit. Kept OUT of `uv`: injecting δ-points
    /// into the shared polyline perturbed every other run consumer
    /// (boundary substitutes, orbit polylines, branch measures) and
    /// tore the propeller hub shell (measured: one closed free-bound
    /// wire in the OCC oracle).
    pub refined_first: Option<(f64, f64)>,
    /// Endpoint-refined direction sample next to the run's last point
    /// (see `refined_first`).
    pub refined_last: Option<(f64, f64)>,
}

/// A sampled run beside its endpoint-refined δ-samples (see
/// [`ChartEmbedding::refined_first`]).
pub(crate) type RunWithRefinedEnds = (Vec<(f64, f64)>, Option<(f64, f64)>, Option<(f64, f64)>);

/// Acceptance cone for endpoint-refined directions: the chord error
/// the refinement corrects is ≤ ~0.06 rad, so a refined direction more
/// than ~0.3 rad off the chord is an artifact, never a correction.
const REFINED_DIR_AGREEMENT: f64 = 0.95;

impl ChartEmbedding {
    /// Chart departure direction at the run start (the endpoint-refined
    /// δ-sample when it AGREES with the coarse chord, else the first
    /// non-degenerate run segment), or `None` when the run is a single
    /// point.
    ///
    /// The refinement corrects the 1/16-run chord's tangent error
    /// (≤ ~0.06 rad measured); a δ-sample whose direction leaves that
    /// cone is a projection artifact — on a closed chart the δ-point
    /// can project across the seam branch and the "tangent" points
    /// across the chart (measured, PK-460: the propeller hub's blade
    /// junction wedges scrambled and the shell tore — one closed
    /// free-bound wire). Such a sample falls back to the chord.
    pub fn start_dir(&self) -> Option<(f64, f64)> {
        let chord = first_dir(self.uv.iter());
        if let (Some(&(u0, v0)), Some((ur, vr)), Some((cu, cv))) =
            (self.uv.first(), self.refined_first, chord)
        {
            let (du, dv) = (ur - u0, vr - v0);
            let len = (du * du + dv * dv).sqrt();
            if len > 1e-13 && (du * cu + dv * cv) / len > REFINED_DIR_AGREEMENT {
                return Some((du / len, dv / len));
            }
        }
        chord
    }

    /// Chart arrival direction at the run end (direction of travel of
    /// the last non-degenerate segment; endpoint-refined when the
    /// δ-sample agrees with the coarse chord — see `start_dir`).
    pub fn end_dir(&self) -> Option<(f64, f64)> {
        let chord = first_dir(self.uv.iter().rev()).map(|(du, dv)| (-du, -dv));
        if let (Some(&(u1, v1)), Some((ur, vr)), Some((cu, cv))) =
            (self.uv.last(), self.refined_last, chord)
        {
            let (du, dv) = (u1 - ur, v1 - vr);
            let len = (du * du + dv * dv).sqrt();
            if len > 1e-13 && (du * cu + dv * cv) / len > REFINED_DIR_AGREEMENT {
                return Some((du / len, dv / len));
            }
        }
        chord
    }

    /// Translate the whole run in u by `delta` (period-branch move).
    pub fn shift_u(&mut self, delta: f64) {
        for uv in &mut self.uv {
            uv.0 += delta;
        }
        if let Some(r) = &mut self.refined_first {
            r.0 += delta;
        }
        if let Some(r) = &mut self.refined_last {
            r.0 += delta;
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

    let mut run: Option<(RunWithRefinedEnds, EmbeddingSource)> = None;

    if let Some(pc) = graph.half_edges[h].pcurve.as_ref()
        && let Some(samples) = pcurve_run(pc, curve, curve_is_own, surface)
    {
        let source = if h == raw {
            EmbeddingSource::OwnPcurve
        } else {
            EmbeddingSource::TwinPcurve
        };
        // IMAGE HONESTY: a stored record pcurve can be authored on
        // the wrong chart BRANCH — its image then sweeps the long
        // way round while the 3D curve takes the short (measured,
        // PK-464: the M5 countersunk's tip-band arc pcurve spanned
        // 4.71 rad through the +x azimuth where the band does not
        // exist; the anchored chart overwrapped 7.85 and the cone
        // declined to 3D-Newell). Verify the run's midpoint image
        // against the 3D curve's midpoint; a mismatch falls through
        // to the projection source, which walks the true samples.
        let image_honest = match curve {
            Some(c) => {
                use popkern_core::geometry::curve::CurveEval as _;
                use popkern_core::geometry::surface::SurfaceEval as _;
                let mid_uv = samples.0[samples.0.len() / 2];
                let s = surface.point_at(mid_uv.0, mid_uv.1);
                let cm = c.geometry.point_at(0.5 * (c.t_start + c.t_end));
                s.distance_to(&cm) < 1.0e-2
            }
            None => true,
        };
        if image_honest {
            run = Some((samples, source));
        }
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
        run = Some((
            (vec![ua, (ub, vb)], None, None),
            EmbeddingSource::EndpointProjection,
        ));
    }

    if run.is_none() && popkern_core::trace::trace_enabled("chart-dump") {
        let a = graph.nodes[graph.half_edges[raw].origin].position;
        let b = graph.nodes[graph.half_edges[twin].origin].position;
        eprintln!(
            "    [embed-none] he {raw} src {:?} ({:.3},{:.3},{:.3})->({:.3},{:.3},{:.3})",
            graph.half_edges[raw].source, a.x, a.y, a.z, b.x, b.y, b.z
        );
    }
    let ((mut uv, mut refined_first, mut refined_last), source) = run?;
    // Period-fold unfolding for pcurve-sourced runs — mechanism and
    // measured rationale in `chart_embedding_unfold`.
    if matches!(
        source,
        EmbeddingSource::OwnPcurve | EmbeddingSource::TwinPcurve
    ) {
        chart_embedding_unfold::unfold_pcurve_run(
            surface,
            &mut uv,
            &mut refined_first,
            &mut refined_last,
        );
    }
    // The resolved half-edge `h` runs opposite the registered `raw`
    // when the twin owned this face's pcurve/curve — reverse into
    // `raw`'s traversal direction (the fragment sampler's measured
    // zigzag lesson).
    if h != raw {
        uv.reverse();
        std::mem::swap(&mut refined_first, &mut refined_last);
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
        std::mem::swap(&mut refined_first, &mut refined_last);
    }
    // A closed record ring meets one node at both ends, so the node
    // test above cannot orient it; the curve holder's traversal
    // convention does (see `chart_embedding_closed`).
    if let Some(c) = curve
        && chart_embedding_closed::closed_run_against_curve(graph, raw, c, surface, &uv)
    {
        if popkern_core::trace::trace_enabled("boolean") {
            eprintln!(
                "  [embed] closed record he {raw} (resolved {h}, {source:?}, own-curve {curve_is_own}) on {face:?}: run reversed onto its curve direction"
            );
        }
        uv.reverse();
        std::mem::swap(&mut refined_first, &mut refined_last);
    }
    super::chart_poles::snap_pole_end_azimuths(surface, &mut uv, refined_first, refined_last);
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
        refined_first,
        refined_last,
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
) -> Option<RunWithRefinedEnds> {
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
    let mut run = Vec::with_capacity(N_SAMPLES + 1);
    for i in 0..=N_SAMPLES {
        let t = s0 + (s1 - s0) * i as f64 / N_SAMPLES as f64;
        let p = pc.point_at(t);
        run.push((p.x, p.y));
    }
    // Endpoint-refined direction samples at δ = span/1000, carried
    // BESIDE the run (see `ChartEmbedding::refined_first`).
    let pd0 = pc.point_at(s0 + (s1 - s0) * 1.0e-3);
    let pd1 = pc.point_at(s0 + (s1 - s0) * (1.0 - 1.0e-3));
    Some((run, Some((pd0.x, pd0.y)), Some((pd1.x, pd1.y))))
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
pub(crate) fn projection_run(
    c: &CurveOnEdge,
    curve_is_own: bool,
    surface: &SurfaceGeometry,
) -> Option<RunWithRefinedEnds> {
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
    // Continuation steps scale with the curve's winding (PK-670): a
    // 13.9-turn root edge at 16 stations stepped most of a turn per
    // projection and continued onto the adjacent turn (embedded
    // diagonally across its strip, uv (0,0.18)->(1,0.82)). Between
    // stations the projection walks `sub` steps of at most two per
    // knot span; `sub = 1` is the historical station-to-station walk.
    let sub = crate::merge::param_basin_search::curve_scan_count(&c.geometry, t0, t1, N_SAMPLES)
        .div_ceil(N_SAMPLES);
    let walk = |from: usize, to: usize, seed: (f64, f64)| -> Option<(f64, f64)> {
        let mut uv = seed;
        for k in 1..=sub {
            let t = if k == sub {
                ts[to]
            } else {
                ts[from] + (ts[to] - ts[from]) * k as f64 / sub as f64
            };
            uv = proj_near(&c.geometry.point_at(t), uv)?;
        }
        Some(uv)
    };
    for i in (mid + 1)..ts.len() {
        run[i] = walk(i - 1, i, run[i - 1])?;
    }
    for i in (0..mid).rev() {
        run[i] = walk(i + 1, i, run[i + 1])?;
    }
    // Endpoint-refined direction samples (same rationale as
    // `pcurve_run`), carried BESIDE the run (see
    // `ChartEmbedding::refined_first`).
    let td0 = t0 + (t1 - t0) * 1.0e-3;
    let td1 = t0 + (t1 - t0) * (1.0 - 1.0e-3);
    let d0 = proj_near(&c.geometry.point_at(td0), run[0])?;
    let d1 = proj_near(&c.geometry.point_at(td1), run[ts.len() - 1])?;
    Some((run, Some(d0), Some(d1)))
}
