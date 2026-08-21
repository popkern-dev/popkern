//! Per-face chart subdivision: branch points, rotation substrate, and
//! the I1–I5 validity checks (PK-104 stage 1 — see
//! docs/architecture/fragment-walk-chart-subdivision.md §3.4–3.5).
//!
//! Period-branch resolution is combinatorial: embeddings are anchored
//! to a common frame by whole-period translations (node-consensus BFS
//! across non-seam nodes; chart seam pairs assigned DETERMINISTICALLY
//! to the border where their travel direction keeps the interior on
//! the material side; only sub-period floating components keep a
//! bounded two-way search — see `chart_branches`). σ is
//! self-calibrated from the stored boundary's chart shoelace, with
//! the orientation flag as fallback. No UV epsilon welding is
//! introduced — branch points are keyed (graph node, cluster index).

use std::collections::BTreeMap;

use super::chart_branches::{
    anchor_periodic,
    clustering::{balance_and_duplicates, boundary_shoelace_sigma},
};

use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};

use super::chart_embedding::{ChartEmbedding, EmbeddingSource, build_embedding};
use crate::merge::graph::{FaceRef, IntersectionGraph};

/// A chart incidence point: a graph node on a specific period branch.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub(crate) struct BranchPoint {
    /// Graph node index.
    pub node: usize,
    /// Period-branch cluster index at that node (ascending u).
    pub branch: u32,
}

/// A structural defect that makes the face's chart subdivision
/// unwalkable; any issue routes the face to the legacy walker.
#[derive(Debug, Clone, PartialEq)]
pub(crate) enum SubdivisionIssue {
    /// No embedding source produced a UV run for this half-edge (I1).
    MissingEmbedding(usize),
    /// The chart is periodic in v — outside the implemented scope.
    VPeriodicChart,
    /// The anchored interior spans more than one period (I2).
    ChartOverwrap { span: f64, period: f64 },
    /// A node away from the chart seam produced >1 branch cluster (I2).
    BranchSplitOffSeam { node: usize },
    /// A seam half-edge could not be aligned to a chart border (I2).
    SeamBorderMisfit { he: usize },
    /// In-degree != out-degree at a branch point (I4).
    DegreeImbalance {
        node: usize,
        branch: u32,
        arrivals: usize,
        departures: usize,
    },
    /// An imbalanced arrival-only point lies on another registered
    /// half-edge's run interior — an uncompleted T-junction (I3).
    UnsplitTJunction { node: usize },
    /// Two departures at one branch point leave along the same chart
    /// segment (I5 — duplicate/overlap records).
    DuplicateChartSegment(usize, usize),
    /// The seam/component balance search exceeded its config budget.
    SeamSearchExhausted,
    /// A zero-extent or graze drop touched a pole node: the face
    /// carries degenerate pole-border topology (full-sphere BReps
    /// hold zero-length pole edges) whose semantics the pole-azimuth
    /// rotation does not model — recorded fallback.
    PoleAdjacentDegeneracy { node: usize },
}

impl std::fmt::Display for SubdivisionIssue {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::MissingEmbedding(h) => write!(f, "no chart embedding for he {h}"),
            Self::VPeriodicChart => write!(f, "v-periodic chart (unimplemented scope)"),
            Self::ChartOverwrap { span, period } => {
                write!(f, "anchored interior spans {span:.4} > period {period:.4}")
            }
            Self::BranchSplitOffSeam { node } => {
                write!(f, "branch split at non-seam node {node}")
            }
            Self::SeamBorderMisfit { he } => write!(f, "seam he {he} misfits chart borders"),
            Self::DegreeImbalance {
                node,
                branch,
                arrivals,
                departures,
            } => write!(
                f,
                "degree imbalance at node {node} branch {branch}: {arrivals} in / {departures} out"
            ),
            Self::UnsplitTJunction { node } => write!(f, "unsplit T-junction at node {node}"),
            Self::DuplicateChartSegment(a, b) => {
                write!(f, "duplicate chart segment: he {a} vs he {b}")
            }
            Self::SeamSearchExhausted => write!(f, "seam balance search exceeded budget"),
            Self::PoleAdjacentDegeneracy { node } => {
                write!(f, "degenerate pole-border topology at node {node}")
            }
        }
    }
}

/// Why a half-edge was excluded from the subdivision (not an error —
/// the record bounds nothing; legacy behavior preserved, but observable).
#[derive(Debug, Clone, PartialEq)]
pub(crate) enum FilteredRecord {
    /// Chart extent below noise — a tangency/touch-point record.
    ZeroExtent(usize),
    /// An intersection record whose chart segment coincides with a
    /// boundary edge or a lower-index record (I5): coplanar overlap
    /// booleans trace segments ON existing edges; the duplicate cover
    /// bounds nothing and unbalances every junction it touches.
    DuplicateCover { dropped: usize, kept: usize },
    /// One direction of an overlap-pair record standing in for a
    /// MISSING boundary piece (the coincident weld path leaves a gap
    /// in the face's boundary registration and the record spans it).
    /// The kept direction is the boundary substitute; the dropped one
    /// pointed material off-face and unbalanced the gap endpoints.
    GapSubstituteDirection { dropped: usize, kept: usize },
    /// One direction of a Shared pair lying ON the face's ORIGINAL
    /// outer boundary, resolved by flow continuation around the
    /// boundary cycle (see `chart_boundary_substitute`): the kept
    /// direction continues the registered boundary flow; the dropped
    /// one bounded the coplanar neighbour's region, not this face's.
    BoundarySubstituteDirection { dropped: usize, kept: usize },
    /// An interior graze: a record pair dead-ending at a node with no
    /// other registered half-edge on this face. Both sides of the
    /// trace bound the same region there — the legacy walk's
    /// interior-graze candidate filter refuses to enter these; the
    /// substrate drops the pair (design doc §5).
    InteriorGraze(usize),
}

/// Structured result of the substrate pass for one face.
#[derive(Debug, Default)]
pub(crate) struct ChartSubdivisionReport {
    /// Embedding source counts: [own pcurve, twin pcurve, projection,
    /// endpoint projection].
    pub source_counts: [u32; 4],
    /// Blocking defects (empty ⇔ the face is permutation-walkable).
    pub issues: Vec<SubdivisionIssue>,
    /// Records excluded as non-subdivision data.
    pub filtered: Vec<FilteredRecord>,
    /// True when the bounded seam/component balance search ran.
    pub seam_search_used: bool,
}

/// The rotation substrate for one face: embeddings, branch points,
/// and departures — everything `next()` needs, and nothing mutable.
#[derive(Debug, Default)]
pub(crate) struct FaceChartSubdivision {
    /// Chart embedding per registered half-edge (post-anchoring).
    pub embeddings: BTreeMap<usize, ChartEmbedding>,
    /// Branch point each half-edge departs from.
    pub start_bp: BTreeMap<usize, BranchPoint>,
    /// Branch point each half-edge arrives at.
    pub end_bp: BTreeMap<usize, BranchPoint>,
    /// Half-edges departing each branch point (ascending he index).
    pub departures: BTreeMap<BranchPoint, Vec<usize>>,
    /// +1 when chart CCW = CCW about the outward normal
    /// (`SurfaceOnFace.orientation == true`), −1 otherwise.
    pub sigma: f64,
    /// Chart u-period (rotational charts), for pole-azimuth rotation.
    pub u_period: Option<f64>,
    /// Pole nodes (chart-degenerate points: sphere poles, cone
    /// apices) → rotation direction for the azimuth successor:
    /// +1 ⇒ smallest positive (u_arrive − u_depart) mod period wins,
    /// −1 ⇒ the mirrored order (σ and pole side folded in).
    pub poles: BTreeMap<usize, f64>,
    /// Substrate diagnostics.
    pub report: ChartSubdivisionReport,
}

impl FaceChartSubdivision {
    /// True when the face can be walked as a pure permutation.
    pub fn complete(&self) -> bool {
        self.report.issues.is_empty()
    }
}

/// Build the chart subdivision for one face's registered half-edges.
///
/// `overlap_pairs`: coplanar face-overlap pairs from coincidence
/// detection; records authored by exactly these pairs get the
/// partial-cover dedup (their traces along existing boundary edges
/// bound nothing).
///
/// `original_brep`: the operand BRep this face belongs to — the
/// boundary-substitute resolution reads the face's ORIGINAL outer
/// loop from it (the graph registration can carry weld gaps).
#[allow(clippy::too_many_arguments)]
pub(crate) fn build_face_subdivision(
    graph: &IntersectionGraph,
    face: &FaceRef,
    he_indices: &[usize],
    surface: &SurfaceGeometry,
    orientation: bool,
    overlap_pairs: &[(
        popkern_core::topology::FaceHandle,
        popkern_core::topology::FaceHandle,
    )],
    original_brep: &popkern_core::topology::BRep,
    tolerance: &popkern_core::tolerance::Tolerance,
) -> FaceChartSubdivision {
    let sub = build_face_subdivision_inner(
        graph,
        face,
        he_indices,
        tolerance,
        surface,
        orientation,
        overlap_pairs,
        original_brep,
    );
    if std::env::var_os("PK_CHART_DUMP").is_some() {
        eprintln!(
            "  [chart] {face:?} sigma={} poles={:?} issues={:?} filtered={:?}",
            sub.sigma, sub.poles, sub.report.issues, sub.report.filtered
        );
        for (&h, emb) in &sub.embeddings {
            let (s0, e0) = (emb.uv[0], emb.uv[emb.uv.len() - 1]);
            let (na, nb) = (
                graph.half_edges[h].origin,
                graph.half_edges[graph.half_edges[h].twin].origin,
            );
            let (pa, pb) = (graph.nodes[na].position, graph.nodes[nb].position);
            eprintln!(
                "    [chart] he {h} tw {} {:?} src={:?} ({:.3},{:.3})->({:.3},{:.3}) n{na}({:.3},{:.3},{:.3})->n{nb}({:.3},{:.3},{:.3}) sbp={:?} ebp={:?} seam={}",
                graph.half_edges[h].twin,
                emb.source,
                graph.half_edges[h].source,
                s0.0,
                s0.1,
                e0.0,
                e0.1,
                pa.x,
                pa.y,
                pa.z,
                pb.x,
                pb.y,
                pb.z,
                sub.start_bp.get(&h),
                sub.end_bp.get(&h),
                emb.seam_like
            );
        }
    }
    sub
}

#[allow(clippy::too_many_arguments)]
fn build_face_subdivision_inner(
    graph: &IntersectionGraph,
    face: &FaceRef,
    he_indices: &[usize],
    tolerance: &popkern_core::tolerance::Tolerance,
    surface: &SurfaceGeometry,
    orientation: bool,
    overlap_pairs: &[(
        popkern_core::topology::FaceHandle,
        popkern_core::topology::FaceHandle,
    )],
    original_brep: &popkern_core::topology::BRep,
) -> FaceChartSubdivision {
    let mut sub = FaceChartSubdivision {
        sigma: if orientation { 1.0 } else { -1.0 },
        ..Default::default()
    };

    let (u_period, v_period) = surface.parameter_periodicity();
    // Doubly-periodic (torus) charts: a SUB-PERIOD face (quarter
    // patch — the popkern torus primitive's native granularity) needs
    // no seam machinery at all. Its chart cell is computed from the
    // stored boundary ring, every embedding is folded into that cell
    // by whole-period shifts, and the rest of the pipeline runs the
    // aperiodic path. Faces that WRAP a period (no resolvable cell)
    // remain out of scope and take the recorded legacy fallback.
    let cell_fold = if v_period.is_some() {
        match crate::clip::chart_cell::face_chart_cell(original_brep, face.face, tolerance) {
            Some(cell) => Some((cell, u_period, v_period)),
            None => {
                sub.report.issues.push(SubdivisionIssue::VPeriodicChart);
                return sub;
            }
        }
    } else {
        None
    };
    let u_period = if cell_fold.is_some() { None } else { u_period };

    // I1: embed every registered half-edge; drop zero-extent records.
    for &h in he_indices {
        match build_embedding(graph, face, h, surface) {
            Some(emb) => {
                let idx = match emb.source {
                    EmbeddingSource::OwnPcurve => 0,
                    EmbeddingSource::TwinPcurve => 1,
                    EmbeddingSource::Projection => 2,
                    EmbeddingSource::EndpointProjection => 3,
                };
                sub.report.source_counts[idx] += 1;
                if emb.start_dir().is_none() {
                    sub.report.filtered.push(FilteredRecord::ZeroExtent(h));
                    continue;
                }
                sub.embeddings.insert(h, emb);
            }
            None => sub
                .report
                .issues
                .push(SubdivisionIssue::MissingEmbedding(h)),
        }
    }
    if !sub.report.issues.is_empty() {
        return sub;
    }
    if let Some((cell, pu, pv)) = &cell_fold {
        // Whole-period fold per RUN (a sub-period face's runs never
        // straddle a seam, so one shift per run is exact): anchor the
        // run's midpoint into the cell on both axes.
        let shift = |x: f64, range: (f64, f64), p: &Option<f64>| match p {
            Some(p) if p.is_finite() && *p > 0.0 => {
                ((0.5 * (range.0 + range.1) - x) / p).round() * p
            }
            _ => 0.0,
        };
        for emb in sub.embeddings.values_mut() {
            let mid = emb.uv[emb.uv.len() / 2];
            let du = shift(mid.0, cell.u_range, pu);
            let dv = shift(mid.1, cell.v_range, pv);
            if du != 0.0 || dv != 0.0 {
                for q in emb.uv.iter_mut() {
                    q.0 += du;
                    q.1 += dv;
                }
            }
        }
    }
    sub.u_period = u_period;

    // POLE NODES (design doc §5): a chart-degenerate endpoint (sphere
    // pole, cone apex) is ONE geometric point whose UV image is the
    // whole u-line at the pole v — every incident half-edge lands
    // there at ITS OWN meridian azimuth. The subdivision stays valid:
    // the pole is a single branch point, and the rotational successor
    // there orders departures by azimuth (the endpoint u), not by
    // local UV direction. The rotation sign folds σ with the pole
    // side (rotation sense reverses between the two poles).
    if u_period.is_some()
        && matches!(
            surface,
            SurfaceGeometry::Sphere { .. } | SurfaceGeometry::Cone { .. }
        )
    {
        let (mut v_lo, mut v_hi) = (f64::MAX, f64::MIN);
        for emb in sub.embeddings.values() {
            for &(_, v) in &emb.uv {
                v_lo = v_lo.min(v);
                v_hi = v_hi.max(v);
            }
        }
        let v_mid = 0.5 * (v_lo + v_hi);
        let hes: Vec<usize> = sub.embeddings.keys().copied().collect();
        for h in hes {
            let emb = &sub.embeddings[&h];
            let (n_start, n_end) = (
                graph.half_edges[h].origin,
                graph.half_edges[graph.half_edges[h].twin].origin,
            );
            for (node, &(u, v)) in [(n_start, &emb.uv[0]), (n_end, &emb.uv[emb.uv.len() - 1])] {
                if sub.poles.contains_key(&node) {
                    continue;
                }
                let p0 = surface.point_at(u, v);
                let p1 = surface.point_at(u + 1e-4, v);
                let scale =
                    1.0 + p0.distance_to(&popkern_core::geometry::Point3::new(0.0, 0.0, 0.0));
                if p0.distance_to(&p1) < 1e-11 * scale {
                    // Pole side decides the rotation sense; σ is not
                    // yet calibrated here, so store the side and fold
                    // σ in at the end of the build.
                    let side = if v >= v_mid { 1.0 } else { -1.0 };
                    sub.poles.insert(node, side);
                }
            }
        }
    }

    // I5 — duplicate-cover dedup: coplanar overlap records trace
    // chart segments that COINCIDE with boundary edges (or with each
    // other). The duplicate cover bounds no region and unbalances
    // every junction it touches (measured on the flush-gusset weld:
    // ±1 in/out at each overlap endpoint, whole faces falling back).
    // Boundary edges outrank records; between records, the lower
    // index wins. Comparison is geometric (endpoint + midpoint UVs,
    // either direction) so it is authoring-order independent.
    {
        let hes: Vec<usize> = sub.embeddings.keys().copied().collect();
        let mut scale = 0.0_f64;
        for emb in sub.embeddings.values() {
            let (s, e) = (emb.uv[0], emb.uv[emb.uv.len() - 1]);
            scale = scale.max((s.0 - e.0).abs() + (s.1 - e.1).abs());
        }
        let tol = 1e-7 * (1.0 + scale);
        let near = |a: (f64, f64), b: (f64, f64)| -> bool {
            (a.0 - b.0).abs() < tol && (a.1 - b.1).abs() < tol
        };
        let seg = |h: usize| -> ((f64, f64), (f64, f64), (f64, f64)) {
            let uv = &sub.embeddings[&h].uv;
            (uv[0], uv[uv.len() / 2], uv[uv.len() - 1])
        };
        let is_record = |h: usize| -> bool {
            matches!(
                graph.half_edges[h].source,
                crate::merge::graph::EdgeSource::Shared { .. }
            )
        };
        let mut drop: Vec<(usize, usize)> = Vec::new();
        for (i, &a) in hes.iter().enumerate() {
            for &b in hes.iter().skip(i + 1) {
                if graph.half_edges[a].twin == b {
                    continue; // twins are the two sides of one edge
                }
                let (sa, ma, ea) = seg(a);
                let (sb, mb, eb) = seg(b);
                let same = near(sa, sb) && near(ea, eb) && near(ma, mb);
                let reversed = near(sa, eb) && near(ea, sb) && near(ma, mb);
                if !(same || reversed) {
                    continue;
                }
                match (is_record(a), is_record(b)) {
                    (true, false) => drop.push((a, b)),
                    (false, true) => drop.push((b, a)),
                    (true, true) => drop.push((b, a)),
                    (false, false) => {} // two boundary edges: not ours to judge
                }
            }
        }
        for (d, k) in drop {
            if sub.embeddings.remove(&d).is_some() {
                sub.report.filtered.push(FilteredRecord::DuplicateCover {
                    dropped: d,
                    kept: k,
                });
            }
        }
        // PARTIAL covers: overlap records are windowed at different
        // nodes than the boundary edges they trace (a weld trace can
        // span half of one boundary edge and all of the next), so
        // segment-identity misses them. SCOPE: only records authored
        // by a known coplanar overlap pair — unconditional
        // sample-on-polyline dropping was measured to also remove
        // records that legitimately bound flush regions (two bool
        // fixtures regressed). Within the overlap pair the semantics
        // are exact: the record traces IN the shared plane, so where
        // it lies on this face's boundary the boundary already
        // bounds, and the record covers nothing new.
        let from_overlap_pair = |h: usize| -> bool {
            match &graph.half_edges[h].source {
                crate::merge::graph::EdgeSource::Shared { face_a, face_b } => overlap_pairs
                    .iter()
                    .any(|&(a, b)| *face_a == a && *face_b == b),
                _ => false,
            }
        };
        let boundary_polys: Vec<Vec<(f64, f64)>> = sub
            .embeddings
            .iter()
            .filter(|(h, _)| !is_record(**h))
            .map(|(_, e)| e.uv.clone())
            .collect();
        let on_boundary = |p: (f64, f64)| -> bool {
            for poly in &boundary_polys {
                for w in poly.windows(2) {
                    let (a, b) = (w[0], w[1]);
                    let (dx, dy) = (b.0 - a.0, b.1 - a.1);
                    let l2 = (dx * dx + dy * dy).max(1e-300);
                    let t = (((p.0 - a.0) * dx + (p.1 - a.1) * dy) / l2).clamp(0.0, 1.0);
                    let d2 = (p.0 - a.0 - dx * t).powi(2) + (p.1 - a.1 - dy * t).powi(2);
                    if d2 < tol * tol {
                        return true;
                    }
                }
            }
            false
        };
        let overlap_hes: Vec<usize> = sub
            .embeddings
            .keys()
            .copied()
            .filter(|&h| from_overlap_pair(h))
            .collect();
        for h in overlap_hes {
            let covered = sub.embeddings[&h].uv.iter().all(|&p| on_boundary(p));
            if covered && sub.embeddings.remove(&h).is_some() {
                sub.report.filtered.push(FilteredRecord::DuplicateCover {
                    dropped: h,
                    kept: usize::MAX,
                });
            }
        }
    }

    // Node endpoints per half-edge (start = origin(h), end = origin(twin)).
    let ends = |h: usize| -> (usize, usize) {
        (
            graph.half_edges[h].origin,
            graph.half_edges[graph.half_edges[h].twin].origin,
        )
    };

    // Zero-extent records whose pcurves wiggle enough to fake a
    // direction slip past the start_dir test — filter by UV bbox
    // extent when both graph endpoints are the same node (measured:
    // a touch-point pair on the flush-gusset arm edge carried 17
    // near-identical samples and survived to unbalance its node).
    {
        let mut scale = 0.0_f64;
        for emb in sub.embeddings.values() {
            let (s0, e0) = (emb.uv[0], emb.uv[emb.uv.len() - 1]);
            scale = scale.max((s0.0 - e0.0).abs() + (s0.1 - e0.1).abs());
        }
        let tol = 1e-9 * (1.0 + scale);
        let degenerate: Vec<usize> = sub
            .embeddings
            .iter()
            .filter(|&(&h, emb)| {
                let (a, b) = ends(h);
                if a != b {
                    return false;
                }
                let (mut lo_u, mut hi_u, mut lo_v, mut hi_v) =
                    (f64::MAX, f64::MIN, f64::MAX, f64::MIN);
                for &(u, v) in &emb.uv {
                    lo_u = lo_u.min(u);
                    hi_u = hi_u.max(u);
                    lo_v = lo_v.min(v);
                    hi_v = hi_v.max(v);
                }
                (hi_u - lo_u) + (hi_v - lo_v) < tol
            })
            .map(|(&h, _)| h)
            .collect();
        for h in degenerate {
            sub.embeddings.remove(&h);
            sub.report.filtered.push(FilteredRecord::ZeroExtent(h));
        }
    }

    // BOUNDARY-SUBSTITUTE resolution (coincident booleans, non-
    // periodic charts): Shared pairs lying ON the face's ORIGINAL
    // outer boundary replace welded-away boundary pieces; exactly one
    // direction bounds this face, decided by agreement with the
    // original loop traversal. Runs BEFORE the balance search below,
    // which cannot decide this class — the full story and measured
    // failures live in `chart_boundary_substitute`.
    if !overlap_pairs.is_empty() && u_period.is_none() {
        super::chart_boundary_substitute::resolve_boundary_substitutes(
            graph,
            face,
            &mut sub,
            surface,
            original_brep,
        );
    }

    // GAP-SUBSTITUTE completion (non-periodic charts): the coincident
    // weld path can leave a face's boundary registration with a GAP
    // (a split middle piece missing) that an overlap-pair record
    // spans — legacy walks the record as the boundary substitute.
    // With BOTH record directions registered the gap endpoints are
    // ±1 imbalanced; exactly one drop choice restores balance at
    // both, and it is the direction agreeing with the adjacent
    // boundary flow. Greedy per pair, deterministic, recorded.
    if u_period.is_none() {
        let degree = |sub: &FaceChartSubdivision| -> std::collections::HashMap<usize, i64> {
            let mut d: std::collections::HashMap<usize, i64> = std::collections::HashMap::new();
            for &h in sub.embeddings.keys() {
                let (a, b) = ends(h);
                *d.entry(a).or_default() += 1; // departure
                *d.entry(b).or_default() -= 1; // arrival
            }
            d
        };
        let imbalance_total = |d: &std::collections::HashMap<usize, i64>| -> i64 {
            d.values().map(|v| v.abs()).sum()
        };
        let deg = degree(&sub);
        if imbalance_total(&deg) > 0 {
            let mut pairs: Vec<(usize, usize)> = Vec::new();
            for &h in sub.embeddings.keys() {
                let tw = graph.half_edges[h].twin;
                if h < tw
                    && sub.embeddings.contains_key(&tw)
                    && matches!(
                        graph.half_edges[h].source,
                        crate::merge::graph::EdgeSource::Shared { .. }
                    )
                {
                    pairs.push((h, tw));
                }
            }
            // Candidate pairs: those touching an imbalanced node. A
            // single drop can be NET-ZERO (it moves the imbalance to
            // its far endpoint) while a combination across pairs
            // resolves everything — greedy per-pair stalls there
            // (measured on the mirrored flush-arm faces). Bounded
            // exhaustive search: 3 states per pair (keep / drop
            // forward / drop twin), ≤4 pairs → ≤81 configs; first
            // config in deterministic order reaching TOTAL BALANCE
            // wins. No full balance → no drops (fallback as before).
            let candidates: Vec<(usize, usize)> = pairs
                .iter()
                .copied()
                .filter(|&(h, _)| {
                    let (a, b) = ends(h);
                    deg.get(&a).copied().unwrap_or(0) != 0 || deg.get(&b).copied().unwrap_or(0) != 0
                })
                .collect();
            if !candidates.is_empty() && candidates.len() <= 4 {
                let n_cfg = 3usize.pow(candidates.len() as u32);
                let mut winning: Option<Vec<(usize, usize)>> = None;
                'cfg: for cfg in 0..n_cfg {
                    let mut d = deg.clone();
                    let mut drops: Vec<(usize, usize)> = Vec::new();
                    let mut c = cfg;
                    for &(h, tw) in &candidates {
                        let choice = c % 3;
                        c /= 3;
                        let dropped = match choice {
                            0 => None,
                            1 => Some((h, tw)),
                            _ => Some((tw, h)),
                        };
                        if let Some((dh, kh)) = dropped {
                            let (da, db) = ends(dh);
                            *d.entry(da).or_default() -= 1;
                            *d.entry(db).or_default() += 1;
                            drops.push((dh, kh));
                        }
                    }
                    if imbalance_total(&d) == 0 {
                        winning = Some(drops);
                        break 'cfg;
                    }
                }
                if let Some(drops) = winning {
                    for (dh, kh) in drops {
                        if sub.embeddings.remove(&dh).is_some() {
                            sub.report
                                .filtered
                                .push(FilteredRecord::GapSubstituteDirection {
                                    dropped: dh,
                                    kept: kh,
                                });
                        }
                    }
                }
            }
        }
    }

    // INTERIOR-GRAZE filter (design doc §5): a record pair whose
    // endpoint node has NO other registered half-edge dead-ends
    // strictly inside the face — both sides of the trace bound the
    // same region there, and walking it weaves an out-and-back spur
    // into a fragment loop (measured on the flush-gusset wall: the
    // gusset's sloped-face trace ends mid-wall at a degree-1 node,
    // and the chart outer loop carried the spur as its two extra
    // half-edges vs the legacy walk, which filters grazes at
    // candidate time). Iterated: removing one graze can expose a
    // chain.
    for _ in 0..8 {
        let mut incident: std::collections::HashMap<usize, Vec<usize>> =
            std::collections::HashMap::new();
        for &h in sub.embeddings.keys() {
            let (a, b) = ends(h);
            incident.entry(a).or_default().push(h);
            incident.entry(b).or_default().push(h);
        }
        let mut grazes: Vec<(usize, usize)> = Vec::new();
        for &h in sub.embeddings.keys() {
            let tw = graph.half_edges[h].twin;
            if h > tw || !sub.embeddings.contains_key(&tw) {
                continue;
            }
            if !matches!(
                graph.half_edges[h].source,
                crate::merge::graph::EdgeSource::Shared { .. }
            ) {
                continue;
            }
            let (a, b) = ends(h);
            if a == b {
                continue; // closed rings are legitimate cuts
            }
            let dead_end = |n: usize| -> bool {
                incident
                    .get(&n)
                    .map(|v| v.iter().all(|&x| x == h || x == tw))
                    .unwrap_or(true)
            };
            if dead_end(a) || dead_end(b) {
                grazes.push((h, tw));
            }
        }
        if grazes.is_empty() {
            break;
        }
        for (h, tw) in grazes {
            sub.embeddings.remove(&h);
            sub.embeddings.remove(&tw);
            sub.report.filtered.push(FilteredRecord::InteriorGraze(h));
        }
    }

    // Pole-adjacency guard: the azimuth rotation models poles as
    // clean fans of meridian endpoints. When a RECORD drop or a graze
    // cascade touched a pole node, the face lost real connectivity
    // there (measured: a full sphere's garbage hull-arc records
    // dead-ended mid-face, graze-filtered into a cascade, and the
    // single surviving orbit broke Euler downstream) — those faces
    // keep the legacy walk, recorded. A ZERO-EXTENT drop of the
    // operand's own degenerate pole BORDER is exempt: a pole is a
    // collapsed border by construction, and dropping its zero-length
    // edge is precisely how the fan model expects to see it.
    if !sub.poles.is_empty() {
        let dropped: Vec<(usize, bool)> = sub
            .report
            .filtered
            .iter()
            .map(|f| match f {
                FilteredRecord::ZeroExtent(h) => (*h, true),
                FilteredRecord::InteriorGraze(h) => (*h, false),
                FilteredRecord::DuplicateCover { dropped, .. } => (*dropped, false),
                FilteredRecord::GapSubstituteDirection { dropped, .. } => (*dropped, false),
                FilteredRecord::BoundarySubstituteDirection { dropped, .. } => (*dropped, false),
            })
            .collect();
        for (h, zero_extent) in dropped {
            let benign_pole_border = zero_extent
                && !matches!(
                    graph.half_edges[h].source,
                    crate::merge::graph::EdgeSource::Shared { .. }
                );
            if benign_pole_border {
                continue;
            }
            let (a, b) = ends(h);
            for n in [a, b] {
                if sub.poles.contains_key(&n) {
                    sub.report
                        .issues
                        .push(SubdivisionIssue::PoleAdjacentDegeneracy { node: n });
                }
            }
        }
        sub.report.issues.dedup();
        if !sub.report.issues.is_empty() {
            return sub;
        }
    }

    match u_period {
        Some(period) => anchor_periodic(graph, &mut sub, period, &ends),
        None => {
            // σ self-calibration: the operand face's stored boundary
            // loops wind CCW about outward (winding census E1), so
            // their chart shoelace sign IS the chart handedness. The
            // stored orientation flag is only the fallback — measured
            // dishonest on extrude bottom caps (flag true, chart
            // left-handed vs outward), which mirror-flipped every
            // junction pick on the breached-bore cap face.
            if let Some(s) = boundary_shoelace_sigma(graph, &sub) {
                sub.sigma = s;
            }
            super::chart_branches::clustering::cluster_single_branch(&mut sub, &ends);
        }
    }
    if !sub.report.issues.is_empty() {
        return sub;
    }

    balance_and_duplicates(&mut sub, u_period);
    // Fold the calibrated σ into the pole rotation directions (σ can
    // change during the periodic anchoring's shoelace calibration).
    let sigma = sub.sigma;
    for dir in sub.poles.values_mut() {
        *dir *= sigma;
    }
    sub
}

/// Resolve a face's surface geometry and orientation flag.
pub(crate) fn face_surface<'a>(
    face: &FaceRef,
    brep_a: &'a popkern_core::topology::BRep,
    brep_b: &'a popkern_core::topology::BRep,
) -> Option<(&'a SurfaceGeometry, bool)> {
    let brep = match face.operand {
        crate::merge::graph::Operand::A => brep_a,
        crate::merge::graph::Operand::B => brep_b,
    };
    let s = brep.faces.get(face.face)?.surface.as_ref()?;
    Some((&s.geometry, s.orientation))
}
