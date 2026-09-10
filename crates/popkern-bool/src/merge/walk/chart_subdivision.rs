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
/// unwalkable; any issue routes the face to the 3D-Newell walker.
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
/// the record bounds nothing; earlier behavior preserved, but observable).
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
    /// trace bound the same region there — the 3D-Newell walk's
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
    if popkern_core::trace::trace_enabled("chart-dump") {
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

    let (u_analytic, v_period) = surface.parameter_periodicity();
    // The walk runs entirely in chart space: `u_period` only drives
    // seam-branch identification (u=0 with u=period), never
    // evaluation beyond the chart — so the GEOMETRIC closure of a
    // NURBS chart counts exactly like an analytic period here
    // (PK-443: a lofted closed profile's wall is one u-closed NURBS
    // face; without the period the wall's chart-wrapping intersection
    // loops walked as interior ISLANDS and the face never split into
    // bands). v-closed NURBS charts stay out of scope — the v-period
    // cell-fold path is torus-specific.
    let (u_chart, v_chart) = crate::clip::classification_chart_periods(surface);
    let u_period = u_analytic.or(u_chart);
    // A chart closed geometrically in v ONLY (a NURBS surface of
    // revolution stores the angle in v) has no seam machinery on
    // this path: its azimuth ring pairs onto itself at the junction
    // node, the loop partition strands the seam pair as a zero-area
    // slit, and the roles coin-flip on zero shoelaces (measured,
    // PK-455: the SECOND same-face sphere-dimple subtract rebuilt the
    // surviving cap as rim-ring + seam-slit and the trim clip shred
    // the mesh). Decline to the recorded 3D-Newell fallback instead.
    if u_period.is_none() && v_period.is_none() && v_chart.is_some() {
        sub.report.issues.push(SubdivisionIssue::VPeriodicChart);
        return sub;
    }
    // Doubly-periodic (torus) charts: a SUB-PERIOD face (quarter
    // patch — the popkern torus primitive's native granularity) needs
    // no seam machinery at all. Its chart cell is computed from the
    // stored boundary ring, every embedding is folded into that cell
    // by whole-period shifts, and the rest of the pipeline runs the
    // aperiodic path. Faces that WRAP a period (no resolvable cell)
    // remain out of scope and take the recorded 3D-Newell fallback.
    let cell_fold = if v_period.is_some() {
        match crate::clip::chart_cell::face_chart_cell(original_brep, face.face, tolerance) {
            Some(cell) => Some((cell, u_period, v_period)),
            // No resolvable cell — the face WRAPS at least one
            // period. A face that wraps u but is v-BOUNDED (a
            // revolved groove arc: full-ring torus band) is, for the
            // walk, exactly a u-closed wall: the v-period of the
            // SURFACE never enters its chart, and the tuned u-seam
            // machinery owns it (measured, PK-459: a radial drill
            // through a two-arc revolved tube's band declined here
            // and the 3D-Newell fallback elected a drill rim as the
            // face's outer — the wall classified Inside and the
            // serve kept two once-used rim circles). Faces wrapping
            // BOTH periods (a closed profile revolved full circle as
            // one face) stay out of scope and decline.
            None if face_v_extent_bounded(original_brep, face.face, v_period) => None,
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
        super::chart_cell_fold::fold_embeddings_into_cell(
            graph, surface, &mut sub, cell, pu, pv, face,
        );
    }
    sub.u_period = u_period;

    // POLE NODES (design doc §5) — detection extracted to
    // `chart_poles.rs`; the pole is a single branch point and the
    // rotational successor orders departures by azimuth.
    super::chart_poles::detect_pole_nodes(graph, surface, &mut sub, u_period);

    // I5 — duplicate-cover dedup, extracted to
    // `chart_duplicate_cover.rs` per the file-size rule.
    super::chart_duplicate_cover::dedup_duplicate_covers(graph, &mut sub, overlap_pairs);

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
    // spans — the 3D-Newell walk takes the record as the boundary substitute.
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

    super::chart_graze_filter::filter_hairpins_and_grazes(graph, &mut sub, tolerance, &ends);

    // Pole-adjacency guard: the azimuth rotation models poles as
    // clean fans of meridian endpoints. When a RECORD drop or a graze
    // cascade touched a pole node, the face lost real connectivity
    // there (measured: a full sphere's garbage hull-arc records
    // dead-ended mid-face, graze-filtered into a cascade, and the
    // single surviving orbit broke Euler downstream) — those faces
    // keep the 3D-Newell walk, recorded. A ZERO-EXTENT drop of the
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
        Some(period) => {
            if popkern_core::trace::trace_enabled("branch-pairing") {
                eprintln!("  [bp-face] {face:?}");
            }
            anchor_periodic(graph, &mut sub, period, &ends)
        }
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

/// True when the face's boundary pcurves span strictly less than the
/// surface's v-period: the face is a v-bounded band of a
/// doubly-periodic surface, and its chart never touches the v-seam.
fn face_v_extent_bounded(
    brep: &popkern_core::topology::BRep,
    face: popkern_core::topology::FaceHandle,
    v_period: Option<f64>,
) -> bool {
    let Some(period) = v_period else {
        return false;
    };
    let Some(boundary) = crate::clip::face_uv_boundary(brep, face) else {
        return false;
    };
    let mut lo = f64::INFINITY;
    let mut hi = f64::NEG_INFINITY;
    for contour in
        std::iter::once(boundary).chain(crate::clip::face_uv_inner_boundaries(brep, face))
    {
        for (_he, pcurve) in &contour {
            let (t0, t1) = pcurve.parameter_range();
            for i in 0..=16 {
                let q = pcurve.point_at(t0 + (t1 - t0) * i as f64 / 16.0);
                lo = lo.min(q.y);
                hi = hi.max(q.y);
            }
        }
    }
    lo.is_finite() && hi - lo < period * (1.0 - 1.0e-6)
}
