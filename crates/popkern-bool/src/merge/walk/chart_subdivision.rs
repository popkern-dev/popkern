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

use super::chart_branches::{anchor_periodic, balance_and_duplicates, boundary_shoelace_sigma};

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
    /// A registered half-edge touches a chart-degenerate point (pole).
    PoleIncidence(usize),
    /// The face has no surface geometry to chart against.
    NoSurface,
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
}

impl std::fmt::Display for SubdivisionIssue {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::MissingEmbedding(h) => write!(f, "no chart embedding for he {h}"),
            Self::VPeriodicChart => write!(f, "v-periodic chart (unimplemented scope)"),
            Self::PoleIncidence(h) => write!(f, "he {h} touches a chart pole"),
            Self::NoSurface => write!(f, "face has no surface geometry"),
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
        }
    }
}

/// Why a half-edge was excluded from the subdivision (not an error —
/// the record bounds nothing; legacy behavior preserved, but observable).
#[derive(Debug, Clone, PartialEq)]
pub(crate) enum FilteredRecord {
    /// Chart extent below noise — a tangency/touch-point record.
    ZeroExtent(usize),
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
pub(crate) fn build_face_subdivision(
    graph: &IntersectionGraph,
    face: &FaceRef,
    he_indices: &[usize],
    surface: &SurfaceGeometry,
    orientation: bool,
) -> FaceChartSubdivision {
    let mut sub = FaceChartSubdivision {
        sigma: if orientation { 1.0 } else { -1.0 },
        ..Default::default()
    };

    let (u_period, v_period) = surface.parameter_periodicity();
    if v_period.is_some() {
        sub.report.issues.push(SubdivisionIssue::VPeriodicChart);
        return sub;
    }

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
                if pole_incident(surface, &emb) {
                    sub.report.issues.push(SubdivisionIssue::PoleIncidence(h));
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

    // Node endpoints per half-edge (start = origin(h), end = origin(twin)).
    let ends = |h: usize| -> (usize, usize) {
        (
            graph.half_edges[h].origin,
            graph.half_edges[graph.half_edges[h].twin].origin,
        )
    };

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
            cluster_single_branch(&mut sub, &ends);
        }
    }
    if !sub.report.issues.is_empty() {
        return sub;
    }

    balance_and_duplicates(&mut sub, u_period);
    sub
}

/// Finite-difference pole test at both run endpoints — only chart
/// types with rotational degeneracies are checked.
fn pole_incident(surface: &SurfaceGeometry, emb: &ChartEmbedding) -> bool {
    if !matches!(
        surface,
        SurfaceGeometry::Sphere { .. } | SurfaceGeometry::Cone { .. }
    ) {
        return false;
    }
    for &(u, v) in [emb.uv[0], emb.uv[emb.uv.len() - 1]].iter() {
        let p0 = surface.point_at(u, v);
        let p1 = surface.point_at(u + 1e-4, v);
        let scale = 1.0 + p0.distance_to(&popkern_core::geometry::Point3::new(0.0, 0.0, 0.0));
        if p0.distance_to(&p1) < 1e-11 * scale {
            return true;
        }
    }
    false
}

/// Non-periodic charts: every node is a single branch point.
fn cluster_single_branch(sub: &mut FaceChartSubdivision, ends: &dyn Fn(usize) -> (usize, usize)) {
    let hes: Vec<usize> = sub.embeddings.keys().copied().collect();
    for h in hes {
        let (s, e) = ends(h);
        let s_bp = BranchPoint { node: s, branch: 0 };
        let e_bp = BranchPoint { node: e, branch: 0 };
        sub.start_bp.insert(h, s_bp);
        sub.end_bp.insert(h, e_bp);
        sub.departures.entry(s_bp).or_default().push(h);
    }
}

/// True when chart point `p` lies on the interior of some OTHER
/// registered half-edge's run (coarse, diagnostics-only).
pub(super) fn node_on_some_run(sub: &FaceChartSubdivision, exclude: usize, p: (f64, f64)) -> bool {
    let (mut lo, mut hi) = ((f64::MAX, f64::MAX), (f64::MIN, f64::MIN));
    for emb in sub.embeddings.values() {
        for &(u, v) in &emb.uv {
            lo = (lo.0.min(u), lo.1.min(v));
            hi = (hi.0.max(u), hi.1.max(v));
        }
    }
    let eps = 1e-6 * (1.0 + (hi.0 - lo.0).abs() + (hi.1 - lo.1).abs());
    let ex_twin = exclude; // run of the arriving he itself is fine to skip only
    for (&h, emb) in &sub.embeddings {
        if h == ex_twin {
            continue;
        }
        for w in emb.uv.windows(2) {
            let (a, b) = (w[0], w[1]);
            let (dx, dy) = (b.0 - a.0, b.1 - a.1);
            let l2 = (dx * dx + dy * dy).max(1e-300);
            let t = (((p.0 - a.0) * dx + (p.1 - a.1) * dy) / l2).clamp(0.0, 1.0);
            let d = ((p.0 - a.0 - dx * t).powi(2) + (p.1 - a.1 - dy * t).powi(2)).sqrt();
            if d < eps && t > 1e-6 && t < 1.0 - 1e-6 {
                return true;
            }
        }
    }
    false
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
