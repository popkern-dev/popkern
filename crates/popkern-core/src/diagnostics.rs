//! Comprehensive BRep diagnostics — a one-shot health card.
//!
//! Where `BRep::validate` answers a binary "is this BRep structurally
//! and geometrically sound?", `Diagnostics::run` answers the harder
//! question: "what's *interesting* about this BRep?" It surfaces
//! near-degenerate geometry, suspicious sliver faces, near-duplicate
//! vertices, T-junctions, NURBS irregularities, pcurve UV-escape, and
//! a surface-kind histogram — the kind of summary you want to print
//! after a build to know whether the output is worth exporting at all.
//!
//! Coverage: topology + geometry summaries; naked, degenerate,
//! sliver, coincident defect detectors.
//!
//! Extended coverage: T-junctions, NURBS irregularities, pcurve
//! UV-escape detectors. NURBS surface-grid validation deferred to
//! later revisions.
//!
//! Detectors run in O(F + E + V·E) — the V·E T-junction check
//! dominates and is fine for showcase-fixture sizes (tens of edges,
//! tens of vertices).

use std::collections::BTreeMap;
use std::fmt;

use crate::geometry::Point3;
use crate::topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle, VertexHandle};

#[path = "diagnostics_detectors.rs"]
mod detectors;

impl Default for GeometryStats {
    fn default() -> Self {
        Self {
            edge_length_min: 0.0,
            edge_length_median: 0.0,
            edge_length_max: 0.0,
            face_area_min: 0.0,
            face_area_median: 0.0,
            face_area_max: 0.0,
            bbox_min: Point3::new(0.0, 0.0, 0.0),
            bbox_max: Point3::new(0.0, 0.0, 0.0),
            bbox_diagonal: 0.0,
        }
    }
}

/// Top-level diagnostic report. Build with [`Diagnostics::run`].
#[derive(Debug, Clone, Default)]
pub struct Diagnostics {
    /// Counts and surface-kind histogram.
    pub topology: TopologyStats,
    /// Edge-length and face-area distributions.
    pub geometry: GeometryStats,
    /// Per-category defect lists (each entry has handle + location).
    pub defects: Defects,
}

/// Topology counts and surface-kind histogram. The Euler characteristic
/// `V - E + F` should be 2 for a closed genus-0 solid; deviations are
/// surfaced for hand-investigation.
#[derive(Debug, Clone, Default)]
pub struct TopologyStats {
    /// Total face count.
    pub faces: usize,
    /// Number of faces that have a surface assigned.
    pub faces_with_surface: usize,
    /// Total edge count.
    pub edges: usize,
    /// Total vertex count.
    pub vertices: usize,
    /// Total half-edge count.
    pub half_edges: usize,
    /// `V - E + F` (genus-0 closed solid → 2).
    pub euler_chi: i64,
    /// Histogram of surface kinds, keyed by short label.
    pub surface_kinds: BTreeMap<&'static str, usize>,
}

/// Edge-length and face-area distribution statistics. All values are
/// in unit-of-length (typically millimetres after the recent STEP
/// unit fix). Bounding-box info is the corner-and-extent of all
/// vertex positions — a sanity-check against the wrong-scale-on-import
/// failure mode.
#[derive(Debug, Clone)]
pub struct GeometryStats {
    /// Smallest edge length observed.
    pub edge_length_min: f64,
    /// Median edge length.
    pub edge_length_median: f64,
    /// Largest edge length.
    pub edge_length_max: f64,
    /// Smallest face area observed.
    pub face_area_min: f64,
    /// Median face area.
    pub face_area_median: f64,
    /// Largest face area.
    pub face_area_max: f64,
    /// Min corner of the vertex-position AABB.
    pub bbox_min: Point3,
    /// Max corner of the vertex-position AABB.
    pub bbox_max: Point3,
    /// AABB diagonal length.
    pub bbox_diagonal: f64,
}

/// Aggregate of all defect categories. Empty vectors mean no defects
/// in that category.
#[derive(Debug, Clone, Default)]
pub struct Defects {
    /// Edges where one or both adjacent faces have no surface.
    pub naked_edges: Vec<NakedEdge>,
    /// Edges shorter than `tolerance.linear()`.
    pub degenerate_edges: Vec<DegenerateEdge>,
    /// Faces with very small area or extreme aspect ratio.
    pub sliver_faces: Vec<SliverFace>,
    /// Pairs of distinct vertex handles within tolerance.
    pub coincident_vertices: Vec<CoincidentPair>,
    /// Vertices lying on an edge's interior.
    pub t_junctions: Vec<TJunction>,
    /// Ill-formed NURBS curves/pcurves.
    pub nurbs_defects: Vec<NurbsDefect>,
    /// Pcurves whose evaluations fall outside their face surface's
    /// parameter range.
    pub pcurve_escapes: Vec<PcurveEscape>,
    /// Pairs of faces with the same analytic surface and overlapping
    /// AABBs — strong indicator of a duplicate / un-merged face left
    /// behind by a boolean / fillet operation.
    pub face_overlaps: Vec<FaceOverlap>,
    /// Count of edges flagged by `find_tessellation_non_manifold_edges`
    /// using a 5×tolerance threshold. Fold-in: the full
    /// per-edge listing is available via `report_non_manifold`; this
    /// is just the headline number for the diagnostic block.
    pub tess_non_manifold_count: usize,
    /// Largest interior-divergence value observed across all edges.
    pub tess_non_manifold_max: f64,
}

/// A naked edge.
#[derive(Debug, Clone, Copy)]
pub struct NakedEdge {
    /// The edge.
    pub edge: EdgeHandle,
    /// Face on the primary half-edge side.
    pub face_a: FaceHandle,
    /// Whether `face_a` has a surface assigned.
    pub face_a_has_surface: bool,
    /// Face on the twin side.
    pub face_b: FaceHandle,
    /// Whether `face_b` has a surface assigned.
    pub face_b_has_surface: bool,
    /// Midpoint of the edge.
    pub midpoint: Point3,
}

/// A degenerate (very short) edge.
#[derive(Debug, Clone, Copy)]
pub struct DegenerateEdge {
    /// The edge.
    pub edge: EdgeHandle,
    /// Computed length.
    pub length: f64,
    /// Midpoint of the edge.
    pub midpoint: Point3,
}

/// A sliver face — small area or bad aspect ratio.
#[derive(Debug, Clone, Copy)]
pub struct SliverFace {
    /// The face.
    pub face: FaceHandle,
    /// Approximate area (fan-triangulation).
    pub area: f64,
    /// Ratio of longest to shortest outer-loop edge.
    pub aspect_ratio: f64,
    /// Centroid of the outer-loop vertex positions.
    pub centroid: Point3,
}

/// A pair of distinct vertex handles within tolerance.
#[derive(Debug, Clone, Copy)]
pub struct CoincidentPair {
    /// First vertex.
    pub a: VertexHandle,
    /// Second vertex.
    pub b: VertexHandle,
    /// 3-D distance between the two positions.
    pub distance: f64,
    /// Midpoint between the two.
    pub midpoint: Point3,
}

/// A vertex lying on the *interior* of an edge it isn't an endpoint
/// of. T-junctions break manifold-edge invariants and confuse
/// OCCT-strict slicers.
#[derive(Debug, Clone, Copy)]
pub struct TJunction {
    /// The stray vertex.
    pub vertex: VertexHandle,
    /// The edge whose interior it touches.
    pub edge: EdgeHandle,
    /// Parameter value on the edge curve at which the vertex lies.
    pub parameter: f64,
    /// 3-D distance from the vertex to the curve.
    pub distance: f64,
    /// 3-D position of the offending vertex.
    pub vertex_position: Point3,
}

/// Where a NURBS irregularity was found.
#[derive(Debug, Clone, Copy)]
pub enum NurbsScope {
    /// On an edge's 3-D curve.
    EdgeCurve(EdgeHandle),
    /// On a half-edge's pcurve.
    HalfEdgePcurve(HalfEdgeHandle),
    /// On a face's NURBS surface.
    FaceSurface(FaceHandle),
}

/// Class of NURBS irregularity.
#[derive(Debug, Clone, Copy)]
pub enum NurbsKind {
    /// `weights.len() != control_points.len()`.
    CpWeightMismatch {
        /// Number of control points.
        cps: usize,
        /// Number of weights.
        weights: usize,
    },
    /// `knots.len() != cps + degree + 1`.
    CpKnotMismatch {
        /// Number of control points.
        cps: usize,
        /// Curve degree.
        degree: usize,
        /// Number of knots.
        knots: usize,
    },
    /// Knot vector contains a non-finite value.
    NonFiniteKnot {
        /// Index of the offending knot.
        index: usize,
    },
    /// Knot vector is not non-decreasing at this index.
    NonMonotonicKnots {
        /// Index where the decrease happens.
        index: usize,
    },
    /// A run of identical knots exceeds `degree + 1`.
    HighKnotMultiplicity {
        /// Curve degree.
        degree: usize,
        /// Length of the offending run.
        run: usize,
    },
    /// A weight is zero, negative, or non-finite.
    InvalidWeight {
        /// Index of the offending weight.
        index: usize,
        /// Its value.
        weight: f64,
    },
    /// A control point coordinate is non-finite.
    NonFiniteControlPoint {
        /// Index of the offending control point.
        index: usize,
    },
}

/// One ill-formed NURBS instance.
#[derive(Debug, Clone, Copy)]
pub struct NurbsDefect {
    /// Where the defect lives.
    pub scope: NurbsScope,
    /// Class of irregularity.
    pub kind: NurbsKind,
}

/// A half-edge pcurve that produces (u, v) values outside the
/// surface's bounded parameter range.
#[derive(Debug, Clone, Copy)]
pub struct PcurveEscape {
    /// The half-edge whose pcurve escapes.
    pub half_edge: HalfEdgeHandle,
    /// The half-edge's parent edge.
    pub edge: EdgeHandle,
    /// The face whose surface UV the pcurve is supposed to live in.
    pub face: FaceHandle,
    /// Pcurve parameter value at which the worst escape occurs.
    pub parameter: f64,
    /// The (u, v, 0) value at that parameter.
    pub uv: Point3,
    /// The surface's `((u_min, u_max), (v_min, v_max))` range.
    pub uv_range: ((f64, f64), (f64, f64)),
    /// Largest absolute UV-axis overshoot found across samples.
    pub max_overshoot: f64,
}

/// A pair of faces with the same surface (within tolerance) and
/// overlapping bounding boxes — strong indicator of a duplicate or
/// un-merged face left behind by a boolean / fillet operation.
#[derive(Debug, Clone, Copy)]
pub struct FaceOverlap {
    /// First face.
    pub face_a: FaceHandle,
    /// Second face.
    pub face_b: FaceHandle,
    /// Whether the two faces share the same analytic surface
    /// definition (within tolerance).
    pub same_surface: bool,
    /// Distance between the two centroids.
    pub centroid_distance: f64,
    /// Surface kind label ("Plane", "Cyl", …).
    pub surface_kind: &'static str,
}

impl Diagnostics {
    /// Run all detectors and return the aggregate report.
    pub fn run(brep: &BRep) -> Self {
        Self {
            topology: detectors::topology_stats(brep),
            geometry: detectors::geometry_stats(brep),
            defects: detectors::detect_defects(brep),
        }
    }
}

impl fmt::Display for Diagnostics {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let t = &self.topology;
        let g = &self.geometry;
        let d = &self.defects;
        let kinds: Vec<String> = t
            .surface_kinds
            .iter()
            .map(|(k, n)| format!("{k}={n}"))
            .collect();
        writeln!(
            f,
            "  Topology: {} faces ({} with surface; {}), {} edges, {} verts; χ=V-E+F={}",
            t.faces,
            t.faces_with_surface,
            kinds.join(", "),
            t.edges,
            t.vertices,
            t.euler_chi
        )?;
        writeln!(
            f,
            "  Geometry: edge len min/med/max = {:.4}/{:.4}/{:.4}; face area min/med/max = {:.4}/{:.4}/{:.4}",
            g.edge_length_min,
            g.edge_length_median,
            g.edge_length_max,
            g.face_area_min,
            g.face_area_median,
            g.face_area_max
        )?;
        writeln!(
            f,
            "  BBox:     [{:+.3},{:+.3},{:+.3}] → [{:+.3},{:+.3},{:+.3}] diag={:.4}",
            g.bbox_min.x,
            g.bbox_min.y,
            g.bbox_min.z,
            g.bbox_max.x,
            g.bbox_max.y,
            g.bbox_max.z,
            g.bbox_diagonal
        )?;
        let n_naked = d.naked_edges.len();
        let n_degen = d.degenerate_edges.len();
        let n_sliver = d.sliver_faces.len();
        let n_coinc = d.coincident_vertices.len();
        let n_tj = d.t_junctions.len();
        let n_nurbs = d.nurbs_defects.len();
        let n_esc = d.pcurve_escapes.len();
        let n_overlap = d.face_overlaps.len();
        let n_tess = d.tess_non_manifold_count;
        let total =
            n_naked + n_degen + n_sliver + n_coinc + n_tj + n_nurbs + n_esc + n_overlap + n_tess;
        if total == 0 {
            writeln!(f, "  Defects:  none ✓")?;
            return Ok(());
        }
        writeln!(
            f,
            "  Defects:  {n_naked} naked, {n_degen} degen-edges, {n_sliver} slivers, {n_coinc} coincidents, {n_tj} T-junctions, {n_nurbs} nurbs, {n_esc} pcurve-escape, {n_overlap} face-overlaps, {n_tess} tess-nm (max {:.3e})",
            d.tess_non_manifold_max
        )?;
        for n in d.naked_edges.iter().take(2) {
            writeln!(
                f,
                "    naked: {:?} faces={:?}|{:?} (surf {} | {}) at ({:+.3},{:+.3},{:+.3})",
                n.edge,
                n.face_a,
                n.face_b,
                if n.face_a_has_surface { "y" } else { "n" },
                if n.face_b_has_surface { "y" } else { "n" },
                n.midpoint.x,
                n.midpoint.y,
                n.midpoint.z
            )?;
        }
        for de in d.degenerate_edges.iter().take(2) {
            writeln!(
                f,
                "    degenerate: {:?} len={:.2e} at ({:+.3},{:+.3},{:+.3})",
                de.edge, de.length, de.midpoint.x, de.midpoint.y, de.midpoint.z
            )?;
        }
        for sf in d.sliver_faces.iter().take(2) {
            writeln!(
                f,
                "    sliver: {:?} area={:.2e} aspect={:.1} at ({:+.3},{:+.3},{:+.3})",
                sf.face, sf.area, sf.aspect_ratio, sf.centroid.x, sf.centroid.y, sf.centroid.z
            )?;
        }
        for cp in d.coincident_vertices.iter().take(2) {
            writeln!(
                f,
                "    coincident: {:?}~{:?} d={:.2e} at ({:+.3},{:+.3},{:+.3})",
                cp.a, cp.b, cp.distance, cp.midpoint.x, cp.midpoint.y, cp.midpoint.z
            )?;
        }
        for tj in d.t_junctions.iter().take(2) {
            writeln!(
                f,
                "    T-junction: {:?} on {:?} t={:.3} d={:.2e} at ({:+.3},{:+.3},{:+.3})",
                tj.vertex,
                tj.edge,
                tj.parameter,
                tj.distance,
                tj.vertex_position.x,
                tj.vertex_position.y,
                tj.vertex_position.z
            )?;
        }
        for nd in d.nurbs_defects.iter().take(2) {
            writeln!(f, "    nurbs: {:?} {:?}", nd.scope, nd.kind)?;
        }
        for pe in d.pcurve_escapes.iter().take(2) {
            writeln!(
                f,
                "    pcurve-escape: {:?} on {:?} uv=({:+.3},{:+.3}) range=u[{:.2},{:.2}] v[{:.2},{:.2}] over={:.2e}",
                pe.edge,
                pe.face,
                pe.uv.x,
                pe.uv.y,
                pe.uv_range.0.0,
                pe.uv_range.0.1,
                pe.uv_range.1.0,
                pe.uv_range.1.1,
                pe.max_overshoot
            )?;
        }
        for fo in d.face_overlaps.iter().take(2) {
            writeln!(
                f,
                "    face-overlap: {:?}~{:?} ({}) cdist={:.4}",
                fo.face_a, fo.face_b, fo.surface_kind, fo.centroid_distance
            )?;
        }
        Ok(())
    }
}

#[cfg(test)]
#[path = "diagnostics_tests.rs"]
mod tests;
