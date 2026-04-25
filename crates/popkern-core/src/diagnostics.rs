//! Comprehensive BRep diagnostics — a one-shot health card.
//!
//! Where `BRep::validate` answers a binary "is this BRep structurally
//! and geometrically sound?", `Diagnostics::run` answers the harder
//! question: "what's *interesting* about this BRep?" It surfaces
//! near-degenerate geometry, suspicious sliver faces, near-duplicate
//! vertices, and surface-kind histograms — the kind of summary you
//! want to print after a build to know whether the output is worth
//! exporting at all.
//!
//! This is Round 1: topology + geometry summaries + naked edges +
//! degenerate edges + sliver faces + coincident vertices. Future
//! rounds add T-junction detection, NURBS irregularities, pcurve
//! UV-escape, face overlap, and so on.
//!
//! The detectors run in O(F + E + V²) — the V² coincident-vertex
//! pairwise check dominates and is fine for the BRep sizes seen in
//! showcase/test fixtures (tens to low hundreds of vertices).

use std::collections::BTreeMap;
use std::fmt;

use crate::geometry::Point3;
use crate::geometry::curve::CurveEval;
use crate::geometry::surface::SurfaceGeometry;
use crate::topology::{BRep, EdgeHandle, FaceHandle, VertexHandle};

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
/// unit fix).
#[derive(Debug, Clone, Default)]
pub struct GeometryStats {
    /// Smallest edge length observed (computed from edge curves where
    /// available, else from vertex-position distance).
    pub edge_length_min: f64,
    /// Median edge length.
    pub edge_length_median: f64,
    /// Largest edge length.
    pub edge_length_max: f64,
    /// Smallest face area observed (planar approximation of outer
    /// loop via fan triangulation).
    pub face_area_min: f64,
    /// Median face area.
    pub face_area_median: f64,
    /// Largest face area.
    pub face_area_max: f64,
}

/// Aggregate of all defect categories. Empty vectors mean no defects
/// in that category.
#[derive(Debug, Clone, Default)]
pub struct Defects {
    /// Edges where one or both adjacent faces have no surface
    /// assigned — the half-edge structure is intact but the
    /// geometry is incomplete. OCCT-strict importers often refuse
    /// these.
    pub naked_edges: Vec<NakedEdge>,
    /// Edges shorter than `tolerance.linear()`. These typically
    /// indicate failed boolean / fillet operations or stale
    /// near-coincident vertices.
    pub degenerate_edges: Vec<DegenerateEdge>,
    /// Faces with area smaller than `tolerance.linear()²`, or with
    /// extreme aspect ratio (longest-edge / shortest-edge > 1000).
    /// Slicers struggle with these and they often indicate sliver
    /// remnants from boolean cuts.
    pub sliver_faces: Vec<SliverFace>,
    /// Pairs of distinct vertex handles whose positions are within
    /// `tolerance.linear()` of each other. Likely a fillet/blend
    /// failed to merge these via Euler ops.
    pub coincident_vertices: Vec<CoincidentPair>,
}

/// A naked edge. `face_a_has_surface == false` or `face_b_has_surface == false`.
#[derive(Debug, Clone, Copy)]
pub struct NakedEdge {
    /// The edge.
    pub edge: EdgeHandle,
    /// Face on the primary half-edge side (always present).
    pub face_a: FaceHandle,
    /// Whether face_a has a surface assigned.
    pub face_a_has_surface: bool,
    /// Face on the twin side.
    pub face_b: FaceHandle,
    /// Whether face_b has a surface assigned.
    pub face_b_has_surface: bool,
    /// Midpoint of the edge — useful for visualisation.
    pub midpoint: Point3,
}

/// A degenerate (very short) edge.
#[derive(Debug, Clone, Copy)]
pub struct DegenerateEdge {
    /// The edge.
    pub edge: EdgeHandle,
    /// Computed length. Less than `tolerance.linear()`.
    pub length: f64,
    /// Midpoint between the edge's two vertex positions.
    pub midpoint: Point3,
}

/// A sliver face — small area or bad aspect ratio.
#[derive(Debug, Clone, Copy)]
pub struct SliverFace {
    /// The face.
    pub face: FaceHandle,
    /// Approximate area (fan-triangulation of outer loop).
    pub area: f64,
    /// Ratio of longest to shortest edge in the outer loop. Values
    /// above ~1000 indicate sliver geometry even at non-trivial area.
    pub aspect_ratio: f64,
    /// Centroid of the outer-loop vertex positions.
    pub centroid: Point3,
}

/// A pair of distinct vertex handles within tolerance of each other.
#[derive(Debug, Clone, Copy)]
pub struct CoincidentPair {
    /// First vertex.
    pub a: VertexHandle,
    /// Second vertex.
    pub b: VertexHandle,
    /// 3-D distance between the two vertex positions.
    pub distance: f64,
    /// Midpoint between the two — useful for visualisation.
    pub midpoint: Point3,
}

impl Diagnostics {
    /// Run all Round-1 detectors and return the aggregate report.
    pub fn run(brep: &BRep) -> Self {
        Self {
            topology: topology_stats(brep),
            geometry: geometry_stats(brep),
            defects: detect_defects(brep),
        }
    }
}

// ---------------------------------------------------------------------------
// Detectors
// ---------------------------------------------------------------------------

fn topology_stats(brep: &BRep) -> TopologyStats {
    let mut surface_kinds: BTreeMap<&'static str, usize> = BTreeMap::new();
    let mut faces_with_surface = 0;
    for (_, face) in &brep.faces {
        if let Some(sof) = face.surface.as_ref() {
            faces_with_surface += 1;
            let key = surface_kind_label(&sof.geometry);
            *surface_kinds.entry(key).or_insert(0) += 1;
        }
    }
    let v = brep.vertices.len() as i64;
    let e = brep.edges.len() as i64;
    let f = brep.faces.len() as i64;
    TopologyStats {
        faces: brep.faces.len(),
        faces_with_surface,
        edges: brep.edges.len(),
        vertices: brep.vertices.len(),
        half_edges: brep.half_edges.len(),
        euler_chi: v - e + f,
        surface_kinds,
    }
}

fn geometry_stats(brep: &BRep) -> GeometryStats {
    let edge_lengths: Vec<f64> = brep.edges.keys().map(|eh| edge_length(brep, eh)).collect();
    let face_areas: Vec<f64> = brep.faces.keys().map(|fh| face_area(brep, fh)).collect();
    let (e_min, e_med, e_max) = min_median_max(&edge_lengths);
    let (f_min, f_med, f_max) = min_median_max(&face_areas);
    GeometryStats {
        edge_length_min: e_min,
        edge_length_median: e_med,
        edge_length_max: e_max,
        face_area_min: f_min,
        face_area_median: f_med,
        face_area_max: f_max,
    }
}

fn detect_defects(brep: &BRep) -> Defects {
    Defects {
        naked_edges: detect_naked_edges(brep),
        degenerate_edges: detect_degenerate_edges(brep),
        sliver_faces: detect_sliver_faces(brep),
        coincident_vertices: detect_coincident_vertices(brep),
    }
}

fn detect_naked_edges(brep: &BRep) -> Vec<NakedEdge> {
    let mut out = Vec::new();
    for (eh, edge) in &brep.edges {
        let Some(he) = brep.half_edges.get(edge.half_edge) else {
            continue;
        };
        let Some(twin) = brep.half_edges.get(he.twin) else {
            continue;
        };
        let face_a = brep.loops[he.loop_].face;
        let face_b = brep.loops[twin.loop_].face;
        let face_a_has = brep.faces.get(face_a).is_some_and(|f| f.surface.is_some());
        let face_b_has = brep.faces.get(face_b).is_some_and(|f| f.surface.is_some());
        if face_a_has && face_b_has {
            continue;
        }
        let midpoint = edge_midpoint(brep, eh);
        out.push(NakedEdge {
            edge: eh,
            face_a,
            face_a_has_surface: face_a_has,
            face_b,
            face_b_has_surface: face_b_has,
            midpoint,
        });
    }
    out
}

fn detect_degenerate_edges(brep: &BRep) -> Vec<DegenerateEdge> {
    let tol = brep.tolerance.linear().max(1e-12);
    let mut out = Vec::new();
    for (eh, _) in &brep.edges {
        let length = edge_length(brep, eh);
        if length < tol {
            out.push(DegenerateEdge {
                edge: eh,
                length,
                midpoint: edge_midpoint(brep, eh),
            });
        }
    }
    out
}

fn detect_sliver_faces(brep: &BRep) -> Vec<SliverFace> {
    let tol = brep.tolerance.linear().max(1e-12);
    let area_floor = tol * tol * 100.0; // 100× tol² gives a generous "sliver" cutoff
    let aspect_ceiling = 1000.0;
    let mut out = Vec::new();
    for (fh, _) in &brep.faces {
        let area = face_area(brep, fh);
        let aspect = face_edge_aspect_ratio(brep, fh);
        let centroid = face_centroid(brep, fh);
        if area < area_floor || aspect > aspect_ceiling {
            out.push(SliverFace {
                face: fh,
                area,
                aspect_ratio: aspect,
                centroid,
            });
        }
    }
    out
}

fn detect_coincident_vertices(brep: &BRep) -> Vec<CoincidentPair> {
    let tol = brep.tolerance.linear().max(1e-12);
    let handles: Vec<VertexHandle> = brep.vertices.keys().collect();
    let mut out = Vec::new();
    for i in 0..handles.len() {
        let a = handles[i];
        let pa = brep.vertices[a].position;
        for &b in &handles[i + 1..] {
            let pb = brep.vertices[b].position;
            let d = pa.distance_to(&pb);
            if d < tol {
                out.push(CoincidentPair {
                    a,
                    b,
                    distance: d,
                    midpoint: Point3::new(
                        0.5 * (pa.x + pb.x),
                        0.5 * (pa.y + pb.y),
                        0.5 * (pa.z + pb.z),
                    ),
                });
            }
        }
    }
    out
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

fn surface_kind_label(sg: &SurfaceGeometry) -> &'static str {
    match sg {
        SurfaceGeometry::Plane { .. } => "Plane",
        SurfaceGeometry::Cylinder { .. } => "Cyl",
        SurfaceGeometry::Sphere { .. } => "Sph",
        SurfaceGeometry::Cone { .. } => "Cone",
        SurfaceGeometry::Torus { .. } => "Torus",
        SurfaceGeometry::Nurbs(_) => "Nurbs",
        SurfaceGeometry::Gregory(_) => "Gregory",
    }
}

fn edge_length(brep: &BRep, eh: EdgeHandle) -> f64 {
    let Some(edge) = brep.edges.get(eh) else {
        return 0.0;
    };
    if let Some(curve_on_edge) = edge.curve.as_ref() {
        let len = curve_on_edge
            .geometry
            .length_between(curve_on_edge.t_start, curve_on_edge.t_end);
        if len.is_finite() && len > 0.0 {
            return len;
        }
    }
    // Fallback: distance between the edge's two vertex positions.
    let Some(he) = brep.half_edges.get(edge.half_edge) else {
        return 0.0;
    };
    let Some(twin) = brep.half_edges.get(he.twin) else {
        return 0.0;
    };
    let pa = brep.vertices[he.origin].position;
    let pb = brep.vertices[twin.origin].position;
    pa.distance_to(&pb)
}

fn edge_midpoint(brep: &BRep, eh: EdgeHandle) -> Point3 {
    let Some(edge) = brep.edges.get(eh) else {
        return Point3::new(0.0, 0.0, 0.0);
    };
    if let Some(curve_on_edge) = edge.curve.as_ref() {
        let t_mid = 0.5 * (curve_on_edge.t_start + curve_on_edge.t_end);
        return curve_on_edge.geometry.point_at(t_mid);
    }
    let Some(he) = brep.half_edges.get(edge.half_edge) else {
        return Point3::new(0.0, 0.0, 0.0);
    };
    let Some(twin) = brep.half_edges.get(he.twin) else {
        return Point3::new(0.0, 0.0, 0.0);
    };
    let pa = brep.vertices[he.origin].position;
    let pb = brep.vertices[twin.origin].position;
    Point3::new(
        0.5 * (pa.x + pb.x),
        0.5 * (pa.y + pb.y),
        0.5 * (pa.z + pb.z),
    )
}

/// Approximate face area by fan-triangulating the outer-loop vertex
/// positions in 3-D from the first vertex and summing the unsigned
/// areas of each triangle. Works for any face shape and ignores
/// surface curvature — it's intentionally a cheap defect detector,
/// not an exact-area computation.
fn face_area(brep: &BRep, fh: FaceHandle) -> f64 {
    let positions = brep.face_outer_loop_positions(fh);
    if positions.len() < 3 {
        return 0.0;
    }
    let p0 = positions[0];
    let mut area = 0.0;
    for i in 1..positions.len() - 1 {
        let p1 = positions[i];
        let p2 = positions[i + 1];
        let ax = p1.x - p0.x;
        let ay = p1.y - p0.y;
        let az = p1.z - p0.z;
        let bx = p2.x - p0.x;
        let by = p2.y - p0.y;
        let bz = p2.z - p0.z;
        // |a × b| / 2
        let cx = ay * bz - az * by;
        let cy = az * bx - ax * bz;
        let cz = ax * by - ay * bx;
        area += 0.5 * (cx * cx + cy * cy + cz * cz).sqrt();
    }
    area
}

fn face_edge_aspect_ratio(brep: &BRep, fh: FaceHandle) -> f64 {
    let positions = brep.face_outer_loop_positions(fh);
    if positions.len() < 2 {
        return 0.0;
    }
    let mut min_len = f64::INFINITY;
    let mut max_len = 0.0_f64;
    for i in 0..positions.len() {
        let pa = positions[i];
        let pb = positions[(i + 1) % positions.len()];
        let len = pa.distance_to(&pb);
        if len > 0.0 && len < min_len {
            min_len = len;
        }
        if len > max_len {
            max_len = len;
        }
    }
    if min_len <= 0.0 || !min_len.is_finite() {
        return f64::INFINITY;
    }
    max_len / min_len
}

fn face_centroid(brep: &BRep, fh: FaceHandle) -> Point3 {
    let positions = brep.face_outer_loop_positions(fh);
    if positions.is_empty() {
        return Point3::new(0.0, 0.0, 0.0);
    }
    let n = positions.len() as f64;
    let (mut x, mut y, mut z) = (0.0_f64, 0.0_f64, 0.0_f64);
    for p in &positions {
        x += p.x;
        y += p.y;
        z += p.z;
    }
    Point3::new(x / n, y / n, z / n)
}

fn min_median_max(values: &[f64]) -> (f64, f64, f64) {
    if values.is_empty() {
        return (0.0, 0.0, 0.0);
    }
    let mut sorted: Vec<f64> = values.iter().copied().filter(|v| v.is_finite()).collect();
    if sorted.is_empty() {
        return (0.0, 0.0, 0.0);
    }
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let min = sorted[0];
    let max = sorted[sorted.len() - 1];
    let median = sorted[sorted.len() / 2];
    (min, median, max)
}

// ---------------------------------------------------------------------------
// Pretty-printing
// ---------------------------------------------------------------------------

impl fmt::Display for Diagnostics {
    /// Produce a multi-line diagnostic block. Indented with two
    /// leading spaces so it nests cleanly inside per-model blocks
    /// in showcase output.
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
        let n_naked = d.naked_edges.len();
        let n_degen = d.degenerate_edges.len();
        let n_sliver = d.sliver_faces.len();
        let n_coinc = d.coincident_vertices.len();
        let total = n_naked + n_degen + n_sliver + n_coinc;
        if total == 0 {
            writeln!(f, "  Defects:  none ✓")?;
            return Ok(());
        }
        writeln!(
            f,
            "  Defects:  {n_naked} naked, {n_degen} degenerate edges, {n_sliver} sliver faces, {n_coinc} coincident vertices"
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
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tolerance::Tolerance;

    #[test]
    fn unit_cube_has_no_defects() {
        let brep = BRep::unit_cube(Tolerance::default()).unwrap();
        let report = Diagnostics::run(&brep);
        assert_eq!(report.topology.faces, 6);
        assert_eq!(report.topology.vertices, 8);
        assert_eq!(report.topology.edges, 12);
        // Genus-0 closed: V - E + F = 8 - 12 + 6 = 2
        assert_eq!(report.topology.euler_chi, 2);
        assert!(report.defects.naked_edges.is_empty());
        assert!(report.defects.degenerate_edges.is_empty());
        assert!(report.defects.coincident_vertices.is_empty());
    }

    #[test]
    fn surface_kinds_histogram_counts_planes() {
        let brep = BRep::unit_cube(Tolerance::default()).unwrap();
        let report = Diagnostics::run(&brep);
        // Cube has 6 planar faces.
        assert_eq!(report.topology.surface_kinds.get("Plane"), Some(&6));
    }

    #[test]
    fn coincident_vertices_detected() {
        // Two vertices placed within tolerance; should show up as a
        // coincident pair.
        let mut brep = BRep::new(Tolerance::default());
        let _v1 = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
        // Manually insert a near-duplicate vertex (no edges/loops —
        // the diagnostic walks brep.vertices directly).
        let v2 = brep.vertices.insert(crate::topology::Vertex {
            position: Point3::new(1e-12, 0.0, 0.0),
            outgoing: Default::default(),
            tolerance: None,
            name: None,
        });
        let report = Diagnostics::run(&brep);
        assert!(
            !report.defects.coincident_vertices.is_empty(),
            "expected coincident pair detection"
        );
        let pair = report.defects.coincident_vertices[0];
        // The detected pair must include the new vertex.
        assert!(pair.a == v2 || pair.b == v2);
    }
}
