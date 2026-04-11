//! Debug visualization overlay for entity imprecision regions.
//!
//! Generates geometric representations (spheres around vertices, tubes along
//! edges) that show where tracked imprecision exists in a BRep model.

use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::CurveEval;
use popkern_core::topology::BRep;

/// Number of sample points along an edge curve for tube polylines.
const EDGE_SAMPLE_COUNT: usize = 20;

/// Debug overlay geometry representing entity imprecision regions.
///
/// Each vertex with explicitly tracked tolerance produces a sphere
/// (center + radius). Each edge with explicitly tracked tolerance
/// and a curve produces a tube (polyline + radius).
///
/// Face bands are omitted in v1 (they require tessellation complexity).
pub struct ImprecisionOverlay {
    /// `(center, radius)` per vertex with explicit imprecision.
    pub spheres: Vec<(Point3, f64)>,
    /// `(polyline along edge, tube radius)` per edge with explicit imprecision.
    pub tubes: Vec<(Vec<Point3>, f64)>,
}

impl ImprecisionOverlay {
    /// Generate overlay for entities with explicitly tracked (non-`None`) imprecision.
    ///
    /// Iterates all vertices and edges in the BRep. Vertices with
    /// `tolerance.is_some()` contribute a sphere entry. Edges with
    /// `tolerance.is_some()` and an attached curve contribute a tube
    /// entry sampled at [`EDGE_SAMPLE_COUNT`] points.
    pub fn from_brep(brep: &BRep) -> Self {
        let mut spheres = Vec::new();
        let mut tubes = Vec::new();

        for (_handle, vertex) in &brep.vertices {
            if let Some(tol) = vertex.tolerance {
                spheres.push((vertex.position, tol));
            }
        }

        for (_handle, edge) in &brep.edges {
            if let Some(tol) = edge.tolerance
                && let Some(ref curve_on_edge) = edge.curve
            {
                let polyline = sample_curve(
                    &curve_on_edge.geometry,
                    curve_on_edge.t_start,
                    curve_on_edge.t_end,
                );
                tubes.push((polyline, tol));
            }
        }

        Self { spheres, tubes }
    }

    /// Only entities using more than `budget_fraction` of `max_imprecision`.
    ///
    /// Builds the full overlay via [`Self::from_brep`], then retains only
    /// entries whose radius exceeds `brep.max_imprecision * budget_fraction`.
    pub fn above_budget_fraction(brep: &BRep, budget_fraction: f64) -> Self {
        let threshold = brep.max_imprecision * budget_fraction;
        let mut overlay = Self::from_brep(brep);
        overlay.spheres.retain(|(_, r)| *r > threshold);
        overlay.tubes.retain(|(_, r)| *r > threshold);
        overlay
    }
}

/// Sample `count` evenly-spaced points along a curve between `t_start` and `t_end`.
fn sample_curve(curve: &dyn CurveEval, t_start: f64, t_end: f64) -> Vec<Point3> {
    let count = EDGE_SAMPLE_COUNT;
    let mut points = Vec::with_capacity(count);
    for i in 0..count {
        let frac = if count <= 1 {
            0.0
        } else {
            i as f64 / (count - 1) as f64
        };
        let t = t_start + (t_end - t_start) * frac;
        points.push(curve.point_at(t));
    }
    points
}

#[cfg(test)]
mod tests {
    use super::*;
    use popkern_core::geometry::Vec3;
    use popkern_core::geometry::curve::{CurveGeometry, CurveOnEdge};
    use popkern_core::tolerance::Tolerance;
    use popkern_core::topology::BRep;

    /// Helper: build a minimal BRep with one vertex that has explicit tolerance.
    fn brep_with_tolerant_vertex() -> BRep {
        let tol = Tolerance::default();
        let mut brep = BRep::new(tol);

        let result = brep
            .mvfs(Point3::new(1.0, 2.0, 3.0))
            .expect("mvfs should succeed in test");
        if let Some(v) = brep.vertices.get_mut(result.vertex) {
            v.tolerance = Some(0.01);
        }

        brep
    }

    /// Helper: build a BRep with an edge that has explicit tolerance and a line curve.
    fn brep_with_tolerant_edge() -> BRep {
        let tol = Tolerance::default();
        let mut brep = BRep::new(tol);

        let r = brep
            .mvfs(Point3::new(0.0, 0.0, 0.0))
            .expect("mvfs should succeed in test");
        let r2 = brep
            .mev_at_vertex(r.vertex, r.loop_, Point3::new(1.0, 0.0, 0.0))
            .expect("mev_at_vertex should succeed in test");

        // Attach a line curve and tolerance to the edge.
        if let Some(edge) = brep.edges.get_mut(r2.edge) {
            edge.tolerance = Some(0.005);
            edge.curve = Some(CurveOnEdge {
                geometry: CurveGeometry::Line {
                    origin: Point3::new(0.0, 0.0, 0.0),
                    direction: Vec3::new(1.0, 0.0, 0.0),
                },
                t_start: 0.0,
                t_end: 1.0,
            });
        }

        brep
    }

    #[test]
    fn from_brep_empty() {
        let brep = BRep::new(Tolerance::default());
        let overlay = ImprecisionOverlay::from_brep(&brep);
        assert!(overlay.spheres.is_empty());
        assert!(overlay.tubes.is_empty());
    }

    #[test]
    fn from_brep_vertex_no_tolerance() {
        let tol = Tolerance::default();
        let mut brep = BRep::new(tol);
        let _r = brep
            .mvfs(Point3::new(0.0, 0.0, 0.0))
            .expect("mvfs should succeed in test");
        // vertex.tolerance is None by default
        let overlay = ImprecisionOverlay::from_brep(&brep);
        assert!(overlay.spheres.is_empty());
    }

    #[test]
    fn from_brep_vertex_with_tolerance() {
        let brep = brep_with_tolerant_vertex();
        let overlay = ImprecisionOverlay::from_brep(&brep);
        assert_eq!(overlay.spheres.len(), 1);
        let (center, radius) = &overlay.spheres[0];
        assert!((center.x - 1.0).abs() < 1e-10);
        assert!((center.y - 2.0).abs() < 1e-10);
        assert!((center.z - 3.0).abs() < 1e-10);
        assert!((radius - 0.01).abs() < 1e-10);
    }

    #[test]
    fn from_brep_edge_with_tolerance_and_curve() {
        let brep = brep_with_tolerant_edge();
        let overlay = ImprecisionOverlay::from_brep(&brep);
        assert_eq!(overlay.tubes.len(), 1);
        let (polyline, radius) = &overlay.tubes[0];
        assert_eq!(polyline.len(), EDGE_SAMPLE_COUNT);
        assert!((radius - 0.005).abs() < 1e-10);
        // First point should be at origin, last at (1,0,0)
        assert!(polyline[0].distance_to(&Point3::new(0.0, 0.0, 0.0)) < 1e-10);
        assert!(polyline[EDGE_SAMPLE_COUNT - 1].distance_to(&Point3::new(1.0, 0.0, 0.0)) < 1e-10);
    }

    #[test]
    fn above_budget_fraction_filters() {
        let tol = Tolerance::default();
        let mut brep = BRep::new(tol);
        // max_imprecision defaults to tolerance.linear() * 100.0
        // With default tolerance (1e-6), max_imprecision = 1e-4

        let r1 = brep
            .mvfs(Point3::new(0.0, 0.0, 0.0))
            .expect("mvfs should succeed in test");
        if let Some(v) = brep.vertices.get_mut(r1.vertex) {
            v.tolerance = Some(0.00008); // 80% of 1e-4
        }

        // Add another vertex with lower tolerance
        let r2 = brep
            .mev_at_vertex(r1.vertex, r1.loop_, Point3::new(1.0, 0.0, 0.0))
            .expect("mev_at_vertex should succeed in test");
        if let Some(v) = brep.vertices.get_mut(r2.vertex) {
            v.tolerance = Some(0.00002); // 20% of 1e-4
        }

        // budget_fraction = 0.5 => threshold = 0.00005
        let overlay = ImprecisionOverlay::above_budget_fraction(&brep, 0.5);
        // Only the first vertex (0.00008 > 0.00005) should remain
        assert_eq!(overlay.spheres.len(), 1);
        assert!((overlay.spheres[0].1 - 0.00008).abs() < 1e-10);
    }

    #[test]
    fn edge_without_curve_skipped() {
        let tol = Tolerance::default();
        let mut brep = BRep::new(tol);

        let r = brep
            .mvfs(Point3::new(0.0, 0.0, 0.0))
            .expect("mvfs should succeed in test");
        let r2 = brep
            .mev_at_vertex(r.vertex, r.loop_, Point3::new(1.0, 0.0, 0.0))
            .expect("mev_at_vertex should succeed in test");

        // Set tolerance but no curve
        if let Some(edge) = brep.edges.get_mut(r2.edge) {
            edge.tolerance = Some(0.01);
            // edge.curve is None by default
        }

        let overlay = ImprecisionOverlay::from_brep(&brep);
        assert!(overlay.tubes.is_empty());
    }

    #[test]
    fn sample_curve_single_point() {
        // Verify the sample_curve function handles edge case gracefully
        let curve = CurveGeometry::Line {
            origin: Point3::new(0.0, 0.0, 0.0),
            direction: Vec3::new(1.0, 0.0, 0.0),
        };
        let points = super::sample_curve(&curve, 0.0, 1.0);
        assert_eq!(points.len(), EDGE_SAMPLE_COUNT);
    }
}
