//! Geometric criteria for stable topological naming.
//!
//! These types capture geometric properties of topological entities at the
//! time they are named. The criteria serve as tie-breakers when structural
//! names collide during feature-tree resolution — i.e. when two entities have
//! the same [`EntityName`] recipe but the modeller needs to pick the closest
//! match by position/shape.

use serde::{Deserialize, Serialize};

use crate::geometry::{Point3, Vec3};
use crate::topology::{BRep, EdgeHandle, FaceHandle, VertexHandle};

// ---------------------------------------------------------------------------
// Surface / curve type tags
// ---------------------------------------------------------------------------

/// The kind of parametric surface backing a face.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum SurfaceType {
    /// Flat plane.
    Plane,
    /// Circular cylinder.
    Cylinder,
    /// Circular cone.
    Cone,
    /// Sphere.
    Sphere,
    /// Torus.
    Torus,
    /// NURBS / B-spline surface.
    Nurbs,
    /// Rational Gregory patch (§5.3 N-sided corner resolution).
    Gregory,
}

/// The kind of parametric curve backing an edge.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum CurveType {
    /// Straight line.
    Line,
    /// Full circle.
    Circle,
    /// Circular arc.
    Arc,
    /// Elliptical arc.
    Ellipse,
    /// NURBS / B-spline curve.
    Nurbs,
}

// ---------------------------------------------------------------------------
// Per-entity criteria
// ---------------------------------------------------------------------------

/// Geometric criteria for a face.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct FaceCriteria {
    /// Type of the backing surface.
    pub surface_type: SurfaceType,
    /// Centroid of the boundary vertices.
    pub centroid: Point3,
    /// Surface normal at the centroid (unit length).
    pub normal: Vec3,
    /// Approximate area (shoelace on boundary polygon).
    pub area: f64,
}

/// Geometric criteria for an edge.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct EdgeCriteria {
    /// Type of the backing curve.
    pub curve_type: CurveType,
    /// Point at the midpoint of the curve parameter range.
    pub midpoint: Point3,
    /// Unit tangent at the midpoint.
    pub tangent: Vec3,
    /// Approximate arc length.
    pub length: f64,
}

/// Geometric criteria for a vertex.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct VertexCriteria {
    /// Position of the vertex.
    pub position: Point3,
}

// ---------------------------------------------------------------------------
// Top-level enum
// ---------------------------------------------------------------------------

/// Geometric criteria for any topological entity.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub enum GeometricCriteria {
    /// Criteria for a face.
    Face(FaceCriteria),
    /// Criteria for an edge.
    Edge(EdgeCriteria),
    /// Criteria for a vertex.
    Vertex(VertexCriteria),
}

// ---------------------------------------------------------------------------
// Capture helpers
// ---------------------------------------------------------------------------

/// Capture geometric criteria for a face.
///
/// Returns `None` if `face` is invalid or the face has no surface geometry.
pub fn capture_face(brep: &BRep, face: FaceHandle) -> Option<FaceCriteria> {
    use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};
    use crate::geometry::uv_projection::project_point_to_uv;

    let face_data = brep.faces.get(face)?;
    let surface_on_face = face_data.surface.as_ref()?;
    let geom = &surface_on_face.geometry;

    let surface_type = match geom {
        SurfaceGeometry::Plane { .. } => SurfaceType::Plane,
        SurfaceGeometry::Cylinder { .. } => SurfaceType::Cylinder,
        SurfaceGeometry::Cone { .. } => SurfaceType::Cone,
        SurfaceGeometry::Sphere { .. } => SurfaceType::Sphere,
        SurfaceGeometry::Torus { .. } => SurfaceType::Torus,
        SurfaceGeometry::Nurbs(_) => SurfaceType::Nurbs,
        SurfaceGeometry::Gregory(_) => SurfaceType::Gregory,
    };

    // Boundary vertex positions.
    let positions = brep.face_outer_loop_positions(face);
    if positions.is_empty() {
        return None;
    }

    // Centroid = average of boundary positions.
    let n = positions.len() as f64;
    let cx = positions.iter().map(|p| p.x).sum::<f64>() / n;
    let cy = positions.iter().map(|p| p.y).sum::<f64>() / n;
    let cz = positions.iter().map(|p| p.z).sum::<f64>() / n;
    let centroid = Point3::new(cx, cy, cz);

    // Normal at the centroid's UV projection.
    let normal = if let SurfaceGeometry::Plane { normal, .. } = geom {
        // Planes have a constant normal; avoid the projection.
        *normal
    } else if let Some((u, v)) = project_point_to_uv(geom, &centroid) {
        let raw = geom.normal_at(u, v);
        raw.normalized().unwrap_or(Vec3::new(0.0, 0.0, 1.0))
    } else {
        // Fallback: use the normal at the surface origin.
        let raw = geom.normal_at(0.0, 0.0);
        raw.normalized().unwrap_or(Vec3::new(0.0, 0.0, 1.0))
    };

    // Flip if surface orientation says reversed.
    let normal = if surface_on_face.orientation {
        normal
    } else {
        Vec3::new(-normal.x, -normal.y, -normal.z)
    };

    // Approximate area: shoelace formula on the boundary polygon projected
    // onto the face plane.  We use two tangent basis vectors derived from
    // the actual (unit) normal.
    let area = polygon_area(&positions, &normal);

    Some(FaceCriteria {
        surface_type,
        centroid,
        normal,
        area,
    })
}

/// Capture geometric criteria for an edge.
///
/// Returns `None` if `edge` is invalid or has no curve geometry.
pub fn capture_edge(brep: &BRep, edge: EdgeHandle) -> Option<EdgeCriteria> {
    use crate::geometry::curve::{CurveEval, CurveGeometry};

    let edge_data = brep.edges.get(edge)?;
    let curve_on_edge = edge_data.curve.as_ref()?;
    let geom = &curve_on_edge.geometry;

    let curve_type = match geom {
        CurveGeometry::Line { .. } => CurveType::Line,
        CurveGeometry::Circle { .. } => CurveType::Circle,
        CurveGeometry::Arc { .. } => CurveType::Arc,
        CurveGeometry::Ellipse { .. } => CurveType::Ellipse,
        CurveGeometry::Nurbs(_) => CurveType::Nurbs,
    };

    let t0 = curve_on_edge.t_start;
    let t1 = curve_on_edge.t_end;
    let t_mid = (t0 + t1) * 0.5;

    let midpoint = geom.point_at(t_mid);
    let tangent = geom.tangent_at(t_mid);
    let length = geom.length_between(t0, t1);

    Some(EdgeCriteria {
        curve_type,
        midpoint,
        tangent,
        length,
    })
}

/// Capture geometric criteria for a vertex.
///
/// Returns `None` if `vertex` is invalid.
pub fn capture_vertex(brep: &BRep, vertex: VertexHandle) -> Option<VertexCriteria> {
    let vertex_data = brep.vertices.get(vertex)?;
    Some(VertexCriteria {
        position: vertex_data.position,
    })
}

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

/// Compute the unsigned area of a planar polygon using the 3D shoelace formula.
///
/// The polygon vertices are projected onto the plane defined by `normal`; the
/// result is always non-negative (absolute value of the cross-product sum).
fn polygon_area(positions: &[Point3], normal: &Vec3) -> f64 {
    if positions.len() < 3 {
        return 0.0;
    }
    // 3D shoelace: area = |Σ (P_i × P_{i+1}) · normal| / 2
    let mut cross_sum = Vec3::zero();
    let n = positions.len();
    for i in 0..n {
        let p0 = positions[i];
        let p1 = positions[(i + 1) % n];
        let v0 = Vec3::new(p0.x, p0.y, p0.z);
        let v1 = Vec3::new(p1.x, p1.y, p1.z);
        cross_sum = Vec3::new(
            cross_sum.x + v0.y * v1.z - v0.z * v1.y,
            cross_sum.y + v0.z * v1.x - v0.x * v1.z,
            cross_sum.z + v0.x * v1.y - v0.y * v1.x,
        );
    }
    cross_sum.dot(normal).abs() * 0.5
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tolerance::Tolerance;
    use crate::topology::BRep;

    #[test]
    fn capture_vertex_correct_position() {
        let brep = BRep::unit_cube(Tolerance::default()).unwrap();
        let vertex_handle = brep.vertices.keys().next().expect("at least one vertex");
        let crit = capture_vertex(&brep, vertex_handle).expect("capture_vertex returned None");
        let stored = brep.vertices[vertex_handle].position;
        assert_eq!(crit.position.x, stored.x);
        assert_eq!(crit.position.y, stored.y);
        assert_eq!(crit.position.z, stored.z);
    }

    #[test]
    fn capture_vertex_invalid_handle_returns_none() {
        let brep = BRep::unit_cube(Tolerance::default()).unwrap();
        let bogus = crate::topology::VertexHandle::default();
        assert!(capture_vertex(&brep, bogus).is_none());
    }

    #[test]
    fn capture_edge_line_on_cube() {
        let brep = BRep::unit_cube(Tolerance::default()).unwrap();
        // Find any edge that has curve geometry.
        let edge_handle = brep
            .edges
            .iter()
            .find(|(_, e)| e.curve.is_some())
            .map(|(h, _)| h)
            .expect("unit cube should have edges with curves");
        let crit = capture_edge(&brep, edge_handle).expect("capture_edge returned None");
        assert_eq!(crit.curve_type, CurveType::Line);
        // Length of a unit-cube edge should be 1.0.
        assert!((crit.length - 1.0).abs() < 1e-9, "length={}", crit.length);
    }

    #[test]
    fn capture_edge_invalid_handle_returns_none() {
        let brep = BRep::unit_cube(Tolerance::default()).unwrap();
        let bogus = crate::topology::EdgeHandle::default();
        assert!(capture_edge(&brep, bogus).is_none());
    }

    #[test]
    fn capture_face_plane_on_cube() {
        let brep = BRep::unit_cube(Tolerance::default()).unwrap();
        let face_handle = brep.faces.keys().next().expect("at least one face");
        let crit = capture_face(&brep, face_handle).expect("capture_face returned None");

        assert_eq!(crit.surface_type, SurfaceType::Plane);

        // Centroid of a unit-cube face should be inside [0,1]^3.
        assert!(
            crit.centroid.x >= -1e-9 && crit.centroid.x <= 1.0 + 1e-9,
            "centroid.x={}",
            crit.centroid.x
        );
        assert!(
            crit.centroid.y >= -1e-9 && crit.centroid.y <= 1.0 + 1e-9,
            "centroid.y={}",
            crit.centroid.y
        );
        assert!(
            crit.centroid.z >= -1e-9 && crit.centroid.z <= 1.0 + 1e-9,
            "centroid.z={}",
            crit.centroid.z
        );

        // Normal should be a unit vector (one of ±X, ±Y, ±Z).
        let len = (crit.normal.x * crit.normal.x
            + crit.normal.y * crit.normal.y
            + crit.normal.z * crit.normal.z)
            .sqrt();
        assert!((len - 1.0).abs() < 1e-9, "normal not unit: len={}", len);

        // Area of a unit face = 1.0.
        assert!((crit.area - 1.0).abs() < 1e-9, "area={}", crit.area);
    }

    #[test]
    fn capture_face_no_surface_returns_none() {
        use crate::topology::FaceHandle;
        // An invalid face handle returns None.
        let brep = BRep::new(Tolerance::default());
        let bogus = FaceHandle::default();
        assert!(capture_face(&brep, bogus).is_none());
    }

    #[test]
    fn capture_edge_midpoint_on_unit_cube_edge() {
        let brep = BRep::unit_cube(Tolerance::default()).unwrap();
        let edge_handle = brep
            .edges
            .iter()
            .find(|(_, e)| e.curve.is_some())
            .map(|(h, _)| h)
            .expect("unit cube should have edges with curves");
        let crit = capture_edge(&brep, edge_handle).unwrap();
        // Midpoint should be within the unit cube bounds.
        assert!(crit.midpoint.x >= -1e-9 && crit.midpoint.x <= 1.0 + 1e-9);
        assert!(crit.midpoint.y >= -1e-9 && crit.midpoint.y <= 1.0 + 1e-9);
        assert!(crit.midpoint.z >= -1e-9 && crit.midpoint.z <= 1.0 + 1e-9);
    }

    #[test]
    fn capture_face_all_cube_faces_have_unit_area() {
        let brep = BRep::unit_cube(Tolerance::default()).unwrap();
        for (handle, _) in &brep.faces {
            let crit = capture_face(&brep, handle).expect("all cube faces have surfaces");
            assert!(
                (crit.area - 1.0).abs() < 1e-9,
                "face area={} expected 1.0",
                crit.area
            );
        }
    }
}
