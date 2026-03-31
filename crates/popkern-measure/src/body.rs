//! Body-level measurement functions (volume, surface area, centroid, bounding box,
//! moment of inertia).

use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use crate::MeasureError;

// ---------------------------------------------------------------------------
// Supporting types
// ---------------------------------------------------------------------------

/// Axis selection for moment-of-inertia queries.
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub enum Axis {
    /// The X axis.
    X,
    /// The Y axis.
    Y,
    /// The Z axis.
    Z,
}

/// Axis-aligned bounding box.
#[derive(Debug, Clone, Copy)]
pub struct BoundingBox {
    /// Minimum corner of the bounding box.
    pub min: Point3,
    /// Maximum corner of the bounding box.
    pub max: Point3,
}

impl BoundingBox {
    /// Extent of the bounding box along the X axis.
    pub fn size_x(&self) -> f64 {
        self.max.x - self.min.x
    }

    /// Extent of the bounding box along the Y axis.
    pub fn size_y(&self) -> f64 {
        self.max.y - self.min.y
    }

    /// Extent of the bounding box along the Z axis.
    pub fn size_z(&self) -> f64 {
        self.max.z - self.min.z
    }
}

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

/// Tessellate the BRep with fine settings suitable for volumetric queries.
fn tessellate_fine(brep: &BRep) -> Result<popkern_mesh::TriMesh, MeasureError> {
    let opts = popkern_mesh::TessellationOptions {
        chord_height: 0.001,
        min_subdivisions: 8,
        max_subdivisions: 128,
    };
    popkern_mesh::tessellate(brep, &opts)
        .map_err(|e| MeasureError::TessellationFailed(e.to_string()))
}

/// Signed volume contribution of the tetrahedron formed by the origin and triangle (v0,v1,v2).
///
/// Based on the divergence theorem: V_signed = v0 · (v1 × v2) / 6
#[inline]
fn signed_tet_volume(v0: Vec3, v1: Vec3, v2: Vec3) -> f64 {
    v0.dot(&v1.cross(&v2)) / 6.0
}

/// Convert a `Point3` to a `Vec3` (position vector from origin).
#[inline]
fn to_vec3(p: Point3) -> Vec3 {
    Vec3::new(p.x, p.y, p.z)
}

// ---------------------------------------------------------------------------
// Public measurement functions
// ---------------------------------------------------------------------------

/// Compute the volume of a closed BRep solid using the divergence theorem.
///
/// If every face is planar the volume is computed exactly from the polygon
/// faces without tessellating. Otherwise the solid is tessellated at high
/// resolution and the signed tetrahedral volume formula
/// `V = (1/6) * Σ v0 · (v1 × v2)` is summed over all triangles.
pub fn body_volume(brep: &BRep, _tol: &Tolerance) -> Result<f64, MeasureError> {
    // Fast path: all-planar bodies (boxes, prisms, polyhedra) are exact.
    let all_planar = brep.faces.iter().all(|(_, face)| {
        face.surface
            .as_ref()
            .map(|s| matches!(s.geometry, SurfaceGeometry::Plane { .. }))
            .unwrap_or(false)
    });

    if all_planar {
        return exact_polyhedral_volume(brep);
    }

    // Fallback: tessellate and apply divergence theorem on triangles.
    let mesh = tessellate_fine(brep)?;

    let volume: f64 = mesh
        .triangles
        .iter()
        .map(|tri| {
            let v0 = to_vec3(mesh.vertices[tri[0]].position);
            let v1 = to_vec3(mesh.vertices[tri[1]].position);
            let v2 = to_vec3(mesh.vertices[tri[2]].position);
            signed_tet_volume(v0, v1, v2)
        })
        .sum::<f64>()
        .abs();

    Ok(volume)
}

/// Exact volume for a polyhedron with only planar faces via the divergence theorem.
///
/// Each face polygon is fan-triangulated from its first vertex and the signed
/// tetrahedral contribution `v0 · (v1 × v2) / 6` is accumulated. Inner loops
/// (holes) are handled by the same accumulation — their opposite winding
/// automatically subtracts their area from the enclosing face contribution.
fn exact_polyhedral_volume(brep: &BRep) -> Result<f64, MeasureError> {
    let mut volume = 0.0_f64;

    for (_, face) in &brep.faces {
        let positions = crate::face_helpers::loop_positions(brep, face.outer_loop);
        if positions.len() < 3 {
            continue;
        }
        // Fan triangulation from the first vertex.
        let v0 = to_vec3(positions[0]);
        for i in 1..positions.len() - 1 {
            let v1 = to_vec3(positions[i]);
            let v2 = to_vec3(positions[i + 1]);
            volume += signed_tet_volume(v0, v1, v2);
        }

        // Inner loops (holes): their opposite winding already gives the correct
        // sign, so we just accumulate their contributions directly.
        for &inner_loop in &face.inner_loops {
            let hole_positions = crate::face_helpers::loop_positions(brep, inner_loop);
            if hole_positions.len() < 3 {
                continue;
            }
            let h0 = to_vec3(hole_positions[0]);
            for i in 1..hole_positions.len() - 1 {
                let h1 = to_vec3(hole_positions[i]);
                let h2 = to_vec3(hole_positions[i + 1]);
                volume += signed_tet_volume(h0, h1, h2);
            }
        }
    }

    Ok(volume.abs())
}

/// Compute the total surface area of a BRep solid.
///
/// Sums `face_area` for every face in the solid.
pub fn body_surface_area(brep: &BRep, tol: &Tolerance) -> Result<f64, MeasureError> {
    let mut total = 0.0_f64;
    for (fh, _) in &brep.faces {
        total += crate::face::face_area(brep, fh, tol)?;
    }
    Ok(total)
}

/// Compute the volume centroid of a closed BRep solid.
///
/// Uses the divergence theorem to compute the first moments of volume, then
/// divides by total volume to obtain the centroid coordinates.
pub fn body_centroid(brep: &BRep, _tol: &Tolerance) -> Result<Point3, MeasureError> {
    let mesh = tessellate_fine(brep)?;

    let mut vol = 0.0_f64;
    let mut cx = 0.0_f64;
    let mut cy = 0.0_f64;
    let mut cz = 0.0_f64;

    for tri in &mesh.triangles {
        let v0 = to_vec3(mesh.vertices[tri[0]].position);
        let v1 = to_vec3(mesh.vertices[tri[1]].position);
        let v2 = to_vec3(mesh.vertices[tri[2]].position);

        let sv = signed_tet_volume(v0, v1, v2);
        vol += sv;
        // First moment: each coordinate's contribution for the signed tetrahedron
        cx += sv * (v0.x + v1.x + v2.x) / 4.0;
        cy += sv * (v0.y + v1.y + v2.y) / 4.0;
        cz += sv * (v0.z + v1.z + v2.z) / 4.0;
    }

    if vol.abs() < f64::EPSILON {
        return Err(MeasureError::DegenerateGeometry(
            "body has zero computed volume".into(),
        ));
    }

    Ok(Point3::new(cx / vol, cy / vol, cz / vol))
}

/// Compute the axis-aligned bounding box of a BRep solid.
///
/// Iterates all vertex positions and records the per-axis min/max values.
pub fn body_bounding_box(brep: &BRep) -> Result<BoundingBox, MeasureError> {
    let mut vertices = brep.vertices.values();

    let first = vertices
        .next()
        .ok_or_else(|| MeasureError::DegenerateGeometry("BRep has no vertices".into()))?;

    let p0 = first.position;
    let mut min_x = p0.x;
    let mut min_y = p0.y;
    let mut min_z = p0.z;
    let mut max_x = p0.x;
    let mut max_y = p0.y;
    let mut max_z = p0.z;

    for v in vertices {
        let p = v.position;
        if p.x < min_x {
            min_x = p.x;
        }
        if p.y < min_y {
            min_y = p.y;
        }
        if p.z < min_z {
            min_z = p.z;
        }
        if p.x > max_x {
            max_x = p.x;
        }
        if p.y > max_y {
            max_y = p.y;
        }
        if p.z > max_z {
            max_z = p.z;
        }
    }

    Ok(BoundingBox {
        min: Point3::new(min_x, min_y, min_z),
        max: Point3::new(max_x, max_y, max_z),
    })
}

/// Compute the second moment of mass (moment of inertia) of a BRep solid about
/// the given axis through the body's centroid.
///
/// The solid is treated as having uniform unit density. The computation
/// approximates the integral `I = ∫ r² dV` by summing over signed tetrahedra,
/// where `r²` is the squared perpendicular distance from each tetrahedron's
/// centroid to the chosen axis through the body centroid.
pub fn body_moment_of_inertia(
    brep: &BRep,
    axis: Axis,
    tol: &Tolerance,
) -> Result<f64, MeasureError> {
    let mesh = tessellate_fine(brep)?;
    let centroid = body_centroid(brep, tol)?;

    let mut moi = 0.0_f64;

    for tri in &mesh.triangles {
        let v0 = to_vec3(mesh.vertices[tri[0]].position);
        let v1 = to_vec3(mesh.vertices[tri[1]].position);
        let v2 = to_vec3(mesh.vertices[tri[2]].position);

        let sv = signed_tet_volume(v0, v1, v2);

        // Centroid of the tetrahedron relative to the body centroid.
        let tx = (v0.x + v1.x + v2.x) / 4.0 - centroid.x;
        let ty = (v0.y + v1.y + v2.y) / 4.0 - centroid.y;
        let tz = (v0.z + v1.z + v2.z) / 4.0 - centroid.z;

        // Squared perpendicular distance to the axis through the centroid.
        let r2 = match axis {
            Axis::X => ty * ty + tz * tz,
            Axis::Y => tx * tx + tz * tz,
            Axis::Z => tx * tx + ty * ty,
        };

        moi += r2 * sv;
    }

    Ok(moi.abs())
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use popkern_core::tolerance::Tolerance;
    use popkern_core::topology::BRep;

    fn make_cube() -> BRep {
        BRep::unit_cube(Tolerance::default()).expect("cube")
    }

    #[test]
    fn unit_cube_volume() {
        let cube = make_cube();
        let tol = Tolerance::default();
        let vol = body_volume(&cube, &tol).unwrap();
        assert!((vol - 1.0).abs() < 0.01, "expected ~1.0, got {vol}");
    }

    #[test]
    fn unit_cube_surface_area() {
        let cube = make_cube();
        let tol = Tolerance::default();
        let sa = body_surface_area(&cube, &tol).unwrap();
        assert!((sa - 6.0).abs() < 0.01, "expected ~6.0, got {sa}");
    }

    #[test]
    fn unit_cube_centroid() {
        let cube = make_cube();
        let tol = Tolerance::default();
        let c = body_centroid(&cube, &tol).unwrap();
        // Unit cube is [0,1]^3, centroid should be (0.5, 0.5, 0.5)
        assert!((c.x - 0.5).abs() < 0.01, "centroid x={}", c.x);
        assert!((c.y - 0.5).abs() < 0.01, "centroid y={}", c.y);
        assert!((c.z - 0.5).abs() < 0.01, "centroid z={}", c.z);
    }

    #[test]
    fn unit_cube_bounding_box() {
        let cube = make_cube();
        let bb = body_bounding_box(&cube).unwrap();
        assert!((bb.min.x - 0.0).abs() < 1e-10);
        assert!((bb.max.x - 1.0).abs() < 1e-10);
        assert!((bb.size_x() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn unit_cube_moment_of_inertia() {
        let cube = make_cube();
        let tol = Tolerance::default();
        // Moment of inertia of unit cube (side=1, mass=volume=1) about an axis through centroid:
        // I = (1/6) for a unit cube ≈ 0.1667
        let moi = body_moment_of_inertia(&cube, Axis::Z, &tol).unwrap();
        assert!((moi - 1.0 / 6.0).abs() < 0.05, "expected ~0.167, got {moi}");
    }
}
