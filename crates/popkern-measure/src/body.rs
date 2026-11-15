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
/// Fine-tessellation chord as a fraction of the body's vertex-box
/// diagonal (PK-843), capped at the historical 0.001. An absolute
/// 0.001 was the whole of a millimetre part (a drilled dome's volume
/// lost 2.4e-7 at k = 1e-3); 5e-4 of the diagonal reproduces 0.001 on
/// the unit-scale fixtures. Large parts keep 0.001 (a kilometre-scale
/// dome meshes 28x slower) until the quadrature is chord-independent.
const FINE_CHORD_PER_DIAGONAL: f64 = 5.0e-4;

pub(crate) fn tessellate_fine(brep: &BRep) -> Result<popkern_mesh::TriMesh, MeasureError> {
    let (lo, hi) = brep.vertices.values().fold(
        (
            Point3::new(f64::MAX, f64::MAX, f64::MAX),
            Point3::new(f64::MIN, f64::MIN, f64::MIN),
        ),
        |(lo, hi), v| {
            let p = v.position;
            (
                Point3::new(lo.x.min(p.x), lo.y.min(p.y), lo.z.min(p.z)),
                Point3::new(hi.x.max(p.x), hi.y.max(p.y), hi.z.max(p.z)),
            )
        },
    );
    let diagonal = if lo.x <= hi.x {
        lo.distance_to(&hi)
    } else {
        0.0
    };
    // Capped at the historical 0.001: the quadrature's trim slivers
    // are not chord-independent yet (PK-818), and a 200 mm tee read
    // 1.3e-4 of its volume at the coarser relative chord.
    let cap = 0.001 * (diagonal / crate::body_origin::TUNED_EXTENT).max(1.0);
    let chord = if diagonal > 0.0 {
        (FINE_CHORD_PER_DIAGONAL * diagonal).min(cap)
    } else {
        0.001
    };
    let opts = popkern_mesh::TessellationOptions {
        chord_height: chord,
        min_subdivisions: 8,
        max_subdivisions: 128,
        ..Default::default()
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
///
/// Gregory corner caps (§5.3 item #9) are integrated on their exact
/// surface like every other curved face, oriented by the mesh.
///
/// A general body's volume is the sum of its solid regions'
/// ([`crate::region_volumes`]).
pub fn body_volume(brep: &BRep, tol: &Tolerance) -> Result<f64, MeasureError> {
    let kind = brep.require_kind(
        "volume",
        &[
            popkern_core::topology::BodyKind::Solid,
            popkern_core::topology::BodyKind::General,
            popkern_core::topology::BodyKind::Empty,
        ],
    )?;
    // The divergence terms r · n are taken about the origin; for a body
    // far from it they are large and cancel to the small volume, and
    // the quadrature's error on them swamps it (a pole cap 1e-4 high at
    // z = 5 read 4.6 times its volume, PK-1135). About a point inside
    // the body's box the sum is the same (∮ c · n dA = 0 on a closed
    // boundary) without the cancellation: measure the body moved there.
    if let Some(offset) = crate::body_origin::far_from_origin(brep) {
        let moved = brep
            .transform(&popkern_core::geometry::Transform3::translation(
                -offset.x, -offset.y, -offset.z,
            ))
            .map_err(|e| MeasureError::DegenerateGeometry(e.to_string()))?;
        return body_volume(&moved, tol);
    }
    if kind == popkern_core::topology::BodyKind::General {
        return Ok(crate::region_volumes(brep, tol)?.iter().map(|r| r.1).sum());
    }
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

    // Curved bodies, Gregory corner caps included: the fine
    // tessellation's triangles integrated on the exact surface, each
    // face's term signed by the mesh's outward winding (see
    // `body_quadrature`). A separate Gregory path signed each face's term
    // by its `orientation` flag, which is a winding contract, not a
    // material side: the G2 boss plate read 1.0043 for a mesh volume of
    // 1.7593 (PK-1193).
    let mesh = tessellate_fine(brep)?;
    let gate = 100.0 * brep.tolerance.linear().max(1e-9);
    Ok(crate::body_quadrature::quadrature_volume(brep, &mesh, gate).abs())
}

/// Exact volume of a body bounded only by planar faces, via the
/// divergence theorem: each face contributes `⅓ (o · A)` from its
/// exact vector area `A` (edge curves integrated, so arcs and
/// circular holes count) and any vertex `o`. Inner loops wind
/// opposite to the outer one, so holes subtract themselves.
fn exact_polyhedral_volume(brep: &BRep) -> Result<f64, MeasureError> {
    let volume: f64 = brep
        .faces
        .keys()
        .map(|face| crate::planar_area::planar_face_divergence(brep, face))
        .sum();
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
    brep.require_kind(
        "centroid",
        &[
            popkern_core::topology::BodyKind::Solid,
            popkern_core::topology::BodyKind::Empty,
        ],
    )?;
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
/// The solid is treated as having uniform unit density. The integral
/// `I = ∫ r² dV` is computed EXACTLY per signed tetrahedron (apex at the
/// body centroid, base at each surface triangle) via the closed-form
/// second-moment integral: for a tetrahedron `(0, A, B, C)` with signed
/// volume `V`,
///
/// ```text
/// ∫ x² dV = (V/10)·(Ax² + Bx² + Cx² + Ax·Bx + Ax·Cx + Bx·Cx)
/// ```
///
/// and likewise for y and z (barycentric monomial integrals
/// `∫u² = 1/60`, `∫uv = 1/120` over the unit simplex, scaled by `6V`).
/// A point-mass-at-tet-centroid approximation is NOT used: it drops each
/// tetrahedron's own second moment and undercounts a unit cube by
/// exactly 2× (1/12 instead of 1/6).
pub fn body_moment_of_inertia(
    brep: &BRep,
    axis: Axis,
    tol: &Tolerance,
) -> Result<f64, MeasureError> {
    brep.require_kind(
        "moment of inertia",
        &[
            popkern_core::topology::BodyKind::Solid,
            popkern_core::topology::BodyKind::Empty,
        ],
    )?;
    let mesh = tessellate_fine(brep)?;
    let centroid = body_centroid(brep, tol)?;

    // Accumulate the exact diagonal second moments ∫x²dV, ∫y²dV, ∫z²dV
    // with the body centroid as the coordinate origin (so the tet apex
    // is the origin and the closed-form applies directly).
    let mut sxx = 0.0_f64;
    let mut syy = 0.0_f64;
    let mut szz = 0.0_f64;

    for tri in &mesh.triangles {
        let a = to_vec3(mesh.vertices[tri[0]].position) - to_vec3(centroid);
        let b = to_vec3(mesh.vertices[tri[1]].position) - to_vec3(centroid);
        let c = to_vec3(mesh.vertices[tri[2]].position) - to_vec3(centroid);

        let sv = a.dot(&b.cross(&c)) / 6.0;

        let second = |pa: f64, pb: f64, pc: f64| -> f64 {
            (sv / 10.0) * (pa * pa + pb * pb + pc * pc + pa * pb + pa * pc + pb * pc)
        };
        sxx += second(a.x, b.x, c.x);
        syy += second(a.y, b.y, c.y);
        szz += second(a.z, b.z, c.z);
    }

    // A consistently inward-wound mesh flips the sign of every signed
    // volume uniformly; the axis moment is the magnitude.
    let moi = match axis {
        Axis::X => syy + szz,
        Axis::Y => sxx + szz,
        Axis::Z => sxx + syy,
    };

    Ok(moi.abs())
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
#[path = "body_tests.rs"]
mod tests;
