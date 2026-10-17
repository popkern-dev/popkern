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

fn tessellate_fine(brep: &BRep) -> Result<popkern_mesh::TriMesh, MeasureError> {
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
    let chord = if diagonal > 0.0 {
        (FINE_CHORD_PER_DIAGONAL * diagonal).min(0.001)
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
/// §5.3 item #9: when the body contains `SurfaceGeometry::Gregory` faces,
/// those faces contribute via native 8×8 Gauss-Legendre quadrature of
/// `(1/3) S · (S_u × S_v) du dv` over the UV rectangle (orders of
/// magnitude tighter than the chord-height tessellation error).
/// Non-Gregory non-planar faces still go through tessellation.
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

    let has_gregory = brep.faces.values().any(|f| {
        f.surface
            .as_ref()
            .is_some_and(|s| matches!(s.geometry, SurfaceGeometry::Gregory(_)))
    });
    if has_gregory {
        return native_divergence_volume(brep);
    }

    // Curved bodies: the fine tessellation's triangles integrated on
    // the exact surface (see `body_quadrature`).
    let mesh = tessellate_fine(brep)?;
    let gate = 100.0 * brep.tolerance.linear().max(1e-9);
    Ok(crate::body_quadrature::quadrature_volume(brep, &mesh, gate).abs())
}

/// Volume via divergence theorem, computing each face's contribution
/// natively by surface type. Planar faces use exact polygon
/// triangulation; Gregory faces use 8×8 Gauss-Legendre quadrature of
/// `(1/3) · S · (S_u × S_v)`; other curved faces fall back to a
/// per-face fine tessellation. Used only when `body_volume` sees at
/// least one Gregory face.
fn native_divergence_volume(brep: &BRep) -> Result<f64, MeasureError> {
    let mut total = 0.0_f64;
    for (face_h, face) in &brep.faces {
        let Some(sof) = face.surface.as_ref() else {
            continue;
        };
        let orient_sign = if sof.orientation { 1.0 } else { -1.0 };
        let contrib = match &sof.geometry {
            SurfaceGeometry::Plane { .. } => {
                crate::planar_area::planar_face_divergence(brep, face_h)
            }
            SurfaceGeometry::Gregory(_) if face.inner_loops.is_empty() => {
                match crate::face_helpers::pcurve_uv_range(brep, face_h) {
                    Some(((u0, u1), (v0, v1)))
                        if crate::face_helpers::outer_loop_matches_uv_rectangle(
                            brep,
                            face_h,
                            ((u0, u1), (v0, v1)),
                        ) =>
                    {
                        gauss_legendre_face_volume(&sof.geometry, (u0, u1), (v0, v1))
                    }
                    _ => tessellated_face_divergence_contribution(brep, face_h)?,
                }
            }
            _ => tessellated_face_divergence_contribution(brep, face_h)?,
        };
        total += orient_sign * contrib;
    }
    Ok(total.abs())
}

/// Divergence contribution from a single curved face via fine
/// tessellation. Used as the fallback branch in
/// `native_divergence_volume` for non-Gregory non-planar faces.
fn tessellated_face_divergence_contribution(
    brep: &BRep,
    face: popkern_core::topology::FaceHandle,
) -> Result<f64, MeasureError> {
    let opts = popkern_mesh::TessellationOptions {
        chord_height: 0.001,
        min_subdivisions: 8,
        max_subdivisions: 128,
        ..Default::default()
    };
    let face_mesh = popkern_mesh::tessellate_face(brep, face, &opts)
        .map_err(|e| MeasureError::TessellationFailed(e.to_string()))?;
    let mut v = 0.0_f64;
    for tri in &face_mesh.triangles {
        let v0 = to_vec3(face_mesh.vertices[tri[0]].position);
        let v1 = to_vec3(face_mesh.vertices[tri[1]].position);
        let v2 = to_vec3(face_mesh.vertices[tri[2]].position);
        v += signed_tet_volume(v0, v1, v2);
    }
    Ok(v)
}

/// 8×8 tensor-product Gauss-Legendre quadrature of
/// `(1/3) · S(u,v) · (S_u × S_v)` over the UV rectangle — the
/// divergence-theorem volume contribution from a single face.
/// `1/3` comes from the factor in `V = (1/3) ∮ r · n dA`.
fn gauss_legendre_face_volume(
    geom: &SurfaceGeometry,
    (u0, u1): (f64, f64),
    (v0, v1): (f64, f64),
) -> f64 {
    use popkern_core::geometry::surface::SurfaceEval;
    const GL_X: [f64; 8] = [
        -0.960_289_856_497_536_2,
        -0.796_666_477_413_626_7,
        -0.525_532_409_916_329,
        -0.183_434_642_495_649_8,
        0.183_434_642_495_649_8,
        0.525_532_409_916_329,
        0.796_666_477_413_626_7,
        0.960_289_856_497_536_2,
    ];
    const GL_W: [f64; 8] = [
        0.101_228_536_290_376_26,
        0.222_381_034_453_374_4,
        0.313_706_645_877_887_3,
        0.362_683_783_378_362,
        0.362_683_783_378_362,
        0.313_706_645_877_887_3,
        0.222_381_034_453_374_4,
        0.101_228_536_290_376_26,
    ];
    let du = 0.5 * (u1 - u0);
    let dv = 0.5 * (v1 - v0);
    let mu = 0.5 * (u1 + u0);
    let mv = 0.5 * (v1 + v0);
    let mut sum = 0.0_f64;
    for i in 0..8 {
        let u = mu + du * GL_X[i];
        for j in 0..8 {
            let v = mv + dv * GL_X[j];
            let p = geom.point_at(u, v);
            let su = geom.partial_u(u, v);
            let sv = geom.partial_v(u, v);
            let n = su.cross(&sv);
            // r · n where r = (p.x, p.y, p.z).
            let r_dot_n = p.x * n.x + p.y * n.y + p.z * n.z;
            sum += GL_W[i] * GL_W[j] * r_dot_n;
        }
    }
    du * dv * sum / 3.0
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
