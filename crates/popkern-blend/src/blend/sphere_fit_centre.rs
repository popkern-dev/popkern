//! Cap-centre selection helpers for the sphere-fit boundary
//! pipeline.
//!
//! Split out of [`super::sphere_fit_boundary`] to keep that file
//! under the 500-line structural cap.  Contains:
//!
//! * [`face_outward_normal_at`] — robust outward-normal lookup
//!   that resolves cube/pyramid sign-convention drift empirically.
//! * [`fillet_nominal_radius`] — extracts radius from Cylinder /
//!   Torus surfaces.
//! * [`compute_cap_centre`] — solves the inscribed-sphere system
//!   over N adjacent face planes.
//! * [`per_fillet_rolling_ball_centre`] — falls back to the
//!   cylinder-axis-tangent centre.
//! * [`select_per_fillet_centre`] — orchestrates the choice
//!   between the inscribed and rolling-ball solutions with a
//!   side-check against the apex's interior direction.

use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::uv_projection::project_point_to_uv;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::topology::{BRep, FaceHandle};

/// Outward face normal at `point` (used here only at the apex).
///
/// Different BRep constructors disagree on the natural-normal sign
/// convention.  `unit_cube` traverses outer loops such that
/// `plane_from_points`' right-hand-rule cross product yields an
/// **inward**-pointing normal; `square_pyramid` traverses the
/// opposite way and yields an **outward** normal.  Trusting
/// `SurfaceOnFace::orientation` alone gives inconsistent results
/// across these constructors.
///
/// To return a consistent outward direction we resolve the
/// stored-normal sign empirically: project the brep's vertex
/// centroid (which approximates the solid's interior for a convex
/// shell) onto the face plane via the candidate outward direction
/// and pick the sign that points AWAY from the centroid.  This is
/// robust across cube, pyramid, n_gonal_pyramid, and any future
/// solid whose loop orientation we don't manually audit.
pub(crate) fn face_outward_normal_at(
    brep: &BRep,
    face: FaceHandle,
    point: &Point3,
) -> Option<Vec3> {
    let face_data = brep.faces.get(face)?;
    let sof = face_data.surface.as_ref()?;
    let (u, v) = project_point_to_uv(&sof.geometry, point)?;
    let n = sof.geometry.normal_at(u, v).normalized()?;
    // Centroid of all brep vertices — for a closed convex solid
    // this lies in the interior; the outward normal at any
    // boundary point therefore points AWAY from the centroid.
    let mut sum = Vec3::new(0.0, 0.0, 0.0);
    let mut count = 0_usize;
    for (_, v) in &brep.vertices {
        sum.x += v.position.x;
        sum.y += v.position.y;
        sum.z += v.position.z;
        count += 1;
    }
    if count == 0 {
        return None;
    }
    let inv = 1.0 / (count as f64);
    let centroid = Point3::new(sum.x * inv, sum.y * inv, sum.z * inv);
    let to_outside = Vec3::new(
        point.x - centroid.x,
        point.y - centroid.y,
        point.z - centroid.z,
    );
    if n.dot(&to_outside) >= 0.0 {
        Some(n)
    } else {
        Some(-n)
    }
}

/// Extract the nominal fillet radius from a fillet's surface.
pub(crate) fn fillet_nominal_radius(surface: &SurfaceGeometry) -> Option<f64> {
    match surface {
        SurfaceGeometry::Cylinder { radius, .. } => Some(*radius),
        SurfaceGeometry::Torus { minor_radius, .. } => Some(*minor_radius),
        _ => None,
    }
}

/// Solve `M (cap - apex) = -r * 1` for `cap`, where `M` has the
/// outward face normals as rows.  For 3 faces this is a direct 3x3
/// solve; for 4+ faces this is the normal-equations least-squares.
pub(super) fn compute_cap_centre(apex: Point3, face_normals: &[Vec3], r: f64) -> Option<Point3> {
    let k = face_normals.len();
    if k < 3 {
        return None;
    }
    // Build A^T A (3x3) and A^T b (3) where A rows are face normals
    // and b = -r * 1.
    let mut ata = [[0.0f64; 3]; 3];
    let mut atb = [0.0f64; 3];
    for n in face_normals {
        let nx = n.x;
        let ny = n.y;
        let nz = n.z;
        ata[0][0] += nx * nx;
        ata[0][1] += nx * ny;
        ata[0][2] += nx * nz;
        ata[1][0] += ny * nx;
        ata[1][1] += ny * ny;
        ata[1][2] += ny * nz;
        ata[2][0] += nz * nx;
        ata[2][1] += nz * ny;
        ata[2][2] += nz * nz;
        atb[0] += nx * (-r);
        atb[1] += ny * (-r);
        atb[2] += nz * (-r);
    }
    // Solve 3x3 system via Cramer's rule (cap - apex).
    let det = det3(&ata);
    if det.abs() < 1e-18 {
        return None;
    }
    let dx = det3(&replace_col(&ata, 0, atb));
    let dy = det3(&replace_col(&ata, 1, atb));
    let dz = det3(&replace_col(&ata, 2, atb));
    Some(Point3::new(
        apex.x + dx / det,
        apex.y + dy / det,
        apex.z + dz / det,
    ))
}

fn det3(m: &[[f64; 3]; 3]) -> f64 {
    m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
}

fn replace_col(m: &[[f64; 3]; 3], col: usize, b: [f64; 3]) -> [[f64; 3]; 3] {
    let mut out = *m;
    for r in 0..3 {
        out[r][col] = b[r];
    }
    out
}

/// Per-fillet rolling-ball cap centre — the point
/// on the fillet's cylinder axis closest to the apex.  When the
/// fillet's sphere centre lies on its cylinder axis, sphere ∩
/// cylinder is a clean great circle, giving exact analytic
/// G1 between cap and cylinder.
///
/// `compute_cap_centre` solves a different problem (sphere
/// tangent to ALL adjacent face planes — the inscribed sphere).
/// For trihedral cube-like apices the inscribed sphere happens
/// to land on every cylinder axis simultaneously — they coincide.
/// For pyramid-like apices (N ≥ 4) the inscribed sphere is NOT
/// on the cylinder axes, and using it as cap centre breaks the
/// cap-cylinder tangency.  In those cases use this per-fillet
/// formula, which guarantees axis-tangency for THIS fillet (and
/// each per-fillet sub-cap will have its own centre).
pub(crate) fn per_fillet_rolling_ball_centre(
    apex: Point3,
    fillet_cylinder: &SurfaceGeometry,
) -> Option<Point3> {
    let SurfaceGeometry::Cylinder { origin, axis, .. } = fillet_cylinder else {
        return None;
    };
    let d = axis.normalized()?;
    let v = Vec3::new(apex.x - origin.x, apex.y - origin.y, apex.z - origin.z);
    let t = v.dot(&d);
    Some(Point3::new(
        origin.x + d.x * t,
        origin.y + d.y * t,
        origin.z + d.z * t,
    ))
}

/// Pick the per-fillet sphere centre that
/// guarantees clean axis-tangency.  Tries the inscribed-sphere
/// solution from `compute_cap_centre` first — if it lies on the
/// fillet's cylinder axis (within tolerance) AND is on the
/// solid's interior side of the apex, use it.  Otherwise fall
/// back to [`per_fillet_rolling_ball_centre`].
///
/// The inscribed-sphere solution is preferred when valid because
/// it agrees across all per-fillet caps in a trihedral apex
/// (cube), so the `FaceEdge` bridges between adjacent fillets'
/// endpoints collapse to single points (clean shared topology).
/// For pyramid-like N ≥ 4 apices the inscribed sphere is off-axis,
/// so we accept distinct per-fillet centres at the cost of
/// non-degenerate face-edge bridges.
pub(crate) fn select_per_fillet_centre(
    apex: Point3,
    face_normals: &[Vec3],
    r: f64,
    fillet_cylinder: &SurfaceGeometry,
) -> Option<Point3> {
    let inscribed = compute_cap_centre(apex, face_normals, r);
    if popkern_core::trace::trace_enabled("cap-centre") {
        eprintln!(
            "[cap-centre] apex=({:.4},{:.4},{:.4}) r={r:.3} normals.len()={} inscribed={:?}",
            apex.x,
            apex.y,
            apex.z,
            face_normals.len(),
            inscribed,
        );
        for (i, n) in face_normals.iter().enumerate() {
            eprintln!(
                "[cap-centre]   normal[{i}] = ({:.4},{:.4},{:.4})",
                n.x, n.y, n.z
            );
        }
    }
    if let (Some(c), SurfaceGeometry::Cylinder { origin, axis, .. }) = (inscribed, fillet_cylinder)
    {
        // Two acceptance conditions for the inscribed-sphere
        // solution:
        //   (1) inscribed centre lies on the fillet's cylinder
        //       axis (within tolerance) — required for clean
        //       sphere-cylinder tangency along a great circle, AND
        //   (2) inscribed centre is on the INTERIOR side of the
        //       apex relative to the surrounding solid.  For
        //       trihedral cube apices the inscribed centre lands
        //       inside the solid; for pyramid-like N≥4 apices it
        //       can land on the far side of the apex from the
        //       solid — a phantom sphere not rounding off any
        //       material.
        //
        // The interior direction at the apex is computed from the
        // adjacent face normals: `interior_dir = -sum(n_k)`.  Each
        // n_k is an OUTWARD-pointing face normal (made consistent
        // by `face_outward_normal_at`), so the negated sum points
        // into the solid.  This replaces an earlier
        // `apex_to_origin`-based check that depended on the
        // unspecified placement of `cylinder.origin` and broke on
        // cube edges where `compute_fillet_geometry` placed
        // `origin` on the apex side instead of past the far end.
        let Some(d) = axis.normalized() else {
            return Some(c);
        };
        let v = Vec3::new(c.x - origin.x, c.y - origin.y, c.z - origin.z);
        let along = v.dot(&d);
        let perp = Vec3::new(v.x - d.x * along, v.y - d.y * along, v.z - d.z * along);
        let dist_to_axis = perp.length();
        let apex_to_c = Vec3::new(c.x - apex.x, c.y - apex.y, c.z - apex.z);
        let mut interior_dir = Vec3::new(0.0, 0.0, 0.0);
        for n in face_normals {
            interior_dir.x -= n.x;
            interior_dir.y -= n.y;
            interior_dir.z -= n.z;
        }
        let side_dot = apex_to_c.dot(&interior_dir);
        if popkern_core::trace::trace_enabled("cap-centre") {
            eprintln!("[cap-centre]   dist_to_axis={dist_to_axis:.4e} side_dot={side_dot:.4e}");
        }
        if dist_to_axis < (r * 1e-6).max(1e-9) && side_dot > 0.0 {
            return Some(c);
        }
    }
    per_fillet_rolling_ball_centre(apex, fillet_cylinder)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Cube corner cap centre: apex (0,0,0), three faces with
    /// outward normals -x, -y, -z. Cap centre at distance r=0.1
    /// from each face (inside solid) is (0.1, 0.1, 0.1).
    #[test]
    fn cap_centre_cube_corner_uniform() {
        let apex = Point3::new(0.0, 0.0, 0.0);
        let normals = [
            Vec3::new(-1.0, 0.0, 0.0),
            Vec3::new(0.0, -1.0, 0.0),
            Vec3::new(0.0, 0.0, -1.0),
        ];
        let cap = compute_cap_centre(apex, &normals, 0.1).unwrap();
        assert!((cap.x - 0.1).abs() < 1e-12);
        assert!((cap.y - 0.1).abs() < 1e-12);
        assert!((cap.z - 0.1).abs() < 1e-12);
    }

    /// Cube apex translated off origin: cap centre still moves the
    /// cap radius along the inward direction (sum of -outward
    /// normals scaled).
    #[test]
    fn cap_centre_translated_apex() {
        let apex = Point3::new(2.0, -3.0, 5.0);
        let normals = [
            Vec3::new(-1.0, 0.0, 0.0),
            Vec3::new(0.0, -1.0, 0.0),
            Vec3::new(0.0, 0.0, -1.0),
        ];
        let cap = compute_cap_centre(apex, &normals, 0.1).unwrap();
        assert!((cap.x - 2.1).abs() < 1e-12);
        assert!((cap.y - (-2.9)).abs() < 1e-12);
        assert!((cap.z - 5.1).abs() < 1e-12);
    }

    /// Singular case: three colinear normals → no unique cap.
    #[test]
    fn cap_centre_singular_returns_none() {
        let apex = Point3::new(0.0, 0.0, 0.0);
        let normals = [
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(2.0, 0.0, 0.0),
            Vec3::new(-1.0, 0.0, 0.0),
        ];
        assert!(compute_cap_centre(apex, &normals, 0.1).is_none());
    }

    /// `fillet_nominal_radius` must accept Cylinder + Torus and
    /// reject other variants (sphere-fit only handles those two
    /// today).
    #[test]
    fn fillet_radius_extraction() {
        let cyl = SurfaceGeometry::Cylinder {
            origin: Point3::new(0.0, 0.0, 0.0),
            axis: Vec3::new(1.0, 0.0, 0.0),
            radius: 0.07,
        };
        assert_eq!(fillet_nominal_radius(&cyl), Some(0.07));

        let tor = SurfaceGeometry::Torus {
            center: Point3::new(0.0, 0.0, 0.0),
            axis: Vec3::new(0.0, 0.0, 1.0),
            major_radius: 1.0,
            minor_radius: 0.05,
        };
        assert_eq!(fillet_nominal_radius(&tor), Some(0.05));

        let sph = SurfaceGeometry::Sphere {
            center: Point3::new(0.0, 0.0, 0.0),
            radius: 0.1,
            axis: popkern_core::geometry::Vec3::Z,
        };
        assert_eq!(fillet_nominal_radius(&sph), None);
    }
}
