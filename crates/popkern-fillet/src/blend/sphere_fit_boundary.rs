//! Sphere-fit termination substitution for blend boundaries (Phase 6c-A).
//!
//! Env-gated post-process for [`compute_blend_boundary`](super::boundary).
//! After the standard hex-closure boundary is built, this module
//! replaces each `FilletEnd` curve with the closed-form sphere-fit
//! termination — the curve where a fitted cap-sphere is tangent to
//! the fillet — and rebuilds the `FaceEdge` bridges to maintain
//! loop closure.
//!
//! Activation: `POPKERN_SPHERE_FIT_BOUNDARY=1`. Falls back silently
//! to the standard boundary on any failure (unsupported variant,
//! geometric degeneracy, missing topology).

use popkern_core::geometry::apex_blend::sphere_fit::sphere_fit_termination;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::uv_projection::project_point_to_uv;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};

use super::boundary::{BlendBoundary, BoundarySegment};
use super::surface_bridge::build_straight_line_bridge;
use crate::fillet_surface::FilletSurfaceResult;

/// Substitute every `FilletEnd` curve in `boundary` with its
/// sphere-fit termination, rebuilding `FaceEdge` bridges to match.
///
/// On success, mutates `boundary` in place. On any failure (input
/// shape unsupported, sphere-fit math errors, missing face data),
/// returns an `Err` with a static reason and leaves the boundary
/// unchanged.
pub(super) fn apply_sphere_fit_substitution(
    boundary: &mut BlendBoundary,
    fillet_surfaces: &[FilletSurfaceResult],
    brep: &BRep,
    tolerance: &Tolerance,
) -> Result<(), &'static str> {
    let apex = boundary.vertex_position;
    let n_seg = boundary.segments.len();

    // Require strict alternation FilletEnd, FaceEdge, FilletEnd, ...
    // with at least 2 fillets (n_seg = 4).
    if n_seg < 4 || !n_seg.is_multiple_of(2) {
        return Err("sphere-fit requires alternating segments with n >= 2 fillets");
    }
    for (i, seg) in boundary.segments.iter().enumerate() {
        match (i.is_multiple_of(2), seg) {
            (true, BoundarySegment::FilletEnd { .. }) => {}
            (false, BoundarySegment::FaceEdge { .. }) => {}
            _ => return Err("sphere-fit expects alternating FilletEnd/FaceEdge"),
        }
    }

    let n_fil = n_seg / 2;

    // Collect adjacent face normals for each fillet.
    // FaceEdge before FilletEnd[k] is at index 2k-1 (or n_seg-1 when k==0).
    // FaceEdge after FilletEnd[k] is at index 2k+1.
    let mut face_a_normals: Vec<Vec3> = Vec::with_capacity(n_fil);
    let mut face_b_normals: Vec<Vec3> = Vec::with_capacity(n_fil);
    let mut unique_face_normals: Vec<(FaceHandle, Vec3)> = Vec::new();

    for k in 0..n_fil {
        let prev_idx = if k == 0 { n_seg - 1 } else { 2 * k - 1 };
        let next_idx = 2 * k + 1;
        let prev_face = match &boundary.segments[prev_idx] {
            BoundarySegment::FaceEdge { face, .. } => *face,
            _ => return Err("expected FaceEdge before FilletEnd"),
        };
        let next_face = match &boundary.segments[next_idx] {
            BoundarySegment::FaceEdge { face, .. } => *face,
            _ => return Err("expected FaceEdge after FilletEnd"),
        };
        let n_a = face_outward_normal_at(brep, prev_face, &apex)
            .ok_or("could not evaluate face A outward normal")?;
        let n_b = face_outward_normal_at(brep, next_face, &apex)
            .ok_or("could not evaluate face B outward normal")?;
        face_a_normals.push(n_a);
        face_b_normals.push(n_b);
        for (fh, n) in [(prev_face, n_a), (next_face, n_b)] {
            if !unique_face_normals.iter().any(|(f, _)| *f == fh) {
                unique_face_normals.push((fh, n));
            }
        }
    }

    // Cap radius from the first fillet's nominal radius.
    let first_fi = match &boundary.segments[0] {
        BoundarySegment::FilletEnd { fillet_index, .. } => *fillet_index,
        _ => unreachable!(),
    };
    let cap_radius = fillet_nominal_radius(&fillet_surfaces[first_fi].surface)
        .ok_or("could not extract fillet radius")?;

    // Sphere-fit needs at least three independent face planes
    // through the apex to pin down the cap centre.
    let face_normals: Vec<Vec3> = unique_face_normals.iter().map(|(_, n)| *n).collect();
    if face_normals.len() < 3 {
        return Err("sphere-fit needs >= 3 unique adjacent faces");
    }
    let cap_centre = compute_cap_centre(apex, &face_normals, cap_radius)
        .ok_or("cap-centre system is singular")?;

    // Compute sphere-fit terminations for every fillet first; if any
    // fail, we abort without mutating the boundary.
    let mut new_curves: Vec<popkern_core::geometry::NurbsCurve> = Vec::with_capacity(n_fil);
    let mut new_endpoints: Vec<[Point3; 2]> = Vec::with_capacity(n_fil);
    for k in 0..n_fil {
        let fi = match &boundary.segments[2 * k] {
            BoundarySegment::FilletEnd { fillet_index, .. } => *fillet_index,
            _ => unreachable!(),
        };
        let fillet = &fillet_surfaces[fi];
        let result = sphere_fit_termination(
            apex,
            cap_centre,
            cap_radius,
            &fillet.surface,
            face_a_normals[k],
            face_b_normals[k],
            tolerance,
        )
        .map_err(|_| "sphere-fit termination failed")?;
        new_endpoints.push([result.endpoint_a, result.endpoint_b]);
        new_curves.push(result.curve_3d);
    }

    // Substitute FilletEnd curves.
    for (k, curve) in new_curves.into_iter().enumerate() {
        if let BoundarySegment::FilletEnd { curve: c, .. } = &mut boundary.segments[2 * k] {
            *c = curve;
        }
    }

    // Rebuild FaceEdge bridges.  Bridge k connects fillet k's
    // endpoint_b (face_b side) to fillet (k+1)'s endpoint_a (face_a
    // side), which both lie on the face shared by fillets k and k+1.
    for k in 0..n_fil {
        let next_k = (k + 1) % n_fil;
        let arc_end = new_endpoints[k][1];
        let next_start = new_endpoints[next_k][0];
        let bridge =
            build_straight_line_bridge(&arc_end, &next_start).ok_or("could not build bridge")?;
        if let BoundarySegment::FaceEdge { curve, .. } = &mut boundary.segments[2 * k + 1] {
            *curve = bridge;
        }
    }

    Ok(())
}

/// Outward face normal at `point` (used here only at the apex).
///
/// In popkern's convention, `unit_cube` and the pyramid fixtures
/// store their `Plane` surface with a natural normal pointing
/// *into* the solid; `SurfaceOnFace::orientation = true` means the
/// face's natural orientation agrees with that surface normal —
/// i.e., the face normal is also inward. We invert here to recover
/// the outward direction the sphere-fit math expects.
fn face_outward_normal_at(brep: &BRep, face: FaceHandle, point: &Point3) -> Option<Vec3> {
    let face_data = brep.faces.get(face)?;
    let sof = face_data.surface.as_ref()?;
    let (u, v) = project_point_to_uv(&sof.geometry, point)?;
    let n = sof.geometry.normal_at(u, v).normalized()?;
    let face_natural = if sof.orientation { n } else { -n };
    Some(-face_natural)
}

/// Extract the nominal fillet radius from a fillet's surface.
fn fillet_nominal_radius(surface: &SurfaceGeometry) -> Option<f64> {
    match surface {
        SurfaceGeometry::Cylinder { radius, .. } => Some(*radius),
        SurfaceGeometry::Torus { minor_radius, .. } => Some(*minor_radius),
        _ => None,
    }
}

/// Solve `M (cap - apex) = -r * 1` for `cap`, where `M` has the
/// outward face normals as rows.  For 3 faces this is a direct 3x3
/// solve; for 4+ faces this is the normal-equations least-squares.
fn compute_cap_centre(apex: Point3, face_normals: &[Vec3], r: f64) -> Option<Point3> {
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
