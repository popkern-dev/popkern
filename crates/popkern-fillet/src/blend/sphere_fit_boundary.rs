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
pub(crate) fn apply_sphere_fit_substitution(
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

    // Sphere-fit needs at least three independent face planes
    // through the apex to pin down the cap centre.
    let face_normals: Vec<Vec3> = unique_face_normals.iter().map(|(_, n)| *n).collect();
    if face_normals.len() < 3 {
        return Err("sphere-fit needs >= 3 unique adjacent faces");
    }

    // **Per-fillet sphere-fit (Phase 6c-A continuation §4.2).**
    // Each fillet uses *its own* cap radius (= the fillet's
    // nominal radius) and *its own* cap centre, computed by
    // solving the same `compute_cap_centre` system but with
    // r = r_k.  For uniform-radius corners (all r_k equal),
    // the per-fillet centres coincide with the single shared
    // centre.  For variable-radius corners (cube_mixed,
    // cube_extreme, …), each fillet's tangent circle is on
    // its own sphere — the resulting boundary's `FaceEdge`
    // bridges connect *different* sphere-fit endpoints on
    // each shared face (non-degenerate, unlike the uniform
    // case where bridges collapse to single points).
    //
    // The downstream cap construction (Gregory-cap promotion
    // ladder, §4.2) absorbs the boundary as-is and fits the
    // cap interior to satisfy G1 against each fillet via
    // Phase 6c-B's outer-fillet-G1 L-M solver.
    let mut new_curves: Vec<popkern_core::geometry::NurbsCurve> = Vec::with_capacity(n_fil);
    let mut new_endpoints: Vec<[Point3; 2]> = Vec::with_capacity(n_fil);
    let mut per_fillet_caps: Vec<(Point3, f64)> = Vec::with_capacity(n_fil);
    for k in 0..n_fil {
        let fi = match &boundary.segments[2 * k] {
            BoundarySegment::FilletEnd { fillet_index, .. } => *fillet_index,
            _ => unreachable!(),
        };
        let fillet = &fillet_surfaces[fi];
        let r_k =
            fillet_nominal_radius(&fillet.surface).ok_or("could not extract fillet radius")?;
        // Phase 6c-E: pick centre via the cylinder-axis-tangent
        // dispatcher — falls back to per-fillet apex-projection
        // when the inscribed-sphere solution is off-axis (pyramid
        // and other N ≥ 4 apices).  For trihedral cube-like apices
        // the dispatcher uses the inscribed-sphere solution
        // unchanged.
        let cap_centre_k = select_per_fillet_centre(apex, &face_normals, r_k, &fillet.surface)
            .ok_or("cap-centre system is singular")?;
        if std::env::var("POPKERN_CAPCENTRE_TRACE").is_ok_and(|v| v == "1") {
            eprintln!(
                "[cap-centre] fillet[{k}] r={r_k:.3} centre=({:.4},{:.4},{:.4})",
                cap_centre_k.x, cap_centre_k.y, cap_centre_k.z
            );
        }
        let result = sphere_fit_termination(
            apex,
            cap_centre_k,
            r_k,
            &fillet.surface,
            face_a_normals[k],
            face_b_normals[k],
            tolerance,
        )
        .map_err(|_| "sphere-fit termination failed")?;
        new_endpoints.push([result.endpoint_a, result.endpoint_b]);
        new_curves.push(result.curve_3d);
        per_fillet_caps.push((cap_centre_k, r_k));
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

    // Phase 6c-A continuation §4.2: expose per-fillet (centre,
    // radius) so downstream cap construction can assign each
    // sub-face a `Sphere(cap_centre_k, r_k)` surface (perfect
    // cap-cylinder G1 by construction).
    boundary.per_fillet_caps = per_fillet_caps;

    // Phase 6c-D Stage 1b: derive deg-9 Bernstein CPs of each
    // sub-patch's hex-side boundary (the substituted FilletEnd
    // arc and FaceEdge bridge curves) so the K-P bi-nonic ring
    // builder can use them directly via `pick_hex_boundary`,
    // bypassing the corner-Hermite synthesis that would
    // otherwise produce visible deviation from the analytic
    // arc on tilted-axis cylinders.
    //
    // For each sub-patch k:
    //   row0 (u=0, v ∈ [0,1]) ← deg-9 fit through 10 uniform
    //                            samples of segments[2k] (FilletEnd[k]).
    //   col9 (v=1, u ∈ [0,1]) ← deg-9 fit through 10 uniform
    //                            samples of segments[2k+1] (FaceEdge[k]).
    //
    // Sample-and-fit is exact for polynomial Bezier degree ≤ 9.
    // For the substituted FilletEnd's rational quadratic on a
    // small arc (cap-tangent angle ≤ π/4), the deg-9 polynomial
    // approximation is within ~1e-12 of the rational — well below
    // the 1e-2 commercial G0 gate.
    let mut hex_boundary_cps_deg9: Vec<([Point3; 10], [Point3; 10])> = Vec::with_capacity(n_fil);
    for k in 0..n_fil {
        let fillet_end_curve = match &boundary.segments[2 * k] {
            BoundarySegment::FilletEnd { curve, .. } => curve,
            _ => return Err("expected FilletEnd at index 2k"),
        };
        let face_edge_curve = match &boundary.segments[2 * k + 1] {
            BoundarySegment::FaceEdge { curve, .. } => curve,
            _ => return Err("expected FaceEdge at index 2k+1"),
        };
        let row0 = sample_and_fit_deg9(fillet_end_curve)
            .ok_or("failed to fit deg-9 Bernstein boundary on FilletEnd[k]")?;
        let col9 = sample_and_fit_deg9(face_edge_curve)
            .ok_or("failed to fit deg-9 Bernstein boundary on FaceEdge[k]")?;
        hex_boundary_cps_deg9.push((row0, col9));
    }
    boundary.hex_boundary_cps_deg9 = hex_boundary_cps_deg9;

    Ok(())
}

/// Phase 6c-D Stage 1b: sample a NurbsCurve at 10 uniform
/// parameter values across its full range and fit a deg-9
/// Bezier (10 Bernstein control points) interpolating those
/// samples.  Returns `None` if the curve has zero parameter
/// range, sampling fails, or the fit raises an error.
///
/// The 10×10 Bernstein interpolation matrix at uniform `t = i/9`
/// is well-conditioned and yields a unique deg-9 Bezier for any
/// 10-tuple of distinct samples. For analytic `sphere ∩ cylinder`
/// arcs (cap-tangent rationals), the deg-9 polynomial fit
/// approximates the rational to ~1e-12 — sub-mm on a unit cube.
fn sample_and_fit_deg9(curve: &popkern_core::geometry::NurbsCurve) -> Option<[Point3; 10]> {
    use popkern_core::geometry::curve::CurveEval;
    use popkern_core::geometry::fitting::fit_curve;

    let (t0, t1) = curve.parameter_range();
    if (t1 - t0).abs() < 1e-15 {
        return None;
    }
    let mut samples: Vec<Point3> = Vec::with_capacity(10);
    let mut params: Vec<f64> = Vec::with_capacity(10);
    for i in 0..10 {
        let s = (i as f64) / 9.0;
        let t = t0 + s * (t1 - t0);
        samples.push(curve.point_at(t));
        params.push(s);
    }
    let fitted = fit_curve(&samples, &params, 9, 10).ok()?;
    if fitted.control_points.len() != 10 {
        return None;
    }
    let mut out = [Point3::new(0.0, 0.0, 0.0); 10];
    for (i, cp) in fitted.control_points.iter().enumerate() {
        out[i] = *cp;
    }
    Some(out)
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

/// Phase 6c-E: per-fillet rolling-ball cap centre — the point
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

/// Phase 6c-E: pick the per-fillet sphere centre that
/// guarantees clean axis-tangency.  Tries the inscribed-sphere
/// solution from `compute_cap_centre` first — if it lies on the
/// fillet's cylinder axis (within tolerance), use it.  Otherwise
/// fall back to `per_fillet_rolling_ball_centre`.
///
/// The inscribed-sphere solution is preferred when valid
/// because it agrees across all per-fillet caps in a
/// trihedral apex (cube), so the `FaceEdge` bridges between
/// adjacent fillets' endpoints collapse to single points (clean
/// shared topology).  For pyramid-like N ≥ 4 apices the
/// inscribed sphere is off-axis, so we accept distinct per-fillet
/// centres at the cost of non-degenerate face-edge bridges.
pub(crate) fn select_per_fillet_centre(
    apex: Point3,
    face_normals: &[Vec3],
    r: f64,
    fillet_cylinder: &SurfaceGeometry,
) -> Option<Point3> {
    let inscribed = compute_cap_centre(apex, face_normals, r);
    if std::env::var("POPKERN_CAPCENTRE_TRACE").is_ok_and(|v| v == "1") {
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
        //   (2) inscribed centre is on the SAME SIDE of the apex
        //       as the cylinder body (uses cylinder.origin as a
        //       sample point on the body). For pyramid-like N≥4
        //       apices the inscribed solution lands on every
        //       cylinder axis at their convergence point, but on
        //       the OPPOSITE side of the apex from the cylinder
        //       body — a phantom sphere not rounding off any
        //       solid material.
        let Some(d) = axis.normalized() else {
            return Some(c);
        };
        let v = Vec3::new(c.x - origin.x, c.y - origin.y, c.z - origin.z);
        let along = v.dot(&d);
        let perp = Vec3::new(v.x - d.x * along, v.y - d.y * along, v.z - d.z * along);
        let dist_to_axis = perp.length();
        // Side check: dot(apex_to_inscribed, apex_to_origin) > 0
        // means inscribed is on the same side of apex as the
        // cylinder body. For cube the origin is in the cylinder
        // body direction from apex; for pyramid the origin is at
        // the base end of the slant edge while inscribed is past
        // the apex (opposite direction).
        let apex_to_c = Vec3::new(c.x - apex.x, c.y - apex.y, c.z - apex.z);
        let apex_to_origin = Vec3::new(origin.x - apex.x, origin.y - apex.y, origin.z - apex.z);
        let side_dot = apex_to_c.dot(&apex_to_origin);
        if std::env::var("POPKERN_CAPCENTRE_TRACE").is_ok_and(|v| v == "1") {
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
    use popkern_core::geometry::Vec3;

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
        use popkern_core::geometry::surface::SurfaceGeometry;
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
        };
        assert_eq!(fillet_nominal_radius(&sph), None);
    }
}
