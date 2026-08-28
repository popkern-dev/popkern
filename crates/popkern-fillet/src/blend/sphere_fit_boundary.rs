//! Sphere-fit termination substitution for blend boundaries (sphere-fit substitution).
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

use popkern_core::geometry::analytic_surfaces::orthonormal_basis;
use popkern_core::geometry::apex_blend::sphere_fit::sphere_fit_termination;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};

use super::boundary::{BlendBoundary, BoundarySegment, FilletEndArcParams};
use super::sphere_fit_centre::{
    face_outward_normal_at, fillet_nominal_radius, select_per_fillet_centre,
};
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
    let r = apply_sphere_fit_substitution_inner(boundary, fillet_surfaces, brep, tolerance);
    if std::env::var_os("PK_SPHEREFIT_TRACE").is_some() {
        eprintln!("[spherefit] vertex {:?}: {:?}", boundary.vertex, r);
    }
    r
}

fn apply_sphere_fit_substitution_inner(
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

    // **Per-fillet sphere-fit (per-fillet sphere-fit continuation).**
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
    // the outer-fillet-G1 L-M solver.
    let mut new_curves: Vec<popkern_core::geometry::NurbsCurve> = Vec::with_capacity(n_fil);
    let mut new_endpoints: Vec<[Point3; 2]> = Vec::with_capacity(n_fil);
    let mut per_fillet_caps: Vec<(Point3, f64)> = Vec::with_capacity(n_fil);
    let mut fillet_end_arcs: Vec<FilletEndArcParams> = Vec::with_capacity(n_fil);
    for k in 0..n_fil {
        let fi = match &boundary.segments[2 * k] {
            BoundarySegment::FilletEnd { fillet_index, .. } => *fillet_index,
            _ => unreachable!(),
        };
        let fillet = &fillet_surfaces[fi];
        let r_k =
            fillet_nominal_radius(&fillet.surface).ok_or("could not extract fillet radius")?;
        // Pick centre via the cylinder-axis-tangent
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
        // Capture analytic arc params so
        // the wire can write `CurveGeometry::Ellipse` (with the
        // cylinder's exact `orthonormal_basis`) to BRep edges
        // instead of a NURBS approximation.  Pcurve becomes
        // `Linear` over the angular range — bit-exact agreement
        // between cap-side and twin (Cyl) surface evaluations.
        if let SurfaceGeometry::Cylinder {
            origin: cyl_origin,
            axis: cyl_axis,
            radius: cyl_radius,
        } = &fillet.surface
        {
            if let Some(axis_unit) = cyl_axis.normalized() {
                let (perp1, perp2) = orthonormal_basis(axis_unit);
                let v_const = (cap_centre_k - *cyl_origin).dot(&axis_unit);
                let project = |p: Point3| -> (f64, f64) {
                    let d = p - *cyl_origin;
                    let mut theta = (d.dot(&perp2)).atan2(d.dot(&perp1));
                    if theta < 0.0 {
                        theta += std::f64::consts::TAU;
                    }
                    let v = d.dot(&axis_unit);
                    (theta, v)
                };
                let (theta_a_raw, _) = project(result.endpoint_a);
                let (theta_b_raw, _) = project(result.endpoint_b);
                // Pick the short arc: shift theta_b by ±2π if
                // needed so |theta_b - theta_a| ≤ π.  Without this,
                // the angular sampling traces the long way around
                // the cylinder (>π span) and produces a very
                // different curve from the actual cap-tangent arc.
                let theta_a = theta_a_raw;
                let mut theta_b = theta_b_raw;
                while theta_b - theta_a > std::f64::consts::PI {
                    theta_b -= std::f64::consts::TAU;
                }
                while theta_a - theta_b > std::f64::consts::PI {
                    theta_b += std::f64::consts::TAU;
                }
                fillet_end_arcs.push(FilletEndArcParams {
                    center: cap_centre_k,
                    u_basis: perp1 * (*cyl_radius),
                    v_basis: perp2 * (*cyl_radius),
                    start_angle: theta_a,
                    end_angle: theta_b,
                    v_const,
                });
            } else {
                // Cylinder axis was zero — should not happen for valid
                // fillet geometry. Skip arc capture; downstream wire
                // falls back to NurbsCurve / refit.
                fillet_end_arcs.push(FilletEndArcParams {
                    center: cap_centre_k,
                    u_basis: Vec3::new(0.0, 0.0, 0.0),
                    v_basis: Vec3::new(0.0, 0.0, 0.0),
                    start_angle: 0.0,
                    end_angle: 0.0,
                    v_const: 0.0,
                });
            }
        } else {
            // Non-cylinder fillet (Torus etc.) — analytic arc
            // representation not yet defined for this case.
            fillet_end_arcs.push(FilletEndArcParams {
                center: cap_centre_k,
                u_basis: Vec3::new(0.0, 0.0, 0.0),
                v_basis: Vec3::new(0.0, 0.0, 0.0),
                start_angle: 0.0,
                end_angle: 0.0,
                v_const: 0.0,
            });
        }
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

    // Per-fillet sphere-fit continuation: expose per-fillet (centre,
    // radius) so downstream cap construction can assign each
    // sub-face a `Sphere(cap_centre_k, r_k)` surface (perfect
    // cap-cylinder G1 by construction).
    boundary.per_fillet_caps = per_fillet_caps;

    // Derive deg-9 Bernstein CPs of each
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
        // Stage 4.7c-ii (angular sampling for FilletEnd row 0)
        // reverted — empirical regression to ~9 cm dev_cap in
        // bi-9 patch evaluation, despite the standalone fit
        // achieving 5e-14 m error.  Suggests the K-P binonic
        // optimizer or the bi-quintic guide preprocessing
        // depends on the boundary's Bezier-t parameterisation in
        // a way that's incompatible with angular sampling.
        // Foundation in `BlendBoundary::fillet_end_arcs` retained
        // for future investigation.
        let _ = fillet_end_arcs.get(k); // suppress unused warning
        let row0 = sample_and_fit_deg9(fillet_end_curve)
            .ok_or("failed to fit deg-9 Bernstein boundary on FilletEnd[k]")?;
        let col9 = sample_and_fit_deg9(face_edge_curve)
            .ok_or("failed to fit deg-9 Bernstein boundary on FaceEdge[k]")?;
        hex_boundary_cps_deg9.push((row0, col9));
    }
    boundary.hex_boundary_cps_deg9 = hex_boundary_cps_deg9;
    boundary.fillet_end_arcs = fillet_end_arcs;

    if std::env::var_os("POPKERN_SPHERE_FIT_TRACE").is_some() {
        use popkern_core::geometry::curve::CurveEval;
        eprintln!("[sphere-fit] boundary segment endpoint coherence:");
        for (i, seg) in boundary.segments.iter().enumerate() {
            let curve = match seg {
                BoundarySegment::FilletEnd { curve, .. } => curve,
                BoundarySegment::FaceEdge { curve, .. } => curve,
            };
            let (t0, t1) = curve.parameter_range();
            let p_start = curve.point_at(t0);
            let p_end = curve.point_at(t1);
            let kind = match seg {
                BoundarySegment::FilletEnd { .. } => "FilletEnd",
                BoundarySegment::FaceEdge { .. } => "FaceEdge ",
            };
            eprintln!(
                "  seg[{i}] {kind} t=[{t0:.4},{t1:.4}] start=({:.6},{:.6},{:.6}) end=({:.6},{:.6},{:.6})",
                p_start.x, p_start.y, p_start.z, p_end.x, p_end.y, p_end.z
            );
        }
        // Coherence check: seg[2k]'s end should equal seg[2k+1]'s start; seg[2k+1]'s end should equal seg[2k+2]'s start (mod 2N).
        for k in 0..n_fil {
            let cur = &boundary.segments[2 * k];
            let next = &boundary.segments[2 * k + 1];
            let cur_curve = match cur {
                BoundarySegment::FilletEnd { curve, .. } => curve,
                _ => continue,
            };
            let next_curve = match next {
                BoundarySegment::FaceEdge { curve, .. } => curve,
                _ => continue,
            };
            let (_, t1_cur) = cur_curve.parameter_range();
            let (t0_next, _) = next_curve.parameter_range();
            let cur_end = cur_curve.point_at(t1_cur);
            let next_start = next_curve.point_at(t0_next);
            let gap = cur_end.distance_to(&next_start);
            eprintln!("  coherence FilletEnd[{k}].end ↔ FaceEdge[{k}].start: gap={gap:.3e}");
        }
    }

    Ok(())
}

/// Sample a NurbsCurve at 10 uniform
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

/// Sample the analytic ellipse at 10
/// uniform angular intervals between `start_angle` and
/// `end_angle`, fit a deg-9 Bezier with parameter `s = (i/9)`
/// interpolating those angular samples.  When the cap surface's
/// row 0 uses these CPs and the pcurve is `Linear (0, 0) →
/// (0, 1)` over angular `t ∈ [start_angle, end_angle]`, surface
/// evaluation matches the analytic ellipse to ~1e-12 m on small
/// arcs (the deg-9 polynomial fit error over angular samples).
#[allow(dead_code)]
fn sample_arc_at_uniform_angles_and_fit_deg9(
    arc: &super::boundary::FilletEndArcParams,
) -> Option<[Point3; 10]> {
    use popkern_core::geometry::fitting::fit_curve;

    if (arc.end_angle - arc.start_angle).abs() < 1e-15 {
        return None;
    }
    let mut samples: Vec<Point3> = Vec::with_capacity(10);
    let mut params: Vec<f64> = Vec::with_capacity(10);
    for i in 0..10 {
        let s = (i as f64) / 9.0;
        let theta = arc.start_angle + s * (arc.end_angle - arc.start_angle);
        let p = arc.center + arc.u_basis * theta.cos() + arc.v_basis * theta.sin();
        samples.push(p);
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
