//! Chamfer curve and surface computation.
//!
//! For each edge to be chamfered, computes offset curves on each adjacent face
//! and builds a ruled surface (flat plane) between them.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{NurbsCurve, NurbsSurface, Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use crate::analysis::EdgeAnalysis;
use crate::center_curve::fit_nurbs_through_points;
use crate::types::ChamferKind;

#[path = "chamfer_offset_rows.rs"]
mod chamfer_offset_rows;
use chamfer_offset_rows::{
    compute_offset_distances, compute_offset_points, face_interior_hint, find_uv_for_point,
};

/// Result of chamfer geometry computation for one edge.
#[derive(Debug)]
pub struct ChamferResult {
    /// The ruled surface connecting the two trim curves.
    pub surface: SurfaceGeometry,
    /// 3D trim curve on the first face.
    pub trim_curve_1: NurbsCurve,
    /// 3D trim curve on the second face.
    pub trim_curve_2: NurbsCurve,
    /// Exact analytic circle behind `trim_curve_1` when the bevel
    /// was recognized as a cone frustum (closed rims) — the
    /// closed-rim surgery builds exact Arc edges from it.
    pub trim_circle_1: Option<(Point3, Vec3, f64)>,
    /// Exact analytic circle behind `trim_curve_2`.
    pub trim_circle_2: Option<(Point3, Vec3, f64)>,
    /// Parameter-space trim curve on the first face (Point3(u, v, 0) convention).
    pub pcurve_1: NurbsCurve,
    /// Parameter-space trim curve on the second face (Point3(u, v, 0) convention).
    pub pcurve_2: NurbsCurve,
}

/// Compute chamfer geometry for an edge.
///
/// Based on the `ChamferKind`, computes the offset distance on each face,
/// finds trim curves at the offset distance, builds a ruled surface between
/// them, and computes parameter-space curves on each face.
pub fn compute_chamfer(
    analysis: &EdgeAnalysis,
    brep: &BRep,
    kind: &ChamferKind,
    tolerance: &Tolerance,
) -> Result<ChamferResult, String> {
    // 1. Determine offset distances on each face.
    let (dist_1, dist_2) = compute_offset_distances(analysis, kind)?;

    // 2. Get the edge curve.
    let edge = brep
        .edges
        .get(analysis.edge)
        .ok_or_else(|| format!("edge {:?} not found", analysis.edge))?;
    let curve = edge
        .curve
        .as_ref()
        .ok_or_else(|| format!("edge {:?} has no curve geometry", analysis.edge))?;

    // 3. Sample points along the edge.
    let n_samples = 20usize;
    let (t_start, t_end) = (curve.t_start, curve.t_end);

    let mut edge_points = Vec::with_capacity(n_samples);
    let mut edge_tangents = Vec::with_capacity(n_samples);
    for i in 0..n_samples {
        let t = t_start + (t_end - t_start) * (i as f64) / ((n_samples - 1).max(1) as f64);
        edge_points.push(curve.geometry.point_at(t));
        // Exact tangents: chord tangents carry ~1e-2 rad noise at 20
        // samples; where back-projection can't absorb it (a plane)
        // the row spread (1.45e-2·distance) crosses the 1e-3·r cone
        // recognition gate past ~2 mm legs → cracked NURBS band.
        edge_tangents.push(curve.geometry.derivative_at(t, 1));
    }

    // 4. Compute offset points on each face. Each face's offset must
    // walk INTO that face's material — `normal × tangent` alone is
    // sign-blind to which side the face lies (measured: the face-2
    // offset row landed OUTSIDE the solid, so every chamfer bevel
    // was authored perpendicular to the true one). The interior
    // side is derived from the winding invariant: interior =
    // outer-loop Newell normal × the face's own half-edge direction
    // along this edge.
    let interior_1 = face_interior_hint(brep, analysis.edge, analysis.face_1);
    let interior_2 = face_interior_hint(brep, analysis.edge, analysis.face_2);
    let offset_points_1 = compute_offset_points(
        &edge_points,
        &edge_tangents,
        &analysis.surface_1,
        dist_1,
        interior_1.as_ref(),
        tolerance,
    )?;
    let offset_points_2 = compute_offset_points(
        &edge_points,
        &edge_tangents,
        &analysis.surface_2,
        dist_2,
        interior_2.as_ref(),
        tolerance,
    )?;

    if std::env::var_os("PK_BEVEL_TRACE").is_some() {
        let a = offset_points_1.first();
        let b = offset_points_2.first();
        eprintln!("  [chamfer-offsets] p1[0]={a:?} p2[0]={b:?}");
    }
    // 4.5. Conical bevel recognition: a chamfer around a circular
    // rim is an exact CONE frustum bounded by two circles (cap∩cone
    // and wall∩cone). Snap the sampled rows onto the exact circles
    // BEFORE fitting trims — the chord-tangent offsets carry ~1e-3·r
    // noise that would otherwise live on in every rail, and the
    // generic closed NURBS ruled ring meshed with cracked seams
    // (measured: ~1800 open mesh edges on a plain rim chamfer).
    let mut offset_points_1 = offset_points_1;
    let mut offset_points_2 = offset_points_2;
    let cone = recognize_conical_bevel(&offset_points_1, &offset_points_2);
    if let Some((_, circle_1, circle_2)) = &cone {
        snap_row_to_circle(&mut offset_points_1, circle_1);
        snap_row_to_circle(&mut offset_points_2, circle_2);
    }

    // 5. Fit trim curves through the offset points. Recognized rim
    // circles become EXACT rational-arc NURBS anchored at the row
    // start (the seam station): the interpolated 20-point fit left
    // ~1e-3·r ripple on every rail and the wall/cap meshed against
    // it with thousands of cracked slivers.
    let exact_rim = |circle: &BevelCircle, anchor: &Point3| -> NurbsCurve {
        let (c, axis, r) = circle;
        // Angles measured in the Arc EVALUATOR's frame (circle_basis)
        // so that point_at(a0) IS the anchor — the surface frame
        // (orthonormal_basis) can sit a phase apart and the authored
        // start then names a different point on the circle.
        let (p1, p2) = popkern_core::geometry::analytic_curves::circle_basis(*axis);
        let d = *anchor - *c;
        let a0 = d.dot(&p2).atan2(d.dot(&p1));
        if std::env::var_os("PK_RIM_TRACE").is_some() {
            use popkern_core::geometry::curve::CurveEval as _;
            let probe = popkern_core::geometry::curve::CurveGeometry::Arc {
                center: *c,
                normal: *axis,
                radius: *r,
                start_angle: a0,
                end_angle: a0 + std::f64::consts::TAU,
            };
            let p = probe.point_at(a0);
            eprintln!(
                "  [rim-anchor] exact_rim anchor=({:.3},{:.3},{:.3}) a0={a0:.4} residual={:.2e}",
                anchor.x,
                anchor.y,
                anchor.z,
                p.distance_to(anchor)
            );
        }
        popkern_core::geometry::curve::CurveGeometry::Arc {
            center: *c,
            normal: *axis,
            radius: *r,
            start_angle: a0,
            end_angle: a0 + std::f64::consts::TAU,
        }
        .to_nurbs_in_range((a0, a0 + std::f64::consts::TAU))
    };
    let (trim_curve_1, trim_curve_2) = match &cone {
        Some((_, circle_1, circle_2)) => (
            exact_rim(circle_1, &offset_points_1[0]),
            exact_rim(circle_2, &offset_points_2[0]),
        ),
        None => (
            fit_nurbs_through_points(&offset_points_1)?,
            fit_nurbs_through_points(&offset_points_2)?,
        ),
    };

    // 6. Compute pcurves by projecting onto face UV spaces.
    let pcurve_1 = compute_pcurve(&offset_points_1, &analysis.surface_1)?;
    let pcurve_2 = compute_pcurve(&offset_points_2, &analysis.surface_2)?;

    // 7. Build the surface: the recognized cone, or the ruled
    // strip between the two trim curves.
    let (trim_circle_1, trim_circle_2) = match &cone {
        Some((_, c1, c2)) => (Some(*c1), Some(*c2)),
        None => (None, None),
    };
    let surface = match cone {
        Some((cone_surface, _, _)) => {
            if std::env::var_os("PK_BEVEL_TRACE").is_some() {
                eprintln!("  [bevel] conical recognition: {cone_surface:?}");
            }
            cone_surface
        }
        None => build_ruled_surface(&offset_points_1, &offset_points_2)?,
    };

    Ok(ChamferResult {
        surface,
        trim_curve_1,
        trim_curve_2,
        trim_circle_1,
        trim_circle_2,
        pcurve_1,
        pcurve_2,
    })
}

/// Compute a parameter-space (pcurve) by projecting 3D points into surface UV.
fn compute_pcurve(points: &[Point3], surface: &SurfaceGeometry) -> Result<NurbsCurve, String> {
    let mut uv_points = Vec::with_capacity(points.len());
    for pt in points {
        let (u, v) = find_uv_for_point(surface, pt);
        uv_points.push(Point3::new(u, v, 0.0));
    }
    fit_nurbs_through_points(&uv_points)
}

/// Build a ruled (degree-1 in v) NURBS surface between two sets of points.
///
/// The surface is degree-1 in both u (along edge) and v (across chamfer).
/// Row 0 = points on face 1, row 1 = points on face 2.
fn build_ruled_surface(
    points_1: &[Point3],
    points_2: &[Point3],
) -> Result<SurfaceGeometry, String> {
    let n = points_1.len();
    if n < 2 || points_2.len() != n {
        return Err("ruled surface requires matching point arrays with >= 2 points".into());
    }

    // Control points: row-major, [row0 = points_1, row1 = points_2].
    // U direction = along the edge (n points), V direction = across chamfer (2 points).
    let mut control_points = Vec::with_capacity(2 * n);
    control_points.extend_from_slice(points_1);
    control_points.extend_from_slice(points_2);

    let degree_u = 1.min(n - 1);
    let degree_v = 1;

    let knots_u = build_clamped_uniform_knots(n, degree_u);
    let knots_v = vec![0.0, 0.0, 1.0, 1.0]; // degree-1 with 2 control points

    let weights = vec![1.0; 2 * n];

    // Flat bevel recognition: a chamfer across a straight edge
    // between planar supports is itself a PLANE. Emitting the
    // analytic chart (like straight-edge fillets emit Cylinder)
    // keeps the blend on the kernel's analytic fast paths — pcurve
    // attachment, tessellation, STEP export. Verified point-wise:
    // every control point must lie on the candidate plane to 1e-9
    // of the bevel scale.
    if let Some(plane) = recognize_planar_bevel(&control_points) {
        return Ok(plane);
    }

    let surface = NurbsSurface {
        degree_u,
        degree_v,
        knots_u,
        knots_v,
        control_points,
        num_u: n,
        num_v: 2,
        weights,
    };

    Ok(SurfaceGeometry::Nurbs(surface))
}

/// Recognize a flat ruled bevel as an analytic plane. The chart runs
/// u along the edge and v across the bevel; `None` when the points
/// are not coplanar (curved-edge chamfers stay NURBS).
fn recognize_planar_bevel(control_points: &[Point3]) -> Option<SurfaceGeometry> {
    let n = control_points.len() / 2;
    if n < 2 {
        return None;
    }
    let origin = control_points[0];
    let u_dir = (control_points[n - 1] - origin).normalized()?;
    let across = control_points[n] - origin;
    let v_dir = (across - u_dir * across.dot(&u_dir)).normalized()?;
    let normal = u_dir.cross(&v_dir).normalized()?;
    let scale = control_points
        .iter()
        .map(|p| (*p - origin).length())
        .fold(0.0_f64, f64::max)
        .max(1e-9);
    for p in control_points {
        let dev = (*p - origin).dot(&normal).abs();
        if dev > 1e-9 * scale {
            if std::env::var_os("PK_BEVEL_TRACE").is_some() {
                eprintln!("  [bevel] planar recognition rejected: dev={dev:.3e} scale={scale:.1}");
            }
            return None;
        }
    }
    Some(SurfaceGeometry::Plane {
        origin,
        normal,
        u_axis: u_dir,
        v_axis: v_dir,
    })
}

/// Recognize a conical (or cylindrical) bevel from its two offset
/// rows: both rows must be concentric circles about a common axis.
/// Returns `None` when either row is not circular to 1e-6 of its
/// radius, the axes disagree, or the frustum degenerates.
type BevelCircle = (Point3, Vec3, f64);

fn recognize_conical_bevel(
    points_1: &[Point3],
    points_2: &[Point3],
) -> Option<(SurfaceGeometry, BevelCircle, BevelCircle)> {
    let fit_circle = |pts: &[Point3]| -> Option<(Point3, Vec3, f64, f64)> {
        // Closed rows duplicate the closure sample — drop it or the
        // centroid biases toward it by r/n and every radius check
        // fails.
        let pts = match (pts.first(), pts.last()) {
            (Some(a), Some(b)) if pts.len() > 4 && a.distance_to(b) < 1e-9 => &pts[..pts.len() - 1],
            _ => pts,
        };
        if pts.len() < 4 {
            return None;
        }
        let n = pts.len() as f64;
        let c = pts.iter().fold(Vec3::new(0.0, 0.0, 0.0), |a, p| {
            Vec3::new(a.x + p.x / n, a.y + p.y / n, a.z + p.z / n)
        });
        let mean = Point3::new(c.x, c.y, c.z);
        // Newell normal over the sampled ring.
        let mut nw = Vec3::new(0.0, 0.0, 0.0);
        for i in 0..pts.len() {
            let p = pts[i];
            let q = pts[(i + 1) % pts.len()];
            nw.x += (p.y - q.y) * (p.z + q.z);
            nw.y += (p.z - q.z) * (p.x + q.x);
            nw.z += (p.x - q.x) * (p.y + q.y);
        }
        let axis = nw.normalized()?;
        // Least-squares (Kåsa) circle centre in the ring's plane —
        // the raw point mean biases toward any angular gap in the
        // sampling (measured: a rim row starting 0.5 short of
        // closure shifted the mean by r/n and every radius check
        // failed).
        let (bu, bv) = {
            let any = if axis.x.abs() < 0.9 {
                Vec3::new(1.0, 0.0, 0.0)
            } else {
                Vec3::new(0.0, 1.0, 0.0)
            };
            let bu = axis.cross(&any).normalized()?;
            (bu, axis.cross(&bu))
        };
        let (mut suu, mut suv, mut svv, mut su3, mut sv3) = (0.0, 0.0, 0.0, 0.0, 0.0);
        let (mut suvv, mut svuu) = (0.0, 0.0);
        for p in pts {
            let d = *p - mean;
            let u = d.dot(&bu);
            let v = d.dot(&bv);
            suu += u * u;
            suv += u * v;
            svv += v * v;
            su3 += u * u * u;
            sv3 += v * v * v;
            suvv += u * v * v;
            svuu += v * u * u;
        }
        let det = suu * svv - suv * suv;
        if det.abs() < 1e-18 {
            return None;
        }
        let rhs_u = 0.5 * (su3 + suvv);
        let rhs_v = 0.5 * (sv3 + svuu);
        let cu = (rhs_u * svv - rhs_v * suv) / det;
        let cv = (rhs_v * suu - rhs_u * suv) / det;
        let centroid = mean + bu * cu + bv * cv;
        if std::env::var_os("PK_BEVEL_TRACE").is_some() {
            for (i, p) in pts.iter().enumerate() {
                let d = *p - centroid;
                let radial = (d - axis * d.dot(&axis)).length();
                eprintln!(
                    "  [bevel-sample] {i}: r={radial:.4} p=({:.2},{:.2},{:.2})",
                    p.x, p.y, p.z
                );
            }
        }
        let mut r_sum = 0.0;
        let mut r_min = f64::MAX;
        let mut r_max: f64 = 0.0;
        for p in pts {
            let d = *p - centroid;
            let radial = d - axis * d.dot(&axis);
            let axial = d.dot(&axis).abs();
            // Chord-tangent offsets carry ~1e-3-scale sampling
            // noise both radially and axially; rows are snapped
            // onto the fitted circle afterwards, so the gates only
            // exclude genuinely non-circular rows.
            if axial > 1e-3 * (1.0 + radial.length()) {
                if std::env::var_os("PK_BEVEL_TRACE").is_some() {
                    eprintln!("  [bevel-circle] non-coplanar axial={axial:.3e}");
                }
                return None; // not coplanar
            }
            let r = radial.length();
            r_sum += r;
            r_min = r_min.min(r);
            r_max = r_max.max(r);
        }
        let r = r_sum / n;
        if std::env::var_os("PK_BEVEL_TRACE").is_some() {
            eprintln!(
                "  [bevel-circle] n={} r={r:.4} spread={:.3e} c=({:.2},{:.2},{:.2})",
                pts.len(),
                r_max - r_min,
                centroid.x,
                centroid.y,
                centroid.z
            );
            for (i, pt) in pts.iter().enumerate().take(6) {
                eprintln!("    [{i}] ({:.3},{:.3},{:.3})", pt.x, pt.y, pt.z);
            }
        }
        // Chord-based offset tangents carry ~1e-3·r sampling noise
        // (one-sided secants at the row ends) — the fit centre and
        // radius are least-squares exact; rows are SNAPPED onto the
        // fitted circle by the caller, so the gate only needs to
        // exclude genuinely non-circular rows.
        if r < 1e-9 || (r_max - r_min) > 1e-3 * r {
            return None; // not circular
        }
        Some((centroid, axis, r, r_max - r_min))
    };
    let trace = std::env::var_os("PK_BEVEL_TRACE").is_some();
    let (c1, a1, r1, spread1) = fit_circle(points_1)?;
    let (c2, a2, r2, spread2) = fit_circle(points_2)?;
    let scale = r1.max(r2);
    // Common axis: parallel normals, centre offset along the axis.
    if a1.cross(&a2).length() > 1e-6 {
        if trace {
            eprintln!("  [bevel-cone] axes disagree: {a1:?} vs {a2:?}");
        }
        return None;
    }
    // RECONCILE the two fits onto one axis: rows carry unequal
    // sampling noise (an analytic-support row fits to 1e-14 while a
    // chord-tangent row wobbles ~1e-2), so the CLEANER circle's axis
    // line is authoritative — re-centre the noisier circle onto it
    // and refit its radius about the corrected centre. Without this
    // the two centres disagree perpendicular to the axis by the
    // noise (measured 3.9e-3) and recognition declines.
    let refit = |c_fixed: Point3, axis: &Vec3, pts: &[Point3]| -> f64 {
        let mut sum = 0.0;
        for p in pts {
            let d = *p - c_fixed;
            sum += (d - *axis * d.dot(axis)).length();
        }
        sum / pts.len() as f64
    };
    let (c1, r1, c2, r2) = if spread1 <= spread2 {
        let c2f = c1 + a1 * (c2 - c1).dot(&a1);
        (c1, r1, c2f, refit(c2f, &a1, points_2))
    } else {
        let c1f = c2 + a2 * (c1 - c2).dot(&a2);
        (c1f, refit(c1f, &a2, points_1), c2, r2)
    };
    let d = c2 - c1;
    let along = d.dot(&a1);
    if (d - a1 * along).length() > 1e-6 * scale || along.abs() < 1e-9 * scale {
        if trace {
            eprintln!(
                "  [bevel-cone] not concentric/coincident: perp={:.3e} along={along:.3e}",
                (d - a1 * along).length()
            );
        }
        return None;
    }
    let dr = r2 - r1;
    if dr.abs() < 1e-9 * scale {
        if trace {
            eprintln!("  [bevel-cone] cylindrical band dr={dr:.3e}");
        }
        return None; // cylindrical band — not this bevel's shape
    }
    // Apex where the radius extrapolates to zero along the axis.
    let t_apex = -r1 * along / dr;
    let apex = c1 + a1 * t_apex;
    // Axis points from the apex INTO the body (toward growing radius).
    let axis = if (along > 0.0) == (dr > 0.0) {
        a1
    } else {
        a1 * -1.0
    };
    let half_angle = (dr / along).abs().atan();
    Some((
        SurfaceGeometry::Cone {
            apex,
            axis,
            half_angle,
        },
        (c1, a1, r1),
        (c2, a2, r2),
    ))
}

/// Snap a sampled offset row onto its recognized exact circle.
fn snap_row_to_circle(row: &mut [Point3], circle: &BevelCircle) {
    let (c, axis, r) = circle;
    for p in row.iter_mut() {
        let d = *p - *c;
        let radial = d - *axis * d.dot(axis);
        if let Some(dir) = radial.normalized() {
            *p = *c + dir * *r;
        }
    }
}

/// Build clamped uniform knot vector for `n` control points at given degree.
fn build_clamped_uniform_knots(n: usize, degree: usize) -> Vec<f64> {
    let n_knots = n + degree + 1;
    let n_interior = n.saturating_sub(degree);
    let mut knots = Vec::with_capacity(n_knots);

    knots.extend(std::iter::repeat_n(0.0, degree + 1));
    for i in 1..n_interior {
        knots.push(i as f64 / n_interior as f64);
    }
    knots.extend(std::iter::repeat_n(1.0, degree + 1));

    knots
}

#[cfg(test)]
#[path = "chamfer_curve_tests.rs"]
mod tests;
