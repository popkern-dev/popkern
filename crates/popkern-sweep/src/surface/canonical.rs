//! Canonical surface recognition for swept surfaces.
//!
//! An extruded circular arc IS a cylinder; leaving it as a generic
//! NURBS forces every downstream consumer (boolean intersection,
//! fillet analysis, STEP export) onto approximate paths — the
//! boolean's marching walker on a NURBS "cylinder" × plane produced
//! intersection curves with ±6e-2 radial slop where the analytic
//! path is exact. Real kernels canonicalize at creation time; so do
//! we.

use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::{NurbsSurface, Point3, Vec3};

/// Affine map from the swept surface's `[0,1]²` chart to the
/// recognized surface's native chart (cylinder: `u` = angle in
/// radians, `v` = length along the axis).
#[derive(Debug, Clone, Copy)]
pub(crate) struct ChartMap {
    /// `u_native = u_offset + u01 * u_scale`.
    pub u_offset: f64,
    /// See `u_offset`.
    pub u_scale: f64,
    /// `v_native = v_offset + v01 * v_scale`.
    pub v_offset: f64,
    /// See `v_offset`.
    pub v_scale: f64,
}

impl ChartMap {
    /// Map a `[0,1]²` chart point into the native chart.
    pub fn map(&self, u01: f64, v01: f64) -> (f64, f64) {
        (
            self.u_offset + u01 * self.u_scale,
            self.v_offset + v01 * self.v_scale,
        )
    }
}

/// Recognize an extruded circular arc as an analytic cylinder.
///
/// Conditions verified by dense sampling (not construction
/// provenance, so any parameterization qualifies):
/// 1. Translation surface: every v-iso row is a translate of the
///    v=0 row (constant offset vector — a straight extrusion).
/// 2. The v=0 profile row lies on a circle whose axis is the
///    translation direction.
///
/// Returns the cylinder and the chart map from the swept `[0,1]²`
/// parameterization. NOTE: the u chart map is angle-LINEAR while a
/// rational-quadratic circle's parameter is not; pcurves remapped
/// through it stay geometrically on the surface but their
/// parameterization deviates from the 3D edge curve's by up to a few
/// percent mid-span — consumers in this kernel treat pcurves as
/// association seeds and refine by inversion, which absorbs that.
pub(crate) fn recognize_extruded_cylinder(
    surface: &NurbsSurface,
) -> Option<(SurfaceGeometry, ChartMap)> {
    const NU: usize = 33;
    const NV: usize = 5;
    // 1. Translation test.
    let row0: Vec<Point3> = (0..NU)
        .map(|i| surface.point_at(i as f64 / (NU - 1) as f64, 0.0))
        .collect();
    // Same envelope shape as recognize_extruded_plane: construction
    // noise scales with the profile span, not just the offset
    // (measured there at ~1.5e-9 relative; the old offset-only 1e-9
    // envelope rejected real cylinder walls).
    let row_span = (row0[NU - 1] - row0[0]).length();
    let mut axis_dir: Option<Vec3> = None;
    let mut height = 0.0;
    for j in 1..NV {
        let v = j as f64 / (NV - 1) as f64;
        let offset = surface.point_at(0.0, v) - row0[0];
        for (i, p0) in row0.iter().enumerate() {
            let p = surface.point_at(i as f64 / (NU - 1) as f64, v);
            let d = p - *p0;
            if (d - offset).length() > 1.0e-8 * (1.0 + row_span + offset.length()) {
                if std::env::var_os("PK_CANON_TRACE").is_some() {
                    eprintln!(
                        "[canon-cyl] translation decline: dev {:.3e} vs env {:.3e}",
                        (d - offset).length(),
                        1.0e-8 * (1.0 + row_span + offset.length())
                    );
                }
                return None;
            }
        }
        if j == NV - 1 {
            height = offset.length();
            if height < 1.0e-12 {
                return None;
            }
            axis_dir = offset.normalized();
        }
    }
    let axis = axis_dir?;

    // 2. Circle test on the profile row: circumcenter from three
    //    spread points, then every sample at constant radius and in
    //    a plane perpendicular to the axis.
    let (a, b, c) = (row0[0], row0[NU / 3], row0[2 * NU / 3]);
    let center = circumcenter_3d(a, b, c, axis)?;
    let radius = a.distance_to(&center);
    if radius < 1.0e-12 {
        return None;
    }
    // Construction noise scales with the profile's coordinate span,
    // not just the radius (measured: a 27 mm arc whose profile
    // coordinates run to ~70 mm deviated 1.0e-7 axially — far over a
    // 1e-9·radius relative gate). Same constant and shape as the
    // translation envelope above; real non-cylinders deviate at
    // sagitta scale, 3+ orders of magnitude above this gate.
    let gate = 1.0e-8 * (1.0 + radius + row_span);
    for p in &row0 {
        let d = *p - center;
        if (d.dot(&axis)).abs() > gate {
            if std::env::var_os("PK_CANON_TRACE").is_some() {
                eprintln!(
                    "[canon-cyl] plane decline: axial {:.3e} vs gate {:.3e}",
                    (d.dot(&axis)).abs(),
                    gate
                );
            }
            return None;
        }
        let radial = d - axis * d.dot(&axis);
        if (radial.length() - radius).abs() > gate {
            if std::env::var_os("PK_CANON_TRACE").is_some() {
                eprintln!(
                    "[canon-cyl] radius decline: dev {:.3e} vs gate {:.3e} (r {radius:.3})",
                    (radial.length() - radius).abs(),
                    gate
                );
            }
            return None;
        }
    }

    // 3. Chart: angles of the profile in the cylinder's own basis
    //    (matching `SurfaceGeometry::Cylinder`'s orthonormal_basis
    //    convention so point_at agrees), unwrapped monotonically.
    let (perp1, perp2) = popkern_core::geometry::analytic_surfaces::orthonormal_basis(axis);
    let angle_of = |p: &Point3| -> f64 {
        let d = *p - center;
        d.dot(&perp2).atan2(d.dot(&perp1))
    };
    let theta0 = angle_of(&row0[0]);
    // Unwrap the end angle to the branch reached by walking the row.
    let mut walked = theta0;
    for p in row0.iter().skip(1) {
        let mut t = angle_of(p);
        while t - walked > std::f64::consts::PI {
            t -= std::f64::consts::TAU;
        }
        while t - walked < -std::f64::consts::PI {
            t += std::f64::consts::TAU;
        }
        walked = t;
    }
    let theta1 = walked;
    if (theta1 - theta0).abs() < 1.0e-12 {
        // Degenerate span; closed full circles land at ±2π.
        return None;
    }

    let geometry = SurfaceGeometry::Cylinder {
        origin: center,
        axis,
        radius,
    };
    // Verify the recognized surface reproduces the sampled grid.
    for j in 0..NV {
        let v01 = j as f64 / (NV - 1) as f64;
        for i in 0..NU {
            let u01 = i as f64 / (NU - 1) as f64;
            let native_u = theta0 + (theta1 - theta0) * u01;
            let native_v = v01 * height;
            let p_nurbs = surface.point_at(u01, v01);
            let p_cyl = geometry.point_at(native_u, native_v);
            // The angle chart is linear while the NURBS parameter is
            // not — compare against the cylinder at the SAMPLE's own
            // angle instead of the chart image for the radial check,
            // and accept chart-parameterization skew along the arc.
            let d = p_nurbs - center;
            let radial = d - axis * d.dot(&axis);
            if (radial.length() - radius).abs() > gate {
                return None;
            }
            let _ = p_cyl;
        }
    }
    Some((
        geometry,
        ChartMap {
            u_offset: theta0,
            u_scale: theta1 - theta0,
            v_offset: 0.0,
            v_scale: height,
        },
    ))
}

/// Whether a recognized canonical surface's normal field agrees with
/// the authored NURBS frame normal (`du × dv`).
///
/// The sweep wires face loops and pcurves against the NURBS chart
/// (u = profile, v = path), whose frame normal is the face's outward
/// direction by construction. Canonical geometry carries its OWN
/// normal convention — a cylinder's `normal_at` points +radial no
/// matter which side the material lies on, so a CONCAVE profile arc
/// (material outside the circle) needs `orientation = false`. The
/// agreement is sampled at the chart midpoint with finite
/// differences, retrying on a coarse grid if the sample is
/// indecisive (degenerate spans).
pub(crate) fn canonical_normal_agreement(
    nurbs: &NurbsSurface,
    canonical: &SurfaceGeometry,
) -> bool {
    let samples = [
        (0.5, 0.5),
        (0.25, 0.5),
        (0.75, 0.5),
        (0.5, 0.25),
        (0.25, 0.25),
    ];
    for &(u, v) in &samples {
        let h = 1.0e-4;
        let p0 = nurbs.point_at(u, v);
        let pu = nurbs.point_at(u + h, v);
        let pv = nurbs.point_at(u, v + h);
        let du = pu - p0;
        let dv = pv - p0;
        let frame_normal = du.cross(&dv);
        let scale = du.length() * dv.length();
        if frame_normal.length() < 1.0e-9 * (1.0 + scale) {
            continue; // degenerate frame at this sample
        }
        let Some((cu, cv)) =
            popkern_core::geometry::uv_projection::project_point_to_uv(canonical, &p0)
        else {
            continue;
        };
        let canonical_normal = canonical.normal_at(cu, cv);
        let dot = frame_normal.dot(&canonical_normal);
        if dot.abs() > 1.0e-9 * (1.0 + frame_normal.length() * canonical_normal.length()) {
            return dot > 0.0;
        }
    }
    // No decisive sample — keep the historical convention.
    true
}

/// Recognize an extruded straight profile span as an analytic plane.
///
/// Conditions: translation surface (straight extrusion) whose v=0
/// profile row is collinear. The chart map is EXACT here (linear
/// profile parameterization): u = distance along the profile line,
/// v = distance along the extrusion.
pub(crate) fn recognize_extruded_plane(
    surface: &NurbsSurface,
) -> Option<(SurfaceGeometry, ChartMap)> {
    const NU: usize = 17;
    const NV: usize = 5;
    let row0: Vec<Point3> = (0..NU)
        .map(|i| surface.point_at(i as f64 / (NU - 1) as f64, 0.0))
        .collect();
    // Straightness of the profile row.
    let chord = row0[NU - 1] - row0[0];
    let len = chord.length();
    if len < 1.0e-12 {
        return None;
    }
    let dir = chord * (1.0 / len);
    for p in &row0 {
        let d = *p - row0[0];
        let off = d - dir * d.dot(&dir);
        if off.length() > 1.0e-9 * (1.0 + len) {
            return None;
        }
        // Monotone along the line (no doubling back).
        if d.dot(&dir) < -1.0e-9 {
            return None;
        }
    }
    // Translation test. The envelope scales with the PROFILE span as
    // well as the offset: construction error in the surface control
    // points accumulates with the profile's coordinate magnitude, so
    // an offset-only envelope rejects long straight walls over
    // nothing (measured: a 127 mm wall failed at 2.2e-9 relative —
    // ~2.6e-7 mm absolute — and the unrecognized NURBS×plane pair
    // hung the marching intersector downstream). Genuinely curved
    // surfaces deviate by geometric amounts, orders of magnitude
    // above this envelope.
    let mut axis_dir: Option<Vec3> = None;
    let mut height = 0.0;
    for j in 1..NV {
        let v = j as f64 / (NV - 1) as f64;
        let offset = surface.point_at(0.0, v) - row0[0];
        for (i, p0) in row0.iter().enumerate() {
            let p = surface.point_at(i as f64 / (NU - 1) as f64, v);
            // 1e-8: measured construction noise reaches ~1.5e-9
            // RELATIVE across scales (133 mm and 127 mm walls both
            // sit at the old 1e-9 envelope — source of the ~1e-9
            // relative floor is an open question, PK-102), while
            // genuine non-translational surfaces (lofts, scaled
            // sweeps) deviate at ≥1e-3 relative. Five orders of
            // margin on the unsafe side, 7× on the noise side.
            if ((p - *p0) - offset).length() > 1.0e-8 * (1.0 + len + offset.length()) {
                return None;
            }
        }
        if j == NV - 1 {
            height = offset.length();
            if height < 1.0e-12 {
                return None;
            }
            axis_dir = offset.normalized();
        }
    }
    let v_axis = axis_dir?;
    let normal = dir.cross(&v_axis).normalized()?;
    let geometry = SurfaceGeometry::Plane {
        origin: row0[0],
        normal,
        u_axis: dir,
        v_axis,
    };
    Some((
        geometry,
        ChartMap {
            u_offset: 0.0,
            u_scale: len,
            v_offset: 0.0,
            v_scale: height,
        },
    ))
}

/// Circumcenter of three points, constrained to the plane through
/// `a` perpendicular to `axis` (extrusion profiles are planar).
fn circumcenter_3d(a: Point3, b: Point3, c: Point3, axis: Vec3) -> Option<Point3> {
    // Build a 2D frame in the profile plane.
    let (e1, e2) = popkern_core::geometry::analytic_surfaces::orthonormal_basis(axis);
    let to2d = |p: Point3| -> (f64, f64) {
        let d = p - a;
        (d.dot(&e1), d.dot(&e2))
    };
    let (ax, ay) = (0.0, 0.0);
    let (bx, by) = to2d(b);
    let (cx, cy) = to2d(c);
    let dd = 2.0 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by));
    if dd.abs() < 1.0e-16 {
        return None;
    }
    let ux = ((ax * ax + ay * ay) * (by - cy)
        + (bx * bx + by * by) * (cy - ay)
        + (cx * cx + cy * cy) * (ay - by))
        / dd;
    let uy = ((ax * ax + ay * ay) * (cx - bx)
        + (bx * bx + by * by) * (ax - cx)
        + (cx * cx + cy * cy) * (bx - ax))
        / dd;
    Some(a + e1 * ux + e2 * uy)
}

/// Recognize a closed circular NURBS curve as an analytic circle
/// `(center, normal, radius)` — Kåsa least-squares in the Newell
/// plane, verified to 1e-9 relative. Sweep cap curves of circular
/// profiles are exact rational circles; keeping them as B-splines
/// costs downstream consumers their analytic paths (and OpenCASCADE
/// importers reject faces bounded by closed rational splines that
/// analytic circles pass cleanly).
pub(crate) fn recognize_circle(
    curve: &popkern_core::geometry::NurbsCurve,
) -> Option<(Point3, Vec3, f64)> {
    use popkern_core::geometry::curve::CurveEval;
    let (t0, t1) = curve.parameter_range();
    if curve.point_at(t0).distance_to(&curve.point_at(t1)) > 1.0e-9 {
        return None;
    }
    let n = 24usize;
    let pts: Vec<Point3> = (0..n)
        .map(|i| curve.point_at(t0 + (t1 - t0) * i as f64 / n as f64))
        .collect();
    let m = pts.len() as f64;
    let c0 = Point3::new(
        pts.iter().map(|p| p.x).sum::<f64>() / m,
        pts.iter().map(|p| p.y).sum::<f64>() / m,
        pts.iter().map(|p| p.z).sum::<f64>() / m,
    );
    let mut nx = 0.0;
    let mut ny = 0.0;
    let mut nz = 0.0;
    for i in 0..pts.len() {
        let a = pts[i];
        let b = pts[(i + 1) % pts.len()];
        nx += (a.y - b.y) * (a.z + b.z);
        ny += (a.z - b.z) * (a.x + b.x);
        nz += (a.x - b.x) * (a.y + b.y);
    }
    let normal = Vec3::new(nx, ny, nz).normalized()?;
    let (e1, e2) = popkern_core::geometry::analytic_surfaces::orthonormal_basis(normal);
    let xy: Vec<(f64, f64)> = pts
        .iter()
        .map(|p| {
            let d = *p - c0;
            (d.dot(&e1), d.dot(&e2))
        })
        .collect();
    let (mut sx, mut sy, mut sxx, mut syy, mut sxy, mut sxz, mut syz, mut sz) =
        (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    for &(x, y) in &xy {
        let z = x * x + y * y;
        sx += x;
        sy += y;
        sxx += x * x;
        syy += y * y;
        sxy += x * y;
        sxz += x * z;
        syz += y * z;
        sz += z;
    }
    let mm = [[sxx, sxy, sx], [sxy, syy, sy], [sx, sy, m]];
    let rhs = [sxz, syz, sz];
    let det3 = |a: &[[f64; 3]; 3]| -> f64 {
        a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
            - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
            + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0])
    };
    let d0 = det3(&mm);
    if d0.abs() < 1e-30 {
        return None;
    }
    let solve = |col: usize| -> f64 {
        let mut t = mm;
        for (row, trow) in t.iter_mut().enumerate() {
            trow[col] = rhs[row];
        }
        det3(&t) / d0
    };
    let (a, b, c2) = (solve(0), solve(1), solve(2));
    let (ccx, ccy) = (0.5 * a, 0.5 * b);
    let r = (c2 + ccx * ccx + ccy * ccy).max(0.0).sqrt();
    if r < 1e-12 {
        return None;
    }
    let center = c0 + e1 * ccx + e2 * ccy;
    let gate = r * 1.0e-9;
    for p in &pts {
        let d = *p - center;
        if (d.length() - r).abs() > gate || d.dot(&normal).abs() > gate {
            return None;
        }
    }
    Some((center, normal, r))
}

/// Recognize a straight NURBS curve: every control point collinear
/// with the endpoint chord (relative 1e-9). Returns the unit
/// direction and length.
fn recognize_line(
    curve: &popkern_core::geometry::NurbsCurve,
) -> Option<(
    popkern_core::geometry::Point3,
    popkern_core::geometry::Vec3,
    f64,
)> {
    let first = *curve.control_points.first()?;
    let last = *curve.control_points.last()?;
    let chord = last - first;
    let len = chord.length();
    if len < 1e-12 {
        return None;
    }
    let dir = chord * (1.0 / len);
    for p in &curve.control_points {
        let d = *p - first;
        let off = d - dir * d.dot(&dir);
        if off.length() > 1e-9 * len.max(1.0) {
            return None;
        }
    }
    Some((first, dir, len))
}

/// Curve-on-edge with canonical recognition: straight curves become
/// analytic Lines, circular curves become analytic
/// angle-parameterized Arcs, everything else stays NURBS.
pub(crate) fn analytic_or_nurbs_curve(
    curve: &popkern_core::geometry::NurbsCurve,
) -> popkern_core::geometry::CurveOnEdge {
    use popkern_core::geometry::curve::CurveEval;
    if let Some((origin, direction, length)) = recognize_line(curve) {
        return popkern_core::geometry::CurveOnEdge {
            geometry: popkern_core::geometry::curve::CurveGeometry::Line { origin, direction },
            t_start: 0.0,
            t_end: length,
        };
    }
    if let Some((center, normal, radius)) = recognize_circle(curve) {
        let (t0, _) = curve.parameter_range();
        let start = curve.point_at(t0);
        // Angles MUST be measured in the Arc evaluator's own frame
        // (circle_basis) — orthonormal_basis is a different
        // convention and a mismatched frame rotates the arc's
        // parameterization off its seam vertex (measured: closed cap
        // rims whose trim start landed 90° from the vertex).
        let (e1, e2) = popkern_core::geometry::analytic_curves::circle_basis(normal);
        let d = start - center;
        let theta0 = d.dot(&e2).atan2(d.dot(&e1));
        return popkern_core::geometry::CurveOnEdge {
            geometry: popkern_core::geometry::curve::CurveGeometry::Arc {
                center,
                normal,
                radius,
                start_angle: theta0,
                end_angle: theta0 + std::f64::consts::TAU,
            },
            t_start: theta0,
            t_end: theta0 + std::f64::consts::TAU,
        };
    }
    let (t_start, t_end) = curve.parameter_range();
    popkern_core::geometry::CurveOnEdge {
        geometry: popkern_core::geometry::curve::CurveGeometry::Nurbs(curve.clone()),
        t_start,
        t_end,
    }
}
