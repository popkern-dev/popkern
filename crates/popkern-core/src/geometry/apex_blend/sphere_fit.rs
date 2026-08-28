//! Sphere-fit termination math (sphere-fit substitution).
//!
//! Closed-form algebraic intersection of a cap sphere (centred at
//! the apex region offset, radius r_cap) with a fillet surface,
//! producing the curve where the sphere is tangent to the fillet —
//! the natural fillet termination for sphere-fit caps.
//!
//! # Surface pairs supported
//!
//! - Cylinder × Sphere — the cube/pyramid common case
//! - Plane × Sphere — for face-bridge curves on planar faces
//! - Torus × Sphere — for variable-radius fillets
//! - Cone × Sphere — for tapered edges
//!
//! All four are closed-form. the sphere-fit substitution work wires these through a
//! unified entry point [`sphere_fit_termination`].

use crate::geometry::analytic_curves::circle_basis;
use crate::geometry::nurbs_curve::NurbsCurve;
use crate::geometry::surface::SurfaceGeometry;
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;
use std::f64::consts::TAU;

/// Result of a sphere-fit termination computation.
#[derive(Clone, Debug)]
pub struct SphereFitTermination {
    /// 3-D curve on the fillet surface from one face-side
    /// endpoint to the other. Always a circle (or arc) for the
    /// surface pairs supported by this module.
    pub curve_3d: NurbsCurve,
    /// Pcurve in the fillet's UV space. Closed-form linear
    /// segment for cylinder; closed-form arc for torus/cone.
    pub pcurve_on_fillet: NurbsCurve,
    /// First endpoint (on adjacent face A).
    pub endpoint_a: Point3,
    /// Second endpoint (on adjacent face B).
    pub endpoint_b: Point3,
}

/// Failure modes from [`sphere_fit_termination`].
#[derive(Debug, thiserror::Error)]
pub enum SphereFitError {
    /// The cap sphere doesn't intersect the fillet surface, or
    /// the intersection is empty after restricting to the
    /// fillet's UV range.
    #[error("sphere does not intersect fillet surface in the fillet's UV range")]
    NoIntersection,
    /// The fillet surface variant isn't supported by the
    /// sphere-fit closed-form path. Caller should fall back to
    /// the iterative Gregory-cap termination math.
    #[error("sphere-fit termination not implemented for fillet variant {variant}")]
    UnsupportedFilletVariant {
        /// Static-string description of the unsupported variant.
        variant: &'static str,
    },
    /// Geometric degeneracy — sphere centre on the fillet
    /// surface, fillet surface degenerate, etc.
    #[error("geometric degeneracy: {reason}")]
    Degenerate {
        /// Static-string reason for downstream display.
        reason: &'static str,
    },
}

/// Compute the sphere-fit termination curve for a fillet.
///
/// `apex` is the original (pre-blend) apex vertex position.
/// `cap_offset` is the offset from the apex to the sphere
/// centre (typically `(r_cap, r_cap, r_cap) - apex` for a cube
/// corner with `r_cap = r`); `cap_radius` is the sphere radius.
///
/// The fillet surface, its UV range, and the two adjacent face
/// "sides" (A and B) bound the termination's endpoints on each
/// adjacent face.
///
/// **Endpoint A / B disambiguation:** the two intersection
/// points where the termination circle meets the fillet's
/// adjacent face boundaries. Caller specifies which face is
/// "A" and which is "B" via `face_a_normal` and `face_b_normal`
/// — endpoint A is the intersection point closer to the
/// half-plane the `face_a_normal` points away from.
#[allow(clippy::too_many_arguments)]
pub fn sphere_fit_termination(
    apex: Point3,
    cap_centre: Point3,
    cap_radius: f64,
    fillet_surf: &SurfaceGeometry,
    face_a_normal: Vec3,
    face_b_normal: Vec3,
    tol: &Tolerance,
) -> Result<SphereFitTermination, SphereFitError> {
    match fillet_surf {
        SurfaceGeometry::Cylinder {
            origin,
            axis,
            radius,
        } => cyl_sphere_termination(
            apex,
            cap_centre,
            cap_radius,
            *origin,
            *axis,
            *radius,
            face_a_normal,
            face_b_normal,
            tol,
        ),
        // Other variants will be implemented in subsequent
        // commits per the sphere-fit substitution work. Cube/pyramid fixtures only need
        // Cylinder×Sphere today.
        SurfaceGeometry::Plane { .. } => Err(SphereFitError::UnsupportedFilletVariant {
            variant: "Plane (the sphere-fit substitution work continuation)",
        }),
        SurfaceGeometry::Sphere { .. } => Err(SphereFitError::UnsupportedFilletVariant {
            variant: "Sphere (the sphere-fit substitution work continuation)",
        }),
        SurfaceGeometry::Torus { .. } => Err(SphereFitError::UnsupportedFilletVariant {
            variant: "Torus (the sphere-fit substitution work continuation)",
        }),
        SurfaceGeometry::Cone { .. } => Err(SphereFitError::UnsupportedFilletVariant {
            variant: "Cone (the sphere-fit substitution work continuation)",
        }),
        SurfaceGeometry::Nurbs(_) => Err(SphereFitError::UnsupportedFilletVariant {
            variant: "Nurbs (use Gregory-cap iterative path)",
        }),
        SurfaceGeometry::Gregory(_) => Err(SphereFitError::UnsupportedFilletVariant {
            variant: "Gregory (cap-of-cap not in sphere-fit scope)",
        }),
        SurfaceGeometry::PolyhedralSpline(_) => Err(SphereFitError::UnsupportedFilletVariant {
            variant: "PolyhedralSpline (cap-of-cap not in sphere-fit scope)",
        }),
    }
}

/// Cylinder × Sphere termination — closed form.
///
/// The termination is a CIRCLE on the cylinder (the cylinder's
/// cross-section at the parametric position where the sphere is
/// tangent). For each face adjacent to the fillet, one endpoint
/// of the termination arc is the intersection of this circle
/// with the face plane.
///
/// Construction:
/// 1. Sphere ∩ Cylinder = circle when sphere centre is on the
///    cylinder axis (axis-aligned case). For our cap-fit setup,
///    the cap sphere centre is constructed to lie on the
///    cylinder fillet's axis (sphere centred at apex offset
///    `(r,r,r)` for cube-corner; for asymmetric or non-axis
///    cases, fall back to the general path — TODO).
/// 2. The intersection circle is at parametric position `v_term`
///    along the cylinder axis, with radius equal to the
///    cylinder's radius.
/// 3. For each face adjacent to the fillet, the face plane
///    intersects the circle at one or two points. The endpoint
///    on face A is the intersection point closer to face A's
///    boundary half-plane.
#[allow(clippy::too_many_arguments)]
fn cyl_sphere_termination(
    _apex: Point3,
    cap_centre: Point3,
    cap_radius: f64,
    cyl_origin: Point3,
    cyl_axis: Vec3,
    cyl_radius: f64,
    face_a_normal: Vec3,
    face_b_normal: Vec3,
    _tol: &Tolerance,
) -> Result<SphereFitTermination, SphereFitError> {
    let cyl_axis_unit = cyl_axis.normalized().ok_or(SphereFitError::Degenerate {
        reason: "cylinder axis is zero-length",
    })?;

    // Find v_term: parametric position along cyl_axis where
    // sphere is tangent to cylinder. The cap sphere centre
    // projected onto the cylinder axis is at:
    //   v_proj = (cap_centre - cyl_origin) · cyl_axis_unit
    let v_proj = (cap_centre - cyl_origin).dot(&cyl_axis_unit);
    // Distance from cap_centre to the cylinder axis:
    let axis_pt = Point3::new(
        cyl_origin.x + cyl_axis_unit.x * v_proj,
        cyl_origin.y + cyl_axis_unit.y * v_proj,
        cyl_origin.z + cyl_axis_unit.z * v_proj,
    );
    let axial_offset = (cap_centre - axis_pt).length();

    // For sphere tangent to cylinder externally, the tangent
    // circle is at axial_offset = sqrt(cap_radius² − cyl_radius²)
    // ... wait, actually the geometry depends on whether the
    // sphere is tangent INTERNALLY or EXTERNALLY. For cube-
    // corner cap-fit:
    //   - Cap sphere of radius r centred at (r,r,r).
    //   - Cylinder fillet axis along x with origin (0,r,r),
    //     radius r.
    //   - Sphere centre to cyl axis: vec is (r,0,0), distance 0
    //     (since (r,0,0) is along x, parallel to axis, so the
    //     centre LIES ON the axis after projection). Wait —
    //     projection onto axis: v_proj = r. axis_pt = (r, r, r).
    //     axial_offset = ‖cap_centre − axis_pt‖ = ‖(0,0,0)‖ = 0.
    //
    // So sphere centre LIES ON the cylinder axis. Sphere ∩
    // cylinder is then the circle at v = v_proj on the cylinder
    // (when r_sphere == r_cyl), or empty (when r_sphere <
    // r_cyl), or two circles (when r_sphere > r_cyl).
    //
    // For r_sphere == r_cyl == r: termination is the cylinder
    // cross-section at v = v_proj.
    if axial_offset > _tol.linear() {
        return Err(SphereFitError::Degenerate {
            reason: "cap centre not on cylinder axis (general case TODO)",
        });
    }
    if (cap_radius - cyl_radius).abs() > _tol.linear() {
        return Err(SphereFitError::Degenerate {
            reason: "cap radius != cylinder radius (general case TODO)",
        });
    }

    // Termination is the cylinder cross-section circle at
    // v_proj. The circle's centre is `axis_pt`, radius
    // `cyl_radius`, in the plane perpendicular to `cyl_axis_unit`.
    let centre = axis_pt;

    // Build a basis for the circle: u_axis × v_axis perpendicular
    // to cyl_axis_unit.
    let (u_axis, v_axis) = circle_basis(cyl_axis_unit);

    // Find the two endpoints where the circle meets the face A
    // and face B planes. Face A plane passes through the
    // cylinder axis (since the fillet is tangent to the face
    // along the axis offset by cyl_radius). Specifically, the
    // face A plane normal is `face_a_normal` and contains the
    // origin of the cylinder.
    //
    // The face plane intersects the circle (centre, radius,
    // basis u/v) at angles theta where:
    //   centre + r·(cos(theta) u_axis + sin(theta) v_axis)
    // is on the face plane (i.e. perpendicular distance from
    // face plane to the point is zero).
    //
    // For face A normal n_a passing through cyl_origin:
    //   (centre + r·(cos·u + sin·v) - cyl_origin) · n_a = 0
    //   ⟹ (centre - cyl_origin)·n_a + r·(cos·u·n_a + sin·v·n_a) = 0
    //
    // Let:
    //   c0 = (centre - cyl_origin)·n_a / r
    //   cu = u·n_a
    //   cv = v·n_a
    //   cu·cos(theta) + cv·sin(theta) = -c0
    //
    // Solve theta via atan2 and amplitude.
    // Face planes pass through the apex (definitional for a vertex
    // blend — all faces incident to the apex contain it).
    let theta_a =
        solve_circle_plane_intersection(centre, _apex, cyl_radius, u_axis, v_axis, face_a_normal)?;
    let theta_b =
        solve_circle_plane_intersection(centre, _apex, cyl_radius, u_axis, v_axis, face_b_normal)?;

    let endpoint_a = circle_eval(centre, cyl_radius, u_axis, v_axis, theta_a);
    let endpoint_b = circle_eval(centre, cyl_radius, u_axis, v_axis, theta_b);

    // Build the 3-D arc as a NURBS rational quadratic Bezier
    // (closed form for a circular arc).
    let curve_3d = arc_as_nurbs(centre, cyl_radius, u_axis, v_axis, theta_a, theta_b);

    // Build the pcurve in cylinder UV. For Cylinder, UV is
    // (angle, height-along-axis). The termination is at fixed
    // height v_proj, with angle going from theta_a to theta_b.
    // It's a horizontal line in UV. As a degree-1 NURBS:
    let pcurve_on_fillet = pcurve_horizontal_line(theta_a, theta_b, v_proj);

    Ok(SphereFitTermination {
        curve_3d,
        pcurve_on_fillet,
        endpoint_a,
        endpoint_b,
    })
}

/// Solve `cu·cos(θ) + cv·sin(θ) = -c0` for `θ`. The angle is
/// returned in `[0, TAU)`. Returns `NoIntersection` when
/// `cu² + cv² < c0²` (circle and plane don't meet).
fn solve_circle_plane_intersection(
    centre: Point3,
    origin: Point3,
    radius: f64,
    u_axis: Vec3,
    v_axis: Vec3,
    plane_normal: Vec3,
) -> Result<f64, SphereFitError> {
    let n = plane_normal
        .normalized()
        .ok_or(SphereFitError::Degenerate {
            reason: "face plane normal is zero-length",
        })?;
    let c0 = (centre - origin).dot(&n) / radius;
    let cu = u_axis.dot(&n);
    let cv = v_axis.dot(&n);
    let amplitude = (cu * cu + cv * cv).sqrt();
    if amplitude < 1e-300 {
        return Err(SphereFitError::Degenerate {
            reason: "plane parallel to circle",
        });
    }
    if c0.abs() > amplitude + 1e-12 {
        // POPKERN_SPHERE_FIT_TRACE: the reachability numbers behind a
        // NoIntersection verdict — |c0| is the plane's normalized
        // offset from the circle centre, amplitude the circle's reach
        // toward the plane.
        if std::env::var_os("POPKERN_SPHERE_FIT_TRACE").is_some() {
            eprintln!(
                "[sphere-fit-trace] circle∩plane MISS: c0={c0:.6} amplitude={amplitude:.6} \
                 centre=({:.4},{:.4},{:.4}) r={radius:.4} n=({:.4},{:.4},{:.4})",
                centre.x, centre.y, centre.z, n.x, n.y, n.z
            );
        }
        return Err(SphereFitError::NoIntersection);
    }
    // Identity: cu·cos(θ) + cv·sin(θ) = R·cos(θ − φ) where
    // R = sqrt(cu² + cv²) and φ = atan2(cv, cu). Solving for θ:
    //   R·cos(θ − φ) = −c0  ⟹  θ = φ ± acos(−c0 / R).
    let phi = cv.atan2(cu);
    let arg = (-c0 / amplitude).clamp(-1.0, 1.0);
    let alpha = arg.acos();
    let theta1 = phi + alpha;
    let theta2 = phi - alpha;
    // Pick the solution in `[0, TAU)`.
    let theta1 = (theta1.rem_euclid(TAU)).rem_euclid(TAU);
    let theta2 = (theta2.rem_euclid(TAU)).rem_euclid(TAU);
    // Disambiguation: pick the smaller angle (closer to theta=0).
    // For cube corner this gives the consistent endpoint.
    Ok(if theta1 <= theta2 { theta1 } else { theta2 })
}

fn circle_eval(centre: Point3, radius: f64, u_axis: Vec3, v_axis: Vec3, theta: f64) -> Point3 {
    let c = theta.cos();
    let s = theta.sin();
    Point3::new(
        centre.x + radius * (c * u_axis.x + s * v_axis.x),
        centre.y + radius * (c * u_axis.y + s * v_axis.y),
        centre.z + radius * (c * u_axis.z + s * v_axis.z),
    )
}

/// Build a circular arc as a degree-2 rational NURBS — the
/// standard exact-on-circle representation. CPs: P0 = arc start,
/// P1 = the "shoulder" point at the intersection of tangent
/// lines at P0 and P2, P2 = arc end. Weights: [1, cos α, 1]
/// where α is the half-angle of the arc.
///
/// For arcs with span ≥ π, splits into multiple pieces (single
/// rational quadratic only spans < π exactly). For our use
/// case (cap-fillet termination), the arc is ≤ π/2 (a single
/// quarter-circle), so this single-piece form suffices.
fn arc_as_nurbs(
    centre: Point3,
    radius: f64,
    u_axis: Vec3,
    v_axis: Vec3,
    theta_a: f64,
    theta_b: f64,
) -> NurbsCurve {
    // Normalise span to [-π, π] so we pick the short arc by default.
    let mut span = theta_b - theta_a;
    while span > std::f64::consts::PI {
        span -= TAU;
    }
    while span < -std::f64::consts::PI {
        span += TAU;
    }
    let theta_m = theta_a + 0.5 * span;
    let half_angle = 0.5 * span.abs();
    let p0 = circle_eval(centre, radius, u_axis, v_axis, theta_a);
    let p2 = circle_eval(centre, radius, u_axis, v_axis, theta_b);
    // P1: intersection of tangent at P0 and tangent at P2. By
    // standard rational-quadratic-on-circle construction, this
    // is the point at angle θ_m on a circle of radius r/cos(α)
    // (the "shoulder" radius).
    let shoulder_radius = if half_angle.cos().abs() > 1e-300 {
        radius / half_angle.cos()
    } else {
        radius
    };
    let p1 = circle_eval(centre, shoulder_radius, u_axis, v_axis, theta_m);
    NurbsCurve {
        degree: 2,
        control_points: vec![p0, p1, p2],
        weights: vec![1.0, half_angle.cos(), 1.0],
        knots: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
    }
}

fn pcurve_horizontal_line(theta_a: f64, theta_b: f64, v: f64) -> NurbsCurve {
    NurbsCurve {
        degree: 1,
        control_points: vec![Point3::new(theta_a, v, 0.0), Point3::new(theta_b, v, 0.0)],
        weights: vec![1.0, 1.0],
        knots: vec![0.0, 0.0, 1.0, 1.0],
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::curve::CurveEval;

    /// Symmetric cube corner with apex at origin, cube extending
    /// into +x +y +z, equal radius r=0.1 on all three apex-
    /// incident edges.
    ///
    /// Test the x-edge fillet. Expected geometry:
    /// - Cylinder axis along x, origin (0, r, r), radius r
    /// - Cap sphere centred at (r, r, r), radius r
    /// - Termination is the cylinder cross-section at x=r
    /// - Endpoints: (r, 0, r) on y-face (y=0) and (r, r, 0) on
    ///   z-face (z=0)
    #[test]
    fn cube_corner_symmetric_x_edge() {
        let r = 0.1;
        let apex = Point3::new(0.0, 0.0, 0.0);
        let cap_centre = Point3::new(r, r, r);
        let cap_radius = r;

        // x-edge fillet
        let cyl = SurfaceGeometry::Cylinder {
            origin: Point3::new(0.0, r, r),
            axis: Vec3::new(1.0, 0.0, 0.0),
            radius: r,
        };
        // Adjacent faces: y-face (y=0, normal +y outward, but
        // for fillet on inside we use -y or just the plane normal)
        let face_a_normal = Vec3::new(0.0, 1.0, 0.0); // y-face
        let face_b_normal = Vec3::new(0.0, 0.0, 1.0); // z-face

        let result = sphere_fit_termination(
            apex,
            cap_centre,
            cap_radius,
            &cyl,
            face_a_normal,
            face_b_normal,
            &Tolerance::default(),
        )
        .expect("sphere fit termination should succeed for symmetric cube");

        // Check endpoints — should be (r, 0, r) and (r, r, 0).
        let expected_a = Point3::new(r, 0.0, r);
        let expected_b = Point3::new(r, r, 0.0);

        let d_a_to_expected_a = result.endpoint_a.distance_to(&expected_a);
        let d_a_to_expected_b = result.endpoint_a.distance_to(&expected_b);
        let d_b_to_expected_a = result.endpoint_b.distance_to(&expected_a);
        let d_b_to_expected_b = result.endpoint_b.distance_to(&expected_b);

        // Either (a, b) or (b, a) — order may differ.
        let ordered = (d_a_to_expected_a < 1e-12 && d_b_to_expected_b < 1e-12)
            || (d_a_to_expected_b < 1e-12 && d_b_to_expected_a < 1e-12);
        assert!(
            ordered,
            "endpoints don't match: got a={:?}, b={:?}; expected one of [(r,0,r), (r,r,0)] in some order",
            result.endpoint_a, result.endpoint_b
        );

        // Curve is on the cylinder: every point should be at
        // distance r from the axis line {(t, r, r) : t ∈ R}.
        for i in 0..=10 {
            let t = (i as f64) / 10.0;
            let p = result.curve_3d.point_at(t);
            // Project p onto cylinder axis: axis is x-direction
            // with origin (0, r, r). Distance from p to axis is
            // sqrt((p.y - r)^2 + (p.z - r)^2).
            let d_to_axis = ((p.y - r).powi(2) + (p.z - r).powi(2)).sqrt();
            assert!(
                (d_to_axis - r).abs() < 1e-12,
                "curve point at t={t} is {d_to_axis} from cylinder axis (expected {r}); \
                 rational quadratic should be exact"
            );
        }
    }
}
