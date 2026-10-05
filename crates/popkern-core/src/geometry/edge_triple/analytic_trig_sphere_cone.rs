//! Sphere × cone axial dispatch.
//!
//! When the cone axis is parallel to the sphere's polar axis
//! `(0, 0, 1)` and the sphere centre lies on the cone axis, the
//! intersection is exactly 1 or 2 circles — solutions of the
//! quadratic `v²·sec²(θ) − 2dv + (d² − r²) = 0` where `d` is the
//! signed distance from cone apex to sphere centre and `θ` is the
//! cone half-angle.
//!
//! Output:
//! - `curve_3d`: exact `CurveGeometry::Arc`.
//! - sphere pcurve: `Pcurve::Linear` (constant colatitude, linear
//!   longitude — guaranteed because the cone axis equals sphere's
//!   polar axis).
//! - cone pcurve: `Pcurve::Linear` (constant v at the intersection
//!   height, linear u).
//! - `ParameterMap::Identity` via shared θ.

use crate::geometry::curve::CurveGeometry;
use crate::geometry::edge_triple::error::EdgeTripleError;
use crate::geometry::edge_triple::inputs::{EdgeAnchors, SurfaceSide};
use crate::geometry::edge_triple::types::{
    AnalyticPair, BranchSelector, EdgeTopology, EdgeTriple, EdgeTripleProvenance,
    EdgeTripleQuality, ParameterMap,
};
use crate::geometry::pcurve::Pcurve;
use crate::geometry::surface::SurfaceGeometry;
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;
use popkern_math::PortableMath;

#[allow(clippy::too_many_arguments)]
pub(super) fn try_sphere_cone_axial(
    sphere_center: Point3,
    sphere_radius: f64,
    sphere_is_b: bool,
    cone_apex: Point3,
    cone_axis: Vec3,
    half_angle: f64,
    anchors: &EdgeAnchors,
    tol: &Tolerance,
) -> Option<Result<EdgeTriple, EdgeTripleError>> {
    let polar = Vec3::new(0.0, 0.0, 1.0);
    // Cone axis must match the sphere's polar axis (within tol).
    if cone_axis.dot(&polar).abs() < 1.0 - tol.angular() {
        return None;
    }
    // Sphere centre must lie on the cone axis.
    let cap = sphere_center - cone_apex;
    let perp = cap - cone_axis * cap.dot(&cone_axis);
    if perp.length() > tol.linear() {
        return None;
    }

    let axis_sign = cone_axis.dot(&polar).signum();
    let d = cap.dot(&cone_axis); // signed offset apex → sphere centre, along cone axis
    let cos_t = half_angle.pcos();
    let sin_t = half_angle.psin();
    let r = sphere_radius;
    // Quadratic: v² − 2d·cos²θ·v + (d² − r²)·cos²θ = 0
    // discriminant / 4 = cos²θ · (r² − d²·sin²θ)
    let disc_inner = r * r - d * d * sin_t * sin_t;
    if disc_inner < tol.linear() {
        return None; // tangent or non-intersecting
    }
    let sqrt_disc = cos_t * disc_inner.sqrt();
    let v_outer = d * cos_t * cos_t + sqrt_disc;
    let v_inner = d * cos_t * cos_t - sqrt_disc;

    // Pick branch by anchor radial distance.
    let (start, end) = match anchors {
        EdgeAnchors::Open { start, end } => (*start, *end),
        _ => return None,
    };
    let radial_at = |p: Point3| -> f64 {
        let dp = p - cone_apex;
        let proj = dp - cone_axis * dp.dot(&cone_axis);
        proj.length()
    };
    let r_outer = v_outer * half_angle.ptan();
    let r_inner = v_inner * half_angle.ptan();
    let r_at_start = radial_at(start);
    let (v_branch, circle_radius) = if (r_at_start - r_outer).abs() <= (r_at_start - r_inner).abs()
    {
        (v_outer, r_outer)
    } else {
        (v_inner, r_inner)
    };
    if v_branch <= tol.linear() || circle_radius <= tol.linear() {
        return None; // branch lies on the negative-v half (other nappe)
    }

    let circle_center = cone_apex + cone_axis * v_branch;
    // colat on sphere: z-component of (circle_center − sphere_center)
    // divided by sphere radius. Result is in [0, π].
    let z_offset = (circle_center.z - sphere_center.z) / sphere_radius;
    let colat = z_offset.clamp(-1.0, 1.0).pacos();
    // axis_sign sanity: surface evaluation only depends on z;
    // non-positive axis_sign means cone is upside-down relative to
    // popkern's sphere v=0 == +z convention. Skip for safety.
    let _ = axis_sign;

    let (cone_perp_a, cone_perp_b) =
        crate::geometry::analytic_surfaces::orthonormal_basis(cone_axis);
    let (arc_basis_u, arc_basis_v) = crate::geometry::analytic_curves::circle_basis(cone_axis);
    let cone_u_offset = arc_basis_u
        .dot(&cone_perp_b)
        .patan2(arc_basis_u.dot(&cone_perp_a));
    // Sphere u-coordinate = atan2(y, x). Sphere uses global (x, y)
    // for u, equivalent to using `(1, 0, 0), (0, 1, 0)` basis.
    let global_x = Vec3::new(1.0, 0.0, 0.0);
    let global_y = Vec3::new(0.0, 1.0, 0.0);
    let sphere_u_offset = arc_basis_u
        .dot(&global_y)
        .patan2(arc_basis_u.dot(&global_x));

    let theta_of = |p: Point3, side: SurfaceSide| -> Result<f64, EdgeTripleError> {
        let dp = p - circle_center;
        let cu = dp.dot(&arc_basis_u);
        let cv = dp.dot(&arc_basis_v);
        let r = (cu * cu + cv * cv).sqrt();
        if (r - circle_radius).abs() > tol.linear() {
            return Err(EdgeTripleError::AnchorOffSurface {
                which: side,
                distance: r - circle_radius,
            });
        }
        Ok(cv.patan2(cu))
    };
    let theta_start = match theta_of(start, SurfaceSide::A) {
        Ok(t) => t,
        Err(e) => return Some(Err(e)),
    };
    let theta_end = match theta_of(end, SurfaceSide::B) {
        Ok(t) => t,
        Err(e) => return Some(Err(e)),
    };
    let theta_end = if theta_end < theta_start {
        theta_end + std::f64::consts::TAU
    } else {
        theta_end
    };

    let curve_3d = CurveGeometry::Arc {
        center: circle_center,
        normal: cone_axis,
        radius: circle_radius,
        start_angle: theta_start,
        end_angle: theta_end,
    };

    let pcurve_sphere = Pcurve::linear(
        (theta_start + sphere_u_offset, colat),
        (theta_end + sphere_u_offset, colat),
        theta_start,
        theta_end,
    );
    let pcurve_cone = Pcurve::linear(
        (theta_start + cone_u_offset, v_branch),
        (theta_end + cone_u_offset, v_branch),
        theta_start,
        theta_end,
    );

    let sphere_surf = SurfaceGeometry::Sphere {
        center: sphere_center,
        radius: sphere_radius,
    };
    let cone_surf = SurfaceGeometry::Cone {
        apex: cone_apex,
        axis: cone_axis,
        half_angle,
    };

    let (pcurve_a, pcurve_b, surf_a, surf_b) = if sphere_is_b {
        (pcurve_cone, pcurve_sphere, &cone_surf, &sphere_surf)
    } else {
        (pcurve_sphere, pcurve_cone, &sphere_surf, &cone_surf)
    };

    let triple = EdgeTriple {
        curve_3d,
        pcurve_a,
        pcurve_b,
        t_start: theta_start,
        t_end: theta_end,
        topology: EdgeTopology::Open,
        parameter_map: ParameterMap::Identity,
        provenance: EdgeTripleProvenance::AnalyticIntersection {
            kind: AnalyticPair::SphereCone,
            branch: BranchSelector::OnlyBranch,
            cardinality: None,
        },
        quality: EdgeTripleQuality::default(),
    };

    Some(super::analytic::finalize_quality_pub(
        triple, surf_a, surf_b, tol,
    ))
}
