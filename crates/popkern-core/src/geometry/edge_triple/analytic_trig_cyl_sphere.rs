//! Cylinder × sphere analytic dispatch.
//!
//! When the sphere centre lies on the cylinder axis, the
//! intersection is up to two circles at heights
//! ±√(r_sphere² − r_cyl²) along the cyl axis (measured from sphere
//! centre projection). Each circle is in a plane perpendicular to
//! the cyl axis, with radius = r_cyl.
//!
//! This case is the workhorse for fillet-vertex-blend topology:
//! corner caps frequently produce sphere blends sitting on
//! cylinder fillet ribbons whose axes pass through the corner.
//!
//! Off-axis sphere centres produce Viviani-like or generic
//! transcendental space curves and fall through.

use crate::geometry::curve::CurveGeometry;
use crate::geometry::edge_triple::error::EdgeTripleError;
use crate::geometry::edge_triple::inputs::{EdgeAnchors, EdgeHint, SurfaceSide};
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
pub(super) fn try_cyl_sphere_axial(
    cyl_origin: Point3,
    cyl_axis: Vec3,
    cyl_radius: f64,
    sphere_is_a: bool,
    sphere_center: Point3,
    sphere_radius: f64,
    anchors: &EdgeAnchors,
    hint: Option<&EdgeHint>,
    tol: &Tolerance,
) -> Option<Result<EdgeTriple, EdgeTripleError>> {
    // Sphere centre's projection onto cyl axis.
    let center_offset = sphere_center - cyl_origin;
    let v_along_axis = center_offset.dot(&cyl_axis);
    let radial_offset = center_offset - cyl_axis * v_along_axis;
    if radial_offset.length() > tol.linear() {
        return None; // off-axis sphere
    }

    // Two circles at v = v_along_axis ± h.
    let h_sq = sphere_radius * sphere_radius - cyl_radius * cyl_radius;
    if h_sq <= tol.linear() {
        return None; // tangent or non-intersecting
    }
    let h = h_sq.sqrt();

    // Branch disambiguation: pick which of the two circles by the
    // anchor's z-projection (along cyl axis from sphere centre).
    let pick_top = match (anchors, hint) {
        (_, Some(EdgeHint::BranchIndex(i))) => *i == 0,
        (EdgeAnchors::Open { start, .. }, _) => (*start - sphere_center).dot(&cyl_axis) >= 0.0,
        _ => true,
    };
    let circle_v = if pick_top {
        v_along_axis + h
    } else {
        v_along_axis - h
    };
    let circle_center = cyl_origin + cyl_axis * circle_v;

    // Sphere v at this height.
    let sphere_h = (circle_center - sphere_center).dot(&cyl_axis);
    // For sphere with polar (0,0,1), the sphere v = colatitude depends
    // on world-z, not on cyl axis. Detect alignment.
    let polar = Vec3::new(0.0, 0.0, 1.0);
    if cyl_axis.dot(&polar).abs() < 1.0 - tol.angular() {
        return None; // cyl axis ≠ sphere polar — sphere pcurve is transcendental
    }

    let (arc_basis_u, arc_basis_v) = crate::geometry::analytic_curves::circle_basis(cyl_axis);
    let (cyl_basis_u, cyl_basis_v) =
        crate::geometry::analytic_surfaces::orthonormal_basis(cyl_axis);
    let cyl_u_offset = arc_basis_u
        .dot(&cyl_basis_v)
        .patan2(arc_basis_u.dot(&cyl_basis_u));
    // Sphere u-axis = global X for points on a constant-z circle
    // (since sphere polar is (0,0,1) and we need cyl_axis aligned).
    let global_x = Vec3::new(1.0, 0.0, 0.0);
    let global_y = Vec3::new(0.0, 1.0, 0.0);
    let sphere_u_offset = arc_basis_u
        .dot(&global_y)
        .patan2(arc_basis_u.dot(&global_x));

    let (start, end) = match anchors {
        EdgeAnchors::Open { start, end } => (*start, *end),
        _ => return None,
    };
    let theta_of = |p: Point3| -> Option<f64> {
        let d = p - circle_center;
        let cu = d.dot(&arc_basis_u);
        let cv = d.dot(&arc_basis_v);
        let r = (cu * cu + cv * cv).sqrt();
        if (r - cyl_radius).abs() > tol.linear() {
            return None;
        }
        Some(cv.patan2(cu))
    };
    let theta_start = match theta_of(start) {
        Some(t) => t,
        None => {
            return Some(Err(EdgeTripleError::AnchorOffSurface {
                which: SurfaceSide::A,
                distance: (start - circle_center).length() - cyl_radius,
            }));
        }
    };
    let theta_end = match theta_of(end) {
        Some(t) => t,
        None => {
            return Some(Err(EdgeTripleError::AnchorOffSurface {
                which: SurfaceSide::B,
                distance: (end - circle_center).length() - cyl_radius,
            }));
        }
    };
    let theta_end = if theta_end < theta_start {
        theta_end + std::f64::consts::TAU
    } else {
        theta_end
    };

    // Cyl pcurve: u = θ + cyl offset, v = circle_v.
    let pcurve_cyl = Pcurve::linear(
        (theta_start + cyl_u_offset, circle_v),
        (theta_end + cyl_u_offset, circle_v),
        theta_start,
        theta_end,
    );

    // Sphere pcurve: u = θ + sphere offset, v = colatitude(sphere_h).
    let v_sphere = (sphere_h / sphere_radius).clamp(-1.0, 1.0).pacos();
    let pcurve_sphere = Pcurve::linear(
        (theta_start + sphere_u_offset, v_sphere),
        (theta_end + sphere_u_offset, v_sphere),
        theta_start,
        theta_end,
    );

    let curve_3d = CurveGeometry::Arc {
        center: circle_center,
        normal: cyl_axis,
        radius: cyl_radius,
        start_angle: theta_start,
        end_angle: theta_end,
    };

    let cyl_surf = SurfaceGeometry::Cylinder {
        origin: cyl_origin,
        axis: cyl_axis,
        radius: cyl_radius,
    };
    let sphere_surf = SurfaceGeometry::Sphere {
        center: sphere_center,
        radius: sphere_radius,
    };

    let (pcurve_a, pcurve_b, surf_a, surf_b) = if sphere_is_a {
        (pcurve_sphere, pcurve_cyl, &sphere_surf, &cyl_surf)
    } else {
        (pcurve_cyl, pcurve_sphere, &cyl_surf, &sphere_surf)
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
            kind: AnalyticPair::CylinderSphere,
            branch: if pick_top {
                BranchSelector::Indexed(0)
            } else {
                BranchSelector::Indexed(1)
            },
            cardinality: None,
        },
        quality: EdgeTripleQuality::default(),
    };

    Some(super::analytic::finalize_quality_pub(
        triple, surf_a, surf_b, tol,
    ))
}
