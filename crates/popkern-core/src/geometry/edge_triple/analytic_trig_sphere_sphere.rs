//! Sphere × sphere analytic dispatch.
//!
//! Two spheres intersect in a circle whenever the distance between
//! centres is between |r₁ − r₂| and r₁ + r₂. The circle lies in
//! the plane equidistant-of-distances from each centre, with normal
//! along the line connecting them.
//!
//! For both pcurves to be `Pcurve::Linear` (constant v, linear u),
//! the circle normal must equal each sphere's UV polar axis. Since
//! popkern's sphere uses a fixed polar axis of `(0, 0, 1)`, this
//! requires the line c₂ − c₁ to be parallel to `(0, 0, 1)`. Other
//! configurations return `None` and fall through to the legacy
//! NURBS path.

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

#[allow(clippy::too_many_arguments)]
pub(super) fn try_sphere_sphere_polar(
    center_a: Point3,
    radius_a: f64,
    center_b: Point3,
    radius_b: f64,
    anchors: &EdgeAnchors,
    tol: &Tolerance,
) -> Option<Result<EdgeTriple, EdgeTripleError>> {
    let polar = Vec3::new(0.0, 0.0, 1.0);
    let delta = center_b - center_a;
    let d = delta.length();
    if d < tol.linear() {
        return None; // coincident or near-coincident
    }
    // Require centre-to-centre direction to align with the polar
    // axis (within angular tolerance).
    let dir = delta * (1.0 / d);
    if dir.dot(&polar).abs() < 1.0 - tol.angular() {
        return None;
    }

    if d >= radius_a + radius_b - tol.linear() || d <= (radius_a - radius_b).abs() + tol.linear() {
        return None;
    }

    // Circle parameters.
    let a = (d * d + radius_a * radius_a - radius_b * radius_b) / (2.0 * d);
    let circle_radius = (radius_a * radius_a - a * a).max(0.0).sqrt();
    if circle_radius < tol.linear() {
        return None;
    }
    let circle_center = center_a + dir * a;

    let (arc_basis_u, arc_basis_v) = crate::geometry::analytic_curves::circle_basis(polar);
    let global_x = Vec3::new(1.0, 0.0, 0.0);
    let global_y = Vec3::new(0.0, 1.0, 0.0);
    let sphere_u_offset = arc_basis_u.dot(&global_y).atan2(arc_basis_u.dot(&global_x));

    let (start, end) = match anchors {
        EdgeAnchors::Open { start, end } => (*start, *end),
        _ => return None,
    };

    let theta_of = |p: Point3| -> Option<f64> {
        let off = p - circle_center;
        let cu = off.dot(&arc_basis_u);
        let cv = off.dot(&arc_basis_v);
        let r = (cu * cu + cv * cv).sqrt();
        if (r - circle_radius).abs() > tol.linear() {
            return None;
        }
        Some(cv.atan2(cu))
    };
    let theta_start = match theta_of(start) {
        Some(t) => t,
        None => {
            return Some(Err(EdgeTripleError::AnchorOffSurface {
                which: SurfaceSide::A,
                distance: (start - circle_center).length() - circle_radius,
            }));
        }
    };
    let theta_end = match theta_of(end) {
        Some(t) => t,
        None => {
            return Some(Err(EdgeTripleError::AnchorOffSurface {
                which: SurfaceSide::B,
                distance: (end - circle_center).length() - circle_radius,
            }));
        }
    };
    let theta_end = if theta_end < theta_start {
        theta_end + std::f64::consts::TAU
    } else {
        theta_end
    };

    // Sphere A pcurve.
    let h_a = (circle_center - center_a).dot(&polar);
    let v_sphere_a = (h_a / radius_a).clamp(-1.0, 1.0).acos();
    let pcurve_a = Pcurve::linear(
        (theta_start + sphere_u_offset, v_sphere_a),
        (theta_end + sphere_u_offset, v_sphere_a),
        theta_start,
        theta_end,
    );
    // Sphere B pcurve.
    let h_b = (circle_center - center_b).dot(&polar);
    let v_sphere_b = (h_b / radius_b).clamp(-1.0, 1.0).acos();
    let pcurve_b = Pcurve::linear(
        (theta_start + sphere_u_offset, v_sphere_b),
        (theta_end + sphere_u_offset, v_sphere_b),
        theta_start,
        theta_end,
    );

    let curve_3d = CurveGeometry::Arc {
        center: circle_center,
        normal: polar,
        radius: circle_radius,
        start_angle: theta_start,
        end_angle: theta_end,
    };

    let sphere_a = SurfaceGeometry::Sphere {
        center: center_a,
        radius: radius_a,
    };
    let sphere_b = SurfaceGeometry::Sphere {
        center: center_b,
        radius: radius_b,
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
            kind: AnalyticPair::SphereSphere,
            branch: BranchSelector::OnlyBranch,
            cardinality: None,
        },
        quality: EdgeTripleQuality::default(),
    };

    Some(super::analytic::finalize_quality_pub(
        triple, &sphere_a, &sphere_b, tol,
    ))
}
