//! Boundary-evidence helpers for edge/surface compatibility repair.

use crate::geometry::curve::CurveEval;
use crate::geometry::edge_triple::{
    BoundarySide, EdgeHint, SurfaceSide, extract_gregory_boundary_curve,
};
use crate::geometry::nurbs_curve::NurbsCurve;
use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};
use crate::topology::{
    BRep, HalfEdgeHandle, SurfaceBoundaryHintDecline, SurfaceBoundaryHintDeclineReason,
};

/// Infer a canonical surface-boundary hint from an owning half-edge pcurve.
pub(super) fn infer_surface_boundary_hint(
    brep: &mut BRep,
    half_edge: HalfEdgeHandle,
    surface: &SurfaceGeometry,
    side: SurfaceSide,
    tolerance: f64,
) -> Option<EdgeHint> {
    if !supports_boundary_owned_reconstruction(surface) {
        return None;
    }
    if let Some(boundary) = brep.surface_boundary_hint(half_edge) {
        match validate_authored_boundary_hint(brep, half_edge, surface, boundary, tolerance) {
            AuthoredBoundaryHintValidation::Accepted => {
                return Some(EdgeHint::SurfaceBoundary { side, boundary });
            }
            AuthoredBoundaryHintValidation::Declined(decline) => {
                brep.set_surface_boundary_decline(half_edge, decline);
                return infer_boundary_from_topology(brep, half_edge, surface, side, tolerance);
            }
        }
    }
    infer_boundary_from_pcurve(brep, half_edge, surface, side, tolerance)
        .or_else(|| infer_boundary_from_topology(brep, half_edge, surface, side, tolerance))
}

/// Infer a canonical surface-boundary side without mutating boundary diagnostics.
pub(super) fn infer_surface_boundary_side_read_only(
    brep: &BRep,
    half_edge: HalfEdgeHandle,
    surface: &SurfaceGeometry,
    tolerance: f64,
) -> Option<BoundarySide> {
    if !supports_boundary_owned_reconstruction(surface) {
        return None;
    }
    if let Some(boundary) = brep.surface_boundary_hint(half_edge)
        && matches!(
            validate_authored_boundary_hint(brep, half_edge, surface, boundary, tolerance),
            AuthoredBoundaryHintValidation::Accepted
        )
    {
        return Some(boundary);
    }
    hint_boundary(infer_boundary_from_pcurve(
        brep,
        half_edge,
        surface,
        SurfaceSide::A,
        tolerance,
    ))
    .or_else(|| {
        hint_boundary(infer_boundary_from_topology(
            brep,
            half_edge,
            surface,
            SurfaceSide::A,
            tolerance,
        ))
    })
}

fn hint_boundary(hint: Option<EdgeHint>) -> Option<BoundarySide> {
    match hint {
        Some(EdgeHint::SurfaceBoundary { boundary, .. }) => Some(boundary),
        _ => None,
    }
}

enum AuthoredBoundaryHintValidation {
    Accepted,
    Declined(SurfaceBoundaryHintDecline),
}

fn validate_authored_boundary_hint(
    brep: &BRep,
    half_edge: HalfEdgeHandle,
    surface: &SurfaceGeometry,
    boundary: BoundarySide,
    tolerance: f64,
) -> AuthoredBoundaryHintValidation {
    let Some(pcurve) = brep
        .half_edges
        .get(half_edge)
        .and_then(|he| he.pcurve.as_ref())
    else {
        return AuthoredBoundaryHintValidation::Declined(SurfaceBoundaryHintDecline {
            reason: SurfaceBoundaryHintDeclineReason::MissingPcurve,
            best_boundary_deviation: f64::INFINITY,
            best_boundary: None,
            endpoint_gap: f64::INFINITY,
        });
    };
    let classification = classify_boundary_pcurve(surface, pcurve);
    let authored_deviation = classification.deviation_for(boundary);
    let threshold = tolerance.max(1.0e-7);
    if authored_deviation > threshold {
        return AuthoredBoundaryHintValidation::Declined(SurfaceBoundaryHintDecline {
            reason: SurfaceBoundaryHintDeclineReason::AuthoredHintBoundaryMismatch,
            best_boundary_deviation: classification.best_deviation,
            best_boundary: classification.best_boundary,
            endpoint_gap: 0.0,
        });
    }
    let endpoint_gap =
        authored_boundary_endpoint_gap(brep, half_edge, surface, boundary).unwrap_or(f64::INFINITY);
    if endpoint_gap > threshold {
        return AuthoredBoundaryHintValidation::Declined(SurfaceBoundaryHintDecline {
            reason: SurfaceBoundaryHintDeclineReason::AuthoredHintEndpointMismatch,
            best_boundary_deviation: classification.best_deviation,
            best_boundary: classification.best_boundary,
            endpoint_gap,
        });
    }
    AuthoredBoundaryHintValidation::Accepted
}

fn infer_boundary_from_pcurve(
    brep: &BRep,
    half_edge: HalfEdgeHandle,
    surface: &SurfaceGeometry,
    side: SurfaceSide,
    tolerance: f64,
) -> Option<EdgeHint> {
    let pcurve = brep.half_edges.get(half_edge)?.pcurve.as_ref()?;
    let classification = classify_boundary_pcurve(surface, pcurve);
    let threshold = tolerance.max(1.0e-8);
    let boundary =
        (classification.best_deviation <= threshold).then_some(classification.best_boundary?)?;
    Some(EdgeHint::SurfaceBoundary { side, boundary })
}

fn infer_boundary_from_topology(
    brep: &BRep,
    half_edge: HalfEdgeHandle,
    surface: &SurfaceGeometry,
    side: SurfaceSide,
    tolerance: f64,
) -> Option<EdgeHint> {
    let he = brep.half_edges.get(half_edge)?;
    let next = brep.half_edges.get(he.next)?;
    let start = brep.vertices.get(he.origin)?.position;
    let end = brep.vertices.get(next.origin)?.position;
    let threshold = tolerance.max(1.0e-7);
    [
        BoundarySide::UMin,
        BoundarySide::UMax,
        BoundarySide::VMin,
        BoundarySide::VMax,
    ]
    .into_iter()
    .filter_map(|boundary| {
        let curve = surface_boundary_curve(surface, boundary)?;
        let (t0, t1) = curve.parameter_range();
        let start_t = closest_curve_parameter(&curve, start, t0, t1)?;
        let end_t = closest_curve_parameter(&curve, end, t0, t1)?;
        if (end_t - start_t).abs() <= 1.0e-12 {
            return None;
        }
        let gap = curve
            .point_at(start_t)
            .distance_to(&start)
            .max(curve.point_at(end_t).distance_to(&end));
        Some((boundary, gap))
    })
    .filter(|(_, gap)| *gap <= threshold)
    .min_by(|a, b| a.1.total_cmp(&b.1))
    .map(|(boundary, _)| EdgeHint::SurfaceBoundary { side, boundary })
}

struct BoundaryPcurveClassification {
    best_boundary: Option<BoundarySide>,
    best_deviation: f64,
    u_min: f64,
    u_max: f64,
    v_min: f64,
    v_max: f64,
}

impl BoundaryPcurveClassification {
    fn deviation_for(&self, boundary: BoundarySide) -> f64 {
        match boundary {
            BoundarySide::UMin => self.u_min,
            BoundarySide::UMax => self.u_max,
            BoundarySide::VMin => self.v_min,
            BoundarySide::VMax => self.v_max,
            BoundarySide::Indexed(_) => f64::INFINITY,
        }
    }
}

fn classify_boundary_pcurve(
    surface: &SurfaceGeometry,
    pcurve: &crate::geometry::pcurve::Pcurve,
) -> BoundaryPcurveClassification {
    let ((u_min, u_max), (v_min, v_max)) = surface.parameter_range();
    if !u_min.is_finite() || !u_max.is_finite() || !v_min.is_finite() || !v_max.is_finite() {
        return BoundaryPcurveClassification {
            best_boundary: None,
            best_deviation: f64::INFINITY,
            u_min: f64::INFINITY,
            u_max: f64::INFINITY,
            v_min: f64::INFINITY,
            v_max: f64::INFINITY,
        };
    }
    let (t0, t1) = pcurve.parameter_range();
    let mut max_u_min = 0.0_f64;
    let mut max_u_max = 0.0_f64;
    let mut max_v_min = 0.0_f64;
    let mut max_v_max = 0.0_f64;
    for index in 0..9 {
        let f = index as f64 / 8.0;
        let p = pcurve.point_at(t0 + (t1 - t0) * f);
        max_u_min = max_u_min.max((p.x - u_min).abs());
        max_u_max = max_u_max.max((p.x - u_max).abs());
        max_v_min = max_v_min.max((p.y - v_min).abs());
        max_v_max = max_v_max.max((p.y - v_max).abs());
    }
    let (best_boundary, best_deviation) = [
        (BoundarySide::UMin, max_u_min),
        (BoundarySide::UMax, max_u_max),
        (BoundarySide::VMin, max_v_min),
        (BoundarySide::VMax, max_v_max),
    ]
    .into_iter()
    .min_by(|a, b| a.1.total_cmp(&b.1))
    .map(|(boundary, deviation)| (Some(boundary), deviation))
    .unwrap_or((None, f64::INFINITY));
    BoundaryPcurveClassification {
        best_boundary,
        best_deviation,
        u_min: max_u_min,
        u_max: max_u_max,
        v_min: max_v_min,
        v_max: max_v_max,
    }
}

fn authored_boundary_endpoint_gap(
    brep: &BRep,
    half_edge: HalfEdgeHandle,
    surface: &SurfaceGeometry,
    boundary: BoundarySide,
) -> Option<f64> {
    let half_edge_data = brep.half_edges.get(half_edge)?;
    let next = brep.half_edges.get(half_edge_data.next)?;
    let start = brep.vertices.get(half_edge_data.origin)?.position;
    let end = brep.vertices.get(next.origin)?.position;
    let curve = surface_boundary_curve(surface, boundary)?;
    let (t0, t1) = curve.parameter_range();
    let start_t = closest_curve_parameter(&curve, start, t0, t1)?;
    let end_t = closest_curve_parameter(&curve, end, t0, t1)?;
    Some(
        curve
            .point_at(start_t)
            .distance_to(&start)
            .max(curve.point_at(end_t).distance_to(&end)),
    )
}

fn supports_boundary_owned_reconstruction(surface: &SurfaceGeometry) -> bool {
    matches!(
        surface,
        SurfaceGeometry::Gregory(_) | SurfaceGeometry::Nurbs(_)
    )
}

pub(super) fn surface_boundary_curve(
    surface: &SurfaceGeometry,
    boundary: BoundarySide,
) -> Option<NurbsCurve> {
    match surface {
        SurfaceGeometry::Gregory(patch) => extract_gregory_boundary_curve(patch, boundary),
        SurfaceGeometry::Nurbs(nurbs) => match boundary {
            BoundarySide::UMin => Some(nurbs_row(nurbs, 0)),
            BoundarySide::UMax => Some(nurbs_row(nurbs, nurbs.num_u.checked_sub(1)?)),
            BoundarySide::VMin => Some(nurbs_col(nurbs, 0)),
            BoundarySide::VMax => Some(nurbs_col(nurbs, nurbs.num_v.checked_sub(1)?)),
            BoundarySide::Indexed(_) => None,
        },
        _ => None,
    }
}

fn nurbs_row(surface: &crate::geometry::nurbs_surface::NurbsSurface, i_u: usize) -> NurbsCurve {
    let control_points = (0..surface.num_v)
        .map(|i_v| surface.control_points[i_u * surface.num_v + i_v])
        .collect();
    let weights = (0..surface.num_v)
        .map(|i_v| surface.weights[i_u * surface.num_v + i_v])
        .collect();
    NurbsCurve {
        degree: surface.degree_v,
        knots: surface.knots_v.clone(),
        control_points,
        weights,
    }
}

fn nurbs_col(surface: &crate::geometry::nurbs_surface::NurbsSurface, i_v: usize) -> NurbsCurve {
    let control_points = (0..surface.num_u)
        .map(|i_u| surface.control_points[i_u * surface.num_v + i_v])
        .collect();
    let weights = (0..surface.num_u)
        .map(|i_u| surface.weights[i_u * surface.num_v + i_v])
        .collect();
    NurbsCurve {
        degree: surface.degree_u,
        knots: surface.knots_u.clone(),
        control_points,
        weights,
    }
}

fn closest_curve_parameter(
    curve: &NurbsCurve,
    target: crate::geometry::Point3,
    t_start: f64,
    t_end: f64,
) -> Option<f64> {
    let samples = 33usize;
    let mut best_index = 0usize;
    let mut best_distance = f64::INFINITY;
    for index in 0..samples {
        let f = index as f64 / (samples - 1) as f64;
        let t = t_start + (t_end - t_start) * f;
        let distance = curve.point_at(t).distance_to(&target);
        if distance < best_distance {
            best_distance = distance;
            best_index = index;
        }
    }
    if !best_distance.is_finite() {
        return None;
    }
    let step = (t_end - t_start) / (samples - 1) as f64;
    let mut lo = t_start + step * best_index.saturating_sub(1) as f64;
    let mut hi = t_start + step * (best_index + 1).min(samples - 1) as f64;
    for _ in 0..32 {
        let m1 = lo + (hi - lo) / 3.0;
        let m2 = hi - (hi - lo) / 3.0;
        if curve.point_at(m1).distance_to(&target) <= curve.point_at(m2).distance_to(&target) {
            hi = m2;
        } else {
            lo = m1;
        }
    }
    Some(0.5 * (lo + hi))
}

/// True when the generic analytic intersection factory has bounded repair coverage.
pub(super) fn is_repairable_surface_pair(a: &SurfaceGeometry, b: &SurfaceGeometry) -> bool {
    use SurfaceGeometry::{Cone, Cylinder, Gregory, Plane, Sphere, Torus};
    matches!(
        (a, b),
        (Plane { .. }, Plane { .. })
            | (Plane { .. }, Cylinder { .. })
            | (Cylinder { .. }, Plane { .. })
            | (Plane { .. }, Cone { .. })
            | (Cone { .. }, Plane { .. })
            | (Plane { .. }, Sphere { .. })
            | (Sphere { .. }, Plane { .. })
            | (Plane { .. }, Torus { .. })
            | (Torus { .. }, Plane { .. })
            | (Cylinder { .. }, Cylinder { .. })
            | (Cylinder { .. }, Cone { .. })
            | (Cone { .. }, Cylinder { .. })
            | (Cylinder { .. }, Sphere { .. })
            | (Sphere { .. }, Cylinder { .. })
            | (Cylinder { .. }, Torus { .. })
            | (Torus { .. }, Cylinder { .. })
            | (Cone { .. }, Cone { .. })
            | (Cone { .. }, Sphere { .. })
            | (Sphere { .. }, Cone { .. })
            | (Cone { .. }, Torus { .. })
            | (Torus { .. }, Cone { .. })
            | (Sphere { .. }, Sphere { .. })
            | (Sphere { .. }, Torus { .. })
            | (Torus { .. }, Sphere { .. })
            | (Torus { .. }, Torus { .. })
            | (Gregory(_), Cylinder { .. })
            | (Cylinder { .. }, Gregory(_))
            | (Gregory(_), Plane { .. })
            | (Plane { .. }, Gregory(_))
            | (Gregory(_), Sphere { .. })
            | (Sphere { .. }, Gregory(_))
            | (Gregory(_), Torus { .. })
            | (Torus { .. }, Gregory(_))
    )
}
