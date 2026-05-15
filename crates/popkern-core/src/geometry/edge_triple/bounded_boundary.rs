//! Bounded-surface boundary dispatch for NURBS-owned edges.

use crate::geometry::Point3;
use crate::geometry::curve::{CurveEval, CurveGeometry};
use crate::geometry::edge_triple::error::EdgeTripleError;
use crate::geometry::edge_triple::gregory_boundary::{
    boundary_side_to_index, build_other_pcurve, reverse_nurbs,
};
use crate::geometry::edge_triple::inputs::{BoundarySide, EdgeAnchors, EdgeHint, SurfaceSide};
use crate::geometry::edge_triple::types::{
    EdgeTopology, EdgeTriple, EdgeTripleProvenance, EdgeTripleQuality, ParameterMap,
};
use crate::geometry::nurbs_curve::NurbsCurve;
use crate::geometry::nurbs_surface::NurbsSurface;
use crate::geometry::pcurve::Pcurve;
use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};
use crate::tolerance::Tolerance;

pub(super) fn try_dispatch_bounded_surface_boundary(
    surf_a: &SurfaceGeometry,
    surf_b: &SurfaceGeometry,
    anchors: &EdgeAnchors,
    hint: Option<&EdgeHint>,
    tol: &Tolerance,
) -> Option<Result<EdgeTriple, EdgeTripleError>> {
    let (owner_side, boundary) = match hint? {
        EdgeHint::SurfaceBoundary { side, boundary } => (*side, *boundary),
        _ => return None,
    };
    let (owner, other, owner_is_a) = match owner_side {
        SurfaceSide::A => (surf_a, surf_b, true),
        SurfaceSide::B => (surf_b, surf_a, false),
    };
    let SurfaceGeometry::Nurbs(surface) = owner else {
        return None;
    };
    let boundary_curve = extract_nurbs_boundary_curve(surface, boundary)?;
    let owner_pcurve = owner_boundary_pcurve(surface, boundary)?;
    let other_pcurve = build_other_pcurve(&boundary_curve, other)?;
    let (t_start, t_end) = boundary_curve.parameter_range();
    let (curve_3d, owner_pcurve, other_pcurve, endpoint_error) = match anchors {
        EdgeAnchors::Open { start, end } => orient_boundary(BoundaryOrientation {
            curve: boundary_curve,
            owner_pcurve,
            other_pcurve,
            start: *start,
            end: *end,
            parameter_range: (t_start, t_end),
            tolerance: tol.linear().max(1.0e-9),
        })?,
        _ => return None,
    };
    let (pcurve_a, pcurve_b) = if owner_is_a {
        (owner_pcurve, Pcurve::from(other_pcurve))
    } else {
        (Pcurve::from(other_pcurve), owner_pcurve)
    };
    let (curve_t_start, curve_t_end) = curve_3d.parameter_range();
    Some(Ok(EdgeTriple {
        curve_3d: CurveGeometry::Nurbs(curve_3d),
        pcurve_a,
        pcurve_b,
        t_start: curve_t_start,
        t_end: curve_t_end,
        topology: EdgeTopology::Open,
        parameter_map: ParameterMap::Identity,
        provenance: EdgeTripleProvenance::SurfaceBoundary {
            side: owner_side,
            boundary_index: boundary_side_to_index(boundary),
        },
        quality: EdgeTripleQuality {
            max_positional_residual: 0.0,
            max_endpoint_error: endpoint_error,
            samples: 0,
            tolerance_used: tol.linear(),
        },
    }))
}

struct BoundaryOrientation {
    curve: NurbsCurve,
    owner_pcurve: Pcurve,
    other_pcurve: NurbsCurve,
    start: Point3,
    end: Point3,
    parameter_range: (f64, f64),
    tolerance: f64,
}

fn orient_boundary(input: BoundaryOrientation) -> Option<(NurbsCurve, Pcurve, NurbsCurve, f64)> {
    let BoundaryOrientation {
        curve,
        owner_pcurve,
        other_pcurve,
        start,
        end,
        parameter_range: (t_start, t_end),
        tolerance,
    } = input;
    let p0 = curve.point_at(t_start);
    let p1 = curve.point_at(t_end);
    let forward = start.distance_to(&p0).max(end.distance_to(&p1));
    let reverse = start.distance_to(&p1).max(end.distance_to(&p0));
    if forward <= tolerance && forward <= reverse {
        Some((curve, owner_pcurve, other_pcurve, forward))
    } else if reverse <= tolerance {
        Some((
            reverse_nurbs(&curve),
            reverse_pcurve(owner_pcurve),
            reverse_nurbs(&other_pcurve),
            reverse,
        ))
    } else {
        None
    }
}

fn extract_nurbs_boundary_curve(
    surface: &NurbsSurface,
    boundary: BoundarySide,
) -> Option<NurbsCurve> {
    match boundary {
        BoundarySide::UMin => Some(nurbs_row(surface, 0)),
        BoundarySide::UMax => Some(nurbs_row(surface, surface.num_u.checked_sub(1)?)),
        BoundarySide::VMin => Some(nurbs_col(surface, 0)),
        BoundarySide::VMax => Some(nurbs_col(surface, surface.num_v.checked_sub(1)?)),
        BoundarySide::Indexed(_) => None,
    }
}

fn nurbs_row(surface: &NurbsSurface, i_u: usize) -> NurbsCurve {
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

fn nurbs_col(surface: &NurbsSurface, i_v: usize) -> NurbsCurve {
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

fn owner_boundary_pcurve(surface: &NurbsSurface, boundary: BoundarySide) -> Option<Pcurve> {
    let ((u_min, u_max), (v_min, v_max)) = surface.parameter_range();
    match boundary {
        BoundarySide::UMin => Some(Pcurve::linear((u_min, v_min), (u_min, v_max), v_min, v_max)),
        BoundarySide::UMax => Some(Pcurve::linear((u_max, v_min), (u_max, v_max), v_min, v_max)),
        BoundarySide::VMin => Some(Pcurve::linear((u_min, v_min), (u_max, v_min), u_min, u_max)),
        BoundarySide::VMax => Some(Pcurve::linear((u_min, v_max), (u_max, v_max), u_min, u_max)),
        BoundarySide::Indexed(_) => None,
    }
}

fn reverse_pcurve(pcurve: Pcurve) -> Pcurve {
    let start = pcurve.point_at(pcurve.t_end);
    let end = pcurve.point_at(pcurve.t_start);
    Pcurve::linear(
        (start.x, start.y),
        (end.x, end.y),
        pcurve.t_start,
        pcurve.t_end,
    )
}
