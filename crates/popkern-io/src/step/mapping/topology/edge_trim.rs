//! An imported edge's trim: the parameter range its STEP vertices bound
//! on the resolved curve.

use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::{CurveEval, CurveGeometry, CurveOnEdge};
use popkern_core::tolerance::Tolerance;

use crate::step::diagnostics::{DiagnosticLog, Phase};
use crate::step::parser::types::EntityId;

pub(super) fn recover_edge_curve_trim(
    edge_curve_id: EntityId,
    mut curve: CurveOnEdge,
    start: Point3,
    end: Point3,
    same_sense: bool,
    tolerance: Tolerance,
    log: &mut DiagnosticLog,
) -> CurveOnEdge {
    if bounded_nurbs_endpoints_match(&curve, start, end, same_sense, tolerance) {
        return curve;
    }

    let Some(mut start_t) = curve.geometry.closest_point(&start) else {
        log.warn(
            Phase::MapTopology,
            Some(edge_curve_id),
            "could not recover EDGE_CURVE start parameter from topological vertex",
        );
        return curve;
    };
    let Some(mut end_t) = curve.geometry.closest_point(&end) else {
        log.warn(
            Phase::MapTopology,
            Some(edge_curve_id),
            "could not recover EDGE_CURVE end parameter from topological vertex",
        );
        return curve;
    };

    let start_gap = curve.geometry.point_at(start_t).distance_to(&start);
    let end_gap = curve.geometry.point_at(end_t).distance_to(&end);
    let limit = tolerance.linear().max(1e-9);
    if start_gap > limit || end_gap > limit {
        log.warn(
            Phase::MapTopology,
            Some(edge_curve_id),
            format!(
                "declined EDGE_CURVE trim recovery: endpoint gaps {:.3e}, {:.3e} exceed {:.3e}",
                start_gap, end_gap, limit
            ),
        );
        return curve;
    }

    normalize_periodic_parameter_order(&curve.geometry, same_sense, &mut start_t, &mut end_t);

    curve.t_start = start_t;
    curve.t_end = end_t;

    curve
}

fn bounded_nurbs_endpoints_match(
    curve: &CurveOnEdge,
    start: Point3,
    end: Point3,
    same_sense: bool,
    tolerance: Tolerance,
) -> bool {
    if !matches!(curve.geometry, CurveGeometry::Nurbs(_)) {
        return false;
    }
    let (curve_start, curve_end) = if same_sense {
        (
            curve.geometry.point_at(curve.t_start),
            curve.geometry.point_at(curve.t_end),
        )
    } else {
        (
            curve.geometry.point_at(curve.t_end),
            curve.geometry.point_at(curve.t_start),
        )
    };
    let limit = tolerance.linear().max(1e-9);
    curve_start.distance_to(&start) <= limit && curve_end.distance_to(&end) <= limit
}

pub(super) fn normalize_periodic_parameter_order(
    geometry: &CurveGeometry,
    same_sense: bool,
    start_t: &mut f64,
    end_t: &mut f64,
) {
    let is_periodic = matches!(
        geometry,
        CurveGeometry::Circle { .. } | CurveGeometry::Arc { .. } | CurveGeometry::Ellipse { .. }
    );
    if !is_periodic {
        return;
    }
    // Projection lands on any period branch (Newton from spread
    // starts converges turns away: measured turbine.stp rims at
    // t = −15.6 → 3.2, a 3-turn trim whose pcurves then spanned 4π and
    // livelocked the mesher). The trim is the vertex angles folded to
    // one turn, then ordered.
    *start_t = start_t.rem_euclid(std::f64::consts::TAU);
    *end_t = end_t.rem_euclid(std::f64::consts::TAU);

    if same_sense {
        while *end_t < *start_t {
            *end_t += std::f64::consts::TAU;
        }
        // CLOSED rim (start vertex == end vertex): both parameters
        // project to the same point and the trim would collapse to
        // zero length — the edge covers the FULL period.
        if (*end_t - *start_t).abs() < 1e-9 {
            *end_t += std::f64::consts::TAU;
        }
    } else {
        while *start_t < *end_t {
            *start_t += std::f64::consts::TAU;
        }
        if (*start_t - *end_t).abs() < 1e-9 {
            *start_t += std::f64::consts::TAU;
        }
    }
}
