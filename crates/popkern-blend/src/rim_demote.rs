//! Temporary containment (PK-417 → PK-419): the G2 construction ladder
//! does not yet serve a closed rim the boolean re-spelled as an exact
//! `Arc`. Measured on the flagship boss plate: the union's junction
//! rim as a NURBS (the pre-PK-417 form) takes the nine-edge G2 fillet
//! watertight; the same rim as an `Arc` — whichever way round, with
//! exact, fitted or the original pcurves — leaves 70 open mesh edges
//! along the strip's rails, while the extrude's own `Arc` rim on the
//! same cylinder takes a G2 fillet cleanly. The mechanism lives in
//! the G2 strip's station origin and is filed as PK-419; until it is
//! fixed, a G2 request runs on a working copy where every closed
//! boolean-made `Arc` rim is spelled as the exact rational circle it
//! also is (`to_nurbs_in_range` over the same parameter range — the
//! same curve, the pre-PK-417 representation). G1 requests and every
//! other consumer keep the analytic arc.

use popkern_core::geometry::curve::{CurveEval, CurveGeometry};
use popkern_core::naming::RoleCategory;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_core::topology::edge_pcurve_reauthor::coparameterised_pcurves;

/// The body a fillet request runs on: under a G2 request, the copy
/// with boolean-made closed `Arc` rims spelled as NURBS (PK-419
/// containment); otherwise the input unchanged.
pub(crate) fn for_request(brep: BRep, wants_g2: bool, tolerance: &Tolerance) -> BRep {
    if !wants_g2 {
        return brep;
    }
    demote_boolean_arc_rims(&brep, tolerance).unwrap_or(brep)
}

/// A copy of `brep` with boolean-made closed `Arc` rims spelled as
/// NURBS, or `None` when there is nothing to demote.
pub(crate) fn demote_boolean_arc_rims(brep: &BRep, tolerance: &Tolerance) -> Option<BRep> {
    let mut demoted = 0usize;
    let mut work = brep.clone();
    for (handle, edge) in brep.edges.iter() {
        let Some(curve) = edge.curve.as_ref() else {
            continue;
        };
        let CurveGeometry::Arc {
            start_angle,
            end_angle,
            ..
        } = &curve.geometry
        else {
            continue;
        };
        let closed = ((end_angle - start_angle).abs() - std::f64::consts::TAU).abs() < 1.0e-9;
        let boolean_made = edge
            .name
            .and_then(|id| brep.name_pool.get(id))
            .is_some_and(|name| name.role.category == RoleCategory::Intersection);
        if !(closed && boolean_made) {
            continue;
        }
        let mut nurbs = curve
            .geometry
            .to_nurbs_in_range((curve.t_start, curve.t_end));
        // The pre-PK-417 rim was parameterised over [0, 1] from its
        // vertex; the G2 strip's station origin follows the rim's
        // parameter origin, so the spelling is re-knotted (affinely,
        // the same curve) onto that range.
        let (k0, k1) = (nurbs.knots[0], nurbs.knots[nurbs.knots.len() - 1]);
        if (k1 - k0).abs() > 0.0 {
            for k in &mut nurbs.knots {
                *k = (*k - k0) / (k1 - k0);
            }
        }
        let (t_start, t_end) = nurbs.parameter_range();
        let spelled = popkern_core::geometry::CurveOnEdge {
            geometry: CurveGeometry::Nurbs(nurbs),
            t_start,
            t_end,
        };
        let Some(pcurves) =
            coparameterised_pcurves(brep, handle, &spelled, tolerance.linear() * 10.0, tolerance)
        else {
            continue;
        };
        if let Some(e) = work.edges.get_mut(handle) {
            e.curve = Some(spelled);
        }
        for (he, pcurve) in pcurves {
            if let Some(h) = work.half_edges.get_mut(he) {
                h.pcurve = Some(pcurve);
            }
        }
        demoted += 1;
    }
    if demoted == 0 {
        return None;
    }
    if popkern_core::trace::trace_enabled("fillet") {
        eprintln!(
            "  [fillet] G2 request: {demoted} boolean arc rim(s) run as NURBS (PK-419 containment)"
        );
    }
    Some(work)
}
