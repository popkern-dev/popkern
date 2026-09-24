//! Blend-territory consistency gate (PK-530).
//!
//! The PK-516 corner-resolution gate measures only corners of two or
//! more fillets. Everything else a fillet authors served unchecked, and
//! the fast `validate()` excludes check 9 by design. Single-edge strips
//! running out onto oblique end faces returned `Ok` with their ends
//! 8.4e-2 off the strip (hexagonal-pyramid base edge), and G2 strips at
//! acute dihedrals served rails 2.8e-2 off their strips, with 5.1e-3
//! bands declared on edges against a 1e-4 Precision ceiling.
//!
//! The gate measures the blend territory: every face the operation
//! created, plus every face adjacent to one, which covers the
//! supports it trimmed and the caps it rounded. Over the edges of
//! those faces it applies:
//! - the Precision ceiling on each declared band
//!   (`BRep::max_imprecision`, the bound `validate()` holds entities
//!   to);
//! - every edge carries a curve (check 9 skips curve-less edges);
//! - check 9 at each edge's own band
//!   (`BRep::validate_curve_surface_edges`).
//!
//! Faces the operation never touched are not measured, so the gate
//! neither pays for nor fails on geometry it did not author.

use std::collections::HashSet;

use popkern_core::geometry::CurveEval;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle};

use crate::error::BlendFailure;

/// The faces the operation created (absent from `source`) and their
/// neighbours across any edge.
fn territory_faces(source: &BRep, served: &BRep) -> Vec<FaceHandle> {
    let created: HashSet<FaceHandle> = served
        .faces
        .keys()
        .filter(|f| !source.faces.contains_key(*f))
        .collect();
    let mut territory: Vec<FaceHandle> = created.iter().copied().collect();
    for he in served.half_edges.values() {
        let Some(own) = served.loops.get(he.loop_).map(|l| l.face) else {
            continue;
        };
        if !created.contains(&own) {
            continue;
        }
        let Some(across) = served
            .half_edges
            .get(he.twin)
            .and_then(|tw| served.loops.get(tw.loop_))
            .map(|l| l.face)
        else {
            continue;
        };
        if !territory.contains(&across) {
            territory.push(across);
        }
    }
    territory.sort();
    territory
}

/// Every edge of `faces` (outer and inner loops), deduplicated.
fn territory_edges(brep: &BRep, faces: &[FaceHandle]) -> Vec<EdgeHandle> {
    let wanted: HashSet<FaceHandle> = faces.iter().copied().collect();
    let mut edges: Vec<EdgeHandle> = brep
        .half_edges
        .values()
        .filter(|he| {
            brep.loops
                .get(he.loop_)
                .is_some_and(|l| wanted.contains(&l.face))
        })
        .map(|he| he.edge)
        .collect();
    edges.sort();
    edges.dedup();
    edges
}

fn midpoint(brep: &BRep, edge: EdgeHandle) -> Option<popkern_core::geometry::Point3> {
    let curve = brep.edges.get(edge)?.curve.as_ref()?;
    Some(curve.geometry.point_at(0.5 * (curve.t_start + curve.t_end)))
}

/// The served blend territory's inconsistency, if any, as the
/// transactional decline.
pub(crate) fn territory_inconsistency(source: &BRep, served: &BRep) -> Option<BlendFailure> {
    let edges = territory_edges(served, &territory_faces(source, served));
    let ceiling = served.max_imprecision;
    let over_ceiling: Vec<(EdgeHandle, f64)> = edges
        .iter()
        .filter_map(|e| {
            let band = served.edges.get(*e)?.tolerance?;
            (band > ceiling).then_some((*e, band))
        })
        .collect();
    if let Some(&(edge, band)) = over_ceiling.iter().max_by(|a, b| a.1.total_cmp(&b.1)) {
        return Some(BlendFailure::InconsistentGeometry {
            reason: format!(
                "{} edge(s) declare a band above the {ceiling:.1e} Precision ceiling, \
                 worst {edge:?} at {band:.3e}",
                over_ceiling.len()
            ),
            location: midpoint(served, edge),
        });
    }
    // An edge with no curve cannot be measured by check 9, which
    // skips it; that is a defect of its own (an open-sheet cap's
    // boundary once served two such edges).
    let curveless: Vec<EdgeHandle> = edges
        .iter()
        .copied()
        .filter(|e| served.edges.get(*e).is_some_and(|d| d.curve.is_none()))
        .collect();
    if let Some(&edge) = curveless.first() {
        let location = served
            .edges
            .get(edge)
            .and_then(|d| served.half_edges.get(d.half_edge))
            .and_then(|h| served.vertices.get(h.origin))
            .map(|v| v.position);
        return Some(BlendFailure::InconsistentGeometry {
            reason: format!(
                "{} edge(s) without a curve, first {edge:?}",
                curveless.len()
            ),
            location,
        });
    }
    let errors = served.validate_curve_surface_edges(&edges);
    let worst = errors
        .iter()
        .filter_map(|error| match error {
            popkern_core::error::ValidationError::CurveSurfaceDeviation {
                edge,
                max_deviation,
            } => Some((*edge, *max_deviation)),
            _ => None,
        })
        .max_by(|a, b| a.1.total_cmp(&b.1))?;
    Some(BlendFailure::InconsistentGeometry {
        reason: format!(
            "{} edge(s) off their faces (curve-surface consistency), worst {:?} at {:.3e}",
            errors.len(),
            worst.0,
            worst.1
        ),
        location: midpoint(served, worst.0),
    })
}
