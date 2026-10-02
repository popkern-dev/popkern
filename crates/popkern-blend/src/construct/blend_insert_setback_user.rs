//! User-set corner setbacks at the setback rung (PK-623,
//! `docs/architecture/fillet-corner-setback.md`): which corners carry
//! one, the up-front range check, and the strict decline that replaces
//! a fallback.

use popkern_core::geometry::Point3;
use popkern_core::topology::{EdgeHandle, VertexHandle};

use crate::construct::kp20_paper_border_extract::setback_station_ranges;
use crate::fillet_surface::FilletSurfaceResult;
use crate::types::{ApexConstructionDiagnostic, CornerConstruction, CornerSetback};

/// What the fillet asks of one corner: its corner-cap preference
/// (PK-603) and its user setbacks (PK-623), if any.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct CornerRequest {
    pub cap: crate::types::CornerCapPreference,
    pub user: Option<UserSetback>,
}

impl CornerRequest {
    /// The request at `vertex` for the rounds `edge_indices`; `input`
    /// is the body before any blend (its edges are the rounds' sharp
    /// edges).
    pub(crate) fn at(
        cap: crate::types::CornerCapPreference,
        setback: &CornerSetback,
        vertex: VertexHandle,
        (edge_indices, edge_geometries): (&[usize], &[crate::construct::EdgeGeometry]),
        input: &popkern_core::topology::BRep,
    ) -> Self {
        CornerRequest {
            cap,
            user: UserSetback::at(setback, vertex, edge_indices, edge_geometries, input),
        }
    }

    /// The preference alone (no user setback).
    #[cfg(test)]
    pub(crate) fn cap_only(cap: crate::types::CornerCapPreference) -> Self {
        CornerRequest { cap, user: None }
    }
}

/// A corner's user setbacks: one station per round (by fillet index;
/// `None` keeps the automatic station) and each round's edge.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct UserSetback {
    pub stations: Vec<Option<f64>>,
    pub edges: Vec<EdgeHandle>,
    /// Each round's room along its edge on the input body: the edge's
    /// length, less the far end's user setback when that end is also a
    /// setback corner.
    pub edge_lengths: Vec<f64>,
    /// Whether an override names this corner (a uniform setback alone
    /// leaves corners the rung does not reach untouched).
    pub named: bool,
}

impl UserSetback {
    /// The user setbacks `setback` puts on `vertex`'s rounds (the edges
    /// of `edge_indices`, in fillet order), or `None` when every one is
    /// automatic.
    pub(crate) fn at(
        setback: &CornerSetback,
        vertex: VertexHandle,
        edge_indices: &[usize],
        edge_geometries: &[crate::construct::EdgeGeometry],
        input: &popkern_core::topology::BRep,
    ) -> Option<Self> {
        let edges: Vec<EdgeHandle> = edge_indices
            .iter()
            .filter_map(|&i| edge_geometries.get(i).map(|g| g.analysis.edge))
            .collect();
        let stations: Vec<Option<f64>> = edges.iter().map(|&e| setback.at(vertex, e)).collect();
        // Each round's room along its edge: the edge's length, less the
        // user setback at its far end when that end is also a corner of
        // three or more rounds (two setbacks on one edge must not overlap).
        let room = |&(i, e): &(usize, EdgeHandle)| {
            let length = edge_length(input, e).unwrap_or(f64::INFINITY);
            let Some(analysis) = edge_geometries.get(i).map(|g| &g.analysis) else {
                return length;
            };
            let far = if analysis.vertex_start == vertex {
                analysis.vertex_end
            } else {
                analysis.vertex_start
            };
            let rounds_at_far = edge_geometries
                .iter()
                .filter(|g| g.analysis.vertex_start == far || g.analysis.vertex_end == far)
                .count();
            match setback.at(far, e) {
                Some(d) if rounds_at_far >= 3 => length - d,
                _ => length,
            }
        };
        stations.iter().any(Option::is_some).then(|| UserSetback {
            stations,
            edge_lengths: edge_indices
                .iter()
                .copied()
                .zip(edges.iter().copied())
                .map(|pair| room(&pair))
                .collect(),
            edges,
            named: setback.names_vertex(vertex),
        })
    }
}

/// The first round whose user station lies outside what the corner
/// admits, as the diagnostic the apply phase turns into
/// [`crate::BlendFailure::SetbackOutOfRange`]; a range the corner
/// cannot compute declines strictly with its reason.
pub(crate) fn out_of_range(
    user: &UserSetback,
    fillet_surfaces: &[FilletSurfaceResult],
    vertex: VertexHandle,
    vertex_position: Point3,
) -> Option<ApexConstructionDiagnostic> {
    let ranges = match setback_station_ranges(fillet_surfaces, vertex_position) {
        Ok(ranges) => ranges,
        Err(reason) => {
            return Some(unavailable(
                vertex,
                format!("the user setback's range could not be measured: {reason}"),
            ));
        }
    };
    user.stations.iter().zip(&ranges).zip(&user.edges).find_map(
        |((station, &(minimum, maximum)), &edge)| {
            let requested = (*station)?;
            // The rails run past the edge's far end (trim overshoot): a
            // setback is also bounded by the sharp edge it sits on.
            let edge_length = user
                .edges
                .iter()
                .position(|&e| e == edge)
                .map(|k| user.edge_lengths[k]);
            let maximum = maximum.min(edge_length.unwrap_or(maximum));
            (!(requested > minimum && requested <= maximum)).then_some(
                ApexConstructionDiagnostic::SetbackOutOfRange {
                    vertex,
                    edge,
                    requested,
                    minimum,
                    maximum,
                },
            )
        },
    )
}

/// The strict decline of a corner carrying a user setback: no fallback
/// may serve it with another setback.
pub(crate) fn unavailable(vertex: VertexHandle, reason: String) -> ApexConstructionDiagnostic {
    ApexConstructionDiagnostic::CornerCapUnavailable {
        vertex,
        construction: CornerConstruction::Fitted,
        reason: format!("user setback: {reason}"),
    }
}

/// The arclength of `edge`'s curve (64 chords), or the chord between its
/// vertices when it has no curve.
fn edge_length(brep: &popkern_core::topology::BRep, edge: EdgeHandle) -> Option<f64> {
    use popkern_core::geometry::curve::CurveEval;
    let data = brep.edges.get(edge)?;
    if let Some(curve) = data.curve.as_ref() {
        let eval: &dyn CurveEval = &curve.geometry;
        let (a, b) = (curve.t_start, curve.t_end);
        const CHORDS: usize = 64;
        let at = |k: usize| eval.point_at(a + (b - a) * k as f64 / CHORDS as f64);
        return Some((0..CHORDS).map(|k| at(k).distance_to(&at(k + 1))).sum());
    }
    let he = brep.half_edges.get(data.half_edge)?;
    let twin = brep.half_edges.get(he.twin)?;
    Some(
        brep.vertices[he.origin]
            .position
            .distance_to(&brep.vertices[twin.origin].position),
    )
}
