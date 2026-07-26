//! G¹ hole filling — the tangent-continuous alternative to the
//! Coons fan.
//!
//! [`super::n_sided_fill::fill_n_sided_hole`] covers a hole with a fan
//! of 4-sided Coons patches radiating from the centroid. That is
//! robust and cheap, but the fan is only G⁰ across its *own* internal
//! spokes: adjacent sub-patches share a boundary curve without sharing
//! a tangent plane, so a filled hole shows creases running out from
//! the centre. The same shape and the same defect appear in the
//! fillet crate's exact-boundary cap, measured at up to 89° of normal
//! misalignment mid-seam.
//!
//! This module fills the same holes with the PnS bi-4 cap
//! ([`crate::geometry::gregory::fill_n_sided_boundary`]), which is G¹
//! across its sector seams by construction — that is precisely what
//! the paper's eq (2) enforces — and meets the surrounding surfaces
//! tangentially rather than merely touching them.
//!
//! # Which to use
//!
//! - Coons fan: any valence, no requirement on the surroundings, G⁰
//!   internally. Fine when the hole is small relative to the
//!   tolerance the consumer cares about.
//! - PnS cap (this module): valence 3–12, needs a usable tangent
//!   plane on the surviving neighbour of every boundary edge, G¹
//!   internally.
//!
//! Neither commits topology. Both return the patch set for a caller
//! that knows how it wants the faces built — see the module docs of
//! [`crate::geometry::gregory::fill_n_sided_boundary`].

use crate::geometry::curve::CurveEval;
use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};
use crate::geometry::{FillBoundaryEdge, NSidedFill, NSidedFillError, fill_n_sided_boundary};
use crate::geometry::{Point3, Vec3};
use crate::topology::{BRep, EdgeHandle};

use super::hole_fill::DetectedHole;

/// Why a G¹ fill could not be produced for a detected hole.
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum G1HoleFillError {
    /// A boundary edge had no curve geometry to sample.
    #[error("G1 hole fill: boundary edge {index} has no curve geometry")]
    MissingCurve {
        /// Position of the edge in the hole's boundary list.
        index: usize,
    },
    /// The face beside a boundary edge had no surface, so no tangent
    /// plane is available there.
    #[error(
        "G1 hole fill: boundary edge {index} has no neighbouring surface; the tangent plane the G1 join needs is unavailable"
    )]
    MissingNeighbourSurface {
        /// Position of the edge in the hole's boundary list.
        index: usize,
    },
    /// The neighbouring surface could not be inverted at a sample.
    #[error(
        "G1 hole fill: boundary edge {index} sample {sample} does not project onto its neighbouring surface"
    )]
    ProjectionFailed {
        /// Position of the edge in the hole's boundary list.
        index: usize,
        /// Which of the three samples (0 = start, 1 = mid, 2 = end).
        sample: usize,
    },
    /// The boundary edges are not in loop order.
    #[error(
        "G1 hole fill: boundary edges are not consecutive — edge {index} does not continue from its predecessor (gap {gap:.3e})"
    )]
    BoundaryOutOfOrder {
        /// Position of the offending edge.
        index: usize,
        /// Measured gap to the predecessor's end.
        gap: f64,
    },
    /// The underlying cap construction declined.
    #[error("G1 hole fill: {0}")]
    Fill(#[from] NSidedFillError),
}

/// Fill a detected hole with a G¹ PnS cap.
///
/// `tolerance` is the linear tolerance used to check that consecutive
/// boundary edges meet, and is passed through to the cap
/// construction. It is explicit, never global state.
///
/// # Errors
///
/// See [`G1HoleFillError`]. A hole whose surroundings cannot supply a
/// tangent plane declines here rather than producing a cap that only
/// looks continuous — callers wanting a result regardless should fall
/// back to [`super::n_sided_fill::fill_n_sided_hole`] and accept G⁰.
pub fn fill_n_sided_hole_g1(
    brep: &BRep,
    hole: &DetectedHole,
    tolerance: f64,
) -> Result<NSidedFill, G1HoleFillError> {
    let boundary = boundary_edges_with_frames(brep, hole, tolerance)?;
    Ok(fill_n_sided_boundary(&boundary, tolerance)?)
}

/// Build the per-side frame data the cap needs from the BRep.
fn boundary_edges_with_frames(
    brep: &BRep,
    hole: &DetectedHole,
    tolerance: f64,
) -> Result<Vec<FillBoundaryEdge>, G1HoleFillError> {
    let mut edges = Vec::with_capacity(hole.boundary_edges.len());
    for (index, &edge_handle) in hole.boundary_edges.iter().enumerate() {
        edges.push(frame_for_edge(brep, edge_handle, index)?);
    }
    orient_loop(&mut edges, tolerance)?;
    Ok(edges)
}

fn frame_for_edge(
    brep: &BRep,
    edge_handle: EdgeHandle,
    index: usize,
) -> Result<FillBoundaryEdge, G1HoleFillError> {
    let edge = brep
        .edges
        .get(edge_handle)
        .ok_or(G1HoleFillError::MissingCurve { index })?;
    let curve = edge
        .curve
        .as_ref()
        .ok_or(G1HoleFillError::MissingCurve { index })?;
    let surface = neighbour_surface(brep, edge_handle)
        .ok_or(G1HoleFillError::MissingNeighbourSurface { index })?;

    let (t0, t1) = (curve.t_start, curve.t_end);
    let eval = |t: f64| (&curve.geometry as &dyn CurveEval).point_at(t);
    let tm = t0 + 0.5 * (t1 - t0);
    let (start, mid, end) = (eval(t0), eval(tm), eval(t1));

    // Tangent by central difference along the curve, which works for
    // every curve kind without assuming a derivative API.
    let dt = 1e-6 * (t1 - t0).abs().max(1e-12);
    let tangent_at = |t: f64| -> Vec3 {
        let a = eval((t - dt).max(t0));
        let b = eval((t + dt).min(t1));
        (b - a).normalized().unwrap_or(Vec3::new(1.0, 0.0, 0.0))
    };

    let mut normals = [Vec3::new(0.0, 0.0, 1.0); 3];
    for (sample, point) in [start, mid, end].into_iter().enumerate() {
        let (u, v) = surface
            .closest_point(&point)
            .ok_or(G1HoleFillError::ProjectionFailed { index, sample })?;
        let n = surface.normal_at(u, v);
        normals[sample] = n
            .normalized()
            .ok_or(G1HoleFillError::ProjectionFailed { index, sample })?;
    }

    // Sample the boundary curve densely so the position row can be
    // fitted rather than pinned to its midpoint — a corner rim is an
    // arc joined to a line, which no quadratic passes through.
    let extra_samples: Vec<(f64, Point3)> = (1..8)
        .map(|k| {
            let f = f64::from(k) / 8.0;
            (f, eval(t0 + f * (t1 - t0)))
        })
        .collect();

    Ok(FillBoundaryEdge {
        start,
        mid,
        end,
        extra_samples,
        normal_start: normals[0],
        normal_mid: normals[1],
        normal_end: normals[2],
        tangent_start: tangent_at(t0),
        tangent_mid: tangent_at(tm),
        tangent_end: tangent_at(t1),
    })
}

/// The surface of the face that survived beside a hole edge.
fn neighbour_surface(brep: &BRep, edge_handle: EdgeHandle) -> Option<&SurfaceGeometry> {
    let edge = brep.edges.get(edge_handle)?;
    let he = brep.half_edges.get(edge.half_edge)?;
    let face = brep.loops.get(he.loop_)?.face;
    Some(&brep.faces.get(face)?.surface.as_ref()?.geometry)
}

/// Put every side in loop order, flipping any edge whose stored
/// direction runs against the loop.
///
/// Hole detection returns edges in adjacency order but says nothing
/// about each edge's own orientation, so roughly half of them can
/// point backwards. The cap needs `end[k] == start[k+1]`; without this
/// the ring reads as open and the fill declines on a boundary that is
/// in fact perfectly closed.
fn orient_loop(edges: &mut [FillBoundaryEdge], tolerance: f64) -> Result<(), G1HoleFillError> {
    let n = edges.len();
    if n < 2 {
        return Ok(());
    }
    // Seed: orient edge 0 so its end meets edge 1 at either endpoint.
    let e1_touches =
        |p: Point3, e: &FillBoundaryEdge| p.distance_to(&e.start).min(p.distance_to(&e.end));
    if e1_touches(edges[0].end, &edges[1]) > e1_touches(edges[0].start, &edges[1]) {
        flip(&mut edges[0]);
    }
    for k in 1..n {
        let prev_end = edges[k - 1].end;
        let to_start = prev_end.distance_to(&edges[k].start);
        let to_end = prev_end.distance_to(&edges[k].end);
        if to_end < to_start {
            flip(&mut edges[k]);
        }
        let gap = prev_end.distance_to(&edges[k].start);
        if gap > tolerance {
            return Err(G1HoleFillError::BoundaryOutOfOrder { index: k, gap });
        }
    }
    Ok(())
}

fn flip(edge: &mut FillBoundaryEdge) {
    std::mem::swap(&mut edge.start, &mut edge.end);
    std::mem::swap(&mut edge.normal_start, &mut edge.normal_end);
    std::mem::swap(&mut edge.tangent_start, &mut edge.tangent_end);
    // Reversing the side reverses its tangents.
    edge.tangent_start = edge.tangent_start * -1.0;
    edge.tangent_mid = edge.tangent_mid * -1.0;
    edge.tangent_end = edge.tangent_end * -1.0;
}

#[cfg(test)]
#[path = "n_sided_fill_g1_tests.rs"]
mod tests;
