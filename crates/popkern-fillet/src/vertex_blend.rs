//! §5.5 Vertex blend — first-class N-sided blend patch at a 3D
//! vertex. Thin public-API wrapper over [`crate::fillet_edges`]
//! that auto-selects every edge incident to the target vertex
//! and routes the call through the batch fillet pipeline.
//!
//! The batch pipeline's corner-resolution phase (see
//! `crate::construct::blend_insert`) produces a rolling-ball
//! sphere patch at the vertex as a side effect: every incident
//! edge becomes a torus strip, and where the torus strips meet
//! the corner gets filled with an N-sided spherical patch
//! tangent to all N incident faces.
//!
//! ## When to use this
//!
//! `vertex_blend_3d` is for the common "round off this corner"
//! use case where the caller knows which vertex they want
//! blended but doesn't want to enumerate edges manually. Use
//! [`crate::fillet_edges`] directly when you need:
//!
//! - Different radii per incident edge (variable profile)
//! - Fillet options other than [`crate::FilletOptions::default`]
//! - Cancellation / progress reporting
//! - Edges that only partially surround the vertex (e.g.
//!   "blend only two of the three edges at this corner")

use popkern_core::history::TrackingContext;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, VertexHandle};

use crate::construct::blend_topo::collect_edges_through_vertex;
use crate::error::BlendError;
use crate::types::{FilletOptions, FilletRadius, FilletSpec};
use crate::{BlendResult, fillet_edges};

/// §5.5 Vertex blend at `vertex` with a constant `radius`.
///
/// Collects every edge incident to `vertex` and fillets them as
/// a single batch via [`fillet_edges`]. The batch pipeline's
/// corner-resolution phase produces the vertex blend patch
/// automatically: for a trihedral vertex (e.g. a cube corner)
/// the output has 3 torus strips plus 1 spherical patch at the
/// corner; for higher-valence vertices the patch is an N-sided
/// spherical region tangent to all N incident faces.
///
/// The far ends of the incident edges become shorter (they're
/// consumed by the fillet's setback); vertices at those far ends
/// don't get their own blend unless the caller also targets them.
///
/// # Errors
///
/// - [`BlendError::InvalidParameter`] if `radius` is non-finite
///   or non-positive, or if `vertex` isn't present in `brep`.
/// - [`BlendError::InvalidParameter`] if `vertex` has no
///   incident edges (isolated vertex — nothing to blend).
/// - Any error surfaced by the underlying [`fillet_edges`] call
///   (propagated unchanged).
pub fn vertex_blend_3d(
    brep: &BRep,
    vertex: VertexHandle,
    radius: f64,
    options: &FilletOptions,
    tolerance: Tolerance,
    tracking: Option<&mut TrackingContext>,
) -> Result<BlendResult, BlendError> {
    if !radius.is_finite() || radius <= 0.0 {
        return Err(BlendError::InvalidParameter {
            parameter: "radius".into(),
            value: format!("{radius}"),
            reason: "vertex blend radius must be a positive finite value".into(),
        });
    }
    if brep.vertices.get(vertex).is_none() {
        return Err(BlendError::InvalidParameter {
            parameter: "vertex".into(),
            value: format!("{vertex:?}"),
            reason: "vertex handle is not present in the BRep".into(),
        });
    }

    let edges = collect_edges_through_vertex(brep, vertex).map_err(|reason| {
        BlendError::InvalidParameter {
            parameter: "vertex".into(),
            value: format!("{vertex:?}"),
            reason: format!("cannot enumerate incident edges: {reason}"),
        }
    })?;
    if edges.is_empty() {
        return Err(BlendError::InvalidParameter {
            parameter: "vertex".into(),
            value: format!("{vertex:?}"),
            reason: "vertex has no incident edges — nothing to blend".into(),
        });
    }

    let specs: Vec<FilletSpec> = edges
        .into_iter()
        .map(|edge| FilletSpec {
            edge,
            radius: FilletRadius::Constant(radius),
        })
        .collect();

    fillet_edges(brep, &specs, options, tolerance, tracking)
}

/// Variable-radius §5.5 vertex blend. Each `(edge, radius)` pair
/// in `edge_radii` specifies the blend size for one incident
/// edge of `vertex`; unspecified incident edges are rejected as
/// an error (to avoid silent under-blending that would leave
/// sharp edge fragments adjacent to the sphere patch).
///
/// This is the vertex-blend analog of §5.7 variable-profile
/// edge-sequence fillets: it lets the caller shape an asymmetric
/// corner blend — e.g. a gentle 0.15 radius on two edges and a
/// tighter 0.05 on the third.
///
/// # Errors
///
/// - [`BlendError::InvalidParameter`] if:
///   - `vertex` isn't present in `brep`
///   - any radius in `edge_radii` is non-finite or non-positive
///   - any edge in `edge_radii` isn't incident to `vertex`
///   - `edge_radii` has duplicate edges
///   - any edge incident to `vertex` is missing from `edge_radii`
/// - Any error surfaced by the underlying [`fillet_edges`] call.
pub fn vertex_blend_3d_variable(
    brep: &BRep,
    vertex: VertexHandle,
    edge_radii: &[(EdgeHandle, f64)],
    options: &FilletOptions,
    tolerance: Tolerance,
    tracking: Option<&mut TrackingContext>,
) -> Result<BlendResult, BlendError> {
    if brep.vertices.get(vertex).is_none() {
        return Err(BlendError::InvalidParameter {
            parameter: "vertex".into(),
            value: format!("{vertex:?}"),
            reason: "vertex handle is not present in the BRep".into(),
        });
    }
    for &(edge, r) in edge_radii {
        if !r.is_finite() || r <= 0.0 {
            return Err(BlendError::InvalidParameter {
                parameter: "edge_radii".into(),
                value: format!("({edge:?}, {r})"),
                reason: "vertex blend radius must be a positive finite value".into(),
            });
        }
    }

    // Dedupe check.
    let mut seen: Vec<EdgeHandle> = Vec::with_capacity(edge_radii.len());
    for &(edge, _) in edge_radii {
        if seen.contains(&edge) {
            return Err(BlendError::InvalidParameter {
                parameter: "edge_radii".into(),
                value: format!("{edge:?}"),
                reason: "duplicate edge in edge_radii list".into(),
            });
        }
        seen.push(edge);
    }

    // Every edge in edge_radii must be incident to `vertex`, and
    // every edge incident to `vertex` must appear in edge_radii.
    let incident = collect_edges_through_vertex(brep, vertex).map_err(|reason| {
        BlendError::InvalidParameter {
            parameter: "vertex".into(),
            value: format!("{vertex:?}"),
            reason: format!("cannot enumerate incident edges: {reason}"),
        }
    })?;
    for &(edge, _) in edge_radii {
        if !incident.contains(&edge) {
            return Err(BlendError::InvalidParameter {
                parameter: "edge_radii".into(),
                value: format!("{edge:?}"),
                reason: "edge is not incident to `vertex`".into(),
            });
        }
    }
    for &edge in &incident {
        if !edge_radii.iter().any(|(e, _)| *e == edge) {
            return Err(BlendError::InvalidParameter {
                parameter: "edge_radii".into(),
                value: format!("{edge:?}"),
                reason: "incident edge missing from edge_radii — all must be specified".into(),
            });
        }
    }

    let specs: Vec<FilletSpec> = edge_radii
        .iter()
        .map(|(edge, r)| FilletSpec {
            edge: *edge,
            radius: FilletRadius::Constant(*r),
        })
        .collect();

    fillet_edges(brep, &specs, options, tolerance, tracking)
}

#[cfg(test)]
#[path = "vertex_blend_tests.rs"]
mod tests;
