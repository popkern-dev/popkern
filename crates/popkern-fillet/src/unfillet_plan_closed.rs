//! Closed-rim unfillet planning (PK-180).
//!
//! `unfillet_face` has executed closed-rim (plane–cylinder annulus)
//! unfillets since Round 5c-4, but `compute_unfillet_plan` declined
//! every curved support, so the preview/plan API could not describe
//! the very blends the surgery restores. This module derives the
//! plan analytically from the same support classification the
//! executor uses: the restored sharp edge is the plane × cylinder
//! intersection circle, typed as
//! [`RestoredEdgeGeometry::Circle`].

use popkern_core::geometry::{Point3, Vec3};
use popkern_core::topology::{BRep, FaceHandle};

use crate::error::BlendError;
use crate::unfillet::UnfilletPlan;
use crate::unfillet_closed::seamed_band_bridge;
use crate::unfillet_closed_helpers::{
    Side, classify_support_pair, compute_rim_vertex_position, extract_cylinder_geometry,
    extract_plane_geometry, plane_cylinder_intersection_circle,
};

/// Geometry of the sharp edge an unfillet would restore, typed by
/// class so a caller (or a GUI preview) can highlight it without
/// re-deriving anything.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum RestoredEdgeGeometry {
    /// Open plane–plane blend: the restored edge is a straight
    /// segment on the support planes' intersection line.
    Line {
        /// Segment start.
        start: Point3,
        /// Segment end.
        end: Point3,
        /// Unit direction from start to end.
        direction: Vec3,
    },
    /// Closed rim blend (plane–cylinder annulus): the restored edge
    /// is the full intersection circle of the two supports.
    Circle {
        /// Circle centre (cylinder axis pierced through the plane).
        center: Point3,
        /// Unit circle axis (the cylinder axis).
        axis: Vec3,
        /// Circle radius (the cylinder radius).
        radius: f64,
    },
}

/// Compute the unfillet plan for a closed-rim (annulus or seamed
/// band) blend face. Mirrors the support classification of
/// [`crate::unfillet_closed::unfillet_closed_annulus`] so the plan
/// always describes what that surgery would restore.
pub(crate) fn compute_closed_rim_plan(
    brep: &BRep,
    blend_face: FaceHandle,
) -> Result<UnfilletPlan, BlendError> {
    // A seamed band charts as one 4-half-edge loop; un-seam VIRTUALLY
    // by reading its two contact rings. A plain annulus reads them
    // from its outer + inner loops.
    let face_data = brep
        .faces
        .get(blend_face)
        .ok_or_else(|| BlendError::InvalidParameter {
            parameter: "blend_face".to_string(),
            value: format!("{blend_face:?}"),
            reason: "face handle does not exist".to_string(),
        })?;
    let mut bridge_edge = None;
    let (ring_a, ring_b) = if let Some(bridge) = seamed_band_bridge(brep, blend_face) {
        bridge_edge = Some(brep.half_edges[bridge].edge);
        // The band loop is bridge → ring → bridge-twin-side → ring.
        let next = brep.half_edges[bridge].next;
        let after = brep.half_edges[next].next;
        let last = brep.half_edges[after].next;
        (next, last)
    } else if face_data.inner_loops.len() == 1 {
        let outer = brep
            .loops
            .get(face_data.outer_loop)
            .ok_or_else(|| BlendError::InvalidParameter {
                parameter: "blend_face.outer_loop".to_string(),
                value: format!("{:?}", face_data.outer_loop),
                reason: "outer loop missing".to_string(),
            })?
            .half_edge;
        let inner = brep
            .loops
            .get(face_data.inner_loops[0])
            .ok_or_else(|| BlendError::InvalidParameter {
                parameter: "blend_face.inner_loop".to_string(),
                value: format!("{:?}", face_data.inner_loops[0]),
                reason: "inner loop missing".to_string(),
            })?
            .half_edge;
        (outer, inner)
    } else {
        return Err(BlendError::InvalidParameter {
            parameter: "blend_face".to_string(),
            value: format!("{blend_face:?}"),
            reason: "closed-rim plan requires an annulus or seamed closed-fillet band".to_string(),
        });
    };

    let face_of = |he: popkern_core::topology::HalfEdgeHandle| -> Result<FaceHandle, BlendError> {
        let twin = brep.half_edges.get(he).map(|d| d.twin).ok_or_else(|| {
            BlendError::InvalidParameter {
                parameter: "blend_face.ring".to_string(),
                value: format!("{he:?}"),
                reason: "ring half-edge missing".to_string(),
            }
        })?;
        brep.half_edges
            .get(twin)
            .and_then(|t| brep.loops.get(t.loop_))
            .map(|l| l.face)
            .ok_or_else(|| BlendError::InvalidParameter {
                parameter: "blend_face.ring.twin".to_string(),
                value: format!("{twin:?}"),
                reason: "ring twin has no live loop".to_string(),
            })
    };
    let support_ring_a = face_of(ring_a)?;
    let support_ring_b = face_of(ring_b)?;

    let (plane_support, cyl_support, plane_ring_he, cyl_ring_he) =
        match classify_support_pair(brep, support_ring_a, support_ring_b) {
            Side::Outer => (support_ring_a, support_ring_b, ring_a, ring_b),
            Side::Inner => (support_ring_b, support_ring_a, ring_b, ring_a),
            Side::Ambiguous => {
                let kind = |f: FaceHandle| -> String {
                    brep.faces
                        .get(f)
                        .and_then(|d| d.surface.as_ref())
                        .map(|s| {
                            format!("{:?}", s.geometry)
                                .split(['(', ' ', '{'])
                                .next()
                                .unwrap_or("?")
                                .to_string()
                        })
                        .unwrap_or_else(|| "NONE".to_string())
                };
                return Err(BlendError::InvalidParameter {
                    parameter: "blend_face.supports".to_string(),
                    value: format!(
                        "{support_ring_a:?} ({}), {support_ring_b:?} ({})",
                        kind(support_ring_a),
                        kind(support_ring_b)
                    ),
                    reason: "closed-rim plan: could not classify which support is the plane \
                             and which is the cylinder"
                        .to_string(),
                });
            }
        };

    let plane_geom = extract_plane_geometry(brep, plane_support)?;
    let cyl_geom = extract_cylinder_geometry(brep, cyl_support)?;
    let (rim_center, rim_radius) = plane_cylinder_intersection_circle(&plane_geom, &cyl_geom)?;
    let axis_len = cyl_geom.axis.length().max(1e-300);
    let axis = cyl_geom.axis * (1.0 / axis_len);

    // Seam point: where the cylinder-side ring vertex lands on the
    // restored rim (the executor keeps that vertex and moves it
    // there) — a true point ON the restored circle, so the legacy
    // line-oriented fields still carry honest geometry.
    let cyl_ring_vertex = brep.half_edges[cyl_ring_he].origin;
    let seam_point = compute_rim_vertex_position(
        &brep.vertices[cyl_ring_vertex].position,
        &rim_center,
        &cyl_geom,
        rim_radius,
    );
    let radial = seam_point - rim_center;
    let mut tangent = axis.cross(&radial);
    let tangent_len = tangent.length();
    if tangent_len > 1e-300 {
        tangent = tangent * (1.0 / tangent_len);
    }

    let mut edges_to_remove = vec![
        brep.half_edges[plane_ring_he].edge,
        brep.half_edges[cyl_ring_he].edge,
    ];
    if let Some(bridge) = bridge_edge {
        edges_to_remove.push(bridge);
    }

    Ok(UnfilletPlan {
        blend_face,
        support_a: plane_support,
        support_b: cyl_support,
        restored_edge_start: seam_point,
        restored_edge_end: seam_point,
        restored_edge_direction: tangent,
        restored_edge: RestoredEdgeGeometry::Circle {
            center: rim_center,
            axis,
            radius: rim_radius,
        },
        long_he_to_a: Some(plane_ring_he),
        long_he_to_b: Some(cyl_ring_he),
        edges_to_remove,
    })
}
