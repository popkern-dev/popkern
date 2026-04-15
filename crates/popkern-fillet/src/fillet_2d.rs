//! §5.7 2D corner fillets for sheet-body planar wires.
//!
//! A **2D corner fillet** replaces a sharp convex corner of a
//! planar polygon with a circular arc of the requested radius.
//! Given a sheet body (produced by e.g. `BRep::sheet_quad`) and a
//! vertex where two straight-line edges meet, [`fillet_corner_2d`]
//! builds two tangent points (one on each adjacent edge) and a
//! new arc edge that spans between them, then deletes the original
//! corner vertex.
//!
//! This is the "sheet body" half of §5.7 Planar Wire Fillets:
//! A commercial kernel's 2D profile fillet does the analogous
//! operation for 2D wires. Since popkern represents sheet bodies
//! as two-face paper-thin shells (front + backside of the same
//! plane, opposite orientations), the surgery runs **twice** —
//! once on the front face's outer loop and once on the backside's
//! — keeping the manifold pairing intact.
//!
//! ## Algorithm
//!
//! 1. Find the two half-edges whose origin is the corner vertex.
//!    They live on distinct edges `edge_a` and `edge_b`.
//! 2. Verify both edges carry `CurveGeometry::Line`. 2D corner
//!    fillets only operate on straight-line edges today; curved
//!    boundary edges are a future extension.
//! 3. Compute unit directions `dir_a`, `dir_b` from the corner
//!    toward each edge's other endpoint.
//! 4. Compute the interior half-angle `θ/2` from
//!    `cos(θ) = dir_a · dir_b`. Reject if `θ ≈ 0` (colinear) or
//!    `θ ≈ π` (near-straight, no meaningful corner).
//! 5. Tangent setback = `r / tan(θ/2)`; arc centre offset = `r /
//!    sin(θ/2)` along the interior bisector.
//! 6. Reject if the setback exceeds either adjacent edge's length
//!    (`RadiusTooLarge`).
//! 7. Create new vertices `T1` (on edge_a) and `T2` (on edge_b)
//!    at the tangent points.
//! 8. Update the two half-edges with origin = corner: their
//!    origins now become `T1` and `T2` respectively, and the
//!    underlying edge curves get their origin-side endpoint
//!    shortened.
//! 9. Create a new arc edge with a `CurveGeometry::Arc` and two
//!    fresh half-edges (one per face loop). Splice both into
//!    their loops between the predecessor of the old corner
//!    half-edge and the updated half-edge.
//! 10. Delete the corner vertex (no longer referenced).
//! 11. Re-attach planar pcurves on the affected faces.
//! 12. Validate the result; on failure return `InvalidTopology`.
//!
//! The caller receives a fresh [`BRep`] on success (the input is
//! never mutated).

use popkern_core::geometry::curve::{CurveGeometry, CurveOnEdge};
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::topology::{BRep, Edge, HalfEdge, HalfEdgeHandle, VertexHandle};

use crate::error::{BlendError, BlendFailure};

/// Fillet a single sharp corner of a sheet body.
///
/// Returns a new BRep with the corner replaced by a circular arc
/// of the requested radius. The input `brep` is unchanged.
///
/// # Errors
///
/// - [`BlendError::InvalidParameter`] — `radius ≤ 0`, `radius` is
///   NaN or infinite, or `vertex` isn't a handle in `brep`.
/// - [`BlendError::PartialFailure`] with a `BlendFailure::InvalidTopology`
///   wrapping one of:
///   - corner has ≠ 2 incident edges
///   - an incident edge is not straight (not a `Line` curve)
///   - edges are colinear / near-parallel (no meaningful corner)
///   - radius is too large for the shorter adjacent edge
///   - post-surgery validation fails
pub fn fillet_corner_2d(
    brep: &BRep,
    vertex: VertexHandle,
    radius: f64,
) -> Result<BRep, BlendError> {
    if !radius.is_finite() || radius <= 0.0 {
        return Err(BlendError::InvalidParameter {
            parameter: "radius".into(),
            value: format!("{radius}"),
            reason: "radius must be a positive finite value".into(),
        });
    }
    if brep.vertices.get(vertex).is_none() {
        return Err(BlendError::InvalidParameter {
            parameter: "vertex".into(),
            value: format!("{vertex:?}"),
            reason: "vertex handle is not present in the BRep".into(),
        });
    }

    let plan =
        corner_plan(brep, vertex, radius).map_err(|reason| topology_failure(vertex, reason))?;

    let mut out = brep.clone();
    apply_corner_plan(&mut out, vertex, &plan)
        .map_err(|reason| topology_failure(vertex, reason))?;

    out.attach_planar_pcurves()
        .map_err(|e| topology_failure(vertex, format!("pcurve refresh failed: {e:?}")))?;

    out.validate()
        .map_err(|e| topology_failure(vertex, format!("validation failed: {e:?}")))?;

    Ok(out)
}

/// Pre-surgery analysis captured before we touch the mutable BRep.
/// Collecting everything up front lets us bail out cleanly on any
/// input validation failure without having to roll back partial
/// state.
#[derive(Debug, Clone, Copy)]
struct CornerPlan {
    /// Outgoing half-edge at the corner whose edge will carry the
    /// first tangent point `t1`. Walked out of the corner vertex
    /// toward the far end of `edge_a`.
    he_out_a: HalfEdgeHandle,
    /// Outgoing half-edge at the corner toward the far end of
    /// `edge_b`.
    he_out_b: HalfEdgeHandle,
    /// Tangent point on `edge_a` (distance `r / tan(θ/2)` from
    /// the corner along `dir_a`).
    t1: Point3,
    /// Tangent point on `edge_b`.
    t2: Point3,
    /// Arc centre, on the interior bisector at distance
    /// `r / sin(θ/2)` from the corner.
    arc_center: Point3,
    /// Unit arc-plane normal, `dir_a × dir_b`.
    arc_normal: Vec3,
    /// Arc start angle at `t1` relative to `arc_center`.
    start_angle: f64,
    /// Arc end angle at `t2`.
    end_angle: f64,
    /// Fillet radius (copied in for convenience).
    radius: f64,
}

fn corner_plan(brep: &BRep, vertex: VertexHandle, radius: f64) -> Result<CornerPlan, String> {
    // Collect all outgoing half-edges at the corner. A sheet-body
    // corner has two: one on the front face's outer loop, one on
    // the backside face's loop. Each points to a different edge.
    let mut outgoing: Vec<HalfEdgeHandle> = brep
        .half_edges
        .iter()
        .filter(|(_, he)| he.origin == vertex)
        .map(|(h, _)| h)
        .collect();
    outgoing.sort();

    if outgoing.len() != 2 {
        return Err(format!(
            "corner vertex must have exactly 2 outgoing half-edges, got {}",
            outgoing.len()
        ));
    }

    let corner_pos = brep.vertices[vertex].position;

    // Both outgoing half-edges should belong to distinct edges.
    let edge_a = brep.half_edges[outgoing[0]].edge;
    let edge_b = brep.half_edges[outgoing[1]].edge;
    if edge_a == edge_b {
        return Err("both outgoing half-edges belong to the same edge (degenerate loop)".into());
    }

    // Both edges must carry straight-line curves.
    require_line(brep, edge_a)?;
    require_line(brep, edge_b)?;

    // Destination vertex (the "far end") of each edge: walk the
    // half-edge's twin → origin.
    let v_a_far = brep.half_edges[brep.half_edges[outgoing[0]].twin].origin;
    let v_b_far = brep.half_edges[brep.half_edges[outgoing[1]].twin].origin;
    let p_a = brep.vertices[v_a_far].position;
    let p_b = brep.vertices[v_b_far].position;

    let dir_a_vec = p_a - corner_pos;
    let dir_b_vec = p_b - corner_pos;
    let len_a = dir_a_vec.length();
    let len_b = dir_b_vec.length();
    if len_a < 1e-12 || len_b < 1e-12 {
        return Err("adjacent edge has near-zero length".into());
    }
    let dir_a = dir_a_vec * (1.0 / len_a);
    let dir_b = dir_b_vec * (1.0 / len_b);

    let cos_angle = dir_a.dot(&dir_b).clamp(-1.0, 1.0);
    let angle = cos_angle.acos();
    if angle < 1.0e-4 {
        return Err("edges are nearly colinear (no meaningful corner)".into());
    }
    if (std::f64::consts::PI - angle).abs() < 1.0e-4 {
        return Err("edges form a near-180° corner (degenerate)".into());
    }

    let half_angle = angle * 0.5;
    let tan_half = half_angle.tan();
    let sin_half = half_angle.sin();
    let setback = radius / tan_half;
    let offset = radius / sin_half;

    if setback >= len_a || setback >= len_b {
        return Err(format!(
            "radius {radius:.6} too large: setback {setback:.6} ≥ adjacent edge length ({len_a:.6} / {len_b:.6})"
        ));
    }

    let t1 = corner_pos + dir_a * setback;
    let t2 = corner_pos + dir_b * setback;

    // Interior bisector: halfway between dir_a and dir_b. If
    // |dir_a + dir_b| is near zero (edges are anti-parallel) the
    // angle check above would already have caught it.
    let bisector = (dir_a + dir_b)
        .normalized()
        .ok_or_else(|| "bisector degenerates — edges near-parallel".to_string())?;
    let arc_center = corner_pos + bisector * offset;

    // Arc plane normal: `dir_a × dir_b` — perpendicular to the
    // corner's plane.
    let raw_normal = dir_a.cross(&dir_b);
    let arc_normal = raw_normal
        .normalized()
        .ok_or_else(|| "arc normal degenerate".to_string())?;

    // `CurveGeometry::Arc::point_at` internally uses the
    // `circle_basis(normal)` frame, which picks `u` by crossing
    // the normal with +x or +y (whichever isn't parallel) and
    // derives `v = normal × u`. We must compute our start and end
    // angles IN THAT FRAME so the evaluated arc endpoints actually
    // land on `t1` and `t2`. Replicating `circle_basis` locally
    // here avoids coupling to a private helper.
    let (u_basis, v_basis) = circle_basis_like_core(arc_normal);
    let t1_rel = t1 - arc_center;
    let t2_rel = t2 - arc_center;
    let start_angle = t1_rel.dot(&v_basis).atan2(t1_rel.dot(&u_basis));
    let end_angle = t2_rel.dot(&v_basis).atan2(t2_rel.dot(&u_basis));

    Ok(CornerPlan {
        he_out_a: outgoing[0],
        he_out_b: outgoing[1],
        t1,
        t2,
        arc_center,
        arc_normal,
        start_angle,
        end_angle,
        radius,
    })
}

/// Rebuild `edge`'s curve as a straight `Line` whose direction
/// runs from the canonical half-edge's origin vertex to its
/// twin's origin. The origin-vertex positions are read from the
/// current BRep, so this picks up any vertex moves that happened
/// earlier in the surgery.
///
/// Shared between `fillet_2d` and `chamfer_2d` — both have to
/// rebuild support-edge curves after moving their endpoints onto
/// new tangent-point vertices.
pub(crate) fn rebuild_line_curve(brep: &mut BRep, edge: popkern_core::topology::EdgeHandle) {
    let canonical = brep.edges[edge].half_edge;
    let twin = brep.half_edges[canonical].twin;
    let v_from = brep.half_edges[canonical].origin;
    let v_to = brep.half_edges[twin].origin;
    let p_from = brep.vertices[v_from].position;
    let p_to = brep.vertices[v_to].position;
    brep.edges[edge].curve = Some(CurveOnEdge {
        geometry: CurveGeometry::Line {
            origin: p_from,
            direction: p_to - p_from,
        },
        t_start: 0.0,
        t_end: 1.0,
    });
}

/// Replicate `popkern_core::geometry::analytic_curves::circle_basis`
/// — which is private — so our computed start/end angles match the
/// frame the core's Arc evaluation actually uses at runtime.
fn circle_basis_like_core(normal: Vec3) -> (Vec3, Vec3) {
    let axis = if normal.x.abs() < 0.9 {
        Vec3::new(1.0, 0.0, 0.0)
    } else {
        Vec3::new(0.0, 1.0, 0.0)
    };
    let perp = normal.cross(&axis);
    let u = perp.normalized().unwrap_or(Vec3::new(1.0, 0.0, 0.0));
    let v = normal
        .cross(&u)
        .normalized()
        .unwrap_or(Vec3::new(0.0, 1.0, 0.0));
    (u, v)
}

/// Validate that the given edge's curve is a straight `Line`.
/// Used by both `fillet_corner_2d` and `chamfer_corner_2d` to
/// guard against curved boundary edges that neither currently
/// handles.
pub(crate) fn require_line(
    brep: &BRep,
    edge_h: popkern_core::topology::EdgeHandle,
) -> Result<(), String> {
    let edge = brep
        .edges
        .get(edge_h)
        .ok_or_else(|| format!("edge {edge_h:?} not found"))?;
    let curve = edge
        .curve
        .as_ref()
        .ok_or_else(|| format!("edge {edge_h:?} has no curve geometry"))?;
    match &curve.geometry {
        CurveGeometry::Line { .. } => Ok(()),
        other => Err(format!(
            "edge {edge_h:?} curve must be a Line, got {other:?}"
        )),
    }
}

fn apply_corner_plan(
    brep: &mut BRep,
    corner: VertexHandle,
    plan: &CornerPlan,
) -> Result<(), String> {
    // 1. Insert the new tangent-point vertices.
    let v_t1 = brep.vertices.insert(popkern_core::topology::Vertex {
        position: plan.t1,
        outgoing: Default::default(),
        name: None,
        tolerance: None,
    });
    let v_t2 = brep.vertices.insert(popkern_core::topology::Vertex {
        position: plan.t2,
        outgoing: Default::default(),
        name: None,
        tolerance: None,
    });

    // 2. Move the origin of the outgoing half-edges onto the new
    //    tangent vertices. `he_out_a` now starts at `t1`, and its
    //    underlying edge's curve is a line from `t1` to `v_a_far`.
    //    Same for `he_out_b` with `t2` and `v_b_far`.
    brep.half_edges[plan.he_out_a].origin = v_t1;
    brep.half_edges[plan.he_out_b].origin = v_t2;

    // Rebuild the two support edges' curves with the shorter
    // Line geometry. The curve must run from the canonical
    // half-edge's origin to its twin's origin — which may or may
    // not be the direction we just updated, depending on which of
    // the edge's two half-edges was selected as canonical when
    // the edge was first created.
    let edge_a_h = brep.half_edges[plan.he_out_a].edge;
    let edge_b_h = brep.half_edges[plan.he_out_b].edge;
    rebuild_line_curve(brep, edge_a_h);
    rebuild_line_curve(brep, edge_b_h);

    // 3. Create the arc edge and its two half-edges (one per face
    //    loop). `new_front` splices into the loop containing
    //    `he_out_a` (just before it) and points from `t2` → `t1`.
    //    `new_back` splices into the loop containing `he_out_b`
    //    and points from `t1` → `t2`.
    let arc_geometry = CurveGeometry::Arc {
        center: plan.arc_center,
        normal: plan.arc_normal,
        radius: plan.radius,
        start_angle: plan.start_angle,
        end_angle: plan.end_angle,
    };
    let arc_edge = brep.edges.insert(Edge {
        half_edge: Default::default(),
        curve: Some(CurveOnEdge {
            geometry: arc_geometry,
            t_start: plan.start_angle,
            t_end: plan.end_angle,
        }),
        name: None,
        tolerance: None,
        canonical: Default::default(),
    });

    // Gather the loop and predecessor data we need BEFORE inserting
    // the new half-edges, because both belong to the same half-edge
    // slot-map and we mustn't borrow while inserting.
    let loop_front = brep.half_edges[plan.he_out_a].loop_;
    let loop_back = brep.half_edges[plan.he_out_b].loop_;
    let prev_front = brep.half_edges[plan.he_out_a].prev;
    let prev_back = brep.half_edges[plan.he_out_b].prev;

    // Front half-edge: origin `t2`, destination `t1`. Inserted
    // before `he_out_a` in the front loop.
    let new_front = brep.half_edges.insert(HalfEdge {
        origin: v_t2,
        twin: Default::default(),
        next: plan.he_out_a,
        prev: prev_front,
        edge: arc_edge,
        loop_: loop_front,
        pcurve: None,
    });
    // Back half-edge: origin `t1`, destination `t2`. Inserted
    // before `he_out_b` in the back loop.
    let new_back = brep.half_edges.insert(HalfEdge {
        origin: v_t1,
        twin: new_front,
        next: plan.he_out_b,
        prev: prev_back,
        edge: arc_edge,
        loop_: loop_back,
        pcurve: None,
    });
    brep.half_edges[new_front].twin = new_back;
    // Canonical half-edge direction matches the arc curve's
    // parameterisation: start_angle → t1, end_angle → t2. The
    // half-edge that starts at t1 is `new_back` (origin = v_t1).
    brep.edges[arc_edge].half_edge = new_back;

    // 4. Rewire next/prev across both splices.
    brep.half_edges[prev_front].next = new_front;
    brep.half_edges[plan.he_out_a].prev = new_front;
    brep.half_edges[prev_back].next = new_back;
    brep.half_edges[plan.he_out_b].prev = new_back;

    // 5. Make sure each vertex's `outgoing` handle points at a
    //    valid half-edge. The new tangent vertices emit the
    //    modified outgoing half-edges; the old corner vertex loses
    //    its references.
    brep.vertices[v_t1].outgoing = plan.he_out_a;
    brep.vertices[v_t2].outgoing = plan.he_out_b;

    // 6. Any loop anchor that pointed at `he_out_a` / `he_out_b`
    //    is still valid (the handles are unchanged), but anchors
    //    that happened to point at the front/back loop via the
    //    predecessor are also valid. Nothing to do.

    // 7. Delete the corner vertex. It is no longer referenced —
    //    any surviving reference would indicate a wiring bug and
    //    validation will catch it.
    brep.vertices.remove(corner);

    Ok(())
}

/// Build a uniform `PartialFailure` wrapping an `InvalidTopology`
/// blend failure for a 2D corner operation. Shared between
/// fillet and chamfer so both emit the same error shape.
pub(crate) fn topology_failure(vertex: VertexHandle, reason: String) -> BlendError {
    BlendError::PartialFailure {
        failed: 1,
        total: 1,
        partial_result: Box::new(BRep::new(popkern_core::tolerance::Tolerance::default())),
        failures: vec![(
            popkern_core::topology::EdgeHandle::default(),
            BlendFailure::InvalidTopology {
                reason: format!("fillet_corner_2d at {vertex:?}: {reason}"),
                location: None,
            },
        )],
    }
}

#[cfg(test)]
#[path = "fillet_2d_tests.rs"]
mod tests;
