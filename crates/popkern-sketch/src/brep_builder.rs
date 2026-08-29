// Items are pub(crate) for use by the pipeline module.
// Suppress dead_code until those consumers are added.
#![allow(dead_code)]

//! BRep construction from classified nesting trees using Euler operators.
//!
//! Converts the nesting tree produced by [`crate::nesting::classify_nesting`]
//! into a valid BRep using Euler operators (mvfs, mev, mef) for the outer
//! loop and direct inner-loop wiring for holes. Each face-flagged node
//! produces a two-sided planar face suitable for extrusion.

use crate::curve_convert;
use crate::error::SketchError;
use crate::loop_tracer::TracedLoop;
use crate::nesting::NestingNode;
use crate::plane::SketchPlane;
use popkern_constraint::Sketch;
use popkern_core::geometry::surface::{SurfaceGeometry, SurfaceOnFace};
use popkern_core::topology::{
    Edge, FaceHandle, HalfEdge, HalfEdgeHandle, Loop, LoopHandle, VertexHandle,
};

/// A constructed face region with its associated inner loops (holes).
#[derive(Debug, Clone)]
pub(crate) struct BuiltRegion {
    /// Front face handle in the BRep.
    pub face: FaceHandle,
    /// Inner loop handles (holes) on the front face.
    pub holes: Vec<LoopHandle>,
}

/// Shared context for building a sketch BRep.
struct BuildContext<'a> {
    loops: &'a [TracedLoop],
    sketch: &'a Sketch,
    plane: &'a SketchPlane,
    /// The sketch plane as a face surface, `orientation: true`. Both
    /// lamina faces share this chart; the back face carries the same
    /// geometry with `orientation: false` (PK-376) so every pcurve is
    /// written once, in one chart, and the flag alone encodes the side.
    surface: SurfaceOnFace,
}

/// Build a BRep from a nesting tree, traced loops, and sketch geometry.
///
/// Uses Euler operators (mvfs, mev, mef) for the outer loop topology
/// and direct inner-loop wiring for holes. Returns the BRep and the
/// list of built regions (one per face-flagged node).
///
/// # Errors
///
/// Returns [`SketchError::InvalidTopology`] if an Euler operation fails
/// or the loop data is inconsistent.
pub(crate) fn build_brep(
    nesting_roots: &[NestingNode],
    loops: &[TracedLoop],
    sketch: &Sketch,
    plane: &SketchPlane,
    tolerance: &popkern_core::tolerance::Tolerance,
) -> Result<(popkern_core::topology::BRep, Vec<BuiltRegion>), SketchError> {
    let mut brep = popkern_core::topology::BRep::new(*tolerance);
    let mut regions = Vec::new();

    let ctx = BuildContext {
        loops,
        sketch,
        plane,
        surface: SurfaceOnFace {
            geometry: SurfaceGeometry::Plane {
                origin: plane.origin,
                normal: plane.normal,
                u_axis: plane.u_axis,
                v_axis: plane.v_axis,
            },
            orientation: true,
        },
    };

    for node in nesting_roots {
        if node.is_face {
            let region = build_face_region(&mut brep, node, &ctx)?;
            regions.push(region);
        }
    }

    Ok((brep, regions))
}

/// Build a single face region from a face-flagged nesting node.
///
/// Creates an outer loop face using Euler operators and adds inner
/// loops for hole-flagged children.
fn build_face_region(
    brep: &mut popkern_core::topology::BRep,
    node: &NestingNode,
    ctx: &BuildContext<'_>,
) -> Result<BuiltRegion, SketchError> {
    let traced = &ctx.loops[node.loop_index];

    // Build outer loop with Euler operators.
    let (mef_face, mvfs_face) = build_outer_loop(brep, traced, ctx)?;

    // The orientation flag is DERIVED state (PK-244): the face whose
    // outer loop winds counter-clockwise about the sketch normal faces
    // `+normal` (`orientation: true`); its twin faces `-normal`. The
    // Euler-op role (mef vs mvfs) says nothing about the winding — it
    // follows the tracer's direction — so assigning by role handed
    // BOTH faces a flag contradicting their loop (measured: a
    // rectangle lamina on XY, both faces `WindingDishonest`, and every
    // sketch→extrude body carried a dishonest base cap; PK-376).
    let (front_face, back_face) = if loop_winding_about(brep, mef_face, ctx.plane.normal) >= 0.0 {
        (mef_face, mvfs_face)
    } else {
        (mvfs_face, mef_face)
    };
    brep.faces[front_face].surface = Some(ctx.surface.clone());
    brep.faces[back_face].surface = Some(SurfaceOnFace {
        geometry: ctx.surface.geometry.clone(),
        orientation: false,
    });

    // Build inner loops (holes) from children that are NOT faces.
    let mut holes = Vec::new();
    for child in &node.children {
        if !child.is_face {
            let hole_traced = &ctx.loops[child.loop_index];
            let hole_loop = build_hole_loop(brep, hole_traced, ctx, front_face, back_face)?;
            if let Some(lh) = hole_loop {
                holes.push(lh);
            }
        }
    }

    Ok(BuiltRegion {
        face: front_face,
        holes,
    })
}

/// Signed winding of `face`'s outer loop about `normal`: the Newell
/// normal of the loop's vertex polygon dotted with `normal` (positive =
/// counter-clockwise seen from `+normal`).
fn loop_winding_about(
    brep: &popkern_core::topology::BRep,
    face: FaceHandle,
    normal: popkern_core::geometry::Vec3,
) -> f64 {
    let pts = brep.face_outer_loop_positions(face);
    let mut newell = popkern_core::geometry::Vec3::new(0.0, 0.0, 0.0);
    for i in 0..pts.len() {
        let a = pts[i];
        let b = pts[(i + 1) % pts.len()];
        newell.x += (a.y - b.y) * (a.z + b.z);
        newell.y += (a.z - b.z) * (a.x + b.x);
        newell.z += (a.x - b.x) * (a.y + b.y);
    }
    newell.dot(&normal)
}

/// Build the outer loop of a face using Euler operators (mvfs, mev, mef).
///
/// Returns `(mef_face, mvfs_face)`: the face `mef` enclosed and the
/// original `mvfs` face. Which one faces the sketch normal is decided
/// by the caller from the measured loop winding.
fn build_outer_loop(
    brep: &mut popkern_core::topology::BRep,
    traced: &TracedLoop,
    ctx: &BuildContext<'_>,
) -> Result<(FaceHandle, FaceHandle), SketchError> {
    let n = traced.edges.len();
    if n < 2 {
        return Err(SketchError::InvalidTopology(
            "loop must have at least 2 edges".into(),
        ));
    }

    let map_err = |e| SketchError::InvalidTopology(format!("Euler op failed: {e}"));

    // 1. mvfs at first vertex.
    let first_pos = ctx
        .plane
        .to_3d(traced.edges[0].start.0, traced.edges[0].start.1);
    let mvfs = brep.mvfs(first_pos).map_err(map_err)?;

    // 2. mev_at_vertex to create the second vertex.
    let second_pos = ctx
        .plane
        .to_3d(traced.edges[1].start.0, traced.edges[1].start.1);
    let mev1 = brep
        .mev_at_vertex(mvfs.vertex, mvfs.loop_, second_pos)
        .map_err(map_err)?;

    // 3. mev chain for remaining vertices.
    let mut last_twin = mev1.half_edge_twin;
    for i in 2..n {
        let pos = ctx
            .plane
            .to_3d(traced.edges[i].start.0, traced.edges[i].start.1);
        let mev = brep.mev(last_twin, pos).map_err(map_err)?;
        last_twin = mev.half_edge_twin;
    }

    // 4. mef to close the loop and create the face.
    let mef = brep.mef(last_twin, mev1.half_edge).map_err(map_err)?;
    let mef_face = mef.face;
    let mvfs_face = mvfs.face;

    // Attach curves and pcurves to the edges (one shared chart, so the
    // side each face ends up on does not matter here).
    attach_loop_geometry(brep, traced, ctx, mef_face)?;
    attach_twin_pcurves(brep, traced, ctx, mvfs_face)?;

    Ok((mef_face, mvfs_face))
}

/// Attach 3D curves and pcurves to the half-edges around a face's outer loop.
fn attach_loop_geometry(
    brep: &mut popkern_core::topology::BRep,
    traced: &TracedLoop,
    ctx: &BuildContext<'_>,
    face: FaceHandle,
) -> Result<(), SketchError> {
    let face_loop = brep.faces[face].outer_loop;
    let start_he = brep.loops[face_loop].half_edge;
    let n = traced.edges.len();

    let mut cursor = start_he;
    loop {
        let origin = brep.half_edges[cursor].origin;
        let next = brep.half_edges[cursor].next;
        let origin_2d = ctx.plane.to_2d(brep.vertices[origin].position);
        let next_origin = brep.half_edges[next].origin;
        let next_2d = ctx.plane.to_2d(brep.vertices[next_origin].position);

        // The half-edge's SEGMENT is origin → next origin; match the
        // traced edge by BOTH endpoints. Start-only matching attached
        // the sketch's OTHER edge leaving the origin whenever a loop's
        // winding ran opposite the traced direction (measured: every
        // curve one segment off on the replay sheet — validate red
        // with a full-edge-length loop gap on both faces, PK-317).
        if let Some((idx, forward)) = find_matching_segment(traced, origin_2d, next_2d, n) {
            let directed = &traced.edges[idx];
            let edge_handle = brep.half_edges[cursor].edge;

            // Attach 3D curve (canonical per edge; consumers resolve
            // the sense through the half-edge).
            if brep.edges[edge_handle].curve.is_none() {
                let curve = curve_convert::convert_to_3d(directed, ctx.sketch, ctx.plane)?;
                brep.edges[edge_handle].curve = Some(curve);
            }

            // Attach pcurve FOLLOWING THE HALF-EDGE: a half-edge that
            // runs opposite the traced direction stores the reversed
            // span (t-swap — the canonical-direction convention,
            // PK-303).
            let pcurve = curve_convert::convert_to_pcurve(directed, ctx.sketch)?;
            let mut pcurve: popkern_core::geometry::pcurve::Pcurve = pcurve.into();
            if !forward {
                std::mem::swap(&mut pcurve.t_start, &mut pcurve.t_end);
            }
            brep.half_edges[cursor].pcurve = Some(pcurve);
        }

        cursor = next;
        if cursor == start_he {
            break;
        }
    }

    Ok(())
}

/// Attach pcurves to twin half-edges (on the back face's outer loop).
fn attach_twin_pcurves(
    brep: &mut popkern_core::topology::BRep,
    traced: &TracedLoop,
    ctx: &BuildContext<'_>,
    back_face: FaceHandle,
) -> Result<(), SketchError> {
    let face_loop = brep.faces[back_face].outer_loop;
    let start_he = brep.loops[face_loop].half_edge;
    let n = traced.edges.len();

    let mut cursor = start_he;
    loop {
        let next = brep.half_edges[cursor].next;
        if brep.half_edges[cursor].pcurve.is_none() {
            let origin = brep.half_edges[cursor].origin;
            let origin_2d = ctx.plane.to_2d(brep.vertices[origin].position);
            let next_origin = brep.half_edges[next].origin;
            let next_2d = ctx.plane.to_2d(brep.vertices[next_origin].position);

            if let Some((idx, forward)) = find_matching_segment(traced, origin_2d, next_2d, n) {
                let directed = &traced.edges[idx];
                let pcurve = curve_convert::convert_to_pcurve(directed, ctx.sketch)?;
                let mut pcurve: popkern_core::geometry::pcurve::Pcurve = pcurve.into();
                if !forward {
                    std::mem::swap(&mut pcurve.t_start, &mut pcurve.t_end);
                }
                brep.half_edges[cursor].pcurve = Some(pcurve);
            }
        }

        cursor = next;
        if cursor == start_he {
            break;
        }
    }

    Ok(())
}

/// Find the traced edge whose start position matches the given 2D point.
/// Match a half-edge SEGMENT (both endpoints) to a traced directed
/// edge. Returns the traced index and whether the half-edge runs the
/// traced (forward) direction. A directed edge's end is the next
/// traced edge's start (the loop is a directed cycle).
fn find_matching_segment(
    traced: &TracedLoop,
    a_2d: (f64, f64),
    b_2d: (f64, f64),
    n: usize,
) -> Option<(usize, bool)> {
    let close = |p: (f64, f64), q: (f64, f64)| -> bool {
        let dx = p.0 - q.0;
        let dy = p.1 - q.1;
        (dx * dx + dy * dy).sqrt() < 1e-6
    };
    for i in 0..n {
        let start = traced.edges[i].start;
        let end = traced.edges[(i + 1) % n].start;
        if close(start, a_2d) && close(end, b_2d) {
            return Some((i, true));
        }
        if close(start, b_2d) && close(end, a_2d) {
            return Some((i, false));
        }
    }
    None
}

/// Build inner loop edges and wire them as a hole on the front and back faces.
///
/// Uses direct topology construction for inner loops (the standard
/// approach for holes that avoids the complexity of Euler-operator-based
/// hole creation via mev+mef+kemr chains).
fn build_hole_loop(
    brep: &mut popkern_core::topology::BRep,
    traced: &TracedLoop,
    ctx: &BuildContext<'_>,
    front_face: FaceHandle,
    back_face: FaceHandle,
) -> Result<Option<LoopHandle>, SketchError> {
    let n = traced.edges.len();
    if n < 2 {
        return Err(SketchError::InvalidTopology(
            "loop must have at least 2 edges".into(),
        ));
    }

    // Create vertices at each directed edge's start position.
    let vertices: Vec<VertexHandle> = traced
        .edges
        .iter()
        .map(|de| {
            brep.vertices.insert(popkern_core::topology::Vertex {
                position: ctx.plane.to_3d(de.start.0, de.start.1),
                outgoing: Default::default(),
                name: None,
                tolerance: None,
            })
        })
        .collect();

    // Create edges and half-edge pairs.
    let mut front_hes: Vec<HalfEdgeHandle> = Vec::with_capacity(n);
    let mut back_hes: Vec<HalfEdgeHandle> = Vec::with_capacity(n);

    for i in 0..n {
        let j = (i + 1) % n;
        let edge_h = brep.edges.insert(Edge {
            half_edge: Default::default(),
            curve: None,
            name: None,
            tolerance: None,
            canonical: Default::default(),
            role: Default::default(),
        });

        let he_fwd = brep.half_edges.insert(HalfEdge {
            origin: vertices[i],
            twin: Default::default(),
            next: Default::default(),
            prev: Default::default(),
            edge: edge_h,
            loop_: Default::default(),
            pcurve: None,
        });

        let he_rev = brep.half_edges.insert(HalfEdge {
            origin: vertices[j],
            twin: he_fwd,
            next: Default::default(),
            prev: Default::default(),
            edge: edge_h,
            loop_: Default::default(),
            pcurve: None,
        });

        brep.half_edges[he_fwd].twin = he_rev;
        brep.edges[edge_h].half_edge = he_fwd;
        if brep.vertices[vertices[i]].outgoing == HalfEdgeHandle::default() {
            brep.vertices[vertices[i]].outgoing = he_fwd;
        }
        if brep.vertices[vertices[j]].outgoing == HalfEdgeHandle::default() {
            brep.vertices[vertices[j]].outgoing = he_rev;
        }

        // Attach 3D curve.
        let curve = curve_convert::convert_to_3d(&traced.edges[i], ctx.sketch, ctx.plane)?;
        brep.edges[edge_h].curve = Some(curve);

        // Attach pcurves.
        let pcurve = curve_convert::convert_to_pcurve(&traced.edges[i], ctx.sketch)?;
        brep.half_edges[he_fwd].pcurve = Some(pcurve.into());

        // The twin runs j → i: the same chart span, reversed
        // (t-swap — the canonical-direction convention, PK-303).
        let back_pcurve = curve_convert::convert_to_pcurve(&traced.edges[i], ctx.sketch)?;
        let mut back_pcurve: popkern_core::geometry::pcurve::Pcurve = back_pcurve.into();
        std::mem::swap(&mut back_pcurve.t_start, &mut back_pcurve.t_end);
        brep.half_edges[he_rev].pcurve = Some(back_pcurve);

        front_hes.push(he_fwd);
        back_hes.push(he_rev);
    }

    // Back face uses reversed twin half-edges.
    back_hes.reverse();

    // Wire front hole as inner loop.
    let front_loop = wire_inner_loop(brep, &front_hes, front_face)?;

    // Wire back hole as inner loop.
    wire_inner_loop(brep, &back_hes, back_face)?;

    Ok(Some(front_loop))
}

/// Wire a sequence of half-edges as an inner loop on a face.
fn wire_inner_loop(
    brep: &mut popkern_core::topology::BRep,
    half_edges: &[HalfEdgeHandle],
    face: FaceHandle,
) -> Result<LoopHandle, SketchError> {
    if half_edges.is_empty() {
        return Err(SketchError::InvalidTopology(
            "inner loop requires at least one half-edge".into(),
        ));
    }

    let lh = brep.loops.insert(Loop {
        half_edge: half_edges[0],
        face,
    });

    let n = half_edges.len();
    for i in 0..n {
        let curr = half_edges[i];
        let next = half_edges[(i + 1) % n];
        let prev = half_edges[(i + n - 1) % n];
        brep.half_edges[curr].next = next;
        brep.half_edges[curr].prev = prev;
        brep.half_edges[curr].loop_ = lh;
    }

    brep.faces[face].inner_loops.push(lh);
    Ok(lh)
}

#[cfg(test)]
#[path = "brep_builder_tests.rs"]
mod tests;

#[cfg(test)]
#[path = "brep_builder_winding_tests.rs"]
mod winding_tests;
