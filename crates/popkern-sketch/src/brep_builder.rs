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
    surface: SurfaceOnFace,
    back_surface: SurfaceOnFace,
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
        back_surface: SurfaceOnFace {
            geometry: SurfaceGeometry::Plane {
                origin: plane.origin,
                normal: -plane.normal,
                u_axis: plane.v_axis,
                v_axis: plane.u_axis,
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
    let (front_face, back_face) = build_outer_loop(brep, traced, ctx)?;

    // Attach surfaces.
    brep.faces[front_face].surface = Some(ctx.surface.clone());
    brep.faces[back_face].surface = Some(ctx.back_surface.clone());

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

/// Build the outer loop of a face using Euler operators (mvfs, mev, mef).
///
/// Returns `(front_face, back_face)` where `front_face` is the mef face
/// (the enclosed polygon) and `back_face` is the original mvfs face.
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
    let front_face = mef.face;
    let back_face = mvfs.face;

    // Attach curves and pcurves to the edges.
    attach_loop_geometry(brep, traced, ctx, front_face)?;
    attach_twin_pcurves(brep, traced, ctx, back_face)?;

    Ok((front_face, back_face))
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
        let origin_pos = brep.vertices[origin].position;
        let origin_2d = ctx.plane.to_2d(origin_pos);

        if let Some(idx) = find_matching_edge(traced, origin_2d, n) {
            let directed = &traced.edges[idx];
            let edge_handle = brep.half_edges[cursor].edge;

            // Attach 3D curve.
            if brep.edges[edge_handle].curve.is_none() {
                let curve = curve_convert::convert_to_3d(directed, ctx.sketch, ctx.plane)?;
                brep.edges[edge_handle].curve = Some(curve);
            }

            // Attach pcurve.
            let pcurve = curve_convert::convert_to_pcurve(directed, ctx.sketch)?;
            brep.half_edges[cursor].pcurve = Some(pcurve);
        }

        cursor = brep.half_edges[cursor].next;
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
        if brep.half_edges[cursor].pcurve.is_none() {
            let origin = brep.half_edges[cursor].origin;
            let origin_pos = brep.vertices[origin].position;
            let origin_2d = ctx.plane.to_2d(origin_pos);

            if let Some(idx) = find_matching_edge(traced, origin_2d, n) {
                let directed = &traced.edges[idx];
                let pcurve = curve_convert::convert_to_pcurve(directed, ctx.sketch)?;
                brep.half_edges[cursor].pcurve = Some(pcurve);
            }
        }

        cursor = brep.half_edges[cursor].next;
        if cursor == start_he {
            break;
        }
    }

    Ok(())
}

/// Find the traced edge whose start position matches the given 2D point.
fn find_matching_edge(traced: &TracedLoop, pos_2d: (f64, f64), n: usize) -> Option<usize> {
    for i in 0..n {
        let dx = traced.edges[i].start.0 - pos_2d.0;
        let dy = traced.edges[i].start.1 - pos_2d.1;
        if (dx * dx + dy * dy).sqrt() < 1e-6 {
            return Some(i);
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
        brep.half_edges[he_fwd].pcurve = Some(pcurve);

        let back_pcurve = curve_convert::convert_to_pcurve(&traced.edges[i], ctx.sketch)?;
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
