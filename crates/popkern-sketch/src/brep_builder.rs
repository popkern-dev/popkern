// Items are pub(crate) for use by the pipeline module.
// Suppress dead_code until those consumers are added.
#![allow(dead_code)]

//! BRep construction from classified nesting trees.
//!
//! Converts the nesting tree produced by [`crate::nesting::classify_nesting`]
//! into a valid BRep sheet body. Each face-flagged node produces a
//! two-sided planar face (front + back) with edge curves and pcurves.
//! Hole-flagged children become inner loops on the front face.

use crate::curve_convert;
use crate::error::SketchError;
use crate::loop_tracer::TracedLoop;
use crate::nesting::NestingNode;
use crate::plane::SketchPlane;
use popkern_constraint::Sketch;
use popkern_core::constructors::BRepBuilder;
use popkern_core::geometry::surface::{SurfaceGeometry, SurfaceOnFace};
use popkern_core::topology::{FaceHandle, HalfEdgeHandle, LoopHandle, VertexHandle};

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
    shell: popkern_core::topology::ShellHandle,
    surface: SurfaceOnFace,
    back_surface: SurfaceOnFace,
}

/// Build a BRep from a nesting tree, traced loops, and sketch geometry.
///
/// Creates a two-sided sheet body for each face-flagged nesting node.
/// All regions share a single BRep. Returns the BRep and the list of
/// built regions (one per face-flagged node).
///
/// # Errors
///
/// Returns [`SketchError::InvalidTopology`] if the builder fails or
/// the loop data is inconsistent.
pub(crate) fn build_brep(
    nesting_roots: &[NestingNode],
    loops: &[TracedLoop],
    sketch: &Sketch,
    plane: &SketchPlane,
    tolerance: &popkern_core::tolerance::Tolerance,
) -> Result<(popkern_core::topology::BRep, Vec<BuiltRegion>), SketchError> {
    let mut builder = BRepBuilder::new(*tolerance);
    let mut regions = Vec::new();

    let solid = builder.add_solid();
    let shell = builder.add_shell(solid);

    let ctx = BuildContext {
        loops,
        sketch,
        plane,
        shell,
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
            let region = build_face_region(&mut builder, node, &ctx)?;
            regions.push(region);
        }
    }

    let brep = builder.build_unchecked();
    Ok((brep, regions))
}

/// Build a single face region from a face-flagged nesting node.
///
/// Creates front/back faces for the outer loop and adds inner loops
/// for hole-flagged children.
fn build_face_region(
    builder: &mut BRepBuilder,
    node: &NestingNode,
    ctx: &BuildContext<'_>,
) -> Result<BuiltRegion, SketchError> {
    let traced = &ctx.loops[node.loop_index];

    // Create front and back faces.
    let front_face = builder.add_face(ctx.shell);
    let back_face = builder.add_face(ctx.shell);

    builder.brep_mut().faces[front_face].surface = Some(ctx.surface.clone());
    builder.brep_mut().faces[back_face].surface = Some(ctx.back_surface.clone());

    // Build outer loop: vertices, edges, half-edge chains.
    build_loop_edges(builder, traced, ctx, front_face, back_face, false)?;

    // Build inner loops (holes) from children that are NOT faces.
    let mut holes = Vec::new();
    for child in &node.children {
        if !child.is_face {
            let hole_traced = &ctx.loops[child.loop_index];
            let hole_loop =
                build_loop_edges(builder, hole_traced, ctx, front_face, back_face, true)?;
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

/// Build edges and half-edge chains for a single traced loop.
///
/// When `is_hole` is false, wires the half-edges as outer loops on the
/// front and back faces. When `is_hole` is true, wires them as inner
/// loops and returns the front face's inner loop handle.
fn build_loop_edges(
    builder: &mut BRepBuilder,
    traced: &TracedLoop,
    ctx: &BuildContext<'_>,
    front_face: FaceHandle,
    back_face: FaceHandle,
    is_hole: bool,
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
        .map(|de| builder.add_vertex(ctx.plane.to_3d(de.start.0, de.start.1)))
        .collect();

    // Create edges and attach geometry.
    let mut front_hes: Vec<HalfEdgeHandle> = Vec::with_capacity(n);
    let mut back_hes: Vec<HalfEdgeHandle> = Vec::with_capacity(n);

    for i in 0..n {
        let j = (i + 1) % n;
        let (edge_h, he_fwd, he_rev) = builder.add_edge(vertices[i], vertices[j]);

        // Attach 3D curve.
        let curve = curve_convert::convert_to_3d(&traced.edges[i], ctx.sketch, ctx.plane)?;
        builder.brep_mut().edges[edge_h].curve = Some(curve);

        // Attach pcurves.
        let pcurve = curve_convert::convert_to_pcurve(&traced.edges[i], ctx.sketch)?;
        builder.brep_mut().half_edges[he_fwd].pcurve = Some(pcurve);

        let back_pcurve = curve_convert::convert_to_pcurve(&traced.edges[i], ctx.sketch)?;
        builder.brep_mut().half_edges[he_rev].pcurve = Some(back_pcurve);

        front_hes.push(he_fwd);
        back_hes.push(he_rev);
    }

    // Back face uses reversed twin half-edges.
    back_hes.reverse();

    if is_hole {
        let hole_loop = builder
            .add_inner_loop(&front_hes, front_face)
            .map_err(|e| SketchError::InvalidTopology(format!("hole loop: {e}")))?;

        builder
            .add_inner_loop(&back_hes, back_face)
            .map_err(|e| SketchError::InvalidTopology(format!("back hole loop: {e}")))?;

        Ok(Some(hole_loop))
    } else {
        builder
            .add_loop(&front_hes, front_face)
            .map_err(|e| SketchError::InvalidTopology(format!("front loop: {e}")))?;

        builder
            .add_loop(&back_hes, back_face)
            .map_err(|e| SketchError::InvalidTopology(format!("back loop: {e}")))?;

        Ok(None)
    }
}

#[cfg(test)]
#[path = "brep_builder_tests.rs"]
mod tests;
