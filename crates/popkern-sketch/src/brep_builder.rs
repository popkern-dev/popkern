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
    let mef_outer = brep.faces[mef_face].outer_loop;
    let (front_face, back_face) = if loop_winding_about(brep, mef_outer, ctx.plane.normal) >= 0.0 {
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

/// Signed winding of loop `lh` about `normal`: the Newell
/// normal of the loop's SAMPLED outline dotted with `normal` (positive
/// = counter-clockwise seen from `+normal`). Every half-edge's curve is
/// sampled in traversal direction — the vertex polygon alone is
/// degenerate for a bigon (a split circle has two vertices and zero
/// polygon area), which assigned the circle lamina's faces by a coin
/// toss (both `WindingDishonest`).
fn loop_winding_about(
    brep: &popkern_core::topology::BRep,
    lh: LoopHandle,
    normal: popkern_core::geometry::Vec3,
) -> f64 {
    use popkern_core::geometry::curve::CurveEval;
    const SAMPLES: usize = 16;
    let mut pts: Vec<popkern_core::geometry::Point3> = Vec::new();
    let start_he = brep.loops[lh].half_edge;
    let mut cursor = start_he;
    loop {
        let he = &brep.half_edges[cursor];
        let origin = brep.vertices[he.origin].position;
        match brep.edges[he.edge].curve.as_ref() {
            Some(c) => {
                // Run the curve from the half-edge's origin.
                // Sample strictly inside the curve's own range: a
                // malformed edge range must degrade to its vertex,
                // never panic.
                let (r0, r1) = c.geometry.parameter_range();
                let (lo, hi) = (r0.min(r1), r0.max(r1));
                let clamp = |t: f64| {
                    if lo.is_finite() && hi.is_finite() {
                        t.clamp(lo, hi)
                    } else {
                        t
                    }
                };
                let (ta, tb) = (clamp(c.t_start), clamp(c.t_end));
                let forward = c.geometry.point_at(ta).distance_to(&origin)
                    <= c.geometry.point_at(tb).distance_to(&origin);
                for i in 0..SAMPLES {
                    let f = i as f64 / SAMPLES as f64;
                    let t = if forward {
                        ta + (tb - ta) * f
                    } else {
                        tb + (ta - tb) * f
                    };
                    pts.push(c.geometry.point_at(t));
                }
            }
            None => pts.push(origin),
        }
        cursor = he.next;
        if cursor == start_he {
            break;
        }
    }
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

    // One verified Euler transaction: mvfs, an mev per further vertex,
    // the closing mef (core `closed_lamina`, PK-709 — per-step checks
    // made a 1000-edge profile quadratic).
    //
    // Vertex k is traced edge k's start; edge k runs vertex k → k+1
    // and IS traced edge k. Recorded here so geometry attaches by
    // construction — never by matching endpoints, which is ambiguous
    // for a bigon (a split circle's two halves share both endpoints;
    // endpoint matching handed both edges the SAME half, validate red
    // with a √2 loop gap on every circle sketch).
    let points: Vec<_> = traced
        .edges
        .iter()
        .map(|e| ctx.plane.to_3d(e.start.0, e.start.1))
        .collect();
    let lamina = brep.closed_lamina(&points).map_err(map_err)?;
    let (vertices, edges) = (lamina.vertices, lamina.edges);

    attach_edge_geometry(brep, traced, ctx, &vertices, &edges)?;

    Ok((lamina.mef_face, lamina.mvfs_face))
}

/// Attach the 3-D curve and both half-edge pcurves of every outer-loop
/// edge from its traced edge (`edges[k]` ↔ `traced.edges[k]`, running
/// `vertices[k]` → `vertices[k + 1]`). The half-edge leaving
/// `vertices[k]` follows the traced direction and takes the pcurve
/// as-is; its twin stores the reversed span (t-swap — the
/// canonical-direction convention, PK-303). Both lamina faces share
/// one chart, so the side each face ends up on does not matter here.
fn attach_edge_geometry(
    brep: &mut popkern_core::topology::BRep,
    traced: &TracedLoop,
    ctx: &BuildContext<'_>,
    vertices: &[VertexHandle],
    edges: &[popkern_core::topology::EdgeHandle],
) -> Result<(), SketchError> {
    for (k, (&edge_handle, directed)) in edges.iter().zip(&traced.edges).enumerate() {
        let curve = curve_convert::convert_to_3d(directed, ctx.sketch, ctx.plane)?;
        brep.edges[edge_handle].curve = Some(curve);

        let he_a = brep.edges[edge_handle].half_edge;
        let he_b = brep.half_edges[he_a].twin;
        let (he_forward, he_backward) = if brep.half_edges[he_a].origin == vertices[k] {
            (he_a, he_b)
        } else {
            (he_b, he_a)
        };
        let pcurve = curve_convert::convert_to_pcurve(directed, ctx.sketch, ctx.plane)?;
        let mut back = pcurve.clone();
        std::mem::swap(&mut back.t_start, &mut back.t_end);
        brep.half_edges[he_forward].pcurve = Some(pcurve);
        brep.half_edges[he_backward].pcurve = Some(back);
    }
    Ok(())
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
        let pcurve = curve_convert::convert_to_pcurve(&traced.edges[i], ctx.sketch, ctx.plane)?;
        brep.half_edges[he_fwd].pcurve = Some(pcurve);

        // The twin runs j → i: the same chart span, reversed
        // (t-swap — the canonical-direction convention, PK-303).
        let mut back_pcurve =
            curve_convert::convert_to_pcurve(&traced.edges[i], ctx.sketch, ctx.plane)?;
        std::mem::swap(&mut back_pcurve.t_start, &mut back_pcurve.t_end);
        brep.half_edges[he_rev].pcurve = Some(back_pcurve);

        front_hes.push(he_fwd);
        back_hes.push(he_rev);
    }

    // Back face uses reversed twin half-edges.
    back_hes.reverse();

    let front_loop = wire_inner_loop(brep, &front_hes, front_face)?;
    let back_loop = wire_inner_loop(brep, &back_hes, back_face)?;

    // A hole keeps the face's material on its left, so it winds
    // clockwise about the front face's normal. The tracer's direction
    // follows how the profile was drawn: a hole drawn counter-clockwise
    // wound like the outer loop on both faces (PK-923: the lamina with
    // a hole read signed UV areas +1 under flag true and -1 under
    // false). Such a hole trades its two sides between the faces.
    if loop_winding_about(brep, front_loop, ctx.plane.normal) > 0.0 {
        brep.loops[front_loop].face = back_face;
        brep.loops[back_loop].face = front_face;
        for (face, from, to) in [
            (front_face, front_loop, back_loop),
            (back_face, back_loop, front_loop),
        ] {
            for l in brep.faces[face].inner_loops.iter_mut() {
                if *l == from {
                    *l = to;
                }
            }
        }
        return Ok(Some(back_loop));
    }
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
