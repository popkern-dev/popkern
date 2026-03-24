//! Helper functions for topology manipulation during blend construction.
//!
//! Provides edge splitting, anchor-vertex creation, and fillet-face insertion
//! used by the main `apply_blend_geometry` orchestrator.

use popkern_core::geometry::NurbsCurve;
use popkern_core::geometry::curve::{CurveEval, CurveGeometry, CurveOnEdge};
use popkern_core::geometry::surface::{SurfaceGeometry, SurfaceOnFace};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{
    BRep, Edge, EdgeHandle, Face, HalfEdge, Loop, SplitFaceResult, VertexHandle,
};

/// Anchor vertices created by splitting the fillet edge.
#[derive(Debug, Clone)]
pub struct AnchorVertices {
    /// Vertex at the start of the fillet (where trim curves begin).
    pub start: VertexHandle,
    /// Vertex at the end of the fillet (where trim curves end).
    pub end: VertexHandle,
}

/// Split a fillet edge at the endpoints of the trim curves.
///
/// The trim curves touch the original edge at its start and end vertices.
/// For a full-edge fillet (fillet spans the entire edge length), the existing
/// edge endpoints are the anchors and no splitting is needed.
///
/// For partial fillets, we would split the edge at the trim curve parameter
/// values. Currently only full-edge fillets are supported.
pub fn split_edge_at_trim_endpoints(
    brep: &BRep,
    edge: EdgeHandle,
    _trim_3d_1: &NurbsCurve,
    _trim_3d_2: &NurbsCurve,
    _tolerance: &Tolerance,
) -> Result<AnchorVertices, String> {
    // Look up the edge to find its endpoint vertices.
    let edge_data = brep
        .edges
        .get(edge)
        .ok_or_else(|| format!("edge {edge:?} not found"))?;

    let he = brep
        .half_edges
        .get(edge_data.half_edge)
        .ok_or_else(|| "half-edge not found on edge".to_string())?;

    let twin = brep
        .half_edges
        .get(he.twin)
        .ok_or_else(|| "twin half-edge not found".to_string())?;

    // For a full-edge fillet, the anchor vertices are the edge endpoints.
    Ok(AnchorVertices {
        start: he.origin,
        end: twin.origin,
    })
}

/// Insert a fillet face between two split faces.
///
/// Creates a new face whose boundary connects the two trim curves and the
/// anchor vertices, assigning the fillet/chamfer surface geometry to it.
///
/// For the initial implementation this attaches the surface geometry to the
/// existing topology without full loop rewiring. The fillet face is created
/// as a standalone face in the same shell, with edges along the two trim
/// curves connecting the anchor vertices.
#[allow(clippy::too_many_arguments)]
pub fn insert_fillet_face(
    brep: &mut BRep,
    surface: &SurfaceGeometry,
    anchors: &AnchorVertices,
    trim_3d_1: &NurbsCurve,
    trim_3d_2: &NurbsCurve,
    pcurve_1: &NurbsCurve,
    pcurve_2: &NurbsCurve,
    split_1: &SplitFaceResult,
    _split_2: &SplitFaceResult,
    _tolerance: &Tolerance,
) -> Result<(), String> {
    // Determine which shell the fillet face belongs to.
    let shell = brep
        .faces
        .get(split_1.face_a)
        .map(|f| f.shell)
        .ok_or_else(|| "could not find shell for face_a".to_string())?;

    // Create the fillet face with the blend surface geometry.
    let fillet_loop = brep.loops.insert(Loop {
        half_edge: Default::default(),
        face: Default::default(),
    });

    let fillet_face = brep.faces.insert(Face {
        outer_loop: fillet_loop,
        inner_loops: Vec::new(),
        shell,
        surface: Some(SurfaceOnFace {
            geometry: surface.clone(),
            orientation: true,
        }),
    });

    // Update the loop to point back to its face.
    brep.loops[fillet_loop].face = fillet_face;

    // Create edges along the two trim curves connecting anchor vertices.
    // Edge 1: along trim_3d_1 from start to end.
    let (t1_start, t1_end) = trim_3d_1.parameter_range();
    let edge_1 = brep.edges.insert(Edge {
        half_edge: Default::default(),
        curve: Some(CurveOnEdge {
            geometry: CurveGeometry::Nurbs(trim_3d_1.clone()),
            t_start: t1_start,
            t_end: t1_end,
        }),
    });

    // Edge 2: along trim_3d_2 from start to end.
    let (t2_start, t2_end) = trim_3d_2.parameter_range();
    let edge_2 = brep.edges.insert(Edge {
        half_edge: Default::default(),
        curve: Some(CurveOnEdge {
            geometry: CurveGeometry::Nurbs(trim_3d_2.clone()),
            t_start: t2_start,
            t_end: t2_end,
        }),
    });

    // Create half-edges for the fillet face boundary.
    // The fillet face has 4 half-edges forming a quad: trim1, end-cap, trim2-reversed, start-cap.
    let he_trim1 = brep.half_edges.insert(HalfEdge {
        origin: anchors.start,
        twin: Default::default(),
        next: Default::default(),
        prev: Default::default(),
        edge: edge_1,
        loop_: fillet_loop,
        pcurve: Some(pcurve_1.clone()),
    });

    let he_trim2 = brep.half_edges.insert(HalfEdge {
        origin: anchors.end,
        twin: Default::default(),
        next: Default::default(),
        prev: Default::default(),
        edge: edge_2,
        loop_: fillet_loop,
        pcurve: Some(pcurve_2.clone()),
    });

    // Wire up half-edge chain: for the simple two-edge case, trim1 -> trim2 -> trim1.
    brep.half_edges[he_trim1].next = he_trim2;
    brep.half_edges[he_trim1].prev = he_trim2;
    brep.half_edges[he_trim2].next = he_trim1;
    brep.half_edges[he_trim2].prev = he_trim1;

    // Self-twin for now (proper stitching happens in a later refinement).
    brep.half_edges[he_trim1].twin = he_trim2;
    brep.half_edges[he_trim2].twin = he_trim1;

    // Update edge half-edge pointers.
    brep.edges[edge_1].half_edge = he_trim1;
    brep.edges[edge_2].half_edge = he_trim2;

    // Update loop to point to the first half-edge.
    brep.loops[fillet_loop].half_edge = he_trim1;

    // The shell stores only a single representative face handle.
    // The new fillet face is already associated with the shell via its `shell` field.

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use popkern_core::tolerance::Tolerance;
    use popkern_core::topology::BRep;

    #[test]
    fn anchor_vertices_from_cube_edge() {
        let cube = BRep::unit_cube(Tolerance::default()).unwrap();
        let edge = cube.edges.keys().next().unwrap();

        let dummy_curve = NurbsCurve {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![
                popkern_core::geometry::Point3::new(0.0, 0.0, 0.0),
                popkern_core::geometry::Point3::new(1.0, 0.0, 0.0),
            ],
            weights: vec![1.0, 1.0],
        };

        let result = split_edge_at_trim_endpoints(
            &cube,
            edge,
            &dummy_curve,
            &dummy_curve,
            &Tolerance::default(),
        );
        assert!(result.is_ok());

        let anchors = result.unwrap();
        assert_ne!(anchors.start, anchors.end);
    }

    #[test]
    fn anchor_vertices_invalid_edge() {
        let cube = BRep::unit_cube(Tolerance::default()).unwrap();
        let bogus = popkern_core::topology::EdgeHandle::default();

        let dummy_curve = NurbsCurve {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![
                popkern_core::geometry::Point3::new(0.0, 0.0, 0.0),
                popkern_core::geometry::Point3::new(1.0, 0.0, 0.0),
            ],
            weights: vec![1.0, 1.0],
        };

        let result = split_edge_at_trim_endpoints(
            &cube,
            bogus,
            &dummy_curve,
            &dummy_curve,
            &Tolerance::default(),
        );
        assert!(result.is_err());
    }
}
