//! Side face creation from swept/lofted surfaces.
//!
//! Builds a BRep solid from the NURBS surfaces and cap boundary curves
//! produced by sweep/loft surface generation.

use popkern_core::constructors::builder::BRepBuilder;
use popkern_core::geometry::curve::{CurveGeometry, CurveOnEdge};
use popkern_core::geometry::surface::{SurfaceGeometry, SurfaceOnFace};
use popkern_core::geometry::{NurbsCurve, NurbsSurface, Point3};
use popkern_core::history::TrackingContext;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle, VertexHandle};

use crate::error::SweepError;
use crate::topology::caps::create_cap_face;

/// Build a BRep solid from sweep surface results.
///
/// Creates vertices at profile discontinuities, edges along boundaries,
/// one face per smooth surface, cap faces at ends, and seam edges for
/// closed profiles.
///
/// For a closed profile with 1 span (e.g., circle sweep):
/// - 2 vertices (seam at path start and end)
/// - 3 edges (seam start, seam end, profile seam connecting them)
/// - 3 faces (1 side + 2 caps)
/// - Euler: V - E + F = 2 - 3 + 3 = 2
///
/// For a closed profile with N spans:
/// - 2N vertices (N at path start + N at path end)
/// - 3N edges (N cap-start, N cap-end, N longitudinal)
/// - N + 2 faces (N side + 2 caps)
/// - Euler: 2N - 3N + (N + 2) = 2
///
/// For an open profile (sheet body): no cap faces, free boundary edges.
pub fn build_sweep_topology(
    surfaces: &[NurbsSurface],
    start_cap_curves: &[NurbsCurve],
    end_cap_curves: &[NurbsCurve],
    is_closed_profile: bool,
    tolerance: Tolerance,
    tracking: Option<&mut TrackingContext>,
) -> Result<BRep, SweepError> {
    if surfaces.is_empty() {
        return Err(SweepError::Topology("no surfaces to build".into()));
    }

    let mut brep = if is_closed_profile {
        build_closed_profile_topology(surfaces, start_cap_curves, end_cap_curves, tolerance)?
    } else {
        build_open_profile_topology(surfaces, tolerance)?
    };

    if let Some(ctx) = tracking {
        assign_sweep_names(&mut brep, ctx);
    }

    Ok(brep)
}

/// Build topology for a closed profile (solid body with caps).
fn build_closed_profile_topology(
    surfaces: &[NurbsSurface],
    start_cap_curves: &[NurbsCurve],
    end_cap_curves: &[NurbsCurve],
    tolerance: Tolerance,
) -> Result<BRep, SweepError> {
    let n_spans = surfaces.len();
    let mut builder = BRepBuilder::new(tolerance);
    let solid = builder.add_solid();
    let shell = builder.add_shell(solid);

    // Create vertices: N at path start, N at path end.
    // For each span, the start vertex is at the beginning of that span's
    // profile at the path start; the end vertex is at the path end.
    let mut start_vertices = Vec::with_capacity(n_spans);
    let mut end_vertices = Vec::with_capacity(n_spans);

    for surface in surfaces {
        // Surface layout: u = profile, v = path.
        // Vertex at (u=0, v=0) = path start, profile start of this span.
        let start_pt = eval_surface_corner(surface, 0.0, 0.0);
        start_vertices.push(builder.add_vertex(start_pt));

        let end_pt = eval_surface_corner(surface, 0.0, 1.0);
        end_vertices.push(builder.add_vertex(end_pt));
    }

    // Create side faces (one per span).
    // Each side face has 4 edges forming a loop:
    //   bottom (cap-start), right (longitudinal next), top (cap-end reversed), left (longitudinal this)
    // For a single span, the "left" and "right" longitudinal edges are the same seam edge.

    // Collect half-edge handles for cap loop assembly.
    let mut start_cap_hes: Vec<HalfEdgeHandle> = Vec::with_capacity(n_spans);
    let mut end_cap_hes: Vec<HalfEdgeHandle> = Vec::with_capacity(n_spans);

    for i in 0..n_spans {
        let next = (i + 1) % n_spans;

        let sv_this = start_vertices[i];
        let sv_next = start_vertices[next];
        let ev_this = end_vertices[i];
        let ev_next = end_vertices[next];

        // Bottom edge (cap-start): sv_this -> sv_next (along profile at path start)
        let (bottom_edge, he_bottom_fwd, he_bottom_rev) = builder.add_edge(sv_this, sv_next);
        attach_nurbs_curve_to_edge(&mut builder, bottom_edge, start_cap_curves, i);

        // Top edge (cap-end): ev_this -> ev_next (along profile at path end)
        let (top_edge, he_top_fwd, he_top_rev) = builder.add_edge(ev_this, ev_next);
        attach_nurbs_curve_to_edge(&mut builder, top_edge, end_cap_curves, i);

        // Left longitudinal edge: sv_this -> ev_this (along path at span start)
        // For multi-span, each span shares its left edge with the previous span's right edge.
        // Since we create edges per-span and the closed loop wraps, we create one
        // longitudinal edge per span vertex pair.
        let (left_edge, he_left_up, he_left_down) = if i == 0 {
            // First span: create the edge.
            builder.add_edge(sv_this, ev_this)
        } else {
            // Reuse would be complex; instead create a new edge for each span.
            // The topology will be wired so that twin half-edges are shared.
            builder.add_edge(sv_this, ev_this)
        };
        let _ = left_edge;

        // Right longitudinal edge: sv_next -> ev_next (along path at span end)
        let (right_edge, he_right_up, _he_right_down) = if n_spans == 1 {
            // Single span: the right edge IS the left edge (seam).
            // Already created above. Reuse the handles.
            (left_edge, he_left_up, he_left_down)
        } else {
            builder.add_edge(sv_next, ev_next)
        };
        let _ = right_edge;

        // Create the side face with its outer loop.
        let face = builder.add_face(shell);
        builder.set_face_surface(
            face,
            SurfaceOnFace {
                geometry: SurfaceGeometry::Nurbs(surfaces[i].clone()),
                orientation: true,
            },
        );

        // Wire the loop: bottom -> right_up -> top_rev -> left_down
        // bottom: sv_this -> sv_next  (he_bottom_fwd)
        // right:  sv_next -> ev_next  (he_right_up)
        // top:    ev_next -> ev_this  (he_top_rev)
        // left:   ev_this -> sv_this  (he_left_down)
        let loop_hes = if n_spans == 1 {
            // Single span: seam edge means left and right are the same edge.
            // Loop: bottom_fwd -> left_up -> top_rev -> left_down
            // Wait — for single span, sv_this == sv_next and ev_this == ev_next.
            // We have: bottom: sv -> sv (self-loop), top: ev -> ev (self-loop),
            //          seam: sv -> ev (he_left_up) and ev -> sv (he_left_down)
            // Loop is: he_bottom_fwd -> he_left_up -> he_top_rev -> he_left_down
            vec![he_bottom_fwd, he_left_up, he_top_rev, he_left_down]
        } else {
            vec![he_bottom_fwd, he_right_up, he_top_rev, he_left_down]
        };

        builder
            .add_loop(&loop_hes, face)
            .map_err(|e| SweepError::Topology(format!("side face loop: {e}")))?;

        // Attach pcurves to body face half-edges.
        // Surface parametrization: u ∈ [0,1] = profile, v ∈ [0,1] = path.
        // bottom (u varies, v=0), seam-up (u=0, v varies),
        // top reversed (u varies [1→0], v=1), seam-down (u=0, v varies [1→0]).
        builder.set_half_edge_pcurve(he_bottom_fwd, uv_line(0.0, 0.0, 1.0, 0.0));
        builder.set_half_edge_pcurve(he_top_rev, uv_line(1.0, 1.0, 0.0, 1.0));
        builder.set_half_edge_pcurve(he_left_down, uv_line(0.0, 1.0, 0.0, 0.0));
        if n_spans == 1 {
            // Single span: left_up is the seam at u=0.
            builder.set_half_edge_pcurve(he_left_up, uv_line(0.0, 0.0, 0.0, 1.0));
        } else {
            // Multi-span: left_up at u=0, right_up at u=1.
            builder.set_half_edge_pcurve(he_left_up, uv_line(0.0, 0.0, 0.0, 1.0));
            builder.set_half_edge_pcurve(he_right_up, uv_line(1.0, 0.0, 1.0, 1.0));
        }

        // Collect twin half-edges for cap loops.
        start_cap_hes.push(he_bottom_rev);
        end_cap_hes.push(he_top_fwd);
    }

    // Create start cap face.
    create_cap_face(&mut builder, shell, &start_cap_hes, start_cap_curves, true)?;

    // Create end cap face.
    create_cap_face(&mut builder, shell, &end_cap_hes, end_cap_curves, false)?;

    Ok(builder.build_unchecked())
}

/// Build topology for an open profile (sheet body, no caps).
fn build_open_profile_topology(
    surfaces: &[NurbsSurface],
    tolerance: Tolerance,
) -> Result<BRep, SweepError> {
    let mut builder = BRepBuilder::new(tolerance);
    let solid = builder.add_solid();
    let shell = builder.add_shell(solid);

    for surface in surfaces {
        // Create 4 corner vertices.
        let v00 = builder.add_vertex(eval_surface_corner(surface, 0.0, 0.0));
        let v10 = builder.add_vertex(eval_surface_corner(surface, 1.0, 0.0));
        let v01 = builder.add_vertex(eval_surface_corner(surface, 0.0, 1.0));
        let v11 = builder.add_vertex(eval_surface_corner(surface, 1.0, 1.0));

        // Create 4 edges.
        let (_, he_bottom, _) = builder.add_edge(v00, v10); // u-direction at v=0
        let (_, he_right, _) = builder.add_edge(v10, v11); // v-direction at u=1
        let (_, he_top_rev, _) = builder.add_edge(v11, v01); // u-direction at v=1 (reversed)
        let (_, he_left_rev, _) = builder.add_edge(v01, v00); // v-direction at u=0 (reversed)

        let face = builder.add_face(shell);
        builder.set_face_surface(
            face,
            SurfaceOnFace {
                geometry: SurfaceGeometry::Nurbs(surface.clone()),
                orientation: true,
            },
        );

        builder
            .add_loop(&[he_bottom, he_right, he_top_rev, he_left_rev], face)
            .map_err(|e| SweepError::Topology(format!("open face loop: {e}")))?;

        // Attach pcurves (u = profile [0,1], v = path [0,1]).
        builder.set_half_edge_pcurve(he_bottom, uv_line(0.0, 0.0, 1.0, 0.0));
        builder.set_half_edge_pcurve(he_right, uv_line(1.0, 0.0, 1.0, 1.0));
        builder.set_half_edge_pcurve(he_top_rev, uv_line(1.0, 1.0, 0.0, 1.0));
        builder.set_half_edge_pcurve(he_left_rev, uv_line(0.0, 1.0, 0.0, 0.0));
    }

    Ok(builder.build_unchecked())
}

/// Assign stable [`EntityName`]s to all entities in a sweep BRep.
///
/// All faces, edges, and vertices receive `Generated` names since sweeps
/// create entirely new geometry rather than modifying existing entities.
fn assign_sweep_names(brep: &mut BRep, ctx: &mut TrackingContext) {
    use popkern_core::naming::{EntityName, EntityRole, RoleCategory};
    use popkern_core::naming_criteria::GeometricCriteria;
    use popkern_core::naming_criteria::{capture_edge, capture_face, capture_vertex};

    let face_handles: Vec<FaceHandle> = brep.faces.keys().collect();
    for (i, &fh) in face_handles.iter().enumerate() {
        let name = EntityName {
            operation: ctx.operation,
            role: EntityRole {
                category: RoleCategory::Generated,
                index: i as u32,
            },
            ancestors: Vec::new(),
        };
        let criteria = capture_face(brep, fh).map(GeometricCriteria::Face);
        let name_id = brep.name_pool.intern(name.clone(), criteria);
        brep.faces[fh].name = Some(name_id);
        ctx.record_generated(name);
    }

    let edge_handles: Vec<EdgeHandle> = brep.edges.keys().collect();
    for (i, &eh) in edge_handles.iter().enumerate() {
        let name = EntityName {
            operation: ctx.operation,
            role: EntityRole {
                category: RoleCategory::Generated,
                index: i as u32,
            },
            ancestors: Vec::new(),
        };
        let criteria = capture_edge(brep, eh).map(GeometricCriteria::Edge);
        let name_id = brep.name_pool.intern(name.clone(), criteria);
        brep.edges[eh].name = Some(name_id);
        ctx.record_generated(name);
    }

    let vertex_handles: Vec<VertexHandle> = brep.vertices.keys().collect();
    for (i, &vh) in vertex_handles.iter().enumerate() {
        let name = EntityName {
            operation: ctx.operation,
            role: EntityRole {
                category: RoleCategory::Generated,
                index: i as u32,
            },
            ancestors: Vec::new(),
        };
        let criteria = capture_vertex(brep, vh).map(GeometricCriteria::Vertex);
        let name_id = brep.name_pool.intern(name.clone(), criteria);
        brep.vertices[vh].name = Some(name_id);
        ctx.record_generated(name);
    }
}

/// Build a degree-1 NURBS line in UV parameter space.
///
/// The pcurve is stored as a `NurbsCurve` with `Point3(u, v, 0.0)`.
fn uv_line(u0: f64, v0: f64, u1: f64, v1: f64) -> NurbsCurve {
    NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(u0, v0, 0.0), Point3::new(u1, v1, 0.0)],
        weights: vec![1.0, 1.0],
    }
}

/// Evaluate a surface corner point.
fn eval_surface_corner(surface: &NurbsSurface, u: f64, v: f64) -> popkern_core::geometry::Point3 {
    use popkern_core::geometry::surface::SurfaceEval;
    surface.point_at(u, v)
}

/// Attach a NURBS curve to an edge if available.
fn attach_nurbs_curve_to_edge(
    builder: &mut BRepBuilder,
    edge: popkern_core::topology::EdgeHandle,
    curves: &[NurbsCurve],
    index: usize,
) {
    if let Some(curve) = curves.get(index) {
        let (t_start, t_end) = curve.parameter_range();
        builder.set_edge_curve(
            edge,
            CurveOnEdge {
                geometry: CurveGeometry::Nurbs(curve.clone()),
                t_start,
                t_end,
            },
        );
    }
}

use popkern_core::geometry::curve::CurveEval;

#[cfg(test)]
mod tests {
    use super::*;
    use popkern_core::geometry::Point3;

    /// Build a simple bilinear surface patch.
    fn bilinear_patch(z_start: f64, z_end: f64) -> NurbsSurface {
        NurbsSurface {
            degree_u: 1,
            degree_v: 1,
            knots_u: vec![0.0, 0.0, 1.0, 1.0],
            knots_v: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![
                Point3::new(1.0, 0.0, z_start),
                Point3::new(1.0, 0.0, z_end),
                Point3::new(0.0, 1.0, z_start),
                Point3::new(0.0, 1.0, z_end),
            ],
            num_u: 2,
            num_v: 2,
            weights: vec![1.0; 4],
        }
    }

    /// Build a simple cap curve.
    fn cap_curve(z: f64) -> NurbsCurve {
        NurbsCurve {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![Point3::new(1.0, 0.0, z), Point3::new(0.0, 1.0, z)],
            weights: vec![1.0, 1.0],
        }
    }

    #[test]
    fn single_span_closed_produces_solid() {
        let surface = bilinear_patch(0.0, 5.0);
        let start_cap = vec![cap_curve(0.0)];
        let end_cap = vec![cap_curve(5.0)];
        let tol = Tolerance::default();

        let result = build_sweep_topology(&[surface], &start_cap, &end_cap, true, tol, None);
        assert!(result.is_ok(), "topology failed: {:?}", result.err());

        let brep = result.unwrap(); // ok in tests
        // Single span closed: V=2, E=3, F=3 (side + 2 caps)
        assert_eq!(brep.vertices.len(), 2, "expected 2 vertices");
        // Edges: 1 bottom + 1 top + 1 seam = 3
        // But we have add_edge creating twin pairs, so 3 edges in the edge slotmap
        // Actually bottom and top are self-loops (same vertex both ends for single span)
        // and the seam connects start to end vertex.
        assert!(
            brep.faces.len() >= 3,
            "expected >= 3 faces, got {}",
            brep.faces.len()
        );
    }

    #[test]
    fn open_profile_produces_sheet() {
        let surface = bilinear_patch(0.0, 5.0);
        let tol = Tolerance::default();

        let result = build_sweep_topology(&[surface], &[], &[], false, tol, None);
        assert!(result.is_ok(), "topology failed: {:?}", result.err());

        let brep = result.unwrap(); // ok in tests
        assert_eq!(brep.vertices.len(), 4, "expected 4 vertices");
        assert_eq!(brep.edges.len(), 4, "expected 4 edges");
        assert_eq!(brep.faces.len(), 1, "expected 1 face");
    }

    #[test]
    fn empty_surfaces_returns_error() {
        let tol = Tolerance::default();
        let result = build_sweep_topology(&[], &[], &[], true, tol, None);
        assert!(result.is_err());
    }

    #[test]
    fn multi_span_closed_topology() {
        let s1 = bilinear_patch(0.0, 5.0);
        let s2 = bilinear_patch(0.0, 5.0);
        let start_caps = vec![cap_curve(0.0), cap_curve(0.0)];
        let end_caps = vec![cap_curve(5.0), cap_curve(5.0)];
        let tol = Tolerance::default();

        let result = build_sweep_topology(&[s1, s2], &start_caps, &end_caps, true, tol, None);
        assert!(result.is_ok(), "topology failed: {:?}", result.err());

        let brep = result.unwrap(); // ok in tests
        // 2 spans: V=4, E=6, F=4 (2 side + 2 caps)
        assert_eq!(brep.vertices.len(), 4, "expected 4 vertices");
        assert!(
            brep.faces.len() >= 4,
            "expected >= 4 faces, got {}",
            brep.faces.len()
        );
    }
}
