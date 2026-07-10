//! Side face creation from swept/lofted surfaces.
//!
//! Builds a BRep solid from the NURBS surfaces and cap boundary curves
//! produced by sweep/loft surface generation.

use popkern_core::constructors::builder::BRepBuilder;
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

    for (i, surface) in surfaces.iter().enumerate() {
        // Surface layout: u = profile, v = path.
        // Vertex at (u=0, v=0) = path start, profile start of this span.
        let start_pt = eval_surface_corner(surface, 0.0, 0.0);
        start_vertices.push(builder.add_vertex(start_pt));

        let end_pt = eval_surface_corner(surface, 0.0, 1.0);
        end_vertices.push(builder.add_vertex(end_pt));
        if std::env::var_os("PK_SWEEP_TRACE").is_some() {
            let u1 = eval_surface_corner(surface, 1.0, 0.0);
            eprintln!(
                "  [sweep] span {i}: u0 ({:.2},{:.2},{:.2}) u1 ({:.2},{:.2},{:.2})",
                start_pt.x, start_pt.y, start_pt.z, u1.x, u1.y, u1.z
            );
        }
    }

    // Create side faces (one per span).
    // Each side face has 4 edges forming a loop:
    //   bottom (cap-start), right (longitudinal next), top (cap-end reversed), left (longitudinal this)
    // For a single span, the "left" and "right" longitudinal edges are the same seam edge.

    // Longitudinal edges: ONE per profile vertex, shared between the
    // two adjacent side faces (span i's right edge = span i+1's left
    // edge). Creating them per-span duplicated every seam — a
    // rectangle extrusion had E=16 instead of 12 and failed
    // Euler–Poincaré at −2.
    let mut long_up: Vec<HalfEdgeHandle> = Vec::with_capacity(n_spans);
    let mut long_down: Vec<HalfEdgeHandle> = Vec::with_capacity(n_spans);
    for j in 0..n_spans {
        let (edge, he_up, he_down) = builder.add_edge(start_vertices[j], end_vertices[j]);
        // Longitudinal edges carry the wall surface's exact u=0
        // boundary iso-curve (v = path direction). Leaving them
        // curveless starved every curve consumer downstream — the
        // boundary-first tessellator has no canonical samples for
        // them, so bare (un-booleaned, un-filleted) sweep solids
        // could not mesh watertight.
        let start = builder.brep().vertices[start_vertices[j]].position;
        let end = builder.brep().vertices[end_vertices[j]].position;
        if let Some(iso) = path_boundary_iso_curve(&surfaces[j], start, end) {
            builder.set_edge_curve(
                edge,
                crate::surface::canonical::analytic_or_nurbs_curve(&iso),
            );
        }
        long_up.push(he_up);
        long_down.push(he_down);
    }

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

        // Longitudinal half-edges from the shared per-vertex edges:
        // span i walks UP its right seam (vertex `next`) and DOWN its
        // left seam (vertex `i`). Each half-edge is consumed by
        // exactly one face loop; the twins belong to the neighbors.
        let he_left_up = long_up[i];
        let he_left_down = long_down[i];
        let he_right_up = long_up[next];

        // Create the side face with its outer loop. Canonical
        // recognition: an extruded circular arc is stored as an
        // analytic Cylinder (exact boolean/fillet/export paths);
        // the pcurves below are then authored in the cylinder's
        // native (angle, length) chart via the chart map.
        let face = builder.add_face(shell);
        let (geometry, chart) = crate::surface::canonical::recognize_extruded_plane(&surfaces[i])
            .or_else(|| crate::surface::canonical::recognize_extruded_cylinder(&surfaces[i]))
            .map(|(g, c)| (g, Some(c)))
            .unwrap_or_else(|| (SurfaceGeometry::Nurbs(surfaces[i].clone()), None));
        // The loop below is authored against the NURBS frame (u =
        // profile, v = path), whose du×dv normal is the face's
        // outward direction. A recognized canonical surface carries
        // its OWN normal convention — a cylinder's normal_at is
        // +radial regardless of which side the material lies on
        // (concave profile arcs put material outside the circle), so
        // the orientation flag must record whether the canonical
        // normal agrees with the authored frame.
        let orientation =
            crate::surface::canonical::canonical_normal_agreement(&surfaces[i], &geometry);
        builder.set_face_surface(
            face,
            SurfaceOnFace {
                geometry,
                orientation,
            },
        );
        let uv_line_chart = |x0: f64, y0: f64, x1: f64, y1: f64| match &chart {
            Some(c) => {
                let (a0, b0) = c.map(x0, y0);
                let (a1, b1) = c.map(x1, y1);
                uv_line(a0, b0, a1, b1)
            }
            None => uv_line(x0, y0, x1, y1),
        };

        // Wire the loop: bottom -> right_up -> top_rev -> left_down
        // bottom: sv_this -> sv_next  (he_bottom_fwd)
        // right:  sv_next -> ev_next  (he_right_up)
        // top:    ev_next -> ev_this  (he_top_rev)
        // left:   ev_this -> sv_this  (he_left_down)
        let loop_hes = if n_spans == 1 {
            // Single span: left and right are the same seam edge
            // (sv_this == sv_next, ev_this == ev_next).
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
        builder.set_half_edge_pcurve(he_bottom_fwd, uv_line_chart(0.0, 0.0, 1.0, 0.0));
        builder.set_half_edge_pcurve(he_top_rev, uv_line_chart(1.0, 1.0, 0.0, 1.0));
        builder.set_half_edge_pcurve(he_left_down, uv_line_chart(0.0, 1.0, 0.0, 0.0));
        if n_spans == 1 {
            // Single span: left_up is the seam at u=0.
            builder.set_half_edge_pcurve(he_left_up, uv_line_chart(0.0, 0.0, 0.0, 1.0));
        } else {
            // Multi-span: this face's upward half-edge on its RIGHT
            // seam sits at u=1 of THIS surface (the twin — the
            // neighbor's left_down — gets its own pcurve from the
            // neighbor's iteration at u=0).
            builder.set_half_edge_pcurve(he_right_up, uv_line_chart(1.0, 0.0, 1.0, 1.0));
        }

        // Collect twin half-edges for cap loops.
        start_cap_hes.push(he_bottom_rev);
        end_cap_hes.push(he_top_fwd);
    }

    // Create start cap face. Its boundary half-edges run
    // sv_next → sv_this, so span order chains only when REVERSED
    // (the start cap winds opposite to the profile; single-span
    // circles never exposed this).
    let start_hes_rev: Vec<HalfEdgeHandle> = start_cap_hes.iter().rev().copied().collect();
    let start_curves_rev: Vec<NurbsCurve> = start_cap_curves.iter().rev().cloned().collect();
    create_cap_face(&mut builder, shell, &start_hes_rev, &start_curves_rev, true)?;

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
/// The path-side boundary iso-curve of a sweep wall surface at the
/// profile start (the longitudinal edge's exact geometry).
///
/// Extraction is VERIFICATION-DRIVEN rather than layout-trusting:
/// both grid interpretations (column stride `num_u` with v knots,
/// and leading row with u knots) are formed and the one whose
/// endpoints coincide with `(start, end)` wins. Grid-layout
/// assumptions have drifted before; endpoints cannot lie.
fn path_boundary_iso_curve(
    surface: &NurbsSurface,
    start: popkern_core::geometry::Point3,
    end: popkern_core::geometry::Point3,
) -> Option<NurbsCurve> {
    let strided = |stride: usize, count: usize, degree: usize, knots: &[f64]| {
        let mut control_points = Vec::with_capacity(count);
        let mut weights = Vec::with_capacity(count);
        for i in 0..count {
            let idx = i * stride;
            if idx >= surface.control_points.len() {
                return None;
            }
            control_points.push(surface.control_points[idx]);
            weights.push(surface.weights[idx]);
        }
        Some(NurbsCurve {
            degree,
            knots: knots.to_vec(),
            control_points,
            weights,
        })
    };
    let candidates = [
        strided(
            surface.num_u,
            surface.num_v,
            surface.degree_v,
            &surface.knots_v,
        ),
        strided(1, surface.num_u, surface.degree_u, &surface.knots_u),
        strided(
            surface.num_v,
            surface.num_u,
            surface.degree_u,
            &surface.knots_u,
        ),
        strided(1, surface.num_v, surface.degree_v, &surface.knots_v),
    ];
    for cand in candidates.into_iter().flatten() {
        let (Some(first), Some(last)) = (
            cand.control_points.first().copied(),
            cand.control_points.last().copied(),
        ) else {
            continue;
        };
        let forward = first.distance_to(&start) < 1e-9 && last.distance_to(&end) < 1e-9;
        let reversed = first.distance_to(&end) < 1e-9 && last.distance_to(&start) < 1e-9;
        if forward {
            return Some(cand);
        }
        if reversed {
            // Reverse to run start → end (edge direction).
            let mut c = cand;
            c.control_points.reverse();
            c.weights.reverse();
            let (lo, hi) = (c.knots[0], *c.knots.last()?);
            c.knots.reverse();
            for k in &mut c.knots {
                *k = lo + hi - *k;
            }
            return Some(c);
        }
    }
    None
}

fn attach_nurbs_curve_to_edge(
    builder: &mut BRepBuilder,
    edge: popkern_core::topology::EdgeHandle,
    curves: &[NurbsCurve],
    index: usize,
) {
    if let Some(curve) = curves.get(index) {
        builder.set_edge_curve(
            edge,
            crate::surface::canonical::analytic_or_nurbs_curve(curve),
        );
    }
}

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
