//! Topology of a PARTIAL revolution of a closed profile with an edge
//! on the axis: a sector of a cylinder, cone or sphere. The axis span
//! sweeps no surface; it stays one edge, shared by the start cap and
//! the end cap, which meet along it. Each profile vertex on the axis
//! sweeps a point: one pole vertex used at both angles, with no
//! longitudinal edge, so a side face touching a pole closes there
//! (a triangle at a cone's apex, a two-edge lune on a sphere).
//!
//! Geometry follows the full revolution's rule (PK-409): when every
//! span is recognised as a canonical surface of revolution the faces
//! carry the exact surface with pcurves in its native chart and the
//! arcs are angle-parameterised; otherwise every face stays NURBS.
//!
//! Counts for an n-span off-axis chain (n + 1 profile vertices, two of
//! them poles): V = 2 + 2(n − 1), E = 1 + (n − 1) + 2n, F = n + 2, so
//! χ = 2.

use popkern_core::constructors::builder::BRepBuilder;
use popkern_core::geometry::surface::{SurfaceGeometry, SurfaceOnFace};
use popkern_core::geometry::{NurbsCurve, NurbsSurface, Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, HalfEdgeHandle};

use crate::error::SweepError;
use crate::surface::canonical::canonical_normal_agreement;
use crate::topology::caps::create_cap_face;
use crate::topology::revolution_canonical::{RevolvedSurface, recognize_revolved_span};
use crate::topology::revolution_faces::{canonical_pcurves, sweep_sense, vertex_arc};
use crate::topology::side_faces::{
    attach_nurbs_curve_to_edge, eval_surface_corner, path_boundary_iso_curve, uv_line,
};

/// A sector to build.
///
/// `spans` run from one pole to the other off the axis (span `i` ends
/// where span `i + 1` starts); `axis_span` closes the profile along
/// the axis from the last pole back to the first. `surfaces[i]` is the
/// revolution of span `i` (u = profile, v = angle normalised to
/// [0, 1]); `end_spans[i]` is span `i` turned through the full angle,
/// the axis span unmoved by the turn closing both caps.
pub(crate) struct SectorInput<'a> {
    /// Off-axis spans, pole to pole.
    pub spans: &'a [NurbsCurve],
    /// The axis run from the last pole back to the first.
    pub axis_span: &'a NurbsCurve,
    /// One revolution surface per span.
    pub surfaces: &'a [NurbsSurface],
    /// The spans turned through `angle`.
    pub end_spans: &'a [NurbsCurve],
    /// A point on the axis.
    pub axis_origin: Point3,
    /// The axis direction.
    pub axis_dir: Vec3,
    /// The sector's angle (radians, below a whole turn).
    pub angle: f64,
}

/// Build the topology of a sector (see [`SectorInput`]).
pub(crate) fn build_sector_topology(
    input: &SectorInput<'_>,
    tolerance: Tolerance,
) -> Result<BRep, SweepError> {
    let SectorInput {
        spans,
        axis_span,
        surfaces,
        end_spans,
        axis_origin,
        axis_dir,
        angle,
    } = *input;
    let n = surfaces.len();
    if n == 0 || spans.len() != n || end_spans.len() != n {
        return Err(SweepError::Topology(
            "a sector needs one surface and one end curve per off-axis span".into(),
        ));
    }
    let axis_dir = axis_dir
        .normalized()
        .ok_or_else(|| SweepError::InvalidProfile("zero-length revolution axis".into()))?;
    let recognized: Option<Vec<RevolvedSurface>> = spans
        .iter()
        .zip(surfaces)
        .map(|(span, surface)| {
            recognize_revolved_span(
                span,
                surface,
                &axis_origin,
                &axis_dir,
                angle,
                tolerance.linear(),
            )
        })
        .collect();
    let sense = sweep_sense(&surfaces[0], &axis_origin, &axis_dir);
    let mut builder = BRepBuilder::new(tolerance);
    let solid = builder.add_solid();
    let shell = builder.add_shell(solid);

    // Profile vertex j at the start and end angles; the poles (j = 0
    // and j = n) are one vertex each.
    let first_pole = builder.add_vertex(eval_surface_corner(&surfaces[0], 0.0, 0.0));
    let last_pole = builder.add_vertex(eval_surface_corner(&surfaces[n - 1], 1.0, 0.0));
    let mut start_vertices = vec![first_pole];
    let mut end_vertices = vec![first_pole];
    for surface in &surfaces[1..] {
        start_vertices.push(builder.add_vertex(eval_surface_corner(surface, 0.0, 0.0)));
        end_vertices.push(builder.add_vertex(eval_surface_corner(surface, 0.0, 1.0)));
    }
    start_vertices.push(last_pole);
    end_vertices.push(last_pole);

    // The arc each interior profile vertex sweeps, shared by the two
    // side faces on either side of it.
    let mut long_up: Vec<Option<HalfEdgeHandle>> = vec![None; n + 1];
    let mut long_down: Vec<Option<HalfEdgeHandle>> = vec![None; n + 1];
    for j in 1..n {
        let (edge, up, down) = builder.add_edge(start_vertices[j], end_vertices[j]);
        let start = builder.brep().vertices[start_vertices[j]].position;
        let end = builder.brep().vertices[end_vertices[j]].position;
        if recognized.is_some() {
            builder.set_edge_curve(
                edge,
                vertex_arc(&start, &axis_origin, &axis_dir, sense, angle),
            );
        } else if let Some(iso) = path_boundary_iso_curve(&surfaces[j], start, end) {
            builder.set_edge_curve(
                edge,
                crate::surface::canonical::analytic_or_nurbs_curve(&iso),
            );
        }
        long_up[j] = Some(up);
        long_down[j] = Some(down);
    }

    let mut start_cap_hes = Vec::with_capacity(n + 1);
    let mut end_cap_hes = Vec::with_capacity(n + 1);
    for i in 0..n {
        let (bottom, bottom_fwd, bottom_rev) =
            builder.add_edge(start_vertices[i], start_vertices[i + 1]);
        attach_nurbs_curve_to_edge(&mut builder, bottom, spans, i);
        let (top, top_fwd, top_rev) = builder.add_edge(end_vertices[i], end_vertices[i + 1]);
        attach_nurbs_curve_to_edge(&mut builder, top, end_spans, i);

        let face = builder.add_face(shell);
        let (surface_on_face, (bottom_pc, up_pc, top_pc, down_pc)) = match &recognized {
            Some(exact) => (
                SurfaceOnFace {
                    orientation: canonical_normal_agreement(&surfaces[i], &exact[i].geometry),
                    geometry: exact[i].geometry.clone(),
                },
                canonical_pcurves(&exact[i].chart),
            ),
            // The face chart is the transpose of the construction
            // surface, angle in u and profile in v, as a full
            // revolution's (PK-456): a pole is then a collapsed v
            // boundary, where the meshers and the measuring quadrature
            // expect one (in the construction chart a fan measured
            // 0.3% over).
            None => (
                SurfaceOnFace {
                    geometry: SurfaceGeometry::Nurbs(surfaces[i].transposed()),
                    orientation: false,
                },
                (
                    uv_line(0.0, 0.0, 0.0, 1.0),
                    uv_line(0.0, 1.0, 1.0, 1.0),
                    uv_line(1.0, 1.0, 1.0, 0.0),
                    uv_line(1.0, 0.0, 0.0, 0.0),
                ),
            ),
        };
        builder.set_face_surface(face, surface_on_face);
        // bottom → right arc up → top back → left arc down; an arc at
        // a pole is absent and the loop closes on the pole vertex.
        let mut loop_hes = vec![bottom_fwd];
        builder.set_half_edge_pcurve(bottom_fwd, bottom_pc);
        if let Some(up) = long_up[i + 1] {
            loop_hes.push(up);
            builder.set_half_edge_pcurve(up, up_pc);
        }
        loop_hes.push(top_rev);
        builder.set_half_edge_pcurve(top_rev, top_pc);
        if let Some(down) = long_down[i] {
            loop_hes.push(down);
            builder.set_half_edge_pcurve(down, down_pc);
        }
        builder
            .add_loop(&loop_hes, face)
            .map_err(|e| SweepError::Topology(format!("sector side face loop: {e}")))?;
        start_cap_hes.push(bottom_rev);
        end_cap_hes.push(top_fwd);
    }

    // The axis edge runs from the last pole to the first, as the
    // profile does; the end cap walks it that way, the start cap back.
    let (_, axis_fwd, axis_rev) = builder.add_edge(last_pole, first_pole);
    start_cap_hes.push(axis_rev);
    end_cap_hes.push(axis_fwd);
    let mut start_curves: Vec<NurbsCurve> = spans.to_vec();
    start_curves.push(axis_span.clone());
    let mut end_curves: Vec<NurbsCurve> = end_spans.to_vec();
    end_curves.push(axis_span.clone());

    // The start cap winds against the profile (see the closed sweep
    // builder): its half-edges chain in reverse span order.
    let start_hes_rev: Vec<HalfEdgeHandle> = start_cap_hes.iter().rev().copied().collect();
    let start_curves_rev: Vec<NurbsCurve> = start_curves.iter().rev().cloned().collect();
    create_cap_face(&mut builder, shell, &start_hes_rev, &start_curves_rev, true)?;
    create_cap_face(&mut builder, shell, &end_cap_hes, &end_curves, false)?;
    Ok(builder.build_unchecked())
}
