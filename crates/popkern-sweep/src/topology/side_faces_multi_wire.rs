//! Closed multi-wire sweep topology (PK-477): one side-face ring per
//! wire, one cap per path end whose outer loop is wire 0 and whose
//! inner loops are the remaining wires.
//!
//! The single-wire builder in `side_faces.rs` rings ALL spans with
//! `next = (i + 1) % n_spans`; fed a flattened multi-wire span list it
//! cross-stitched the outer and inner wires' walls into one bogus ring
//! (measured: an outer+inner square tube swept Ok with 10 faces and
//! failed `validate()` with a 1.414 loop gap). This module is a
//! separate path that only activates for `wires.len() > 1`, so the
//! single-wire builder stays byte-identical.
//!
//! Invariants (stated per the topology-change contract):
//! - Winding: `analyze_multi_wire` normalizes the outer wire CCW and
//!   inner wires CW in the profile plane, so `du×dv` of every side
//!   surface — hole walls included — is the face's outward direction,
//!   and the per-wire ring wiring below is identical for both roles.
//! - Postconditions: the result must pass `validate()` (any-genus
//!   Euler–Poincaré: a W-wire tube body has genus W−1) — callers gate
//!   through the sweep pipeline's existing validation and fixtures.
//! - Rollback: pure construction; errors return before any partial
//!   body escapes.
//! - Diagnostics: wiring failures name the wire and stage in the
//!   `SweepError::Topology` message.
//! - Not fixture-specific: the ring code is span-count agnostic per
//!   wire (single-span self-seam wires included) and wire-count
//!   agnostic beyond the outer/inner role split.
//! - Tolerance: geometric tolerance is the builder's, passed
//!   explicitly; no additional epsilons are introduced.

use popkern_core::constructors::builder::BRepBuilder;
use popkern_core::geometry::surface::{SurfaceGeometry, SurfaceOnFace};
use popkern_core::geometry::{NurbsCurve, NurbsSurface};
use popkern_core::history::TrackingContext;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, HalfEdgeHandle};

use super::side_faces::{
    SideFaceGeometry, attach_nurbs_curve_to_edge, eval_surface_corner, path_boundary_iso_curve,
    uv_line,
};
use crate::error::SweepError;
use crate::topology::caps::create_cap_face_with_inners;

/// [`build_sweep_topology`] with per-wire span grouping: a closed
/// profile with more than one wire routes to the multi-wire builder
/// (one side-face ring per wire, caps with inner loops — PK-477);
/// everything else takes the single-wire path unchanged.
#[allow(clippy::too_many_arguments)]
pub fn build_sweep_topology_wired(
    surfaces: &[NurbsSurface],
    wires: &[usize],
    start_cap_curves: &[NurbsCurve],
    end_cap_curves: &[NurbsCurve],
    is_closed_profile: bool,
    tolerance: Tolerance,
    tracking: Option<&mut TrackingContext>,
) -> Result<BRep, SweepError> {
    if is_closed_profile && wires.len() > 1 {
        let mut brep = build_closed_multi_wire_topology(
            surfaces,
            wires,
            start_cap_curves,
            end_cap_curves,
            SideFaceGeometry::ExtrudedCanonical,
            tolerance,
        )?;
        if let Some(ctx) = tracking {
            super::side_faces::assign_sweep_names(&mut brep, ctx);
        }
        return Ok(brep);
    }
    super::side_faces::build_sweep_topology(
        surfaces,
        start_cap_curves,
        end_cap_curves,
        is_closed_profile,
        tolerance,
        tracking,
    )
}

/// Per-wire cap boundary: the cap-side half-edges and their curves,
/// in span order.
struct WireCapRing {
    half_edges: Vec<HalfEdgeHandle>,
    curves: Vec<NurbsCurve>,
}

/// Build a closed multi-wire sweep solid: `wires` gives the span
/// count per wire (outer first) over the flat `surfaces` /
/// `*_cap_curves` lists.
pub(crate) fn build_closed_multi_wire_topology(
    surfaces: &[NurbsSurface],
    wires: &[usize],
    start_cap_curves: &[NurbsCurve],
    end_cap_curves: &[NurbsCurve],
    geometry: SideFaceGeometry,
    tolerance: Tolerance,
) -> Result<BRep, SweepError> {
    let total: usize = wires.iter().sum();
    if total != surfaces.len() {
        return Err(SweepError::Topology(format!(
            "multi-wire span census mismatch: {} spans across wires, {} surfaces",
            total,
            surfaces.len()
        )));
    }
    let mut builder = BRepBuilder::new(tolerance);
    let solid = builder.add_solid();
    let shell = builder.add_shell(solid);

    let mut start_rings: Vec<WireCapRing> = Vec::with_capacity(wires.len());
    let mut end_rings: Vec<WireCapRing> = Vec::with_capacity(wires.len());
    let mut offset = 0usize;
    for (wi, &n_spans) in wires.iter().enumerate() {
        if n_spans == 0 {
            return Err(SweepError::Topology(format!(
                "multi-wire profile: wire {wi} has no spans"
            )));
        }
        let wire_surfaces = &surfaces[offset..offset + n_spans];
        let wire_start_curves = &start_cap_curves[offset..offset + n_spans];
        let wire_end_curves = &end_cap_curves[offset..offset + n_spans];
        let (start_ring, end_ring) = build_wire_ring(
            &mut builder,
            shell,
            wire_surfaces,
            wire_start_curves,
            wire_end_curves,
            geometry,
            wi,
        )?;
        start_rings.push(start_ring);
        end_rings.push(end_ring);
        offset += n_spans;
    }

    // Start cap: boundary half-edges run sv_next → sv_this per wire,
    // so span order chains only when reversed (same rule as the
    // single-wire builder; the winding roles come out opposed because
    // the inner wires' profiles are CW where the outer is CCW).
    let rev = |ring: &WireCapRing| -> (Vec<HalfEdgeHandle>, Vec<NurbsCurve>) {
        (
            ring.half_edges.iter().rev().copied().collect(),
            ring.curves.iter().rev().cloned().collect(),
        )
    };
    let (outer_hes, outer_curves) = rev(&start_rings[0]);
    let start_inners: Vec<(Vec<HalfEdgeHandle>, Vec<NurbsCurve>)> =
        start_rings[1..].iter().map(&rev).collect();
    create_cap_face_with_inners(
        &mut builder,
        shell,
        &outer_hes,
        &outer_curves,
        &start_inners,
        true,
    )?;

    let end_inners: Vec<(Vec<HalfEdgeHandle>, Vec<NurbsCurve>)> = end_rings[1..]
        .iter()
        .map(|r| (r.half_edges.clone(), r.curves.clone()))
        .collect();
    create_cap_face_with_inners(
        &mut builder,
        shell,
        &end_rings[0].half_edges,
        &end_rings[0].curves,
        &end_inners,
        false,
    )?;

    Ok(builder.build_unchecked())
}

/// Build one wire's side-face ring; returns its start/end cap rings.
fn build_wire_ring(
    builder: &mut BRepBuilder,
    shell: popkern_core::topology::ShellHandle,
    surfaces: &[NurbsSurface],
    start_cap_curves: &[NurbsCurve],
    end_cap_curves: &[NurbsCurve],
    geometry: SideFaceGeometry,
    wire_index: usize,
) -> Result<(WireCapRing, WireCapRing), SweepError> {
    let n_spans = surfaces.len();
    let mut start_vertices = Vec::with_capacity(n_spans);
    let mut end_vertices = Vec::with_capacity(n_spans);
    for surface in surfaces {
        start_vertices.push(builder.add_vertex(eval_surface_corner(surface, 0.0, 0.0)));
        end_vertices.push(builder.add_vertex(eval_surface_corner(surface, 0.0, 1.0)));
    }

    // One longitudinal edge per profile vertex, shared between the
    // adjacent side faces of THIS wire.
    let mut long_up: Vec<HalfEdgeHandle> = Vec::with_capacity(n_spans);
    let mut long_down: Vec<HalfEdgeHandle> = Vec::with_capacity(n_spans);
    for j in 0..n_spans {
        let (edge, he_up, he_down) = builder.add_edge(start_vertices[j], end_vertices[j]);
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

    let mut start_ring = WireCapRing {
        half_edges: Vec::with_capacity(n_spans),
        curves: start_cap_curves.to_vec(),
    };
    let mut end_ring = WireCapRing {
        half_edges: Vec::with_capacity(n_spans),
        curves: end_cap_curves.to_vec(),
    };

    for i in 0..n_spans {
        let next = (i + 1) % n_spans;
        let sv_this = start_vertices[i];
        let sv_next = start_vertices[next];
        let ev_this = end_vertices[i];
        let ev_next = end_vertices[next];

        let (bottom_edge, he_bottom_fwd, he_bottom_rev) = builder.add_edge(sv_this, sv_next);
        attach_nurbs_curve_to_edge(builder, bottom_edge, start_cap_curves, i);
        let (top_edge, he_top_fwd, he_top_rev) = builder.add_edge(ev_this, ev_next);
        attach_nurbs_curve_to_edge(builder, top_edge, end_cap_curves, i);

        let he_left_up = long_up[i];
        let he_left_down = long_down[i];
        let he_right_up = long_up[next];

        let face = builder.add_face(shell);
        let recognized = match geometry {
            SideFaceGeometry::ExtrudedCanonical => {
                crate::surface::canonical::recognize_extruded_plane(&surfaces[i]).or_else(|| {
                    crate::surface::canonical::recognize_extruded_cylinder(&surfaces[i])
                })
            }
            SideFaceGeometry::Nurbs => None,
        };
        let (surface_geometry, chart) = recognized
            .map(|(g, c)| (g, Some(c)))
            .unwrap_or_else(|| (SurfaceGeometry::Nurbs(surfaces[i].clone()), None));
        let orientation =
            crate::surface::canonical::canonical_normal_agreement(&surfaces[i], &surface_geometry);
        builder.set_face_surface(
            face,
            SurfaceOnFace {
                geometry: surface_geometry,
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

        let loop_hes = if n_spans == 1 {
            vec![he_bottom_fwd, he_left_up, he_top_rev, he_left_down]
        } else {
            vec![he_bottom_fwd, he_right_up, he_top_rev, he_left_down]
        };
        builder
            .add_loop(&loop_hes, face)
            .map_err(|e| SweepError::Topology(format!("wire {wire_index} side face loop: {e}")))?;

        builder.set_half_edge_pcurve(he_bottom_fwd, uv_line_chart(0.0, 0.0, 1.0, 0.0));
        builder.set_half_edge_pcurve(he_top_rev, uv_line_chart(1.0, 1.0, 0.0, 1.0));
        builder.set_half_edge_pcurve(he_left_down, uv_line_chart(0.0, 1.0, 0.0, 0.0));
        if n_spans == 1 {
            // Single-span wire: both longitudinal boundaries are one
            // seam edge; NURBS walls need the chart-continuous u=1
            // half (see the single-wire builder, PK-441).
            let seam_u = if chart.is_some() { 0.0 } else { 1.0 };
            builder.set_half_edge_pcurve(he_left_up, uv_line_chart(seam_u, 0.0, seam_u, 1.0));
        } else {
            builder.set_half_edge_pcurve(he_right_up, uv_line_chart(1.0, 0.0, 1.0, 1.0));
        }

        start_ring.half_edges.push(he_bottom_rev);
        end_ring.half_edges.push(he_top_fwd);
    }

    Ok((start_ring, end_ring))
}
