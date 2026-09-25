//! The saddle setback ring (PK-527; design
//! `docs/architecture/saddle-corner-setback.md`).
//!
//! A saddle's bridged setback border is a ring of strip arcs meeting
//! at the convex faces' rail crossings, plus a rail-flat bridge on each
//! reflex face. The pre-carve leaves a mitre complex at a saddle: rim
//! chords to the sharp edges, and rails replaced by off-surface
//! chords. The setback surgery's crossing carve assumes clean strip
//! faces, so it cannot place arcs there. The ring carve
//! (`boundary_ring`) resolves exactly that complex: a rail split, or a
//! rail-spur recut where a chord replaced the rail. So the hybrid ring
//! is carved there, side by side:
//! - strips carry their exact section rims;
//! - bridges carry their rail-flat curves;
//! - consecutive strips share their corner vertex.
//!
//! The cap fan then commits on the merged blend face exactly as on a
//! convex corner.

use popkern_core::geometry::NurbsCurve;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, VertexHandle};

use super::Kp20SetbackArc;
use super::boundary_ring::{RingSide, carve_ring_sides};
use crate::fillet_surface::FilletSurfaceResult;

/// Carve the bridged setback ring at `vertex`. Returns the blend face
/// and each side's ring edge, in side order.
pub(super) fn carve_saddle_ring(
    brep: &mut BRep,
    vertex: VertexHandle,
    vertex_position: popkern_core::geometry::Point3,
    sides: &[Kp20SetbackArc],
    fillet_surfaces: &[FilletSurfaceResult],
    tolerance: &Tolerance,
) -> Result<(FaceHandle, Vec<EdgeHandle>), String> {
    let curves: Vec<NurbsCurve> = sides
        .iter()
        .map(|side| match &side.bridge {
            Some(bridge) => bridge.curve.clone(),
            None => side.section_curve(),
        })
        .collect();
    let ring_sides: Vec<RingSide<'_>> = sides
        .iter()
        .zip(&curves)
        .map(|(side, curve)| match side.bridge {
            Some(_) => RingSide::Bridge { curve },
            None => RingSide::Arc {
                fillet_index: side.fillet_index,
                curve,
            },
        })
        .collect();
    let ring = carve_ring_sides(
        brep,
        vertex,
        vertex_position,
        &ring_sides,
        fillet_surfaces,
        tolerance,
    )?;
    Ok((ring.blend_face, ring.side_edges))
}
