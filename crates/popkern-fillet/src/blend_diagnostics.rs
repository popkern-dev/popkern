//! UI-forward, one-call diagnostic helpers for a completed blend.
//!
//! Each function here takes a [`BRep`] and a blend-face [`FaceHandle`]
//! and returns a diagnostic ready for a UI layer to render without
//! the caller having to:
//!
//! - resolve the surface geometry from the face
//! - compute a finite parameter range from the face's vertices
//! - locate the adjacent support faces topologically
//!
//! The primitives live in [`crate::curvature_comb`] and
//! [`crate::highlight_line`]; this module glues them to the BRep so
//! a GUI can call one function with one face handle.

use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{NurbsCurve, Vec3};
use popkern_core::topology::{BRep, FaceHandle};

use crate::curvature_comb::{CombAlong, CurvatureComb, surface_isoline_curvature_comb};
use crate::highlight_line::highlight_boundary_kink;

/// Curvature comb along the mid-`u` iso-line of a blend face (sweeps
/// `v` with `u = u_mid`).
///
/// For a [`SurfaceGeometry::Nurbs`] blend built by the G2 sweep,
/// `v` IS the cross-section direction, so this comb walks across
/// the fillet width — showing the cross-section's curvature
/// profile. For an analytic [`SurfaceGeometry::Cylinder`]
/// (straight-edge plane-plane fillet), `v` is the cylinder axis
/// direction, so the comb reads κ = 0 along the spine.
///
/// Use [`blend_curvature_comb_v`] for the same comb shape on a
/// cylinder-output fillet (or on any surface where the cross-
/// section is parameterised by `u`).
///
/// Returns `None` if `face` has no surface geometry.
pub fn blend_curvature_comb_v(
    brep: &BRep,
    face: FaceHandle,
    n_samples: usize,
) -> Option<CurvatureComb> {
    let surface_on_face = brep.faces.get(face)?.surface.as_ref()?;
    let (u_range, v_range) =
        crate::geometry_compute::face_param_range(brep, face, &surface_on_face.geometry);
    let u_mid = 0.5 * (u_range.0 + u_range.1);
    Some(surface_isoline_curvature_comb(
        &surface_on_face.geometry,
        CombAlong::V,
        u_mid,
        v_range,
        n_samples,
    ))
}

/// Curvature comb along the mid-`v` iso-line of a blend face (sweeps
/// `u` with `v = v_mid`).
///
/// For analytic [`SurfaceGeometry::Cylinder`] blend surfaces this
/// is the angular direction — the cross-section — so the comb reads
/// κ = 1/r. For G2 NURBS sweeps `u` is the along-spine direction;
/// the comb reads the rolling-ball centre trace's curvature.
///
/// Companion to [`blend_curvature_comb_v`]; use whichever matches
/// the blend's parameterisation. The G2 NURBS fingerprint is
/// `degree_v == 5` — see
/// [`crate::quality::BlendQuality::achieved_continuity`].
///
/// Returns `None` if `face` has no surface geometry.
pub fn blend_curvature_comb_u(
    brep: &BRep,
    face: FaceHandle,
    n_samples: usize,
) -> Option<CurvatureComb> {
    let surface_on_face = brep.faces.get(face)?.surface.as_ref()?;
    let (u_range, v_range) =
        crate::geometry_compute::face_param_range(brep, face, &surface_on_face.geometry);
    let v_mid = 0.5 * (v_range.0 + v_range.1);
    Some(surface_isoline_curvature_comb(
        &surface_on_face.geometry,
        CombAlong::U,
        v_mid,
        u_range,
        n_samples,
    ))
}

/// §15 highlight-line boundary kink between a blend face and each
/// of its adjacent support faces.
///
/// Walks every edge on `blend_face`'s outer loop, picks the
/// other-side face handle, and calls [`highlight_boundary_kink`]
/// for each unique support face. Returns the **maximum** kink angle
/// (in radians) found across all boundary edges — the single
/// "worst-case" Class-A metric for this blend.
///
/// - G2-continuous blends produce `≈ 0` (the tangent direction of
///   the highlight line is continuous across the boundary).
/// - G1 blends produce small non-zero values (a few degrees).
/// - G0 (positional-only) blends produce large angles up to π/2.
///
/// Returns `None` if `blend_face` has no adjacent faces with surface
/// geometry, or if the kink analysis produces no matching sample
/// pairs at any boundary.
pub fn blend_highlight_kink_max(
    brep: &BRep,
    blend_face: FaceHandle,
    light_dir: Vec3,
    view_dir: Vec3,
    grid_density: usize,
    boundary_tol: f64,
) -> Option<f64> {
    let face = brep.faces.get(blend_face)?;
    let outer_loop = brep.loops.get(face.outer_loop)?;
    let start_he = outer_loop.half_edge;

    // Collect unique adjacent support face handles by walking the
    // outer loop's half-edges and reading their twin's face.
    let mut supports: Vec<FaceHandle> = Vec::new();
    let mut cursor = start_he;
    let bound = brep.half_edges.len().max(1);
    for _ in 0..=bound {
        let he = brep.half_edges.get(cursor)?;
        let twin_h = he.twin;
        if twin_h != popkern_core::topology::HalfEdgeHandle::default()
            && let Some(twin) = brep.half_edges.get(twin_h)
            && let Some(loop_) = brep.loops.get(twin.loop_)
        {
            let f = loop_.face;
            if f != blend_face
                && !supports.contains(&f)
                && brep
                    .faces
                    .get(f)
                    .and_then(|ff| ff.surface.as_ref())
                    .is_some()
            {
                supports.push(f);
            }
        }
        cursor = he.next;
        if cursor == start_he {
            break;
        }
    }

    let mut worst: Option<f64> = None;
    for support in supports {
        if let Some(angle) = highlight_boundary_kink(
            brep,
            blend_face,
            support,
            light_dir,
            view_dir,
            grid_density,
            boundary_tol,
        ) {
            worst = Some(match worst {
                Some(m) => m.max(angle),
                None => angle,
            });
        }
    }
    worst
}

/// Best-effort support pair lookup helper used by tests. Returns the
/// two neighbour faces of `blend_face` that carry surface geometry
/// distinct from the blend itself, suitable for passing to the
/// `g2_cross_section` primitives.
#[allow(dead_code)]
fn blend_support_pair(brep: &BRep, blend_face: FaceHandle) -> Option<(FaceHandle, FaceHandle)> {
    let face = brep.faces.get(blend_face)?;
    let outer_loop = brep.loops.get(face.outer_loop)?;
    let start_he = outer_loop.half_edge;
    let mut supports: Vec<FaceHandle> = Vec::new();
    let mut cursor = start_he;
    let bound = brep.half_edges.len().max(1);
    for _ in 0..=bound {
        let he = brep.half_edges.get(cursor)?;
        let twin_h = he.twin;
        if twin_h != popkern_core::topology::HalfEdgeHandle::default()
            && let Some(twin) = brep.half_edges.get(twin_h)
            && let Some(loop_) = brep.loops.get(twin.loop_)
        {
            let f = loop_.face;
            if f != blend_face && !supports.contains(&f) && brep.faces.contains_key(f) {
                supports.push(f);
                if supports.len() >= 2 {
                    break;
                }
            }
        }
        cursor = he.next;
        if cursor == start_he {
            break;
        }
    }
    match supports.as_slice() {
        [a, b] => Some((*a, *b)),
        _ => None,
    }
}

/// Silence unused-NurbsCurve import warning — `curvature_comb` re-
/// exports it transitively for surface isolines.
#[allow(dead_code)]
fn _unused_import_pacifier(_c: NurbsCurve, _v: SurfaceGeometry) {}

#[cfg(test)]
#[path = "blend_diagnostics_tests.rs"]
mod tests;
