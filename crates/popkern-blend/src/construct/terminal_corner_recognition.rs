//! Stale chain-end corner RECOGNITION for the terminal settle:
//! iso-station stub gates, chart spans, and the off-surface
//! signature per strip surface kind.
//!
//! Split from `terminal_arc_trim.rs` to honour the structural cap.

use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle, HalfEdgeHandle, VertexHandle};
use popkern_math::PortableMath;

use super::terminal_arc_trim::{CornerCandidate, PlanarStubRecognition};

/// Walk one face's outer loop and recognize stale chain-end corners:
/// two consecutive half-edges whose strip pcurves sit on one iso-v
/// station while the shared vertex lies off the strip surface.
pub(super) fn recognize_corners(
    brep: &BRep,
    face: FaceHandle,
    tolerance: &Tolerance,
    planar_mode: PlanarStubRecognition,
) -> Vec<CornerCandidate> {
    let mut found = Vec::new();
    let Some(f) = brep.faces.get(face) else {
        return found;
    };
    let Some(surface) = f.surface.as_ref() else {
        return found;
    };
    let Some(l) = brep.loops.get(f.outer_loop) else {
        return found;
    };
    let start = l.half_edge;
    let mut he_in = start;
    // Bounded walk: loops are cycles; cap at half-edge count.
    for _ in 0..brep.half_edges.len() {
        let Some(a) = brep.half_edges.get(he_in) else {
            break;
        };
        let he_out = a.next;
        let Some(b) = brep.half_edges.get(he_out) else {
            break;
        };
        if let Some(candidate) = recognize_corner_pair(
            brep,
            face,
            &surface.geometry,
            he_in,
            he_out,
            b.origin,
            tolerance,
            planar_mode,
        ) {
            found.push(candidate);
        }
        he_in = he_out;
        if he_in == start {
            break;
        }
    }
    found
}

/// UV endpoints of a half-edge's strip-side pcurve, oriented along
/// the half-edge (origin → destination).
pub(super) fn strip_pcurve_span(
    brep: &BRep,
    face: FaceHandle,
    he: HalfEdgeHandle,
) -> Option<((f64, f64), (f64, f64))> {
    let pc = brep
        .face_use_pcurve(face, he)
        .or(brep.half_edges.get(he)?.pcurve.as_ref())?;
    let (t0, t1) = pc.parameter_range();
    let a = pc.point_at(t0);
    let b = pc.point_at(t1);
    Some(((a.x, a.y), (b.x, b.y)))
}

#[allow(clippy::too_many_arguments)]
fn recognize_corner_pair(
    brep: &BRep,
    face: FaceHandle,
    surface: &SurfaceGeometry,
    he_in: HalfEdgeHandle,
    he_out: HalfEdgeHandle,
    corner: VertexHandle,
    tolerance: &Tolerance,
    planar_mode: PlanarStubRecognition,
) -> Option<CornerCandidate> {
    // A planar strip (chamfer band) needs no chart gate at all: its
    // loop must lie in its plane, so an off-plane corner is a defect
    // by definition regardless of stub chart direction. The iso
    // gates below reject SKEWED stubs — a bevel meeting its end face
    // obliquely on a yawed+tilted boolean-cut rest plane (PK-168:
    // the strip loop folded through the uncut corner by d/sqrt(2)
    // and every mesh gate stayed green; armed check 9 exposed it) —
    // and the downstream guards (off-surface signature, common
    // planar cap, vertex isolation, anchors) fully determine the
    // candidate on planes.
    let surface_scoped_plane = (matches!(surface, SurfaceGeometry::Plane { .. })
        && planar_mode == PlanarStubRecognition::SurfaceScoped)
        // CONE strips (circular chamfers): the stub pcurves in the
        // cone chart don't satisfy the cylinder iso gates (and are
        // often absent on freshly-cut interrupted-rim terminals), so
        // chart gating starves recognition entirely. The downstream
        // guards — off-surface signature, vertex isolation, common
        // planar cap, anchors-on-cone — fully determine the
        // candidate, exactly as for surface-scoped planes (PK-171:
        // the notched flange's stale rim-end corner was never even
        // considered).
        || matches!(surface, SurfaceGeometry::Cone { .. });
    let mut shared_chart_corner = None;
    if !surface_scoped_plane {
        let (a_start, a_end) = strip_pcurve_span(brep, face, he_in)?;
        let (b_start, b_end) = strip_pcurve_span(brep, face, he_out)?;
        shared_chart_corner = Some(a_end);
        // Both stubs on one iso-v station, meeting at the same u. The
        // chart values may be stale — they gate recognition only; all
        // repair positions are computed from trusted geometry.
        // Recognition-only gate, SCALE-FREE: swept-strip stub pcurves
        // projected from a stale corner drift off the exact iso by the
        // corner's chart shadow, and the shadow grows with the end
        // radius (measured 1.4e-3 in u at r=0.10, 6.5e-3 at r=0.20), so
        // an absolute station tolerance cannot serve every radius. Gate
        // on the SLOPE instead — each stub's variation must be dominated
        // by one shared chart direction — plus a meeting check relative
        // to the spans. Cylinder strips author iso-v stubs, swept
        // strips near-iso with shadow drift; every ACCEPT is still
        // guarded downstream by off-surface, isolation, common-cap,
        // circularity, and anchor checks.
        let d_a = (a_end.0 - a_start.0, a_end.1 - a_start.1);
        let d_b = (b_end.0 - b_start.0, b_end.1 - b_start.1);
        let iso_u = d_a.0.abs() <= 0.2 * d_a.1.abs() && d_b.0.abs() <= 0.2 * d_b.1.abs();
        let iso_v = d_a.1.abs() <= 0.2 * d_a.0.abs() && d_b.1.abs() <= 0.2 * d_b.0.abs();
        let meet_eps = 1e-6 + 0.25 * (d_a.0.abs() + d_a.1.abs() + d_b.0.abs() + d_b.1.abs());
        let meets =
            (a_end.0 - b_start.0).abs() <= meet_eps && (a_end.1 - b_start.1).abs() <= meet_eps;
        if !((iso_v || iso_u) && meets) {
            return None;
        }
    }
    // The vertex must be OFF the strip surface — that is the stale
    // signature. Coherent corners (already settled) are left alone,
    // which keeps the pass idempotent. Measured against the true
    // support geometry, not the possibly-stale chart.
    let corner_pos = brep.vertices.get(corner)?.position;
    if popkern_core::trace::trace_enabled("terminal-arc")
        && matches!(surface, SurfaceGeometry::Cone { .. })
    {
        eprintln!(
            "[tat]   cone corner probe at ({:.3},{:.3},{:.3})",
            corner_pos.x, corner_pos.y, corner_pos.z
        );
    }
    let off_surface = match surface {
        SurfaceGeometry::Cylinder {
            origin,
            axis,
            radius,
        } => {
            let d = corner_pos - *origin;
            let radial = d - *axis * d.dot(axis);
            (radial.length() - radius).abs() > tolerance.linear() * 10.0
        }
        SurfaceGeometry::Plane { origin, normal, .. } => {
            (corner_pos - *origin).dot(normal).abs() > tolerance.linear() * 10.0
        }
        SurfaceGeometry::Cone {
            apex,
            axis,
            half_angle,
        } => {
            let an = axis.normalized()?;
            let d = corner_pos - *apex;
            let station = d.dot(&an);
            let radial = (d - an * station).length();
            (radial - station.abs() * half_angle.ptan().abs()).abs() > tolerance.linear() * 10.0
        }
        SurfaceGeometry::Torus {
            center,
            axis,
            major_radius,
            minor_radius,
        } => {
            let a = axis.normalized()?;
            super::terminal_torus_section::torus_distance(
                &corner_pos,
                center,
                &a,
                *major_radius,
                *minor_radius,
            ) > tolerance.linear() * 10.0
        }
        SurfaceGeometry::Nurbs(_) => {
            // Swept (variable-radius) strips: the chart + surface are
            // the construction authority; the stale signature is the
            // corner vertex sitting off the surface point at the
            // stubs' shared chart corner.
            use popkern_core::geometry::SurfaceEval;
            let shared = shared_chart_corner?;
            let on_surface = surface.point_at(shared.0, shared.1);
            corner_pos.distance_to(&on_surface) > tolerance.linear() * 10.0
        }
        _ => return None,
    };
    if !off_surface {
        return None;
    }
    Some(CornerCandidate {
        strip_face: face,
        he_in,
        he_out,
        corner,
    })
}
