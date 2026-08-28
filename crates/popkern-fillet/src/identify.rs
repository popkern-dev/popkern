//! Blend identification and query API (§14.1, §14.3).
//!
//! Given a [`BRep`] face, this module answers the commercially
//! essential question "is this face part of a blend, and if so,
//! what kind and what parameters?". The engine uses purely
//! geometric analysis — no metadata from the create-time
//! operation is required — so identification works on imported
//! models (STEP, IGES, stray BReps from unknown sources) as well
//! as models built inside popkern.
//!
//! ## How it works
//!
//! Identification combines two diagnostics:
//!
//! 1. **Boundary tangency.** For each boundary edge of the
//!    candidate face, the angle between the face normal and the
//!    neighbour face normal is computed. A rolling-ball fillet is
//!    G1-tangent along its support-face boundaries, so at least
//!    two boundary edges will report an angle near zero. A chamfer
//!    face is G0 on both support boundaries, so no boundary edge
//!    is tangent. An original un-blended face has every boundary
//!    at a sharp corner.
//!
//! 2. **Surface type.** The analytic surface kind names the blend:
//!    `Torus` → rolling-ball fillet on a curved edge between two
//!    planar supports, `Cylinder` → rolling-ball fillet on a
//!    straight edge (plane-plane) or on a plane + cylindrical
//!    support, `Sphere` → vertex blend, `Plane` with non-tangent
//!    boundaries → chamfer.
//!
//! Combining both gives per-identification confidence in `[0, 1]`:
//! 0.9 for an analytic torus section (strong match), 0.75 for a
//! cylinder-based fillet, 0.85 for a spherical vertex blend, 0.7
//! for a chamfer, 0.5 for a NURBS sweep with G1 boundaries.
//! Identifications with confidence below 0.5 return `None`.
//!
//! ## Limitations
//!
//! This round uses **only** curvature/tangency analysis —
//! metadata-tagged identification (where the creating operation
//! records "I made this face as a radius-0.1 rolling-ball fillet")
//! is a future extension. Curvature analysis is less
//! reliable on heavily approximated NURBS fillets and on edge cases
//! like variable-radius blends where the cross-section curvature
//! changes along the spine.

use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::topology::{BRep, FaceHandle};

use crate::identify_helpers::{
    all_neighbours_are_planes, boundary_tangency, chamfer_support_pair, distinct_neighbours,
    estimate_chamfer_half_width, tangent_support_pair,
};

/// Classification of an identified blend.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BlendKind {
    /// Rolling-ball fillet with a circular cross-section.
    Fillet,
    /// Chamfer with a straight-line cross-section.
    Chamfer,
    /// Vertex blend — an N-sided patch at a multi-fillet junction.
    VertexBlend,
}

/// Where the identification came from.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IdentificationSource {
    /// The engine recognised the blend from metadata attached at
    /// creation time. Authoritative — the parameters are exactly
    /// what the caller specified. *(Not populated in this round;
    /// future work.)*
    Metadata,
    /// The engine recognised the blend from curvature and
    /// tangency analysis on the face's geometry. Best-effort
    /// heuristic — see the [`BlendInfo::confidence`] field.
    CurvatureAnalysis,
}

/// A complete description of a blend recognised on a face.
#[derive(Debug, Clone)]
pub struct BlendInfo {
    /// The face the identification applies to.
    pub face: FaceHandle,
    /// Classification of the blend.
    pub kind: BlendKind,
    /// Nominal radius of the rolling ball (for `Fillet` /
    /// `VertexBlend`) or half-width of the chamfer face (for
    /// `Chamfer`). `None` when the engine could not extract a
    /// reliable value (typically NURBS approximations).
    pub radius: Option<f64>,
    /// The two support faces the blend bridges. For `Fillet` and
    /// `Chamfer` this is the primary support pair. For
    /// `VertexBlend` this is the pair of **adjacent** blend
    /// faces that most obviously touch the vertex blend — the
    /// full neighbour set can be retrieved via
    /// [`blend_neighbours`].
    pub support_faces: Option<(FaceHandle, FaceHandle)>,
    /// Where the identification came from.
    pub source: IdentificationSource,
    /// Confidence in `[0, 1]`. 1.0 for metadata-backed
    /// identifications; curvature-analysis identifications range
    /// from 0.5 (weak) to 0.95 (strong analytic surface match).
    pub confidence: f64,
    /// Cap provenance when the face ships as a first-class
    /// polyhedral-net spline sector (PK-200): which cap (valence)
    /// and which sector of it this face is. `None` for every other
    /// blend kind and for caps that shipped as anonymous surfaces.
    pub cap_provenance: Option<popkern_core::geometry::PolyhedralSplineProvenance>,
}

/// Identify a blend on the given face.
///
/// Returns `Some(BlendInfo)` when the engine recognises the face
/// as a blend with confidence ≥ 0.5, or `None` otherwise.
///
/// # Example
///
/// ```ignore
/// let brep_before = BRep::unit_cube(tol).unwrap();
/// let edge = brep_before.edges.keys().next().unwrap();
/// let result = popkern_fillet::fillet_edge(&brep_before, edge, 0.1, tol, None).unwrap();
///
/// for (fh, _) in result.brep.faces.iter() {
///     if let Some(info) = popkern_fillet::identify_blend(&result.brep, fh) {
///         println!("Face {fh:?} is a {:?} with radius {:?}", info.kind, info.radius);
///     }
/// }
/// ```
pub fn identify_blend(brep: &BRep, face: FaceHandle) -> Option<BlendInfo> {
    // Unknown face handles return None immediately.
    let face_data = brep.faces.get(face)?;
    let surface_on_face = face_data.surface.as_ref()?;
    let geometry = &surface_on_face.geometry;

    // Collect boundary data: (neighbour face, edge, midpoint, dihedral angle).
    let boundary = boundary_tangency(brep, face)?;
    if boundary.is_empty() {
        return None;
    }

    // Count tangent and sharp boundary edges. The angles come from
    // 3D-projected sample points; topologically-noisy boundaries
    // (split contact curves, NURBS approximations) can show as
    // "almost tangent" rather than "perfectly tangent", so we use a
    // generous tolerance.
    let tangent_tol = 0.05; // radians ≈ 2.9°
    if std::env::var_os("PK_IDENTIFY_TRACE").is_some() {
        for b in &boundary {
            let kind = brep
                .faces
                .get(b.neighbour)
                .and_then(|f| f.surface.as_ref())
                .map(|s| format!("{:?}", std::mem::discriminant(&s.geometry)));
            eprintln!(
                "[identify] face {face:?} nb {:?} kind {kind:?} angle {:.4}",
                b.neighbour, b.angle
            );
        }
    }
    let tangent_count = boundary
        .iter()
        .filter(|b| b.angle.abs() < tangent_tol)
        .count();
    let perp_count = boundary
        .iter()
        .filter(|b| (b.angle - std::f64::consts::FRAC_PI_2).abs() < tangent_tol)
        .count();
    let mixed_count = boundary.len() - tangent_count - perp_count;
    let distinct_neighbours = distinct_neighbours(&boundary);

    // -----------------------------------------------------------------
    // Classification rules — surface type is the primary signal.
    // -----------------------------------------------------------------
    //
    // Strategy: analytic blend surfaces (Torus, Sphere) are
    // diagnostic on their own. For a `Cylinder` or `Plane` face we
    // need at least one tangent boundary (cylinder) or a non-trivial
    // mixed-angle bridge pattern (plane → chamfer) to distinguish
    // it from a primitive that happens to share the same surface
    // type.
    match geometry {
        SurfaceGeometry::Torus { minor_radius, .. } => {
            // A torus face in a BRep is almost always a fillet
            // strip — toroidal primitives are vanishingly rare in
            // CAD models. High confidence.
            let support_faces = tangent_support_pair(&boundary, tangent_tol)
                .or_else(|| chamfer_support_pair(&boundary, &distinct_neighbours));
            Some(BlendInfo {
                face,
                kind: BlendKind::Fillet,
                radius: Some(*minor_radius),
                support_faces,
                source: IdentificationSource::CurvatureAnalysis,
                confidence: 0.9,
                cap_provenance: None,
            })
        }
        SurfaceGeometry::Sphere { radius, .. } => {
            // A spherical face with 3+ neighbours is a vertex blend
            // (smooth N-sided cap). Otherwise it's a fillet variant
            // or a sphere primitive — call it a fillet with lower
            // confidence so the caller can decide.
            let kind = if distinct_neighbours.len() >= 3 {
                BlendKind::VertexBlend
            } else {
                BlendKind::Fillet
            };
            let confidence = if matches!(kind, BlendKind::VertexBlend) {
                0.85
            } else {
                0.65
            };
            Some(BlendInfo {
                face,
                kind,
                radius: Some(*radius),
                support_faces: tangent_support_pair(&boundary, tangent_tol)
                    .or_else(|| chamfer_support_pair(&boundary, &distinct_neighbours)),
                source: IdentificationSource::CurvatureAnalysis,
                confidence,
                cap_provenance: None,
            })
        }
        SurfaceGeometry::Cylinder { radius, .. } => {
            // A cylinder face is a fillet in two scenarios:
            //
            // 1. It has at least one G1-tangent boundary to a
            //    neighbour (the classic cylinder-on-curved-support
            //    case, e.g. plane + cylindrical hole).
            //
            // 2. All distinct neighbours are planes. This is the
            //    straight-edge plane-plane fillet (§11.1): the
            //    rolling ball sweeps along a straight dihedral edge
            //    and the resulting surface is an exact cylinder
            //    whose radius equals the fillet radius. The
            //    boundaries are *not* G1-tangent (the cylinder
            //    meets each plane at the rolling-ball contact
            //    angle), so the tangent_count check alone would
            //    miss it. Confidence is the same as the tangent
            //    path (0.75): a cylinder between two planes is a
            //    strong geometric signal.
            //
            // Otherwise the cylinder is a primitive (hole, shaft)
            // and we return None.
            if tangent_count >= 1 {
                Some(BlendInfo {
                    face,
                    kind: BlendKind::Fillet,
                    radius: Some(*radius),
                    support_faces: tangent_support_pair(&boundary, tangent_tol),
                    source: IdentificationSource::CurvatureAnalysis,
                    confidence: 0.75,
                    cap_provenance: None,
                })
            } else if all_neighbours_are_planes(brep, &boundary) {
                let support_faces = chamfer_support_pair(&boundary, &distinct_neighbours);
                Some(BlendInfo {
                    face,
                    kind: BlendKind::Fillet,
                    radius: Some(*radius),
                    support_faces,
                    source: IdentificationSource::CurvatureAnalysis,
                    confidence: 0.75,
                    cap_provenance: None,
                })
            } else {
                None
            }
        }
        SurfaceGeometry::Plane { .. } => {
            // A plane face is a chamfer only if it has at least
            // one mixed-angle boundary (the chamfer-angle edge to
            // a support face) AND no tangent boundaries (which
            // would indicate a fillet-modified original face).
            // Original cube-style faces have only perpendicular
            // boundaries and are correctly rejected.
            if tangent_count == 0 && mixed_count >= 1 {
                let support_faces = chamfer_support_pair(&boundary, &distinct_neighbours);
                let radius = estimate_chamfer_half_width(brep, face);
                Some(BlendInfo {
                    face,
                    kind: BlendKind::Chamfer,
                    radius,
                    support_faces,
                    source: IdentificationSource::CurvatureAnalysis,
                    confidence: 0.7,
                    cap_provenance: None,
                })
            } else {
                None
            }
        }
        SurfaceGeometry::Nurbs(_) => {
            // NURBS faces can be either fillet sweeps (when at
            // least two boundaries are G1-tangent to their
            // neighbours — the rolling-ball signature) or chamfer
            // ruled surfaces (when at least one boundary shows the
            // chamfer-angle pattern and no boundary is tangent).
            if tangent_count >= 2 {
                Some(BlendInfo {
                    face,
                    kind: BlendKind::Fillet,
                    radius: None,
                    support_faces: tangent_support_pair(&boundary, tangent_tol),
                    source: IdentificationSource::CurvatureAnalysis,
                    confidence: 0.5,
                    cap_provenance: None,
                })
            } else if tangent_count == 0 && mixed_count >= 1 {
                let support_faces = chamfer_support_pair(&boundary, &distinct_neighbours);
                let radius = estimate_chamfer_half_width(brep, face);
                Some(BlendInfo {
                    face,
                    kind: BlendKind::Chamfer,
                    radius,
                    support_faces,
                    source: IdentificationSource::CurvatureAnalysis,
                    confidence: 0.6,
                    cap_provenance: None,
                })
            } else {
                None
            }
        }
        SurfaceGeometry::Cone { .. } => {
            // Cones can be fillet-derived (offset of a planar edge
            // adjacent to a cone), but support-face detection is
            // ambiguous. Defer to a future round.
            let _ = perp_count;
            None
        }
        SurfaceGeometry::Gregory(_) | SurfaceGeometry::PolyhedralSpline(_) => {
            // A Gregory or polyhedral-net-spline patch at a face
            // implies a §5.3 N-sided vertex blend — the corner patch
            // constructed when multiple fillets meet. High confidence
            // because the cap surface type itself is the signal.
            //
            // This is exactly what the PnS variant was added for: as
            // anonymous NURBS a committed cap sector was invisible
            // here, so `identify_blend` could not see the corner it
            // was looking at. The patch also carries its valence and
            // sector index, so a caller can tell which piece of which
            // cap a face is rather than inferring it from adjacency.
            Some(BlendInfo {
                face,
                kind: BlendKind::VertexBlend,
                radius: None,
                support_faces: tangent_support_pair(&boundary, tangent_tol)
                    .or_else(|| chamfer_support_pair(&boundary, &distinct_neighbours)),
                source: IdentificationSource::CurvatureAnalysis,
                confidence: 0.95,
                cap_provenance: match &surface_on_face.geometry {
                    SurfaceGeometry::PolyhedralSpline(patch) => Some(patch.provenance),
                    _ => None,
                },
            })
        }
    }
}

/// Identify all blends in a [`BRep`].
///
/// Walks every face, calling [`identify_blend`] on each. Returns
/// one entry per recognised blend face. The order is the face-map
/// iteration order.
///
/// # Performance
///
/// Cost is O(F × E) where F is face count and E is average boundary
/// edges per face — the per-face walk does N projection calls. For
/// typical models this is fast enough for interactive use (sub-ms
/// on a unit cube). For heavy models, cache the result.
pub fn identify_all_blends(brep: &BRep) -> Vec<BlendInfo> {
    brep.faces
        .keys()
        .filter_map(|fh| identify_blend(brep, fh))
        .collect()
}

/// Convenience wrapper: return the nominal radius of the blend on
/// `face`, or `None` if the face is not recognised as a blend or
/// if the radius could not be extracted.
pub fn blend_radius(brep: &BRep, face: FaceHandle) -> Option<f64> {
    identify_blend(brep, face).and_then(|info| info.radius)
}

/// Convenience wrapper: return the support-face pair behind the
/// blend on `face`, or `None` if the face is not a blend or the
/// supports could not be determined.
pub fn blend_support_faces(brep: &BRep, face: FaceHandle) -> Option<(FaceHandle, FaceHandle)> {
    identify_blend(brep, face).and_then(|info| info.support_faces)
}

/// Return every distinct neighbour face of `face`, in outer-loop
/// traversal order. Useful for inspecting vertex blends that may
/// touch more than two adjacent fillets.
pub fn blend_neighbours(brep: &BRep, face: FaceHandle) -> Vec<FaceHandle> {
    let Some(boundary) = boundary_tangency(brep, face) else {
        return Vec::new();
    };
    distinct_neighbours(&boundary)
}

// Internal helpers live in `identify_helpers` — kept separate so
// each source file stays under the 500-line structural limit.

#[cfg(test)]
#[path = "identify_tests.rs"]
mod tests;
