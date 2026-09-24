//! Per-vertex surface-curvature fields over a tessellation (PK-520).
//!
//! Surface-quality evaluation for humans and gates alike: Gaussian
//! and mean curvature are the standard class-A inspection scalars
//! (a curvature-discontinuous seam shows as a value step across the
//! seam; an inflection band shows as a sign change), and the field is
//! the substrate for false-color rendering in popkern-viz and the
//! playground viewer.
//!
//! Curvature is evaluated ANALYTICALLY at each mesh vertex's recorded
//! surface parameters through [`SurfaceEval`] first and second
//! partials (first/second fundamental forms), never estimated from
//! the triangles — mesh density cannot alias the reported values.

use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::topology::{BRep, FaceHandle};
use popkern_mesh::TriMesh;

use crate::MeasureError;

/// Which curvature scalar a field carries.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CurvatureKind {
    /// Gaussian curvature `K = κ₁·κ₂` — sign separates elliptic,
    /// parabolic (0), and hyperbolic regions.
    Gaussian,
    /// Mean curvature `H = (κ₁ + κ₂) / 2` — the fillet inspection
    /// scalar (a constant-radius rolling-ball strip reads `±1/(2r)`
    /// plus the spine term; steps across seams are G2 defects).
    Mean,
    /// Largest absolute principal curvature `max(|κ₁|, |κ₂|)` —
    /// a radius-of-curvature check (`1/κ`).
    MaxPrincipalAbs,
}

/// Per-vertex curvature field over a [`TriMesh`], parallel to
/// `mesh.vertices`.
#[derive(Debug, Clone)]
pub struct CurvatureField {
    /// The scalar the field carries.
    pub kind: CurvatureKind,
    /// One value per mesh vertex. A vertex whose chart is too
    /// degenerate to evaluate (collapsed patch corner) or that no
    /// triangle references carries `0.0` and is counted in
    /// [`Self::unevaluated`].
    pub per_vertex: Vec<f64>,
    /// Smallest evaluated value.
    pub min: f64,
    /// Largest evaluated value.
    pub max: f64,
    /// Vertices that could not be evaluated (degenerate chart or
    /// unreferenced); their entries are `0.0`.
    pub unevaluated: usize,
}

/// Evaluate the curvature scalar `kind` at every vertex of `mesh`
/// against the analytic face surfaces of `brep`.
///
/// Each vertex is evaluated on the face whose triangles reference it
/// (mesh vertices are per-face; zipper-added vertices inherit the
/// face of their referencing triangle), at its recorded `(u, v)`.
pub fn curvature_field(
    brep: &BRep,
    mesh: &TriMesh,
    kind: CurvatureKind,
) -> Result<CurvatureField, MeasureError> {
    // Vertex → owning face via the triangle/face map.
    let mut vertex_face: Vec<Option<FaceHandle>> = vec![None; mesh.vertices.len()];
    for (tri, face) in mesh.triangles.iter().zip(mesh.face_map.iter()) {
        for &vi in tri {
            if vi < vertex_face.len() && vertex_face[vi].is_none() {
                vertex_face[vi] = Some(*face);
            }
        }
    }
    let mut per_vertex = vec![0.0_f64; mesh.vertices.len()];
    let mut unevaluated = 0usize;
    let (mut min, mut max) = (f64::INFINITY, f64::NEG_INFINITY);
    for (vi, vertex) in mesh.vertices.iter().enumerate() {
        let value = vertex_face[vi]
            .and_then(|face| brep.faces.get(face))
            .and_then(|f| f.surface.as_ref())
            .and_then(|s| curvature_at(&s.geometry, vertex.uv.0, vertex.uv.1, kind));
        match value {
            Some(v) => {
                per_vertex[vi] = v;
                min = min.min(v);
                max = max.max(v);
            }
            None => unevaluated += 1,
        }
    }
    if !min.is_finite() {
        (min, max) = (0.0, 0.0);
    }
    Ok(CurvatureField {
        kind,
        per_vertex,
        min,
        max,
        unevaluated,
    })
}

/// Analytic curvature scalar at surface parameters `(u, v)`, or
/// `None` where the chart is too degenerate to trust.
pub fn curvature_at(
    surface: &impl SurfaceEval,
    u: f64,
    v: f64,
    kind: CurvatureKind,
) -> Option<f64> {
    let su = surface.partial_u(u, v);
    let sv = surface.partial_v(u, v);
    let (lu, lv) = (su.length(), sv.length());
    if !lu.is_finite() || !lv.is_finite() || lu.min(lv) < 1e-9 * lu.max(lv).max(1e-12) {
        return None;
    }
    let n = su.cross(&sv).normalized()?;
    let suu = surface.second_partial_uu(u, v);
    let suv = surface.second_partial_uv(u, v);
    let svv = surface.second_partial_vv(u, v);
    let (e, f, g) = (su.dot(&su), su.dot(&sv), sv.dot(&sv));
    let (l, m, o) = (n.dot(&suu), n.dot(&suv), n.dot(&svv));
    let det = e * g - f * f;
    if !det.is_finite() || det.abs() < 1e-30 {
        return None;
    }
    let gaussian = (l * o - m * m) / det;
    let mean = (e * o + g * l - 2.0 * f * m) / (2.0 * det);
    let value = match kind {
        CurvatureKind::Gaussian => gaussian,
        CurvatureKind::Mean => mean,
        CurvatureKind::MaxPrincipalAbs => {
            // κ = H ± sqrt(H² − K); the discriminant is clamped at 0
            // (it is non-negative analytically; FD noise can dip).
            let disc = (mean * mean - gaussian).max(0.0).sqrt();
            (mean + disc).abs().max((mean - disc).abs())
        }
    };
    value.is_finite().then_some(value)
}

#[cfg(test)]
#[path = "surface_curvature_tests.rs"]
mod tests;
