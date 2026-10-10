//! Edges lying in a face's surface (PK-750), split from
//! `intersect_edge_face.rs` under the file-size budget.

use popkern_core::geometry::CurveEval;
use popkern_core::geometry::curve::CurveGeometry;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::tolerance::Tolerance;

/// Samples along the edge's own window.
const SAMPLES: usize = 9;

/// Whether the edge curve over `t_range` lies on `surface`: every sample
/// (ends included) within LINEAR tolerance of its closest surface point.
/// Coincidence is the model tolerance's decision; the 10 × linear
/// honesty band is too loose for it (measured: at a coarse 1e-2
/// tolerance every unit-cube edge stays within 0.1 of an r = 0.8 sphere
/// it crosses, and skipping them dropped every crossing).
///
/// Such an edge (a rim circle in a coplanar cap) has distance
/// identically zero to the surface: any root the solver returns is an
/// arbitrary point of a coincidence, not a crossing, and the pair is
/// skipped. Measured (PK-750): a rod's half-rim in a flush box's bottom
/// face hit at a point 4.2e-5 off the true junction, the spurious node
/// chopped the rim record into zero-length pieces, and the union
/// declined on curveless edges. Coincident contact is the overlap
/// records' business.
pub(crate) fn edge_lies_on_surface(
    curve: &CurveGeometry,
    t_range: (f64, f64),
    surface: &SurfaceGeometry,
    tolerance: &Tolerance,
) -> bool {
    let band = tolerance.linear();
    (0..SAMPLES).all(|k| {
        let t = t_range.0 + (t_range.1 - t_range.0) * k as f64 / (SAMPLES - 1) as f64;
        let p = curve.point_at(t);
        surface
            .closest_point(&p)
            .is_some_and(|(u, v)| surface.point_at(u, v).distance_to(&p) <= band)
    })
}

/// Whether the hit at `t` is the contact at the nearer endpoint of the
/// window: the curve stays on `surface` (within the honesty band) all
/// the way from the hit to that endpoint. A curve TANGENT to the
/// surface at its own endpoint has a double root there, and a double
/// root's parameter is ill-conditioned: the solver's hit slides along
/// the curve by ~sqrt(solver residual) and passes a parameter endpoint
/// test (measured, PK-750: a rod's half-rim touching a flush box's
/// tangent face at its own vertex hit 1e-5 rad inside the window; the
/// node 4.2e-5 off the vertex chopped the rim into zero-length pieces).
pub(crate) fn touches_at_endpoint(
    curve: &CurveGeometry,
    (t_lo, t_hi): (f64, f64),
    t: f64,
    surface: &SurfaceGeometry,
    tolerance: &Tolerance,
) -> bool {
    let end = if t - t_lo <= t_hi - t { t_lo } else { t_hi };
    edge_lies_on_surface(curve, (t, end), surface, tolerance)
}

#[cfg(test)]
#[path = "intersect_edge_face_coincident_tests.rs"]
mod tests;
