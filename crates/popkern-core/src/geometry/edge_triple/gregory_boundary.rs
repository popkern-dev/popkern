//! Gregory-boundary dispatch (ADR §5, stage 5).
//!
//! When the caller knows an edge corresponds to a canonical Gregory
//! UV boundary (passing
//! [`EdgeHint::SurfaceBoundary { side, boundary }`]), this dispatcher
//! constructs all three pieces of edge geometry from one source — the
//! Gregory's symbolic boundary curve. That eliminates the cap-to-
//! fillet seam drift documented in Task #144 / journal entry
//! "Commercial-grade at N=3,4 claim was hidden-gap-coverage".
//!
//! # Construction
//!
//! 1. **`curve_3d`** = the Gregory's UV-boundary curve as a NURBS,
//!    extracted symbolically from `boundary_cps`. Exact, not a fit.
//! 2. **`pcurve` on the Gregory side** = the canonical UV-axis line
//!    (e.g. `(0, t)` for the `u=0` boundary). Hits the boundary by
//!    construction: `gregory(0, t) ≡ curve_3d(t)` for all `t`.
//! 3. **`pcurve` on the other side** = projection of `curve_3d` onto
//!    the other surface's UV. For plane/cylinder/sphere/torus, the
//!    projection is closed-form via [`project_point_to_uv`].
//!
//! All three reference the same canonical parameter `t`. G3
//! (parameter map) is satisfied as `Identity`.
//!
//! # Coverage (stage 5)
//!
//! Implemented:
//! - Gregory × Cylinder
//! - Gregory × Plane
//! - Gregory × Sphere
//! - Gregory × Torus
//!
//! Not yet:
//! - Gregory × Cone (rare in fillet pipeline; defer)
//! - Gregory × NURBS (marching; defer to stage 7 catch-all)
//! - Gregory × Gregory (corner-of-corner; mention in comments only —
//!   needs both boundaries simultaneously)

use crate::geometry::curve::{CurveEval, CurveGeometry};
use crate::geometry::edge_triple::error::EdgeTripleError;
use crate::geometry::edge_triple::inputs::{BoundarySide, EdgeAnchors, EdgeHint, SurfaceSide};
use crate::geometry::edge_triple::types::{
    EdgeTopology, EdgeTriple, EdgeTripleProvenance, EdgeTripleQuality, GregorySide, ParameterMap,
};
use crate::geometry::gregory::GregoryPatch;
use crate::geometry::nurbs_curve::NurbsCurve;
use crate::geometry::pcurve::Pcurve;
use crate::geometry::surface::SurfaceGeometry;
use crate::geometry::uv_projection::{project_point_to_uv, project_point_to_uv_near};
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;

#[path = "gregory_boundary_trim.rs"]
mod gregory_boundary_trim;
use gregory_boundary_trim::orient_boundary_to_anchors;

/// Sampling density for projecting `curve_3d` onto the non-Gregory
/// surface's UV. Higher than the legacy `reproject_loop_pcurves`
/// fitter (12 samples) — the cap-boundary edge is the user-visible
/// seam, so we accept slightly more compute to push G2 below 1e-6.
const PROJECTION_SAMPLES: usize = 32;

/// Dispatcher entry. Returns:
/// * `Some(Ok(triple))` — Gregory boundary applies and was constructed
/// * `Some(Err(e))` — Gregory boundary applies but failed (e.g. anchor off surface)
/// * `None` — not a Gregory-boundary case; let the next dispatcher try
pub(super) fn try_dispatch_gregory_boundary(
    surf_a: &SurfaceGeometry,
    surf_b: &SurfaceGeometry,
    anchors: &EdgeAnchors,
    hint: Option<&EdgeHint>,
    tol: &Tolerance,
) -> Option<Result<EdgeTriple, EdgeTripleError>> {
    let (cap_side, boundary_side) = match hint? {
        EdgeHint::SurfaceBoundary { side, boundary } => (*side, *boundary),
        _ => return None,
    };
    let (cap_surf_geo, other_surf_geo, gregory_side, cap_is_a) = match cap_side {
        SurfaceSide::A => match surf_a {
            SurfaceGeometry::Gregory(g) => (g, surf_b, GregorySide::A, true),
            _ => return None,
        },
        SurfaceSide::B => match surf_b {
            SurfaceGeometry::Gregory(g) => (g, surf_a, GregorySide::B, false),
            _ => return None,
        },
    };

    let boundary_curve = match extract_gregory_boundary_curve(cap_surf_geo, boundary_side) {
        Some(nc) => nc,
        None => {
            return Some(Err(EdgeTripleError::DegenerateSurfacePair {
                reason: "gregory boundary extraction not implemented for this patch variant",
            }));
        }
    };
    let cap_pcurve_uv = canonical_uv_line(boundary_side);

    // Project curve_3d onto the other surface's UV. For analytic
    // surfaces this is closed-form; for NURBS we leave it to a
    // later stage.
    let other_pcurve = match build_other_pcurve(&boundary_curve, other_surf_geo) {
        Some(pc) => pc,
        None => {
            // Other surface is non-analytic (NURBS or Gregory). Caller
            // falls through to a later dispatcher (marching catch-all).
            return None;
        }
    };

    let (curve_3d_oriented, cap_pcurve_oriented, other_pcurve_oriented, max_endpoint_err) =
        match orient_boundary_to_anchors(
            boundary_curve,
            cap_pcurve_uv,
            other_pcurve,
            anchors,
            cap_side,
            tol,
        )? {
            Ok(oriented) => oriented,
            Err(err) => return Some(Err(err)),
        };

    let (curve_t_start, curve_t_end) = curve_3d_oriented.parameter_range();

    // Map cap/other roles back to (a, b) per call argument order.
    let (pcurve_a, pcurve_b) = if cap_is_a {
        (cap_pcurve_oriented, other_pcurve_oriented)
    } else {
        (other_pcurve_oriented, cap_pcurve_oriented)
    };

    let triple = EdgeTriple {
        curve_3d: CurveGeometry::Nurbs(curve_3d_oriented),
        pcurve_a: Pcurve::from(pcurve_a),
        pcurve_b: Pcurve::from(pcurve_b),
        t_start: curve_t_start,
        t_end: curve_t_end,
        topology: EdgeTopology::Open,
        parameter_map: ParameterMap::Identity,
        provenance: EdgeTripleProvenance::GregoryBoundary {
            side: gregory_side,
            boundary_index: boundary_side_to_index(boundary_side),
        },
        quality: EdgeTripleQuality {
            max_positional_residual: 0.0,
            max_endpoint_error: max_endpoint_err,
            samples: 0,
            tolerance_used: tol.linear(),
        },
    };
    Some(Ok(triple))
}

/// Extract the Gregory's UV boundary as a NurbsCurve. Stage 5
/// supports the rectangular variants used by the fillet pipeline
/// (Rectangular = cubic, Deg5 = quintic, BiquinticFull = quintic
/// grid, Biseptic = septic, TwoStripe = bi-9 with two-stripe layout,
/// BiNonic = bi-9). Triangular variants return `None`.
pub(crate) fn extract_gregory_boundary_curve(
    patch: &GregoryPatch,
    boundary: BoundarySide,
) -> Option<NurbsCurve> {
    match patch {
        GregoryPatch::Rectangular(r) => {
            let cps = boundary_cps_4x4(&r.boundary_cps, boundary);
            Some(bezier_to_nurbs(&cps, 3))
        }
        GregoryPatch::RectangularDeg5(r) => {
            let cps = boundary_cps_6x6(&r.boundary_cps, boundary);
            Some(bezier_to_nurbs(&cps, 5))
        }
        GregoryPatch::RectangularBiquinticFull(r) => {
            let cps = boundary_cps_6x6(&r.grid, boundary);
            Some(bezier_to_nurbs(&cps, 5))
        }
        GregoryPatch::RectangularBiseptic(r) => {
            let cps = boundary_cps_8x8(&r.grid, boundary);
            Some(bezier_to_nurbs(&cps, 7))
        }
        GregoryPatch::RectangularTwoStripe(r) => {
            let cps = boundary_cps_8x8(&r.grid, boundary);
            Some(bezier_to_nurbs(&cps, 7))
        }
        GregoryPatch::RectangularBiNonic(r) => {
            let cps = boundary_cps_10x10(&r.grid, boundary);
            Some(bezier_to_nurbs(&cps, 9))
        }
        // Boundary rows/columns carry plain controls — the
        // rational dual-twist blend affects interior CPs only.
        GregoryPatch::RectangularBiNonicDualTwist(r) => {
            let cps = boundary_cps_10x10(&r.grid, boundary);
            Some(bezier_to_nurbs(&cps, 9))
        }
        GregoryPatch::Triangular(_) => None,
    }
}

/// Map ADR §1's [`BoundarySide`] to a grid row/column slice.
/// Mirrors the boundary-orientation convention in
/// `popkern_core::geometry::gregory::build`:
/// `boundary_cps[i][j]` is at Bezier index `(i, j)` where `i` is
/// the u-index and `j` is the v-index. Therefore:
/// * `u=0` boundary = `grid[0][0..N]` (first row)
/// * `u=1` boundary = `grid[N-1][0..N]` (last row)
/// * `v=0` boundary = `grid[0..N][0]` (first column)
/// * `v=1` boundary = `grid[0..N][N-1]` (last column)
fn boundary_cps_4x4(grid: &[[Point3; 4]; 4], boundary: BoundarySide) -> Vec<Point3> {
    match boundary {
        BoundarySide::UMin => grid[0].to_vec(),
        BoundarySide::UMax => grid[3].to_vec(),
        BoundarySide::VMin => grid.iter().map(|row| row[0]).collect(),
        BoundarySide::VMax => grid.iter().map(|row| row[3]).collect(),
        BoundarySide::Indexed(_) => Vec::new(),
    }
}

fn boundary_cps_6x6(grid: &[[Point3; 6]; 6], boundary: BoundarySide) -> Vec<Point3> {
    match boundary {
        BoundarySide::UMin => grid[0].to_vec(),
        BoundarySide::UMax => grid[5].to_vec(),
        BoundarySide::VMin => grid.iter().map(|row| row[0]).collect(),
        BoundarySide::VMax => grid.iter().map(|row| row[5]).collect(),
        BoundarySide::Indexed(_) => Vec::new(),
    }
}

fn boundary_cps_8x8(grid: &[[Point3; 8]; 8], boundary: BoundarySide) -> Vec<Point3> {
    match boundary {
        BoundarySide::UMin => grid[0].to_vec(),
        BoundarySide::UMax => grid[7].to_vec(),
        BoundarySide::VMin => grid.iter().map(|row| row[0]).collect(),
        BoundarySide::VMax => grid.iter().map(|row| row[7]).collect(),
        BoundarySide::Indexed(_) => Vec::new(),
    }
}

fn boundary_cps_10x10(grid: &[[Point3; 10]; 10], boundary: BoundarySide) -> Vec<Point3> {
    match boundary {
        BoundarySide::UMin => grid[0].to_vec(),
        BoundarySide::UMax => grid[9].to_vec(),
        BoundarySide::VMin => grid.iter().map(|row| row[0]).collect(),
        BoundarySide::VMax => grid.iter().map(|row| row[9]).collect(),
        BoundarySide::Indexed(_) => Vec::new(),
    }
}

/// Build a NurbsCurve from `n+1` Bezier CPs of degree `degree` over
/// the parameter range `[0, 1]`. Knot vector is `[0]*(degree+1) ++
/// [1]*(degree+1)`, weights all 1. CP ordering matches input order.
fn bezier_to_nurbs(cps: &[Point3], degree: usize) -> NurbsCurve {
    debug_assert_eq!(cps.len(), degree + 1);
    let knots: Vec<f64> = std::iter::repeat_n(0.0, degree + 1)
        .chain(std::iter::repeat_n(1.0, degree + 1))
        .collect();
    let weights: Vec<f64> = vec![1.0; cps.len()];
    NurbsCurve {
        degree,
        control_points: cps.to_vec(),
        weights,
        knots,
    }
}

/// Canonical UV line for the given Gregory boundary as a degree-1
/// NurbsCurve over `[0, 1]`. Exact by construction:
/// `gregory(line(t)) ≡ boundary_curve(t)` for `t ∈ [0, 1]`.
fn canonical_uv_line(boundary: BoundarySide) -> NurbsCurve {
    let (a, b) = match boundary {
        BoundarySide::UMin => (Point3::new(0.0, 0.0, 0.0), Point3::new(0.0, 1.0, 0.0)),
        BoundarySide::UMax => (Point3::new(1.0, 0.0, 0.0), Point3::new(1.0, 1.0, 0.0)),
        BoundarySide::VMin => (Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)),
        BoundarySide::VMax => (Point3::new(0.0, 1.0, 0.0), Point3::new(1.0, 1.0, 0.0)),
        BoundarySide::Indexed(_) => (Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)),
    };
    NurbsCurve {
        degree: 1,
        control_points: vec![a, b],
        weights: vec![1.0, 1.0],
        knots: vec![0.0, 0.0, 1.0, 1.0],
    }
}

/// Project a 3-D NurbsCurve onto another surface's UV space by
/// sampling, projecting each sample, and fitting a NurbsCurve through
/// the (u, v, 0) samples.
///
/// Analytic surfaces use their closed-form inverse. NURBS surfaces use
/// continuation from the previous UV sample, which avoids branch jumping
/// along a boundary curve. Gregory projection currently uses the generic
/// closest-point inverse per sample; the caller still validates the final
/// 3D curve against both surfaces before committing topology repair.
pub(super) fn build_other_pcurve(
    curve_3d: &NurbsCurve,
    other: &SurfaceGeometry,
) -> Option<NurbsCurve> {
    let (t_lo, t_hi) = curve_3d.parameter_range();
    let mut params = Vec::with_capacity(PROJECTION_SAMPLES);
    let mut uv_samples = Vec::with_capacity(PROJECTION_SAMPLES);
    let mut prev_uv: Option<(f64, f64)> = None;
    for i in 0..PROJECTION_SAMPLES {
        let s = (i as f64) / ((PROJECTION_SAMPLES - 1) as f64);
        let t = t_lo + s * (t_hi - t_lo);
        let p = curve_3d.point_at(t);
        let (mut u, mut v) = if let Some((pu, pv)) = prev_uv {
            project_point_to_uv_near(other, &p, pu, pv)?
        } else {
            project_point_to_uv(other, &p)?
        };
        // Periodic-axis unwrapping for cyl/cone/torus u-axis.
        if let Some((pu, pv)) = prev_uv {
            if let (Some(period_u), _) = other.periodicity() {
                while u - pu > period_u * 0.5 {
                    u -= period_u;
                }
                while pu - u > period_u * 0.5 {
                    u += period_u;
                }
            }
            if let (_, Some(period_v)) = other.periodicity() {
                while v - pv > period_v * 0.5 {
                    v -= period_v;
                }
                while pv - v > period_v * 0.5 {
                    v += period_v;
                }
            }
        }
        prev_uv = Some((u, v));
        params.push(s);
        uv_samples.push(Point3::new(u, v, 0.0));
    }
    // Interpolate every projected UV sample. Boundary-owned shared
    // trims are committed through strict curve/surface validation, so
    // least-squares compression here can incorrectly decline an
    // otherwise valid exact-boundary candidate on curved neighbors.
    let degree = 5usize.min(uv_samples.len().saturating_sub(1));
    crate::geometry::fitting::interpolate_curve(&uv_samples, &params, degree).ok()
}

/// Reverse a NurbsCurve so `forward(t)` becomes `reverse(1 - t')`.
/// Used when the caller's anchors are in the opposite direction
/// from the boundary's natural parameter.
pub(super) fn reverse_nurbs(curve: &NurbsCurve) -> NurbsCurve {
    let n = curve.control_points.len();
    let mut cps = curve.control_points.clone();
    cps.reverse();
    let mut weights = curve.weights.clone();
    weights.reverse();
    let Some(&k_lo) = curve.knots.first() else {
        return curve.clone();
    };
    let Some(&k_hi) = curve.knots.last() else {
        return curve.clone();
    };
    let knots: Vec<f64> = curve
        .knots
        .iter()
        .rev()
        .map(|k| k_lo + (k_hi - k))
        .collect();
    let _ = n;
    NurbsCurve {
        degree: curve.degree,
        control_points: cps,
        weights,
        knots,
    }
}

pub(super) fn boundary_side_to_index(boundary: BoundarySide) -> u8 {
    match boundary {
        BoundarySide::VMin => 0,
        BoundarySide::UMax => 1,
        BoundarySide::VMax => 2,
        BoundarySide::UMin => 3,
        BoundarySide::Indexed(i) => i,
    }
}

// Suppress unused-warning for `Vec3` if cleanup removes the only
// reference; placeholder for future cone dispatch.
#[allow(dead_code)]
fn _vec3_placeholder(v: Vec3) -> Vec3 {
    v
}

#[cfg(test)]
#[path = "gregory_boundary_tests.rs"]
mod tests;
