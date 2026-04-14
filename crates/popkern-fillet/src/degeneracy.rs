//! §13.1 degeneracy detectors for fillet input edges.
//!
//! Runs before `compute_fillet_geometry` does any real work so that
//! degenerate input produces a structured diagnostic (`FilletIssue`
//! in the pre-flight path, `BlendFailure::InvalidTopology` with a
//! `degenerate: …` prefix in the pipeline path) rather than a
//! cryptic downstream failure.
//!
//! This module covers all six §13.1 sub-categories:
//!
//! - **(1) Cusp at edge endpoints** — edge curve tangent magnitude
//!   below `tolerance.linear()` at `t_start` or `t_end`.
//! - **(2) Zero-area support patches** — one of the adjacent face's
//!   polygonal outer-loop area below `tolerance.linear()²`.
//! - **(3) Nearly-tangent supports** — dihedral angle at the edge
//!   midpoint below `2 * tolerance.angular()`, meaning the two
//!   adjacent face normals are nearly anti-parallel (popkern face
//!   normals point into the solid, so anti-parallel == knife-edge).
//! - **(4) Parameterization singularity** — `point_at`, `partial_u`,
//!   or `partial_v` returns a non-finite value (NaN or ∞) at some
//!   point in the support surface's sample grid. Distinct from (6):
//!   zero-magnitude normals are handled there, while this catches
//!   genuine NaN/Inf from pathological NURBS control polygons or
//!   bad knot vectors.
//! - **(5) Collapsed edge** — 3D length between the start and end
//!   vertex positions below `tolerance.linear()`. Skipped for closed
//!   edges (`v_start == v_end`), which are the intentional closed-
//!   rim case, not degeneracy.
//! - **(6) Degenerate support normal** — the surface normal at the
//!   projected edge midpoint has magnitude below `tolerance.linear()`
//!   on one or both support faces.
//!
//! Detection short-circuits on the first issue found — a caller
//! receives one degeneracy at a time, which is both simpler to
//! report and consistent with the §16.1 "one structured error per
//! failure" pattern.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use crate::analysis::EdgeAnalysis;

/// Which endpoint of an edge a cusp was detected at.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EdgeEndpoint {
    /// The start vertex (parameter `t_start` on the edge curve).
    Start,
    /// The end vertex (parameter `t_end` on the edge curve).
    End,
}

/// Which support face a surface-level degeneracy was detected on.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SupportSide {
    /// `analysis.surface_1`, the first adjacent face.
    Face1,
    /// `analysis.surface_2`, the second adjacent face.
    Face2,
}

/// A structured §13.1 degeneracy diagnostic. Each variant names the
/// sub-category and carries the measured value so callers can render
/// an actionable message without string parsing.
#[derive(Debug, Clone)]
pub enum DegeneracyKind {
    /// §13.1(5): the edge is shorter than `tolerance.linear()`.
    CollapsedEdge {
        /// The measured 3D length between the start and end vertices.
        length: f64,
    },
    /// §13.1(1): the edge curve's tangent vanishes at one endpoint.
    DegenerateEdgeTangent {
        /// Which endpoint exhibits the vanishing tangent.
        at: EdgeEndpoint,
        /// The measured tangent magnitude.
        magnitude: f64,
    },
    /// §13.1(3): the two support faces are nearly tangent at the
    /// edge midpoint — their inward-pointing normals are nearly
    /// anti-parallel, meaning the solid has a knife-edge there.
    NearlyTangentSupports {
        /// The measured dihedral angle in radians. Smaller means
        /// more tangent; 0 is a perfect knife-edge.
        angle_rad: f64,
    },
    /// §13.1(6): one of the support face normals has magnitude
    /// below `tolerance.linear()` at the projected edge midpoint.
    DegenerateSupportNormal {
        /// Which support face exhibits the degenerate normal.
        side: SupportSide,
        /// The measured normal magnitude.
        magnitude: f64,
    },
    /// §13.1(2): a support face has polygonal outer-loop area
    /// below `tolerance.linear()²`. Indicates a face whose vertices
    /// have collapsed to a line or point — filleting an edge whose
    /// support is zero-area produces no meaningful geometry.
    ZeroAreaSupport {
        /// Which support face is zero-area.
        side: SupportSide,
        /// The computed polygonal area.
        area: f64,
    },
    /// §13.1(4): evaluating the support surface produces a non-
    /// finite (NaN or ∞) coordinate at some sampled (u, v).
    /// Typically indicates a pathological NURBS control polygon
    /// or a bad knot vector. Distinct from `DegenerateSupportNormal`:
    /// that catches zero-magnitude normals (well-defined point but
    /// undefined tangent plane), while this catches the point
    /// itself being non-evaluable.
    ParameterizationSingularity {
        /// Which support face exhibits the singularity.
        side: SupportSide,
        /// The (u, v) parameter where the non-finite value was
        /// first observed during the sample sweep.
        uv: (f64, f64),
    },
}

impl DegeneracyKind {
    /// Render a short, user-facing description suitable for a
    /// `BlendFailure::InvalidTopology.reason` or a `FilletIssue` text.
    pub fn describe(&self) -> String {
        match self {
            DegeneracyKind::CollapsedEdge { length } => format!(
                "degenerate: collapsed edge (length {length:.6e} below \
                 geometric tolerance)"
            ),
            DegeneracyKind::DegenerateEdgeTangent { at, magnitude } => {
                let at_str = match at {
                    EdgeEndpoint::Start => "start",
                    EdgeEndpoint::End => "end",
                };
                format!(
                    "degenerate: cusp at edge {at_str} (tangent magnitude \
                     {magnitude:.6e} below tolerance)"
                )
            }
            DegeneracyKind::NearlyTangentSupports { angle_rad } => format!(
                "degenerate: nearly-tangent support faces (dihedral \
                 {angle_rad:.6e} rad — below 2·tol_angular)"
            ),
            DegeneracyKind::DegenerateSupportNormal { side, magnitude } => {
                let side_str = match side {
                    SupportSide::Face1 => "face_1",
                    SupportSide::Face2 => "face_2",
                };
                format!(
                    "degenerate: vanishing {side_str} normal at edge \
                     midpoint (magnitude {magnitude:.6e})"
                )
            }
            DegeneracyKind::ZeroAreaSupport { side, area } => {
                let side_str = match side {
                    SupportSide::Face1 => "face_1",
                    SupportSide::Face2 => "face_2",
                };
                format!(
                    "degenerate: zero-area {side_str} support (polygonal area \
                     {area:.6e} below tolerance²)"
                )
            }
            DegeneracyKind::ParameterizationSingularity { side, uv } => {
                let side_str = match side {
                    SupportSide::Face1 => "face_1",
                    SupportSide::Face2 => "face_2",
                };
                format!(
                    "degenerate: parameterization singularity on {side_str} at \
                     ({:.6}, {:.6}) (NaN or ∞ in surface evaluation)",
                    uv.0, uv.1
                )
            }
        }
    }
}

/// Run all §13.1 degeneracy detectors against an analysed edge.
///
/// Short-circuits on the first issue; returns `Ok(())` when all
/// six detectors pass. Intended to be called both from
/// `precheck::precheck_fillet` (as a non-destructive query) and
/// from `compute_fillet_geometry` (as a pipeline pre-flight).
///
/// Detector order is chosen so each check runs on input that has
/// already passed the cheaper structural gates — structure first
/// (collapsed / zero-area / cusp), then parameterisation
/// (parameterization singularity), then local surface properties
/// (degenerate normal / nearly-tangent).
pub fn check_edge_degeneracy(
    brep: &BRep,
    analysis: &EdgeAnalysis,
    tolerance: &Tolerance,
) -> Result<(), DegeneracyKind> {
    check_collapsed_edge(brep, analysis, tolerance)?;
    check_zero_area_support(brep, analysis, tolerance)?;
    check_edge_tangent_cusp(brep, analysis, tolerance)?;
    check_parameterization_singularity(analysis)?;
    check_degenerate_support_normal(brep, analysis, tolerance)?;
    check_nearly_tangent_supports(brep, analysis, tolerance)?;
    Ok(())
}

/// §13.1(5): collapsed edge detection.
fn check_collapsed_edge(
    brep: &BRep,
    analysis: &EdgeAnalysis,
    tolerance: &Tolerance,
) -> Result<(), DegeneracyKind> {
    // Closed edges (v_start == v_end) are the intentional closed-
    // rim case — skip the length check for them.
    if analysis.vertex_start == analysis.vertex_end {
        return Ok(());
    }
    let p_start = brep.vertices[analysis.vertex_start].position;
    let p_end = brep.vertices[analysis.vertex_end].position;
    let length = (p_end - p_start).length();
    if length < tolerance.linear() {
        return Err(DegeneracyKind::CollapsedEdge { length });
    }
    Ok(())
}

/// §13.1(1): vanishing edge tangent at start or end.
fn check_edge_tangent_cusp(
    brep: &BRep,
    analysis: &EdgeAnalysis,
    tolerance: &Tolerance,
) -> Result<(), DegeneracyKind> {
    let Some(edge) = brep.edges.get(analysis.edge) else {
        return Ok(()); // Edge missing — caught by other checks.
    };
    let Some(curve_on_edge) = &edge.curve else {
        return Ok(()); // No curve to inspect — can't evaluate tangent.
    };
    let curve = &curve_on_edge.geometry;
    let tan_start = curve.tangent_at(curve_on_edge.t_start);
    let mag_start = tan_start.length();
    if mag_start < tolerance.linear() {
        return Err(DegeneracyKind::DegenerateEdgeTangent {
            at: EdgeEndpoint::Start,
            magnitude: mag_start,
        });
    }
    let tan_end = curve.tangent_at(curve_on_edge.t_end);
    let mag_end = tan_end.length();
    if mag_end < tolerance.linear() {
        return Err(DegeneracyKind::DegenerateEdgeTangent {
            at: EdgeEndpoint::End,
            magnitude: mag_end,
        });
    }
    Ok(())
}

/// §13.1(6): vanishing surface normal at the edge midpoint, on
/// either support face.
fn check_degenerate_support_normal(
    brep: &BRep,
    analysis: &EdgeAnalysis,
    tolerance: &Tolerance,
) -> Result<(), DegeneracyKind> {
    let midpoint = match super::geometry_compute::get_edge_midpoint(brep, analysis) {
        Ok(p) => p,
        Err(_) => return Ok(()), // Can't compute midpoint — let the caller handle it.
    };

    let check_side = |surface: &popkern_core::geometry::surface::SurfaceGeometry,
                      side: SupportSide|
     -> Result<(), DegeneracyKind> {
        let Some((u, v)) = surface.closest_point(&midpoint) else {
            return Ok(());
        };
        let normal = surface.normal_at(u, v);
        let mag = normal.length();
        if mag < tolerance.linear() {
            return Err(DegeneracyKind::DegenerateSupportNormal {
                side,
                magnitude: mag,
            });
        }
        Ok(())
    };

    check_side(&analysis.surface_1, SupportSide::Face1)?;
    check_side(&analysis.surface_2, SupportSide::Face2)?;
    Ok(())
}

/// §13.1(3): nearly-tangent supports — the two face normals are
/// anti-parallel within `2 * tolerance.angular()`.
///
/// popkern face normals point INTO the solid. A dihedral angle of
/// π (flat, no edge) corresponds to `n1 · n2 = 1` (parallel); a
/// dihedral angle of 0 (knife-edge) corresponds to `n1 · n2 = -1`
/// (anti-parallel). We flag the knife-edge case because it is what
/// makes a rolling-ball fillet impossible to place: as the dihedral
/// shrinks, the rolling-ball centre escapes to infinity.
///
/// The flat case (dihedral ≈ π, `n1 · n2 ≈ 1`) is handled by the
/// existing `Convexity::Mixed` classification in
/// `analysis::analyze_edges` and is not re-checked here to avoid
/// double-reporting.
fn check_nearly_tangent_supports(
    brep: &BRep,
    analysis: &EdgeAnalysis,
    tolerance: &Tolerance,
) -> Result<(), DegeneracyKind> {
    let midpoint = match super::geometry_compute::get_edge_midpoint(brep, analysis) {
        Ok(p) => p,
        Err(_) => return Ok(()),
    };
    let Some((u1, v1)) = analysis.surface_1.closest_point(&midpoint) else {
        return Ok(());
    };
    let Some((u2, v2)) = analysis.surface_2.closest_point(&midpoint) else {
        return Ok(());
    };
    let n1 = analysis.surface_1.normal_at(u1, v1);
    let n2 = analysis.surface_2.normal_at(u2, v2);
    let m1 = n1.length();
    let m2 = n2.length();
    if m1 < tolerance.linear() || m2 < tolerance.linear() {
        return Ok(()); // already caught upstream by check_degenerate_support_normal
    }
    let n1u = n1 * (1.0 / m1);
    let n2u = n2 * (1.0 / m2);
    let dot = n1u.dot(&n2u);

    // Apply popkern's orientation convention: if orientation is
    // `true`, the surface's natural normal equals the face normal;
    // otherwise the face normal is `-natural`. Face normals point
    // into the solid.
    let sign_1 = if analysis.orientation_1 { 1.0 } else { -1.0 };
    let sign_2 = if analysis.orientation_2 { 1.0 } else { -1.0 };
    let face_dot = sign_1 * sign_2 * dot;

    // Knife-edge: face_dot near -1 → dihedral near 0.
    // Threshold: angle between face normals > π - 2·tol_angular
    //   ⇔ face_dot < -cos(2·tol_angular) ≈ -(1 - 2·tol_angular²)
    let threshold = -(1.0 - 2.0 * tolerance.angular() * tolerance.angular());
    if face_dot < threshold {
        // Convert dot to angle for the diagnostic.
        let angle_between_face_normals = face_dot.clamp(-1.0, 1.0).acos();
        let dihedral = std::f64::consts::PI - angle_between_face_normals;
        return Err(DegeneracyKind::NearlyTangentSupports {
            angle_rad: dihedral,
        });
    }
    Ok(())
}

/// §13.1(2): zero-area polygonal outer loop on either support face.
///
/// Computes the magnitude of the sum of edge cross products around
/// the outer loop. For a planar polygon this is exactly twice the
/// polygonal area; for a non-planar loop it's a weighted measure
/// that still scales with the loop's overall extent. A loop whose
/// 3D cross-product sum is below `tolerance.linear()²` is degenerate
/// regardless of planarity.
///
/// Inner loops are not considered — a face's "degenerate or not"
/// classification is structural and comes from the outer loop
/// alone. A face with a meaningful outer loop but an inner loop
/// that spans most of it is a different problem (zero effective
/// area) that this detector does not try to catch.
fn check_zero_area_support(
    brep: &BRep,
    analysis: &EdgeAnalysis,
    tolerance: &Tolerance,
) -> Result<(), DegeneracyKind> {
    let tol_sq = tolerance.linear() * tolerance.linear();
    let check_side = |face: popkern_core::topology::FaceHandle,
                      side: SupportSide|
     -> Result<(), DegeneracyKind> {
        let positions = brep.face_outer_loop_positions(face);
        if positions.len() < 3 {
            return Err(DegeneracyKind::ZeroAreaSupport { side, area: 0.0 });
        }
        // Shoelace in 3D via sum of cross products from a reference
        // vertex. Uses positions[0] as the anchor; for each triangle
        // (p0, p_i, p_{i+1}) add the signed cross product.
        let p0 = positions[0];
        let mut total = popkern_core::geometry::Vec3::new(0.0, 0.0, 0.0);
        for i in 1..positions.len() - 1 {
            let v1 = positions[i] - p0;
            let v2 = positions[i + 1] - p0;
            total = total + v1.cross(&v2);
        }
        // |total| / 2 is the (vector-sum) polygonal area.
        let area = total.length() * 0.5;
        if area < tol_sq {
            return Err(DegeneracyKind::ZeroAreaSupport { side, area });
        }
        Ok(())
    };

    check_side(analysis.face_1, SupportSide::Face1)?;
    check_side(analysis.face_2, SupportSide::Face2)?;
    Ok(())
}

/// §13.1(4): NaN/∞ in surface evaluation at sampled (u, v) points.
///
/// Sweeps each support surface on an 8×8 UV grid over the natural
/// parameter range (finite dimensions only; infinite dimensions
/// clamp to `[-10, 10]` so the sweep terminates). At each sample
/// checks `point_at`, `partial_u`, and `partial_v` for non-finite
/// coordinates. Fires on the first non-finite sample found, so the
/// cost is bounded at 8×8×6 = 384 float checks in the clean case.
///
/// This detector catches genuinely broken surfaces — pathological
/// NURBS control polygons, bad knot vectors, arithmetic singularities
/// in derived formulas. popkern's analytic primitives return finite
/// values across their natural parameter range, so they will only
/// fire on the NURBS case.
fn check_parameterization_singularity(analysis: &EdgeAnalysis) -> Result<(), DegeneracyKind> {
    check_surface_singularity(&analysis.surface_1, SupportSide::Face1)?;
    check_surface_singularity(&analysis.surface_2, SupportSide::Face2)?;
    Ok(())
}

fn check_surface_singularity(
    surface: &popkern_core::geometry::surface::SurfaceGeometry,
    side: SupportSide,
) -> Result<(), DegeneracyKind> {
    let ((u_min, u_max), (v_min, v_max)) = surface.parameter_range();
    // Clamp infinite dimensions to a bounded sweep range so the grid
    // evaluation terminates. The specific bound is arbitrary — any
    // finite value works because the detector is looking for NaN/∞
    // from the surface evaluator, not from the range values.
    let u_lo = if u_min.is_finite() { u_min } else { -10.0 };
    let u_hi = if u_max.is_finite() { u_max } else { 10.0 };
    let v_lo = if v_min.is_finite() { v_min } else { -10.0 };
    let v_hi = if v_max.is_finite() { v_max } else { 10.0 };

    const N: usize = 8;
    let denom = (N - 1) as f64;
    for i in 0..N {
        let u = u_lo + (u_hi - u_lo) * (i as f64) / denom;
        for j in 0..N {
            let v = v_lo + (v_hi - v_lo) * (j as f64) / denom;
            let p = surface.point_at(u, v);
            let pu = surface.partial_u(u, v);
            let pv = surface.partial_v(u, v);
            if !p.x.is_finite()
                || !p.y.is_finite()
                || !p.z.is_finite()
                || !pu.x.is_finite()
                || !pu.y.is_finite()
                || !pu.z.is_finite()
                || !pv.x.is_finite()
                || !pv.y.is_finite()
                || !pv.z.is_finite()
            {
                return Err(DegeneracyKind::ParameterizationSingularity { side, uv: (u, v) });
            }
        }
    }
    Ok(())
}

#[cfg(test)]
#[path = "degeneracy_tests.rs"]
mod tests;
