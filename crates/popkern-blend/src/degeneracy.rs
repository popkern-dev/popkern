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
    use crate::degeneracy_surface::{
        check_degenerate_support_normal, check_nearly_tangent_supports,
        check_parameterization_singularity,
    };
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
///
/// Uses [`CurveEval::derivative_at(t, 1)`] rather than
/// [`CurveEval::tangent_at`]: popkern's `tangent_at` returns a
/// **normalized** tangent (or exactly zero when the derivative
/// magnitude is below `normalized()`'s floor), so checking its
/// length against a small tolerance only catches exact-zero
/// derivatives. The un-normalized first derivative preserves
/// magnitude information and is the correct quantity for "tangent
/// vanishing" detection.
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
    let d_start = curve.derivative_at(curve_on_edge.t_start, 1);
    let mag_start = d_start.length();
    if mag_start < tolerance.linear() {
        return Err(DegeneracyKind::DegenerateEdgeTangent {
            at: EdgeEndpoint::Start,
            magnitude: mag_start,
        });
    }
    let d_end = curve.derivative_at(curve_on_edge.t_end, 1);
    let mag_end = d_end.length();
    if mag_end < tolerance.linear() {
        return Err(DegeneracyKind::DegenerateEdgeTangent {
            at: EdgeEndpoint::End,
            magnitude: mag_end,
        });
    }
    Ok(())
}

// Surface-level detectors for §13.1(3), (4), (6) live in
// `degeneracy_surface.rs` — split out to keep this file under the
// 500-line structural cap. They are called from
// `check_edge_degeneracy` above.

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
        // Boundary polygon from EDGE-CURVE samples, not vertices: a
        // curved face bounded by closed curves (a cylinder wall
        // between two circles) has as few as two vertices, and the
        // vertex polygon is degenerate even though the face is not.
        let mut positions: Vec<popkern_core::geometry::Point3> = Vec::new();
        if let Some(face_data) = brep.faces.get(face) {
            let start = brep.loops[face_data.outer_loop].half_edge;
            let mut he = start;
            for _ in 0..10_000 {
                let data = &brep.half_edges[he];
                if let Some(curve) = brep.edges[data.edge].curve.as_ref() {
                    use popkern_core::geometry::curve::CurveEval;
                    let ceval: &dyn CurveEval = &curve.geometry;
                    let origin = brep.vertices[data.origin].position;
                    let (t0, t1) = (curve.t_start, curve.t_end);
                    let forward = ceval.point_at(t0).distance_to(&origin)
                        <= ceval.point_at(t1).distance_to(&origin);
                    for k in 0..8 {
                        let f = k as f64 / 8.0;
                        let t = if forward {
                            t0 + (t1 - t0) * f
                        } else {
                            t1 - (t1 - t0) * f
                        };
                        positions.push(ceval.point_at(t));
                    }
                } else {
                    positions.push(brep.vertices[data.origin].position);
                }
                he = data.next;
                if he == start {
                    break;
                }
            }
        }
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

#[cfg(test)]
#[path = "degeneracy_tests.rs"]
mod tests;
