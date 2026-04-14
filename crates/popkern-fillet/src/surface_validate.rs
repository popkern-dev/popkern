//! Post-construction self-intersection detection for fillet output
//! surfaces (§13.2 "Self-Intersection Detection").
//!
//! Called by `compute_fillet_geometry` after `build_fillet_surface`
//! returns, before the topology stitch-up phase. Rejects surfaces
//! that would leave the output BRep in an invalid state:
//!
//! - **Spindle torus** (`minor_radius > major_radius`) — inner tube
//!   passes through itself. Caught at the primitive level.
//! - **Twisted NURBS sweeps** — cross-section rotation > π causes
//!   adjacent UV samples to have opposing normals. Caught by a
//!   sign-flip check on the normal field.
//! - **Self-intersecting NURBS sweeps** — two non-neighbouring UV
//!   samples land at the same 3D point (within a multiple of the
//!   geometric tolerance). Caught by an all-pairs 3D-distance scan
//!   on a small UV grid.
//!
//! The §13.2 false-positive acceptance is < 1 %. Thresholds here
//! are deliberately conservative — the detector only fires on
//! clear failures, and valid fillet surfaces that happen to have
//! tight concavities are not flagged. Any threshold-tuning that
//! increases sensitivity must be accompanied by regression runs
//! of the full fillet test suite.
//!
//! The pre-construction undercut check in
//! [`crate::offset::self_intersect`] catches most of the "radius
//! exceeds local curvature" cases before the output surface is ever
//! built; this module is the safety net that catches the rest and
//! the specifically-analytic spindle-torus case that the pre-check
//! doesn't model.

use popkern_core::geometry::Point3;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::tolerance::Tolerance;

/// Outcome of [`validate_fillet_surface`] on a single surface.
#[derive(Debug, Clone)]
pub enum SurfaceValidationIssue {
    /// A torus with `minor_radius > major_radius` is a spindle torus
    /// whose inner tube self-intersects. Produced by the primitive-
    /// level check on `SurfaceGeometry::Torus`.
    SpindleTorus {
        /// The offending major radius.
        major_radius: f64,
        /// The offending minor radius.
        minor_radius: f64,
    },
    /// Two non-neighbouring UV samples on a NURBS surface landed at
    /// the same 3D point (within threshold). Typically indicates the
    /// rolling-ball sweep has wrapped back on itself.
    SelfIntersecting {
        /// First UV parameter of the coincident pair.
        uv_a: (f64, f64),
        /// Second UV parameter of the coincident pair.
        uv_b: (f64, f64),
        /// 3D distance between the two samples.
        distance: f64,
    },
    /// A NURBS surface's normal field flipped sign between adjacent
    /// UV samples. Indicates a twist exceeding π in the cross-section
    /// rotation.
    Twisted {
        /// UV parameter of one adjacent sample.
        uv_a: (f64, f64),
        /// UV parameter of the twisted neighbour.
        uv_b: (f64, f64),
    },
}

impl SurfaceValidationIssue {
    /// Produce a short, user-facing description of the issue.
    /// Used as the `reason` field on `BlendFailure::SelfIntersectingBlend`.
    pub fn describe(&self) -> String {
        match self {
            SurfaceValidationIssue::SpindleTorus {
                major_radius,
                minor_radius,
            } => format!(
                "spindle torus detected: minor_radius {minor_radius:.6} > \
                 major_radius {major_radius:.6}; the resulting surface \
                 self-intersects along the inner tube"
            ),
            SurfaceValidationIssue::SelfIntersecting {
                uv_a,
                uv_b,
                distance,
            } => format!(
                "NURBS surface self-intersects: samples at ({:.3}, {:.3}) \
                 and ({:.3}, {:.3}) coincide in 3D within {:.6}",
                uv_a.0, uv_a.1, uv_b.0, uv_b.1, distance
            ),
            SurfaceValidationIssue::Twisted { uv_a, uv_b } => format!(
                "NURBS surface twisted: normal field flipped between \
                 adjacent samples at ({:.3}, {:.3}) and ({:.3}, {:.3})",
                uv_a.0, uv_a.1, uv_b.0, uv_b.1
            ),
        }
    }
}

/// Check a fillet output surface for self-intersection, twist, and
/// other post-construction §13.2 failure modes.
///
/// Returns `Ok(())` when the surface is clean, `Err(issue)` with a
/// structured description on the first detected problem. Early-
/// returns on the first issue — does not enumerate all problems.
pub fn validate_fillet_surface(
    surface: &SurfaceGeometry,
    u_range: (f64, f64),
    v_range: (f64, f64),
    tolerance: &Tolerance,
) -> Result<(), SurfaceValidationIssue> {
    match surface {
        SurfaceGeometry::Torus {
            major_radius,
            minor_radius,
            ..
        } => {
            if *minor_radius > *major_radius {
                return Err(SurfaceValidationIssue::SpindleTorus {
                    major_radius: *major_radius,
                    minor_radius: *minor_radius,
                });
            }
            Ok(())
        }
        // Convex analytic primitives (infinite planes, cylinders,
        // spheres, cones) cannot self-intersect by construction.
        SurfaceGeometry::Plane { .. }
        | SurfaceGeometry::Cylinder { .. }
        | SurfaceGeometry::Sphere { .. }
        | SurfaceGeometry::Cone { .. } => Ok(()),
        SurfaceGeometry::Nurbs(_) => check_nurbs_surface(surface, u_range, v_range, tolerance),
    }
}

/// Run the grid-sampling checks on a NURBS fillet surface.
fn check_nurbs_surface(
    surface: &SurfaceGeometry,
    u_range: (f64, f64),
    v_range: (f64, f64),
    tolerance: &Tolerance,
) -> Result<(), SurfaceValidationIssue> {
    // 6×6 grid = 36 samples. Dense enough to catch clear failures
    // without being expensive on every fillet call. If the grid
    // size is tuned up for tighter detection, re-measure the
    // false-positive rate against the full fillet test suite.
    const GRID_N: usize = 6;
    let samples = sample_uv_grid(surface, u_range, v_range, GRID_N);

    check_twist(surface, u_range, v_range, GRID_N)?;
    check_self_intersection(&samples, GRID_N, tolerance)?;

    Ok(())
}

/// Sample the surface on a regular `n × n` UV grid.
fn sample_uv_grid(
    surface: &SurfaceGeometry,
    (u_lo, u_hi): (f64, f64),
    (v_lo, v_hi): (f64, f64),
    n: usize,
) -> Vec<(f64, f64, Point3)> {
    let mut out = Vec::with_capacity(n * n);
    for i in 0..n {
        let u = u_lo + (u_hi - u_lo) * (i as f64) / ((n - 1).max(1) as f64);
        for j in 0..n {
            let v = v_lo + (v_hi - v_lo) * (j as f64) / ((n - 1).max(1) as f64);
            out.push((u, v, surface.point_at(u, v)));
        }
    }
    out
}

/// Detect a normal-field sign flip between adjacent UV samples.
///
/// Walks the UV grid, computing the surface normal at each sample.
/// Between any two neighbours (u-step or v-step), if the normals have
/// `n_i · n_j < −0.1` the surface has twisted by more than roughly
/// 95°. Valid fillet surfaces (even degree-5 G2 cross-sections)
/// never twist by this much between a pair of neighbouring samples
/// on a 6×6 grid, so the threshold is well below the false-positive
/// floor.
fn check_twist(
    surface: &SurfaceGeometry,
    (u_lo, u_hi): (f64, f64),
    (v_lo, v_hi): (f64, f64),
    n: usize,
) -> Result<(), SurfaceValidationIssue> {
    let denom = (n - 1).max(1) as f64;
    let u_at = |i: usize| u_lo + (u_hi - u_lo) * (i as f64) / denom;
    let v_at = |j: usize| v_lo + (v_hi - v_lo) * (j as f64) / denom;

    for i in 0..n {
        for j in 0..n {
            let u0 = u_at(i);
            let v0 = v_at(j);
            let n0 = surface.normal_at(u0, v0);
            let n0_len = n0.length();
            if n0_len < 1e-12 {
                continue; // degenerate point; skip silently
            }
            let n0_unit = n0 * (1.0 / n0_len);

            // u-neighbour
            if i + 1 < n {
                let u1 = u_at(i + 1);
                let n1 = surface.normal_at(u1, v0);
                let n1_len = n1.length();
                if n1_len >= 1e-12 {
                    let n1_unit = n1 * (1.0 / n1_len);
                    if n0_unit.dot(&n1_unit) < -0.1 {
                        return Err(SurfaceValidationIssue::Twisted {
                            uv_a: (u0, v0),
                            uv_b: (u1, v0),
                        });
                    }
                }
            }
            // v-neighbour
            if j + 1 < n {
                let v1 = v_at(j + 1);
                let n1 = surface.normal_at(u0, v1);
                let n1_len = n1.length();
                if n1_len >= 1e-12 {
                    let n1_unit = n1 * (1.0 / n1_len);
                    if n0_unit.dot(&n1_unit) < -0.1 {
                        return Err(SurfaceValidationIssue::Twisted {
                            uv_a: (u0, v0),
                            uv_b: (u0, v1),
                        });
                    }
                }
            }
        }
    }
    Ok(())
}

/// Detect coincident 3D points from non-neighbouring UV samples.
///
/// All-pairs O(n^4) scan on an `n × n` grid (36 samples → 630 pairs).
/// For each pair, require:
///
/// 1. **UV distance > 1/3 of the diagonal** — excludes immediate
///    neighbours whose 3D separation is naturally small on dense
///    samplings of thin fillet surfaces.
/// 2. **3D distance < 5 × `tolerance.linear()`** — flags coincidence
///    well above numerical noise but below any separation a valid
///    fillet surface would have across a 1/3-range UV gap.
///
/// Both conditions must hold. The two-gate design is what keeps
/// the false-positive rate below the §13.2 1 % acceptance: tight
/// but valid concavities satisfy at most one of the conditions.
fn check_self_intersection(
    samples: &[(f64, f64, Point3)],
    n: usize,
    tolerance: &Tolerance,
) -> Result<(), SurfaceValidationIssue> {
    let uv_gap_threshold = (n as f64 - 1.0) / 3.0;
    let distance_threshold = tolerance.linear() * 5.0;

    let index = |i: usize, j: usize| i * n + j;
    for i1 in 0..n {
        for j1 in 0..n {
            let (u1, v1, p1) = samples[index(i1, j1)];
            for i2 in 0..n {
                for j2 in 0..n {
                    if (i2, j2) <= (i1, j1) {
                        continue; // only unordered pairs
                    }
                    // UV grid distance, in integer-cell units.
                    let du = (i2 as i64 - i1 as i64).unsigned_abs() as f64;
                    let dv = (j2 as i64 - j1 as i64).unsigned_abs() as f64;
                    let uv_gap = (du * du + dv * dv).sqrt();
                    if uv_gap < uv_gap_threshold {
                        continue;
                    }
                    let (u2, v2, p2) = samples[index(i2, j2)];
                    let distance = p1.distance_to(&p2);
                    if distance < distance_threshold {
                        return Err(SurfaceValidationIssue::SelfIntersecting {
                            uv_a: (u1, v1),
                            uv_b: (u2, v2),
                            distance,
                        });
                    }
                }
            }
        }
    }
    Ok(())
}

#[cfg(test)]
#[path = "surface_validate_tests.rs"]
mod tests;
