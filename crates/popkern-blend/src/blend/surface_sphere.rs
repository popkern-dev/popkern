//! Rolling-ball sphere corner blend: least-squares sphere fit over
//! all boundary samples (PK-233), a radius-scaled G0 acceptance
//! window, and a G1 tangency certification against the incident
//! fillet strips — boundary-on-sphere alone is necessary but not
//! sufficient (a corner's arc+bridge ring can fit a non-tangent
//! sphere; see the displaced-ring impostor test below).
//!
//! Split from `blend/surface.rs` to honour the structural line cap.

use popkern_core::geometry::Point3;
use popkern_core::geometry::nurbs_curve::NurbsCurve;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::SurfaceEval;

use super::boundary::{BlendBoundary, BoundarySegment};
use super::surface::BlendSurfaceResult;
use super::surface_math::fit_sphere_least_squares;
use crate::fillet_surface::FilletSurfaceResult;
use crate::types::{ContinuityLevel, CornerMethod};

/// Angular gate for the tangency certification, as the sine of the
/// angle between the strip normal and the sphere's radial direction.
/// Exact rolling-ball corners measure ~1e-6 here and boundary-fit
/// noise (~100·tol on the sample positions) perturbs the angle by
/// ~1e-3 at typical radii; the measured non-tangent impostor (the
/// mitred-ring sphere coincidence on the uniform cube corner) sits
/// far above it.
const TANGENCY_SIN_GATE: f64 = 2.0e-2;

/// Certify that the fitted sphere is tangent (G1) to every incident
/// fillet strip along its FilletEnd boundary segments: the strip's
/// surface normal at each sample must be parallel to the sphere's
/// radial direction there, and the sample must sit at sphere-radius
/// distance from the center.
fn sphere_tangent_to_strips(
    boundary: &BlendBoundary,
    fillet_surfaces: &[FilletSurfaceResult],
    center: &Point3,
    radius: f64,
    trace: bool,
) -> bool {
    const SAMPLES_PER_SEGMENT: usize = 5;
    let mut certified_segments = 0usize;
    let mut worst_sin = 0.0f64;
    for segment in &boundary.segments {
        let BoundarySegment::FilletEnd {
            fillet_index,
            curve,
        } = segment
        else {
            continue;
        };
        let Some(strip) = fillet_surfaces.get(*fillet_index) else {
            return false;
        };
        let (t0, t1) = curve.parameter_range();
        for k in 0..SAMPLES_PER_SEGMENT {
            let t = t0 + (t1 - t0) * (k as f64 + 0.5) / SAMPLES_PER_SEGMENT as f64;
            let p = curve.point_at(t);
            let radial = p - *center;
            let dist = radial.length();
            if dist < radius * 1.0e-6 {
                return false;
            }
            let Some((u, v)) = strip.surface.closest_point(&p) else {
                return false;
            };
            let Some(n) = strip.surface.normal_at(u, v).normalized() else {
                return false;
            };
            let sin_angle = n.cross(&(radial * (1.0 / dist))).length();
            worst_sin = worst_sin.max(sin_angle);
            if sin_angle > TANGENCY_SIN_GATE {
                if trace {
                    eprintln!(
                        "[spherefit] tangency REJECT fillet={fillet_index} sin={sin_angle:.4e} gate={TANGENCY_SIN_GATE:.1e}"
                    );
                }
                return false;
            }
        }
        certified_segments += 1;
    }
    if trace {
        eprintln!(
            "[spherefit] tangency certified segments={certified_segments} worst_sin={worst_sin:.4e}"
        );
    }
    // A sphere that touched no FilletEnd segment certified nothing.
    certified_segments >= 2
}

/// Try to build a spherical cap blend.
///
/// Fits the least-squares sphere through ALL boundary samples
/// (PK-233; see `fit_sphere_least_squares`). If the points lie on
/// that sphere (max deviation within the scale-relative gate),
/// returns it. This handles cube corners with equal-radius
/// perpendicular fillets exactly, and — since the PK-233 rewrite —
/// regular pyramid apexes whose incident cylinder axes are
/// concurrent (the fit's predecessor solved a 3×3 built from
/// co-circular first-arc samples, which is singular even when the
/// full sample set determines the sphere exactly).
///
/// Returns `None` if the samples cannot determine a unique sphere
/// (coplanar/co-circular data), the fit is geometrically
/// implausible for a corner cap, the points don't lie on the
/// fitted sphere, or the sphere is not tangent to the incident
/// fillet strips (the G1 certification below).
pub(super) fn try_spherical_blend(
    boundary: &BlendBoundary,
    fillet_surfaces: &[FilletSurfaceResult],
    points: &[Point3],
    boundary_curves: &[NurbsCurve],
    tolerance: &Tolerance,
) -> Option<BlendSurfaceResult> {
    let (center, avg_radius) = fit_sphere_least_squares(points)?;

    if avg_radius < tolerance.linear() {
        return None;
    }

    // Verify ALL boundary points are equidistant from center.
    let max_deviation: f64 = points
        .iter()
        .map(|p| (p.distance_to(&center) - avg_radius).abs())
        .fold(0.0_f64, f64::max);

    // Loosened from 10× to 100× after §5.3: boundary sample
    // points come off NURBS fits with ~1e-5 drift on unit-scale
    // geometry, and the 10× window was rejecting legitimate cube-
    // corner spheres. The sphere is an exact analytic surface —
    // if its max_deviation fits within 100× tol_linear, it's the
    // right answer and beats the polar-NURBS fallback on both
    // representation clarity and allocation cost.
    // SCALE-RELATIVE gate: the boundary samples come off NURBS fits
    // whose drift scales with the corner radius, so an absolute
    // window rejects exact cube-corner spheres at radii above ~1
    // (measured: uniform r=3 corner — sphere declined, the swept
    // fallback shipped a boundary lying 0.0712·r off the strips, and
    // check 9 caught it at serving; at r ≤ 1 the same corner is
    // machine-exact). The sphere is still verified as the analytic
    // answer, just against a fit-noise window that grows with r.
    let boundary_diag = {
        let (mut lo, mut hi) = (
            Point3::new(f64::MAX, f64::MAX, f64::MAX),
            Point3::new(f64::MIN, f64::MIN, f64::MIN),
        );
        for p in points {
            lo = Point3::new(lo.x.min(p.x), lo.y.min(p.y), lo.z.min(p.z));
            hi = Point3::new(hi.x.max(p.x), hi.y.max(p.y), hi.z.max(p.z));
        }
        lo.distance_to(&hi)
    };
    // Plausibility guard (PK-233): co-circular samples are exactly
    // equidistant from EVERY center on the circle axis, so a
    // near-singular fit on near-coplanar data can push the center
    // arbitrarily far away and still pass the deviation gate. A
    // corner cap's sphere radius is O(boundary size) by
    // construction (rolling-ball center sits one fillet radius off
    // each support); a radius far beyond the boundary diagonal
    // means the samples did not determine the sphere.
    if avg_radius > boundary_diag * 100.0 {
        return None;
    }

    let sphere_gate = (tolerance.linear() * 100.0).max(boundary_diag * 1.0e-3);
    if popkern_core::trace::trace_enabled("sphere-fit") {
        eprintln!(
            "[spherefit] rolling-ball sphere candidate r={avg_radius:.4} max_dev={max_deviation:.4e} gate={sphere_gate:.4e}"
        );
    }
    if max_deviation > sphere_gate {
        return None;
    }

    // G1 tangency certification (PK-233). Boundary-on-sphere is
    // necessary but NOT sufficient: the uniform circular cube
    // corner's hexagonal arc+bridge ring fits a sphere at 1.4e-4
    // (inside the G0 gate above) whose radius is 1.48× the strip
    // radius and which is tangent to nothing (measured: fitted
    // r = 0.2225 vs strip r = 0.15; the committed rim edge then
    // validated 6.2e-2 off the cap). A serving rolling-ball sphere
    // must be tangent to every incident fillet strip: at samples of
    // each FilletEnd segment the strip's surface normal must be
    // radial to the fitted center.
    if !sphere_tangent_to_strips(
        boundary,
        fillet_surfaces,
        &center,
        avg_radius,
        popkern_core::trace::trace_enabled("sphere-fit"),
    ) {
        return None;
    }

    let _ = CornerMethod::RollingBallSphere; // classification tag
    let surface = SurfaceGeometry::Sphere {
        center,
        radius: avg_radius,
        axis: popkern_core::geometry::Vec3::Z,
    };

    Some(BlendSurfaceResult {
        surface,
        boundary_curves: boundary_curves.to_vec(),
        method: {
            if popkern_core::trace::trace_enabled("sphere-fit") {
                eprintln!("[corner-method] rolling-ball");
            }
            CornerMethod::RollingBallSphere
        },
        continuity: ContinuityLevel::G1,
    })
}

#[cfg(test)]
mod tests {
    use popkern_core::geometry::curve::CurveGeometry;
    use popkern_core::geometry::{Point3, Vec3};
    use popkern_core::topology::VertexHandle;

    use super::super::boundary::{BlendBoundary, BoundarySegment};
    use super::sphere_tangent_to_strips;
    use crate::fillet_surface::FilletSurfaceResult;

    const R: f64 = 0.15;

    /// The three fillet cylinders of a unit-cube origin corner with
    /// uniform radius R: axes parallel to x/y/z through the
    /// concurrency point (R, R, R).
    fn corner_strips() -> Vec<FilletSurfaceResult> {
        let axes = [
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
        ];
        axes.iter()
            .map(|axis| {
                let arc = end_circle_arc(*axis, R);
                FilletSurfaceResult {
                    surface: popkern_core::geometry::surface::SurfaceGeometry::Cylinder {
                        origin: Point3::new(R, R, R),
                        axis: *axis,
                        radius: R,
                    },
                    trim_curve_1: arc.clone(),
                    trim_curve_2: arc,
                    continuity_downgrade: None,
                }
            })
            .collect()
    }

    /// Circle on the cylinder with the given axis, in the plane at
    /// `offset` along the axis from the concurrency point (R, R, R).
    /// `offset = 0` is the tangency great circle of the rolling-ball
    /// sphere; `offset != 0` is a displaced ring that still lies on
    /// A sphere (centered at the concurrency point, radius
    /// sqrt(offset² + R²)) but is tangent to nothing.
    fn end_circle_arc(axis: Vec3, offset: f64) -> popkern_core::geometry::NurbsCurve {
        // Circle center = concurrency point shifted by `offset` along the axis.
        let center = Point3::new(
            R + axis.x * offset,
            R + axis.y * offset,
            R + axis.z * offset,
        );
        CurveGeometry::Arc {
            center,
            normal: axis,
            radius: R,
            start_angle: 0.4,
            end_angle: 1.9,
        }
        .to_nurbs()
    }

    fn boundary_from_arcs(arcs: Vec<popkern_core::geometry::NurbsCurve>) -> BlendBoundary {
        BlendBoundary {
            vertex: VertexHandle::default(),
            vertex_position: Point3::new(0.0, 0.0, 0.0),
            segments: arcs
                .into_iter()
                .enumerate()
                .map(|(k, curve)| BoundarySegment::FilletEnd {
                    fillet_index: k,
                    curve,
                })
                .collect(),
            bridge_fallbacks: Vec::new(),
            sphere_fit_substituted: false,
            per_fillet_caps: Vec::new(),
            hex_boundary_cps_deg9: Vec::new(),
            fillet_end_arcs: Vec::new(),
        }
    }

    /// The true rolling-ball sphere — center at the axis concurrency
    /// point, radius = strip radius — certifies tangent against all
    /// three strips (exact: strip normal is radial by construction).
    #[test]
    fn rolling_ball_sphere_certifies_tangent() {
        let strips = corner_strips();
        let arcs = strips
            .iter()
            .map(|s| match s.surface {
                popkern_core::geometry::surface::SurfaceGeometry::Cylinder { axis, .. } => {
                    end_circle_arc(axis, 0.0)
                }
                _ => unreachable!(),
            })
            .collect();
        let boundary = boundary_from_arcs(arcs);
        assert!(sphere_tangent_to_strips(
            &boundary,
            &strips,
            &Point3::new(R, R, R),
            R,
            false,
        ));
    }

    /// The impostor class that motivated the certification: rings at
    /// axial offset d on each cylinder lie EXACTLY on the sphere of
    /// radius sqrt(d² + R²) at the concurrency point (zero G0
    /// deviation — the boundary-fit gate alone accepts it), but that
    /// sphere is tangent to none of the strips and must be rejected.
    #[test]
    fn displaced_ring_sphere_rejected_despite_exact_boundary_fit() {
        let offset = 0.15;
        let strips = corner_strips();
        let arcs: Vec<_> = strips
            .iter()
            .map(|s| match s.surface {
                popkern_core::geometry::surface::SurfaceGeometry::Cylinder { axis, .. } => {
                    end_circle_arc(axis, offset)
                }
                _ => unreachable!(),
            })
            .collect();
        let impostor_center = Point3::new(R, R, R);
        let impostor_radius = (offset * offset + R * R).sqrt();
        // Document the trap: the boundary lies on the impostor
        // sphere to machine precision.
        use popkern_core::geometry::curve::CurveEval;
        for arc in &arcs {
            let (t0, t1) = arc.parameter_range();
            for k in 0..7 {
                let p = arc.point_at(t0 + (t1 - t0) * k as f64 / 6.0);
                assert!(
                    (p.distance_to(&impostor_center) - impostor_radius).abs() < 1e-9,
                    "impostor construction must be an exact boundary fit"
                );
            }
        }
        let boundary = boundary_from_arcs(arcs);
        assert!(!sphere_tangent_to_strips(
            &boundary,
            &strips,
            &impostor_center,
            impostor_radius,
            false,
        ));
    }
}
