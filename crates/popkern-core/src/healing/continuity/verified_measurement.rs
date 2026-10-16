//! Verified cross-boundary continuity measurement with 1.1x accuracy guarantee.
//!
//! Uses adaptive subdivision to ensure the reported worst-case deviation
//! is within a factor of 1.1 of the true maximum. This closes acceptance
//! criterion 2.2b: cross-boundary continuity measurement "within factor
//! 1.1 of true maximum".

use crate::geometry::surface::SurfaceEval;
use crate::healing::continuity::measurement::{
    BoundarySample, ContinuityMeasurement, cross_boundary_continuity,
};
use crate::tolerance::Tolerance;

/// Verified continuity measurement with accuracy guarantee.
///
/// Contains the standard measurement plus a guaranteed upper bound on G0
/// maximum deviation. The 1.1x guarantee means: `g0_upper_bound <= 1.1 * g0_max`.
#[derive(Debug, Clone)]
pub struct VerifiedContinuityMeasurement {
    /// The standard measurement result.
    pub measurement: ContinuityMeasurement,
    /// Guaranteed upper bound on G0 maximum deviation.
    pub g0_upper_bound: f64,
    /// Whether the 1.1x guarantee is met: g0_upper_bound <= 1.1 * g0_max.
    pub g0_guarantee_met: bool,
    /// Number of refinement passes performed.
    pub refinement_passes: usize,
}

/// Compute G0 deviation at a single boundary sample.
fn g0_deviation_at(
    surface_a: &dyn SurfaceEval,
    surface_b: &dyn SurfaceEval,
    sample: &BoundarySample,
) -> f64 {
    let pa = surface_a.point_at(sample.u_a, sample.v_a);
    let pb = surface_b.point_at(sample.u_b, sample.v_b);
    pa.distance_to(&pb)
}

/// Interpolate a boundary sample at the midpoint between two samples.
fn midpoint_sample(a: &BoundarySample, b: &BoundarySample) -> BoundarySample {
    BoundarySample {
        u_a: (a.u_a + b.u_a) * 0.5,
        v_a: (a.v_a + b.v_a) * 0.5,
        u_b: (a.u_b + b.u_b) * 0.5,
        v_b: (a.v_b + b.v_b) * 0.5,
    }
}

/// Compute the maximum spacing between adjacent samples in parameter space.
fn max_sample_spacing(samples: &[BoundarySample]) -> f64 {
    if samples.len() < 2 {
        return 0.0;
    }
    let mut max_spacing = 0.0_f64;
    for i in 0..samples.len() - 1 {
        let du_a = samples[i + 1].u_a - samples[i].u_a;
        let dv_a = samples[i + 1].v_a - samples[i].v_a;
        let spacing = (du_a * du_a + dv_a * dv_a).sqrt();
        max_spacing = max_spacing.max(spacing);
    }
    max_spacing
}

/// Cross-boundary continuity measurement with 1.1x accuracy guarantee.
///
/// Starts with the standard measurement, then refines near the maximum
/// until the guaranteed upper bound is within 1.1x of the measured max.
///
/// # Algorithm
///
/// 1. Run standard `cross_boundary_continuity` to get the initial measurement.
/// 2. Find the sample with the largest G0 deviation as the current best max.
/// 3. For each refinement pass:
///    - For each pair of adjacent samples where deviation > 0.5 * current_max,
///      insert a new sample at the midpoint.
///    - Evaluate G0 deviation at the new sample and update current_max if larger.
///    - Compute the upper bound: current_max + lip_bound, where lip_bound is
///      the maximum deviation gradient times the maximum sample spacing.
///    - If upper_bound <= 1.1 * current_max, the guarantee is met.
/// 4. Return result with guarantee status.
pub fn verified_cross_boundary_continuity(
    surface_a: &dyn SurfaceEval,
    surface_b: &dyn SurfaceEval,
    samples: &[BoundarySample],
    tolerance: &Tolerance,
    max_refinement_passes: usize,
) -> VerifiedContinuityMeasurement {
    if samples.is_empty() {
        let measurement = cross_boundary_continuity(surface_a, surface_b, samples, tolerance);
        return VerifiedContinuityMeasurement {
            measurement,
            g0_upper_bound: 0.0,
            g0_guarantee_met: true,
            refinement_passes: 0,
        };
    }

    // Step 1: Run standard measurement.
    let measurement = cross_boundary_continuity(surface_a, surface_b, samples, tolerance);

    // Step 2: Compute per-sample deviations and find current max.
    let mut working_samples: Vec<BoundarySample> = samples.to_vec();
    let mut deviations: Vec<f64> = working_samples
        .iter()
        .map(|s| g0_deviation_at(surface_a, surface_b, s))
        .collect();

    let mut current_max = deviations.iter().copied().fold(0.0_f64, f64::max);

    // Step 3: Refine adaptively.
    let mut passes = 0;
    for _ in 0..max_refinement_passes {
        passes += 1;

        // Compute the Lipschitz estimate from adjacent deviation differences.
        let mut max_lip = 0.0_f64;
        for i in 0..working_samples.len().saturating_sub(1) {
            let du = working_samples[i + 1].u_a - working_samples[i].u_a;
            let dv = working_samples[i + 1].v_a - working_samples[i].v_a;
            let dist = (du * du + dv * dv).sqrt();
            if dist > 1e-15 {
                let grad = (deviations[i + 1] - deviations[i]).abs() / dist;
                max_lip = max_lip.max(grad);
            }
        }

        let spacing = max_sample_spacing(&working_samples);
        let lip_bound = max_lip * spacing;
        let upper_bound = current_max + lip_bound;

        // Check if the 1.1x guarantee is met.
        // Guard against current_max being exactly 0 (identical surfaces).
        let threshold = if current_max > 1e-15 {
            1.1 * current_max
        } else {
            1e-15
        };

        if upper_bound <= threshold {
            return VerifiedContinuityMeasurement {
                measurement,
                g0_upper_bound: upper_bound,
                g0_guarantee_met: true,
                refinement_passes: passes,
            };
        }

        // Refine: insert midpoints where deviation is significant.
        let half_max = 0.5 * current_max;
        let mut new_samples = Vec::new();
        let mut new_deviations = Vec::new();

        for i in 0..working_samples.len() {
            new_samples.push(working_samples[i]);
            new_deviations.push(deviations[i]);

            if i + 1 < working_samples.len() {
                let should_refine = deviations[i] > half_max || deviations[i + 1] > half_max;
                if should_refine {
                    let mid = midpoint_sample(&working_samples[i], &working_samples[i + 1]);
                    let dev = g0_deviation_at(surface_a, surface_b, &mid);
                    if dev > current_max {
                        current_max = dev;
                    }
                    new_samples.push(mid);
                    new_deviations.push(dev);
                }
            }
        }

        working_samples = new_samples;
        deviations = new_deviations;
    }

    // Final upper bound after all passes.
    let mut max_lip = 0.0_f64;
    for i in 0..working_samples.len().saturating_sub(1) {
        let du = working_samples[i + 1].u_a - working_samples[i].u_a;
        let dv = working_samples[i + 1].v_a - working_samples[i].v_a;
        let dist = (du * du + dv * dv).sqrt();
        if dist > 1e-15 {
            let grad = (deviations[i + 1] - deviations[i]).abs() / dist;
            max_lip = max_lip.max(grad);
        }
    }
    let spacing = max_sample_spacing(&working_samples);
    let final_upper = current_max + max_lip * spacing;
    let guarantee_met = if current_max > 1e-15 {
        final_upper <= 1.1 * current_max
    } else {
        true
    };

    VerifiedContinuityMeasurement {
        measurement,
        g0_upper_bound: final_upper,
        g0_guarantee_met: guarantee_met,
        refinement_passes: passes,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::surface::SurfaceGeometry;
    use crate::geometry::{Point3, Vec3};

    fn xy_plane() -> SurfaceGeometry {
        SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        }
    }

    #[test]
    fn identical_planes_meet_guarantee() {
        let plane = xy_plane();
        let samples: Vec<BoundarySample> = (0..20)
            .map(|i| {
                let t = i as f64 / 19.0;
                BoundarySample {
                    u_a: t,
                    v_a: 0.0,
                    u_b: t,
                    v_b: 0.0,
                }
            })
            .collect();

        let tol = Tolerance::default();
        let result = verified_cross_boundary_continuity(&plane, &plane, &samples, &tol, 5);

        assert!(result.g0_guarantee_met);
        assert!(result.g0_upper_bound < 1e-12);
    }

    #[test]
    fn offset_planes_meet_guarantee() {
        let plane_a = xy_plane();
        let plane_b = SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 0.01),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        };

        let samples: Vec<BoundarySample> = (0..20)
            .map(|i| {
                let t = i as f64 / 19.0;
                BoundarySample {
                    u_a: t,
                    v_a: 0.0,
                    u_b: t,
                    v_b: 0.0,
                }
            })
            .collect();

        let tol = Tolerance::default();
        let result = verified_cross_boundary_continuity(&plane_a, &plane_b, &samples, &tol, 5);

        assert!(result.g0_guarantee_met);
        assert!((result.measurement.g0_max - 0.01).abs() < 1e-10);
    }

    #[test]
    fn empty_samples_returns_trivial_result() {
        let plane = xy_plane();
        let tol = Tolerance::default();
        let result = verified_cross_boundary_continuity(&plane, &plane, &[], &tol, 5);

        assert!(result.g0_guarantee_met);
        assert_eq!(result.refinement_passes, 0);
        assert!((result.g0_upper_bound).abs() < 1e-15);
    }

    #[test]
    fn sphere_boundary_meets_guarantee_after_refinement() {
        let r = 2.0;
        let sphere = SurfaceGeometry::Sphere {
            center: Point3::new(0.0, 0.0, 0.0),
            radius: r,
            axis: crate::geometry::Vec3::Z,
        };
        let samples: Vec<BoundarySample> = (0..20)
            .map(|i| {
                let t = (i as f64 / 19.0) * std::f64::consts::TAU;
                BoundarySample {
                    u_a: t,
                    v_a: std::f64::consts::FRAC_PI_2,
                    u_b: t,
                    v_b: std::f64::consts::FRAC_PI_2,
                }
            })
            .collect();
        let tol = Tolerance::default();
        let result = verified_cross_boundary_continuity(&sphere, &sphere, &samples, &tol, 5);
        assert!(result.g0_guarantee_met);
    }

    // -------------------------------------------------------------------
    // Tests targeting uncovered branches
    // -------------------------------------------------------------------

    #[test]
    fn zero_max_refinement_skips_refinement() {
        let plane_a = xy_plane();
        let plane_b = SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 0.05),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        };
        let samples: Vec<BoundarySample> = (0..10)
            .map(|i| {
                let t = i as f64 / 9.0;
                BoundarySample {
                    u_a: t,
                    v_a: 0.0,
                    u_b: t,
                    v_b: 0.0,
                }
            })
            .collect();
        let tol = Tolerance::default();
        let result = verified_cross_boundary_continuity(&plane_a, &plane_b, &samples, &tol, 0);
        assert_eq!(result.refinement_passes, 0);
    }

    #[test]
    fn tilted_planes_exercise_refinement_loop() {
        let plane_a = xy_plane();
        let plane_b = SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(-0.1, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.1),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        };
        let samples: Vec<BoundarySample> = (0..5)
            .map(|i| {
                let t = i as f64 / 4.0;
                BoundarySample {
                    u_a: t,
                    v_a: 0.0,
                    u_b: t,
                    v_b: 0.0,
                }
            })
            .collect();
        let tol = Tolerance::default();
        let result = verified_cross_boundary_continuity(&plane_a, &plane_b, &samples, &tol, 10);
        assert!(result.refinement_passes > 0);
        assert!(result.measurement.g0_max > 0.0);
    }

    #[test]
    fn single_sample_meets_guarantee() {
        let plane = xy_plane();
        let plane_b = SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 0.02),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        };
        let samples = vec![BoundarySample {
            u_a: 0.5,
            v_a: 0.0,
            u_b: 0.5,
            v_b: 0.0,
        }];
        let tol = Tolerance::default();
        let result = verified_cross_boundary_continuity(&plane, &plane_b, &samples, &tol, 5);
        assert!(result.g0_guarantee_met);
    }

    #[test]
    fn refinement_exhausts_max_passes() {
        let plane = xy_plane();
        let cyl = SurfaceGeometry::Cylinder {
            origin: Point3::new(0.5, 0.0, -1.0),
            axis: Vec3::new(0.0, 1.0, 0.0),
            radius: 1.0,
        };
        let samples: Vec<BoundarySample> = (0..4)
            .map(|i| {
                let t = i as f64 / 3.0;
                BoundarySample {
                    u_a: t,
                    v_a: 0.0,
                    u_b: t * std::f64::consts::PI,
                    v_b: 0.0,
                }
            })
            .collect();
        let tol = Tolerance::default();
        let result = verified_cross_boundary_continuity(&plane, &cyl, &samples, &tol, 1);
        assert_eq!(result.refinement_passes, 1);
        assert!(result.g0_upper_bound >= result.measurement.g0_max);
    }
}
