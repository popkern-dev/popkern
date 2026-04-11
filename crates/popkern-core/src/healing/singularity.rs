//! Singularity detection on parametric surfaces.
//!
//! Detects degenerate points where the surface parameterization collapses
//! (e.g., sphere poles, cone apex). Uses both numerical sampling of the
//! cross product |S_u x S_v| and analytical knowledge of known surface types.

use crate::geometry::surface::SurfaceEval;
use serde::{Deserialize, Serialize};

/// Classification of a surface singularity.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum SingularityKind {
    /// A true pole where the surface collapses to a point (e.g., sphere pole).
    Pole,
    /// A parameter line that collapses to zero length.
    CollapsedParameterLine,
    /// A near-degenerate point where the cross product magnitude is very small
    /// but not quite zero.
    NearSingular,
}

/// A detected singularity at a specific parameter location.
#[derive(Debug, Clone)]
pub struct SingularityRecord {
    /// The type of singularity.
    pub kind: SingularityKind,
    /// The u parameter where the singularity was detected.
    pub u: f64,
    /// The v parameter where the singularity was detected.
    pub v: f64,
    /// The magnitude of |S_u x S_v| at this point.
    pub cross_product_magnitude: f64,
}

/// Detect singularities on a surface by sampling the cross product magnitude.
///
/// Samples the surface on a `sample_density x sample_density` grid within the
/// parameter range. At each sample point, computes |S_u x S_v| and classifies:
///
/// - `< epsilon_singular` -> [`SingularityKind::Pole`]
/// - `< 100 * epsilon_singular` -> [`SingularityKind::NearSingular`]
///
/// For analytical surface types (sphere, cone), also checks known pole
/// locations via [`detect_analytical_singularities`].
///
/// Returns all detected singularities, with duplicates possible when both
/// the grid and analytical checks find the same point.
pub fn detect_singularities(
    surface: &dyn SurfaceEval,
    epsilon_singular: f64,
    sample_density: usize,
) -> Vec<SingularityRecord> {
    let mut records = Vec::new();

    let ((u_min, u_max), (v_min, v_max)) = surface.parameter_range();

    // Clamp infinite ranges for sampling.
    let u_lo = u_min.max(-100.0);
    let u_hi = u_max.min(100.0);
    let v_lo = v_min.max(-100.0);
    let v_hi = v_max.min(100.0);

    let density = sample_density.max(2);

    for i in 0..density {
        let u = u_lo + (u_hi - u_lo) * (i as f64) / ((density - 1) as f64);
        for j in 0..density {
            let v = v_lo + (v_hi - v_lo) * (j as f64) / ((density - 1) as f64);

            let su = surface.partial_u(u, v);
            let sv = surface.partial_v(u, v);
            let cross = su.cross(&sv);
            let mag = cross.length();

            if mag < epsilon_singular {
                records.push(SingularityRecord {
                    kind: SingularityKind::Pole,
                    u,
                    v,
                    cross_product_magnitude: mag,
                });
            } else if mag < 100.0 * epsilon_singular {
                records.push(SingularityRecord {
                    kind: SingularityKind::NearSingular,
                    u,
                    v,
                    cross_product_magnitude: mag,
                });
            }
        }
    }

    // Supplement with analytical checks for known surface types.
    let analytical = detect_analytical_singularities(surface, epsilon_singular);
    records.extend(analytical);

    records
}

/// Detect singularities at analytically known locations.
///
/// Checks specific parameter locations for known surface types:
///
/// - **Sphere**: poles at v = 0 and v = PI (where longitude collapses)
/// - **Cone**: apex at v = 0 (where radius collapses to zero)
///
/// Evaluates the normal at each candidate point and checks if its magnitude
/// is below `epsilon_singular`.
pub fn detect_analytical_singularities(
    surface: &dyn SurfaceEval,
    epsilon_singular: f64,
) -> Vec<SingularityRecord> {
    let mut records = Vec::new();

    let ((u_min, u_max), (v_min, v_max)) = surface.parameter_range();

    // Sphere-like: check poles at v = v_min and v = v_max if range spans [0, PI].
    let pi = std::f64::consts::PI;
    let is_sphere_range = (v_min - 0.0).abs() < 1e-10 && (v_max - pi).abs() < 1e-10;

    if is_sphere_range {
        let u_mid = (u_min.max(-100.0) + u_max.min(100.0)) * 0.5;

        // North pole: v = 0
        check_analytical_point(surface, u_mid, 0.0, epsilon_singular, &mut records);
        // South pole: v = PI
        check_analytical_point(surface, u_mid, pi, epsilon_singular, &mut records);
    }

    // Cone-like: check apex at v = 0 if v_min is at 0.
    let tau = std::f64::consts::TAU;
    let is_cone_range = (v_min - 0.0).abs() < 1e-10
        && v_max > 1.0
        && (u_min - 0.0).abs() < 1e-10
        && (u_max - tau).abs() < 1e-10;

    if is_cone_range && !is_sphere_range {
        let u_mid = (u_min + u_max) * 0.5;
        check_analytical_point(surface, u_mid, 0.0, epsilon_singular, &mut records);
    }

    records
}

/// Check a single parameter point for singularity and push a record if found.
fn check_analytical_point(
    surface: &dyn SurfaceEval,
    u: f64,
    v: f64,
    epsilon_singular: f64,
    records: &mut Vec<SingularityRecord>,
) {
    let normal = surface.normal_at(u, v);
    let mag = normal.length();

    if mag < epsilon_singular {
        records.push(SingularityRecord {
            kind: SingularityKind::Pole,
            u,
            v,
            cross_product_magnitude: mag,
        });
    } else if mag < 100.0 * epsilon_singular {
        records.push(SingularityRecord {
            kind: SingularityKind::NearSingular,
            u,
            v,
            cross_product_magnitude: mag,
        });
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::surface::SurfaceGeometry;
    use crate::geometry::{Point3, Vec3};

    fn unit_sphere() -> SurfaceGeometry {
        SurfaceGeometry::Sphere {
            center: Point3::new(0.0, 0.0, 0.0),
            radius: 1.0,
        }
    }

    fn unit_cone() -> SurfaceGeometry {
        SurfaceGeometry::Cone {
            apex: Point3::new(0.0, 0.0, 0.0),
            axis: Vec3::new(0.0, 0.0, 1.0),
            half_angle: std::f64::consts::FRAC_PI_4,
        }
    }

    fn xy_plane() -> SurfaceGeometry {
        SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        }
    }

    #[test]
    fn sphere_poles_detected() {
        let sphere = unit_sphere();
        let eps = 1e-6;
        let results = detect_singularities(&sphere, eps, 20);

        // Should find singularities at v=0 (north pole) and v=PI (south pole).
        let poles: Vec<_> = results
            .iter()
            .filter(|r| r.kind == SingularityKind::Pole)
            .collect();
        assert!(
            poles.len() >= 2,
            "expected at least 2 poles on a sphere, got {}",
            poles.len()
        );

        // At least one pole near v=0 and one near v=PI.
        let has_north = poles.iter().any(|r| r.v.abs() < 0.1);
        let has_south = poles
            .iter()
            .any(|r| (r.v - std::f64::consts::PI).abs() < 0.1);
        assert!(has_north, "expected a pole near v=0 (north)");
        assert!(has_south, "expected a pole near v=PI (south)");
    }

    #[test]
    fn cone_apex_detected() {
        let cone = unit_cone();
        let eps = 1e-6;
        let results = detect_singularities(&cone, eps, 20);

        // Should find singularity at v=0 (apex).
        let at_apex: Vec<_> = results.iter().filter(|r| r.v.abs() < 0.1).collect();
        assert!(
            !at_apex.is_empty(),
            "expected a singularity at cone apex (v=0)"
        );
        assert!(
            at_apex.iter().any(|r| r.kind == SingularityKind::Pole),
            "expected a Pole classification at cone apex"
        );
    }

    #[test]
    fn plane_has_no_singularities() {
        let plane = xy_plane();
        let eps = 1e-6;
        let results = detect_singularities(&plane, eps, 20);

        assert!(
            results.is_empty(),
            "a plane should have no singularities, got {}",
            results.len()
        );
    }

    #[test]
    fn analytical_sphere_poles() {
        let sphere = unit_sphere();
        let eps = 1e-6;
        let results = detect_analytical_singularities(&sphere, eps);

        assert!(
            results.len() >= 2,
            "analytical check should find at least 2 sphere poles, got {}",
            results.len()
        );
    }

    #[test]
    fn analytical_cone_apex() {
        let cone = unit_cone();
        let eps = 1e-6;
        let results = detect_analytical_singularities(&cone, eps);

        assert!(
            !results.is_empty(),
            "analytical check should find cone apex singularity"
        );
        assert_eq!(results[0].kind, SingularityKind::Pole);
    }

    #[test]
    fn singularity_kind_equality() {
        assert_eq!(SingularityKind::Pole, SingularityKind::Pole);
        assert_ne!(SingularityKind::Pole, SingularityKind::NearSingular);
        assert_ne!(
            SingularityKind::CollapsedParameterLine,
            SingularityKind::NearSingular
        );
    }
}
