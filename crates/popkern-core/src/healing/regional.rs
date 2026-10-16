//! Regional surface analysis: curvature bounds, planarity detection,
//! and degeneracy detection over parameter-space regions.
//!
//! These functions sample a surface on a uniform grid within a given
//! parameter rectangle and compute aggregate statistics. They are useful
//! for healing workflows that need to classify surface patches before
//! deciding on a repair strategy.

use crate::geometry::Vec3;
use crate::geometry::surface::SurfaceEval;
use crate::healing::curvature::fundamental_forms::first_fundamental_form;
use crate::healing::curvature::surface_curvature::surface_curvature_at;
use crate::healing::interval_curvature::GuaranteedCurvatureBounds;

/// Curvature bounds over a parameter region.
#[derive(Debug, Clone, Copy)]
pub struct CurvatureBounds {
    /// Minimum Gaussian curvature sampled in the region.
    pub gaussian_min: f64,
    /// Maximum Gaussian curvature sampled in the region.
    pub gaussian_max: f64,
    /// Minimum mean curvature sampled in the region.
    pub mean_min: f64,
    /// Maximum mean curvature sampled in the region.
    pub mean_max: f64,
    /// Minimum radius of curvature in the region (smallest among samples).
    pub min_radius: f64,
}

/// First and second order derivative bounds over a parameter region.
#[derive(Debug, Clone, Copy)]
pub struct DerivativeBounds {
    /// Maximum magnitude of dS/du over the region.
    pub du_max: f64,
    /// Maximum magnitude of dS/dv over the region.
    pub dv_max: f64,
    /// Maximum magnitude of d²S/du² over the region.
    pub duu_max: f64,
    /// Maximum magnitude of d²S/dudv over the region.
    pub duv_max: f64,
    /// Maximum magnitude of d²S/dv² over the region.
    pub dvv_max: f64,
}

/// Regional analysis result combining curvature bounds, planarity, and degeneracy.
#[derive(Debug, Clone)]
pub struct RegionalAnalysis {
    /// Curvature bounds, or `None` if all sample points were singular.
    pub curvature_bounds: Option<CurvatureBounds>,
    /// Whether the region is effectively planar within the given tolerance.
    pub is_planar: bool,
    /// Maximum deviation from the best-fit plane (model units).
    pub planarity_deviation: f64,
    /// Whether the region contains degenerate areas (near-zero area element).
    pub has_degeneracy: bool,
    /// Minimum area element found across the sampled region.
    pub min_area_element: f64,
    /// First and second order derivative bounds, or `None` if not computed.
    pub derivative_bounds: Option<DerivativeBounds>,
    /// Rigorous curvature enclosures (certified interval arithmetic, see
    /// [`crate::healing::interval_curvature`]), or `None` if not
    /// computed. They are the whole line when a cell cannot be enclosed
    /// (`degenerate_cells > 0`).
    pub guaranteed_curvature_bounds: Option<GuaranteedCurvatureBounds>,
}

/// Generate uniformly spaced sample parameters over a 2D range.
///
/// Returns a `Vec` of `(u, v)` pairs on a `density x density` grid.
fn sample_grid(u_range: (f64, f64), v_range: (f64, f64), density: usize) -> Vec<(f64, f64)> {
    let n = density.max(2);
    let mut samples = Vec::with_capacity(n * n);
    for i in 0..n {
        let u = u_range.0 + (u_range.1 - u_range.0) * (i as f64 / (n - 1) as f64);
        for j in 0..n {
            let v = v_range.0 + (v_range.1 - v_range.0) * (j as f64 / (n - 1) as f64);
            samples.push((u, v));
        }
    }
    samples
}

/// Analyze curvature bounds over a parameter region.
///
/// Samples the surface on a `sample_density x sample_density` grid,
/// computes [`surface_curvature_at`] at each point, and tracks
/// min/max of Gaussian curvature, mean curvature, and minimum radius.
///
/// Returns `None` if all sample points are singular (degenerate).
pub fn curvature_bounds(
    surface: &dyn SurfaceEval,
    u_range: (f64, f64),
    v_range: (f64, f64),
    sample_density: usize,
) -> Option<CurvatureBounds> {
    let samples = sample_grid(u_range, v_range, sample_density);
    let mut g_min = f64::INFINITY;
    let mut g_max = f64::NEG_INFINITY;
    let mut m_min = f64::INFINITY;
    let mut m_max = f64::NEG_INFINITY;
    let mut r_min = f64::INFINITY;
    let mut count = 0usize;

    for (u, v) in &samples {
        if let Some(sc) = surface_curvature_at(surface, *u, *v) {
            g_min = g_min.min(sc.gaussian);
            g_max = g_max.max(sc.gaussian);
            m_min = m_min.min(sc.mean);
            m_max = m_max.max(sc.mean);
            r_min = r_min.min(sc.min_radius);
            count += 1;
        }
    }

    if count == 0 {
        return None;
    }

    Some(CurvatureBounds {
        gaussian_min: g_min,
        gaussian_max: g_max,
        mean_min: m_min,
        mean_max: m_max,
        min_radius: r_min,
    })
}

/// Detect if a surface region is effectively planar.
///
/// Samples points on a `sample_density x sample_density` grid, fits a
/// plane via centroid and average normal, and computes the maximum
/// deviation of any sample point from that plane.
///
/// Returns `(is_planar, max_deviation)` where `is_planar` is true when
/// `max_deviation < tolerance`.
pub fn detect_planarity(
    surface: &dyn SurfaceEval,
    u_range: (f64, f64),
    v_range: (f64, f64),
    tolerance: f64,
    sample_density: usize,
) -> (bool, f64) {
    let grid = sample_grid(u_range, v_range, sample_density);
    if grid.is_empty() {
        return (true, 0.0);
    }

    // Collect 3D sample points and normals.
    let mut points = Vec::with_capacity(grid.len());
    let mut normals = Vec::with_capacity(grid.len());
    for (u, v) in &grid {
        points.push(surface.point_at(*u, *v));
        normals.push(surface.normal_at(*u, *v));
    }

    // Compute centroid.
    let n = points.len() as f64;
    let cx = points.iter().map(|p| p.x).sum::<f64>() / n;
    let cy = points.iter().map(|p| p.y).sum::<f64>() / n;
    let cz = points.iter().map(|p| p.z).sum::<f64>() / n;

    // Average normal (sum and normalize).
    let nx = normals.iter().map(|v| v.x).sum::<f64>();
    let ny = normals.iter().map(|v| v.y).sum::<f64>();
    let nz = normals.iter().map(|v| v.z).sum::<f64>();
    let avg_normal = Vec3::new(nx, ny, nz);
    let norm_len = avg_normal.length();

    // If the average normal is degenerate, we cannot determine planarity.
    if norm_len < 1e-15 {
        return (false, f64::INFINITY);
    }
    let unit_normal = avg_normal * (1.0 / norm_len);

    // Compute max deviation: |dot(point - centroid, unit_normal)|.
    let mut max_dev = 0.0_f64;
    for p in &points {
        let dx = p.x - cx;
        let dy = p.y - cy;
        let dz = p.z - cz;
        let dev = (dx * unit_normal.x + dy * unit_normal.y + dz * unit_normal.z).abs();
        max_dev = max_dev.max(dev);
    }

    (max_dev < tolerance, max_dev)
}

/// Detect degenerate regions where the area element approaches zero.
///
/// Samples the first fundamental form on a `sample_density x sample_density`
/// grid and computes the area element at each point. A region is considered
/// degenerate if any sample's area element falls below `threshold`.
///
/// Returns `(has_degeneracy, min_area_element)`.
pub fn detect_degeneracy(
    surface: &dyn SurfaceEval,
    u_range: (f64, f64),
    v_range: (f64, f64),
    threshold: f64,
    sample_density: usize,
) -> (bool, f64) {
    let grid = sample_grid(u_range, v_range, sample_density);
    let mut min_area = f64::INFINITY;

    for (u, v) in &grid {
        let ff1 = first_fundamental_form(surface, *u, *v);
        let area = ff1.area_element();
        min_area = min_area.min(area);
    }

    // If no samples were taken, report no degeneracy with infinite area.
    if min_area == f64::INFINITY {
        return (false, f64::INFINITY);
    }

    (min_area < threshold, min_area)
}

/// Compute first and second order derivative bounds over a parameter region.
///
/// Samples the surface on a `sample_density x sample_density` grid and
/// evaluates partial derivatives at each point, tracking the maximum
/// magnitude of each derivative component.
pub fn derivative_bounds(
    surface: &dyn SurfaceEval,
    u_range: (f64, f64),
    v_range: (f64, f64),
    sample_density: usize,
) -> DerivativeBounds {
    let grid = sample_grid(u_range, v_range, sample_density);
    let mut du_max = 0.0_f64;
    let mut dv_max = 0.0_f64;
    let mut duu_max = 0.0_f64;
    let mut duv_max = 0.0_f64;
    let mut dvv_max = 0.0_f64;

    for (u, v) in &grid {
        du_max = du_max.max(surface.partial_u(*u, *v).length());
        dv_max = dv_max.max(surface.partial_v(*u, *v).length());
        duu_max = duu_max.max(surface.second_partial_uu(*u, *v).length());
        duv_max = duv_max.max(surface.second_partial_uv(*u, *v).length());
        dvv_max = dvv_max.max(surface.second_partial_vv(*u, *v).length());
    }

    DerivativeBounds {
        du_max,
        dv_max,
        duu_max,
        duv_max,
        dvv_max,
    }
}

/// Full regional analysis combining curvature bounds, planarity, and degeneracy.
///
/// Runs [`curvature_bounds`], [`detect_planarity`], and [`detect_degeneracy`]
/// over the same parameter region and returns a combined [`RegionalAnalysis`].
pub fn analyze_region(
    surface: &crate::geometry::surface::SurfaceGeometry,
    u_range: (f64, f64),
    v_range: (f64, f64),
    planarity_tolerance: f64,
    degeneracy_threshold: f64,
    sample_density: usize,
) -> RegionalAnalysis {
    let bounds = curvature_bounds(surface, u_range, v_range, sample_density);
    let (is_planar, planarity_deviation) = detect_planarity(
        surface,
        u_range,
        v_range,
        planarity_tolerance,
        sample_density,
    );
    let (has_degeneracy, min_area_element) = detect_degeneracy(
        surface,
        u_range,
        v_range,
        degeneracy_threshold,
        sample_density,
    );
    let deriv_bounds = derivative_bounds(surface, u_range, v_range, sample_density);
    let guaranteed = crate::healing::interval_curvature::guaranteed_curvature_bounds(
        surface, u_range, v_range, 0.1, 3,
    );

    RegionalAnalysis {
        curvature_bounds: bounds,
        is_planar,
        planarity_deviation,
        has_degeneracy,
        min_area_element,
        derivative_bounds: Some(deriv_bounds),
        guaranteed_curvature_bounds: Some(guaranteed),
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
    fn sphere_curvature_bounds_near_1_over_r_squared() {
        let r = 3.0;
        let sphere = SurfaceGeometry::Sphere {
            center: Point3::new(0.0, 0.0, 0.0),
            radius: r,
            axis: crate::geometry::Vec3::Z,
        };
        let expected_k = 1.0 / (r * r);

        // Sample away from poles: u in [0, 2pi], v in [pi/4, 3pi/4].
        let bounds = curvature_bounds(
            &sphere,
            (0.0, std::f64::consts::TAU),
            (
                std::f64::consts::FRAC_PI_4,
                3.0 * std::f64::consts::FRAC_PI_4,
            ),
            5,
        );
        let b = bounds.expect("sphere should have well-defined curvature");

        assert!(
            (b.gaussian_min - expected_k).abs() < 0.05,
            "gaussian_min={}, expected ~{}",
            b.gaussian_min,
            expected_k
        );
        assert!(
            (b.gaussian_max - expected_k).abs() < 0.05,
            "gaussian_max={}, expected ~{}",
            b.gaussian_max,
            expected_k
        );
        assert!(
            (b.min_radius - r).abs() < 0.5,
            "min_radius={}, expected ~{}",
            b.min_radius,
            r
        );
    }

    #[test]
    fn plane_curvature_bounds_near_zero() {
        let plane = xy_plane();
        let bounds = curvature_bounds(&plane, (0.0, 1.0), (0.0, 1.0), 5);
        let b = bounds.expect("plane should have well-defined curvature");
        assert!(
            b.gaussian_max.abs() < 1e-10,
            "plane gaussian_max should be ~0, got {}",
            b.gaussian_max
        );
        assert!(
            b.mean_max.abs() < 1e-10,
            "plane mean_max should be ~0, got {}",
            b.mean_max
        );
    }

    #[test]
    fn detect_planarity_on_plane() {
        let plane = xy_plane();
        let (is_planar, deviation) = detect_planarity(&plane, (0.0, 1.0), (0.0, 1.0), 1e-6, 5);
        assert!(is_planar, "plane should be detected as planar");
        assert!(
            deviation < 1e-10,
            "plane deviation should be ~0, got {}",
            deviation
        );
    }

    #[test]
    fn detect_planarity_on_sphere_is_not_planar() {
        let sphere = SurfaceGeometry::Sphere {
            center: Point3::new(0.0, 0.0, 0.0),
            radius: 1.0,
            axis: crate::geometry::Vec3::Z,
        };
        // Sphere over a wide parameter range is definitely not planar.
        let (is_planar, deviation) = detect_planarity(
            &sphere,
            (0.0, std::f64::consts::TAU),
            (
                std::f64::consts::FRAC_PI_4,
                3.0 * std::f64::consts::FRAC_PI_4,
            ),
            1e-3,
            5,
        );
        assert!(!is_planar, "sphere should not be planar");
        assert!(
            deviation > 0.01,
            "sphere deviation should be significant, got {}",
            deviation
        );
    }

    #[test]
    fn detect_degeneracy_on_plane_no_degeneracy() {
        let plane = xy_plane();
        let (has_deg, min_area) = detect_degeneracy(&plane, (0.0, 1.0), (0.0, 1.0), 1e-6, 5);
        assert!(!has_deg, "plane should have no degeneracy");
        assert!(
            (min_area - 1.0).abs() < 1e-10,
            "plane area element should be 1.0, got {}",
            min_area
        );
    }

    #[test]
    fn analyze_region_plane() {
        let plane = xy_plane();
        let result = analyze_region(&plane, (0.0, 1.0), (0.0, 1.0), 1e-6, 1e-6, 5);
        assert!(result.is_planar);
        assert!(!result.has_degeneracy);
        assert!(result.curvature_bounds.is_some());
    }

    #[test]
    fn analyze_region_sphere() {
        let sphere = SurfaceGeometry::Sphere {
            center: Point3::new(0.0, 0.0, 0.0),
            radius: 2.0,
            axis: crate::geometry::Vec3::Z,
        };
        let result = analyze_region(
            &sphere,
            (0.0, std::f64::consts::TAU),
            (
                std::f64::consts::FRAC_PI_4,
                3.0 * std::f64::consts::FRAC_PI_4,
            ),
            1e-3,
            1e-6,
            5,
        );
        assert!(!result.is_planar);
        assert!(!result.has_degeneracy);
        assert!(result.curvature_bounds.is_some());
    }

    #[test]
    fn sample_grid_produces_correct_count() {
        let grid = sample_grid((0.0, 1.0), (0.0, 1.0), 4);
        assert_eq!(grid.len(), 16);
    }

    #[test]
    fn sample_grid_minimum_density() {
        // Density of 1 should be clamped to 2.
        let grid = sample_grid((0.0, 1.0), (0.0, 1.0), 1);
        assert_eq!(grid.len(), 4);
    }
}
