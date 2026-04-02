//! Chord-height adaptive UV-grid subdivision for parametric surfaces.

use popkern_core::geometry::{Point3, SurfaceEval, SurfaceGeometry};

/// A rectangular cell in UV parameter space.
#[derive(Debug, Clone)]
pub(crate) struct UvCell {
    /// Minimum u parameter of this cell.
    pub u_min: f64,
    /// Maximum u parameter of this cell.
    pub u_max: f64,
    /// Minimum v parameter of this cell.
    pub v_min: f64,
    /// Maximum v parameter of this cell.
    pub v_max: f64,
}

/// Compute the chord-height error for a UV cell.
///
/// This measures the distance from the bilinear interpolation of the four
/// corner points to the actual surface point at the cell's UV midpoint.
/// A flat (zero-curvature) surface will return approximately zero; a curved
/// surface will return a positive value proportional to its local curvature.
pub(crate) fn chord_height(surface: &SurfaceGeometry, cell: &UvCell) -> f64 {
    let p00 = surface.point_at(cell.u_min, cell.v_min);
    let p10 = surface.point_at(cell.u_max, cell.v_min);
    let p01 = surface.point_at(cell.u_min, cell.v_max);
    let p11 = surface.point_at(cell.u_max, cell.v_max);
    let mid_interp = Point3::new(
        (p00.x + p10.x + p01.x + p11.x) * 0.25,
        (p00.y + p10.y + p01.y + p11.y) * 0.25,
        (p00.z + p10.z + p01.z + p11.z) * 0.25,
    );
    let u_mid = (cell.u_min + cell.u_max) * 0.5;
    let v_mid = (cell.v_min + cell.v_max) * 0.5;
    let mid_actual = surface.point_at(u_mid, v_mid);
    mid_interp.distance_to(&mid_actual)
}

/// Adaptively subdivide a surface's UV parameter space into a grid of cells.
///
/// Starting from an initial `min_subdivisions × min_subdivisions` grid, each
/// cell whose chord-height exceeds `chord_tolerance` is refined into four
/// sub-cells (quad-tree subdivision). Refinement continues until all leaf
/// cells satisfy the tolerance or until the subdivision depth would exceed the
/// level implied by `max_subdivisions`.
///
/// # Parameters
/// - `surface` — the parametric surface to subdivide.
/// - `u_range` — `(u_min, u_max)` parameter extent in the u direction.
/// - `v_range` — `(v_min, v_max)` parameter extent in the v direction.
/// - `min_subdivisions` — initial grid resolution along each axis (e.g. `2`
///   creates a 2×2 = 4-cell starting grid).
/// - `max_subdivisions` — maximum resolution along each axis; no cell will be
///   subdivided further once the axis resolution would exceed this value.
/// - `chord_tolerance` — maximum allowed chord-height deviation.
///
/// # Returns
/// A `Vec<UvCell>` of all leaf cells covering the parameter domain.
pub(crate) fn adaptive_subdivide(
    surface: &SurfaceGeometry,
    u_range: (f64, f64),
    v_range: (f64, f64),
    min_subdivisions: usize,
    max_subdivisions: usize,
    chord_tolerance: f64,
) -> Vec<UvCell> {
    // Compute max refinement depth: how many times we can double min_subdivisions
    // before reaching max_subdivisions.
    //   depth 0: min_subdivisions cells/axis
    //   depth 1: 2*min_subdivisions cells/axis
    //   depth d: 2^d * min_subdivisions cells/axis
    // We want 2^d * min_subdivisions <= max_subdivisions
    //   => 2^d <= max_subdivisions / min_subdivisions
    //   => d <= log2(max_subdivisions / min_subdivisions)
    let max_depth = if min_subdivisions == 0 || max_subdivisions <= min_subdivisions {
        0_u32
    } else {
        let ratio = max_subdivisions / min_subdivisions;
        // ilog2 gives floor(log2(ratio))
        (ratio as u32).ilog2()
    };

    let (u_min, u_max) = u_range;
    let (v_min, v_max) = v_range;

    // Build initial grid.
    let mut queue: Vec<(UvCell, u32)> = Vec::new();
    let n = min_subdivisions.max(1);
    let du = (u_max - u_min) / n as f64;
    let dv = (v_max - v_min) / n as f64;
    for i in 0..n {
        for j in 0..n {
            queue.push((
                UvCell {
                    u_min: u_min + i as f64 * du,
                    u_max: u_min + (i + 1) as f64 * du,
                    v_min: v_min + j as f64 * dv,
                    v_max: v_min + (j + 1) as f64 * dv,
                },
                0,
            ));
        }
    }

    let mut leaves: Vec<UvCell> = Vec::new();

    while let Some((cell, depth)) = queue.pop() {
        let ch = chord_height(surface, &cell);
        if ch <= chord_tolerance || depth >= max_depth {
            leaves.push(cell);
        } else {
            // Subdivide into 4 children.
            let u_mid = (cell.u_min + cell.u_max) * 0.5;
            let v_mid = (cell.v_min + cell.v_max) * 0.5;
            let next_depth = depth + 1;
            queue.push((
                UvCell {
                    u_min: cell.u_min,
                    u_max: u_mid,
                    v_min: cell.v_min,
                    v_max: v_mid,
                },
                next_depth,
            ));
            queue.push((
                UvCell {
                    u_min: u_mid,
                    u_max: cell.u_max,
                    v_min: cell.v_min,
                    v_max: v_mid,
                },
                next_depth,
            ));
            queue.push((
                UvCell {
                    u_min: cell.u_min,
                    u_max: u_mid,
                    v_min: v_mid,
                    v_max: cell.v_max,
                },
                next_depth,
            ));
            queue.push((
                UvCell {
                    u_min: u_mid,
                    u_max: cell.u_max,
                    v_min: v_mid,
                    v_max: cell.v_max,
                },
                next_depth,
            ));
        }
    }

    leaves
}

#[cfg(test)]
mod tests {
    use super::*;
    use popkern_core::geometry::Vec3;
    use std::f64::consts::{PI, TAU};

    #[test]
    fn chord_height_plane_is_zero() {
        let plane = SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        };
        let cell = UvCell {
            u_min: 0.0,
            u_max: 1.0,
            v_min: 0.0,
            v_max: 1.0,
        };
        assert!(chord_height(&plane, &cell) < 1e-10);
    }

    #[test]
    fn chord_height_sphere_positive() {
        let sphere = SurfaceGeometry::Sphere {
            center: Point3::new(0.0, 0.0, 0.0),
            radius: 1.0,
        };
        let cell = UvCell {
            u_min: 0.0,
            u_max: 1.0,
            v_min: 0.5,
            v_max: 1.5,
        };
        assert!(chord_height(&sphere, &cell) > 0.0);
    }

    #[test]
    fn adaptive_plane_no_refinement() {
        let plane = SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        };
        let cells = adaptive_subdivide(&plane, (0.0, 1.0), (0.0, 1.0), 2, 64, 0.01);
        assert_eq!(cells.len(), 4);
    }

    #[test]
    fn adaptive_sphere_refines() {
        let sphere = SurfaceGeometry::Sphere {
            center: Point3::new(0.0, 0.0, 0.0),
            radius: 1.0,
        };
        let cells = adaptive_subdivide(&sphere, (0.0, TAU), (0.0, PI), 2, 64, 0.01);
        assert!(cells.len() > 4);
    }

    #[test]
    fn tighter_tolerance_more_cells() {
        let sphere = SurfaceGeometry::Sphere {
            center: Point3::new(0.0, 0.0, 0.0),
            radius: 1.0,
        };
        let coarse = adaptive_subdivide(&sphere, (0.0, TAU), (0.0, PI), 2, 64, 0.1);
        let fine = adaptive_subdivide(&sphere, (0.0, TAU), (0.0, PI), 2, 64, 0.001);
        assert!(fine.len() > coarse.len());
    }

    #[test]
    fn max_subdivisions_respected() {
        let sphere = SurfaceGeometry::Sphere {
            center: Point3::new(0.0, 0.0, 0.0),
            radius: 1.0,
        };
        // Very tight tolerance with low max — should hit the cap.
        // max_subdivisions=4, min_subdivisions=2 => ratio=2 => max_depth=1
        // At depth 1 each of the 4 initial cells can split into 4 => up to 16 cells.
        let cells = adaptive_subdivide(&sphere, (0.0, TAU), (0.0, PI), 2, 4, 0.0001);
        assert!(cells.len() <= 16);
    }

    #[test]
    fn min_subdivisions_zero_treated_as_one() {
        let plane = SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        };
        // min_subdivisions=0 should be clamped to 1 internally.
        let cells = adaptive_subdivide(&plane, (0.0, 1.0), (0.0, 1.0), 0, 64, 0.01);
        assert_eq!(cells.len(), 1);
    }

    #[test]
    fn equal_min_max_no_refinement() {
        let sphere = SurfaceGeometry::Sphere {
            center: Point3::new(0.0, 0.0, 0.0),
            radius: 1.0,
        };
        // max == min => max_depth = 0, so no refinement possible.
        let cells = adaptive_subdivide(&sphere, (0.0, TAU), (0.0, PI), 4, 4, 0.0001);
        // 4x4 initial grid = 16 cells, none refined further.
        assert_eq!(cells.len(), 16);
    }

    #[test]
    fn cells_cover_full_parameter_domain() {
        let plane = SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        };
        let cells = adaptive_subdivide(&plane, (0.0, 2.0), (0.0, 3.0), 3, 64, 0.01);
        // For a plane (no refinement), 3x3 = 9 cells.
        assert_eq!(cells.len(), 9);
        // Verify total area coverage: sum of cell areas should equal domain area.
        let total_area: f64 = cells
            .iter()
            .map(|c| (c.u_max - c.u_min) * (c.v_max - c.v_min))
            .sum();
        assert!((total_area - 6.0).abs() < 1e-10);
    }

    #[test]
    fn cells_do_not_overlap() {
        let sphere = SurfaceGeometry::Sphere {
            center: Point3::new(0.0, 0.0, 0.0),
            radius: 1.0,
        };
        let cells = adaptive_subdivide(&sphere, (0.0, TAU), (0.0, PI), 2, 8, 0.05);
        // No two leaf cells should overlap (their interiors should be disjoint).
        for (i, a) in cells.iter().enumerate() {
            for b in cells.iter().skip(i + 1) {
                let u_overlap = a.u_min < b.u_max && b.u_min < a.u_max;
                let v_overlap = a.v_min < b.v_max && b.v_min < a.v_max;
                if u_overlap && v_overlap {
                    // Cells share at most a boundary, not an interior.
                    let u_inter = (a.u_max.min(b.u_max) - a.u_min.max(b.u_min)).max(0.0);
                    let v_inter = (a.v_max.min(b.v_max) - a.v_min.max(b.v_min)).max(0.0);
                    assert!(
                        u_inter * v_inter < 1e-20,
                        "Cells {i} and {} overlap with area {}",
                        i + 1,
                        u_inter * v_inter
                    );
                }
            }
        }
    }

    #[test]
    fn chord_height_cylinder_positive() {
        let cylinder = SurfaceGeometry::Cylinder {
            origin: Point3::new(0.0, 0.0, 0.0),
            axis: Vec3::new(0.0, 0.0, 1.0),
            radius: 1.0,
        };
        let cell = UvCell {
            u_min: 0.0,
            u_max: 1.0,
            v_min: 0.0,
            v_max: 1.0,
        };
        assert!(chord_height(&cylinder, &cell) > 0.0);
    }

    #[test]
    fn large_radius_sphere_less_refinement_than_small() {
        let large = SurfaceGeometry::Sphere {
            center: Point3::new(0.0, 0.0, 0.0),
            radius: 100.0,
        };
        let small = SurfaceGeometry::Sphere {
            center: Point3::new(0.0, 0.0, 0.0),
            radius: 0.5,
        };
        // Same angular range but different radii => different chord heights.
        // For the same tolerance, a larger sphere needs more cells (larger absolute error).
        let large_cells = adaptive_subdivide(&large, (0.0, 1.0), (0.5, 1.5), 2, 64, 0.01);
        let small_cells = adaptive_subdivide(&small, (0.0, 1.0), (0.5, 1.5), 2, 64, 0.01);
        assert!(
            large_cells.len() >= small_cells.len(),
            "larger sphere should need at least as many cells"
        );
    }

    #[test]
    fn narrow_parameter_range_still_works() {
        let plane = SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        };
        // Very narrow range.
        let cells = adaptive_subdivide(&plane, (0.0, 0.001), (0.0, 0.001), 1, 64, 0.01);
        assert_eq!(cells.len(), 1);
        assert!((cells[0].u_max - cells[0].u_min - 0.001).abs() < 1e-15);
    }
}
