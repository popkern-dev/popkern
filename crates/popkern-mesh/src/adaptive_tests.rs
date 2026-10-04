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
