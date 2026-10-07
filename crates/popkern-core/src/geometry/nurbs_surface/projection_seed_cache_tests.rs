//! Seed-field memoization tests (PK-513).

use crate::geometry::Point3;
use crate::geometry::PortableMath;
use crate::geometry::nurbs_surface::NurbsSurface;
use crate::geometry::surface::{closest_point_from_seed_points, closest_point_from_seeds};

/// A degree-1 helical ribbon: 2 u rows (inner/outer radius), `n` v
/// stations over `turns` turns — a coiled multi-sheet chart.
fn coil(n: usize, turns: f64) -> NurbsSurface {
    let mut knots_v = vec![0.0];
    for i in 0..n {
        knots_v.push(i as f64 / (n - 1) as f64);
    }
    knots_v.push(1.0);
    let mut control_points = Vec::new();
    for r in [1.0, 1.2] {
        for i in 0..n {
            let t = i as f64 / (n - 1) as f64;
            let a = std::f64::consts::TAU * turns * t;
            control_points.push(Point3::new(r * a.pcos(), r * a.psin(), 0.3 * turns * t));
        }
    }
    NurbsSurface {
        degree_u: 1,
        degree_v: 1,
        knots_u: vec![0.0, 0.0, 1.0, 1.0],
        knots_v,
        control_points,
        num_u: 2,
        num_v: n,
        weights: vec![1.0; 2 * n],
    }
}

#[test]
fn the_memoized_field_ranks_exactly_like_the_evaluating_projection() {
    let surface = coil(181, 6.0);
    let Some(field) = surface.projection_seed_field() else {
        panic!("a coiled chart has seeds");
    };
    for k in 0..40 {
        let t = (k as f64 + 0.37) / 40.0;
        let a = std::f64::consts::TAU * 6.0 * t;
        let query = Point3::new(1.1 * a.pcos(), 1.1 * a.psin(), 1.8 * t + 0.01);
        let evaluated = closest_point_from_seeds(&surface, &query, &field.seeds, (false, false), 6);
        let cached = closest_point_from_seed_points(
            &surface,
            &query,
            &field.seeds,
            &field.points,
            (false, false),
            6,
        );
        let bits = |r: Option<(f64, f64)>| r.map(|(u, v)| (u.to_bits(), v.to_bits()));
        assert_eq!(bits(evaluated), bits(cached), "query {k}");
    }
}

#[test]
fn a_repeat_lookup_hits_and_a_changed_surface_misses() {
    let surface = coil(121, 4.0);
    let (Some(first), Some(second)) = (
        surface.projection_seed_field(),
        surface.projection_seed_field(),
    ) else {
        panic!("seed fields");
    };
    assert!(std::rc::Rc::ptr_eq(&first, &second));
    let mut moved = surface.clone();
    moved.control_points[5].z += 1e-9;
    let Some(third) = moved.projection_seed_field() else {
        panic!("seed field");
    };
    assert!(!std::rc::Rc::ptr_eq(&first, &third));
}
