// This file is only compiled under #[cfg(test)] (see mod.rs: `#[cfg(test)] mod tests;`).
#[cfg(test)]
use super::*;
use crate::geometry::{Point3, SurfaceEval};

/// Helper to build a bilinear patch for testing.
fn bilinear_patch() -> NurbsSurface {
    NurbsSurface {
        degree_u: 1,
        degree_v: 1,
        knots_u: vec![0.0, 0.0, 1.0, 1.0],
        knots_v: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![
            // Row u=0: (0,0,0), (0,1,0)
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
            // Row u=1: (1,0,0), (1,1,0)
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
        ],
        num_u: 2,
        num_v: 2,
        weights: vec![1.0, 1.0, 1.0, 1.0],
    }
}

#[test]
fn bilinear_patch_center() {
    let surface = bilinear_patch();
    let p = surface.point_at(0.5, 0.5);
    assert!((p.x - 0.5).abs() < 1e-10);
    assert!((p.y - 0.5).abs() < 1e-10);
    assert!((p.z).abs() < 1e-10);
}

#[test]
fn bilinear_patch_corners() {
    let surface = bilinear_patch();
    let p00 = surface.point_at(0.0, 0.0);
    assert!((p00.x).abs() < 1e-10);
    assert!((p00.y).abs() < 1e-10);
    // Use slightly less than 1.0 to stay in valid knot span
    let p11 = surface.point_at(1.0 - 1e-12, 1.0 - 1e-12);
    assert!((p11.x - 1.0).abs() < 1e-6);
    assert!((p11.y - 1.0).abs() < 1e-6);
}

#[test]
fn bilinear_patch_normal_is_z() {
    let surface = bilinear_patch();
    let n = surface.normal_at(0.5, 0.5);
    // Normal should point in +z or -z direction for a flat XY patch
    let nlen = n.length();
    assert!(nlen > 1e-10, "Normal should be non-zero");
    let nz = n.z.abs() / nlen;
    assert!((nz - 1.0).abs() < 1e-6, "Normal should be ±z, got {:?}", n);
}

#[test]
fn parameter_range() {
    let surface = bilinear_patch();
    let ((u0, u1), (v0, v1)) = surface.parameter_range();
    assert_eq!(u0, 0.0);
    assert_eq!(u1, 1.0);
    assert_eq!(v0, 0.0);
    assert_eq!(v1, 1.0);
}

#[test]
fn bilinear_patch_partial_u() {
    let surface = bilinear_patch();
    let du = surface.partial_u(0.5, 0.5);
    // For a bilinear patch, d/du should be (1, 0, 0)
    assert!((du.x - 1.0).abs() < 1e-4);
    assert!((du.y).abs() < 1e-4);
}

#[test]
fn bilinear_patch_partial_v() {
    let surface = bilinear_patch();
    let dv = surface.partial_v(0.5, 0.5);
    // For a bilinear patch, d/dv should be (0, 1, 0)
    assert!((dv.x).abs() < 1e-4);
    assert!((dv.y - 1.0).abs() < 1e-4);
}

#[test]
fn bilinear_patch_at_end() {
    let surface = bilinear_patch();
    // Evaluate at exactly t=1.0 (clamped end)
    let p = surface.point_at(1.0, 1.0);
    assert!((p.x - 1.0).abs() < 1e-6);
    assert!((p.y - 1.0).abs() < 1e-6);
}

#[test]
fn bilinear_patch_at_before_start() {
    let surface = bilinear_patch();
    // Evaluate before start — clamped to degree span, extrapolates
    // Just verify it doesn't panic and returns finite values
    let p = surface.point_at(-0.1, -0.1);
    assert!(p.x.is_finite());
    assert!(p.y.is_finite());
}

#[test]
fn surface_insert_knot_u_preserves_shape() {
    let surface = bilinear_patch();
    let new_surface = surface.insert_knot_u(0.5).unwrap();
    // num_u grows by 1
    assert_eq!(new_surface.num_u, surface.num_u + 1);
    assert_eq!(new_surface.num_v, surface.num_v);
    // Shape preserved at multiple sample points
    for &u in &[0.0, 0.25, 0.5, 0.75, 0.99] {
        for &v in &[0.0, 0.25, 0.5, 0.75, 0.99] {
            let p1 = surface.point_at(u, v);
            let p2 = new_surface.point_at(u, v);
            assert!(
                (p1.x - p2.x).abs() < 1e-10,
                "x mismatch at u={u},v={v}: {p1:?} vs {p2:?}"
            );
            assert!(
                (p1.y - p2.y).abs() < 1e-10,
                "y mismatch at u={u},v={v}: {p1:?} vs {p2:?}"
            );
        }
    }
}

#[test]
fn surface_insert_knot_v_preserves_shape() {
    let surface = bilinear_patch();
    let new_surface = surface.insert_knot_v(0.5).unwrap();
    assert_eq!(new_surface.num_u, surface.num_u);
    assert_eq!(new_surface.num_v, surface.num_v + 1);
    for &u in &[0.0, 0.25, 0.5, 0.75, 0.99] {
        for &v in &[0.0, 0.25, 0.5, 0.75, 0.99] {
            let p1 = surface.point_at(u, v);
            let p2 = new_surface.point_at(u, v);
            assert!(
                (p1.x - p2.x).abs() < 1e-10,
                "x mismatch at u={u},v={v}: {p1:?} vs {p2:?}"
            );
            assert!(
                (p1.y - p2.y).abs() < 1e-10,
                "y mismatch at u={u},v={v}: {p1:?} vs {p2:?}"
            );
        }
    }
}

#[test]
fn surface_insert_knot_u_out_of_range() {
    let surface = bilinear_patch();
    assert!(surface.insert_knot_u(2.0).is_err());
    assert!(surface.insert_knot_u(-0.5).is_err());
}

#[test]
fn surface_insert_knot_v_out_of_range() {
    let surface = bilinear_patch();
    assert!(surface.insert_knot_v(2.0).is_err());
    assert!(surface.insert_knot_v(-0.5).is_err());
}

#[test]
fn surface_insert_knots_u_preserves_shape() {
    let surface = bilinear_patch();
    let new_surface = surface.insert_knots_u(&[0.25, 0.75]).unwrap();
    assert_eq!(new_surface.num_u, surface.num_u + 2);
    let p1 = surface.point_at(0.5, 0.5);
    let p2 = new_surface.point_at(0.5, 0.5);
    assert!((p1.x - p2.x).abs() < 1e-10);
    assert!((p1.y - p2.y).abs() < 1e-10);
}

#[test]
fn surface_insert_knots_v_preserves_shape() {
    let surface = bilinear_patch();
    let new_surface = surface.insert_knots_v(&[0.25, 0.75]).unwrap();
    assert_eq!(new_surface.num_v, surface.num_v + 2);
    let p1 = surface.point_at(0.5, 0.5);
    let p2 = new_surface.point_at(0.5, 0.5);
    assert!((p1.x - p2.x).abs() < 1e-10);
    assert!((p1.y - p2.y).abs() < 1e-10);
}

#[test]
fn surface_remove_knot_u_roundtrip() {
    let surface = bilinear_patch();
    // Insert a knot in u, then remove it
    let inserted = surface.insert_knot_u(0.5).unwrap();
    assert_eq!(inserted.num_u, surface.num_u + 1);
    assert_eq!(inserted.knots_u.len(), surface.knots_u.len() + 1);

    // The inserted knot is at index 2 in knots_u = [0,0, 0.5, 1,1]
    let removed = inserted.remove_knot_u(2, 0.0);
    assert!(removed.is_some(), "Should be able to remove inserted knot");
    let removed = removed.unwrap();
    assert_eq!(removed.num_u, surface.num_u);
    assert_eq!(removed.knots_u.len(), surface.knots_u.len());

    // Shape should be preserved
    for &u in &[0.0, 0.25, 0.5, 0.75, 0.99] {
        for &v in &[0.0, 0.25, 0.5, 0.75, 0.99] {
            let p1 = surface.point_at(u, v);
            let p2 = removed.point_at(u, v);
            assert!(
                (p1.x - p2.x).abs() < 1e-9,
                "x mismatch at u={u},v={v}: {p1:?} vs {p2:?}"
            );
            assert!(
                (p1.y - p2.y).abs() < 1e-9,
                "y mismatch at u={u},v={v}: {p1:?} vs {p2:?}"
            );
        }
    }
}

#[test]
fn surface_remove_knot_v_roundtrip() {
    let surface = bilinear_patch();
    // Insert a knot in v, then remove it
    let inserted = surface.insert_knot_v(0.5).unwrap();
    assert_eq!(inserted.num_v, surface.num_v + 1);

    // The inserted knot is at index 2 in knots_v = [0,0, 0.5, 1,1]
    let removed = inserted.remove_knot_v(2, 0.0);
    assert!(removed.is_some(), "Should be able to remove inserted knot");
    let removed = removed.unwrap();
    assert_eq!(removed.num_v, surface.num_v);
    assert_eq!(removed.knots_v.len(), surface.knots_v.len());
}

#[test]
fn surface_remove_knot_u_none_on_boundary() {
    let surface = bilinear_patch();
    // Cannot remove clamped boundary knots
    assert!(surface.remove_knot_u(0, 0.0).is_none());
    assert!(surface.remove_knot_u(3, 0.0).is_none());
}

#[test]
fn surface_remove_knot_v_none_on_boundary() {
    let surface = bilinear_patch();
    assert!(surface.remove_knot_v(0, 0.0).is_none());
    assert!(surface.remove_knot_v(3, 0.0).is_none());
}

#[test]
fn surface_remove_knot_u_times() {
    let surface = bilinear_patch();
    // Insert knot 0.5 twice in u
    let inserted = surface.insert_knots_u(&[0.5, 0.5]).unwrap();
    // knots_u = [0,0, 0.5, 0.5, 1,1], index 2 is the first inserted 0.5
    let (result, count) = inserted.remove_knot_u_times(2, 3, 0.0);
    assert_eq!(count, 2, "Should have removed 2 knots");
    assert_eq!(result.knots_u.len(), inserted.knots_u.len() - 2);
}

#[test]
fn surface_remove_knot_v_times() {
    let surface = bilinear_patch();
    // Insert knot 0.5 twice in v
    let inserted = surface.insert_knots_v(&[0.5, 0.5]).unwrap();
    // knots_v = [0,0, 0.5, 0.5, 1,1], index 2 is the first inserted 0.5
    let (result, count) = inserted.remove_knot_v_times(2, 3, 0.0);
    assert_eq!(count, 2, "Should have removed 2 knots");
    assert_eq!(result.knots_v.len(), inserted.knots_v.len() - 2);
}

#[test]
fn surface_elevate_degree_u_preserves_shape() {
    let surface = bilinear_patch();
    let elevated = surface.elevate_degree_u();
    assert_eq!(elevated.degree_u, 2);
    for u in [0.0f64, 0.25, 0.5, 0.75, 0.99] {
        for v in [0.0f64, 0.25, 0.5, 0.75, 0.99] {
            let p1 = surface.point_at(u, v);
            let p2 = elevated.point_at(u, v);
            assert!(
                (p1.x - p2.x).abs() < 1e-10,
                "x mismatch at u={u},v={v}: {p1:?} vs {p2:?}"
            );
            assert!(
                (p1.y - p2.y).abs() < 1e-10,
                "y mismatch at u={u},v={v}: {p1:?} vs {p2:?}"
            );
        }
    }
}

#[test]
fn surface_elevate_degree_v_preserves_shape() {
    let surface = bilinear_patch();
    let elevated = surface.elevate_degree_v();
    assert_eq!(elevated.degree_v, 2);
    for u in [0.0f64, 0.25, 0.5, 0.75, 0.99] {
        for v in [0.0f64, 0.25, 0.5, 0.75, 0.99] {
            let p1 = surface.point_at(u, v);
            let p2 = elevated.point_at(u, v);
            assert!(
                (p1.x - p2.x).abs() < 1e-10,
                "x mismatch at u={u},v={v}: {p1:?} vs {p2:?}"
            );
            assert!(
                (p1.y - p2.y).abs() < 1e-10,
                "y mismatch at u={u},v={v}: {p1:?} vs {p2:?}"
            );
        }
    }
}

#[test]
fn surface_elevate_degree_u_to_target() {
    let surface = bilinear_patch();
    let elevated = surface.elevate_degree_u_to(3).unwrap();
    assert_eq!(elevated.degree_u, 3);
}

#[test]
fn surface_elevate_degree_v_to_target() {
    let surface = bilinear_patch();
    let elevated = surface.elevate_degree_v_to(3).unwrap();
    assert_eq!(elevated.degree_v, 3);
}

#[test]
fn surface_elevate_degree_u_same_fails() {
    let surface = bilinear_patch();
    assert!(surface.elevate_degree_u_to(1).is_err());
    assert!(surface.elevate_degree_u_to(0).is_err());
}

#[test]
fn surface_elevate_degree_v_same_fails() {
    let surface = bilinear_patch();
    assert!(surface.elevate_degree_v_to(1).is_err());
    assert!(surface.elevate_degree_v_to(0).is_err());
}

#[test]
fn surface_split_at_u_preserves_shape() {
    let surface = bilinear_patch();
    let (left, right) = surface.split_at_u(0.5).unwrap();

    for u in [0.1, 0.2, 0.3, 0.4] {
        for v in [0.0, 0.5, 0.99] {
            let po = surface.point_at(u, v);
            let pl = left.point_at(u, v);
            assert!(
                (po.x - pl.x).abs() < 1e-8,
                "left x mismatch at u={u},v={v}: {po:?} vs {pl:?}"
            );
            assert!(
                (po.y - pl.y).abs() < 1e-8,
                "left y mismatch at u={u},v={v}: {po:?} vs {pl:?}"
            );
        }
    }
    for u in [0.6, 0.7, 0.8, 0.9] {
        for v in [0.0, 0.5, 0.99] {
            let po = surface.point_at(u, v);
            let pr = right.point_at(u, v);
            assert!(
                (po.x - pr.x).abs() < 1e-8,
                "right x mismatch at u={u},v={v}: {po:?} vs {pr:?}"
            );
            assert!(
                (po.y - pr.y).abs() < 1e-8,
                "right y mismatch at u={u},v={v}: {po:?} vs {pr:?}"
            );
        }
    }
}

#[test]
fn surface_split_at_v_preserves_shape() {
    let surface = bilinear_patch();
    let (bottom, top) = surface.split_at_v(0.5).unwrap();

    for u in [0.0, 0.5, 0.99] {
        for v in [0.1, 0.2, 0.3, 0.4] {
            let po = surface.point_at(u, v);
            let pb = bottom.point_at(u, v);
            assert!(
                (po.x - pb.x).abs() < 1e-8,
                "bottom x mismatch at u={u},v={v}: {po:?} vs {pb:?}"
            );
            assert!(
                (po.y - pb.y).abs() < 1e-8,
                "bottom y mismatch at u={u},v={v}: {po:?} vs {pb:?}"
            );
        }
    }
    for u in [0.0, 0.5, 0.99] {
        for v in [0.6, 0.7, 0.8, 0.9] {
            let po = surface.point_at(u, v);
            let pt = top.point_at(u, v);
            assert!(
                (po.x - pt.x).abs() < 1e-8,
                "top x mismatch at u={u},v={v}: {po:?} vs {pt:?}"
            );
            assert!(
                (po.y - pt.y).abs() < 1e-8,
                "top y mismatch at u={u},v={v}: {po:?} vs {pt:?}"
            );
        }
    }
}

#[test]
fn surface_split_at_u_endpoint_fails() {
    let surface = bilinear_patch();
    assert!(surface.split_at_u(0.0).is_err());
    assert!(surface.split_at_u(1.0).is_err());
}

#[test]
fn surface_split_at_v_endpoint_fails() {
    let surface = bilinear_patch();
    assert!(surface.split_at_v(0.0).is_err());
    assert!(surface.split_at_v(1.0).is_err());
}

#[test]
fn surface_insert_knots_u_out_of_range() {
    let surface = bilinear_patch();
    assert!(surface.insert_knots_u(&[0.5, 2.0]).is_err());
    assert!(surface.insert_knots_u(&[-0.5, 0.5]).is_err());
}

#[test]
fn surface_insert_knots_v_out_of_range() {
    let surface = bilinear_patch();
    assert!(surface.insert_knots_v(&[0.5, 2.0]).is_err());
    assert!(surface.insert_knots_v(&[-0.5, 0.5]).is_err());
}

#[test]
fn surface_remove_knot_u_times_zero_removals() {
    let surface = bilinear_patch();
    // No interior knots to remove -> 0 removals
    let (result, count) = surface.remove_knot_u_times(0, 1, 0.0);
    assert_eq!(count, 0);
    assert_eq!(result.num_u, surface.num_u);
}

#[test]
fn surface_remove_knot_v_times_zero_removals() {
    let surface = bilinear_patch();
    let (result, count) = surface.remove_knot_v_times(0, 1, 0.0);
    assert_eq!(count, 0);
    assert_eq!(result.num_v, surface.num_v);
}

/// Test that remove_knot_u_times stops when knot_index >= knots_u.len() (line 271).
/// We insert a knot at a high index and try to remove it multiple times.
#[test]
fn surface_remove_knot_u_times_index_out_of_range_break() {
    let surface = bilinear_patch();
    // Insert a single knot at 0.5
    let inserted = surface.insert_knot_u(0.5).unwrap();
    // knots_u = [0, 0, 0.5, 1, 1]. Index 2 is the inserted knot.
    // Try to remove at index 2 ten times — after the first removal the
    // index will be past the end, triggering the break.
    let (result, count) = inserted.remove_knot_u_times(2, 10, 0.0);
    assert_eq!(count, 1);
    assert_eq!(result.num_u, surface.num_u);
}

/// Same test for v direction (line 297).
#[test]
fn surface_remove_knot_v_times_index_out_of_range_break() {
    let surface = bilinear_patch();
    let inserted = surface.insert_knot_v(0.5).unwrap();
    let (result, count) = inserted.remove_knot_v_times(2, 10, 0.0);
    assert_eq!(count, 1);
    assert_eq!(result.num_v, surface.num_v);
}
