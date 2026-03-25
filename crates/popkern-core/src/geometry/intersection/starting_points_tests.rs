use super::*;
use crate::geometry::Point3;

fn xy_plane() -> SurfaceGeometry {
    SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    }
}

fn xz_plane() -> SurfaceGeometry {
    SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 1.0, 0.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 0.0, 1.0),
    }
}

#[test]
fn find_starting_points_two_planes() {
    let s1 = xy_plane();
    let s2 = xz_plane();
    let points = find_starting_points(
        &s1,
        (-1.0, 1.0),
        (-1.0, 1.0),
        &s2,
        (-1.0, 1.0),
        (-1.0, 1.0),
        &Tolerance::default(),
    );
    assert!(!points.is_empty(), "Should find starting points");
    assert!(
        !points[0].is_tangent,
        "XY/XZ intersection should be transversal"
    );
}

#[test]
fn no_starting_points_disjoint() {
    let s1 = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    let s2 = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 10.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    let points = find_starting_points(
        &s1,
        (0.0, 1.0),
        (0.0, 1.0),
        &s2,
        (0.0, 1.0),
        (0.0, 1.0),
        &Tolerance::default(),
    );
    assert!(
        points.is_empty(),
        "Disjoint parallel planes should yield no starting points"
    );
}

#[test]
fn tangent_classification() {
    // Two identical planes should produce tangent starting points (or not panic).
    let s = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    // For identical surfaces Newton may not converge (over-determined / under-constrained).
    // The key requirement is that this does not panic.
    let _ = find_starting_points(
        &s,
        (0.0, 1.0),
        (0.0, 1.0),
        &s,
        (0.0, 1.0),
        (0.0, 1.0),
        &Tolerance::default(),
    );
}

#[test]
fn classify_tangent_parallel_planes() {
    // Two parallel planes at the same height — normals are identical.
    let s = xy_plane();
    let sample = MarchSample {
        u1: 0.5,
        v1: 0.5,
        u2: 0.5,
        v2: 0.5,
        point: Point3::new(0.5, 0.5, 0.0),
    };
    let tol = Tolerance::default();
    assert!(
        classify_tangent(&s, &s, &sample, &tol),
        "Identical-normal planes should be classified tangent"
    );
}

#[test]
fn classify_transversal_perpendicular_planes() {
    let s1 = xy_plane();
    let s2 = xz_plane();
    let sample = MarchSample {
        u1: 0.0,
        v1: 0.0,
        u2: 0.0,
        v2: 0.0,
        point: Point3::new(0.0, 0.0, 0.0),
    };
    let tol = Tolerance::default();
    assert!(
        !classify_tangent(&s1, &s2, &sample, &tol),
        "Perpendicular planes should be classified transversal"
    );
}

#[test]
fn find_branch_points_does_not_panic() {
    // Basic smoke test: the function must return without panicking regardless of whether
    // any branches are found.
    let s1 = xy_plane();
    let s2 = xz_plane();
    let start = MarchSample {
        u1: 0.0,
        v1: 0.0,
        u2: 0.0,
        v2: 0.0,
        point: Point3::new(0.0, 0.0, 0.0),
    };
    let tol = Tolerance::default();
    let _ = find_branch_points(
        &s1,
        (-1.0, 1.0),
        (-1.0, 1.0),
        &s2,
        (-1.0, 1.0),
        (-1.0, 1.0),
        &[],
        &start,
        &tol,
    );
}

#[test]
fn find_branch_points_excludes_existing_curve() {
    // Candidates that are on the existing traced curve must be filtered out.
    let s1 = xy_plane();
    let s2 = xz_plane();
    let start = MarchSample {
        u1: 0.0,
        v1: 0.0,
        u2: 0.0,
        v2: 0.0,
        point: Point3::new(0.0, 0.0, 0.0),
    };
    let tol = Tolerance::default();

    // Run once with no traced samples to collect candidates.
    let raw_branches = find_branch_points(
        &s1,
        (-1.0, 1.0),
        (-1.0, 1.0),
        &s2,
        (-1.0, 1.0),
        (-1.0, 1.0),
        &[],
        &start,
        &tol,
    );

    // Now supply those same candidates as already-traced; result must be empty.
    let filtered = find_branch_points(
        &s1,
        (-1.0, 1.0),
        (-1.0, 1.0),
        &s2,
        (-1.0, 1.0),
        (-1.0, 1.0),
        &raw_branches,
        &start,
        &tol,
    );
    assert!(
        filtered.is_empty(),
        "All previously-found branches should be suppressed when supplied as traced_samples"
    );
}

#[test]
fn find_branch_points_no_branches_for_disjoint_surfaces() {
    // Disjoint parallel planes should produce no branch candidates.
    let s1 = xy_plane();
    let s2 = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 10.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    };
    let start = MarchSample {
        u1: 0.5,
        v1: 0.5,
        u2: 0.5,
        v2: 0.5,
        point: Point3::new(0.5, 0.5, 0.0),
    };
    let tol = Tolerance::default();
    let branches = find_branch_points(
        &s1,
        (0.0, 1.0),
        (0.0, 1.0),
        &s2,
        (0.0, 1.0),
        (0.0, 1.0),
        &[],
        &start,
        &tol,
    );
    // Newton projection should fail for disjoint surfaces, yielding no branches.
    assert!(
        branches.is_empty(),
        "Disjoint surfaces should yield no branch candidates"
    );
}

#[test]
fn queue_cap_prevents_oom_on_overlapping_surfaces() {
    use crate::geometry::nurbs_surface::NurbsSurface;
    let nurbs = SurfaceGeometry::Nurbs(NurbsSurface {
        degree_u: 1,
        degree_v: 1,
        knots_u: vec![0.0, 0.0, 1.0, 1.0],
        knots_v: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
        ],
        num_u: 2,
        num_v: 2,
        weights: vec![1.0; 4],
    });
    let r = find_starting_points(
        &nurbs,
        (0.0, 1.0),
        (0.0, 1.0),
        &nurbs,
        (0.0, 1.0),
        (0.0, 1.0),
        &Tolerance::default(),
    );
    assert!(r.len() <= MAX_CANDIDATE_PAIRS);
}
