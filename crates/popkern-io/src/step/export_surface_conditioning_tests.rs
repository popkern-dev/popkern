//! Conditioned Gregory export: node clustering, boundary curvature
//! measurement, and the curvature acceptance on an ill-conditioned
//! rational patch (PK-551).

use std::f64::consts::SQRT_2;

use super::*;

/// A quarter cylinder of radius `r` (u around, v along +z) as the exact
/// bi-quadratic rational patch.
fn quarter_cylinder(r: f64) -> SurfaceGeometry {
    let w = SQRT_2 / 2.0;
    let ring = [(r, 0.0, 1.0), (r, r, w), (0.0, r, 1.0)];
    let mut control_points = Vec::new();
    let mut weights = Vec::new();
    // u-major storage (index = i_u * num_v + j_v): u runs around the
    // rim, v along the ruling.
    for &(x, y, wu) in &ring {
        for z in [0.0, 1.0, 2.0] {
            control_points.push(Point3::new(x, y, z));
            weights.push(wu);
        }
    }
    SurfaceGeometry::Nurbs(NurbsSurface {
        degree_u: 2,
        degree_v: 2,
        knots_u: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
        knots_v: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
        control_points,
        num_u: 3,
        num_v: 3,
        weights,
    })
}

#[test]
fn chebyshev_lobatto_nodes_include_the_ends_and_cluster_there() {
    let nodes = clustered(0.0, 2.0, 33);
    assert_eq!(nodes.len(), 33);
    assert!(nodes[0].abs() < 1e-15 && (nodes[32] - 2.0).abs() < 1e-15);
    assert!(
        (nodes[16] - 1.0).abs() < 1e-12,
        "symmetric about the middle"
    );
    assert!(
        nodes[1] - nodes[0] < 0.2 * (nodes[17] - nodes[16]),
        "clustered at the ends"
    );
}

#[test]
fn boundary_curvatures_read_the_cylinder() {
    let cylinder = quarter_cylinder(2.0);
    // Along the rim (u) the normal curvature is ±1/r; along the ruling
    // (v) it is 0.
    let [along, across] =
        boundary_curvatures(&cylinder, 0.3, 0.0, true).unwrap_or_else(|| panic!("regular point"));
    assert!((along.abs() - 0.5).abs() < 1e-12, "{along}");
    assert!(across.abs() < 1e-12, "{across}");
    let [along_v, across_v] =
        boundary_curvatures(&cylinder, 0.0, 0.4, false).unwrap_or_else(|| panic!("regular point"));
    assert!(along_v.abs() < 1e-12 && (across_v.abs() - 0.5).abs() < 1e-12);
}

#[test]
fn the_conditioned_fit_holds_boundary_curvature_where_a_uniform_fit_does_not() {
    // A rational quarter cylinder: its angle parameterization is not
    // polynomial, so a non-rational fit carries real curvature error.
    let patch = quarter_cylinder(2.0);
    let fit = conditioned_fit(&patch, (0.0, 1.0), (0.0, 1.0), 1e-7)
        .unwrap_or_else(|| panic!("a position-valid fit exists"));
    assert!(fit.weights.iter().all(|&w| w == 1.0), "non-rational");
    let boundary = boundary_samples(&patch, (0.0, 1.0), (0.0, 1.0));
    let bar = CURVATURE_FIDELITY * boundary.curvature_scale;
    let error = boundary.curvature_error(&fit);
    assert!(error <= bar, "boundary curvature error {error:e} > {bar:e}");
    // Negative control: the uniform-node, position-only fit this
    // replaced (17 nodes met the position bar) misses the curvature bar.
    let uniform: Vec<f64> = (0..17).map(|k| k as f64 / 16.0).collect();
    let grid: Vec<Vec<Point3>> = uniform
        .iter()
        .map(|&v| uniform.iter().map(|&u| patch.point_at(u, v)).collect())
        .collect();
    let coarse = popkern_core::geometry::fitting::interpolate_surface_grid(
        &grid, &uniform, &uniform, FIT_DEGREE, FIT_DEGREE,
    )
    .unwrap_or_else(|e| panic!("{e}"));
    assert!(boundary.curvature_error(&coarse) > bar);
}

#[test]
fn a_well_conditioned_form_is_written_exactly() {
    let patch = quarter_cylinder(2.0);
    let written = conditioned_gregory_nurbs(&patch, 1e-7);
    assert_eq!(
        written.weights,
        patch.to_nurbs().weights,
        "the exact form is kept"
    );
}
