//! Pins for the row shading contract (PK-914): a constant-curvature
//! round departs by nothing (to third order) in any chart, a cubic
//! profile's rail row departs as its curvature rises (`h²`), tilted
//! vertex normals read their tilt, and the display pass narrows an
//! over-bound rail row until it keeps the bound.

use super::*;
use popkern_core::geometry::Point3;
use popkern_core::geometry::nurbs_surface::NurbsSurface;

use crate::mesh::FaceMesh;

/// A radius-1 cylinder of height 1 as a rational NURBS (u the angle,
/// not arc length within each quarter: the chart is non-uniform).
fn rational_tube() -> SurfaceGeometry {
    let w = std::f64::consts::FRAC_1_SQRT_2;
    let square = [(1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (-1.0, 1.0), (-1.0, 0.0)];
    let (mut control_points, mut weights) = (Vec::new(), Vec::new());
    for (k, (x, y)) in square.iter().enumerate() {
        for z in [0.0, 1.0] {
            control_points.push(Point3::new(*x, *y, z));
            weights.push(if k % 2 == 0 { 1.0 } else { w });
        }
    }
    use std::f64::consts::{FRAC_PI_2, PI};
    SurfaceGeometry::Nurbs(NurbsSurface {
        degree_u: 2,
        degree_v: 1,
        knots_u: vec![0.0, 0.0, 0.0, FRAC_PI_2, FRAC_PI_2, PI, PI, PI],
        knots_v: vec![0.0, 0.0, 1.0, 1.0],
        control_points,
        num_u: 5,
        num_v: 2,
        weights,
    })
}

/// The profile `y = x³` over `x` in [0, 1] extruded along z over
/// [0, 1]: flat at `x = 0` (a G2 rail), its curvature rising as `6x`.
fn cubic_strip() -> SurfaceGeometry {
    let (mut control_points, mut weights) = (Vec::new(), Vec::new());
    for (x, y) in [(0.0, 0.0), (1.0 / 3.0, 0.0), (2.0 / 3.0, 0.0), (1.0, 1.0)] {
        for z in [0.0, 1.0] {
            control_points.push(Point3::new(x, y, z));
            weights.push(1.0);
        }
    }
    SurfaceGeometry::Nurbs(NurbsSurface {
        degree_u: 3,
        degree_v: 1,
        knots_u: vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0],
        knots_v: vec![0.0, 0.0, 1.0, 1.0],
        control_points,
        num_u: 4,
        num_v: 2,
        weights,
    })
}

fn vertex(surface: &SurfaceGeometry, u: f64, v: f64) -> MeshVertex {
    let n = surface.normal_at(u, v);
    let l = n.length();
    MeshVertex {
        position: surface.point_at(u, v),
        normal: Vec3::new(n.x / l, n.y / l, n.z / l),
        uv: (u, v),
    }
}

/// A needle on the base `u = u0`, `v` in [0, 1], apex at `u0 + du`,
/// `v = 0.5`, at the display's default chord height.
fn needle_error(surface: &SurfaceGeometry, u0: f64, du: f64) -> Option<f64> {
    let (a, b, c) = (
        vertex(surface, u0, 0.0),
        vertex(surface, u0, 1.0),
        vertex(surface, u0 + du, 0.5),
    );
    row_shading_error(surface, &a, &b, &c, 0.005)
}

#[test]
fn a_constant_curvature_round_departs_by_nothing_in_a_rational_chart() {
    let tube = rational_tube();
    // Rows 0.05 to 0.2 rad across, near a knot and mid-span, where
    // the rational chart's angle-per-u differs.
    for (u0, du) in [(0.1, 0.05), (0.7, 0.2), (1.4, 0.15), (2.2, 0.1)] {
        let e = needle_error(&tube, u0, du);
        let e = e.unwrap_or_else(|| panic!("u0 {u0} du {du}: not measured"));
        // Third order in the row's turn (3.9e-4 read at du 0.2),
        // far under the 6.5e-3 bound.
        assert!(e <= du * du * du / 16.0, "u0 {u0} du {du}: {e:.3e}");
    }
}

#[test]
fn a_rail_row_departs_as_its_curvature_rises() {
    let strip = cubic_strip();
    // The normal turns by atan(3x²) from a flat rail: a row of height h
    // departs by about 3h², a quarter per halving.
    let wide = needle_error(&strip, 0.0, 0.08).expect("wide row");
    let narrow = needle_error(&strip, 0.0, 0.04).expect("narrow row");
    assert!((wide / 1.92e-2 - 1.0).abs() < 0.1, "wide {wide:.3e}");
    assert!(
        (wide / narrow - 4.0).abs() < 0.4,
        "{wide:.3e} / {narrow:.3e}"
    );
    // The contract holds the narrow row and not the wide one.
    let bound = row_shading_bound(crate::DISPLAY_MAX_FACET_ANGLE);
    assert!(wide > bound && narrow < bound, "bound {bound:.3e}");
}

#[test]
fn a_tilted_apex_normal_reads_its_tilt() {
    // The negative control: the same round with its apex normal turned
    // 0.01 rad about the ruling reads the turn, so the zero above is
    // the surface's, not the metric's blindness.
    let tube = rational_tube();
    let (a, b, mut c) = (
        vertex(&tube, 0.7, 0.0),
        vertex(&tube, 0.7, 1.0),
        vertex(&tube, 0.9, 0.5),
    );
    let (n, t) = (c.normal, Vec3::new(0.0, 0.0, 1.0).cross(&c.normal));
    let (s, k) = (0.01_f64.psin(), 0.01_f64.pcos());
    c.normal = n * k + t * s;
    let e = row_shading_error(&tube, &a, &b, &c, 0.005).expect("measurable");
    assert!((e - 0.01).abs() < 5e-4, "{e:.3e}");
}

#[test]
fn compact_facets_apexes_beyond_the_base_and_sub_resolution_rows_are_not_held() {
    let strip = cubic_strip();
    let v = |u, w| vertex(&strip, u, w);
    // Base 0.2 under a 0.1 height: compact.
    assert_eq!(
        row_shading_error(&strip, &v(0.0, 0.4), &v(0.0, 0.6), &v(0.1, 0.5), 0.005),
        None
    );
    // Apex beyond the base's end: it lies along the base.
    assert_eq!(
        row_shading_error(&strip, &v(0.0, 0.1), &v(0.0, 0.5), &v(0.02, 0.8), 0.005),
        None
    );
    // A row 1e-4 high, under a sixteenth of the chord height.
    assert_eq!(needle_error(&strip, 0.0, 1e-4), None);
}

#[test]
fn the_display_pass_narrows_an_over_bound_rail_row() {
    // One row of needles 0.1 across on the cubic strip's flat rail:
    // stations z in {0, 0.5, 1} at x = 0 and x = 0.1.
    let strip = cubic_strip();
    let vertices: Vec<MeshVertex> = [(0.0, 0.0), (0.0, 0.5), (0.0, 1.0)]
        .iter()
        .chain([(0.1, 0.0), (0.1, 0.5), (0.1, 1.0)].iter())
        .map(|&(u, v)| vertex(&strip, u, v))
        .collect();
    let triangles = vec![[0, 3, 4], [0, 4, 1], [1, 4, 5], [1, 5, 2]];
    let mut mesh = FaceMesh {
        face: Default::default(),
        vertices,
        triangles,
    };
    let max_angle = crate::DISPLAY_MAX_FACET_ANGLE;
    let left = crate::angular_refine::refine_face_shading(&mut mesh, &strip, max_angle, 0.005);
    assert_eq!(left, 0);
    // Every facet on the rail now stands at most 0.05 off it (two
    // halvings of 0.1 meet the bound: 3h² <= 6.5e-3 at h <= 0.047),
    // and the pass stayed local.
    let mut rail_rows = 0;
    for t in &mesh.triangles {
        if t.iter().filter(|&&i| mesh.vertices[i].uv.0 == 0.0).count() == 2 {
            let apex = t.iter().map(|&i| mesh.vertices[i].uv.0).fold(0.0, f64::max);
            assert!(apex <= 0.05 + 1e-12, "rail row {apex}");
            rail_rows += 1;
        }
    }
    assert_eq!(rail_rows, 2);
    assert!(
        mesh.triangles.len() <= 64,
        "{} facets",
        mesh.triangles.len()
    );
}
