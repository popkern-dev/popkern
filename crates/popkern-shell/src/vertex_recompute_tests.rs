use super::*;
use crate::classify::classify_faces;
use crate::inner_skin::build_inner_skin;
use crate::surface_offset::compute_offset_surfaces;
use crate::types::ShellDirection;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

fn tol() -> Tolerance {
    Tolerance::default()
}

#[test]
fn cube_vertices_shift_by_thickness() {
    let brep = BRep::unit_cube(tol()).unwrap();
    let thickness = 0.1;
    let classification =
        classify_faces(&brep, thickness, &[], ShellDirection::Inward, &[]).unwrap();
    let offset_result = compute_offset_surfaces(&brep, &classification, true, &tol()).unwrap();
    let (mut inner, skin_map) = build_inner_skin(&brep, &[]).unwrap();

    recompute_vertices(
        &mut inner,
        &offset_result.surfaces,
        &skin_map,
        &brep,
        &tol(),
    )
    .unwrap();

    // Each vertex of the inner cube should have moved from its original
    // position. Verify all vertices actually changed.
    let reverse_vertex = skin_map.reverse_vertex_map();
    let mut moved_count = 0;
    for (inner_vh, v) in &inner.vertices {
        if let Some(&outer_vh) = reverse_vertex.get(&inner_vh) {
            let outer_pos = brep.vertices[outer_vh].position;
            let dist = outer_pos.distance_to(&v.position);
            if dist > 1e-8 {
                moved_count += 1;
            }
            // The offset distance along the vertex diagonal of a cube corner
            // (where 3 orthogonal planes meet) should be thickness * sqrt(3).
            // This holds regardless of normal direction.
            let expected_diag = thickness * 3.0_f64.sqrt();
            assert!(
                (dist - expected_diag).abs() < 0.05,
                "vertex moved by {dist}, expected ~{expected_diag}"
            );
        }
    }
    assert_eq!(moved_count, 8, "all 8 vertices should have moved");
}

#[test]
fn solve_3_planes_basic() {
    // Three coordinate planes offset by 1.
    let p0 = (Vec3::new(1.0, 0.0, 0.0), 1.0); // x = 1
    let p1 = (Vec3::new(0.0, 1.0, 0.0), 2.0); // y = 2
    let p2 = (Vec3::new(0.0, 0.0, 1.0), 3.0); // z = 3
    let result = solve_3_planes(&p0, &p1, &p2).unwrap();
    assert!((result.x - 1.0).abs() < 1e-12);
    assert!((result.y - 2.0).abs() < 1e-12);
    assert!((result.z - 3.0).abs() < 1e-12);
}

#[test]
fn solve_3_planes_singular() {
    // Two parallel planes + one perpendicular → singular.
    let p0 = (Vec3::new(0.0, 0.0, 1.0), 1.0);
    let p1 = (Vec3::new(0.0, 0.0, 1.0), 2.0);
    let p2 = (Vec3::new(1.0, 0.0, 0.0), 0.0);
    let result = solve_3_planes(&p0, &p1, &p2);
    assert!(result.is_none());
}

#[test]
fn two_plane_vertex_solve_is_exact() {
    let planes = vec![
        (Vec3::new(1.0, 0.0, 0.0), 1.0), // x = 1
        (Vec3::new(0.0, 1.0, 0.0), 2.0), // y = 2
    ];
    let seed = Point3::new(0.0, 0.0, 0.7);
    let result = solve_underdetermined_planes(&planes, seed).unwrap();
    assert!((result.x - 1.0).abs() < 1e-14);
    assert!((result.y - 2.0).abs() < 1e-14);
    assert!((result.z - 0.7).abs() < 1e-14, "seed's line station kept");
}

#[test]
fn near_parallel_planes_collapse_to_one_constraint() {
    let planes = vec![
        (Vec3::new(0.0, 0.0, 1.0), 1.0),
        (Vec3::new(0.0, 0.0, 1.0), 1.0),
    ];
    let seed = Point3::new(0.3, 0.4, 0.0);
    let result = solve_underdetermined_planes(&planes, seed).unwrap();
    assert!((result.z - 1.0).abs() < 1e-14);
    assert!((result.x - 0.3).abs() < 1e-14);
}
