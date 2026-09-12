//! Unit gates for the triple-point junction refinement (PK-480).

use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle};

use super::refine_to_triple_point;

/// The unit cube's vertical edge at (1,1): both endpoints carry
/// x = y = 1 and the two adjacent faces are the x=1 and y=1 planes.
fn vertical_edge_at_1_1(brep: &BRep) -> EdgeHandle {
    for (eh, e) in brep.edges.iter() {
        let Some(he) = brep.half_edges.get(e.half_edge) else {
            continue;
        };
        let a = brep.vertices[he.origin].position;
        let b = brep.vertices[brep.half_edges[he.twin].origin].position;
        if (a.x - 1.0).abs() < 1e-9
            && (b.x - 1.0).abs() < 1e-9
            && (a.y - 1.0).abs() < 1e-9
            && (b.y - 1.0).abs() < 1e-9
        {
            return eh;
        }
    }
    panic!("unit cube lacks the vertical edge at (1,1)");
}

/// An oblique plane crossing that edge: the refinement must pull a
/// slid seed onto the exact three-plane point, and keep a seed
/// perturbed past the move cap unchanged.
#[test]
fn refines_onto_the_three_plane_point_and_respects_the_cap() {
    let tolerance = Tolerance::default();
    let brep = BRep::unit_cube(tolerance).expect("cube");
    let eh = vertical_edge_at_1_1(&brep);
    // Plane z = 0.4 + 0.05·x + 0.03·y: the true triple point with the
    // x=1 and y=1 faces is (1, 1, 0.48).
    let normal = Vec3::new(-0.05, -0.03, 1.0).normalized().expect("n");
    let u_axis = normal
        .cross(&Vec3::new(0.0, 1.0, 0.0))
        .normalized()
        .expect("u");
    let v_axis = normal.cross(&u_axis);
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.4),
        normal,
        u_axis,
        v_axis,
    };
    let exact = Point3::new(1.0, 1.0, 0.48);
    // A seed slid 5e-3 along the edge — the amplified near-tangent
    // crossing class this refinement retires.
    let seed = Point3::new(1.0, 1.0, 0.48 + 5.0e-3);
    let refined = refine_to_triple_point(&brep, eh, &plane, seed, &tolerance);
    assert!(
        refined.distance_to(&exact) < 1e-9,
        "refined {refined:?} should land on {exact:?}"
    );
    // A seed past the move cap must come back unchanged.
    let far = Point3::new(1.0, 1.0, 0.98);
    let kept = refine_to_triple_point(&brep, eh, &plane, far, &tolerance);
    assert_eq!(kept, far, "past-cap seed must be kept");
}

/// A plane PARALLEL to one adjacent face makes the normal triple
/// singular: the seed must be kept, never a wild solve.
#[test]
fn degenerate_normal_triple_keeps_the_seed() {
    let tolerance = Tolerance::default();
    let brep = BRep::unit_cube(tolerance).expect("cube");
    let eh = vertical_edge_at_1_1(&brep);
    let plane = SurfaceGeometry::Plane {
        origin: Point3::new(0.3, 0.0, 0.0),
        normal: Vec3::new(1.0, 0.0, 0.0),
        u_axis: Vec3::new(0.0, 1.0, 0.0),
        v_axis: Vec3::new(0.0, 0.0, 1.0),
    };
    let seed = Point3::new(1.0, 1.0, 0.5);
    let kept = refine_to_triple_point(&brep, eh, &plane, seed, &tolerance);
    assert_eq!(kept, seed);
}
