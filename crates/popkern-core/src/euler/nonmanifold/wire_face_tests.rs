use super::*;
use crate::geometry::surface::{SurfaceGeometry, SurfaceOnFace};
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;
use crate::topology::{BodyKind, HalfEdgeHandle, RegionHandle};

fn tol() -> Tolerance {
    Tolerance::default()
}

fn assert_valid(body: &BRep, what: &str) {
    let errors = body.validate_all();
    assert!(errors.is_empty(), "{what}: {errors:?}");
}

/// A square wire cycle at height `z` over [x0, x0 + s]², counter-
/// clockwise about +z, in `region`. Returns the cycle's fins.
fn square(body: &mut BRep, region: RegionHandle, x0: f64, s: f64, z: f64) -> Vec<HalfEdgeHandle> {
    let corners = [(x0, x0), (x0 + s, x0), (x0 + s, x0 + s), (x0, x0 + s)];
    let v: Vec<_> = corners
        .iter()
        .map(|&(x, y)| {
            body.make_acorn(region, Point3::new(x, y, z))
                .expect("acorn")
                .vertex
        })
        .collect();
    (0..4)
        .map(|i| {
            body.make_wire_edge(region, v[i], v[(i + 1) % 4])
                .expect("wire")
                .fin
        })
        .collect()
}

fn plane_at(z: f64) -> SurfaceOnFace {
    SurfaceOnFace {
        geometry: SurfaceGeometry::Plane {
            origin: Point3::new(0.0, 0.0, z),
            normal: Vec3::new(0.0, 0.0, 1.0),
            u_axis: Vec3::new(1.0, 0.0, 0.0),
            v_axis: Vec3::new(0.0, 1.0, 0.0),
        },
        orientation: true,
    }
}

/// A face spanning a wire square in the outside is a one-face sheet on
/// free edges; turning it back gives the wire square again.
#[test]
fn a_wire_square_spans_a_sheet_and_back() {
    let mut body = BRep::new(tol());
    let out = body.outside;
    let cycle = square(&mut body, out, 0.0, 1.0, 0.0);
    assert_eq!(body.body_kind(), BodyKind::Wire);
    let wires = NonManifoldCounts::of(&body);
    let made = body
        .make_face_from_wire_cycle(&cycle, Some(plane_at(0.0)))
        .expect("face");
    assert_eq!(
        wires.delta(NonManifoldCounts::of(&body)),
        [0, 0, -4, 1, 1, 0, 0, -4, 0]
    );
    body.recompute_face_pcurves(made.face).expect("pcurves");
    assert_eq!(body.body_kind(), BodyKind::Sheet);
    assert!(body.edges.keys().all(|e| body.is_free_edge(e)));
    assert_valid(&body, "sheet");
    let edges = body.kill_face_to_wire_cycle(made.face).expect("back");
    assert_eq!(edges.len(), 4);
    assert_eq!(NonManifoldCounts::of(&body), wires);
    assert_eq!(body.body_kind(), BodyKind::Wire);
    assert_valid(&body, "wires again");
}

/// A wire square inside a cube's solid region spans a lamina in the
/// solid: both sides of the face in the lamina's own shell of that
/// region, a general body.
#[test]
fn a_lamina_inside_a_solid_is_a_general_body() {
    let mut cube = BRep::unit_cube(tol()).expect("cube");
    let solid = cube.solid_regions()[0];
    let cycle = square(&mut cube, solid, 0.25, 0.5, 0.5);
    let made = cube
        .make_face_from_wire_cycle(&cycle, Some(plane_at(0.5)))
        .expect("lamina");
    cube.recompute_face_pcurves(made.face).expect("pcurves");
    let sides = cube.face_sides(made.face).expect("sides");
    assert_eq!(sides.back, sides.front);
    assert_eq!(cube.shell_region(sides.back), Some(solid));
    assert_eq!(cube.body_kind(), BodyKind::General);
    assert_valid(&cube, "lamina in a cube");
}

/// Declines leave the body unchanged: a cycle that does not close, a
/// fin listed twice, a cube's face (not a lamina on free edges).
#[test]
fn wire_faces_decline_without_changing_the_body() {
    let mut body = BRep::new(tol());
    let out = body.outside;
    let cycle = square(&mut body, out, 0.0, 1.0, 0.0);
    let before = NonManifoldCounts::of(&body);
    assert!(matches!(
        body.make_face_from_wire_cycle(&cycle[..3], None),
        Err(EulerError::NotAWireCycle { .. })
    ));
    assert!(matches!(
        body.make_face_from_wire_cycle(&[cycle[0], cycle[1], cycle[0]], None),
        Err(EulerError::NotAWireCycle { .. })
    ));
    assert_eq!(NonManifoldCounts::of(&body), before);
    let mut cube = BRep::unit_cube(tol()).expect("cube");
    let face = cube.faces.keys().next().expect("face");
    assert_eq!(
        cube.kill_face_to_wire_cycle(face),
        Err(EulerError::NotALaminaFace { face })
    );
}
