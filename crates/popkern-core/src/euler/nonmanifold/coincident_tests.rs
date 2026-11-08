use super::super::NonManifoldCounts;
use super::*;
use crate::geometry::transform::Transform3;
use crate::tolerance::Tolerance;
use crate::topology::BodyKind;

fn tol() -> Tolerance {
    Tolerance::default()
}

fn cube_at(x: f64, y: f64, z: f64) -> BRep {
    BRep::unit_cube(tol())
        .and_then(|c| c.transform(&Transform3::translation(x, y, z)))
        .expect("cube")
}

/// Two cubes touching along an edge, nothing shared: the glue makes one
/// vertex at each end and one ring of four, a general body that
/// validates.
#[test]
fn cubes_touching_along_an_edge_glue_into_one_body() {
    let mut body = cube_at(0.0, 0.0, 0.0);
    body.append_body(&cube_at(1.0, 1.0, 0.0));
    let glued = body.glue_coincident().expect("glue");
    assert_eq!(
        glued,
        CoincidentGlue {
            vertices: 2,
            edges: 1
        }
    );
    assert_eq!(body.body_kind(), BodyKind::General);
    assert!(body.edges.keys().any(|e| body.radial_fins(e).len() == 4));
    let errors = body.validate_all();
    assert!(errors.is_empty(), "{errors:?}");
}

/// Negative control: a valid body has nothing coincident to glue.
#[test]
fn a_valid_body_has_nothing_to_glue() {
    let mut cube = cube_at(0.0, 0.0, 0.0);
    let before = NonManifoldCounts::of(&cube);
    assert_eq!(cube.glue_coincident(), Ok(CoincidentGlue::default()));
    assert_eq!(NonManifoldCounts::of(&cube), before);
}

/// A one-face planar sheet in the void over `points` (a convex polygon):
/// a wire cycle of lines spanned by a face.
fn polygon(points: &[(f64, f64, f64)]) -> BRep {
    use crate::geometry::curve::{CurveGeometry, CurveOnEdge};
    use crate::geometry::surface::{SurfaceGeometry, SurfaceOnFace};
    use crate::geometry::{Point3, Vec3};
    let mut body = BRep::new(tol());
    let outside = body.outside;
    let pts: Vec<Point3> = points
        .iter()
        .map(|&(x, y, z)| Point3::new(x, y, z))
        .collect();
    let v: Vec<_> = pts
        .iter()
        .map(|&p| body.make_acorn(outside, p).expect("acorn").vertex)
        .collect();
    let n = pts.len();
    let fins: Vec<_> = (0..n)
        .map(|i| {
            let made = body
                .make_wire_edge(outside, v[i], v[(i + 1) % n])
                .expect("wire");
            body.edges[made.edge].curve = Some(CurveOnEdge {
                geometry: CurveGeometry::Line {
                    origin: pts[i],
                    direction: pts[(i + 1) % n] - pts[i],
                },
                t_start: 0.0,
                t_end: 1.0,
            });
            made.fin
        })
        .collect();
    let u: Vec3 = pts[1] - pts[0];
    let normal = u.cross(&(pts[2] - pts[0])).normalized().expect("plane");
    let u = u.normalized().expect("axis");
    let plane = SurfaceOnFace {
        geometry: SurfaceGeometry::Plane {
            origin: pts[0],
            normal,
            u_axis: u,
            v_axis: normal.cross(&u),
        },
        orientation: true,
    };
    let made = body
        .make_face_from_wire_cycle(&fins, Some(plane))
        .expect("face");
    body.recompute_face_pcurves(made.face).expect("pcurves");
    body
}

fn region_of_sheet(body: &BRep, sheet_faces: usize) -> Option<crate::topology::RegionHandle> {
    let face = body.faces.values().nth(sheet_faces)?;
    body.shell_region(face.back)
}

/// A sheet glued onto a cube's edge on the material side is a lamina of
/// the cube: the radial order puts its fins in the solid wedge, and the
/// glue moves its shell into the solid region. Glued on the other side
/// it stays a sheet in the void. Both validate.
#[test]
fn a_sheet_glued_onto_a_solid_edge_takes_the_side_it_lies_in() {
    for (corner, inside) in [((0.5, 0.5), true), ((-0.5, -0.5), false)] {
        let mut body = cube_at(0.0, 0.0, 0.0);
        let (x, y) = corner;
        // A trapezoid in the plane x = y through the cube's edge on the
        // z axis, touching the cube nowhere else.
        let sheet = polygon(&[(0.0, 0.0, 0.0), (x, y, 0.3), (x, y, 0.7), (0.0, 0.0, 1.0)]);
        assert_eq!(sheet.faces.len(), 1, "a one-face sheet");
        body.append_body(&sheet);
        let glued = body.glue_coincident().expect("glue");
        assert_eq!((glued.vertices, glued.edges), (2, 1), "inside {inside}");
        let solid = body.solid_regions()[0];
        let sheet_region = region_of_sheet(&body, 6);
        assert_eq!(sheet_region == Some(solid), inside, "inside {inside}");
        assert_eq!(body.body_kind(), BodyKind::General);
        let errors = body.validate_all();
        assert!(errors.is_empty(), "inside {inside}: {errors:?}");
    }
}

/// Free shells lying in a cell join it when the cells are filled: a
/// lamina floating in a cube and an acorn at its centre become the
/// solid's; an acorn outside stays in the void.
#[test]
fn free_shells_join_the_cell_they_lie_in() {
    let mut body = cube_at(0.0, 0.0, 0.0);
    body.append_body(&polygon(&[
        (0.25, 0.25, 0.4),
        (0.75, 0.25, 0.4),
        (0.75, 0.75, 0.4),
        (0.25, 0.75, 0.4),
    ]));
    let outside = body.outside;
    let inner = body
        .make_acorn(outside, crate::geometry::Point3::new(0.5, 0.5, 0.6))
        .expect("acorn");
    let outer = body
        .make_acorn(outside, crate::geometry::Point3::new(2.0, 0.5, 0.5))
        .expect("acorn");
    assert!(body.fill_enclosed_cells().is_empty(), "no cells made");
    let solid = body.solid_regions()[0];
    assert_eq!(region_of_sheet(&body, 6), Some(solid));
    assert_eq!(body.shells[inner.shell].region, solid);
    assert_eq!(body.shells[outer.shell].region, outside);
    let errors = body.validate_all();
    assert!(errors.is_empty(), "{errors:?}");
}
