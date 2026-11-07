use super::*;
use popkern_core::geometry::Point3;
use popkern_core::geometry::Vec3;
use popkern_core::geometry::surface::{SurfaceGeometry, SurfaceOnFace};
use popkern_core::geometry::transform::Transform3;

fn tol() -> Tolerance {
    Tolerance::default()
}

fn cube_at(x: f64, y: f64, z: f64) -> BRep {
    BRep::unit_cube(tol())
        .and_then(|c| c.transform(&Transform3::translation(x, y, z)))
        .expect("cube")
}

/// Glue every pair of coincident vertices.
fn glue_all_vertices(body: &mut BRep) {
    loop {
        let keys: Vec<_> = body.vertices.keys().collect();
        let pair = keys.iter().enumerate().find_map(|(i, &a)| {
            keys[i + 1..]
                .iter()
                .copied()
                .find(|&b| {
                    body.vertices[a]
                        .position
                        .distance_to(&body.vertices[b].position)
                        < 1e-9
                })
                .map(|b| (a, b))
        });
        let Some((a, b)) = pair else { break };
        body.glue_vertices(a, b).expect("vertex");
    }
}

fn volumes(body: &BRep) -> Vec<f64> {
    let mut v: Vec<f64> = region_volumes(body, &tol())
        .expect("volumes")
        .into_iter()
        .map(|(_, v)| v)
        .collect();
    v.sort_by(f64::total_cmp);
    v
}

fn close(a: &[f64], b: &[f64]) {
    assert_eq!(a.len(), b.len(), "{a:?} vs {b:?}");
    for (x, y) in a.iter().zip(b) {
        assert!((x - y).abs() < 1e-9, "{a:?} vs {b:?}");
    }
}

#[test]
fn a_cube_has_one_region_of_its_volume() {
    close(&volumes(&cube_at(0.0, 0.0, 0.0)), &[1.0]);
}

/// Two cubes glued along an edge: two regions of one each, and the
/// body's volume their sum.
#[test]
fn cubes_glued_along_an_edge_have_a_region_each() {
    let mut body = cube_at(0.0, 0.0, 0.0);
    body.append_body(&cube_at(1.0, 1.0, 0.0));
    glue_all_vertices(&mut body);
    let edges: Vec<_> = body
        .edges
        .keys()
        .filter(|&e| {
            body.edge_vertices(e).is_some_and(|(a, b)| {
                let (p, q) = (body.vertices[a].position, body.vertices[b].position);
                (p.x - 1.0).abs() < 1e-9
                    && (p.y - 1.0).abs() < 1e-9
                    && (q.x - 1.0).abs() < 1e-9
                    && (q.y - 1.0).abs() < 1e-9
            })
        })
        .collect();
    body.glue_edges(edges[0], edges[1]).expect("edge");
    assert_eq!(body.body_kind(), BodyKind::General);
    close(&volumes(&body), &[1.0, 1.0]);
    let total = crate::body_volume(&body, &tol()).expect("volume");
    assert!((total - 2.0).abs() < 1e-9, "{total}");
}

/// A box standing on a box with one internal face between them: the
/// face counts once for each solid.
#[test]
fn an_internal_face_bounds_both_solids() {
    let mut body = cube_at(0.0, 0.0, 0.0);
    body.append_body(&cube_at(0.0, 0.0, 1.0));
    glue_all_vertices(&mut body);
    let at_one: Vec<_> = body
        .faces
        .keys()
        .filter(|&f| {
            let outer = body.faces[f].outer_loop;
            let start = body.loops[outer].half_edge;
            let mut h = start;
            let mut all = true;
            loop {
                all &= (body.vertices[body.half_edges[h].origin].position.z - 1.0).abs() < 1e-9;
                h = body.half_edges[h].next;
                if h == start {
                    break;
                }
            }
            all
        })
        .collect();
    body.glue_faces(at_one[0], at_one[1]).expect("faces");
    close(&volumes(&body), &[1.0, 1.0]);
}

/// A lamina inside a cube bounds no volume of its own: both its sides
/// face the cube's region and its terms cancel.
#[test]
fn a_lamina_inside_a_solid_adds_no_volume() {
    let mut cube = cube_at(0.0, 0.0, 0.0);
    let solid = cube.solid_regions()[0];
    let corners = [(0.25, 0.25), (0.75, 0.25), (0.75, 0.75), (0.25, 0.75)];
    let v: Vec<_> = corners
        .iter()
        .map(|&(x, y)| {
            cube.make_acorn(solid, Point3::new(x, y, 0.5))
                .expect("acorn")
                .vertex
        })
        .collect();
    let cycle: Vec<_> = (0..4)
        .map(|i| {
            cube.make_wire_edge(solid, v[i], v[(i + 1) % 4])
                .expect("wire")
                .fin
        })
        .collect();
    let face = cube
        .make_face_from_wire_cycle(
            &cycle,
            Some(SurfaceOnFace {
                geometry: SurfaceGeometry::Plane {
                    origin: Point3::new(0.0, 0.0, 0.5),
                    normal: Vec3::new(0.0, 0.0, 1.0),
                    u_axis: Vec3::new(1.0, 0.0, 0.0),
                    v_axis: Vec3::new(0.0, 1.0, 0.0),
                },
                orientation: true,
            }),
        )
        .expect("lamina")
        .face;
    cube.recompute_face_pcurves(face).expect("pcurves");
    assert_eq!(cube.body_kind(), BodyKind::General);
    close(&volumes(&cube), &[1.0]);
}
