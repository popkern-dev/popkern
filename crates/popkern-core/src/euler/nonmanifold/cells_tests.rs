use super::*;
use crate::geometry::transform::Transform3;
use crate::tolerance::Tolerance;
use crate::topology::{BodyKind, FaceHandle};

fn tol() -> Tolerance {
    Tolerance::default()
}

fn cube_at(x: f64, y: f64, z: f64) -> BRep {
    BRep::unit_cube(tol())
        .and_then(|c| c.transform(&Transform3::translation(x, y, z)))
        .expect("cube")
}

/// A cube's region killed leaves its surface a closed sheet: refilling
/// the enclosed cells gives the solid back.
#[test]
fn a_closed_sheet_encloses_one_cell() {
    let mut cube = cube_at(0.0, 0.0, 0.0);
    let solid = cube.solid_regions()[0];
    cube.kill_face_region(solid).expect("kill");
    let made = cube.fill_enclosed_cells();
    assert_eq!(made.len(), 1);
    assert_eq!(cube.body_kind(), BodyKind::Solid);
    let errors = cube.validate_all();
    assert!(errors.is_empty(), "{errors:?}");
}

/// A box standing on a box with one internal face, its two solid
/// regions killed: the faces bound two cells, the internal face both,
/// and refilling gives two solids again.
#[test]
fn two_cells_share_an_internal_face() {
    let mut body = cube_at(0.0, 0.0, 0.0);
    body.append_body(&cube_at(0.0, 0.0, 1.0));
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
    let at_one: Vec<FaceHandle> = body
        .faces
        .keys()
        .filter(|&f| {
            body.loop_cycle(body.faces[f].outer_loop)
                .iter()
                .all(|&h| (body.vertices[body.half_edges[h].origin].position.z - 1.0).abs() < 1e-9)
        })
        .collect();
    body.glue_faces(at_one[0], at_one[1]).expect("faces");
    // Dissolve both solids into the outside: every side in the void.
    for region in body.solid_regions() {
        for shell in body.region_shells(region) {
            body.shells[shell].region = body.outside;
        }
        body.regions.remove(region);
    }
    let made = body.fill_enclosed_cells();
    assert_eq!(made.len(), 2);
    assert_eq!(body.solid_regions().len(), 2);
    let errors = body.validate_all();
    assert!(errors.is_empty(), "{errors:?}");
}

/// A hollow cube with its regions dissolved: the outer faces' material
/// side encloses positive volume (a solid), the cavity's material side
/// negative (an inner shell of that solid, found by winding), the
/// cavity's own side positive with no back side (a bounded void).
#[test]
fn a_cavity_nests_inside_its_solid() {
    let mut host = cube_at(0.0, 0.0, 0.0);
    let mut tool = BRep::unit_cube(tol())
        .and_then(|c| c.transform(&Transform3::scale_uniform(0.5)))
        .and_then(|c| c.transform(&Transform3::translation(0.25, 0.25, 0.25)))
        .expect("tool");
    tool.reverse_body_orientation();
    let solid = host.solid_regions()[0];
    let map = host.append_body(&tool);
    host.invert_regions_into(&[map.regions[&tool.solid_regions()[0]]], solid);
    let outside = host.outside;
    let regions: Vec<_> = host.regions.keys().filter(|&r| r != outside).collect();
    for region in regions {
        for shell in host.region_shells(region) {
            host.shells[shell].region = outside;
        }
        host.regions.remove(region);
    }
    let made = host.fill_enclosed_cells();
    assert_eq!(made.len(), 2, "a solid and the cavity's void");
    let solids = host.solid_regions();
    assert_eq!(solids.len(), 1);
    assert_eq!(
        host.region_shells(solids[0]).len(),
        2,
        "outer and inner shells"
    );
    assert_eq!(host.cavity_regions().len(), 1);
    let errors = host.validate_all();
    assert!(errors.is_empty(), "{errors:?}");
}
