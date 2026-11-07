use super::*;
use popkern_core::geometry::transform::Transform3;
use popkern_core::tolerance::Tolerance;

fn tol() -> Tolerance {
    Tolerance::default()
}

fn cube_at(x: f64, y: f64, z: f64) -> BRep {
    BRep::unit_cube(tol())
        .and_then(|c| c.transform(&Transform3::translation(x, y, z)))
        .expect("cube")
}

fn glue_coincident_vertices(body: &mut BRep) {
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

/// The volume `region`'s triangles enclose, each turned to face out.
fn region_volume(m: &GeneralMesh, region: RegionHandle) -> f64 {
    m.region_triangles(region)
        .into_iter()
        .map(|(i, reversed)| {
            let [a, b, c] = m.mesh.triangles[i].map(|k| m.mesh.vertices[k].position);
            let (b, c) = if reversed { (c, b) } else { (b, c) };
            let v = |p: Point3| popkern_core::geometry::Vec3::new(p.x, p.y, p.z);
            v(a).dot(&v(b).cross(&v(c))) / 6.0
        })
        .sum()
}

/// Two cubes glued along an edge: every face bounds one solid, and each
/// solid's triangles, faced out of it, enclose its volume.
#[test]
fn glued_cubes_mesh_with_a_closed_boundary_per_solid() {
    let mut body = cube_at(0.0, 0.0, 0.0);
    body.append_body(&cube_at(1.0, 1.0, 0.0));
    glue_coincident_vertices(&mut body);
    let shared: Vec<_> = body
        .edges
        .keys()
        .filter(|&e| {
            body.edge_vertices(e).is_some_and(|(a, b)| {
                [a, b].iter().all(|&v| {
                    let p = body.vertices[v].position;
                    (p.x - 1.0).abs() < 1e-9 && (p.y - 1.0).abs() < 1e-9
                })
            })
        })
        .collect();
    body.glue_edges(shared[0], shared[1]).expect("edge");
    let m = tessellate_general(&body, &TessellationOptions::default()).expect("mesh");
    assert_eq!(m.faces.len(), 12);
    assert!(m.faces.iter().all(|f| !f.two_sided && !f.internal));
    for region in body.solid_regions() {
        let v = region_volume(&m, region);
        assert!((v - 1.0).abs() < 1e-9, "{v}");
    }
}

/// A box standing on a box: the face between them is internal, and it
/// closes both solids' boundaries.
#[test]
fn an_internal_face_is_tagged_and_bounds_both_solids() {
    let mut body = cube_at(0.0, 0.0, 0.0);
    body.append_body(&cube_at(0.0, 0.0, 1.0));
    glue_coincident_vertices(&mut body);
    let at_one: Vec<FaceHandle> = body
        .faces
        .keys()
        .filter(|&f| {
            body.face_outer_loop_positions(f)
                .iter()
                .all(|p| (p.z - 1.0).abs() < 1e-9)
        })
        .collect();
    body.glue_faces(at_one[0], at_one[1]).expect("faces");
    let m = tessellate_general(&body, &TessellationOptions::default()).expect("mesh");
    assert_eq!(m.faces.iter().filter(|f| f.internal).count(), 1);
    for region in body.solid_regions() {
        let v = region_volume(&m, region);
        assert!((v - 1.0).abs() < 1e-9, "{v}");
    }
}

/// A wire between two points and an acorn mesh as a polyline and a
/// point; a face-less body has no triangles.
#[test]
fn wires_and_acorns_mesh_as_polylines_and_points() {
    let mut body = BRep::new(tol());
    let out = body.outside;
    let a = body
        .make_acorn(out, Point3::new(0.0, 0.0, 0.0))
        .expect("a")
        .vertex;
    let b = body
        .make_acorn(out, Point3::new(1.0, 0.0, 0.0))
        .expect("b")
        .vertex;
    body.make_wire_edge(out, a, b).expect("wire");
    let lone = body
        .make_acorn(out, Point3::new(5.0, 5.0, 5.0))
        .expect("lone")
        .vertex;
    let m = tessellate_general(&body, &TessellationOptions::default()).expect("mesh");
    assert!(m.mesh.triangles.is_empty());
    assert_eq!(m.wires.len(), 1);
    assert_eq!(
        m.wires[0].points,
        vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)]
    );
    assert_eq!(m.acorns, vec![(lone, Point3::new(5.0, 5.0, 5.0))]);
}

/// A sheet's faces are seen from both sides.
#[test]
fn a_sheet_is_two_sided() {
    let mut open = cube_at(0.0, 0.0, 0.0);
    let face = open.faces.keys().next().expect("face");
    open.retire_faces_to_boundary(&[face]).expect("retire");
    open.settle_open_shells();
    let m = tessellate_general(&open, &TessellationOptions::default()).expect("mesh");
    assert_eq!(m.faces.len(), 5);
    assert!(m.faces.iter().all(|f| f.two_sided && !f.internal));
}
