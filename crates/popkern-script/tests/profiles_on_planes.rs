//! PK-416: an exact profile extrudes to the same body whichever way it
//! is spelled (`path` or `polygon`) and whichever plane it sits on. The
//! sweep kernel's wire contract is the path frame's 2-D plane; the
//! recorder used to lift the wire onto the sketch plane first, which
//! collapsed every span along the plane's v axis — a rectangle on YZ
//! became two zero-normal faces (the `path` spelling tripped the
//! self-intersection gate, the `polygon` spelling served the
//! degenerate body).

use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::topology::BRep;
use popkern_script::prelude::*;
use popkern_script::{Model, Plane};

fn l_bracket_path(m: &mut Model, plane: Plane, label: &str) -> BRep {
    let body = m
        .profile(plane)
        .path(0.0, 0.0)
        .line_to(1.0, 0.0)
        .line_to(1.0, 0.15)
        .line_to(0.15, 0.15)
        .line_to(0.15, 0.8)
        .line_to(0.0, 0.8)
        .close()
        .sketch(label)
        .expect("sketch")
        .extrude(m, End::Blind(0.6.into()))
        .expect("extrude");
    m.brep_of(body).expect("brep")
}

fn l_bracket_polygon(m: &mut Model, plane: Plane, label: &str) -> BRep {
    let body = m
        .profile(plane)
        .polygon(&[
            (0.0, 0.0),
            (1.0, 0.0),
            (1.0, 0.15),
            (0.15, 0.15),
            (0.15, 0.8),
            (0.0, 0.8),
        ])
        .sketch(label)
        .expect("sketch")
        .extrude(m, End::Blind(0.6.into()))
        .expect("extrude");
    m.brep_of(body).expect("brep")
}

fn sorted_vertices(brep: &BRep) -> Vec<[i64; 3]> {
    let mut v: Vec<[i64; 3]> = brep
        .vertices
        .values()
        .map(|v| {
            [
                (v.position.x * 1e9).round() as i64,
                (v.position.y * 1e9).round() as i64,
                (v.position.z * 1e9).round() as i64,
            ]
        })
        .collect();
    v.sort_unstable();
    v
}

fn all_planes(brep: &BRep) -> bool {
    brep.faces.values().all(|f| {
        matches!(
            f.surface.as_ref().map(|s| &s.geometry),
            Some(SurfaceGeometry::Plane { .. })
        )
    })
}

#[test]
fn path_and_polygon_spellings_extrude_to_the_same_planar_body_on_every_plane() {
    let planes = [
        ("xy", Plane::XY),
        ("xz", Plane::XZ),
        ("yz", Plane::YZ),
        ("xy_at", Plane::xy_at(0.3)),
        ("xz_at", Plane::xz_at(-0.2)),
        ("yz_at", Plane::yz_at(0.5)),
    ];
    for (name, plane) in planes {
        let mut m = Model::new();
        let by_path = l_bracket_path(&mut m, plane.clone(), "l_path");
        let mut m2 = Model::new();
        let by_polygon = l_bracket_polygon(&mut m2, plane, "l_polygon");
        assert_eq!(by_path.faces.len(), 8, "[{name}] an L prism has 8 faces");
        assert_eq!(by_polygon.faces.len(), 8, "[{name}]");
        assert!(
            all_planes(&by_path),
            "[{name}] every face of a prism is a plane"
        );
        assert!(all_planes(&by_polygon), "[{name}]");
        assert_eq!(
            sorted_vertices(&by_path),
            sorted_vertices(&by_polygon),
            "[{name}] the two spellings place the same vertices"
        );
        by_path.validate().expect("valid");
        by_path
            .validate_curve_surface()
            .expect("curves on surfaces");
        // Every lateral has a live normal (no collapsed span).
        for face in by_path.faces.values() {
            let Some(SurfaceGeometry::Plane { normal, .. }) =
                face.surface.as_ref().map(|s| &s.geometry)
            else {
                panic!("[{name}] plane");
            };
            assert!((normal.length() - 1.0).abs() < 1e-9, "[{name}] unit normal");
        }
        m.assert_watertight().expect("watertight");
    }
}

#[test]
fn a_rectangle_on_yz_is_a_box_with_its_vertices_on_the_plane() {
    let mut m = Model::new();
    let body = m
        .profile(Plane::YZ)
        .path(0.0, 0.0)
        .line_to(1.0, 0.0)
        .line_to(1.0, 0.5)
        .line_to(0.0, 0.5)
        .close()
        .sketch("r")
        .expect("sketch")
        .extrude(&mut m, End::Blind(0.6.into()))
        .expect("extrude");
    let brep = m.brep_of(body).expect("brep");
    assert_eq!(brep.faces.len(), 6);
    assert_eq!(brep.vertices.len(), 8);
    let xs: Vec<f64> = brep.vertices.values().map(|v| v.position.x).collect();
    assert!(
        xs.iter()
            .all(|x| x.abs() < 1e-12 || (x - 0.6).abs() < 1e-12),
        "extruded along +x from the YZ plane: {xs:?}"
    );
}
