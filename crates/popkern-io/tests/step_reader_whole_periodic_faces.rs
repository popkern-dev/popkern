//! PK-654: whole periodic faces mesh to the chord height. A STEP torus
//! bounded only by a `VERTEX_LOOP` (moon_buggy writes its tyres this
//! way, and its spheres) imports as a face with no boundary edges, and
//! meshes on the uniform-grid path, whose quadtree was capped at 64 cells a side and
//! tested each cell's bilinear interior only: the r 254 sphere's facets
//! sat 0.34 off and the torus's saddle (inner equator) 0.54 off at a
//! 0.005 chord. The bound here is twice the chord height, the facet
//! gate the boundary-first census uses.

use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::topology::BRep;
use popkern_io::{StepImportOptions, import_step};

/// Worst facet-centroid distance to its face's surface, and the
/// welded open-edge count.
fn facet_census(brep: &BRep) -> (f64, usize) {
    let mesh = popkern_mesh::tessellate(brep, &popkern_mesh::TessellationOptions::default())
        .expect("mesh");
    let mut worst: f64 = 0.0;
    for (ti, t) in mesh.triangles.iter().enumerate() {
        let surface = &brep.faces[mesh.face_map[ti]]
            .surface
            .as_ref()
            .expect("surface")
            .geometry;
        let (a, b, c) = (
            mesh.vertices[t[0]].position,
            mesh.vertices[t[1]].position,
            mesh.vertices[t[2]].position,
        );
        let g = popkern_core::geometry::Point3::new(
            (a.x + b.x + c.x) / 3.0,
            (a.y + b.y + c.y) / 3.0,
            (a.z + b.z + c.z) / 3.0,
        );
        let (u, v) = surface.closest_point(&g).expect("projection");
        worst = worst.max(surface.point_at(u, v).distance_to(&g));
    }
    (worst, popkern_mesh::welded_open_edges(&mesh))
}

/// One face on `surface` (a STEP entity line numbered #5) bounded only
/// by a vertex loop at `vertex`.
fn vertex_loop_face_step(surface: &str, vertex: (f64, f64, f64)) -> String {
    format!(
        "ISO-10303-21;\nHEADER;\nFILE_DESCRIPTION(('torus'),'2;1');\nFILE_NAME('torus.step','',(''),(''),'','','');\nFILE_SCHEMA(('AUTOMOTIVE_DESIGN'));\nENDSEC;\nDATA;\n\
#1=CARTESIAN_POINT('',(0.,0.,0.));\n#2=DIRECTION('',(0.,0.,1.));\n#3=DIRECTION('',(1.,0.,0.));\n#4=AXIS2_PLACEMENT_3D('',#1,#2,#3);\n\
#5={surface};\n\
#6=CARTESIAN_POINT('',({x:.17e},{y:.17e},{z:.17e}));\n#7=VERTEX_POINT('',#6);\n#8=VERTEX_LOOP('',#7);\n#9=FACE_OUTER_BOUND('',#8,.T.);\n\
#10=ADVANCED_FACE('',(#9),#5,.T.);\n#11=CLOSED_SHELL('',(#10));\n#12=MANIFOLD_SOLID_BREP('',#11);\n\
#13=(LENGTH_UNIT() NAMED_UNIT(*) SI_UNIT(.MILLI.,.METRE.));\n#14=(NAMED_UNIT(*) PLANE_ANGLE_UNIT() SI_UNIT($,.RADIAN.));\n\
#15=(GEOMETRIC_REPRESENTATION_CONTEXT(3) GLOBAL_UNIT_ASSIGNED_CONTEXT((#13,#14)) REPRESENTATION_CONTEXT('',''));\n\
#16=ADVANCED_BREP_SHAPE_REPRESENTATION('',(#12),#15);\nENDSEC;\nEND-ISO-10303-21;\n",
        x = vertex.0,
        y = vertex.1,
        z = vertex.2,
    )
}

fn import_one(text: &str) -> BRep {
    let imported = import_step(text.as_bytes(), &StepImportOptions::default()).expect("import");
    imported
        .document
        .products
        .iter()
        .flat_map(|p| p.shapes.iter())
        .flat_map(|s| s.solids.iter().cloned())
        .next()
        .expect("solid")
}

#[test]
fn a_vertex_loop_torus_meshes_its_saddle_within_twice_the_chord() {
    let torus = import_one(&vertex_loop_face_step(
        "TOROIDAL_SURFACE('',#4,254.,127.)",
        (381.0, 0.0, 0.0),
    ));
    let chord = popkern_mesh::TessellationOptions::default().chord_height;
    let (worst, open) = facet_census(&torus);
    assert_eq!(open, 0, "watertight");
    assert!(worst <= 2.0 * chord, "torus facet {worst} off the surface");
}

#[test]
fn a_vertex_loop_sphere_meshes_within_twice_the_chord() {
    let sphere = import_one(&vertex_loop_face_step(
        "SPHERICAL_SURFACE('',#4,254.)",
        (0.0, 0.0, 254.0),
    ));
    let chord = popkern_mesh::TessellationOptions::default().chord_height;
    let (worst, open) = facet_census(&sphere);
    assert_eq!(open, 0, "watertight");
    assert!(worst <= 2.0 * chord, "sphere facet {worst} off the surface");
}
