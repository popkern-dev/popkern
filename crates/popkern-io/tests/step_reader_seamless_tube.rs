//! PK-641: a closed band bounded by its rims alone — an extruded
//! ellipse tube written as `SURFACE_OF_LINEAR_EXTRUSION` with one
//! ellipse loop per end and no seam edge — gets the seam its chart
//! needs at import, meshes watertight, and encloses π·a·b·h.

use popkern_io::{StepImportOptions, import_step};

/// An elliptic cylinder: semi-axes `a`, `b`, along +Z from 0 to `h`.
/// `start` places each rim's vertex at that ellipse parameter (the
/// seam must be found off the rims' own vertices when it is not 0).
fn tube_step(a: f64, b: f64, h: f64, start: f64) -> String {
    let rim = |z: f64| {
        format!(
            "({:.17e},{:.17e},{z:.17e})",
            a * start.cos(),
            b * start.sin()
        )
    };
    format!(
        "ISO-10303-21;\nHEADER;\nFILE_DESCRIPTION(('tube'),'2;1');\nFILE_NAME('tube.step','',(''),(''),'','','');\nFILE_SCHEMA(('AUTOMOTIVE_DESIGN'));\nENDSEC;\nDATA;\n\
#1=CARTESIAN_POINT('',(0.,0.,0.));\n#2=DIRECTION('',(0.,0.,1.));\n#3=DIRECTION('',(1.,0.,0.));\n#4=AXIS2_PLACEMENT_3D('',#1,#2,#3);\n\
#5=ELLIPSE('',#4,{a:.17e},{b:.17e});\n#6=CARTESIAN_POINT('',(0.,0.,{h:.17e}));\n#7=AXIS2_PLACEMENT_3D('',#6,#2,#3);\n#8=ELLIPSE('',#7,{a:.17e},{b:.17e});\n\
#9=VECTOR('',#2,1.);\n#10=SURFACE_OF_LINEAR_EXTRUSION('',#5,#9);\n\
#11=CARTESIAN_POINT('',{r0});\n#12=VERTEX_POINT('',#11);\n#13=EDGE_CURVE('',#12,#12,#5,.T.);\n\
#14=CARTESIAN_POINT('',{r1});\n#15=VERTEX_POINT('',#14);\n#16=EDGE_CURVE('',#15,#15,#8,.T.);\n\
#20=ORIENTED_EDGE('',*,*,#13,.T.);\n#21=EDGE_LOOP('',(#20));\n#22=FACE_BOUND('',#21,.T.);\n\
#23=ORIENTED_EDGE('',*,*,#16,.F.);\n#24=EDGE_LOOP('',(#23));\n#25=FACE_BOUND('',#24,.T.);\n\
#26=ADVANCED_FACE('',(#22,#25),#10,.T.);\n\
#30=DIRECTION('',(0.,0.,-1.));\n#31=AXIS2_PLACEMENT_3D('',#1,#30,#3);\n#32=PLANE('',#31);\n#33=ORIENTED_EDGE('',*,*,#13,.F.);\n#34=EDGE_LOOP('',(#33));\n#35=FACE_OUTER_BOUND('',#34,.T.);\n#36=ADVANCED_FACE('',(#35),#32,.T.);\n\
#40=PLANE('',#7);\n#41=ORIENTED_EDGE('',*,*,#16,.T.);\n#42=EDGE_LOOP('',(#41));\n#43=FACE_OUTER_BOUND('',#42,.T.);\n#44=ADVANCED_FACE('',(#43),#40,.T.);\n\
#50=CLOSED_SHELL('',(#26,#36,#44));\n#51=MANIFOLD_SOLID_BREP('',#50);\n\
#52=(LENGTH_UNIT() NAMED_UNIT(*) SI_UNIT(.MILLI.,.METRE.));\n#53=(NAMED_UNIT(*) PLANE_ANGLE_UNIT() SI_UNIT($,.RADIAN.));\n\
#54=(GEOMETRIC_REPRESENTATION_CONTEXT(3) GLOBAL_UNIT_ASSIGNED_CONTEXT((#52,#53)) REPRESENTATION_CONTEXT('',''));\n\
#55=ADVANCED_BREP_SHAPE_REPRESENTATION('',(#51),#54);\nENDSEC;\nEND-ISO-10303-21;\n",
        r0 = rim(0.0),
        r1 = rim(h),
    )
}

fn check_tube(start: f64) {
    let (a, b, h) = (3.0, 1.5, 4.0);
    let text = tube_step(a, b, h, start);
    let imported = import_step(text.as_bytes(), &StepImportOptions::default()).expect("import");
    assert!(
        imported
            .diagnostics
            .entries()
            .iter()
            .any(|e| e.message.contains("added the missing seam")),
        "start {start}: the seam repair ran: {:?}",
        imported.diagnostics.entries()
    );
    let tube = imported
        .document
        .products
        .iter()
        .flat_map(|p| p.shapes.iter())
        .flat_map(|s| s.solids.iter().cloned())
        .next()
        .expect("solid");
    tube.validate().expect("valid");
    let mesh = popkern_mesh::tessellate(&tube, &popkern_mesh::TessellationOptions::default())
        .expect("mesh");
    assert_eq!(
        popkern_mesh::welded_open_edges(&mesh),
        0,
        "start {start}: watertight"
    );
    let got: f64 = mesh
        .triangles
        .iter()
        .map(|t| {
            let (p, q, r) = (
                mesh.vertices[t[0]].position,
                mesh.vertices[t[1]].position,
                mesh.vertices[t[2]].position,
            );
            (p.x * (q.y * r.z - q.z * r.y) - p.y * (q.x * r.z - q.z * r.x)
                + p.z * (q.x * r.y - q.y * r.x))
                / 6.0
        })
        .sum::<f64>()
        .abs();
    let exact = std::f64::consts::PI * a * b * h;
    assert!(
        (got - exact).abs() < 5e-3 * exact,
        "start {start}: volume {got} vs {exact}"
    );
}

#[test]
fn a_seamless_ellipse_tube_gets_its_seam_at_a_rim_vertex() {
    check_tube(0.0);
}

#[test]
fn a_seamless_ellipse_tube_whose_rim_vertices_are_off_the_seam() {
    check_tube(1.0);
}
