//! PK-644: a cone written as a rational B-spline whose first control
//! row is its apex (one point), split into two half sectors bounded by
//! rulings through the apex, meshes watertight on the boundary-first
//! path. The rulings reach the apex in different chart columns; the
//! ring runs along the collapsed row between them (measured: f1-db's
//! cone sectors fell back to the uniform grid, 562 open mesh edges).
//! Each sector is half a turn, so the row run is the half-period tie
//! the chain lift must resolve by closing the ring.

use popkern_io::{StepImportOptions, import_step};

fn p(x: f64, y: f64, z: f64) -> String {
    format!("({x:.17e},{y:.17e},{z:.17e})")
}

/// Cone of base radius `r` (z = 0) and apex height `h`, its side one
/// closed B-spline surface (u: apex row 0 to rim row 1, v: the angle).
fn cone_step(r: f64, h: f64, (rot_a, rot_b): (usize, usize)) -> String {
    let w = std::f64::consts::FRAC_1_SQRT_2;
    let square = [
        (1.0, 0.0),
        (1.0, 1.0),
        (0.0, 1.0),
        (-1.0, 1.0),
        (-1.0, 0.0),
        (-1.0, -1.0),
        (0.0, -1.0),
        (1.0, -1.0),
        (1.0, 0.0),
    ];
    let mut points = String::new();
    let mut rows = Vec::new();
    for row in 0..2 {
        let mut ids = Vec::new();
        for (k, (x, y)) in square.iter().enumerate() {
            let id = 200 + row * 10 + k;
            let pt = if row == 0 {
                p(0.0, 0.0, h)
            } else {
                p(r * x, r * y, 0.0)
            };
            points.push_str(&format!("#{id}=CARTESIAN_POINT('',{pt});\n"));
            ids.push(format!("#{id}"));
        }
        rows.push(format!("({})", ids.join(",")));
    }
    let weights = (0..9)
        .map(|k| {
            if k % 2 == 0 {
                "1.".to_string()
            } else {
                format!("{w:.17e}")
            }
        })
        .collect::<Vec<_>>()
        .join(",");
    let len = (r * r + h * h).sqrt();
    // Each sector's loop starts at a chosen edge: where the walk meets
    // the apex decides which way the half-period row run first lifts.
    let rotate = |edges: [&str; 3], k: usize| -> String {
        (0..3)
            .map(|i| edges[(i + k) % 3])
            .collect::<Vec<_>>()
            .join(",")
    };
    let loop_a = rotate(["#40", "#41", "#42"], rot_a);
    let loop_b = rotate(["#50", "#51", "#52"], rot_b);
    format!(
        "ISO-10303-21;\nHEADER;\nFILE_DESCRIPTION(('cone'),'2;1');\nFILE_NAME('cone.step','',(''),(''),'','','');\nFILE_SCHEMA(('AUTOMOTIVE_DESIGN'));\nENDSEC;\nDATA;\n\
#1=CARTESIAN_POINT('',(0.,0.,0.));\n#2=DIRECTION('',(0.,0.,1.));\n#3=DIRECTION('',(1.,0.,0.));\n#4=AXIS2_PLACEMENT_3D('',#1,#2,#3);\n#5=CIRCLE('',#4,{r:.17e});\n\
#6=CARTESIAN_POINT('',{apex});\n#7=VERTEX_POINT('',#6);\n#8=CARTESIAN_POINT('',{right});\n#9=VERTEX_POINT('',#8);\n#10=CARTESIAN_POINT('',{left});\n#11=VERTEX_POINT('',#10);\n\
#12=EDGE_CURVE('',#9,#11,#5,.T.);\n#13=EDGE_CURVE('',#11,#9,#5,.T.);\n\
#14=DIRECTION('',{d0});\n#15=VECTOR('',#14,{len:.17e});\n#16=LINE('',#6,#15);\n#17=EDGE_CURVE('',#7,#9,#16,.T.);\n\
#18=DIRECTION('',{d1});\n#19=VECTOR('',#18,{len:.17e});\n#20=LINE('',#6,#19);\n#21=EDGE_CURVE('',#7,#11,#20,.T.);\n\
{points}\
#30=(BOUNDED_SURFACE() B_SPLINE_SURFACE(1,2,({rows}),.UNSPECIFIED.,.F.,.T.,.F.) \
B_SPLINE_SURFACE_WITH_KNOTS((2,2),(3,2,2,2,3),(0.,1.),(0.,1.5707963267948966,3.141592653589793,4.71238898038469,6.283185307179586),.UNSPECIFIED.) \
GEOMETRIC_REPRESENTATION_ITEM() RATIONAL_B_SPLINE_SURFACE((({weights}),({weights}))) REPRESENTATION_ITEM('') SURFACE());\n\
#40=ORIENTED_EDGE('',*,*,#17,.T.);\n#41=ORIENTED_EDGE('',*,*,#12,.T.);\n#42=ORIENTED_EDGE('',*,*,#21,.F.);\n#43=EDGE_LOOP('',({loop_a}));\n#44=FACE_OUTER_BOUND('',#43,.T.);\n#45=ADVANCED_FACE('',(#44),#30,.T.);\n\
#50=ORIENTED_EDGE('',*,*,#21,.T.);\n#51=ORIENTED_EDGE('',*,*,#13,.T.);\n#52=ORIENTED_EDGE('',*,*,#17,.F.);\n#53=EDGE_LOOP('',({loop_b}));\n#54=FACE_OUTER_BOUND('',#53,.T.);\n#55=ADVANCED_FACE('',(#54),#30,.T.);\n\
#60=DIRECTION('',(0.,0.,-1.));\n#61=AXIS2_PLACEMENT_3D('',#1,#60,#3);\n#62=PLANE('',#61);\n\
#63=ORIENTED_EDGE('',*,*,#12,.F.);\n#64=ORIENTED_EDGE('',*,*,#13,.F.);\n#65=EDGE_LOOP('',(#63,#64));\n#66=FACE_OUTER_BOUND('',#65,.T.);\n#67=ADVANCED_FACE('',(#66),#62,.T.);\n\
#70=CLOSED_SHELL('',(#45,#55,#67));\n#71=MANIFOLD_SOLID_BREP('',#70);\n\
#72=(LENGTH_UNIT() NAMED_UNIT(*) SI_UNIT(.MILLI.,.METRE.));\n#73=(NAMED_UNIT(*) PLANE_ANGLE_UNIT() SI_UNIT($,.RADIAN.));\n\
#74=(GEOMETRIC_REPRESENTATION_CONTEXT(3) GLOBAL_UNIT_ASSIGNED_CONTEXT((#72,#73)) REPRESENTATION_CONTEXT('',''));\n\
#75=ADVANCED_BREP_SHAPE_REPRESENTATION('',(#71),#74);\nENDSEC;\nEND-ISO-10303-21;\n",
        apex = p(0.0, 0.0, h),
        right = p(r, 0.0, 0.0),
        left = p(-r, 0.0, 0.0),
        d0 = p(r / len, 0.0, -h / len),
        d1 = p(-r / len, 0.0, -h / len),
        rows = rows.join(","),
    )
}

#[test]
fn a_bspline_cone_in_two_half_sectors_meshes_through_its_apex_row() {
    for rotations in (0..3).flat_map(|a| (0..3).map(move |b| (a, b))) {
        check_cone(rotations);
    }
}

fn check_cone(rotations: (usize, usize)) {
    let (r, h) = (30.0, 40.0);
    let imported = import_step(
        cone_step(r, h, rotations).as_bytes(),
        &StepImportOptions::default(),
    )
    .expect("import");
    let cone = imported
        .document
        .products
        .iter()
        .flat_map(|p| p.shapes.iter())
        .flat_map(|s| s.solids.iter().cloned())
        .next()
        .expect("solid");
    cone.validate().expect("valid");
    let mesh = popkern_mesh::tessellate(&cone, &popkern_mesh::TessellationOptions::default())
        .expect("mesh");
    assert_eq!(
        popkern_mesh::welded_open_edges(&mesh),
        0,
        "{rotations:?}: watertight"
    );
    let volume: f64 = mesh
        .triangles
        .iter()
        .map(|t| {
            let (a, b, c) = (
                mesh.vertices[t[0]].position,
                mesh.vertices[t[1]].position,
                mesh.vertices[t[2]].position,
            );
            (a.x * (b.y * c.z - b.z * c.y) - a.y * (b.x * c.z - b.z * c.x)
                + a.z * (b.x * c.y - b.y * c.x))
                / 6.0
        })
        .sum();
    let exact = std::f64::consts::PI * r * r * h / 3.0;
    assert!(
        (volume - exact).abs() < 2e-3 * exact,
        "{rotations:?}: signed volume {volume} vs {exact}"
    );
}
