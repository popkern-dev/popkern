//! PK-646: a closed NURBS band whose chart closes in V — a cylinder
//! written as a rational B-spline revolve (u along the axis, v the
//! angle), bounded by its two rim circles with no seam edge — gets a
//! seam along its exact iso-curve at import, its rims' pcurves are
//! re-projected honestly (the closed rims' pre-seam pcurves had
//! collapsed to points), and it meshes watertight enclosing π·r²·h
//! (measured: jack_in_the_box body 27's revolved collar, open by 1024
//! mesh edges with a 1.26 pcurve deviation).

use popkern_io::{StepImportOptions, import_step};

fn p(x: f64, y: f64, z: f64) -> String {
    format!("({x:.17e},{y:.17e},{z:.17e})")
}

/// A cylinder of radius `r` along +Z from 0 to `h`, its side a
/// rational B-spline closed in v; each rim's vertex at angle `start`.
fn band_step(r: f64, h: f64, start: f64) -> String {
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
    for (row, z) in [0.0, h].into_iter().enumerate() {
        let mut ids = Vec::new();
        for (k, (x, y)) in square.iter().enumerate() {
            let id = 100 + row * 10 + k;
            points.push_str(&format!(
                "#{id}=CARTESIAN_POINT('',{});\n",
                p(r * x, r * y, z)
            ));
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
    let rim = |z: f64| p(r * start.cos(), r * start.sin(), z);
    format!(
        "ISO-10303-21;\nHEADER;\nFILE_DESCRIPTION(('band'),'2;1');\nFILE_NAME('band.step','',(''),(''),'','','');\nFILE_SCHEMA(('AUTOMOTIVE_DESIGN'));\nENDSEC;\nDATA;\n\
#1=CARTESIAN_POINT('',(0.,0.,0.));\n#2=DIRECTION('',(0.,0.,1.));\n#3=DIRECTION('',(1.,0.,0.));\n#4=AXIS2_PLACEMENT_3D('',#1,#2,#3);\n\
#5=CIRCLE('',#4,{r:.17e});\n#6=CARTESIAN_POINT('',(0.,0.,{h:.17e}));\n#7=AXIS2_PLACEMENT_3D('',#6,#2,#3);\n#8=CIRCLE('',#7,{r:.17e});\n\
{points}\
#10=(BOUNDED_SURFACE() B_SPLINE_SURFACE(1,2,({rows}),.UNSPECIFIED.,.F.,.T.,.F.) \
B_SPLINE_SURFACE_WITH_KNOTS((2,2),(3,2,2,2,3),(0.,{h:.17e}),(0.,1.5707963267948966,3.141592653589793,4.71238898038469,6.283185307179586),.UNSPECIFIED.) \
GEOMETRIC_REPRESENTATION_ITEM() RATIONAL_B_SPLINE_SURFACE((({weights}),({weights}))) REPRESENTATION_ITEM('') SURFACE());\n\
#11=CARTESIAN_POINT('',{r0});\n#12=VERTEX_POINT('',#11);\n#13=EDGE_CURVE('',#12,#12,#5,.T.);\n\
#14=CARTESIAN_POINT('',{r1});\n#15=VERTEX_POINT('',#14);\n#16=EDGE_CURVE('',#15,#15,#8,.T.);\n\
#20=ORIENTED_EDGE('',*,*,#13,.T.);\n#21=EDGE_LOOP('',(#20));\n#22=FACE_BOUND('',#21,.T.);\n\
#23=ORIENTED_EDGE('',*,*,#16,.F.);\n#24=EDGE_LOOP('',(#23));\n#25=FACE_BOUND('',#24,.T.);\n\
#26=ADVANCED_FACE('',(#22,#25),#10,.F.);\n\
#30=DIRECTION('',(0.,0.,-1.));\n#31=AXIS2_PLACEMENT_3D('',#1,#30,#3);\n#32=PLANE('',#31);\n#33=ORIENTED_EDGE('',*,*,#13,.F.);\n#34=EDGE_LOOP('',(#33));\n#35=FACE_OUTER_BOUND('',#34,.T.);\n#36=ADVANCED_FACE('',(#35),#32,.T.);\n\
#40=PLANE('',#7);\n#41=ORIENTED_EDGE('',*,*,#16,.T.);\n#42=EDGE_LOOP('',(#41));\n#43=FACE_OUTER_BOUND('',#42,.T.);\n#44=ADVANCED_FACE('',(#43),#40,.T.);\n\
#50=CLOSED_SHELL('',(#26,#36,#44));\n#51=MANIFOLD_SOLID_BREP('',#50);\n\
#52=(LENGTH_UNIT() NAMED_UNIT(*) SI_UNIT(.MILLI.,.METRE.));\n#53=(NAMED_UNIT(*) PLANE_ANGLE_UNIT() SI_UNIT($,.RADIAN.));\n\
#54=(GEOMETRIC_REPRESENTATION_CONTEXT(3) GLOBAL_UNIT_ASSIGNED_CONTEXT((#52,#53)) REPRESENTATION_CONTEXT('',''));\n\
#55=ADVANCED_BREP_SHAPE_REPRESENTATION('',(#51),#54);\nENDSEC;\nEND-ISO-10303-21;\n",
        rows = rows.join(","),
        r0 = rim(0.0),
        r1 = rim(h),
    )
}

fn check_band(start: f64) {
    let (r, h) = (2.0, 3.0);
    let text = band_step(r, h, start);
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
    let band = imported
        .document
        .products
        .iter()
        .flat_map(|p| p.shapes.iter())
        .flat_map(|s| s.solids.iter().cloned())
        .next()
        .expect("solid");
    band.validate().expect("valid");
    // Every pcurve is an honest projection of its edge: no recorded
    // deviation beyond the kernel's own precision.
    for (eh, _) in band.edges.iter() {
        let dev = band.measure_edge_deviation(eh);
        assert!(
            dev < 1e-5,
            "start {start}: edge {eh:?} pcurve deviation {dev}"
        );
    }
    let mesh = popkern_mesh::tessellate(&band, &popkern_mesh::TessellationOptions::default())
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
            let (a, b, c) = (
                mesh.vertices[t[0]].position,
                mesh.vertices[t[1]].position,
                mesh.vertices[t[2]].position,
            );
            (a.x * (b.y * c.z - b.z * c.y) - a.y * (b.x * c.z - b.z * c.x)
                + a.z * (b.x * c.y - b.y * c.x))
                / 6.0
        })
        .sum::<f64>()
        .abs();
    let exact = std::f64::consts::PI * r * r * h;
    assert!(
        (got - exact).abs() < 5e-3 * exact,
        "start {start}: volume {got} vs {exact}"
    );
}

#[test]
fn a_band_closed_in_v_gets_its_seam_at_the_rim_vertices() {
    check_band(0.0);
}

#[test]
fn a_band_closed_in_v_whose_rim_vertices_are_off_the_seam() {
    check_band(2.0);
}
