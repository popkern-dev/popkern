//! PK-636: a domed solid as third-party writers spell it — a sphere
//! face bounded by ONE circle edge (no seam, the pole inside the face)
//! plus the cutting disk — imports valid, meshes watertight, and
//! encloses the analytic spherical-cap volume πh²(3R − h)/3. The cap
//! is cut by planes on both sides of the centre and at a tilt, so it
//! contains the sphere chart's north pole, its south pole, or neither.

use popkern_io::{StepImportOptions, import_step};

fn unit(v: [f64; 3]) -> [f64; 3] {
    let l = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    [v[0] / l, v[1] / l, v[2] / l]
}

fn perpendicular(n: [f64; 3]) -> [f64; 3] {
    let a = if n[0].abs() < 0.9 {
        [1.0, 0.0, 0.0]
    } else {
        [0.0, 1.0, 0.0]
    };
    let d = a[0] * n[0] + a[1] * n[1] + a[2] * n[2];
    unit([a[0] - d * n[0], a[1] - d * n[1], a[2] - d * n[2]])
}

/// The part of the sphere (centre origin, radius `r`) on the `+n` side
/// of the plane at signed distance `d` along unit `n`.
fn dome_step(r: f64, n: [f64; 3], d: f64) -> String {
    let n = unit(n);
    let x = perpendicular(n);
    let rim = (r * r - d * d).sqrt();
    let c = [d * n[0], d * n[1], d * n[2]];
    let p0 = [c[0] + rim * x[0], c[1] + rim * x[1], c[2] + rim * x[2]];
    let f = |v: [f64; 3]| format!("({:.17e},{:.17e},{:.17e})", v[0], v[1], v[2]);
    format!(
        "ISO-10303-21;\nHEADER;\nFILE_DESCRIPTION(('dome'),'2;1');\nFILE_NAME('dome.step','',(''),(''),'','','');\nFILE_SCHEMA(('AUTOMOTIVE_DESIGN'));\nENDSEC;\nDATA;\n\
#1=CARTESIAN_POINT('',{c});\n#2=DIRECTION('',{n});\n#3=DIRECTION('',{x});\n#4=AXIS2_PLACEMENT_3D('',#1,#2,#3);\n\
#5=CIRCLE('',#4,{rim:.17e});\n#6=CARTESIAN_POINT('',{p0});\n#7=VERTEX_POINT('',#6);\n#8=EDGE_CURVE('',#7,#7,#5,.T.);\n\
#9=CARTESIAN_POINT('',(0.,0.,0.));\n#10=DIRECTION('',(0.,0.,1.));\n#11=DIRECTION('',(1.,0.,0.));\n#12=AXIS2_PLACEMENT_3D('',#9,#10,#11);\n\
#13=SPHERICAL_SURFACE('',#12,{r:.17e});\n#14=ORIENTED_EDGE('',*,*,#8,.T.);\n#15=EDGE_LOOP('',(#14));\n#16=FACE_OUTER_BOUND('',#15,.T.);\n\
#17=ADVANCED_FACE('',(#16),#13,.T.);\n#18=DIRECTION('',{m});\n#19=AXIS2_PLACEMENT_3D('',#1,#18,#3);\n#20=PLANE('',#19);\n\
#21=ORIENTED_EDGE('',*,*,#8,.F.);\n#22=EDGE_LOOP('',(#21));\n#23=FACE_OUTER_BOUND('',#22,.T.);\n#24=ADVANCED_FACE('',(#23),#20,.T.);\n\
#25=CLOSED_SHELL('',(#17,#24));\n#26=MANIFOLD_SOLID_BREP('',#25);\n\
#27=(LENGTH_UNIT() NAMED_UNIT(*) SI_UNIT(.MILLI.,.METRE.));\n#28=(NAMED_UNIT(*) PLANE_ANGLE_UNIT() SI_UNIT($,.RADIAN.));\n\
#29=(GEOMETRIC_REPRESENTATION_CONTEXT(3) GLOBAL_UNIT_ASSIGNED_CONTEXT((#27,#28)) REPRESENTATION_CONTEXT('',''));\n\
#30=ADVANCED_BREP_SHAPE_REPRESENTATION('',(#26),#29);\nENDSEC;\nEND-ISO-10303-21;\n",
        c = f(c),
        n = f(n),
        x = f(x),
        p0 = f(p0),
        m = f([-n[0], -n[1], -n[2]]),
    )
}

fn mesh_volume(mesh: &popkern_mesh::TriMesh) -> f64 {
    mesh.triangles
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
        .abs()
}

fn check_dome(n: [f64; 3], d: f64) {
    let r = 10.0;
    let text = dome_step(r, n, d);
    let imported = import_step(text.as_bytes(), &StepImportOptions::default()).expect("import");
    let solids: Vec<_> = imported
        .document
        .products
        .iter()
        .flat_map(|p| p.shapes.iter())
        .flat_map(|s| s.solids.iter().cloned())
        .collect();
    assert_eq!(
        solids.len(),
        1,
        "n {n:?} d {d}: one solid ({:?})",
        imported.diagnostics.entries()
    );
    let dome = &solids[0];
    dome.validate().expect("valid");
    let mesh = popkern_mesh::tessellate(dome, &popkern_mesh::TessellationOptions::default())
        .expect("mesh");
    assert_eq!(
        popkern_mesh::welded_open_edges(&mesh),
        0,
        "n {n:?} d {d}: watertight"
    );
    let h = r - d;
    let exact = std::f64::consts::PI * h * h * (3.0 * r - h) / 3.0;
    let got = mesh_volume(&mesh);
    assert!(
        (got - exact).abs() < 2e-3 * exact,
        "n {n:?} d {d}: volume {got} vs {exact}"
    );
}

#[test]
fn a_cap_around_the_north_pole() {
    check_dome([0.0, 0.0, 1.0], 4.0);
}

#[test]
fn a_cap_larger_than_a_hemisphere() {
    check_dome([0.0, 0.0, 1.0], -4.0);
}

#[test]
fn a_cap_around_the_south_pole() {
    check_dome([0.0, 0.0, -1.0], 4.0);
}

#[test]
fn a_tilted_cap_whose_rim_is_no_latitude() {
    check_dome([1.0, 0.0, 2.0], 3.0);
    check_dome([1.0, 0.5, -2.0], -2.0);
}

#[test]
fn a_cap_clear_of_both_poles() {
    check_dome([1.0, 0.0, 0.0], 4.0);
}
