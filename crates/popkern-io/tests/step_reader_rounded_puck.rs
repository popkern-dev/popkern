//! PK-648: torus fillet bands bounded by two major circles (rounded
//! disk edges, a cup's floor fillet, a thin washer's rounded bore)
//! mesh with the winding their face senses name.
//! The band between the torus's outer equator and its bottom crosses
//! the chart's v wrap, so the band slit lifts a ring a period; the
//! lifted band must still wind with its neighbours (measured: turbine
//! .stp's two fillet bands wound against their planes, every face
//! area right, the solid 66 % over OpenCASCADE's volume).

use popkern_io::{StepImportOptions, import_step};

fn f(v: [f64; 3]) -> String {
    format!("({:.17e},{:.17e},{:.17e})", v[0], v[1], v[2])
}

/// A disk of radius `big + small` and height `h`, both circular edges
/// rounded by radius `small` (two convex torus quarter bands).
fn puck_step(big: f64, small: f64, h: f64) -> String {
    let circle = |id: usize, rho: f64, z: f64| {
        format!(
            "#{id}=CARTESIAN_POINT('',(0.,0.,{z:.17e}));\n#{a}=AXIS2_PLACEMENT_3D('',#{id},#2,#3);\n#{c}=CIRCLE('',#{a},{rho:.17e});\n#{p}=CARTESIAN_POINT('',{pt});\n#{v}=VERTEX_POINT('',#{p});\n#{e}=EDGE_CURVE('',#{v},#{v},#{c},.T.);\n",
            a = id + 1,
            c = id + 2,
            p = id + 3,
            v = id + 4,
            e = id + 5,
            pt = f([rho, 0.0, z]),
        )
    };
    // Edges: top rim #15 (rho big, z h), top equator #25 (rho big+small,
    // z h-small), bottom equator #35 (rho big+small, z small), bottom
    // rim #45 (rho big, z 0).
    let mut text = String::from(
        "ISO-10303-21;\nHEADER;\nFILE_DESCRIPTION(('puck'),'2;1');\nFILE_NAME('puck.step','',(''),(''),'','','');\nFILE_SCHEMA(('AUTOMOTIVE_DESIGN'));\nENDSEC;\nDATA;\n\
#1=CARTESIAN_POINT('',(0.,0.,0.));\n#2=DIRECTION('',(0.,0.,1.));\n#3=DIRECTION('',(1.,0.,0.));\n#4=DIRECTION('',(0.,0.,-1.));\n",
    );
    text += &circle(10, big, h);
    text += &circle(20, big + small, h - small);
    text += &circle(30, big + small, small);
    text += &circle(40, big, 0.0);
    let face = |id: usize, surface: &str, bounds: &[(usize, &str)], sense: &str| {
        let mut out = String::new();
        let mut refs = Vec::new();
        for (k, (edge, dir)) in bounds.iter().enumerate() {
            let (oe, lp, fb) = (id + 3 * k + 1, id + 3 * k + 2, id + 3 * k + 3);
            out += &format!(
                "#{oe}=ORIENTED_EDGE('',*,*,#{edge},{dir});\n#{lp}=EDGE_LOOP('',(#{oe}));\n#{fb}=FACE_BOUND('',#{lp},.T.);\n"
            );
            refs.push(format!("#{fb}"));
        }
        out += &format!(
            "#{id}=ADVANCED_FACE('',({}),{surface},{sense});\n",
            refs.join(",")
        );
        out
    };
    // Surfaces.
    text += &format!(
        "#50=PLANE('',#11);\n#51=CARTESIAN_POINT('',(0.,0.,{tz:.17e}));\n#52=AXIS2_PLACEMENT_3D('',#51,#2,#3);\n#53=TOROIDAL_SURFACE('',#52,{big:.17e},{small:.17e});\n\
#54=AXIS2_PLACEMENT_3D('',#1,#2,#3);\n#55=CYLINDRICAL_SURFACE('',#54,{rc:.17e});\n\
#56=CARTESIAN_POINT('',(0.,0.,{small:.17e}));\n#57=AXIS2_PLACEMENT_3D('',#56,#2,#3);\n#58=TOROIDAL_SURFACE('',#57,{big:.17e},{small:.17e});\n\
#59=AXIS2_PLACEMENT_3D('',#1,#4,#3);\n#60=PLANE('',#59);\n",
        tz = h - small,
        rc = big + small,
    );
    text += &face(100, "#50", &[(15, ".T.")], ".T.");
    text += &face(110, "#53", &[(15, ".F."), (25, ".T.")], ".T.");
    text += &face(120, "#55", &[(25, ".F."), (35, ".T.")], ".T.");
    text += &face(130, "#58", &[(35, ".F."), (45, ".T.")], ".T.");
    text += &face(140, "#60", &[(45, ".F.")], ".T.");
    text += "#150=CLOSED_SHELL('',(#100,#110,#120,#130,#140));\n#151=MANIFOLD_SOLID_BREP('',#150);\n\
#152=(LENGTH_UNIT() NAMED_UNIT(*) SI_UNIT(.MILLI.,.METRE.));\n#153=(NAMED_UNIT(*) PLANE_ANGLE_UNIT() SI_UNIT($,.RADIAN.));\n\
#154=(GEOMETRIC_REPRESENTATION_CONTEXT(3) GLOBAL_UNIT_ASSIGNED_CONTEXT((#152,#153)) REPRESENTATION_CONTEXT('',''));\n\
#155=ADVANCED_BREP_SHAPE_REPRESENTATION('',(#151),#154);\nENDSEC;\nEND-ISO-10303-21;\n";
    text
}

/// A cup: outer radius `outer`, height `h`, floor at `t`, bore radius
/// `a + r` whose bottom corner is a CONCAVE fillet of radius `r` (a
/// torus band from the tube's bottom to its outer equator, sense .F.:
/// turbine.stp's shroud fillets).
fn cup_step(outer: f64, h: f64, t: f64, a: f64, r: f64) -> String {
    let circle = |id: usize, rho: f64, z: f64| {
        format!(
            "#{id}=CARTESIAN_POINT('',(0.,0.,{z:.17e}));\n#{a}=AXIS2_PLACEMENT_3D('',#{id},#2,#3);\n#{c}=CIRCLE('',#{a},{rho:.17e});\n#{p}=CARTESIAN_POINT('',{pt});\n#{v}=VERTEX_POINT('',#{p});\n#{e}=EDGE_CURVE('',#{v},#{v},#{c},.T.);\n",
            a = id + 1,
            c = id + 2,
            p = id + 3,
            v = id + 4,
            e = id + 5,
            pt = f([rho, 0.0, z]),
        )
    };
    let mut text = String::from(
        "ISO-10303-21;\nHEADER;\nFILE_DESCRIPTION(('cup'),'2;1');\nFILE_NAME('cup.step','',(''),(''),'','','');\nFILE_SCHEMA(('AUTOMOTIVE_DESIGN'));\nENDSEC;\nDATA;\n\
#1=CARTESIAN_POINT('',(0.,0.,0.));\n#2=DIRECTION('',(0.,0.,1.));\n#3=DIRECTION('',(1.,0.,0.));\n#4=DIRECTION('',(0.,0.,-1.));\n",
    );
    // C1 bottom outer #15, C2 top outer #25, C3 top inner #35, C4 bore
    // bottom (fillet top) #45, C5 floor rim (fillet bottom) #55.
    text += &circle(10, outer, 0.0);
    text += &circle(20, outer, h);
    text += &circle(30, a + r, h);
    text += &circle(40, a + r, t + r);
    text += &circle(50, a, t);
    let face = |id: usize, surface: &str, bounds: &[(usize, &str)], sense: &str| {
        let mut out = String::new();
        let mut refs = Vec::new();
        for (k, (edge, dir)) in bounds.iter().enumerate() {
            let (oe, lp, fb) = (id + 3 * k + 1, id + 3 * k + 2, id + 3 * k + 3);
            out += &format!(
                "#{oe}=ORIENTED_EDGE('',*,*,#{edge},{dir});\n#{lp}=EDGE_LOOP('',(#{oe}));\n#{fb}=FACE_BOUND('',#{lp},.T.);\n"
            );
            refs.push(format!("#{fb}"));
        }
        out += &format!(
            "#{id}=ADVANCED_FACE('',({}),{surface},{sense});\n",
            refs.join(",")
        );
        out
    };
    text += &format!(
        "#60=AXIS2_PLACEMENT_3D('',#1,#4,#3);\n#61=PLANE('',#60);\n\
#62=AXIS2_PLACEMENT_3D('',#1,#2,#3);\n#63=CYLINDRICAL_SURFACE('',#62,{outer:.17e});\n\
#64=PLANE('',#21);\n#65=CYLINDRICAL_SURFACE('',#62,{bore:.17e});\n\
#66=CARTESIAN_POINT('',(0.,0.,{tz:.17e}));\n#67=AXIS2_PLACEMENT_3D('',#66,#2,#3);\n#68=TOROIDAL_SURFACE('',#67,{a:.17e},{r:.17e});\n\
#69=PLANE('',#51);\n",
        bore = a + r,
        tz = t + r,
    );
    text += &face(100, "#61", &[(15, ".F.")], ".T.");
    text += &face(110, "#63", &[(15, ".T."), (25, ".F.")], ".T.");
    text += &face(120, "#64", &[(25, ".T."), (35, ".F.")], ".T.");
    text += &face(130, "#65", &[(35, ".T."), (45, ".F.")], ".F.");
    text += &face(140, "#68", &[(45, ".T."), (55, ".F.")], ".F.");
    text += &face(150, "#69", &[(55, ".T.")], ".T.");
    text += "#160=CLOSED_SHELL('',(#100,#110,#120,#130,#140,#150));\n#161=MANIFOLD_SOLID_BREP('',#160);\n\
#162=(LENGTH_UNIT() NAMED_UNIT(*) SI_UNIT(.MILLI.,.METRE.));\n#163=(NAMED_UNIT(*) PLANE_ANGLE_UNIT() SI_UNIT($,.RADIAN.));\n\
#164=(GEOMETRIC_REPRESENTATION_CONTEXT(3) GLOBAL_UNIT_ASSIGNED_CONTEXT((#162,#163)) REPRESENTATION_CONTEXT('',''));\n\
#165=ADVANCED_BREP_SHAPE_REPRESENTATION('',(#161),#164);\nENDSEC;\nEND-ISO-10303-21;\n";
    text
}

/// A thin washer: outer radius `a + w`, height `h`, its bore of radius
/// `a - r` rounded at both ends by convex fillets of radius `r` that
/// meet the top and bottom annuli (width `w`) at radius `a`. With `w`
/// under a thousandth of the rim's circumference the annuli's two rims
/// sit inside the old orientation vote's polyline tolerance (turbine).
fn washer_step(a: f64, w: f64, r: f64, h: f64) -> String {
    let circle = |id: usize, rho: f64, z: f64| {
        format!(
            "#{id}=CARTESIAN_POINT('',(0.,0.,{z:.17e}));\n#{a}=AXIS2_PLACEMENT_3D('',#{id},#2,#3);\n#{c}=CIRCLE('',#{a},{rho:.17e});\n#{p}=CARTESIAN_POINT('',{pt});\n#{v}=VERTEX_POINT('',#{p});\n#{e}=EDGE_CURVE('',#{v},#{v},#{c},.T.);\n",
            a = id + 1,
            c = id + 2,
            p = id + 3,
            v = id + 4,
            e = id + 5,
            pt = f([rho, 0.0, z]),
        )
    };
    let mut text = String::from(
        "ISO-10303-21;\nHEADER;\nFILE_DESCRIPTION(('washer'),'2;1');\nFILE_NAME('washer.step','',(''),(''),'','','');\nFILE_SCHEMA(('AUTOMOTIVE_DESIGN'));\nENDSEC;\nDATA;\n\
#1=CARTESIAN_POINT('',(0.,0.,0.));\n#2=DIRECTION('',(0.,0.,1.));\n#3=DIRECTION('',(1.,0.,0.));\n#4=DIRECTION('',(0.,0.,-1.));\n",
    );
    // Top outer #15, top inner #25, bore top #35, bore bottom #45,
    // bottom inner #55, bottom outer #65.
    text += &circle(10, a + w, h);
    text += &circle(20, a, h);
    text += &circle(30, a - r, h - r);
    text += &circle(40, a - r, r);
    text += &circle(50, a, 0.0);
    text += &circle(60, a + w, 0.0);
    let face = |id: usize, surface: &str, bounds: &[(usize, &str)], sense: &str| {
        let mut out = String::new();
        let mut refs = Vec::new();
        for (k, (edge, dir)) in bounds.iter().enumerate() {
            let (oe, lp, fb) = (id + 3 * k + 1, id + 3 * k + 2, id + 3 * k + 3);
            out += &format!(
                "#{oe}=ORIENTED_EDGE('',*,*,#{edge},{dir});\n#{lp}=EDGE_LOOP('',(#{oe}));\n#{fb}=FACE_BOUND('',#{lp},.T.);\n"
            );
            refs.push(format!("#{fb}"));
        }
        out += &format!(
            "#{id}=ADVANCED_FACE('',({}),{surface},{sense});\n",
            refs.join(",")
        );
        out
    };
    text += &format!(
        "#70=PLANE('',#11);\n#71=AXIS2_PLACEMENT_3D('',#1,#2,#3);\n#72=CYLINDRICAL_SURFACE('',#71,{outer:.17e});\n\
#73=AXIS2_PLACEMENT_3D('',#1,#4,#3);\n#74=PLANE('',#73);\n\
#75=CARTESIAN_POINT('',(0.,0.,{tz:.17e}));\n#76=AXIS2_PLACEMENT_3D('',#75,#2,#3);\n#77=TOROIDAL_SURFACE('',#76,{a:.17e},{r:.17e});\n\
#78=CYLINDRICAL_SURFACE('',#71,{bore:.17e});\n\
#79=CARTESIAN_POINT('',(0.,0.,{r:.17e}));\n#80=AXIS2_PLACEMENT_3D('',#79,#2,#3);\n#81=TOROIDAL_SURFACE('',#80,{a:.17e},{r:.17e});\n",
        outer = a + w,
        tz = h - r,
        bore = a - r,
    );
    text += &face(100, "#70", &[(15, ".T."), (25, ".F.")], ".T.");
    text += &face(110, "#72", &[(15, ".F."), (65, ".T.")], ".T.");
    text += &face(120, "#74", &[(65, ".F."), (55, ".T.")], ".T.");
    text += &face(130, "#77", &[(25, ".T."), (35, ".F.")], ".T.");
    text += &face(140, "#78", &[(35, ".T."), (45, ".F.")], ".F.");
    text += &face(150, "#81", &[(45, ".T."), (55, ".F.")], ".T.");
    text += "#160=CLOSED_SHELL('',(#100,#110,#120,#130,#140,#150));\n#161=MANIFOLD_SOLID_BREP('',#160);\n\
#162=(LENGTH_UNIT() NAMED_UNIT(*) SI_UNIT(.MILLI.,.METRE.));\n#163=(NAMED_UNIT(*) PLANE_ANGLE_UNIT() SI_UNIT($,.RADIAN.));\n\
#164=(GEOMETRIC_REPRESENTATION_CONTEXT(3) GLOBAL_UNIT_ASSIGNED_CONTEXT((#162,#163)) REPRESENTATION_CONTEXT('',''));\n\
#165=ADVANCED_BREP_SHAPE_REPRESENTATION('',(#161),#164);\nENDSEC;\nEND-ISO-10303-21;\n";
    text
}

#[test]
fn a_thin_washer_with_a_rounded_bore_keeps_its_bore_strip_outward() {
    use std::f64::consts::PI;
    let (a, w, r, h) = (100.0, 0.5, 1.5, 6.0);
    let got = signed_meshed_volume(&washer_step(a, w, r, h));
    let rect = (w + r) * h * 0.5 * (2.0 * a - r + w);
    let spandrel = r * r * (a - 0.5 * r) - (PI * r * r / 4.0) * (a - 4.0 * r / (3.0 * PI));
    let exact = 2.0 * PI * (rect - 2.0 * spandrel);
    assert!(
        (got - exact).abs() < 2e-3 * exact,
        "signed volume {got} vs {exact} (outward winding)"
    );
}

/// Import, mesh, require watertight and every shared mesh edge walked
/// once each way, and return the signed (outward-positive) volume.
fn signed_meshed_volume(text: &str) -> f64 {
    let imported = import_step(text.as_bytes(), &StepImportOptions::default()).expect("import");
    let solid = imported
        .document
        .products
        .iter()
        .flat_map(|p| p.shapes.iter())
        .flat_map(|s| s.solids.iter().cloned())
        .next()
        .expect("solid");
    solid.validate().expect("valid");
    let mesh = popkern_mesh::tessellate(&solid, &popkern_mesh::TessellationOptions::default())
        .expect("mesh");
    assert_eq!(popkern_mesh::welded_open_edges(&mesh), 0, "watertight");
    let key = |p: popkern_core::geometry::Point3| {
        (
            (p.x * 1e7).round() as i64,
            (p.y * 1e7).round() as i64,
            (p.z * 1e7).round() as i64,
        )
    };
    let mut directed = std::collections::HashMap::new();
    for t in &mesh.triangles {
        let c = t.map(|i| key(mesh.vertices[i].position));
        for (a, b) in [(c[0], c[1]), (c[1], c[2]), (c[2], c[0])] {
            if a != b {
                *directed.entry((a, b)).or_insert(0usize) += 1;
            }
        }
    }
    let same_way = directed.values().filter(|&&n| n > 1).count();
    assert_eq!(
        same_way, 0,
        "{same_way} mesh edges walked the same way twice"
    );
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
        .sum()
}

#[test]
fn a_cup_whose_floor_meets_its_bore_in_a_concave_fillet() {
    use std::f64::consts::PI;
    let (outer, h, t, a, r) = (14.0, 8.0, 2.0, 10.0, 1.5);
    let got = signed_meshed_volume(&cup_step(outer, h, t, a, r));
    let square = 2.0 * PI * (a + 0.5 * r) * r * r;
    let quarter = 2.0 * PI * (a + 4.0 * r / (3.0 * PI)) * (PI * r * r / 4.0);
    let exact = PI * outer * outer * h - PI * (a + r).powi(2) * (h - t) + (square - quarter);
    assert!(
        (got - exact).abs() < 2e-3 * exact,
        "signed volume {got} vs {exact} (outward winding)"
    );
}

#[test]
fn a_puck_with_both_edges_rounded_winds_consistently() {
    let (big, small, h) = (10.0, 1.5, 6.0);
    let imported = import_step(
        puck_step(big, small, h).as_bytes(),
        &StepImportOptions::default(),
    )
    .expect("import");
    let puck = imported
        .document
        .products
        .iter()
        .flat_map(|p| p.shapes.iter())
        .flat_map(|s| s.solids.iter().cloned())
        .next()
        .expect("solid");
    puck.validate().expect("valid");
    let mesh = popkern_mesh::tessellate(&puck, &popkern_mesh::TessellationOptions::default())
        .expect("mesh");
    assert_eq!(popkern_mesh::welded_open_edges(&mesh), 0, "watertight");
    // Every shared mesh edge is walked once each way.
    let key = |p: popkern_core::geometry::Point3| {
        (
            (p.x * 1e7).round() as i64,
            (p.y * 1e7).round() as i64,
            (p.z * 1e7).round() as i64,
        )
    };
    let mut directed = std::collections::HashMap::new();
    for t in &mesh.triangles {
        let c = t.map(|i| key(mesh.vertices[i].position));
        for (a, b) in [(c[0], c[1]), (c[1], c[2]), (c[2], c[0])] {
            if a != b {
                *directed.entry((a, b)).or_insert(0usize) += 1;
            }
        }
    }
    let same_way = directed.values().filter(|&&n| n > 1).count();
    assert_eq!(
        same_way, 0,
        "{same_way} mesh edges walked the same way twice"
    );
    let signed: f64 = mesh
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
    use std::f64::consts::PI;
    let rounded_layer = PI * big * big * small
        + 2.0 * PI * (big + 4.0 * small / (3.0 * PI)) * (PI * small * small / 4.0);
    let exact = PI * (big + small).powi(2) * (h - 2.0 * small) + 2.0 * rounded_layer;
    assert!(
        (signed - exact).abs() < 2e-3 * exact,
        "signed volume {signed} vs {exact} (outward winding)"
    );
}

/// A shell written inside out (every face sense and loop reversed) is
/// still a closed body: the mesh faces outward. The coherence pass
/// flips a closed component on its signed volume; judging closure by
/// the face meshes' boundary segments read such bodies as open (their
/// meshes still carry T-junctions there) and left them inside out
/// (jack_in_the_box's pins, io1-tc; PK-648's corpus volume gate).
#[test]
fn an_inside_out_puck_meshes_outward() {
    let (big, small, h) = (10.0, 1.5, 6.0);
    let text: String = puck_step(big, small, h)
        .lines()
        .map(|line| {
            if line.contains("ORIENTED_EDGE") || line.contains("ADVANCED_FACE") {
                line.replace(".T.", ".X.")
                    .replace(".F.", ".T.")
                    .replace(".X.", ".F.")
            } else {
                line.to_string()
            }
        })
        .collect::<Vec<_>>()
        .join("\n");
    let got = signed_meshed_volume(&text);
    use std::f64::consts::PI;
    let rounded_layer = PI * big * big * small
        + 2.0 * PI * (big + 4.0 * small / (3.0 * PI)) * (PI * small * small / 4.0);
    let exact = PI * (big + small).powi(2) * (h - 2.0 * small) + 2.0 * rounded_layer;
    assert!(
        (got - exact).abs() < 2e-3 * exact,
        "signed volume {got} vs {exact} (outward winding)"
    );
}
