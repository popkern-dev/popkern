//! PK-641: a pipe elbow as third-party writers spell it — a torus face
//! bounded by two tube circles (no seam edge) plus two end disks —
//! imports valid, meshes watertight, and encloses θ·R·π·r². The band
//! between two tube circles is one of TWO bands on a torus; the 270°
//! elbow's complementary band is the 90° one, so the side matters.

use popkern_io::{StepImportOptions, import_step};

fn f(v: [f64; 3]) -> String {
    format!("({:.17e},{:.17e},{:.17e})", v[0], v[1], v[2])
}

/// The torus (centre origin, axis +Z, radii `big`, `small`) swept from
/// azimuth `start` through `theta`, capped by disks.
fn elbow_step(big: f64, small: f64, start: f64, theta: f64) -> String {
    let at = |phi: f64| [big * phi.cos(), big * phi.sin(), 0.0];
    let tangent = |phi: f64| [-phi.sin(), phi.cos(), 0.0];
    let radial = |phi: f64| [phi.cos(), phi.sin(), 0.0];
    let rim = |phi: f64| [(big + small) * phi.cos(), (big + small) * phi.sin(), 0.0];
    let neg = |v: [f64; 3]| [-v[0], -v[1], -v[2]];
    format!(
        "ISO-10303-21;\nHEADER;\nFILE_DESCRIPTION(('elbow'),'2;1');\nFILE_NAME('elbow.step','',(''),(''),'','','');\nFILE_SCHEMA(('AUTOMOTIVE_DESIGN'));\nENDSEC;\nDATA;\n\
#1=CARTESIAN_POINT('',(0.,0.,0.));\n#2=DIRECTION('',(0.,0.,1.));\n#3=DIRECTION('',(1.,0.,0.));\n#4=AXIS2_PLACEMENT_3D('',#1,#2,#3);\n\
#5=TOROIDAL_SURFACE('',#4,{big:.17e},{small:.17e});\n\
#10=CARTESIAN_POINT('',{c0});\n#11=DIRECTION('',{t0});\n#12=DIRECTION('',{r0});\n#13=AXIS2_PLACEMENT_3D('',#10,#11,#12);\n#14=CIRCLE('',#13,{small:.17e});\n\
#15=CARTESIAN_POINT('',{p0});\n#16=VERTEX_POINT('',#15);\n#17=EDGE_CURVE('',#16,#16,#14,.T.);\n\
#20=CARTESIAN_POINT('',{c1});\n#21=DIRECTION('',{t1});\n#22=DIRECTION('',{r1});\n#23=AXIS2_PLACEMENT_3D('',#20,#21,#22);\n#24=CIRCLE('',#23,{small:.17e});\n\
#25=CARTESIAN_POINT('',{p1});\n#26=VERTEX_POINT('',#25);\n#27=EDGE_CURVE('',#26,#26,#24,.T.);\n\
#30=ORIENTED_EDGE('',*,*,#17,.T.);\n#31=EDGE_LOOP('',(#30));\n#32=FACE_BOUND('',#31,.T.);\n\
#33=ORIENTED_EDGE('',*,*,#27,.F.);\n#34=EDGE_LOOP('',(#33));\n#35=FACE_BOUND('',#34,.T.);\n\
#36=ADVANCED_FACE('',(#32,#35),#5,.T.);\n\
#40=DIRECTION('',{n0});\n#41=AXIS2_PLACEMENT_3D('',#10,#40,#12);\n#42=PLANE('',#41);\n#43=ORIENTED_EDGE('',*,*,#17,.F.);\n#44=EDGE_LOOP('',(#43));\n#45=FACE_OUTER_BOUND('',#44,.T.);\n#46=ADVANCED_FACE('',(#45),#42,.T.);\n\
#50=AXIS2_PLACEMENT_3D('',#20,#21,#22);\n#51=PLANE('',#50);\n#52=ORIENTED_EDGE('',*,*,#27,.T.);\n#53=EDGE_LOOP('',(#52));\n#54=FACE_OUTER_BOUND('',#53,.T.);\n#55=ADVANCED_FACE('',(#54),#51,.T.);\n\
#60=CLOSED_SHELL('',(#36,#46,#55));\n#61=MANIFOLD_SOLID_BREP('',#60);\n\
#62=(LENGTH_UNIT() NAMED_UNIT(*) SI_UNIT(.MILLI.,.METRE.));\n#63=(NAMED_UNIT(*) PLANE_ANGLE_UNIT() SI_UNIT($,.RADIAN.));\n\
#64=(GEOMETRIC_REPRESENTATION_CONTEXT(3) GLOBAL_UNIT_ASSIGNED_CONTEXT((#62,#63)) REPRESENTATION_CONTEXT('',''));\n\
#65=ADVANCED_BREP_SHAPE_REPRESENTATION('',(#61),#64);\nENDSEC;\nEND-ISO-10303-21;\n",
        c0 = f(at(start)),
        t0 = f(tangent(start)),
        r0 = f(radial(start)),
        p0 = f(rim(start)),
        c1 = f(at(start + theta)),
        t1 = f(tangent(start + theta)),
        r1 = f(radial(start + theta)),
        p1 = f(rim(start + theta)),
        n0 = f(neg(tangent(start))),
    )
}

fn check_elbow(start: f64, theta: f64) {
    let (big, small) = (5.0, 1.0);
    let text = elbow_step(big, small, start, theta);
    let imported = import_step(text.as_bytes(), &StepImportOptions::default()).expect("import");
    let elbow = imported
        .document
        .products
        .iter()
        .flat_map(|p| p.shapes.iter())
        .flat_map(|s| s.solids.iter().cloned())
        .next()
        .expect("solid");
    let mesh = popkern_mesh::tessellate(&elbow, &popkern_mesh::TessellationOptions::default())
        .expect("mesh");
    assert_eq!(
        popkern_mesh::welded_open_edges(&mesh),
        0,
        "theta {theta}: watertight"
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
    // The meshed tube's cross-section is the rim's chord polygon: at
    // chord height 0.005 on a radius-1 tube (32 sides) it encloses
    // 0.64 % less than the circle, so 1 % is the honest band; the
    // complementary band would miss by far more (e.g. 270° vs 90°).
    let exact = theta * big * std::f64::consts::PI * small * small;
    assert!(
        (got - exact).abs() < 1e-2 * exact,
        "theta {theta}: volume {got} vs {exact}"
    );
}

#[test]
fn a_quarter_elbow() {
    check_elbow(0.0, std::f64::consts::FRAC_PI_2);
}

#[test]
fn a_half_elbow() {
    check_elbow(0.0, std::f64::consts::PI);
}

#[test]
fn a_three_quarter_elbow_is_not_its_complement() {
    check_elbow(0.0, 1.5 * std::f64::consts::PI);
}

#[test]
fn elbows_across_the_chart_seam_keep_their_side() {
    // Starting past the torus chart's u seam, one tube circle's major
    // angle wraps below the other's: the naive band is the complement.
    for start in [4.0, 5.5] {
        for theta in [std::f64::consts::FRAC_PI_2, 1.5 * std::f64::consts::PI] {
            check_elbow(start, theta);
        }
    }
}
