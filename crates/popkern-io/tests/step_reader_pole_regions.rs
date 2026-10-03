//! PK-646: faces that reach a chart pole, bounded the way third-party
//! writers bound them, mesh as the region they name.
//!
//! - A sphere LUNE (the solid between two half-planes through the
//!   sphere's axis): its sphere face is bounded by two meridians from
//!   pole to pole, and the chart strip between them is the face or
//!   its complement as the flag says (measured: jack body 33, a ball
//!   with a 2.85° wedge cut out, meshed the wedge).
//! - A cone tip cut by a TILTED plane: the cone face's rim is no
//!   latitude, and the pole cap's apex wall must not evacuate the
//!   chart's interior (measured: jack body 28 meshed 33.2 of 42.1).

use popkern_io::{StepImportOptions, import_step};

fn f(v: [f64; 3]) -> String {
    format!("({:.17e},{:.17e},{:.17e})", v[0], v[1], v[2])
}

const HEADER: &str = "ISO-10303-21;\nHEADER;\nFILE_DESCRIPTION(('pole'),'2;1');\nFILE_NAME('pole.step','',(''),(''),'','','');\nFILE_SCHEMA(('AUTOMOTIVE_DESIGN'));\nENDSEC;\nDATA;\n";
const FOOTER: &str = "#90=(LENGTH_UNIT() NAMED_UNIT(*) SI_UNIT(.MILLI.,.METRE.));\n#91=(NAMED_UNIT(*) PLANE_ANGLE_UNIT() SI_UNIT($,.RADIAN.));\n\
#92=(GEOMETRIC_REPRESENTATION_CONTEXT(3) GLOBAL_UNIT_ASSIGNED_CONTEXT((#90,#91)) REPRESENTATION_CONTEXT('',''));\n\
#93=ADVANCED_BREP_SHAPE_REPRESENTATION('',(#81),#92);\nENDSEC;\nEND-ISO-10303-21;\n";

/// The solid of the sphere (centre origin, axis +Z, radius `r`) between
/// azimuths `a` and `a + phi`: one sphere face bounded by the two
/// meridian half-circles, two half-disks sharing the axis segment.
fn lune_step(r: f64, a: f64, phi: f64) -> String {
    let b = a + phi;
    let d = |t: f64| [t.cos(), t.sin(), 0.0];
    let axis_of = |t: f64| [-t.sin(), t.cos(), 0.0]; // z × d(t)
    format!(
        "{HEADER}\
#1=CARTESIAN_POINT('',(0.,0.,0.));\n#2=DIRECTION('',(0.,0.,1.));\n#3=DIRECTION('',(1.,0.,0.));\n#4=AXIS2_PLACEMENT_3D('',#1,#2,#3);\n\
#5=SPHERICAL_SURFACE('',#4,{r:.17e});\n\
#6=CARTESIAN_POINT('',(0.,0.,{r:.17e}));\n#7=VERTEX_POINT('',#6);\n#8=CARTESIAN_POINT('',(0.,0.,{nr:.17e}));\n#9=VERTEX_POINT('',#8);\n\
#10=DIRECTION('',{aa});\n#11=AXIS2_PLACEMENT_3D('',#1,#10,#2);\n#12=CIRCLE('',#11,{r:.17e});\n#13=EDGE_CURVE('',#7,#9,#12,.T.);\n\
#20=DIRECTION('',{ab});\n#21=AXIS2_PLACEMENT_3D('',#1,#20,#2);\n#22=CIRCLE('',#21,{r:.17e});\n#23=EDGE_CURVE('',#7,#9,#22,.T.);\n\
#30=VECTOR('',#2,1.);\n#31=LINE('',#8,#30);\n#32=EDGE_CURVE('',#9,#7,#31,.T.);\n\
#40=ORIENTED_EDGE('',*,*,#13,.T.);\n#41=ORIENTED_EDGE('',*,*,#23,.F.);\n#42=EDGE_LOOP('',(#40,#41));\n#43=FACE_OUTER_BOUND('',#42,.T.);\n#44=ADVANCED_FACE('',(#43),#5,.T.);\n\
#50=DIRECTION('',{na});\n#51=DIRECTION('',{da});\n#52=AXIS2_PLACEMENT_3D('',#1,#50,#51);\n#53=PLANE('',#52);\n\
#54=ORIENTED_EDGE('',*,*,#13,.F.);\n#55=ORIENTED_EDGE('',*,*,#32,.F.);\n#56=EDGE_LOOP('',(#54,#55));\n#57=FACE_OUTER_BOUND('',#56,.T.);\n#58=ADVANCED_FACE('',(#57),#53,.T.);\n\
#60=DIRECTION('',{nb});\n#61=DIRECTION('',{db});\n#62=AXIS2_PLACEMENT_3D('',#1,#60,#61);\n#63=PLANE('',#62);\n\
#64=ORIENTED_EDGE('',*,*,#23,.T.);\n#65=ORIENTED_EDGE('',*,*,#32,.T.);\n#66=EDGE_LOOP('',(#64,#65));\n#67=FACE_OUTER_BOUND('',#66,.T.);\n#68=ADVANCED_FACE('',(#67),#63,.T.);\n\
#80=CLOSED_SHELL('',(#44,#58,#68));\n#81=MANIFOLD_SOLID_BREP('',#80);\n{FOOTER}",
        nr = -r,
        aa = f(axis_of(a)),
        ab = f(axis_of(b)),
        na = f([a.sin(), -a.cos(), 0.0]),
        da = f(d(a)),
        nb = f([-b.sin(), b.cos(), 0.0]),
        db = f(d(b)),
    )
}

/// The tip of the cone (apex origin, axis +Z, half-angle `alpha`) under
/// the plane through (0, 0, `h`) tilted by `beta` about Y. Returns the
/// file and the exact volume.
fn tilted_tip_step(alpha: f64, h: f64, beta: f64) -> (String, f64) {
    let n = [beta.sin(), 0.0, beta.cos()];
    let generator = |t: f64| [alpha.sin() * t.cos(), alpha.sin() * t.sin(), alpha.cos()];
    let hit = |t: f64| {
        let g = generator(t);
        let s = h * beta.cos() / (n[0] * g[0] + n[2] * g[2]);
        [s * g[0], s * g[1], s * g[2]]
    };
    let (p0, pi) = (hit(0.0), hit(std::f64::consts::PI));
    let c = [0.5 * (p0[0] + pi[0]), 0.0, 0.5 * (p0[2] + pi[2])];
    let major = 0.5 * ((p0[0] - pi[0]).powi(2) + (p0[2] - pi[2]).powi(2)).sqrt();
    let dir = [
        (p0[0] - pi[0]) / (2.0 * major),
        0.0,
        (p0[2] - pi[2]) / (2.0 * major),
    ];
    let minor = ((c[2] * alpha.tan()).powi(2) - c[0] * c[0]).sqrt();
    let radius = h * alpha.tan();
    let text = format!(
        "{HEADER}\
#1=CARTESIAN_POINT('',(0.,0.,{h:.17e}));\n#2=DIRECTION('',(0.,0.,1.));\n#3=DIRECTION('',(1.,0.,0.));\n#4=AXIS2_PLACEMENT_3D('',#1,#2,#3);\n\
#5=CONICAL_SURFACE('',#4,{radius:.17e},{alpha:.17e});\n\
#10=CARTESIAN_POINT('',{c});\n#11=DIRECTION('',{n});\n#12=DIRECTION('',{dir});\n#13=AXIS2_PLACEMENT_3D('',#10,#11,#12);\n\
#14=ELLIPSE('',#13,{major:.17e},{minor:.17e});\n#15=CARTESIAN_POINT('',{p0});\n#16=VERTEX_POINT('',#15);\n#17=EDGE_CURVE('',#16,#16,#14,.T.);\n\
#20=ORIENTED_EDGE('',*,*,#17,.F.);\n#21=EDGE_LOOP('',(#20));\n#22=FACE_OUTER_BOUND('',#21,.T.);\n#23=ADVANCED_FACE('',(#22),#5,.T.);\n\
#30=PLANE('',#13);\n#31=ORIENTED_EDGE('',*,*,#17,.T.);\n#32=EDGE_LOOP('',(#31));\n#33=FACE_OUTER_BOUND('',#32,.T.);\n#34=ADVANCED_FACE('',(#33),#30,.T.);\n\
#80=CLOSED_SHELL('',(#23,#34));\n#81=MANIFOLD_SOLID_BREP('',#80);\n{FOOTER}",
        c = f(c),
        n = f(n),
        dir = f(dir),
        p0 = f(p0),
    );
    let exact = std::f64::consts::PI * major * minor * (h * beta.cos()) / 3.0;
    (text, exact)
}

/// Import the file's one solid, require it valid and watertight, and
/// return its meshed volume and centroid.
fn meshed_volume(text: &str, label: &str) -> (f64, [f64; 3]) {
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
    assert_eq!(
        popkern_mesh::welded_open_edges(&mesh),
        0,
        "{label}: watertight"
    );
    // Signed tetrahedra from the origin: volume and first moment.
    let (mut volume, mut moment) = (0.0, [0.0; 3]);
    for t in &mesh.triangles {
        let (a, b, c) = (
            mesh.vertices[t[0]].position,
            mesh.vertices[t[1]].position,
            mesh.vertices[t[2]].position,
        );
        let v = (a.x * (b.y * c.z - b.z * c.y) - a.y * (b.x * c.z - b.z * c.x)
            + a.z * (b.x * c.y - b.y * c.x))
            / 6.0;
        volume += v;
        moment[0] += v * (a.x + b.x + c.x) / 4.0;
        moment[1] += v * (a.y + b.y + c.y) / 4.0;
        moment[2] += v * (a.z + b.z + c.z) / 4.0;
    }
    (volume.abs(), moment.map(|m| m / volume))
}

/// The sphere chart's interior lattice is capped at `max_subdivisions`
/// (PK-647): faces covering most of a radius-10 sphere facet past the
/// chord height and mesh up to 2e-3 short.
const LARGE_SPHERE_FACE_BOUND: f64 = 3e-3;

fn check_lune(a: f64, phi: f64) {
    let r = 10.0;
    let label = format!("lune a {a} phi {phi}");
    let (got, centroid) = meshed_volume(&lune_step(r, a, phi), &label);
    let exact = phi * 2.0 / 3.0 * r * r * r;
    assert!(
        (got - exact).abs() < LARGE_SPHERE_FACE_BOUND * exact,
        "{label}: volume {got} vs {exact}"
    );
    // The solid's centroid lies on the bisecting half-plane: this names
    // the side even when the two sides have equal volume (φ = π).
    let azimuth = centroid[1].atan2(centroid[0]);
    let off = (azimuth - (a + 0.5 * phi)).rem_euclid(std::f64::consts::TAU);
    assert!(
        off.min(std::f64::consts::TAU - off) < 1e-2,
        "{label}: centroid azimuth {azimuth} vs bisector {}",
        a + 0.5 * phi
    );
}

#[test]
fn a_thin_lune_is_the_strip_between_its_meridians() {
    check_lune(1.0, 0.3);
}

#[test]
fn a_ball_with_a_thin_wedge_cut_out_is_the_complement_of_that_strip() {
    // Azimuths 1.0 → 2π + 0.7: the face spans the chart's u seam and
    // its meridians bound a 0.3-wide chart strip that is NOT the face.
    check_lune(1.0, std::f64::consts::TAU - 0.3);
}

#[test]
fn a_half_ball_cut_through_its_axis() {
    check_lune(0.2, std::f64::consts::PI);
}

#[test]
fn a_cone_tip_under_a_tilted_cut() {
    for beta in [0.0, 0.35, 0.8] {
        let (text, exact) = tilted_tip_step(0.3, 10.0, beta);
        let (got, _) = meshed_volume(&text, &format!("tilt {beta}"));
        assert!(
            (got - exact).abs() < 2e-3 * exact,
            "tilt {beta}: volume {got} vs {exact}"
        );
    }
}
