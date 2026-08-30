//! Phase 2c (PK-383): the parity matrix. Every substrate row is built
//! once through `popkern-part` and once through the recorder, and the
//! two bodies must agree in census (faces/edges/vertices), volume, and
//! sorted mesh vertices. Plus D3 rename-safety: recorded selections are
//! name snapshots, so renaming features or entity labels changes
//! nothing about what replays.
use popkern_part::Plane as PartPlane;
use popkern_part::prelude::*;
use popkern_script::Plane;
use popkern_script::prelude::*;

fn signature(part: &Part) -> (usize, usize, usize, i64, Vec<[i64; 3]>) {
    let b = part.brep();
    let mesh = part
        .mesh(&popkern_mesh::TessellationOptions::default())
        .expect("mesh");
    let mut v6 = 0.0;
    for t in &mesh.triangles {
        let (a, b, c) = (
            mesh.vertices[t[0]].position,
            mesh.vertices[t[1]].position,
            mesh.vertices[t[2]].position,
        );
        v6 += a.x * (b.y * c.z - b.z * c.y) - a.y * (b.x * c.z - b.z * c.x)
            + a.z * (b.x * c.y - b.y * c.x);
    }
    let mut verts: Vec<[i64; 3]> = mesh
        .vertices
        .iter()
        .map(|v| {
            let p = v.position;
            [
                (p.x * 1e9).round() as i64,
                (p.y * 1e9).round() as i64,
                (p.z * 1e9).round() as i64,
            ]
        })
        .collect();
    verts.sort();
    (
        b.faces.len(),
        b.edges.len(),
        b.vertices.len(),
        (v6 / 6.0 * 1e9).round() as i64,
        verts,
    )
}

fn assert_parity(row: &str, direct: Part, mut recorded: Model) {
    let r = recorded
        .to_part()
        .unwrap_or_else(|e| panic!("{row}: recorder part: {e}"));
    let (a, b) = (signature(&direct), signature(&r));
    assert_eq!(
        (a.0, a.1, a.2),
        (b.0, b.1, b.2),
        "{row}: census (faces, edges, vertices)"
    );
    assert_eq!(a.3, b.3, "{row}: volume");
    assert_eq!(a.4, b.4, "{row}: mesh vertices");
    assert_eq!(
        recorded.report().failures().count(),
        0,
        "{row}: replay failures"
    );
}

fn cube_direct() -> Part {
    Profile::on(PartPlane::XY)
        .rect(1.0, 1.0)
        .extrude(1.0)
        .expect("cube")
}

fn cube_recorded(m: &mut Model) -> Body {
    m.profile(Plane::XY)
        .rect(1.0, 1.0)
        .sketch("cube")
        .expect("sketch")
        .extrude(m, End::Blind(1.0.into()))
        .expect("extrude")
}

#[test]
fn profiles_extrude_in_parity() {
    // rect
    let mut m = Model::new();
    cube_recorded(&mut m);
    assert_parity("rect extrude", cube_direct(), m);
    // circle, placed
    let direct = Profile::on(PartPlane::xy_at(0.3))
        .circle(0.35)
        .at(0.4, -0.2)
        .extrude(0.8)
        .expect("circle");
    let mut m = Model::new();
    m.profile(Plane::xy_at(0.3))
        .circle(0.35)
        .at(0.4, -0.2)
        .sketch("peg")
        .expect("sketch")
        .extrude(&mut m, End::Blind(0.8.into()))
        .expect("extrude");
    assert_parity("circle extrude", direct, m);
    // polygon (hexagon)
    let hex: Vec<(f64, f64)> = (0..6)
        .map(|i| {
            let a = i as f64 * std::f64::consts::TAU / 6.0;
            (a.cos(), a.sin())
        })
        .collect();
    let direct = Profile::on(PartPlane::XY)
        .polygon(&hex)
        .extrude(0.5)
        .expect("hex");
    let mut m = Model::new();
    m.profile(Plane::XY)
        .polygon(&hex)
        .sketch("hex")
        .expect("sketch")
        .extrude(&mut m, End::Blind(0.5.into()))
        .expect("extrude");
    assert_parity("polygon extrude", direct, m);
    // path: line/arc/line
    let direct = Profile::on(PartPlane::XY)
        .path(0.0, 0.0)
        .line_to(2.0, 0.0)
        .arc_to((2.5, 0.5), (2.0, 1.0))
        .line_to(0.0, 1.0)
        .close()
        .extrude(0.4)
        .expect("path");
    let mut m = Model::new();
    m.profile(Plane::XY)
        .path(0.0, 0.0)
        .line_to(2.0, 0.0)
        .arc_to((2.5, 0.5), (2.0, 1.0))
        .line_to(0.0, 1.0)
        .sketch("outline")
        .expect("sketch")
        .extrude(&mut m, End::Blind(0.4.into()))
        .expect("extrude");
    assert_parity("path extrude", direct, m);
}

#[test]
fn dressups_in_parity() {
    // fillet, vertical edges
    let direct = cube_direct()
        .fillet(0.1, |e| e.straight().direction(Dir::Z))
        .expect("fillet");
    let mut m = Model::new();
    let cube = cube_recorded(&mut m);
    cube.fillet(&mut m, edges().straight().direction(Dir::Z), 0.1)
        .expect("fillet");
    assert_parity("vertical fillet", direct, m);
    // chamfer, top rim
    let direct = cube_direct()
        .chamfer(0.1, |e| e.straight().at_height(1.0))
        .expect("chamfer");
    let mut m = Model::new();
    let cube = cube_recorded(&mut m);
    cube.chamfer(&mut m, edges().straight().at_height(1.0), 0.1)
        .expect("chamfer");
    assert_parity("rim chamfer", direct, m);
    // shell, open top
    let direct = cube_direct()
        .shell(0.1, |f| f.facing(Dir::Z))
        .expect("shell");
    let mut m = Model::new();
    let cube = cube_recorded(&mut m);
    cube.shell(
        &mut m,
        faces().facing(Dir::Z),
        0.1,
        ShellDirectionKind::Inward,
    )
    .expect("shell");
    assert_parity("open-top shell", direct, m);
}

#[test]
fn booleans_and_transforms_in_parity() {
    // union, half stagger
    let direct = cube_direct()
        .union(&cube_direct().translate(0.5, 0.0, 0.0).expect("t"))
        .expect("union");
    let mut m = Model::new();
    let a = cube_recorded(&mut m);
    let b = m
        .tool(|m| {
            let b = cube_recorded(m);
            b.translate(m, 0.5, 0.0, 0.0)
        })
        .expect("tool");
    a.union(&mut m, b).expect("union");
    assert_parity("stagger union", direct, m);
    // subtract, bore
    let direct = cube_direct()
        .subtract(
            &Profile::on(PartPlane::xy_at(-0.5))
                .circle(0.2)
                .extrude(2.0)
                .expect("peg"),
        )
        .expect("subtract");
    let mut m = Model::new();
    let a = cube_recorded(&mut m);
    let peg = m
        .tool(|m| {
            m.profile(Plane::xy_at(-0.5))
                .circle(0.2)
                .sketch("peg")?
                .extrude(m, End::Blind(2.0.into()))
        })
        .expect("tool");
    a.subtract(&mut m, peg).expect("subtract");
    assert_parity("bore subtract", direct, m);
    // rotate about z by 30°
    let direct = cube_direct()
        .rotate_z(30.0_f64.to_radians())
        .expect("rotate");
    let mut m = Model::new();
    let cube = cube_recorded(&mut m);
    cube.rotate(
        &mut m,
        AxisDesc::Literal {
            origin: [0.0, 0.0, 0.0],
            direction: [0.0, 0.0, 1.0],
        },
        30.0,
    )
    .expect("rotate");
    assert_parity("rotate z", direct, m);
}

#[test]
fn labels_bind_selections_and_renames_change_nothing() {
    let mut m = Model::new();
    let cube = cube_recorded(&mut m);
    // Label the top face and one vertical edge, then select through the labels.
    let top = m
        .label_face("lid", &faces().facing(Dir::Z))
        .expect("label face");
    let edge = m
        .label_edge(
            "post",
            &edges()
                .straight()
                .direction(Dir::Z)
                .chord_length_between(0.0, 2.0),
        )
        .expect_err("four vertical edges: cardinality");
    assert!(matches!(
        edge,
        ScriptError::SelectionCardinality { got: 4, .. }
    ));
    let post = m
        .label_edge(
            "post",
            &edges()
                .named(top.clone())
                .allow_empty()
                .or(edges().straight().direction(Dir::Z)),
        )
        .err();
    assert!(post.is_some(), "still four edges");
    // Fillet through the label of the lid's first vertical... use a face label in a query.
    let filleted = cube
        .fillet(&mut m, edges().straight().direction(Dir::Z), 0.1)
        .expect("fillet");
    let before = m.to_part().expect("part");
    let sig = signature(&before);
    // Rename the feature and the entity label: name snapshots keep replaying.
    let rounds = m.feature("fillet-1").expect("label");
    m.rename(&rounds, "rounds").expect("rename");
    m.rename_label("lid", "top").expect("rename label");
    assert!(m.feature("fillet-1").is_err());
    assert!(m.feature("rounds").is_ok());
    assert_eq!(m.labels().len(), 1);
    let report = m.replay().expect("replay");
    assert_eq!(report.failures().count(), 0);
    assert_eq!(
        signature(&m.to_part().expect("part")),
        sig,
        "rename changed nothing"
    );
    // Duplicate label is strict (a different entity under a taken name).
    assert!(matches!(
        m.label_face("top", &faces().facing_neg(Dir::Z)),
        Err(ScriptError::History(_))
    ));
    let _ = filleted;
    // A label-bound face selection resolves after replay.
    let by_label = faces().labeled("top");
    let mut again = Model::new();
    again
        .profile(Plane::XY)
        .rect(1.0, 1.0)
        .sketch("c")
        .expect("s")
        .extrude(&mut again, End::Blind(1.0.into()))
        .expect("e");
    again
        .label_face("top", &faces().facing(Dir::Z))
        .expect("label");
    let shelled = again
        .active()
        .expect("active")
        .shell(&mut again, by_label, 0.1, ShellDirectionKind::Inward)
        .expect("shell through label");
    let expected = cube_direct()
        .shell(0.1, |f| f.facing(Dir::Z))
        .expect("direct shell")
        .face_count();
    assert_eq!(again.to_part().expect("part").face_count(), expected);
    let _ = shelled;
}
