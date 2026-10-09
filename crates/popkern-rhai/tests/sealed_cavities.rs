//! Unions that seal a cavity (PK-745): the result is one solid with an
//! outer shell and a cavity shell. A cavity bounded by a cylinder and
//! two discs has only seam vertices, so the shell partition must box
//! components by their faces, not their vertices — it used to find no
//! interior point, keep one shell, and decline as an Euler violation.

use popkern_rhai::ScriptHost;

fn served(source: &str) -> popkern_core::topology::BRep {
    let mut host = ScriptHost::new();
    let outcome = host
        .run(source, "cavity.rhai")
        .unwrap_or_else(|e| panic!("{}", e.render()));
    let mut model = outcome.model.expect("model");
    model.to_part().expect("part").brep().clone()
}

#[test]
fn a_plug_sealing_a_cylindrical_cup_leaves_a_cavity_shell() {
    let brep = served(
        r#"let m = model();
let cup = m.profile(xy()).circle(10.0).sketch("c").extrude(blind(20.0));
cup = cup.subtract(m.tool(|t| t.profile(xy_at(2.0)).circle(8.0).sketch("h").extrude(blind(20.0))));
cup.union(m.tool(|t| t.profile(xy_at(8.0)).circle(9.0).sketch("plug").extrude(blind(4.0))), #{ label: "plug" });
m"#,
    );
    assert_eq!(brep.solids.len(), 1);
    assert_eq!(brep.shells.len(), 2, "outer shell + cavity shell");
    assert!(brep.validate().is_ok());
    // Cup (r 10 × 20 less the r 8 pocket above z 2) plus the plug's
    // r 8 × 4 core: π(2000 − 1152 + 256) = 1104π.
    let volume = popkern_measure::body_volume(&brep, &brep.tolerance).expect("volume");
    let exact = 1104.0 * std::f64::consts::PI;
    assert!(
        ((volume - exact) / exact).abs() < 1e-3,
        "volume {volume} vs {exact}"
    );
}

#[test]
fn a_square_lid_on_a_cylindrical_cup_leaves_a_cavity_shell() {
    let brep = served(
        r#"let m = model();
let cup = m.profile(xy()).circle(10.0).sketch("c").extrude(blind(20.0));
cup = cup.subtract(m.tool(|t| t.profile(xy_at(2.0)).circle(8.0).sketch("h").extrude(blind(20.0))));
cup.union(m.tool(|t| t.profile(xy_at(8.0)).rect(18.0, 18.0).sketch("lid").extrude(blind(4.0))), #{ label: "lid" });
m"#,
    );
    assert_eq!(brep.solids.len(), 1);
    assert_eq!(brep.shells.len(), 2, "outer shell + cavity shell");
    assert!(brep.validate().is_ok());
}
