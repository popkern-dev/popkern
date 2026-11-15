//! PK-1191: the boss plate turned about z before its rounds serves the
//! unturned body. The corner sphere's centre was built by stepping the
//! ball's radius along the global axes from the plate corner, gated on
//! axis-aligned fillets: turned 23 or 137 degrees every corner kept its
//! cube-corner outlier and the rounds declined, 36 edges 3.3e-2 off
//! their faces. The centre now steps along the fillets' own axes
//! (`surface_fix_sphere::corner_frame`).

use popkern_core::tolerance::Tolerance;

fn rounds(bx: f64, by: f64, turn: f64) -> String {
    format!(
        r#"let m = model();
let plate = m.profile(xy()).rect(2.4, 1.6).sketch("plate_sketch").extrude(blind(0.4));
let boss = m.tool(|m| m.profile(xy_at(0.2)).circle(0.35).sketch("boss_sketch")
             .extrude(blind(0.8))).translate({bx}, {by}, 0);
plate = plate.union(boss);
plate = plate.rotate(axis_at([0.0, 0.0, 0.0], [0.0, 0.0, 1.0]), {turn});
plate.fillet(
    edges().between("cylinder", "plane").at_height(0.4)
        .or(edges().straight().at_height(0.4))
        .or(edges().straight().direction("z")),
    0.08, #{{ continuity: "g1", label: "rounds" }});
m.assert_watertight();
m"#
    )
}

/// Faces and volume of the served body, validated (self-intersection
/// included) and watertight at the display chord.
fn serve(source: &str) -> (usize, f64) {
    let mut host = popkern_rhai::ScriptHost::with_options(Default::default());
    let run = host
        .run(source, "turned_boss_plate")
        .unwrap_or_else(|e| panic!("{e}"));
    let b = run.model.expect("model").to_parts_all().expect("parts")[0]
        .brep()
        .clone();
    let errors = b.validate_all();
    assert!(errors.is_empty(), "validate_all: {errors:?}");
    let mesh = popkern_mesh::tessellate(&b, &popkern_mesh::TessellationOptions::display()).unwrap();
    assert_eq!(
        popkern_mesh::welded_open_edges(&mesh),
        0,
        "display mesh open"
    );
    (
        b.faces.len(),
        popkern_measure::body_volume(&b, &Tolerance::default()).unwrap(),
    )
}

#[test]
fn every_turn_serves_the_unturned_boss_plate() {
    // The unturned plate reads OCC's volume, 1.757239352.
    for (bx, by) in [(0.0, 0.0), (0.6, 0.1)] {
        let (faces, unturned) = serve(&rounds(bx, by, 0.0));
        assert!((unturned - 1.757239352).abs() < 1e-8, "{unturned}");
        for turn in [23.0, 45.0, 137.0, -17.5] {
            let (f, v) = serve(&rounds(bx, by, turn));
            assert_eq!(f, faces, "boss ({bx},{by}) turned {turn}: faces");
            assert!(
                (v - unturned).abs() < 2e-9,
                "boss ({bx},{by}) turned {turn}: {v} vs {unturned}"
            );
        }
    }
}
