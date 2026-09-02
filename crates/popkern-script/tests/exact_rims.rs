//! PK-417: a bore cut through a block leaves exact circular rims — the
//! boolean re-spells the plane × cylinder intersection curves as the
//! `Arc`s they are — so the rims select as `circular()`, carry the
//! bore's radius, take a fillet by label, and the body stays
//! watertight and STEP-clean. Union rims (a boss meeting a plate) are
//! the same class.

use std::path::Path;

use popkern_core::geometry::curve::CurveGeometry;
use popkern_script::prelude::*;
use popkern_script::{Model, Plane};

fn bored_block(m: &mut Model) -> Body {
    let block = m
        .profile(Plane::XY)
        .rect(1.2, 0.8)
        .sketch("block_sketch")
        .expect("sketch")
        .extrude(m, End::Blind(0.5.into()))
        .expect("extrude");
    let bore = m
        .tool(|m| {
            m.profile(Plane::xy_at(-0.1))
                .circle(0.18)
                .sketch("bore_sketch")?
                .extrude(m, End::Blind(0.7.into()))
        })
        .expect("bore");
    block.subtract(m, bore).expect("subtract")
}

#[test]
fn a_bore_leaves_two_exact_circular_rims_that_select_and_fillet_by_label() {
    let mut m = Model::new();
    let mut block = bored_block(&mut m);
    let brep = m.brep_of(block).expect("brep");
    let arcs: Vec<f64> = brep
        .edges
        .values()
        .filter_map(|e| match e.curve.as_ref().map(|c| &c.geometry) {
            Some(CurveGeometry::Arc { radius, .. })
            | Some(CurveGeometry::Circle { radius, .. }) => Some(*radius),
            _ => None,
        })
        .collect();
    assert_eq!(arcs.len(), 2, "two rims, both exact: {arcs:?}");
    assert!(arcs.iter().all(|r| (r - 0.18).abs() < 1e-9), "{arcs:?}");
    let nurbs = brep
        .edges
        .values()
        .filter(|e| {
            matches!(
                e.curve.as_ref().map(|c| &c.geometry),
                Some(CurveGeometry::Nurbs(_))
            )
        })
        .count();
    assert_eq!(nurbs, 0, "no NURBS edge is left on a block with a bore");
    let issues = brep.validate_all();
    assert!(issues.is_empty(), "valid after promotion: {issues:?}");
    // The rims select by kind and by height, and take a fillet by label.
    m.label_edge("bore_rim", &edges().circular().at_height(0.5))
        .expect("the top rim is circular");
    block = block
        .fillet(&mut m, edges().labeled("bore_rim"), 0.05)
        .expect("fillet on the labeled rim");
    let _ = block;
    assert_eq!(m.report().failures().count(), 0);
    m.assert_watertight().expect("watertight");
    let step = Path::new(env!("CARGO_TARGET_TMPDIR")).join("bored_block.step");
    m.assert_step_roundtrip(&step).expect("STEP round trip");
    let oracle = Path::new(&std::env::var("HOME").unwrap_or_default()).join("venvs/ocp/bin/python");
    let check = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tools/occ-oracle/occ_check.py");
    if oracle.exists() && check.exists() {
        let output = std::process::Command::new(&oracle)
            .arg(&check)
            .arg(&step)
            .output()
            .expect("oracle");
        let report = String::from_utf8_lossy(&output.stdout);
        assert!(report.contains("\"passed\": true"), "OCC oracle: {report}");
    }
}

/// Containment pin (PK-419): a union's junction rim stays NURBS until
/// the G2 closed-fillet ladder takes an arc there — promoting it left
/// the flagship boss plate's nine-edge G2 fillet with 70 open mesh
/// edges. The boss's own top rim is the extrude's exact arc as ever.
#[test]
fn a_union_junction_rim_stays_nurbs_until_pk419() {
    let mut m = Model::new();
    let plate = m
        .profile(Plane::XY)
        .rect(2.4, 1.6)
        .sketch("plate_sketch")
        .expect("sketch")
        .extrude(&mut m, End::Blind(0.4.into()))
        .expect("extrude");
    let boss = m
        .tool(|m| {
            m.profile(Plane::xy_at(0.2))
                .circle(0.35)
                .sketch("boss_sketch")?
                .extrude(m, End::Blind(0.8.into()))
        })
        .expect("boss");
    let plate = plate.union(&mut m, boss).expect("union");
    let brep = m.brep_of(plate).expect("brep");
    let circular = brep
        .edges
        .values()
        .filter(|e| {
            matches!(
                e.curve.as_ref().map(|c| &c.geometry),
                Some(CurveGeometry::Arc { .. }) | Some(CurveGeometry::Circle { .. })
            )
        })
        .count();
    assert_eq!(circular, 1, "only the boss top rim is an arc (PK-419)");
    let junction = m.label_edge("junction", &edges().circular().at_height(0.4));
    assert!(
        junction.is_err(),
        "the junction rim is not circular yet: {junction:?}"
    );
    m.label_edge("junction", &edges().nurbs().at_height(0.4))
        .expect("the junction rim selects as NURBS");
}
