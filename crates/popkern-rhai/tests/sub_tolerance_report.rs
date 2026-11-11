//! A boolean whose tool is thinner than the linear tolerance resolves
//! it as a coincidence and says so in the model's report (PK-1131): a
//! cube cut by a saw 1e-6 wide serves whole, the boolean's row warning
//! that the tool was under the tolerance and contributed nothing.

use popkern_rhai::ScriptHost;

#[test]
fn a_saw_under_the_tolerance_cuts_nothing_and_says_so() {
    let source = std::fs::read_to_string(
        std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/tangency_bench/hairline_slot.rhai"),
    )
    .expect("script");
    let mut host = ScriptHost::new();
    let mut outcome = host
        .run(&source, "hairline_slot.rhai")
        .unwrap_or_else(|e| panic!("{}", e.render()));
    let model = outcome.model.as_mut().expect("model");
    let warned: Vec<String> = model
        .report()
        .warnings()
        .map(|(_, w)| w.message.clone())
        .filter(|m| m.contains("under the linear tolerance"))
        .collect();
    assert_eq!(warned.len(), 1, "{warned:?}");
    let part = model.to_part().expect("part");
    let v = popkern_measure::body_volume(part.brep(), &part.brep().tolerance).expect("volume");
    assert!((v - 1000.0).abs() < 1e-9 * 1000.0, "{v}");
}

/// Three blocks a 5e-7 gap apart in a row, fused by one N-ary union:
/// the model's report warns of each gap the union closed (PK-1162), as
/// a two-operand union's row does.
#[test]
fn an_nary_union_closing_gaps_under_the_tolerance_says_so() {
    let source = r#"let m = model();
let a = m.profile(xy()).rect(10.0, 10.0).at(5.0, 5.0).sketch("a").extrude(blind(10.0));
let b = m.tool(|m| m.profile(xy()).rect(10.0, 10.0).at(15.0000005, 5.0).sketch("b").extrude(blind(10.0)));
let c = m.tool(|m| m.profile(xy()).rect(10.0, 10.0).at(25.000001, 5.0).sketch("c").extrude(blind(10.0)));
a.nary_boolean([b, c], any_body());
m"#;
    let mut host = ScriptHost::new();
    let mut outcome = host
        .run(source, "three_blocks.rhai")
        .unwrap_or_else(|e| panic!("{}", e.render()));
    let model = outcome.model.as_mut().expect("model");
    let warned: Vec<String> = model
        .report()
        .warnings()
        .map(|(_, w)| w.message.clone())
        .filter(|m| m.contains("the union closed the gap"))
        .collect();
    assert_eq!(warned.len(), 2, "{warned:?}");
    let part = model.to_part().expect("part");
    assert_eq!(part.brep().faces.len(), 6);
}
