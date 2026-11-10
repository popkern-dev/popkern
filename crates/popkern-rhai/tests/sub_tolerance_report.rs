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
