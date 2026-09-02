//! PK-410: a saved document is byte-deterministic — two trees built
//! the same way save to the same bytes, in one process and across
//! processes. Found when the browser build's `.pkn` differed from the
//! command's for the boss plate: the parameter name index (and three
//! other maps) serialised in `HashMap` order, which differs per map
//! instance, so even two native saves of one script disagreed.

use popkern_core::tolerance::Tolerance;
use popkern_history::FeatureTree;
use popkern_history::unit::Unit;

fn tree() -> FeatureTree {
    let mut tree = FeatureTree::new();
    tree.set_build_tolerance(Some(Tolerance::default()));
    // A dozen globals: enough keys that hash order almost surely
    // differs between two independent maps.
    for (i, name) in [
        "wall", "height", "radius", "depth", "pitch", "count", "gap", "lip", "rim", "boss", "web",
        "rib",
    ]
    .iter()
    .enumerate()
    {
        tree.add_global_parameter_display(name, i as f64 + 1.0, Unit::Meters)
            .expect("global");
    }
    tree
}

#[test]
fn two_trees_built_the_same_way_save_to_identical_bytes() {
    let a = tree().save().expect("save a");
    let b = tree().save().expect("save b");
    assert_eq!(a.len(), b.len());
    assert!(a == b, "documents differ (map order leaked into the bytes)");
    // The name index is sorted in the bytes.
    let text = String::from_utf8(a).expect("utf8");
    let start = text.find("\"name_index\"").expect("name_index present");
    let keys: Vec<&str> = text[start..]
        .split('"')
        .filter(|k| ["boss", "count", "depth", "gap", "wall", "web"].contains(k))
        .collect();
    let mut sorted = keys.clone();
    sorted.sort_unstable();
    assert_eq!(keys, sorted, "name index keys are sorted");
}
