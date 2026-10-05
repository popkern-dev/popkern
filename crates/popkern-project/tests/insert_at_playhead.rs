//! PK-666 "insert a feature here": the insert-after migration template
//! at a history playhead builds — the new feature lands right after the
//! anchor, the later features replay after it, and the result round-trips
//! through `.pkn`; an anchor by position matches an anchor by label; the
//! unedited template declines saying what to do.

mod support;

use std::rc::Rc;

use popkern_assembly::Document;
use popkern_project::scaffold::{insert_after_name, insert_after_part_source, insert_after_source};
use popkern_project::{MemoryFs, ProjectError, ProjectFs};
use support::{build, feature_count, memory, runner, statuses, volume};

/// Labels of the part document's tree, in timeline order.
fn labels(document: &mut Document) -> Vec<String> {
    let Document::Part(model) = document else {
        panic!("part document expected");
    };
    model.features().into_iter().map(|f| f.label).collect()
}

/// The bore subtraction of `labelled_rim` (the feature before the
/// rim fillet): its label and 1-based position.
fn anchor(fs: &Rc<MemoryFs>) -> (String, usize) {
    let mut built = build(&mut runner(fs), fs);
    let Some(Document::Part(model)) = built.document.as_mut() else {
        panic!("part document expected");
    };
    let features = model.features();
    let position = features
        .iter()
        .position(|f| f.operation == "boolean")
        .expect("the bore subtraction");
    (features[position].label.clone(), position + 1)
}

/// The template with its closure filled: a chamfer on the block's
/// four straight top edges (not the bore rim the later fillet rounds).
fn filled(source: &str) -> String {
    let filled = source.replace(
        "    body\n",
        "    body.chamfer(edges().straight().at_height(1.0), 0.05)\n",
    );
    assert_ne!(filled, source, "the template has a closure body to fill");
    filled
}

fn insert(fs: &MemoryFs, label: &str, source: &str) {
    let path = format!("migrations/003_{}.rhai", insert_after_name(label));
    fs.write(&path, source.as_bytes()).expect("write");
}

#[test]
fn an_insert_at_the_playhead_builds_after_its_anchor_and_round_trips() {
    let fs_rc = memory("labelled_rim");
    let (label, position) = anchor(&fs_rc);
    let mut runner = runner(&fs_rc);
    let mut before = build(&mut runner, &fs_rc);
    let before_labels = labels(before.document.as_mut().expect("doc"));
    let before_volume = volume(before.document.as_mut().expect("doc"));

    insert(
        &fs_rc,
        &label,
        &filled(&insert_after_source(&label, position, true)),
    );
    let mut after = build(&mut runner, &fs_rc);
    assert!(after.error.is_none(), "{:?}", after.error);
    assert_eq!(statuses(&after), ["cached", "cached", "built"]);
    let document = after.document.as_mut().expect("doc");
    let after_labels = labels(document);
    assert_eq!(feature_count(document), before_labels.len() + 1);
    assert_eq!(after_labels[position - 1], label, "the anchor stays put");
    assert!(
        !before_labels.contains(&after_labels[position]),
        "the new feature sits right after the anchor: {after_labels:?}"
    );
    assert_eq!(
        after_labels[position + 1..],
        before_labels[position..],
        "the later features replay after it"
    );
    let chamfered = volume(document);
    assert!(
        chamfered < before_volume - 1e-6,
        "{chamfered} vs {before_volume}"
    );

    // Saved and opened, the tree and the part are the same.
    let bytes = document.to_pkn_bytes().expect("save");
    let mut opened = Document::from_pkn_bytes(&bytes).expect("open");
    assert_eq!(labels(&mut opened), after_labels);
    assert!((volume(&mut opened) - chamfered).abs() <= 1e-9 * chamfered);
}

#[test]
fn an_anchor_by_position_inserts_where_the_label_would() {
    let by_label = memory("labelled_rim");
    let (label, position) = anchor(&by_label);
    let by_position = memory("labelled_rim");
    insert(
        &by_label,
        &label,
        &filled(&insert_after_source(&label, position, true)),
    );
    insert(
        &by_position,
        &label,
        &filled(&insert_after_source(&label, position, false)),
    );
    let mut a = build(&mut runner(&by_label), &by_label);
    let mut b = build(&mut runner(&by_position), &by_position);
    assert!(a.error.is_none() && b.error.is_none());
    assert_eq!(
        labels(a.document.as_mut().expect("doc")),
        labels(b.document.as_mut().expect("doc"))
    );
}

#[test]
fn the_unedited_template_declines_naming_the_closure_rule() {
    let fs = memory("labelled_rim");
    let (label, position) = anchor(&fs);
    insert(&fs, &label, &insert_after_source(&label, position, true));
    let built = build(&mut runner(&fs), &fs);
    match &built.error {
        Some(ProjectError::MigrationFailed { migration, .. }) => {
            assert!(migration.contains("insert_after"), "{migration}");
            let detail = built
                .error
                .as_ref()
                .map(ToString::to_string)
                .unwrap_or_default();
            assert!(detail.contains("exactly one feature"), "{detail}");
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn labels_are_quoted_and_names_are_file_safe() {
    let source = insert_after_source("rim \"round\" \\ 2", 4, true);
    assert!(
        source.contains(r#"m.feature("rim \"round\" \\ 2")"#),
        "{source}"
    );
    assert!(insert_after_source("x", 4, false).contains("m.features()[3]"));
    assert_eq!(insert_after_name("Rim round!"), "insert_after_rim_round");
    assert_eq!(insert_after_name("***"), "insert_feature");
}

#[test]
fn an_assembly_part_takes_an_insert_from_a_migration() {
    let files = vec![
        (
            "popkern.toml".to_string(),
            b"[project]\nname = \"Pins\"\n".to_vec(),
        ),
        (
            "migrations/001_pins.rhai".to_string(),
            b"let a = assembly(\"pins\");\nlet pin = a.part(\"pin\", |m| { m.profile(xy()).circle(0.1).sketch(\"pin\").extrude(blind(1.0)); });\na.place(pin, \"left\");\na.place(pin, \"right\").translate(0.5, 0.0, 0.0);\na\n".to_vec(),
        ),
    ];
    let fs = Rc::new(MemoryFs::new(files).expect("files"));
    let mut runner = runner(&fs);
    assert!(build(&mut runner, &fs).error.is_none());
    let source = insert_after_part_source("pin", "extrude-1", 2, true).replace(
        "        body\n",
        "        body.chamfer(edges().circular().at_height(1.0), 0.02)\n",
    );
    fs.write(
        &format!("migrations/010_{}.rhai", insert_after_name("pin extrude-1")),
        source.as_bytes(),
    )
    .expect("write");
    let mut built = build(&mut runner, &fs);
    assert!(built.error.is_none(), "{:?}", built.error);
    assert_eq!(statuses(&built), ["cached", "built"]);
    let Some(Document::Assembly(assembly)) = built.document.as_mut() else {
        panic!("assembly document expected");
    };
    let id = assembly.definition_named("pin").expect("pin");
    let labels: Vec<String> = assembly
        .part_model_mut(id)
        .expect("part")
        .features()
        .into_iter()
        .map(|f| f.label)
        .collect();
    assert_eq!(labels, ["pin", "extrude-1", "chamfer-1"]);
    assert_eq!(assembly.leaves().len(), 2, "both instances still placed");
}

/// PK-677: a chamfer inserted before the bore subtraction renumbers the
/// subtraction's names; `rim_round`'s reference rebinds by provenance
/// to the same rim edge, and the part equals the one recorded with the
/// chamfer there from the start.
#[test]
fn an_insert_upstream_of_the_boolean_keeps_the_rim_reference() {
    let fs = memory("labelled_rim");
    let mut runner = runner(&fs);
    assert!(build(&mut runner, &fs).error.is_none());
    let chamfer = "    body.chamfer(edges().at_height(1.0), 0.05)\n";
    let source = insert_after_source("extrude-1", 2, true).replace("    body\n", chamfer);
    insert(&fs, "extrude-1", &source);
    let mut inserted = build(&mut runner, &fs);
    assert!(inserted.error.is_none(), "{:?}", inserted.error);
    let document = inserted.document.as_mut().expect("doc");
    assert_eq!(
        labels(document),
        [
            "block_sketch",
            "extrude-1",
            "chamfer-1",
            "bore_sketch",
            "extrude-2",
            "boolean-1",
            "rim_round"
        ]
    );
    let Document::Part(model) = &*document else {
        panic!("part document expected");
    };
    let report = model.tree().last_report();
    let rim = report
        .features
        .iter()
        .find(|r| r.label == "rim_round")
        .expect("rim_round row");
    assert!(
        rim.fallbacks
            .iter()
            .any(|f| f.mechanism.contains("provenance rebind")),
        "the rebind is recorded on the row: {:?}",
        rim.fallbacks
    );
    let rebound_volume = volume(document);

    // The same part recorded with the chamfer in place from the start.
    let fresh = popkern_rhai::ScriptHost::new()
        .run(
            "let m = model();\nlet block = m.profile(xy()).rect(2.0, 2.0).sketch(\"block_sketch\").extrude(blind(1.0));\nblock.chamfer(edges().at_height(1.0), 0.05);\nlet bore = m.tool(|m| m.profile(xy()).circle(0.4).sketch(\"bore_sketch\").extrude(blind(1.0)));\nm.active().subtract(bore);\nm.active().fillet(edges().circular().at_height(1.0), 0.05);\nm\n",
            "fresh.rhai",
        )
        .map(|mut out| out.take_document())
        .expect("fresh model runs");
    let mut fresh = fresh.expect("a document");
    let fresh_volume = volume(&mut fresh);
    assert!(
        (rebound_volume - fresh_volume).abs() <= 1e-9 * fresh_volume,
        "the rebound fillet is the rim fillet: {rebound_volume} vs {fresh_volume}"
    );
}
