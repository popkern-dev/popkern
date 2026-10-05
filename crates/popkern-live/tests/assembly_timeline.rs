//! PK-666: an assembly scrubs per construction step — each place, move
//! and mate is a timeline row with its script line and migration, and
//! the occurrences after any step are the instances placed by then at
//! their placements then; the last step is the run's assembly.

use std::path::Path;

use popkern_live::envelope::InstanceView;
use popkern_live::{HostOptions, Session};

const TWO_BLOCKS: &str = r#"let a = assembly("two");
let block = a.part("block", |m| { m.profile(xy()).rect(1.0, 0.5).sketch("block").extrude(blind(0.2)); });
let base = a.place(block, "base");
let moving = a.place(block, "moving").translate(3.0, -1.0, 0.7);
a.ground(base);
a.coincident(moving.datum("xy"), base.datum("xy"), #{ aligned: true });
a
"#;

fn tracked() -> Session {
    Session::new(HostOptions {
        track_sources: true,
        ..HostOptions::default()
    })
}

fn z_of(row: &InstanceView) -> f32 {
    row.world[14]
}

fn same(a: &[InstanceView], b: &[InstanceView]) -> bool {
    a.len() == b.len()
        && a.iter()
            .zip(b)
            .all(|(x, y)| x.path == y.path && x.world == y.world && x.definition == y.definition)
}

#[test]
fn a_mate_is_a_step_that_moves_what_it_mates() {
    let mut session = tracked();
    let envelope = session.run(TWO_BLOCKS, "two.rhai");
    assert!(
        matches!(envelope.outcome, popkern_rhai::Outcome::Served),
        "{:?}",
        envelope.outcome
    );
    let timeline = session
        .assembly_timeline()
        .expect("timeline")
        .expect("an assembly");
    assert!(timeline.recorded);
    let kinds: Vec<&str> = timeline.steps.iter().map(|s| s.kind).collect();
    assert_eq!(kinds, ["place", "place", "mate"]);
    assert_eq!(timeline.steps[1].label, "place moving (block)");
    assert_eq!(timeline.steps[2].instances, ["moving"]);
    let lines: Vec<Option<u32>> = timeline
        .steps
        .iter()
        .map(|s| s.source.as_ref().map(|s| s.line))
        .collect();
    assert_eq!(lines, [Some(3), Some(4), Some(6)], "the authoring lines");
    assert!(timeline.steps.iter().all(|s| s.migration.is_none()));

    assert!(session.assembly_step(0).expect("step 0").is_empty());
    assert_eq!(session.assembly_step(1).expect("step 1").len(), 1);
    let before = session.assembly_step(2).expect("step 2");
    assert!((z_of(&before[1]) - 0.7).abs() < 1e-6, "authored height");
    let after = session.assembly_step(3).expect("step 3");
    assert!(z_of(&after[1]).abs() < 1e-6, "mated flush with base");
    let run = envelope.assembly.expect("assembly view").instances;
    assert!(same(&after, &run), "the last step is the run's assembly");
}

#[test]
fn a_project_attributes_each_step_to_its_migration() {
    let mut session = tracked();
    let first: Vec<&str> = TWO_BLOCKS.lines().take(4).collect();
    let files = vec![
        (
            "popkern.toml".to_string(),
            b"[project]\nname = \"Two\"\n".to_vec(),
        ),
        (
            "migrations/001_base.rhai".to_string(),
            format!("{}\na\n", first.join("\n")).into_bytes(),
        ),
        (
            "migrations/002_mate.rhai".to_string(),
            b"let a = previous();\na.coincident(a.occurrence(\"moving\").datum(\"xy\"), a.occurrence(\"base\").datum(\"xy\"), #{ aligned: true });\na\n".to_vec(),
        ),
    ];
    let envelope = session.run_project(files);
    assert!(
        matches!(envelope.outcome, popkern_rhai::Outcome::Served),
        "{:?}",
        envelope.project
    );
    let timeline = session
        .assembly_timeline()
        .expect("timeline")
        .expect("an assembly");
    let by: Vec<Option<&str>> = timeline
        .steps
        .iter()
        .map(|s| s.migration.as_deref())
        .collect();
    assert_eq!(by, [Some("001_base"), Some("001_base"), Some("002_mate")]);
}

#[test]
fn a_part_document_has_no_assembly_steps() {
    let mut session = tracked();
    session.run(
        "let m = model(); m.profile(xy()).rect(1.0, 1.0).sketch(\"s\").extrude(blind(0.1)); m",
        "part.rhai",
    );
    assert!(session.assembly_timeline().expect("timeline").is_none());
    assert!(session.assembly_step(1).is_err());
}

#[test]
fn the_gearbox_gains_one_occurrence_per_place_and_ends_at_the_run() {
    let source = std::fs::read_to_string(
        Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../../web/live/examples/planetary_gearbox.rhai"),
    )
    .expect("gallery script");
    let mut session = tracked();
    let envelope = session.run(&source, "planetary_gearbox.rhai");
    assert!(
        matches!(envelope.outcome, popkern_rhai::Outcome::Served),
        "{:?}",
        envelope.outcome
    );
    let timeline = session
        .assembly_timeline()
        .expect("timeline")
        .expect("an assembly");
    let mut occurrences = 0;
    for step in &timeline.steps {
        let rows = session.assembly_step(step.index).expect("step");
        if step.kind == "place" {
            assert!(rows.len() > occurrences, "{} adds occurrences", step.label);
        }
        occurrences = rows.len();
    }
    let last = session
        .assembly_step(timeline.steps.len())
        .expect("last step");
    let run = envelope.assembly.expect("assembly view").instances;
    assert!(same(&last, &run), "the last step is the run's assembly");
}
