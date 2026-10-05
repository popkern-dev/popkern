//! PK-660 fixture projects: building a folder of migrations, reusing
//! checkpoints, attributing features, failing loudly where an upstream
//! edit broke a downstream migration, and refusing paths outside the
//! project.

mod support;

use popkern_assembly::Document;
use popkern_project::{NativeFs, Project, ProjectError, ProjectFs, VersionStatus};
use popkern_rhai::ScriptHost;
use support::{build, feature_count, fixture, memory, runner, statuses, volume};

#[test]
fn a_project_builds_the_same_part_as_one_script() {
    let fs = memory("boss_plate");
    let mut runner = runner(&fs);
    let mut built = build(&mut runner, &fs);
    assert!(built.complete(), "{:?}", built.error);
    assert_eq!(statuses(&built), ["built", "built", "built"]);
    let names: Vec<&str> = built
        .record
        .migrations
        .iter()
        .map(|m| m.name.as_str())
        .collect();
    assert_eq!(names, ["001_plate", "002_boss", "003_rounds"]);
    let per_migration: Vec<usize> = built.versions.iter().map(|v| v.features.len()).collect();
    assert_eq!(
        per_migration,
        [2, 3, 1],
        "plate sketch+extrude, boss sketch+extrude+union, rounds"
    );

    // The same part as the single gallery script.
    let gallery = std::fs::read_to_string(
        std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../../web/live/examples/boss_plate.rhai"),
    )
    .expect("gallery script");
    let mut single = ScriptHost::new()
        .run(&gallery, "boss_plate.rhai")
        .expect("gallery");
    let mut single = Document::Part(Box::new(single.model.take().expect("model")));
    let mut document = built.document.take().expect("document");
    assert_eq!(feature_count(&mut document), feature_count(&mut single));
    let (a, b) = (volume(&mut document), volume(&mut single));
    assert!((a - b).abs() < 1e-9 * b.abs().max(1.0), "volume {a} vs {b}");
}

#[test]
fn an_unchanged_project_reuses_every_checkpoint() {
    let fs = memory("boss_plate");
    let mut runner = runner(&fs);
    build(&mut runner, &fs);
    let again = build(&mut runner, &fs);
    assert_eq!(statuses(&again), ["cached", "cached", "cached"]);
}

#[test]
fn editing_the_last_migration_reruns_only_it() {
    let fs = memory("boss_plate");
    let mut runner = runner(&fs);
    build(&mut runner, &fs);
    let path = "migrations/003_rounds.rhai";
    let source = String::from_utf8(fs.read(path).expect("read")).expect("utf8");
    fs.write(path, source.replace("0.08", "0.06").as_bytes())
        .expect("edit");
    let again = build(&mut runner, &fs);
    assert_eq!(statuses(&again), ["cached", "cached", "built"]);
}

#[test]
fn a_renamed_or_re_ended_migration_is_the_same_migration() {
    let fs = memory("boss_plate");
    let mut runner = runner(&fs);
    build(&mut runner, &fs);
    // Rename 003 and give 002 Windows line endings with a BOM.
    let rounds = fs.read("migrations/003_rounds.rhai").expect("read");
    let mut files = fs.files();
    files.remove("migrations/003_rounds.rhai");
    files.insert("migrations/003_round_edges.rhai".into(), rounds);
    let boss = String::from_utf8(files["migrations/002_boss.rhai"].clone()).expect("utf8");
    files.insert(
        "migrations/002_boss.rhai".into(),
        format!("\u{feff}{}", boss.replace('\n', "\r\n")).into_bytes(),
    );
    let renamed = std::rc::Rc::new(popkern_project::MemoryFs::new(files).expect("files"));
    let project = Project::load(renamed.as_ref()).expect("load");
    let again = runner.build(&project).expect("build");
    assert_eq!(statuses(&again), ["cached", "cached", "cached"]);
    assert_eq!(again.record.migrations[2].name, "003_round_edges");
}

#[test]
fn an_upstream_parameter_edit_is_attributed_as_an_edit() {
    let fs = memory("boss_plate");
    fs.write(
        "migrations/004_thicker.rhai",
        b"let m = previous();\nm.set_param(\"extrude-1.distance\", 0.5);\nm\n",
    )
    .expect("add");
    let mut runner = runner(&fs);
    let mut built = build(&mut runner, &fs);
    assert!(built.complete(), "{:?}", built.error);
    let thicker = &built.versions[3];
    assert!(
        thicker.features.is_empty(),
        "a parameter edit records no feature"
    );
    assert_eq!(
        thicker
            .edited
            .iter()
            .map(|f| f.label.as_str())
            .collect::<Vec<_>>(),
        ["extrude-1"],
        "the plate extrude is the edited feature"
    );
    // The edit reached the geometry (plate 0.5 thick instead of 0.4).
    let mut document = built.document.take().expect("document");
    let thick = volume(&mut document);
    let fs = memory("boss_plate");
    let mut thin = build(&mut support::runner(&fs), &fs)
        .document
        .expect("document");
    assert!(thick > volume(&mut thin) + 0.3, "thicker plate");
}

#[test]
fn an_upstream_insertion_records_only_the_inserted_feature() {
    let fs = memory("boss_plate");
    fs.write(
        "migrations/004_foot_break.rhai",
        b"let m = previous();\nm.insert_before(m.feature(\"rounds\"), |m, b| b.chamfer(edges().straight().at_height(0.0), 0.05, #{ label: \"foot_break\" }));\nm\n",
    )
    .expect("add");
    let mut runner = runner(&fs);
    let built = build(&mut runner, &fs);
    assert!(built.complete(), "{:?}", built.error);
    let inserted = &built.versions[3];
    assert_eq!(
        inserted
            .features
            .iter()
            .map(|f| f.label.as_str())
            .collect::<Vec<_>>(),
        ["foot_break"]
    );
    assert!(
        inserted.edited.is_empty(),
        "inserting before rounds does not edit rounds: {:?}",
        inserted.edited
    );
}

#[test]
fn an_upstream_edit_that_drops_a_label_fails_downstream_with_attribution() {
    let fs = memory("labelled_rim");
    let mut runner = runner(&fs);
    let built = build(&mut runner, &fs);
    assert!(built.complete(), "{:?}", built.error);

    let path = "migrations/001_block.rhai";
    let source = String::from_utf8(fs.read(path).expect("read")).expect("utf8");
    let without = source
        .lines()
        .filter(|l| !l.contains("label_edge"))
        .collect::<Vec<_>>()
        .join("\n");
    fs.write(path, without.as_bytes()).expect("edit");
    let broken = build(&mut runner, &fs);
    assert_eq!(statuses(&broken), ["built", "failed"]);
    match &broken.error {
        Some(ProjectError::MigrationFailed {
            migration,
            error,
            upstream_changed,
        }) => {
            assert_eq!(migration, "002_round.rhai");
            assert_eq!(error.line, Some(2), "located at the fillet call");
            assert_eq!(upstream_changed, &["001_block".to_string()]);
        }
        other => panic!("expected MigrationFailed, got {other:?}"),
    }
    assert!(broken.document.is_some(), "the last good document is kept");
}

#[test]
fn a_step_base_starts_the_document() {
    let fs = memory("step_base");
    let mut runner = runner(&fs);
    let built = build(&mut runner, &fs);
    assert!(built.complete(), "{:?}", built.error);
    let record = built.record.base.as_ref().expect("base recorded");
    assert_eq!(record.path, "imports/cube.step");
    match built.document.expect("document") {
        Document::Assembly(assembly) => {
            assert_eq!(assembly.definitions().len(), 2, "the cube and the pin");
        }
        Document::Part(_) => panic!("a STEP base is an assembly"),
    }
}

#[test]
fn paths_outside_the_project_are_refused() {
    let fs = NativeFs::new(fixture("outside_base")).expect("fixture");
    let error = Project::load(&fs).expect_err("escaping base");
    assert!(
        matches!(error, ProjectError::PathOutsideProject { .. }),
        "{error}"
    );

    let fs = memory("boss_plate");
    fs.write(
        "migrations/004_peek.rhai",
        b"let m = previous();\nlet other = open(\"../secret.pkn\");\nm\n",
    )
    .expect("add");
    let mut runner = runner(&fs);
    let built = build(&mut runner, &fs);
    match &built.error {
        Some(ProjectError::MigrationFailed { error, .. }) => {
            assert_eq!(error.payload["reason"], "path_outside_project");
            assert_eq!(error.line, Some(2));
        }
        other => panic!("expected the open() to be refused, got {other:?}"),
    }
}

#[test]
fn misnamed_and_duplicate_migrations_are_refused() {
    let fs = NativeFs::new(fixture("misnamed")).expect("fixture");
    assert!(matches!(
        Project::load(&fs),
        Err(ProjectError::MigrationName { .. })
    ));
    let fs = memory("boss_plate");
    fs.write("migrations/002_again.rhai", b"let m = previous();\nm\n")
        .expect("add");
    assert!(matches!(
        Project::load(fs.as_ref()),
        Err(ProjectError::MigrationOrder { number: 2, .. })
    ));
}

#[test]
fn a_native_project_builds_from_disk() {
    let fs = std::rc::Rc::new(NativeFs::new(fixture("boss_plate")).expect("fixture"));
    let project = Project::load(fs.as_ref()).expect("load");
    let mut runner =
        popkern_project::ProjectRunner::new(fs.clone(), popkern_rhai::HostOptions::default());
    let built = runner.build(&project).expect("build");
    assert!(built.complete(), "{:?}", built.error);
    assert!(
        built
            .versions
            .iter()
            .all(|v| v.status == VersionStatus::Built)
    );
}

#[test]
fn attribution_survives_renaming_the_failing_migration() {
    let fs = memory("labelled_rim");
    let mut runner = runner(&fs);
    assert!(build(&mut runner, &fs).complete());
    // Rename the downstream migration (same content), then drop the
    // label upstream.
    let mut files = fs.files();
    let round = files.remove("migrations/002_round.rhai").expect("002");
    files.insert("migrations/002_fillet_the_rim.rhai".into(), round);
    let block = String::from_utf8(files["migrations/001_block.rhai"].clone()).expect("utf8");
    files.insert(
        "migrations/001_block.rhai".into(),
        block
            .lines()
            .filter(|l| !l.contains("label_edge"))
            .collect::<Vec<_>>()
            .join("\n")
            .into_bytes(),
    );
    let edited = std::rc::Rc::new(popkern_project::MemoryFs::new(files).expect("files"));
    let project = Project::load(edited.as_ref()).expect("load");
    let broken = runner.build(&project).expect("build");
    match &broken.error {
        Some(ProjectError::MigrationFailed {
            migration,
            upstream_changed,
            ..
        }) => {
            assert_eq!(migration, "002_fillet_the_rim.rhai");
            assert_eq!(upstream_changed, &["001_block".to_string()]);
        }
        other => panic!("expected MigrationFailed, got {other:?}"),
    }
}
