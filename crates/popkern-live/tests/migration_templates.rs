//! PK-666: the page's migration templates (web/live/migration-templates.js)
//! are the text popkern-project's `scaffold` writes, so a migration the
//! page adds ("add migration", "insert a feature here") is the one the
//! CLI tools would write.

use std::path::{Path, PathBuf};
use std::process::Command;

use popkern_project::scaffold::{
    first_migration_source, insert_after_name, insert_after_part_source, insert_after_source,
    next_migration_source,
};
use serde_json::{Value, json};

fn root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../..")
}

/// The JS module's answers for `cases` (`[label, position, unique]`).
fn javascript(cases: &Value) -> Value {
    let module = root().join("web/live/migration-templates.js");
    let script = format!(
        "import({:?}).then((m) => console.log(JSON.stringify({{ next: m.NEXT_MIGRATION, \
         first: m.firstMigration(\"Bracket \\\"B\\\"\"), inserts: {cases}.map(([l, p, u]) => [m.insertAfterMigration(l, p, u), m.insertAfterName(l), m.insertAfterPartMigration(\"pin \\\"x\\\"\", l, p, u)]) }})))",
        format!("file://{}", module.display())
    );
    let out = Command::new("node")
        .args(["--input-type=module", "-e", &script])
        .output()
        .expect("node");
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    serde_json::from_slice(&out.stdout).expect("json")
}

#[test]
fn the_page_writes_the_templates_popkern_project_writes() {
    let cases = json!([
        ["rim_round", 4, true],
        ["extrude-1", 2, false],
        ["rim \"round\" \\ 2", 7, true],
        ["Top Face  Fillet!", 3, true],
        ["***", 1, true],
        ["line\nbreak", 5, true],
    ]);
    let js = javascript(&cases);
    assert_eq!(js["next"], next_migration_source());
    assert_eq!(js["first"], first_migration_source("Bracket \"B\""));
    for (i, case) in cases.as_array().expect("cases").iter().enumerate() {
        let label = case[0].as_str().expect("label");
        let position = case[1].as_u64().expect("position") as usize;
        let unique = case[2].as_bool().expect("unique");
        assert_eq!(
            js["inserts"][i][0],
            insert_after_source(label, position, unique),
            "source for {label:?}"
        );
        assert_eq!(
            js["inserts"][i][1],
            insert_after_name(label),
            "name for {label:?}"
        );
        assert_eq!(
            js["inserts"][i][2],
            insert_after_part_source("pin \"x\"", label, position, unique),
            "part source for {label:?}"
        );
    }
}
