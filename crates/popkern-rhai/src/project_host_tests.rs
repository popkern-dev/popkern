// Structural lint marker: this file is test code, included via #[path] under #[cfg(test)].

use std::cell::RefCell;
use std::collections::BTreeMap;
use std::path::PathBuf;
use std::rc::Rc;

use popkern_assembly::Document;
use popkern_script::{ProjectReason, ScriptError};

use super::ScriptFiles;
use crate::{HostOptions, ScriptHost};

const BODY: &str = r#"
let m = model();
m.profile(xy()).rect(1.0, 1.0).sketch("plate").extrude(blind(0.2));
"#;

const BASE: &str = r#"
let m = model();
m.profile(xy()).rect(1.0, 1.0).sketch("plate").extrude(blind(0.2));
m
"#;

fn features(document: &Document) -> usize {
    match document {
        Document::Part(model) => model.features().len(),
        Document::Assembly(_) => panic!("part expected"),
    }
}

fn first() -> Document {
    ScriptHost::new()
        .run_migration(BASE, "001_base.rhai", None)
        .expect("first migration")
        .document
}

#[test]
fn the_first_migration_authors_the_document() {
    assert_eq!(features(&first()), 2);
}

#[test]
fn previous_in_a_first_migration_is_refused() {
    let error = ScriptHost::new()
        .run_migration("let m = previous(); m", "001_base.rhai", None)
        .expect_err("no previous document");
    assert_eq!(error.kind, "DocumentError");
    assert_eq!(error.payload["reason"], "no_previous_document");
    assert_eq!(error.line, Some(1), "located at the previous() call");
}

#[test]
fn previous_outside_a_project_is_refused() {
    let error = ScriptHost::new()
        .run("let m = previous(); m", "plain.rhai")
        .expect_err("not a project");
    assert_eq!(error.payload["reason"], "no_previous_document");
}

#[test]
fn a_later_migration_edits_the_previous_document() {
    // Chained without a variable; a tool model left in scope is a tool.
    let source = r#"
let t = model();
previous().profile(xy()).circle(0.1).sketch("pin").extrude(blind(0.5));
"#;
    let out = ScriptHost::new()
        .run_migration(source, "002_pin.rhai", Some(first()))
        .expect("second migration");
    assert_eq!(features(&out.document), 4);
}

#[test]
fn a_later_migration_may_not_replace_the_document() {
    let error = ScriptHost::new()
        .run_migration(BASE, "002_again.rhai", Some(first()))
        .expect_err("replaced");
    assert_eq!(error.payload["reason"], "document_replaced");
}

#[test]
fn a_first_migration_without_a_document_is_refused() {
    let error = ScriptHost::new()
        .run_migration("let x = 1;", "001_empty.rhai", None)
        .expect_err("no document");
    assert_eq!(error.payload["reason"], "no_document");
}

/// An in-memory store that refuses `..`.
#[derive(Default)]
struct Memory(RefCell<BTreeMap<String, Vec<u8>>>);

impl ScriptFiles for Memory {
    fn read(&self, path: &str) -> Result<Vec<u8>, ScriptError> {
        if path.contains("..") {
            return Err(ScriptError::Project {
                reason: ProjectReason::PathOutsideProject,
                path: Some(path.into()),
                detail: format!("{path} is outside the project"),
            });
        }
        self.0
            .borrow()
            .get(path)
            .cloned()
            .ok_or(ScriptError::Project {
                reason: ProjectReason::FileNotFound,
                path: Some(path.into()),
                detail: format!("{path} does not exist"),
            })
    }
    fn write(&self, path: &str, bytes: &[u8]) -> Result<(), ScriptError> {
        self.0.borrow_mut().insert(path.into(), bytes.to_vec());
        Ok(())
    }
    fn native_path(&self, _: &str) -> Result<Option<PathBuf>, ScriptError> {
        Ok(None)
    }
}

#[test]
fn file_verbs_go_through_the_project_store() {
    let files = Rc::new(Memory::default());
    let mut host = ScriptHost::with_files(HostOptions::default(), files.clone());
    let saved = format!("{BODY}m.save(\"out/plate.pkn\");\nm.to_step(\"out/plate.step\");\nm");
    host.run(&saved, "save.rhai").expect("save");
    assert!(files.0.borrow().contains_key("out/plate.pkn"));
    assert!(files.0.borrow().contains_key("out/plate.step"));

    let reopened = host
        .run("let m = open(\"out/plate.pkn\"); m", "open.rhai")
        .expect("open from the store");
    assert_eq!(reopened.model.expect("model").features().len(), 2);

    let escape = host
        .run("let m = open(\"../secret.pkn\"); m", "escape.rhai")
        .expect_err("sandboxed");
    assert_eq!(escape.payload["reason"], "path_outside_project");

    let gate = host
        .run(
            &format!("{BODY}m.assert_step_roundtrip(\"x.step\");"),
            "gate.rhai",
        )
        .expect_err("needs a filesystem");
    assert_eq!(gate.kind, "IoFailed");
}
