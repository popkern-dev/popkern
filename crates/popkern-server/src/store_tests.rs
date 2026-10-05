// Structural lint marker: this file is test code, included via #[path] under #[cfg(test)].

use std::collections::BTreeMap;

use super::{ProjectsRoot, StoreError, rfc3339};

fn temp(name: &str) -> std::path::PathBuf {
    let dir = std::env::temp_dir().join(format!(
        "popkern-server-store-{name}-{}",
        std::process::id()
    ));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).expect("dir");
    dir
}

fn files(source: &str) -> BTreeMap<String, String> {
    BTreeMap::from([
        (
            "popkern.toml".to_string(),
            "[project]\nname = \"Bracket\"\n".to_string(),
        ),
        ("migrations/001_base.rhai".to_string(), source.to_string()),
    ])
}

#[test]
fn projects_round_trip_with_revision_checked_saves() {
    let root = temp("roundtrip");
    let store = ProjectsRoot::new(&root);
    let created = store
        .create("Bracket", &files("let m = model();\nm\n"))
        .expect("create");
    assert_eq!(created.id, "bracket");
    assert_eq!(created.name, "Bracket");
    assert!(root.join("bracket/.gitignore").is_file());
    // A second project with the same name gets its own folder.
    assert_eq!(
        store.create("Bracket", &files("")).expect("again").id,
        "bracket-2"
    );

    let mut next = created.files.clone();
    next.insert(
        "migrations/002_pockets.rhai".into(),
        "let m = previous();\nm\n".into(),
    );
    let saved = store
        .save("bracket", &next, &created.revision)
        .expect("save");
    assert_ne!(saved.revision, created.revision);
    assert!(root.join("bracket/migrations/002_pockets.rhai").is_file());

    // A stale save (the revision it edited is gone) is refused with
    // the stored project; nothing is written.
    match store.save("bracket", &created.files, &created.revision) {
        Err(StoreError::Conflict(current)) => assert_eq!(current.revision, saved.revision),
        other => panic!("expected a conflict, got {other:?}"),
    }
    assert!(root.join("bracket/migrations/002_pockets.rhai").is_file());

    // Removing a migration in a save deletes its file; other files stay.
    std::fs::write(root.join("bracket/notes.txt"), "mine").expect("notes");
    let trimmed = store
        .save("bracket", &created.files, &saved.revision)
        .expect("trim");
    assert!(!root.join("bracket/migrations/002_pockets.rhai").exists());
    assert!(root.join("bracket/notes.txt").is_file());
    assert_eq!(trimmed.files.len(), 2);
    assert_eq!(store.list().expect("list").len(), 2);

    // Delete moves the folder aside.
    store.remove("bracket").expect("remove");
    assert!(store.get("bracket").is_err());
    assert_eq!(
        std::fs::read_dir(root.join(".trash"))
            .expect("trash")
            .count(),
        1
    );
}

#[test]
fn paths_outside_the_exchanged_set_are_refused() {
    let root = temp("paths");
    let store = ProjectsRoot::new(&root);
    let created = store.create("P", &files("")).expect("create");
    for bad in [
        "../escape.rhai",
        "migrations/../../x.rhai",
        "imports/a.step",
        "migrations/sub/a.rhai",
    ] {
        let mut f = created.files.clone();
        f.insert(bad.into(), "x".into());
        assert!(
            matches!(
                store.save("p", &f, &created.revision),
                Err(StoreError::BadPath(_))
            ),
            "{bad}"
        );
    }
    assert!(matches!(store.get("../etc"), Err(StoreError::NotFound(_))));
}

#[test]
fn times_are_rfc3339() {
    assert_eq!(rfc3339(std::time::UNIX_EPOCH), "1970-01-01T00:00:00Z");
    let t = std::time::UNIX_EPOCH + std::time::Duration::from_secs(1_759_622_400);
    assert_eq!(rfc3339(t), "2025-10-05T00:00:00Z");
}
