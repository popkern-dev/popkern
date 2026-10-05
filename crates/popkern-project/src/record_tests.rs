// Structural lint marker: this file is test code, included via #[path] under #[cfg(test)].

use super::{MigrationRecord, ProjectRecord, Reconcile, read_record, reconcile};
use crate::manifest::Manifest;
use crate::project::{Migration, Project, normalize_source, sha256_hex};

fn migration(number: u64, name: &str, source: &str) -> Migration {
    let source = normalize_source(source);
    Migration {
        number,
        name: name.to_string(),
        path: format!("migrations/{name}.rhai"),
        sha256: sha256_hex(source.as_bytes()),
        source,
    }
}

fn project(migrations: Vec<Migration>) -> Project {
    Project {
        manifest: Manifest::new("t"),
        base: None,
        migrations,
    }
}

fn applied(project: &Project, n: usize) -> ProjectRecord {
    let mut record = ProjectRecord::new(project);
    record.migrations = project.migrations[..n]
        .iter()
        .map(|m| MigrationRecord {
            name: m.name.clone(),
            sha256: m.sha256.clone(),
            features: Vec::new(),
            edited: Vec::new(),
        })
        .collect();
    record
}

#[test]
fn the_reconcile_table() {
    let full = project(vec![
        migration(1, "001_a", "a"),
        migration(2, "002_b", "b"),
        migration(3, "003_c", "c"),
    ]);
    assert!(matches!(reconcile(None, &full), Reconcile::Unrecorded));
    assert!(matches!(
        reconcile(Some(&applied(&full, 3)), &full),
        Reconcile::UpToDate
    ));
    assert!(matches!(
        reconcile(Some(&applied(&full, 1)), &full),
        Reconcile::Pending { from: 1 }
    ));

    let edited = project(vec![
        migration(1, "001_a", "a2"),
        migration(2, "002_b", "b"),
        migration(3, "003_c", "c"),
    ]);
    assert!(matches!(
        reconcile(Some(&applied(&full, 2)), &edited),
        Reconcile::Diverged(_)
    ));

    let renamed = project(vec![
        migration(1, "001_first", "a"),
        migration(2, "002_b", "b"),
        migration(3, "003_c", "c"),
    ]);
    assert!(
        matches!(
            reconcile(Some(&applied(&full, 3)), &renamed),
            Reconcile::UpToDate
        ),
        "a rename keeps the content: still applied"
    );

    let deleted = project(vec![migration(1, "001_a", "a"), migration(3, "003_c", "c")]);
    assert!(matches!(
        reconcile(Some(&applied(&full, 3)), &deleted),
        Reconcile::Missing(_)
    ));
}

#[test]
fn line_endings_and_a_bom_do_not_change_the_hash() {
    assert_eq!(
        migration(1, "001_a", "let m = model();\nm\n").sha256,
        migration(1, "001_a", "\u{feff}let m = model();\r\nm\r\n").sha256
    );
}

#[test]
fn a_record_round_trips_through_the_envelope() {
    let full = project(vec![migration(1, "001_a", "a")]);
    let record = applied(&full, 1);
    let mut model = popkern_script::Model::new();
    model.replay().expect("empty replay");
    let document = popkern_assembly::Document::Part(Box::new(model));
    let bytes = super::to_pkn_bytes(&document, &record).expect("bytes");
    assert_eq!(read_record(&bytes).expect("read"), Some(record));
    popkern_assembly::Document::from_pkn_bytes(&bytes).expect("still a document");
    let plain = document.to_pkn_bytes().expect("plain");
    assert_eq!(read_record(&plain).expect("read"), None);
}
