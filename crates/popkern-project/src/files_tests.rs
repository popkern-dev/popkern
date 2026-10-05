// Structural lint marker: this file is test code, included via #[path] under #[cfg(test)].

use super::{MemoryFs, ProjectFs, normalize};
use crate::error::ProjectError;

#[test]
fn paths_normalise_and_never_leave_the_project() {
    assert_eq!(
        normalize("./migrations/001_a.rhai").expect("ok"),
        "migrations/001_a.rhai"
    );
    for outside in ["../x", "a/../../x", "/etc/passwd", "a\\..\\b", "", "."] {
        assert!(
            matches!(
                normalize(outside),
                Err(ProjectError::PathOutsideProject { .. })
            ),
            "{outside:?} must be refused"
        );
    }
}

#[test]
fn memory_lists_only_direct_children() {
    let fs = MemoryFs::new([
        ("migrations/001_a.rhai".to_string(), b"a".to_vec()),
        ("migrations/sub/002_b.rhai".to_string(), b"b".to_vec()),
        ("popkern.toml".to_string(), b"t".to_vec()),
    ])
    .expect("files");
    assert_eq!(fs.list("migrations").expect("list"), ["001_a.rhai"]);
    assert!(fs.list("nowhere").expect("list").is_empty());
    assert!(matches!(
        fs.read("missing.rhai"),
        Err(ProjectError::FileNotFound { .. })
    ));
}

#[test]
fn a_symlink_out_of_a_native_project_is_refused() {
    let root = std::env::temp_dir().join(format!("popkern-project-symlink-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&root);
    std::fs::create_dir_all(root.join("project")).expect("mkdir");
    std::fs::write(root.join("secret.pkn"), b"{}").expect("secret");
    #[cfg(unix)]
    {
        std::os::unix::fs::symlink(root.join("secret.pkn"), root.join("project/link.pkn"))
            .expect("symlink");
        let fs = super::NativeFs::new(root.join("project")).expect("fs");
        assert!(matches!(
            fs.read("link.pkn"),
            Err(ProjectError::PathOutsideProject { .. })
        ));
    }
    let _ = std::fs::remove_dir_all(&root);
}
