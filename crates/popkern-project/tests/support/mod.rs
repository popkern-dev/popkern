//! Shared helpers for the fixture-project tests.

use std::path::{Path, PathBuf};
use std::rc::Rc;

use popkern_assembly::Document;
use popkern_project::{MemoryFs, Project, ProjectBuild, ProjectFs, ProjectRunner};
use popkern_rhai::HostOptions;

/// The fixture folder `name`.
pub fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures")
        .join(name)
}

fn collect(root: &Path, dir: &Path, out: &mut Vec<(String, Vec<u8>)>) {
    for entry in std::fs::read_dir(dir).expect("fixture folder") {
        let path = entry.expect("entry").path();
        if path.is_dir() {
            collect(root, &path, out);
        } else {
            let relative = path
                .strip_prefix(root)
                .expect("inside the fixture")
                .to_string_lossy()
                .replace('\\', "/");
            out.push((relative, std::fs::read(&path).expect("read fixture")));
        }
    }
}

/// Fixture `name` copied into memory (tests edit it freely).
pub fn memory(name: &str) -> Rc<MemoryFs> {
    let root = fixture(name);
    let mut files = Vec::new();
    collect(&root, &root, &mut files);
    Rc::new(MemoryFs::new(files).expect("fixture paths"))
}

/// A runner over `fs`.
pub fn runner(fs: &Rc<MemoryFs>) -> ProjectRunner {
    let store: Rc<dyn ProjectFs> = fs.clone();
    ProjectRunner::new(store, HostOptions::default())
}

/// Load and build `fs` with `runner`, panicking on a structural error.
pub fn build(runner: &mut ProjectRunner, fs: &MemoryFs) -> ProjectBuild {
    let project = Project::load(fs).expect("project loads");
    runner.build(&project).expect("build runs")
}

/// Statuses of a build, in order.
pub fn statuses(build: &ProjectBuild) -> Vec<&'static str> {
    build.versions.iter().map(|v| v.status.as_str()).collect()
}

/// Total feature count over every part of `document`.
pub fn feature_count(document: &mut Document) -> usize {
    document
        .models_mut()
        .into_iter()
        .map(|(_, m)| m.features().len())
        .sum()
}

/// The volume of every alive body of a part document.
pub fn volume(document: &mut Document) -> f64 {
    let Document::Part(model) = document else {
        panic!("part document expected");
    };
    model
        .to_parts_all()
        .expect("parts")
        .iter()
        .map(|p| {
            popkern_measure::body_volume(p.brep(), &popkern_core::tolerance::Tolerance::default())
                .expect("volume")
        })
        .sum()
}
