//! Detects the optional ACM ToMS 1032 package (gitignored — license;
//! see PK-314/PK-315). When present, `cfg(toms_package)` embeds the
//! PnS mask tables and enables the table-dependent tests; when absent
//! (fresh clones, CI), the crate still compiles and table consumers
//! get a structured error instead of a build failure.

fn main() {
    println!("cargo::rustc-check-cfg=cfg(toms_package)");
    let manifest = std::env::var("CARGO_MANIFEST_DIR").unwrap_or_default();
    let package = std::path::Path::new(&manifest)
        .join("../../docs/code/ACM-ToMS-PnS-1032/1032/Src/Source/src/Patch/Table");
    if package.is_dir() {
        println!("cargo::rustc-cfg=toms_package");
    }
    // Track BOTH the package path (content changes) and its parent
    // (appearance/removal): a package dropped into place by mv/untar
    // keeps old mtimes, so only the parent directory's mtime records
    // the transition (measured: a restore invisible to the
    // package-path entry left a stale absent-config build).
    println!("cargo::rerun-if-changed=../../docs/code/ACM-ToMS-PnS-1032");
    println!("cargo::rerun-if-changed=../../docs/code");
}
