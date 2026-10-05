//! Structural test: no faer numeric kernels in library code.
//!
//! faer and gemm choose their SIMD kernels by runtime CPU detection
//! (FMA, AVX2, AVX-512 natively; simd128 in wasm), so a product, LU,
//! QR, SVD or sparse factorisation gave different last bits on
//! different machines — measured under PK-676: the gallery's G2 corner
//! caps differed natively and in the browser by 1e-12..7e-10 while every
//! step before them was bit-identical. Library code computes with
//! `popkern_math::linalg` (fixed operation order, no runtime dispatch);
//! faer `Mat` stays usable as a container, and faer's AMD ordering
//! (integer work only) is allowed. Test code is exempt.
//!
//! REMEDIATION: use `popkern_math::linalg` (`matmul`, `matmul_tn`,
//! `syrk_upper`, `lu_solve`, `cholesky`, the triangular solves,
//! `HouseholderQr`, `ColPivQr`, `jacobi_svd`, `LdltSymbolic`).

use std::path::Path;

const FORBIDDEN: &[&str] = &[
    "faer::linalg::",
    "faer::sparse::linalg::lu",
    "faer::sparse::linalg::cholesky",
    "faer::sparse::linalg::qr",
    "SpSolver",
    ".partial_piv_lu(",
    ".full_piv_lu(",
    ".col_piv_qr(",
    ".qr()",
    ".cholesky(faer",
    ".thin_svd()",
    ".svd()",
    ".singular_values()",
    ".selfadjoint_eigendecomposition(",
    ".eigendecomposition",
];

fn workspace_root() -> std::path::PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .canonicalize()
        .expect("Could not resolve workspace root")
}

fn is_test_file(file: &Path) -> bool {
    if file.components().any(|c| {
        c.as_os_str() == "tests" || c.as_os_str() == "benches" || c.as_os_str() == "examples"
    }) {
        return true;
    }
    let name = file
        .file_name()
        .map(|n| n.to_string_lossy().into_owned())
        .unwrap_or_default();
    name == "tests.rs" || name.ends_with("_tests.rs") || name.contains("_tests_part")
}

/// Which lines belong to a `#[cfg(test)]` item. The item runs from
/// the attribute to the brace that closes its body; a braceless item
/// (`mod tests;`, a gated `use`) covers only its own lines — it must
/// not swallow the library code after it.
fn cfg_test_lines(lines: &[&str]) -> Vec<bool> {
    let mut test_only = vec![false; lines.len()];
    let mut i = 0;
    while i < lines.len() {
        if !lines[i].trim().starts_with("#[cfg(test)]") {
            i += 1;
            continue;
        }
        let mut k = i + 1;
        while k < lines.len() && !lines[k].contains('{') && !lines[k].trim_end().ends_with(';') {
            k += 1;
        }
        if k < lines.len() && lines[k].contains('{') {
            let mut depth = 0i64;
            while k < lines.len() {
                depth += lines[k].matches('{').count() as i64;
                depth -= lines[k].matches('}').count() as i64;
                if depth <= 0 {
                    break;
                }
                k += 1;
            }
        }
        let end = k.min(lines.len() - 1);
        for flag in &mut test_only[i..=end] {
            *flag = true;
        }
        i = end + 1;
    }
    test_only
}

#[test]
fn no_faer_numeric_kernels_in_library_code() {
    let files = popkern_lints::find_rs_files(&workspace_root().join("crates"));
    let mut violations = Vec::new();
    for file in &files {
        let path = file.to_string_lossy();
        if path.contains("popkern-examples") || is_test_file(file) {
            continue;
        }
        let Ok(content) = std::fs::read_to_string(file) else {
            continue;
        };
        let lines: Vec<&str> = content.lines().collect();
        let test_only = cfg_test_lines(&lines);
        for (idx, line) in lines.iter().enumerate() {
            if test_only[idx] || line.trim_start().starts_with("//") {
                continue;
            }
            if let Some(pattern) = FORBIDDEN.iter().find(|p| line.contains(**p)) {
                violations.push(format!("{path}:{}: {pattern}: {}", idx + 1, line.trim()));
            }
        }
    }
    assert!(
        violations.is_empty(),
        "faer numeric kernels in library code (their SIMD dispatch depends on the CPU; \
         use popkern_math::linalg):\n{}",
        violations.join("\n")
    );
}
