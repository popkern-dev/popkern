//! The generated corpus's full tier (PK-1227): one test per kind pair,
//! run by the nextest profile `foundations` (not the per-commit suite).
//! `CORPUS_SEED` and `CORPUS_RANDOM` set the seed and the random scenes
//! per cell; each failure prints its replay descriptor.

#[path = "generated_corpus/mod.rs"]
mod generated_corpus;

use generated_corpus::operand::Kind::{self, *};
use generated_corpus::report::{random_count, run_pair};

fn pair(a: Kind, b: Kind) {
    // The full tier fails, not skips, without its OCC oracle (design §7).
    assert!(
        generated_corpus::occ::python().is_some(),
        "the OCC oracle (~/venvs/ocp) is required by the full tier"
    );
    let (unexpected, fixed) = run_pair(a, b, random_count(300), false);
    let mut problems = unexpected;
    problems.extend(fixed);
    assert!(
        problems.is_empty(),
        "{} problems:\n{}",
        problems.len(),
        problems.join("\n")
    );
}

macro_rules! pairs {
    ($($name:ident: $a:ident $b:ident),* $(,)?) => {
        $(#[test] fn $name() { pair($a, $b); })*
    };
}

pairs!(
    p_p: P P, p_y: P Y, p_c: P C, p_f: P F, p_s: P S, p_t: P T, p_r: P R, p_q: P Q, p_b: P B, p_w: P W, p_l: P L, p_g: P G,
    y_y: Y Y, y_c: Y C, y_f: Y F, y_s: Y S, y_t: Y T, y_r: Y R, y_q: Y Q, y_b: Y B, y_w: Y W, y_l: Y L, y_g: Y G,
    c_c: C C, c_f: C F, c_s: C S, c_t: C T, c_r: C R, c_q: C Q, c_b: C B, c_w: C W, c_l: C L, c_g: C G,
    f_f: F F, f_s: F S, f_t: F T, f_r: F R, f_q: F Q, f_b: F B, f_w: F W, f_l: F L, f_g: F G,
    s_s: S S, s_t: S T, s_r: S R, s_q: S Q, s_b: S B, s_w: S W, s_l: S L, s_g: S G,
    t_t: T T, t_r: T R, t_q: T Q, t_b: T B, t_w: T W, t_l: T L, t_g: T G,
    r_r: R R, r_q: R Q, r_b: R B, r_w: R W, r_l: R L, r_g: R G,
    q_q: Q Q, q_b: Q B, q_w: Q W, q_l: Q L, q_g: Q G,
    b_b: B B, b_w: B W, b_l: B L, b_g: B G,
    w_w: W W, w_l: W L, w_g: W G,
    l_l: L L, l_g: L G,
    g_g: G G,
);

fn sheets(solid: Kind) {
    let unexpected = generated_corpus::report::run_sheets(solid, random_count(300) / 10);
    assert!(
        unexpected.is_empty(),
        "{} problems:\n{}",
        unexpected.len(),
        unexpected.join("\n")
    );
}

macro_rules! sheet_tests {
    ($($name:ident: $k:ident),* $(,)?) => {
        $(#[test] fn $name() { sheets($k); })*
    };
}

sheet_tests!(
    h_p: P, h_y: Y, h_c: C, h_f: F, h_s: S, h_t: T, h_r: R, h_q: Q, h_b: B, h_w: W, h_l: L, h_g: G,
);

/// Every curated tangency corpus case under random placements and
/// self-turns (template cells, design §4).
#[test]
fn templates() {
    use generated_corpus::known::KNOWN;
    let (failures, runs, passed) = generated_corpus::templates::run_templates(
        generated_corpus::report::seed(),
        random_count(300) / 10,
    );
    eprintln!("templates: {runs} runs, {passed} passed");
    let unexpected: Vec<String> = failures
        .iter()
        .filter(|(key, _)| !KNOWN.iter().any(|(k, _)| k == key))
        .map(|(key, detail)| format!("{key}: {detail}"))
        .collect();
    assert!(
        unexpected.is_empty(),
        "{} problems:\n{}",
        unexpected.len(),
        unexpected.join("\n")
    );
}
