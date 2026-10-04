//! Step-state equivalence (PK-659): the document state the step cache
//! holds after every timeline step equals, bit for bit, the state a
//! replay from scratch produces — on the gallery's multi-body
//! documents, after recording (where every record resumed from the
//! cache), after a mid-tree parameter edit (where the replay resumed
//! at the edited feature), and under a memory budget that evicted
//! states (where `state_at` re-executed from the nearest kept one).

use std::path::Path;

use popkern_assembly::Document;
use popkern_history::{CacheStrategy, DocumentState, FeatureTree};
use popkern_rhai::ScriptHost;
use popkern_script::Model;

fn gallery(name: &str) -> String {
    let path = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../web/live/examples")
        .join(format!("{name}.rhai"));
    std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("{}: {e}", path.display()))
}

fn document(name: &str, source: &str) -> Document {
    let mut host = ScriptHost::new();
    let mut outcome = host
        .run(source, &format!("{name}.rhai"))
        .unwrap_or_else(|e| panic!("{name}: {e}\n{}", e.payload));
    outcome.take_document().expect("a document")
}

/// The source of gallery fixture `name`.
fn source(name: &str) -> String {
    gallery(name)
}

/// The feature tree of every part model fixture `name` authors.
fn trees(name: &str) -> Vec<(String, FeatureTree)> {
    match document(name, &source(name)) {
        Document::Part(model) => vec![(name.to_string(), model.tree().clone())],
        Document::Assembly(mut assembly) => assembly
            .part_models_mut()
            .into_iter()
            .map(|(part, model)| {
                (
                    format!("{name}:{}", part.unwrap_or_default()),
                    model.tree().clone(),
                )
            })
            .collect(),
    }
}

/// The model of the part document `name`.
fn part_model(name: &str) -> Model {
    match document(name, &source(name)) {
        Document::Part(model) => *model,
        Document::Assembly(_) => panic!("{name} is an assembly"),
    }
}

/// A state as comparable text: each body's exact serialization (f64
/// round-trips through JSON exactly) and the step that last changed it.
fn fingerprint(state: &DocumentState) -> Vec<(u32, Option<String>, String)> {
    state
        .bodies()
        .map(|(id, brep)| {
            (
                id.0,
                state.changed_by(id).map(|f| format!("{f:?}")),
                serde_json::to_string(brep).expect("BRep serializes"),
            )
        })
        .collect()
}

/// The cached state after every step of `tree`'s last replay.
fn chain(label: &str, tree: &FeatureTree) -> Vec<Vec<(u32, Option<String>, String)>> {
    (0..=tree.timeline().len())
        .map(|i| {
            fingerprint(
                tree.cached_state_at(i)
                    .unwrap_or_else(|| panic!("{label}: no cached state after step {i}")),
            )
        })
        .collect()
}

/// `tree`'s cached chain equals a from-scratch replay's, step by step.
fn assert_matches_scratch(label: &str, tree: &FeatureTree) {
    let cached = chain(label, tree);
    let mut fresh = tree.clone();
    fresh.clear_cache();
    let tolerance = fresh.effective_build_tolerance();
    let result = fresh
        .replay_from_scratch(tolerance)
        .expect("scratch replay");
    assert_eq!(
        result.cached_steps, 0,
        "{label}: scratch replay executes all"
    );
    let scratch = chain(label, &fresh);
    assert_eq!(cached.len(), scratch.len());
    for (i, (a, b)) in cached.iter().zip(&scratch).enumerate() {
        assert!(
            a == b,
            "{label}: the cached state after step {i} differs from a scratch replay"
        );
    }
}

/// Multi-step, multi-body documents: the boss plate's tool bodies and
/// union, and the gearbox's part definitions. Fastener parts are one
/// step each (the fastener feature), and building their threads takes
/// minutes in the test profile, so the fastener kit runs on demand.
const MULTI_BODY: [&str; 2] = ["boss_plate", "planetary_gearbox"];

#[test]
fn recorded_states_equal_a_scratch_replay() {
    for name in MULTI_BODY {
        for (label, tree) in trees(name) {
            assert_matches_scratch(&label, &tree);
        }
    }
}

/// The full gallery fastener kit (M6/M4/M5 threads): about twelve
/// minutes in the test profile, so it runs on demand
/// (`cargo nextest run -p popkern-rhai --run-ignored only`). Passed
/// when PK-659 landed.
#[test]
#[ignore = "heavy: builds the whole gallery fastener kit (~12 min in the test profile); run on demand — the default gate covers multi-body equivalence on the boss plate and gearbox"]
fn fastener_kit_states_equal_a_scratch_replay() {
    for (label, tree) in trees("fastener_kit") {
        assert_matches_scratch(&label, &tree);
    }
}

#[test]
fn shared_bodies_are_stored_once() {
    for name in MULTI_BODY {
        for (label, tree) in trees(name) {
            let unshared: usize = (1..=tree.timeline().len())
                .filter_map(|i| tree.cached_state_at(i))
                .flat_map(|s| s.bodies().map(|(_, b)| estimate(b)).collect::<Vec<_>>())
                .sum();
            let held = tree.cache_memory_bytes();
            assert!(
                held <= unshared,
                "{label}: the cache ({held} B) never exceeds unshared copies ({unshared} B)"
            );
            let multi = (1..=tree.timeline().len())
                .filter_map(|i| tree.cached_state_at(i))
                .any(|s| s.body_count() > 1);
            if multi {
                assert!(
                    held < unshared,
                    "{label}: multi-body steps share untouched bodies ({held} B vs {unshared} B)"
                );
            }
        }
    }
}

fn estimate(brep: &popkern_core::topology::BRep) -> usize {
    brep.vertices.len() * 64 + brep.edges.len() * 96 + brep.faces.len() * 256
}

#[test]
fn a_mid_tree_edit_resumes_at_the_edited_feature() {
    let mut model = part_model("boss_plate");
    let tree = model.tree();
    let rounds = tree
        .timeline()
        .iter()
        .position(|&id| tree.feature(id).is_some_and(|n| n.label == "rounds"))
        .expect("the rounds fillet");
    model.set_param("rounds.radius", 0.06).expect("edit radius");
    let tolerance = model.tree().effective_build_tolerance();
    let result = model.tree_mut().replay(tolerance).expect("replay");
    assert_eq!(
        result.cached_steps, rounds,
        "every step before the edited fillet is served from the cache"
    );
    assert_eq!(
        result.report.features.len(),
        model.tree().timeline().len(),
        "a whole-document replay reports every feature"
    );
    assert_eq!(
        model.tree().last_report().features.len(),
        model.tree().timeline().len(),
        "the tree's report still covers every feature"
    );
    assert_matches_scratch("boss_plate after the edit", model.tree());

    // The original radius is still cached: switching back executes nothing.
    model
        .set_param("rounds.radius", 0.08)
        .expect("restore radius");
    let back = model.tree_mut().replay(tolerance).expect("replay back");
    assert_eq!(back.cached_steps, model.tree().timeline().len());
    assert_matches_scratch("boss_plate restored", model.tree());
}

#[test]
fn evicted_states_re_execute_to_the_same_state() {
    for name in ["boss_plate", "planetary_gearbox"] {
        for (label, tree) in trees(name) {
            let reference = chain(&label, &tree);
            let mut tree = tree.clone();
            // A budget below one body: nearly every state is dropped.
            tree.set_cache_strategy(CacheStrategy::MemoryBudget { max_bytes: 1 });
            let tolerance = tree.effective_build_tolerance();
            tree.replay(tolerance).expect("replay under budget");
            for (i, expected) in reference.iter().enumerate() {
                let state = tree.state_at(i, tolerance).expect("state_at");
                assert!(
                    &fingerprint(&state) == expected,
                    "{label}: state_at({i}) under a memory budget differs"
                );
            }
        }
    }
}
