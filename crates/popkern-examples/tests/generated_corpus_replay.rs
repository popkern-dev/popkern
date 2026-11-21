//! Replays one generated-corpus scene from its descriptor (PK-1227): the
//! command every failure prints. Without `CORPUS_SCENE` it replays
//! nothing and passes.
//!
//!   CORPUS_SCENE='<json>' cargo test -p popkern-examples --test generated_corpus_replay -- --nocapture

#[path = "generated_corpus/mod.rs"]
mod generated_corpus;

use generated_corpus::judge::judge;
use generated_corpus::scene::Scene;

#[test]
fn replay() {
    let Ok(text) = std::env::var("CORPUS_SCENE") else {
        return;
    };
    let scene: Scene = serde_json::from_str(&text).expect("a scene descriptor");
    // A stored scene the generator drew before a fix may no longer meet
    // its independent check: its expectations are void, not the kernel's.
    if let Err(reject) = generated_corpus::scene_check::check(&scene) {
        eprintln!("scene {reject:?}: the descriptor fails its independent check");
        panic!("an invalid scene");
    }
    let verdict = judge(&scene);
    for f in &verdict.failures {
        eprintln!("{} {:?}: {}", f.oracle, f.op, f.detail);
    }
    assert!(
        verdict.failures.is_empty(),
        "{} oracle failures",
        verdict.failures.len()
    );
}
