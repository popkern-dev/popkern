//! `.pkn` load-migration gate (design row 1h, E4): a document written by
//! the version-3 writer loads under the version-4 reader, replays to the
//! body the same recipe builds today, and re-saves at version 4.

mod support;

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_history::{
    CURRENT_VERSION, FeatureTree, HistoryError, MIN_SUPPORTED_VERSION, ReplayOutcome, ReplayResult,
};

/// Frozen bytes from the version-3 writer (see `support::author_reference_document`).
const VERSION_3_DOCUMENT: &[u8] = include_bytes!("fixtures/reference-document-v3.pkn");

fn sorted_vertices(brep: &BRep) -> Vec<[f64; 3]> {
    let mut pts: Vec<[f64; 3]> = brep
        .vertices
        .values()
        .map(|v| [v.position.x, v.position.y, v.position.z])
        .collect();
    pts.sort_by(|a, b| a.partial_cmp(b).expect("finite"));
    pts
}

fn outcomes(tree: &FeatureTree, result: &ReplayResult) -> Vec<(String, String)> {
    tree.feature_labels()
        .into_iter()
        .map(|(label, id)| {
            let outcome = match &result.report.feature(id).expect("report").outcome {
                ReplayOutcome::Served => "served".to_string(),
                ReplayOutcome::Failed { error } => format!("failed: {error}"),
                ReplayOutcome::Suppressed { reason } => format!("suppressed: {reason:?}"),
            };
            (label, outcome)
        })
        .collect()
}

#[test]
fn fixture_is_a_version_3_document() {
    let text = String::from_utf8_lossy(VERSION_3_DOCUMENT);
    assert!(text.starts_with("{\"version\":3,"), "fixture must stay v3");
    assert!(
        !text.contains("naming_scheme"),
        "v3 had no naming-scheme field"
    );
    assert!(!text.contains("\"labels\""), "v3 had no label table");
    assert!(
        !text.contains("build_tolerance"),
        "v3 had no tolerance attribute"
    );
    assert_eq!(MIN_SUPPORTED_VERSION, 3);
    assert_eq!(CURRENT_VERSION, 4);
}

#[test]
fn version_3_document_loads_and_replays_clean_under_the_version_4_reader() {
    let mut loaded = FeatureTree::load(VERSION_3_DOCUMENT).expect("v3 loads");

    // Defaults for every field added since version 3.
    assert_eq!(loaded.build_tolerance(), None);
    assert!(loaded.labels().entities.is_empty());
    let labels: Vec<String> = loaded
        .feature_labels()
        .into_iter()
        .map(|(l, _)| l)
        .collect();
    assert_eq!(
        labels,
        ["Base Sketch", "Block", "Vertical Rounds", "Top Break"]
    );

    let result = loaded.replay(Tolerance::default()).expect("replay");
    assert!(
        result.report.failures().next().is_none() && result.report.fallbacks().next().is_none(),
        "{:?}",
        result.report
    );
    let body = loaded.current_brep().expect("body").clone();
    body.validate().expect("replayed v3 document validates");
    assert_eq!(body.faces.len(), 10, "block + four G2 rounds");

    // Same body the recipe builds today.
    let mut fresh = support::author_reference_document();
    let fresh_result = fresh.replay(Tolerance::default()).expect("replay");
    let fresh_body = fresh.current_brep().expect("body");
    assert_eq!(sorted_vertices(&body), sorted_vertices(fresh_body));
    assert_eq!(body.edges.len(), fresh_body.edges.len());
    assert_eq!(outcomes(&loaded, &result), outcomes(&fresh, &fresh_result));
    assert!(
        outcomes(&loaded, &result)
            .iter()
            .any(|(l, o)| l == "Top Break" && o.starts_with("suppressed")),
        "suppression state survives migration"
    );
    let active = loaded.active_config();
    assert_eq!(
        active.map(|c| loaded.config(c).map(|c| c.name.clone())),
        None
    );
    assert_eq!(
        loaded.configs().len(),
        1,
        "the `tall` configuration survives"
    );
}

#[test]
fn migrated_document_resaves_at_the_current_version_and_round_trips() {
    let mut loaded = FeatureTree::load(VERSION_3_DOCUMENT).expect("v3 loads");
    loaded.replay(Tolerance::default()).expect("replay");
    let body = loaded.current_brep().expect("body").clone();

    let bytes = loaded.save().expect("save");
    let text = String::from_utf8_lossy(&bytes);
    assert!(text.starts_with("{\"version\":4,"), "{}", &text[..40]);
    assert!(text.contains("\"naming_scheme\":1"));

    let mut again = FeatureTree::load(&bytes).expect("v4 loads");
    again.replay(Tolerance::default()).expect("replay");
    let again_body = again.current_brep().expect("body");
    assert_eq!(sorted_vertices(again_body), sorted_vertices(&body));
    assert_eq!(again.feature_labels(), loaded.feature_labels());
}

#[test]
fn versions_outside_the_supported_range_are_refused() {
    let text = String::from_utf8_lossy(VERSION_3_DOCUMENT);
    for (version, bumped) in [
        (2, text.replacen("\"version\":3", "\"version\":2", 1)),
        (5, text.replacen("\"version\":3", "\"version\":5", 1)),
    ] {
        match FeatureTree::load(bumped.as_bytes()) {
            Err(HistoryError::UnsupportedVersion {
                file_version,
                min_supported,
                max_supported,
            }) => {
                assert_eq!(
                    (file_version, min_supported, max_supported),
                    (version, 3, 4)
                );
            }
            other => panic!("version {version}: expected UnsupportedVersion, got {other:?}"),
        }
    }
}
