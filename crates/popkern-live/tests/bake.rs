//! Pre-built gallery examples: the bake format round-trips a served,
//! certified envelope, its buffers start 4-byte aligned, the lookup key
//! is FNV-1a 64 (the page computes the same), and only served models
//! bake.

use std::path::Path;
use std::time::Duration;

use popkern_live::bake::{self, BakeTimings, MAGIC};
use popkern_live::{ContinuityCertification, HostOptions, Session};

fn boss_plate() -> String {
    let path =
        Path::new(env!("CARGO_MANIFEST_DIR")).join("../popkern-rhai/examples/boss_plate.rhai");
    std::fs::read_to_string(path).expect("boss plate example")
}

const TIMINGS: BakeTimings = BakeTimings {
    run_ms: 1.0,
    certify_ms: 2.0,
};

#[test]
fn source_key_is_fnv1a_64_over_utf8() {
    // Published FNV-1a 64 vectors; web/live/bakes.js must agree.
    assert_eq!(bake::source_key(""), "cbf29ce484222325");
    assert_eq!(bake::source_key("a"), "af63dc4c8601ec8c");
    assert_eq!(bake::source_key("foobar"), "85944171f73967e8");
    // Over UTF-8 bytes, not UTF-16 units: a script with non-ASCII text
    // keys the same on both sides.
    assert_eq!(bake::source_key("é∑"), "ba361110a50f20d6");
}

#[test]
fn a_certified_boss_plate_bakes_and_decodes_with_aligned_buffers() {
    let source = boss_plate();
    let mut session = Session::new(HostOptions {
        continuity_certification: ContinuityCertification::Deferred,
        ..HostOptions::default()
    });
    let envelope = session.run(&source, "boss_plate.rhai");
    assert!(envelope.outcome.is_served());
    let report = loop {
        let step = session.certify(Duration::from_secs(60));
        assert!(step.failure.is_none(), "{:?}", step.failure);
        if step.remaining == 0 {
            break step.report.expect("the final step carries the report");
        }
    };
    let bytes = bake::encode(&source, &envelope, report.clone(), TIMINGS, None).expect("bake");
    assert_eq!(&bytes[..8], MAGIC);
    let (header, start) = bake::decode_header(&bytes).expect("header");
    assert_eq!(start % 4, 0, "buffers start aligned");
    assert_eq!(header["source"], source.as_str());
    assert_eq!(header["source_key"], bake::source_key(&source).as_str());
    assert_eq!(header["envelope"]["report"], report, "the certified report");
    assert_eq!(header["envelope"]["outcome"], "served");
    assert_eq!(header["native"]["certify_ms"], 2.0);
    let body = &envelope.bodies[0];
    let counts = &header["buffers"][0];
    assert_eq!(counts["positions"], body.positions.len());
    assert_eq!(counts["indices"], body.indices.len());
    let expected = 4
        * (body.positions.len()
            + body.normals.len()
            + body.indices.len()
            + body.face_ids.len()
            + body.gaussian.len()
            + body.mean.len());
    assert_eq!(bytes.len() - start, expected, "every buffer, nothing else");
    // The first index reads back where the layout puts it.
    let at = start + 4 * (body.positions.len() + body.normals.len());
    let first = u32::from_le_bytes(bytes[at..at + 4].try_into().expect("4 bytes"));
    assert_eq!(first, body.indices[0]);
}

#[test]
fn only_served_models_bake() {
    let mut session = Session::new(HostOptions::default());
    let envelope = session.run("let m = model(\n", "broken.rhai");
    let error = bake::encode(
        "let m = model(\n",
        &envelope,
        serde_json::Value::Null,
        TIMINGS,
        None,
    );
    assert!(matches!(error, Err(bake::BakeError::NotServed(_))));
    assert!(bake::decode_header(b"not a bake").is_none());
}
