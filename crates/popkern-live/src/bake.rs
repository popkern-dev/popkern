//! Pre-built gallery examples: a served, fully certified envelope and
//! its display buffers, written once at build time so the page shows
//! an example the moment it is chosen instead of waiting on the
//! kernel.
//!
//! A bake answers for exactly one script text under exactly one
//! kernel build. The page looks a script up by [`source_key`] and then
//! compares the stored source byte for byte, so a key collision can
//! never show the wrong model. The kernel identity (the wasm bundle's
//! SHA-256, supplied by the build script) is recorded in the index;
//! `tools/build-live.sh` deletes bakes whose identity no longer
//! matches the bundle, so a stale bake is never served.
//!
//! File layout (little-endian): the 8-byte magic [`MAGIC`], a `u32`
//! header length, the header JSON padded with spaces to a multiple of
//! four bytes, then each body's buffers in order — positions, normals
//! (f32), indices, face ids (u32), Gaussian and mean curvature (f32) —
//! with the element counts in the header. Every buffer therefore
//! starts 4-byte aligned, so the page views them in place.

use serde_json::{Value, json};

use crate::envelope::{BodyMesh, Envelope};

/// File magic and format version.
pub const MAGIC: &[u8; 8] = b"PKBAKE01";

/// Index and bake format version.
pub const FORMAT: u32 = 1;

/// The lookup key of a script text: FNV-1a 64 over its UTF-8 bytes,
/// as 16 lowercase hex digits. The page computes the same key; a match
/// is confirmed by comparing the stored source.
pub fn source_key(source: &str) -> String {
    let mut hash: u64 = 0xcbf2_9ce4_8422_2325;
    for byte in source.as_bytes() {
        hash ^= u64::from(*byte);
        hash = hash.wrapping_mul(0x0100_0000_01b3);
    }
    format!("{hash:016x}")
}

/// Native timings recorded beside a bake (what the example costs the
/// kernel, so the page can say so).
#[derive(Debug, Clone, Copy)]
pub struct BakeTimings {
    /// The script run, display meshes included.
    pub run_ms: f64,
    /// Continuity certification of every deferred seam.
    pub certify_ms: f64,
}

/// Errors building a bake.
#[derive(Debug, thiserror::Error)]
pub enum BakeError {
    /// Only a served envelope with a report is worth baking.
    #[error("the envelope is not a served model with a report (outcome {0})")]
    NotServed(String),
    /// The header did not serialise.
    #[error("bake header: {0}")]
    Header(#[from] serde_json::Error),
    /// The header is larger than the format can address.
    #[error("bake header of {0} bytes exceeds the u32 length field")]
    HeaderTooLarge(usize),
}

/// Encode a bake of `source` from its served `envelope`, whose report
/// is `certified_report` (the report once every deferred seam has been
/// certified). `generative` is the run's solved generative designs
/// (`popkern_generative::export_touched`), handed to the page's worker
/// so it serves them without solving (PK-956).
pub fn encode(
    source: &str,
    envelope: &Envelope,
    certified_report: Value,
    timings: BakeTimings,
    generative: Option<&str>,
) -> Result<Vec<u8>, BakeError> {
    if !envelope.outcome.is_served() || envelope.report.is_none() {
        let outcome = serde_json::to_value(&envelope.outcome).unwrap_or(Value::Null);
        return Err(BakeError::NotServed(outcome.to_string()));
    }
    let mut described = envelope.to_json();
    described["report"] = certified_report;
    let counts: Vec<Value> = envelope.bodies.iter().map(buffer_counts).collect();
    let mut header = json!({
        "format": FORMAT,
        "source_key": source_key(source),
        "source": source,
        "envelope": described,
        "buffers": counts,
        "native": { "run_ms": timings.run_ms, "certify_ms": timings.certify_ms },
    });
    if let Some(designs) = generative {
        header["generative"] = serde_json::from_str(designs)?;
    }
    let mut text = serde_json::to_vec(&header)?;
    while text.len() % 4 != 0 {
        text.push(b' ');
    }
    let length = u32::try_from(text.len()).map_err(|_| BakeError::HeaderTooLarge(text.len()))?;
    let mut out = Vec::with_capacity(12 + text.len() + buffer_bytes(&envelope.bodies));
    out.extend_from_slice(MAGIC);
    out.extend_from_slice(&length.to_le_bytes());
    out.extend_from_slice(&text);
    for body in &envelope.bodies {
        for values in [&body.positions, &body.normals] {
            values
                .iter()
                .for_each(|v| out.extend_from_slice(&v.to_le_bytes()));
        }
        for values in [&body.indices, &body.face_ids] {
            values
                .iter()
                .for_each(|v| out.extend_from_slice(&v.to_le_bytes()));
        }
        for values in [&body.gaussian, &body.mean] {
            values
                .iter()
                .for_each(|v| out.extend_from_slice(&v.to_le_bytes()));
        }
    }
    Ok(out)
}

/// The header's element counts for one body, in file order.
fn buffer_counts(body: &BodyMesh) -> Value {
    json!({
        "positions": body.positions.len(),
        "normals": body.normals.len(),
        "indices": body.indices.len(),
        "face_ids": body.face_ids.len(),
        "gaussian": body.gaussian.len(),
        "mean": body.mean.len(),
    })
}

/// Bytes of every body's buffers.
fn buffer_bytes(bodies: &[BodyMesh]) -> usize {
    bodies
        .iter()
        .map(|b| {
            4 * (b.positions.len()
                + b.normals.len()
                + b.indices.len()
                + b.face_ids.len()
                + b.gaussian.len()
                + b.mean.len())
        })
        .sum()
}

/// Split a bake back into its header JSON and the byte offset where
/// the buffers begin (the reader the tests and tools share with the
/// page's `bakes.js`).
pub fn decode_header(bytes: &[u8]) -> Option<(Value, usize)> {
    if bytes.len() < 12 || &bytes[..8] != MAGIC {
        return None;
    }
    let length = u32::from_le_bytes(bytes[8..12].try_into().ok()?) as usize;
    let end = 12usize.checked_add(length)?;
    let header = serde_json::from_slice(bytes.get(12..end)?).ok()?;
    Some((header, end))
}
