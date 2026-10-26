//! Design histories carried between processes (a gallery bake hands the
//! page the histories of the designs it solved, PK-959), as bytes.
//!
//! Layout (little-endian): the 8-byte magic [`MAGIC`], a `u32` header
//! length, the header JSON padded with spaces to a multiple of four
//! bytes, then per design in header order: the space mesh's positions
//! and normals (`f32`) and indices (`u32`), the active cells (`u32`),
//! and each frame's material (`u8`, one per active cell), padded with
//! zeros to a multiple of four bytes. Counts are in the header.

use serde::{Deserialize, Serialize};

use crate::history::{DesignHistory, DisplayMesh, HistoryFrame};

/// File magic and format version.
pub const MAGIC: &[u8; 8] = b"PKDHST01";

#[derive(Serialize, Deserialize)]
struct FrameHeader {
    update: usize,
    compliance_nmm: f64,
    volume_fraction: f64,
}

#[derive(Serialize, Deserialize)]
struct DesignHeader {
    design: String,
    cell: f64,
    corner: [f64; 3],
    dims: [usize; 3],
    updates: usize,
    active: usize,
    space_vertices: usize,
    space_indices: usize,
    frames: Vec<FrameHeader>,
}

#[derive(Serialize, Deserialize)]
struct Header {
    designs: Vec<DesignHeader>,
}

fn pad(out: &mut Vec<u8>, with: u8) {
    while !out.len().is_multiple_of(4) {
        out.push(with);
    }
}

/// `histories` (design id, history) as bytes.
pub(crate) fn encode(histories: &[(String, &DesignHistory)]) -> Vec<u8> {
    let header = Header {
        designs: histories
            .iter()
            .map(|(design, h)| DesignHeader {
                design: design.clone(),
                cell: h.cell,
                corner: h.corner,
                dims: h.dims,
                updates: h.updates,
                active: h.active.len(),
                space_vertices: h.space.positions.len() / 3,
                space_indices: h.space.indices.len(),
                frames: h
                    .frames
                    .iter()
                    .map(|f| FrameHeader {
                        update: f.update,
                        compliance_nmm: f.compliance_nmm,
                        volume_fraction: f.volume_fraction,
                    })
                    .collect(),
            })
            .collect(),
    };
    let mut text = serde_json::to_vec(&header).unwrap_or_default();
    while !text.len().is_multiple_of(4) {
        text.push(b' ');
    }
    let mut out = Vec::new();
    out.extend_from_slice(MAGIC);
    out.extend_from_slice(&(text.len() as u32).to_le_bytes());
    out.extend_from_slice(&text);
    for (_, h) in histories {
        for v in h.space.positions.iter().chain(&h.space.normals) {
            out.extend_from_slice(&v.to_le_bytes());
        }
        for v in h.space.indices.iter().chain(&h.active) {
            out.extend_from_slice(&v.to_le_bytes());
        }
        for f in &h.frames {
            out.extend_from_slice(&f.material);
        }
        pad(&mut out, 0);
    }
    out
}

/// Reads `bytes[at..at + n]`, moving `at` on.
fn take<'a>(bytes: &'a [u8], at: &mut usize, n: usize) -> Result<&'a [u8], String> {
    let end = at
        .checked_add(n)
        .filter(|&end| end <= bytes.len())
        .ok_or_else(|| format!("design history truncated at byte {at}"))?;
    let out = &bytes[*at..end];
    *at = end;
    Ok(out)
}

fn f32s(bytes: &[u8], at: &mut usize, n: usize) -> Result<Vec<f32>, String> {
    Ok(take(bytes, at, 4 * n)?
        .as_chunks::<4>()
        .0
        .iter()
        .map(|&c| f32::from_le_bytes(c))
        .collect())
}

fn u32s(bytes: &[u8], at: &mut usize, n: usize) -> Result<Vec<u32>, String> {
    Ok(take(bytes, at, 4 * n)?
        .as_chunks::<4>()
        .0
        .iter()
        .map(|&c| u32::from_le_bytes(c))
        .collect())
}

/// The histories `bytes` hold (an [`encode`]), each with its design id.
///
/// # Errors
///
/// A message for bytes of another format, truncated, or describing a
/// history that does not hold together (an active cell outside the
/// grid, an index past the space mesh's vertices).
pub(crate) fn decode(bytes: &[u8]) -> Result<Vec<(String, DesignHistory)>, String> {
    if bytes.len() < 12 || &bytes[..8] != MAGIC {
        return Err("not a design history (magic PKDHST01)".into());
    }
    let length = u32::from_le_bytes([bytes[8], bytes[9], bytes[10], bytes[11]]) as usize;
    let mut at = 12;
    let header: Header =
        serde_json::from_slice(take(bytes, &mut at, length)?).map_err(|e| e.to_string())?;
    let mut out = Vec::with_capacity(header.designs.len());
    for d in header.designs {
        let cells = d
            .dims
            .iter()
            .try_fold(1usize, |n, &k| n.checked_mul(k))
            .ok_or("design history grid too large")?;
        let positions = f32s(bytes, &mut at, 3 * d.space_vertices)?;
        let normals = f32s(bytes, &mut at, 3 * d.space_vertices)?;
        let indices = u32s(bytes, &mut at, d.space_indices)?;
        let active = u32s(bytes, &mut at, d.active)?;
        if indices.iter().any(|&i| i as usize >= d.space_vertices)
            || active.iter().any(|&c| c as usize >= cells)
        {
            return Err(format!("design history {} is inconsistent", d.design));
        }
        let mut frames = Vec::with_capacity(d.frames.len());
        for f in d.frames {
            frames.push(HistoryFrame {
                update: f.update,
                compliance_nmm: f.compliance_nmm,
                volume_fraction: f.volume_fraction,
                material: take(bytes, &mut at, d.active)?.to_vec(),
            });
        }
        at = at.next_multiple_of(4).min(bytes.len());
        out.push((
            d.design,
            DesignHistory {
                cell: d.cell,
                corner: d.corner,
                dims: d.dims,
                active,
                updates: d.updates,
                frames,
                space: DisplayMesh {
                    positions,
                    normals,
                    indices,
                },
            },
        ));
    }
    Ok(out)
}
