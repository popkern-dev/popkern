//! `--stl` (docs/architecture/cli-render-export.md §2): one binary STL
//! of every visible body at its world placement, in document units.
//!
//! STL wants outward, counter-clockwise facets with their own normals.
//! The tessellator's winding is not promised per body, so each body is
//! oriented by its signed volume (a closed body encloses positive
//! volume when its facets face out) and a mirrored placement flips it
//! back; facet normals are the triangles' own, in world space. The
//! welded open-edge census (0 = watertight) and the enclosed volume are
//! reported so a print farm can gate on them.

use popkern_assembly::Document;
use popkern_core::geometry::{Point3, Vec3};
use popkern_mesh::{TessellationOptions, TriMesh};
use serde::Serialize;

/// What an STL export wrote.
#[derive(Debug, Clone, Serialize)]
pub struct StlSummary {
    /// Bodies written (placed occurrences for an assembly).
    pub bodies: usize,
    /// Facets written.
    pub triangles: usize,
    /// Open edges of the welded mesh (0 = watertight).
    pub open_edges: usize,
    /// Volume the facets enclose, in cubic document units.
    pub volume: f64,
}

/// A row-major world placement.
type Placement = [[f64; 4]; 4];

/// A body's mesh and its row-major world placement.
struct Placed<'a> {
    mesh: &'a TriMesh,
    placement: [[f64; 4]; 4],
}

fn apply(m: &[[f64; 4]; 4], p: Point3) -> Point3 {
    Point3::new(
        m[0][0] * p.x + m[0][1] * p.y + m[0][2] * p.z + m[0][3],
        m[1][0] * p.x + m[1][1] * p.y + m[1][2] * p.z + m[1][3],
        m[2][0] * p.x + m[2][1] * p.y + m[2][2] * p.z + m[2][3],
    )
}

fn determinant(m: &[[f64; 4]; 4]) -> f64 {
    m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
}

/// Signed volume of `mesh` in its own frame (divergence theorem).
fn signed_volume(mesh: &TriMesh) -> f64 {
    let o = Point3::new(0.0, 0.0, 0.0);
    mesh.triangles
        .iter()
        .map(|&[a, b, c]| {
            let (p, q, r) = (
                mesh.vertices[a].position - o,
                mesh.vertices[b].position - o,
                mesh.vertices[c].position - o,
            );
            p.dot(&q.cross(&r)) / 6.0
        })
        .sum()
}

/// The world-space facets of `bodies`, each oriented outward.
fn facets(bodies: &[Placed<'_>]) -> Vec<[Point3; 3]> {
    let mut out = Vec::new();
    for body in bodies {
        let inward = signed_volume(body.mesh) < 0.0;
        let mirrored = determinant(&body.placement) < 0.0;
        let flip = inward != mirrored;
        for &[a, b, c] in &body.mesh.triangles {
            let p = |i: usize| apply(&body.placement, body.mesh.vertices[i].position);
            out.push(if flip {
                [p(a), p(c), p(b)]
            } else {
                [p(a), p(b), p(c)]
            });
        }
    }
    out
}

/// Binary STL bytes of `facets` (80-byte header naming popkern, count,
/// then normal and corners as little-endian f32 per facet).
fn encode(facets: &[[Point3; 3]]) -> Vec<u8> {
    let mut bytes = Vec::with_capacity(84 + 50 * facets.len());
    let mut header = [0u8; 80];
    let tag = b"popkern binary STL";
    header[..tag.len()].copy_from_slice(tag);
    bytes.extend_from_slice(&header);
    bytes.extend_from_slice(&(facets.len() as u32).to_le_bytes());
    let push = |bytes: &mut Vec<u8>, v: [f64; 3]| {
        for c in v {
            bytes.extend_from_slice(&(c as f32).to_le_bytes());
        }
    };
    for &[a, b, c] in facets {
        let n = (b - a)
            .cross(&(c - a))
            .normalized()
            .unwrap_or(Vec3::new(0.0, 0.0, 0.0));
        push(&mut bytes, [n.x, n.y, n.z]);
        for p in [a, b, c] {
            push(&mut bytes, [p.x, p.y, p.z]);
        }
        bytes.extend_from_slice(&[0, 0]);
    }
    bytes
}

/// The census of the written facets as one mesh: open edges after the
/// shared 1e-7 weld, and the enclosed volume.
fn census(facets: &[[Point3; 3]]) -> (usize, f64) {
    let mut mesh = TriMesh::default();
    let o = Point3::new(0.0, 0.0, 0.0);
    let mut volume = 0.0;
    for &[a, b, c] in facets {
        let base = mesh.vertices.len();
        for p in [a, b, c] {
            mesh.vertices.push(popkern_mesh::MeshVertex {
                position: p,
                normal: Vec3::new(0.0, 0.0, 0.0),
                uv: (0.0, 0.0),
            });
        }
        mesh.triangles.push([base, base + 1, base + 2]);
        volume += (a - o).dot(&(b - o).cross(&(c - o))) / 6.0;
    }
    (popkern_mesh::welded_open_edges(&mesh), volume)
}

/// The tessellation an `--stl` export uses: the display preset unless
/// `chord` (document units) or `angle_deg` refine it.
pub fn options(chord: Option<f64>, angle_deg: Option<f64>) -> TessellationOptions {
    let mut options = TessellationOptions::display();
    if let Some(chord) = chord {
        options.chord_height = chord;
    }
    if let Some(angle) = angle_deg {
        options.max_facet_angle = angle.to_radians();
    }
    options
}

/// Binary STL of every visible body of `document`.
///
/// # Errors
///
/// A body that does not tessellate, or a document with nothing visible.
pub fn document_stl(
    document: &mut Document,
    options: &TessellationOptions,
) -> Result<(Vec<u8>, StlSummary), String> {
    let identity = [
        [1.0, 0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0, 1.0],
    ];
    let tessellate = |brep: &popkern_core::topology::BRep| {
        popkern_mesh::tessellate(brep, options).map_err(|e| format!("tessellation failed: {e}"))
    };
    let (meshes, rows): (Vec<TriMesh>, Vec<(usize, Placement)>) = match document {
        Document::Part(_) => {
            let part = document.to_part_all().map_err(|e| e.to_string())?;
            (vec![tessellate(part.brep())?], vec![(0, identity)])
        }
        Document::Assembly(assembly) => {
            let parts = assembly.definition_parts().map_err(|e| e.to_string())?;
            let mut meshes = Vec::new();
            let mut index = std::collections::HashMap::new();
            let mut rows = Vec::new();
            for row in assembly
                .leaves()
                .into_iter()
                .filter(|r| r.visible && r.is_part)
            {
                let Some(bodies) = parts.get(&row.definition) else {
                    continue;
                };
                for (i, part) in bodies.iter().enumerate() {
                    let k = match index.get(&(row.definition, i)) {
                        Some(&k) => k,
                        None => {
                            meshes.push(tessellate(part.brep())?);
                            index.insert((row.definition, i), meshes.len() - 1);
                            meshes.len() - 1
                        }
                    };
                    rows.push((k, row.world.transform().matrix));
                }
            }
            (meshes, rows)
        }
    };
    if rows.is_empty() {
        return Err("nothing to export: no visible body".into());
    }
    let placed: Vec<Placed<'_>> = rows
        .iter()
        .map(|&(k, placement)| Placed {
            mesh: &meshes[k],
            placement,
        })
        .collect();
    let facets = facets(&placed);
    let (open_edges, volume) = census(&facets);
    Ok((
        encode(&facets),
        StlSummary {
            bodies: placed.len(),
            triangles: facets.len(),
            open_edges,
            volume,
        },
    ))
}

#[cfg(test)]
#[path = "stl_tests.rs"]
mod tests;
