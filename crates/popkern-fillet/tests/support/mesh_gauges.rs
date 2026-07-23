//! Gauges for the hardening matrices.
//!
//! `brep_open_defects` is the WATERTIGHTNESS oracle here: structural
//! plus open-loop/shell-gap counts from `validate_all`. The welded MESH
//! gauge is a measurement only — it counts T-vertex non-conformity
//! along shared edges as "open" (measured: a pristine BRep::prism box
//! reads 358, a pristine unit cube 0, a plain single-edge fillet
//! 26–34 — tessellation conformity, not holes). Product pipelines
//! keep conforming meshes; raw constructor/fillet outputs need not.

use popkern_core::topology::BRep;
use std::collections::HashMap;

pub(crate) fn welded_open_edges(brep: &BRep) -> usize {
    let mesh = popkern_mesh::tessellate(brep, &popkern_mesh::TessellationOptions::default())
        .expect("tessellate");
    let mut key: HashMap<(i64, i64, i64), usize> = HashMap::new();
    let mut id = |p: popkern_core::geometry::Point3| -> usize {
        let k = (
            (p.x * 1.0e7).round() as i64,
            (p.y * 1.0e7).round() as i64,
            (p.z * 1.0e7).round() as i64,
        );
        let next = key.len();
        *key.entry(k).or_insert(next)
    };
    let mut counts: HashMap<(usize, usize), usize> = HashMap::new();
    for t in &mesh.triangles {
        let w = [
            id(mesh.vertices[t[0]].position),
            id(mesh.vertices[t[1]].position),
            id(mesh.vertices[t[2]].position),
        ];
        if w[0] == w[1] || w[1] == w[2] || w[0] == w[2] {
            continue;
        }
        for (a, b) in [(w[0], w[1]), (w[1], w[2]), (w[2], w[0])] {
            *counts.entry((a.min(b), a.max(b))).or_insert(0) += 1;
        }
    }
    counts.values().filter(|&&c| c == 1).count()
}

#[allow(dead_code)]
pub(crate) fn mesh_volume(brep: &BRep) -> f64 {
    let mesh = popkern_mesh::tessellate(brep, &popkern_mesh::TessellationOptions::default())
        .expect("tessellate");
    let mut vol = 0.0_f64;
    for t in &mesh.triangles {
        let (a, b, c) = (
            mesh.vertices[t[0]].position,
            mesh.vertices[t[1]].position,
            mesh.vertices[t[2]].position,
        );
        vol += a.x * (b.y * c.z - b.z * c.y) - a.y * (b.x * c.z - b.z * c.x)
            + a.z * (b.x * c.y - b.y * c.x);
    }
    (vol / 6.0).abs()
}

/// BRep-level defects from `validate_all`: (structural, open-loop/
/// shell-gap, other — incl. self-intersection) error counts. All
/// three must be zero for a body to count as clean; `other` is kept
/// separate so gates can name which class failed.
#[allow(dead_code)] // not every matrix binary uses every gauge
pub(crate) fn brep_open_defects(brep: &BRep) -> (usize, usize, usize) {
    use popkern_core::error::ValidationError;
    let mut structural = 0usize;
    let mut open = 0usize;
    let mut other = 0usize;
    for err in brep.validate_all() {
        if std::env::var_os("PK_GAUGE_TRACE").is_some() {
            eprintln!("[gauge] {err:?}");
        }
        match err {
            ValidationError::EulerPoincare { .. }
            | ValidationError::DanglingReference { .. }
            | ValidationError::BrokenTwin { .. }
            | ValidationError::BrokenCycle { .. }
            | ValidationError::InconsistentVertex { .. }
            | ValidationError::NonManifoldEdge { .. }
            | ValidationError::LowValence { .. } => structural += 1,
            ValidationError::OpenLoop { .. } | ValidationError::ShellGap { .. } => open += 1,
            _ => other += 1,
        }
    }
    (structural, open, other)
}
