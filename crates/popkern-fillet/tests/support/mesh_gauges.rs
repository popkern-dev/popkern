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
            if let popkern_core::error::ValidationError::CurveSurfaceDeviation { edge, .. } = &err {
                let identity = brep.edges.get(*edge).map(|ed| {
                    let he = &brep.half_edges[ed.half_edge];
                    let tw = &brep.half_edges[he.twin];
                    let face_kind = |f: popkern_core::topology::FaceHandle| {
                        brep.faces
                            .get(f)
                            .and_then(|fd| fd.surface.as_ref())
                            .map(|s| format!("{:?}", std::mem::discriminant(&s.geometry)))
                            .unwrap_or_else(|| "none".into())
                    };
                    let a = brep.vertices[he.origin].position;
                    let b = brep.vertices[tw.origin].position;
                    format!(
                        "faces {:?}({})/{:?}({}) span ({:.4},{:.4},{:.4})->({:.4},{:.4},{:.4})",
                        brep.loops[he.loop_].face,
                        face_kind(brep.loops[he.loop_].face),
                        brep.loops[tw.loop_].face,
                        face_kind(brep.loops[tw.loop_].face),
                        a.x,
                        a.y,
                        a.z,
                        b.x,
                        b.y,
                        b.z
                    )
                });
                eprintln!(
                    "[gauge] {err:?} (edge tolerance {:?}; {})",
                    brep.edges.get(*edge).and_then(|e| e.tolerance),
                    identity.unwrap_or_default()
                );
                // Midpoint witness dump: where does each face's chart
                // place the defect curve's midpoint?
                if let Some(ed) = brep.edges.get(*edge)
                    && let Some(cur) = ed.curve.as_ref()
                {
                    use popkern_core::geometry::curve::CurveEval;
                    use popkern_core::geometry::surface::SurfaceEval;
                    for frac in [0.0, 0.25, 0.5, 0.75, 1.0] {
                        let t = cur.t_start + (cur.t_end - cur.t_start) * frac;
                        let pt = cur.geometry.point_at(t);
                        let he0 = &brep.half_edges[ed.half_edge];
                        let tw0 = &brep.half_edges[he0.twin];
                        let mut ds = String::new();
                        for f in [brep.loops[he0.loop_].face, brep.loops[tw0.loop_].face] {
                            if let Some(sof) = brep.faces.get(f).and_then(|fd| fd.surface.as_ref())
                            {
                                let d = popkern_core::geometry::uv_projection::project_point_to_uv(
                                    &sof.geometry,
                                    &pt,
                                )
                                .map(|(u, v)| pt.distance_to(&sof.geometry.point_at(u, v)))
                                .unwrap_or(-1.0);
                                ds.push_str(&format!(" d({f:?})={d:.3e}"));
                            }
                        }
                        eprintln!(
                            "[gauge]   curve f={frac:.2} p=({:.4},{:.4},{:.4}){ds}",
                            pt.x, pt.y, pt.z
                        );
                    }
                    let mid = cur.geometry.point_at(0.5 * (cur.t_start + cur.t_end));
                    let he = &brep.half_edges[ed.half_edge];
                    let tw = &brep.half_edges[he.twin];
                    for f in [brep.loops[he.loop_].face, brep.loops[tw.loop_].face] {
                        if let Some(sof) = brep.faces.get(f).and_then(|fd| fd.surface.as_ref()) {
                            let g = &sof.geometry;
                            let ((u0, u1), (v0, v1)) = g.parameter_range();
                            let w =
                                popkern_core::geometry::uv_projection::project_point_to_uv(g, &mid)
                                    .map(|(u, v)| {
                                        let q = g.point_at(u, v);
                                        format!("uv=({u:.4},{v:.4}) d={:.3e}", mid.distance_to(&q))
                                    })
                                    .unwrap_or_else(|| "proj FAILED".into());
                            eprintln!(
                                "[gauge]   mid on {f:?}: {w} (chart u[{u0:.3},{u1:.3}]                                  v[{v0:.3},{v1:.3}])"
                            );
                        }
                    }
                }
            } else {
                eprintln!("[gauge] {err:?}");
            }
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
