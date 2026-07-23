//! PK-141 stale-vertex settle: relocate boundary vertices that sit
//! off an adjacent blend-strip surface onto the intersection of
//! their incident surfaces, then re-author the incident curves.
//!
//! The region-merged strips (curved-support apply) inherit boundary
//! remnants whose VERTICES predate the strip surface (measured on
//! the boss chain: seam and corner pieces 2.1e-2 off the strip
//! NURBS; the curve-level reconcile cannot close them because the
//! fit endpoints themselves are stale). This pass alternating-
//! projects each stale vertex onto ALL of its incident surfaces
//! (residual-gated — the linear-crawl lesson), moves it with the
//! chord-rebuilding snap, and re-runs the edge reconcile so the
//! incident curves refit between now-on-surface endpoints.
//!
//! Transactional: the whole pass runs on a candidate clone and
//! commits only when `validate()` holds and every moved vertex's
//! residual improved. Scope: vertices adjacent to a NURBS blend
//! face, off it beyond the serving gate, moved at most one
//! fit-scale cap.

use std::collections::{BTreeMap, BTreeSet};

use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{Point3, SurfaceEval};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle, VertexHandle};

/// Settle stale strip-boundary vertices. Returns how many moved.
pub(crate) fn settle_stale_strip_vertices(brep: &mut BRep, tolerance: &Tolerance) -> usize {
    let trace = std::env::var_os("PK_RECONCILE_TRACE").is_some();
    // Candidate gate at the genuinely-STALE scale (linear x 1000):
    // NURBS fit precision leaves boundary vertices ~2e-4 off their
    // strips, and settling those micro-offsets breaks welded
    // coincidence with topology the refit does not touch (measured:
    // both swept-corner hardening gates shipped welded opens under a
    // linear x 100 gate). The stale class this pass exists for is
    // 1e-2-scale.
    let gate = tolerance.linear().max(1e-9) * 5000.0;
    // Incident faces per vertex.
    let mut incident: BTreeMap<VertexHandle, BTreeSet<FaceHandle>> = BTreeMap::new();
    for (_, he) in brep.half_edges.iter() {
        if let Some(l) = brep.loops.get(he.loop_) {
            incident.entry(he.origin).or_default().insert(l.face);
        }
    }
    let off_surface = |b: &BRep, f: FaceHandle, p: Point3| -> Option<f64> {
        let sof = b.faces.get(f)?.surface.as_ref()?;
        let (u, v) = popkern_core::geometry::uv_projection::project_point_to_uv(&sof.geometry, &p)?;
        Some(p.distance_to(&sof.geometry.point_at(u, v)))
    };
    // Candidates: off an incident NURBS blend face beyond the gate.
    let mut candidates: Vec<(VertexHandle, Point3, f64)> = Vec::new();
    for (v, faces) in &incident {
        let Some(pos) = brep.vertices.get(*v).map(|d| d.position) else {
            continue;
        };
        let mut worst = 0.0_f64;
        let mut has_nurbs = false;
        for &f in faces {
            let is_nurbs = matches!(
                brep.faces.get(f).and_then(|fd| fd.surface.as_ref()),
                Some(s) if matches!(s.geometry, SurfaceGeometry::Nurbs(_))
            );
            // Staleness is measured against the NURBS strips ONLY.
            // Widening to all incident surfaces was FALSIFIED on the
            // boss chain (2026-07-23): corner-territory faces ride the
            // merged strips, so corner vertices are loop-incident to
            // surfaces they geometrically don't lie on (measured 0.32
            // "off" the sphere cap), and tangential strip x support
            // pairs stall the alternating projection at 1.49e-2 (the
            // linear-crawl class). Honest incident sets need the
            // per-corner caps first.
            if let Some(d) = off_surface(brep, f, pos)
                && is_nurbs
            {
                has_nurbs = true;
                worst = worst.max(d);
            }
        }
        if has_nurbs && worst > gate {
            if trace {
                let mut parts: Vec<String> = Vec::new();
                for &f in faces {
                    let kind = match brep.faces.get(f).and_then(|fd| fd.surface.as_ref()) {
                        Some(s) => match &s.geometry {
                            SurfaceGeometry::Plane { .. } => "plane",
                            SurfaceGeometry::Cylinder { .. } => "cyl",
                            SurfaceGeometry::Sphere { .. } => "sphere",
                            SurfaceGeometry::Nurbs(_) => "nurbs",
                            SurfaceGeometry::Torus { .. } => "torus",
                            _ => "other",
                        },
                        None => "none",
                    };
                    let d = off_surface(brep, f, pos);
                    parts.push(format!(
                        "{kind}:{}",
                        d.map_or("-".into(), |d| format!("{d:.2e}"))
                    ));
                }
                eprintln!(
                    "  [vsettle] cand {v:?} ({:.4},{:.4},{:.4}) worst {worst:.3e} [{}]",
                    pos.x,
                    pos.y,
                    pos.z,
                    parts.join(" ")
                );
            }
            candidates.push((*v, pos, worst));
        }
    }
    if candidates.is_empty() {
        return 0;
    }
    let mut cand_brep = brep.clone();
    let mut moved = 0usize;
    for (v, pos, before) in &candidates {
        let Some(faces) = incident.get(v) else {
            continue;
        };
        // Alternating projection across ALL incident surfaces.
        let mut p = *pos;
        let surfaces: Vec<SurfaceGeometry> = faces
            .iter()
            .filter_map(|&f| {
                cand_brep
                    .faces
                    .get(f)
                    .and_then(|fd| fd.surface.as_ref())
                    .map(|s| s.geometry.clone())
            })
            .collect();
        if surfaces.is_empty() {
            continue;
        }
        for _ in 0..200 {
            let mut max_step = 0.0_f64;
            for geom in &surfaces {
                if let Some((u, vv)) =
                    popkern_core::geometry::uv_projection::project_point_to_uv(geom, &p)
                {
                    let q = geom.point_at(u, vv);
                    max_step = max_step.max(p.distance_to(&q));
                    p = q;
                }
            }
            if max_step < gate * 1e-3 {
                break;
            }
        }
        // Residual + travel gates.
        let after = surfaces
            .iter()
            .filter_map(|g| {
                popkern_core::geometry::uv_projection::project_point_to_uv(g, &p)
                    .map(|(u, vv)| p.distance_to(&g.point_at(u, vv)))
            })
            .fold(0.0_f64, f64::max);
        let travel = p.distance_to(pos);
        if after >= *before || travel > before * 4.0 + gate {
            if trace {
                eprintln!(
                    "  [vsettle] {v:?} declined: before {before:.3e} after {after:.3e} travel {travel:.3e}"
                );
            }
            continue;
        }
        super::edge_reconcile::snap_vertex_and_adjacent_chords(
            &mut cand_brep,
            *v,
            p,
            popkern_core::topology::EdgeHandle::default(),
        );
        if trace {
            eprintln!("  [vsettle] {v:?} moved {travel:.3e}: off {before:.3e} -> {after:.3e}");
        }
        moved += 1;
    }
    if moved == 0 {
        return 0;
    }
    // Refit the incident curves between the settled endpoints.
    super::edge_reconcile::reconcile_edge_curves(&mut cand_brep);
    if cand_brep.validate().is_ok() {
        *brep = cand_brep;
        moved
    } else {
        if trace {
            eprintln!("  [vsettle] pass declined: candidate failed validation");
        }
        0
    }
}
