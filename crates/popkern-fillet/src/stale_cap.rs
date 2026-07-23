//! PK-143 — stale-cap detection.
//!
//! An in-place edit (radius modify, support drag, imported repair)
//! can leave a corner cap's SURFACE inconsistent with its boundary:
//! the strips move, the cap keeps its old geometry, and the body
//! ships a G0 crack that only surfaces later as mesh opens or STEP
//! sewing failures. This diagnostic measures every cap-like blend
//! face's boundary edges against its own surface and reports the
//! faces whose worst boundary deviation exceeds the gate — the guard
//! callers run BEFORE trusting a cap for further edits (the modify
//! machinery's precondition, and a repair-report input).
//!
//! Pure measurement: no mutation, no repair — detection feeds the
//! [`crate::modify`] re-dispatch decision (a stale cap must go
//! through the record replay, never the in-place path).

use popkern_core::geometry::SurfaceEval;
use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::uv_projection::project_point_to_uv;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};

/// One stale-cap finding.
#[derive(Debug, Clone)]
pub struct StaleCap {
    /// The cap face whose surface disagrees with its boundary.
    pub face: FaceHandle,
    /// Worst boundary-sample distance from the cap surface.
    pub max_deviation: f64,
    /// Number of boundary edges measured.
    pub boundary_edges: usize,
}

/// Detect cap faces whose surface no longer carries their boundary.
///
/// A face is cap-like when its surface is a Sphere or NURBS patch and
/// every boundary edge's far side is a blend strip (Cylinder / Cone /
/// NURBS) or another cap — the corner-cap adjacency signature. Each
/// boundary edge's curve is sampled and projected onto the cap
/// surface; the face is reported when the worst sample exceeds
/// `tolerance.linear() × 100` (the same class of gate the serving
/// validator uses for curve–surface deviation).
pub fn detect_stale_caps(brep: &BRep, tolerance: &Tolerance) -> Vec<StaleCap> {
    let gate = tolerance.linear().max(1e-9) * 100.0;
    let mut out = Vec::new();
    for (fh, f) in brep.faces.iter() {
        let Some(sof) = f.surface.as_ref() else {
            continue;
        };
        if !matches!(
            sof.geometry,
            SurfaceGeometry::Sphere { .. } | SurfaceGeometry::Nurbs(_)
        ) {
            continue;
        }
        // Adjacency signature: every neighbour is blend-family.
        let mut cap_like = true;
        let mut edges = Vec::new();
        let mut loops = vec![f.outer_loop];
        loops.extend(f.inner_loops.iter().copied());
        for lh in loops {
            let Some(l) = brep.loops.get(lh) else {
                continue;
            };
            let start = l.half_edge;
            let mut cur = start;
            for _ in 0..brep.half_edges.len() + 1 {
                let Some(he) = brep.half_edges.get(cur) else {
                    break;
                };
                let Some(tw) = brep.half_edges.get(he.twin) else {
                    break;
                };
                let far = brep.loops.get(tw.loop_).map(|l| l.face);
                let far_blend = far
                    .and_then(|ff| brep.faces.get(ff))
                    .and_then(|fd| fd.surface.as_ref())
                    .map(|s| {
                        matches!(
                            s.geometry,
                            SurfaceGeometry::Cylinder { .. }
                                | SurfaceGeometry::Cone { .. }
                                | SurfaceGeometry::Sphere { .. }
                                | SurfaceGeometry::Torus { .. }
                                | SurfaceGeometry::Nurbs(_)
                        )
                    })
                    .unwrap_or(false);
                if !far_blend {
                    cap_like = false;
                }
                edges.push(he.edge);
                cur = he.next;
                if cur == start {
                    break;
                }
            }
        }
        if !cap_like || edges.is_empty() {
            continue;
        }
        let mut worst = 0.0_f64;
        for eh in &edges {
            let Some(curve) = brep.edges.get(*eh).and_then(|e| e.curve.as_ref()) else {
                continue;
            };
            for k in 0..=8 {
                let t = curve.t_start + (curve.t_end - curve.t_start) * f64::from(k) / 8.0;
                let p = curve.geometry.point_at(t);
                let d = match project_point_to_uv(&sof.geometry, &p) {
                    Some((u, v)) => p.distance_to(&sof.geometry.point_at(u, v)),
                    None => f64::INFINITY,
                };
                worst = worst.max(d);
            }
        }
        if worst > gate {
            out.push(StaleCap {
                face: fh,
                max_deviation: worst,
                boundary_edges: edges.len(),
            });
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{FilletOptions, FilletRadius, FilletSpec, fillet_edges};
    use popkern_core::geometry::Point3;

    fn capped_corner() -> BRep {
        let tol = Tolerance::default();
        let brep = BRep::unit_cube(tol).unwrap();
        let corner = brep
            .vertices
            .iter()
            .find(|(_, v)| {
                (v.position.x - 1.0).abs() < 1e-9
                    && (v.position.y - 1.0).abs() < 1e-9
                    && (v.position.z - 1.0).abs() < 1e-9
            })
            .map(|(vh, _)| vh)
            .unwrap();
        let mut edges = Vec::new();
        for (e, d) in brep.edges.iter() {
            let he = &brep.half_edges[d.half_edge];
            if he.origin == corner || brep.half_edges[he.twin].origin == corner {
                edges.push(e);
            }
        }
        let specs: Vec<FilletSpec> = edges
            .iter()
            .map(|&edge| FilletSpec {
                edge,
                radius: FilletRadius::Constant(0.1),
            })
            .collect();
        fillet_edges(&brep, &specs, &FilletOptions::default(), tol, None)
            .unwrap()
            .brep
    }

    #[test]
    fn fresh_caps_are_not_stale() {
        let brep = capped_corner();
        let stale = detect_stale_caps(&brep, &Tolerance::default());
        assert!(
            stale.is_empty(),
            "freshly built corner caps must measure clean: {stale:?}"
        );
    }

    #[test]
    fn perturbed_cap_is_detected() {
        let mut brep = capped_corner();
        // Stale-ify one cap-like face: displace its sphere centre.
        let mut hit = None;
        for (fh, f) in brep.faces.iter() {
            if let Some(sof) = f.surface.as_ref()
                && matches!(sof.geometry, SurfaceGeometry::Sphere { .. })
            {
                hit = Some(fh);
                break;
            }
        }
        // Fall back to a NURBS cap when the corner dispatched away
        // from the sphere family.
        if hit.is_none() {
            for (fh, f) in brep.faces.iter() {
                if let Some(sof) = f.surface.as_ref()
                    && matches!(sof.geometry, SurfaceGeometry::Nurbs(_))
                {
                    hit = Some(fh);
                    break;
                }
            }
        }
        let fh = hit.expect("capped corner has a cap face");
        if let Some(f) = brep.faces.get_mut(fh)
            && let Some(s) = f.surface.as_mut()
        {
            match &mut s.geometry {
                SurfaceGeometry::Sphere { center, .. } => {
                    *center = Point3::new(center.x + 0.02, center.y, center.z);
                }
                SurfaceGeometry::Nurbs(n) => {
                    for cp in &mut n.control_points {
                        cp.x += 0.02;
                    }
                }
                _ => {}
            }
        }
        let stale = detect_stale_caps(&brep, &Tolerance::default());
        assert!(
            stale.iter().any(|s| s.face == fh && s.max_deviation > 1e-3),
            "displaced cap surface must be detected as stale: {stale:?}"
        );
    }
}
