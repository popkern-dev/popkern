//! Flap and sliver detection for the vestigial-shell
//! consolidation: the aggregate bridge path's duplicate-coverage
//! signatures (same-plane double-cover flaps, strip far-end sliver
//! patches).

use super::*;

/// A flap: Production face whose plane hosts an interior-overlapping
/// Production partner and whose sides all twin that partner or
/// ConstructionResidual faces.
pub(super) fn detect_flaps(brep: &BRep) -> Vec<FaceHandle> {
    let mut by_plane: BTreeMap<(i64, i64, i64, i64), Vec<FaceHandle>> = BTreeMap::new();
    for (face, data) in &brep.faces {
        if data.role != TopologyRole::Production {
            continue;
        }
        if let Some(key) = plane_key(brep, face) {
            by_plane.entry(key).or_default().push(face);
        }
    }
    let mut flaps = Vec::new();
    for faces in by_plane.values() {
        if faces.len() < 2 {
            continue;
        }
        for &candidate in faces {
            let sides = face_sides(brep, candidate);
            if sides.is_empty() {
                continue;
            }
            // Partner: the single Production face this candidate twins.
            let mut partner = None;
            let mut qualifies = true;
            for side in &sides {
                match twin_face(brep, *side) {
                    Some(f) if f == candidate => {}
                    Some(f)
                        if brep
                            .faces
                            .get(f)
                            .is_some_and(|d| d.role == TopologyRole::ConstructionResidual) => {}
                    Some(f) => match partner {
                        None => partner = Some(f),
                        Some(p) if p == f => {}
                        Some(_) => {
                            qualifies = false;
                            break;
                        }
                    },
                    None => {
                        qualifies = false;
                        break;
                    }
                }
            }
            let Some(partner) = partner else {
                continue;
            };
            if !qualifies
                || plane_key(brep, partner) != plane_key(brep, candidate)
                || !faces_share_plane_interior(brep, candidate, partner)
            {
                continue;
            }
            flaps.push(candidate);
        }
    }
    flaps
}

/// Interior-overlap test on the common plane (grid sampling, strict
/// interior with edge clearance) — the same instrument as the
/// duplicate-coverage serving gate.
fn faces_share_plane_interior(brep: &BRep, a: FaceHandle, b: FaceHandle) -> bool {
    let Some(sof) = brep.faces.get(a).and_then(|d| d.surface.as_ref()) else {
        return false;
    };
    let SurfaceGeometry::Plane { origin, normal, .. } = &sof.geometry else {
        return false;
    };
    let Some(n) = normal.normalized() else {
        return false;
    };
    let e1 = if n.x.abs() < 0.9 {
        popkern_core::geometry::Vec3::new(1.0, 0.0, 0.0)
    } else {
        popkern_core::geometry::Vec3::new(0.0, 1.0, 0.0)
    };
    let Some(e1) = (e1 - n * e1.dot(&n)).normalized() else {
        return false;
    };
    let e2 = n.cross(&e1);
    let ring = |face: FaceHandle| -> Vec<(f64, f64)> {
        brep.face_outer_loop_positions(face)
            .iter()
            .map(|p| {
                let r = *p - *origin;
                (r.dot(&e1), r.dot(&e2))
            })
            .collect()
    };
    let ra = ring(a);
    let rb = ring(b);
    if ra.len() < 3 || rb.len() < 3 {
        return false;
    }
    let inside = |poly: &[(f64, f64)], p: (f64, f64)| -> bool {
        let mut winding = false;
        let mut min_edge = f64::MAX;
        for i in 0..poly.len() {
            let (x0, y0) = poly[i];
            let (x1, y1) = poly[(i + 1) % poly.len()];
            if ((y0 > p.1) != (y1 > p.1)) && (p.0 < (x1 - x0) * (p.1 - y0) / (y1 - y0) + x0) {
                winding = !winding;
            }
            let (dx, dy) = (x1 - x0, y1 - y0);
            let len2 = dx * dx + dy * dy;
            let t = if len2 > 0.0 {
                (((p.0 - x0) * dx + (p.1 - y0) * dy) / len2).clamp(0.0, 1.0)
            } else {
                0.0
            };
            let (cx, cy) = (x0 + t * dx, y0 + t * dy);
            min_edge = min_edge.min(((p.0 - cx).powi(2) + (p.1 - cy).powi(2)).sqrt());
        }
        winding && min_edge > 1.0e-6
    };
    let (mut x0, mut y0, mut x1, mut y1) = (f64::MAX, f64::MAX, f64::MIN, f64::MIN);
    for &(x, y) in &ra {
        x0 = x0.min(x);
        y0 = y0.min(y);
        x1 = x1.max(x);
        y1 = y1.max(y);
    }
    for gx in 1..24 {
        for gy in 1..24 {
            let p = (
                x0 + (x1 - x0) * f64::from(gx) / 24.0,
                y0 + (y1 - y0) * f64::from(gy) / 24.0,
            );
            if inside(&ra, p) && inside(&rb, p) {
                return true;
            }
        }
    }
    false
}

/// A sliver: Production PLANE face with >=2 residual-twinned sides
/// whose remaining twins are all non-plane Production faces (the
/// strip far-end removed-material patch).
pub(super) fn detect_slivers(brep: &BRep) -> Vec<FaceHandle> {
    let mut slivers = Vec::new();
    for (face, data) in &brep.faces {
        if data.role != TopologyRole::Production {
            continue;
        }
        if plane_key(brep, face).is_none() {
            continue;
        }
        let sides = face_sides(brep, face);
        if sides.is_empty() {
            continue;
        }
        let mut residual_sides = 0usize;
        let mut qualifies = true;
        for side in &sides {
            match twin_face(brep, *side) {
                Some(f)
                    if brep
                        .faces
                        .get(f)
                        .is_some_and(|d| d.role == TopologyRole::ConstructionResidual) =>
                {
                    residual_sides += 1;
                }
                Some(f)
                    if brep.faces.get(f).is_some_and(|d| {
                        d.role == TopologyRole::Production
                            && d.surface.as_ref().is_some_and(|s| {
                                !matches!(s.geometry, SurfaceGeometry::Plane { .. })
                            })
                    }) => {}
                _ => {
                    qualifies = false;
                    break;
                }
            }
        }
        if qualifies && residual_sides >= 2 && residual_sides < sides.len() {
            slivers.push(face);
        }
    }
    slivers
}
