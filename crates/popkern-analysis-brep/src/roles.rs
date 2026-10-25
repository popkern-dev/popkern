//! Blend roles inferred from how a solid's faces join (the rules are
//! documented in [`crate::features`]): fillets on cylinders, tori and
//! G1 freeform faces, chamfers on narrow oblique planar strips (no two
//! sharing an edge).

use std::collections::BTreeMap;

use popkern_analysis::{FaceId, FaceInfo, FaceRole, RoleSource, Surface};

use crate::features::{FeatureTolerances, NormalAt, PairJoin, angle};
use crate::geometry::{norm, sub};

/// Largest turning angle of a chamfer's join to either neighbour, rad.
const CHAMFER_MAX_TURN: f64 = 85.0 * std::f64::consts::PI / 180.0;

/// How closely a chamfer's two turnings must add up to the angle between
/// its neighbours, rad.
const CHAMFER_SUM_TOL: f64 = 1.0 * std::f64::consts::PI / 180.0;

/// A chamfer strip is at most this wide relative to its length.
const CHAMFER_ASPECT: f64 = 0.5;

/// One solid's welded mesh.
pub(crate) struct Mesh<'a> {
    pub vertices: &'a [[f64; 3]],
    pub triangles: &'a [[u32; 3]],
    pub triangle_face: &'a [FaceId],
}

/// The smallest radius of curvature of freeform face `face` seen along
/// its tessellation's edges, from exact normals at their ends.
fn freeform_radius(mesh: &Mesh<'_>, face: FaceId, freeform: &NormalAt<'_>) -> Option<f64> {
    let mut normals: BTreeMap<u32, Option<[f64; 3]>> = BTreeMap::new();
    let mut kappa = 0.0f64;
    for (t, f) in mesh.triangles.iter().zip(mesh.triangle_face) {
        if *f != face {
            continue;
        }
        for k in 0..3 {
            let (a, b) = (t[k], t[(k + 1) % 3]);
            let mut normal = |v: u32| {
                *normals
                    .entry(v)
                    .or_insert_with(|| freeform(face, mesh.vertices[v as usize]))
            };
            let (Some(na), Some(nb)) = (normal(a), normal(b)) else {
                continue;
            };
            let length = norm(sub(mesh.vertices[b as usize], mesh.vertices[a as usize]));
            // Normals from the evaluator may disagree in sense.
            let turn = angle(na, nb).min(angle(na, nb.map(|c| -c)));
            if length > 0.0 {
                kappa = kappa.max(turn / length);
            }
        }
    }
    (kappa > 0.0).then(|| 1.0 / kappa)
}

/// Sets the inferred fillet and chamfer roles of `faces`.
pub(crate) fn infer_roles(
    mesh: &Mesh<'_>,
    faces: &mut [FaceInfo],
    pairs: &BTreeMap<(FaceId, FaceId), PairJoin>,
    tolerances: &FeatureTolerances,
    freeform: &NormalAt<'_>,
) {
    let mut tangent_neighbours: BTreeMap<FaceId, usize> = BTreeMap::new();
    for (&(fa, fb), join) in pairs {
        if join.tangent {
            *tangent_neighbours.entry(fa).or_default() += 1;
            *tangent_neighbours.entry(fb).or_default() += 1;
        }
    }
    for face in faces.iter_mut() {
        if tangent_neighbours.get(&face.id).copied().unwrap_or(0) < 2 {
            continue;
        }
        let radius = match face.surface {
            Surface::Cylinder { radius, .. } => Some(radius),
            Surface::Torus { minor_radius, .. } => Some(minor_radius),
            Surface::Freeform => freeform_radius(mesh, face.id, freeform),
            _ => None,
        };
        if let Some(radius) = radius {
            face.role = Some(FaceRole::Fillet { radius });
            face.role_source = RoleSource::Inferred;
        }
    }
    for f in chamfers(faces, pairs, tolerances) {
        if let Some(face) = faces.iter_mut().find(|x| x.id == f) {
            face.role = Some(FaceRole::Chamfer);
            face.role_source = RoleSource::Inferred;
        }
    }
}

/// Most chamfer candidates in one connected group searched exhaustively.
const EXHAUSTIVE: usize = 16;

/// Chamfer faces: planar strips that meet exactly two planar faces
/// obliquely with the same sense (sharp, all convex or all concave,
/// under `CHAMFER_MAX_TURN`), whose turnings add up to the angle between
/// those faces, and that are narrow. A face between two chamfers passes
/// these tests too, but chamfers bevel edges and never share one, so
/// each group of edge-adjacent candidates keeps its largest set of
/// non-adjacent faces (lowest ids on ties).
fn chamfers(
    faces: &[FaceInfo],
    pairs: &BTreeMap<(FaceId, FaceId), PairJoin>,
    tolerances: &FeatureTolerances,
) -> Vec<FaceId> {
    let plane_normal = |f: FaceId| {
        faces
            .iter()
            .find(|x| x.id == f)
            .and_then(|x| match x.surface {
                Surface::Plane { normal, .. } => Some(normal),
                _ => None,
            })
    };
    let mut candidates: Vec<FaceId> = Vec::new();
    for face in faces {
        if face.role.is_some() || plane_normal(face.id).is_none() {
            continue;
        }
        let oblique: Vec<(FaceId, f64, f64, bool)> = pairs
            .iter()
            .filter_map(|(&(fa, fb), join)| {
                let other = if fa == face.id {
                    fb
                } else if fb == face.id {
                    fa
                } else {
                    return None;
                };
                let turn = join.all_turning / join.segments.max(1) as f64;
                let sharp = turn >= tolerances.sharp_angle && turn <= CHAMFER_MAX_TURN;
                let uniform = join.sharp_convex || join.concave.len() == join.segments;
                (sharp && uniform).then_some((other, turn, join.length, join.sharp_convex))
            })
            .collect();
        let [(a, ta, la, ca), (b, tb, lb, cb)] = oblique.as_slice() else {
            continue;
        };
        let (Some(na), Some(nb)) = (plane_normal(*a), plane_normal(*b)) else {
            continue;
        };
        let length = la.max(*lb);
        if ca == cb
            && (angle(na, nb) - (ta + tb)).abs() <= CHAMFER_SUM_TOL
            && length > 0.0
            && face.area / length <= CHAMFER_ASPECT * length
        {
            candidates.push(face.id);
        }
    }
    let adjacent = |x: FaceId, y: FaceId| pairs.contains_key(&(x.min(y), x.max(y)));
    let mut out = Vec::new();
    let mut done = vec![false; candidates.len()];
    for start in 0..candidates.len() {
        if done[start] {
            continue;
        }
        // The connected group of candidates around `start`.
        let mut group = vec![start];
        done[start] = true;
        let mut k = 0;
        while k < group.len() {
            for j in 0..candidates.len() {
                if !done[j] && adjacent(candidates[group[k]], candidates[j]) {
                    done[j] = true;
                    group.push(j);
                }
            }
            k += 1;
        }
        group.sort_unstable();
        let n = group.len();
        let independent = |mask: u32| {
            (0..n).all(|i| {
                mask & (1 << i) == 0
                    || (i + 1..n).all(|j| {
                        mask & (1 << j) == 0
                            || !adjacent(candidates[group[i]], candidates[group[j]])
                    })
            })
        };
        let best = if n <= EXHAUSTIVE {
            // Most faces; on ties the set whose lowest ids come first.
            (0..(1u32 << n))
                .filter(|&m| independent(m))
                .max_by_key(|&m| (m.count_ones(), m.reverse_bits()))
                .unwrap_or(0)
        } else {
            // Greedy in id order beyond the exhaustive size.
            let mut mask = 0u32;
            for i in 0..n.min(32) {
                if independent(mask | (1 << i)) {
                    mask |= 1 << i;
                }
            }
            mask
        };
        out.extend(
            (0..n)
                .filter(|i| best & (1 << i) != 0)
                .map(|i| candidates[group[i]]),
        );
    }
    out.sort_unstable();
    out
}
