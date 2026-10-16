//! Corner plans for the paired chord split (see the parent module):
//! a chord whose other copy is a corner triangle, or a triangle with
//! its apex in the chord's cap zone, is split with that copy re-fanned
//! as a pentagon (PK-859).

use std::collections::HashMap;

use popkern_core::geometry::Point3;
use popkern_core::topology::{BRep, FaceHandle};

use super::super::curve::{CurveSamples, edge_curve_midpoint};
use super::super::{EdgeKey, STEP_SLACK, chart_area, edge_key, triangle_spread, triangle_tilt};
use super::{
    Corner, EdgeMap, EdgeMaps, MeshPlan, edge_triangles, plan_copy, plan_tilt, surface_of,
    vertex_at,
};
use crate::mesh::{FaceMesh, MeshVertex};

/// The plan for copy `(mi, ti, a, b)` of a chord split at `m` when its
/// plain halves lean (its apex sits in the chord's cap zone). The new
/// vertex and the copy's triangle become a pentagon with one more
/// corner, fanned from the corner that leans least:
/// - a second boundary edge of the triangle (a corner, where a seam
///   meets a trim) is split too, at its edge-curve midpoint, and every
///   other copy of it takes the plain split;
/// - else the neighbour across one of the apex's edges joins in, its
///   apex the fifth corner (the flip the cap zone needs, with the split).
///
/// PK-859: the coil's seam caps were refused by their seam copies'
/// corner triangles, 6 -> 35 deg.
#[allow(clippy::too_many_arguments)]
fn plan_corner(
    brep: &BRep,
    meshes: &[FaceMesh],
    copies: &HashMap<EdgeKey, Vec<(usize, usize, usize, usize)>>,
    copy: (usize, usize, usize, usize),
    m: Point3,
    gate: f64,
    bound: f64,
    cache: &mut CurveSamples,
    maps: &mut EdgeMaps,
) -> Option<Vec<MeshPlan>> {
    let (mi, ti, a, b) = copy;
    let fm = &meshes[mi];
    let surface = surface_of(brep, fm)?;
    let t = fm.triangles[ti];
    let r = (0..3).find(|&k| t[k] == a && t[(k + 1) % 3] == b)?;
    let c = t[(r + 2) % 3];
    let edges: &EdgeMap = maps.entry(mi).or_insert_with(|| edge_triangles(fm));
    let old = triangle_tilt(fm, surface, t).unwrap_or(f64::INFINITY);
    let sense = chart_area(fm, t).signum();
    let mid = vertex_at(fm, surface, a, b, m, gate)?;
    // The pentagon's corners in winding order, `Corner::New(1)` the
    // second new vertex when there is one.
    struct Form {
        ring: [Corner; 5],
        remove: Vec<usize>,
        second: Option<((usize, usize), Point3, MeshVertex)>,
        old: f64,
        spread_cap: f64,
    }
    let mut forms: Vec<Form> = Vec::new();
    for (x, y) in [(b, c), (c, a)] {
        let tris = edges
            .get(&(x.min(y), x.max(y)))
            .map(Vec::as_slice)
            .unwrap_or(&[]);
        let ring = |fifth: Corner| {
            if (x, y) == (b, c) {
                [
                    Corner::Old(a),
                    Corner::New(0),
                    Corner::Old(b),
                    fifth,
                    Corner::Old(c),
                ]
            } else {
                [
                    Corner::Old(a),
                    Corner::New(0),
                    Corner::Old(b),
                    Corner::Old(c),
                    fifth,
                ]
            }
        };
        match tris {
            [_] => {
                // A boundary edge: split it too.
                let key2 = edge_key(fm.vertices[x].position, fm.vertices[y].position);
                let Some(list2) = copies.get(&key2) else {
                    continue;
                };
                let faces: Vec<FaceHandle> =
                    list2.iter().map(|&(m2, ..)| meshes[m2].face).collect();
                let Some(m2) = edge_curve_midpoint(
                    brep,
                    &faces,
                    fm.vertices[x].position,
                    fm.vertices[y].position,
                    gate,
                    cache,
                    None,
                ) else {
                    continue;
                };
                let Some(v2) = vertex_at(fm, surface, x, y, m2, gate) else {
                    continue;
                };
                forms.push(Form {
                    ring: ring(Corner::New(1)),
                    remove: vec![ti],
                    second: Some(((x, y), m2, v2)),
                    old,
                    spread_cap: bound.max(triangle_spread(fm, t)),
                });
            }
            [p, q] => {
                // An interior edge: the neighbour's apex joins the fan.
                let tn = if *p == ti { *q } else { *p };
                let Some(d) = fm.triangles[tn].iter().copied().find(|&v| v != x && v != y) else {
                    continue;
                };
                let Some(old_n) = triangle_tilt(fm, surface, fm.triangles[tn]) else {
                    continue;
                };
                forms.push(Form {
                    ring: ring(Corner::Old(d)),
                    remove: vec![ti, tn],
                    second: None,
                    old: old.max(old_n),
                    spread_cap: bound
                        .max(triangle_spread(fm, t))
                        .max(triangle_spread(fm, fm.triangles[tn])),
                });
            }
            _ => {}
        }
    }
    let mut best: Option<(f64, Form, Vec<[Corner; 3]>)> = None;
    for form in forms {
        let mut vertices = vec![mid.clone()];
        if let Some((_, _, v2)) = &form.second {
            vertices.push(v2.clone());
        }
        for apex in 0..5 {
            let tris: Vec<[Corner; 3]> = (1..4)
                .map(|k| {
                    [
                        form.ring[apex],
                        form.ring[(apex + k) % 5],
                        form.ring[(apex + k + 1) % 5],
                    ]
                })
                .collect();
            if let Some(w) = plan_tilt(fm, &vertices, &tris, sense, form.spread_cap)
                && w <= form.old.max(bound) + STEP_SLACK
                && best.as_ref().is_none_or(|x| w < x.0)
            {
                best = Some((
                    w,
                    Form {
                        ring: form.ring,
                        remove: form.remove.clone(),
                        second: form.second.clone(),
                        old: form.old,
                        spread_cap: form.spread_cap,
                    },
                    tris,
                ));
            }
        }
    }
    let (_, form, tris) = best?;
    let mut vertices = vec![mid];
    let mut plans = Vec::new();
    if let Some(((x, y), m2, v2)) = form.second {
        vertices.push(v2);
        let key2 = edge_key(fm.vertices[x].position, fm.vertices[y].position);
        for &cp in copies
            .get(&key2)?
            .iter()
            .filter(|&&(m3, t3, ..)| !(m3 == mi && t3 == ti))
        {
            plans.push(plan_copy(brep, meshes, cp, m2, gate, bound)?);
        }
    }
    plans.insert(
        0,
        MeshPlan {
            mesh: mi,
            remove: form.remove,
            vertices,
            triangles: tris,
        },
    );
    Some(plans)
}

/// Plan splitting leaning cap `(mi, ti)`'s chord `a -> b` at its
/// edge-curve midpoint in every copy, where a copy whose plain halves
/// lean takes the corner plan instead. The cap's own halves must lean
/// less than it does.
#[allow(clippy::too_many_arguments)]
pub(super) fn plan_chord_with_corners(
    brep: &BRep,
    meshes: &[FaceMesh],
    copies: &HashMap<EdgeKey, Vec<(usize, usize, usize, usize)>>,
    (mi, ti): (usize, usize),
    (a, b): (usize, usize),
    gate: f64,
    bound: f64,
    cache: &mut CurveSamples,
    maps: &mut EdgeMaps,
) -> Option<Vec<MeshPlan>> {
    let fm = &meshes[mi];
    let surface = surface_of(brep, fm)?;
    let key = edge_key(fm.vertices[a].position, fm.vertices[b].position);
    let list = copies.get(&key)?;
    let faces: Vec<FaceHandle> = list.iter().map(|&(m2, ..)| meshes[m2].face).collect();
    let m = edge_curve_midpoint(
        brep,
        &faces,
        fm.vertices[a].position,
        fm.vertices[b].position,
        gate,
        cache,
        None,
    )?;
    let old = triangle_tilt(fm, surface, fm.triangles[ti])?;
    let mut plans = Vec::new();
    let mut corner_used = false;
    for &cp in list {
        match plan_copy(brep, meshes, cp, m, gate, bound) {
            Some(p) => plans.push(p),
            None => {
                plans.extend(plan_corner(
                    brep, meshes, copies, cp, m, gate, bound, cache, maps,
                )?);
                corner_used = true;
            }
        }
    }
    // Only a corner plan is new here: a plain split the main pass already tried.
    if !corner_used {
        return None;
    }
    // The leaning cap itself must improve.
    let own = plans
        .iter()
        .find(|p| p.mesh == mi && p.remove.contains(&ti))?;
    let worst = plan_tilt(
        &meshes[mi],
        &own.vertices,
        &own.triangles,
        chart_area(fm, fm.triangles[ti]).signum(),
        f64::INFINITY,
    )?;
    (worst < old).then_some(plans)
}
