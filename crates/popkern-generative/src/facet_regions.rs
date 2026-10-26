//! Planar faces from triangles: coplanar patches with their holes,
//! disk-preserving growth where a patch's boundary pinches, and the
//! collinear-vertex pass.

use std::collections::{HashMap, VecDeque};

use super::{V3, cross, dot, sub};

/// One planar face: its owner, its outer cycle (counter-clockwise from
/// outside) and its holes (clockwise).
pub(super) struct Region {
    pub owner: u32,
    pub cycle: Vec<u32>,
    pub holes: Vec<Vec<u32>>,
}

const NONE: u32 = u32::MAX;

/// The planar faces of a closed mesh. Triangles group into coplanar
/// patches, grown from seeds in index order (deterministic): a
/// neighbour joins when it has the seed's owner, its corners lie within
/// `plane_tol` of the seed's plane and its normal agrees. A patch whose
/// boundary is simple cycles, one of them outer, is one face with
/// holes; a patch whose boundary pinches (a vertex it passes twice) is
/// grown again as disks.
pub(super) fn grow(
    vertices: &[V3],
    triangles: &[[u32; 3]],
    owners: &[u32],
    planes: &[(V3, f64)],
    edges: &HashMap<(u32, u32), u32>,
    plane_tol: f64,
) -> Vec<Region> {
    let neighbour = |t: usize, k: usize| -> u32 {
        let tri = triangles[t];
        edges
            .get(&(tri[(k + 1) % 3], tri[k]))
            .copied()
            .unwrap_or(NONE)
    };
    let mut patch = vec![NONE; triangles.len()];
    let mut patches: Vec<(u32, Vec<u32>, V3)> = Vec::new();
    let mut queue: VecDeque<u32> = VecDeque::new();
    for seed in 0..triangles.len() {
        if patch[seed] != NONE {
            continue;
        }
        let p = patches.len() as u32;
        let (n, d) = planes[seed];
        let owner = owners[seed];
        patch[seed] = p;
        let mut members = vec![seed as u32];
        queue.clear();
        queue.extend((0..3).map(|k| neighbour(seed, k)));
        while let Some(t) = queue.pop_front() {
            if t == NONE || patch[t as usize] != NONE {
                continue;
            }
            let ti = t as usize;
            if owners[ti] != owner
                || dot(planes[ti].0, n) < 0.5
                || triangles[ti]
                    .iter()
                    .any(|&v| (dot(n, vertices[v as usize]) - d).abs() > plane_tol)
            {
                continue;
            }
            patch[ti] = p;
            members.push(t);
            queue.extend((0..3).map(|k| neighbour(ti, k)));
        }
        patches.push((owner, members, n));
    }
    let mut regions = Vec::new();
    for (p, (owner, members, n)) in patches.iter().enumerate() {
        match loops(vertices, triangles, members, &patch, p as u32, edges, *n) {
            Some((cycle, holes)) => regions.push(Region {
                owner: *owner,
                cycle,
                holes,
            }),
            None => grow_disks(
                triangles,
                vertices.len(),
                &patch,
                p as u32,
                edges,
                &mut regions,
                *owner,
            ),
        }
    }
    regions
}

/// A patch's boundary as one outer cycle and its holes, or `None` when
/// it pinches (a vertex with two boundary edges leaving it) or has not
/// exactly one cycle winding positively about the patch normal `n`.
fn loops(
    vertices: &[V3],
    triangles: &[[u32; 3]],
    members: &[u32],
    patch: &[u32],
    p: u32,
    edges: &HashMap<(u32, u32), u32>,
    n: V3,
) -> Option<(Vec<u32>, Vec<Vec<u32>>)> {
    let mut next: HashMap<u32, u32> = HashMap::new();
    for &t in members {
        let tri = triangles[t as usize];
        for k in 0..3 {
            let (a, b) = (tri[k], tri[(k + 1) % 3]);
            let inside = edges.get(&(b, a)).is_some_and(|&o| patch[o as usize] == p);
            if !inside && next.insert(a, b).is_some() {
                return None;
            }
        }
    }
    let mut starts: Vec<u32> = next.keys().copied().collect();
    starts.sort_unstable();
    let mut seen: std::collections::HashSet<u32> = std::collections::HashSet::new();
    let mut outer = None;
    let mut holes = Vec::new();
    for start in starts {
        if !seen.insert(start) {
            continue;
        }
        let mut cycle = vec![start];
        let mut v = start;
        while let Some(&w) = next.get(&v) {
            if w == start {
                break;
            }
            if !seen.insert(w) {
                return None;
            }
            cycle.push(w);
            v = w;
        }
        // Newell's area vector against the patch normal.
        let mut area = [0.0; 3];
        for i in 0..cycle.len() {
            let a = vertices[cycle[i] as usize];
            let b = vertices[cycle[(i + 1) % cycle.len()] as usize];
            let c = cross(a, b);
            for d in 0..3 {
                area[d] += c[d];
            }
        }
        if dot(area, n) > 0.0 {
            if outer.is_some() {
                return None;
            }
            outer = Some(cycle);
        } else {
            holes.push(cycle);
        }
    }
    outer.map(|o| (o, holes))
}

/// Grows planar disk regions over the triangles whose `patch` is
/// `only`, in index order (each unassigned triangle seeds one, so the
/// result is deterministic): a neighbour in the patch joins while the
/// region stays a disk, across one shared edge only when its third
/// corner is new to the region, across two (a notch) always, never
/// across three.
fn grow_disks(
    triangles: &[[u32; 3]],
    vertex_count: usize,
    patch: &[u32],
    only: u32,
    edges: &HashMap<(u32, u32), u32>,
    regions: &mut Vec<Region>,
    owner: u32,
) {
    let mut region_of = vec![NONE; triangles.len()];
    // The region a vertex was last seen in (membership of the region
    // being grown only).
    let mut vertex_region = vec![NONE; vertex_count];
    let neighbour = |t: usize, k: usize| -> u32 {
        let tri = triangles[t];
        edges
            .get(&(tri[(k + 1) % 3], tri[k]))
            .copied()
            .unwrap_or(NONE)
    };
    let mut queue: VecDeque<u32> = VecDeque::new();
    for seed in 0..triangles.len() {
        if patch[seed] != only || region_of[seed] != NONE {
            continue;
        }
        let r = regions.len() as u32;
        let mut members = vec![seed as u32];
        region_of[seed] = r;
        for &v in &triangles[seed] {
            vertex_region[v as usize] = r;
        }
        queue.clear();
        queue.extend((0..3).map(|k| neighbour(seed, k)));
        while let Some(t) = queue.pop_front() {
            if t == NONE || region_of[t as usize] != NONE {
                continue;
            }
            let ti = t as usize;
            let tri = triangles[ti];
            if patch[ti] != only {
                continue;
            }
            let shared: Vec<usize> = (0..3)
                .filter(|&k| {
                    let nb = neighbour(ti, k);
                    nb != NONE && region_of[nb as usize] == r
                })
                .collect();
            let joins = match shared.as_slice() {
                // Edge k is (tri[k], tri[k+1]); the third corner is
                // tri[k+2].
                [k] => vertex_region[tri[(k + 2) % 3] as usize] != r,
                [_, _] => true,
                _ => false,
            };
            if !joins {
                continue;
            }
            region_of[ti] = r;
            members.push(t);
            for &v in &tri {
                vertex_region[v as usize] = r;
            }
            queue.extend((0..3).map(|k| neighbour(ti, k)));
        }
        let cycle = boundary(triangles, &members, &region_of, r, edges);
        regions.push(Region {
            owner,
            cycle,
            holes: Vec::new(),
        });
    }
}

/// The boundary cycle of a disk region: its directed edges whose twin
/// lies outside it, chained from the smallest vertex.
fn boundary(
    triangles: &[[u32; 3]],
    members: &[u32],
    region_of: &[u32],
    r: u32,
    edges: &HashMap<(u32, u32), u32>,
) -> Vec<u32> {
    let mut next: HashMap<u32, u32> = HashMap::new();
    for &t in members {
        let tri = triangles[t as usize];
        for k in 0..3 {
            let (a, b) = (tri[k], tri[(k + 1) % 3]);
            let inside = edges
                .get(&(b, a))
                .is_some_and(|&o| region_of[o as usize] == r);
            if !inside {
                next.insert(a, b);
            }
        }
    }
    let Some(&start) = next.keys().min() else {
        return Vec::new();
    };
    let mut cycle = vec![start];
    let mut v = start;
    while let Some(&w) = next.get(&v) {
        if w == start || cycle.len() > next.len() {
            break;
        }
        cycle.push(w);
        v = w;
    }
    cycle
}

/// Drops every vertex with exactly two neighbours in the faces' cycles
/// that lies within `tol` of the segment between them (it splits a
/// straight edge between the same two faces). Returns the count.
pub(super) fn drop_collinear(vertices: &[V3], cycles: &mut [Vec<u32>], tol: f64) -> usize {
    let mut neighbours: HashMap<u32, Vec<u32>> = HashMap::new();
    for c in cycles.iter() {
        let n = c.len();
        for i in 0..n {
            let (a, b) = (c[i], c[(i + 1) % n]);
            for (x, y) in [(a, b), (b, a)] {
                let list = neighbours.entry(x).or_default();
                if !list.contains(&y) {
                    list.push(y);
                }
            }
        }
    }
    let dropped: std::collections::HashSet<u32> = neighbours
        .iter()
        .filter(|(_, list)| list.len() == 2)
        .filter(|&(&v, list)| {
            let (p, a, b) = (
                vertices[v as usize],
                vertices[list[0] as usize],
                vertices[list[1] as usize],
            );
            let ab = sub(b, a);
            let len2 = dot(ab, ab);
            if len2 <= 0.0 {
                return false;
            }
            let ap = sub(p, a);
            let t = dot(ap, ab) / len2;
            let off = cross(ab, ap);
            t > 0.0 && t < 1.0 && dot(off, off) <= tol * tol * len2
        })
        .map(|(&v, _)| v)
        .collect();
    for c in cycles.iter_mut() {
        c.retain(|v| !dropped.contains(v));
    }
    dropped.len()
}
