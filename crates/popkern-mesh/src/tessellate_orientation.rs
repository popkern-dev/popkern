//! Coherent orientation of per-face meshes (split per the
//! file-size rule; see [`orient_face_meshes_coherently`]).

use super::{PairWindingVotes, PositionBitKey};
use popkern_core::geometry::{Point3, SurfaceGeometry};
use popkern_core::topology::BRep;

pub(super) fn orient_face_meshes_coherently(
    brep: &BRep,
    face_meshes: &mut [crate::mesh::FaceMesh],
) {
    use std::collections::HashMap;

    let n = face_meshes.len();
    if popkern_core::trace::trace_enabled("mesh-orientation") {
        for (fi, fm) in face_meshes.iter().enumerate() {
            let class = brep
                .faces
                .get(fm.face)
                .and_then(|f| f.surface.as_ref())
                .map(|s| match &s.geometry {
                    SurfaceGeometry::Plane { origin, .. } => format!("Plane z={:.3}", origin.z),
                    SurfaceGeometry::Cylinder { radius, .. } => format!("Cyl r={radius:.3}"),
                    SurfaceGeometry::Cone { apex, .. } => format!("Cone apex_z={:.3}", apex.z),
                    _ => "other".to_string(),
                });
            let mut v6 = 0.0;
            for t in &fm.triangles {
                let (a, b, c) = (
                    fm.vertices[t[0]].position,
                    fm.vertices[t[1]].position,
                    fm.vertices[t[2]].position,
                );
                v6 += a.x * (b.y * c.z - b.z * c.y) - a.y * (b.x * c.z - b.z * c.x)
                    + a.z * (b.x * c.y - b.y * c.x);
            }
            eprintln!(
                "  [orient] mesh {fi} = {:?} {class:?} tris={} verts={} wedge={:.5}",
                fm.face,
                fm.triangles.len(),
                fm.vertices.len(),
                v6 / 6.0
            );
        }
    }
    let mesh_index_of: HashMap<popkern_core::topology::FaceHandle, usize> = face_meshes
        .iter()
        .enumerate()
        .map(|(fi, fm)| (fm.face, fi))
        .collect();

    let boundary_segments: Vec<Vec<(Point3, Point3)>> =
        face_meshes.iter().map(boundary_segments).collect();

    // Each face's segment start points, indexed: a segment can only
    // vote on an edge when its start lies within the projection
    // tolerance of the edge's polyline, so each edge visits only the
    // segments starting inside its windows' grown boxes. Testing every
    // segment of both faces against every edge was quadratic on a
    // thread rod, whose wall has both more edges and more segments per
    // turn (orientation 1.0 s -> 3.6 s at 40 -> 80 mm, PK-794).
    let segment_starts: Vec<crate::point_box_index::PointBoxIndex> = boundary_segments
        .iter()
        .map(|segs| crate::point_box_index::PointBoxIndex::build(segs.iter().map(|s| s.0)))
        .collect();

    // Signed arc-length advance of a face's boundary along a
    // canonical edge polyline. Positive = traverses the edge in
    // increasing parameter, negative = decreasing; magnitude = the
    // summed length of matched segments (the vote's weight).
    let advance_along = |segs: &[(Point3, Point3)],
                         starts: &crate::point_box_index::PointBoxIndex,
                         polyline: &[(f64, Point3)],
                         curve: &dyn popkern_core::geometry::curve::CurveEval|
     -> f64 {
        let mut length = 0.0_f64;
        for w in polyline.windows(2) {
            length += w[0].1.distance_to(&w[1].1);
        }
        if length <= 0.0 {
            return 0.0;
        }
        let tol = (length * 1.0e-3).max(1.0e-9);
        let exact_tol = (length * 1.0e-4).max(1.0e-9);
        // Arc-length parameter of the nearest polyline point and the
        // point's distance from the edge CURVE, or None if the point is
        // farther than tol from every polyline segment.
        let project = |p: &Point3| -> Option<(f64, f64)> {
            let mut best: Option<(f64, f64, f64)> = None; // (dist, s, t)
            let mut s0 = 0.0_f64;
            for w in polyline.windows(2) {
                let (a, b) = (w[0].1, w[1].1);
                let ab = b - a;
                let len = ab.length();
                if len <= 0.0 {
                    continue;
                }
                let t = ((*p - a).dot(&ab) / (len * len)).clamp(0.0, 1.0);
                let q = a + ab * t;
                let d = p.distance_to(&q);
                if best.is_none_or(|(bd, _, _)| d < bd) {
                    best = Some((d, s0 + len * t, w[0].0 + (w[1].0 - w[0].0) * t));
                }
                s0 += len;
            }
            // The polyline is only a coarse stand-in (its chords sag
            // up to ~R·(π/32)²/2 off a circle), so the distance that
            // decides membership is to the curve itself (below).
            let (d, s, t_seed) = best?;
            if d > tol {
                return None;
            }
            let mut t = t_seed;
            for _ in 0..4 {
                let c = curve.point_at(t);
                let d1 = curve.derivative_at(t, 1);
                let d2 = curve.derivative_at(t, 2);
                let r = c - *p;
                let fp = d1.dot(&d1) + r.dot(&d2);
                if fp.abs() <= 1e-300 {
                    break;
                }
                t -= r.dot(&d1) / fp;
            }
            Some((s, curve.point_at(t).distance_to(p)))
        };
        // Closed polylines (rim circles): a raw arc-length delta
        // telescopes to zero around the loop — the wraparound
        // segment's jump (−circumference + step) cancels every other
        // step and the cap↔wall constraint silently vanishes
        // (measured: a cylinder cap stayed inside-out, mesh wedge
        // volume exactly ⅔ of 2πr²h). Wrap each per-segment delta to
        // the short way around.
        let closed = polyline
            .first()
            .zip(polyline.last())
            .is_some_and(|(a, b)| a.1.distance_to(&b.1) <= tol);
        // Candidates: segments starting inside a window's box grown by
        // the tolerance (plus rounding slack), a superset of those
        // whose start `project` accepts. Visiting them in segment
        // order adds the same terms in the same order as a full scan.
        let mut candidates = Vec::new();
        for w in polyline.windows(2) {
            let (a, b) = (w[0].1, w[1].1);
            let grow = tol * (1.0 + 1e-9)
                + 1e-12 * (1.0 + a.x.abs().max(a.y.abs()).max(a.z.abs()))
                + 1e-12 * b.x.abs().max(b.y.abs()).max(b.z.abs());
            let lo = Point3::new(
                a.x.min(b.x) - grow,
                a.y.min(b.y) - grow,
                a.z.min(b.z) - grow,
            );
            let hi = Point3::new(
                a.x.max(b.x) + grow,
                a.y.max(b.y) + grow,
                a.z.max(b.z) + grow,
            );
            starts.query(lo, hi, &mut candidates);
        }
        candidates.sort_unstable();
        candidates.dedup();
        let mut advance = 0.0_f64;
        for (p, q) in candidates.into_iter().map(|k| &segs[k]) {
            if let (Some((sp, dp)), Some((sq, dq))) = (project(p), project(q)) {
                // A segment is ON the edge when the curve passes within
                // a quarter of the segment's own length of both ends (a
                // fallback-grid boundary sags off by about a chord
                // height) or a ten-thousandth of the edge. The polyline
                // tolerance alone admitted a parallel boundary of the
                // same face (a thin annulus's other rim 0.6 away under
                // 1.3-long segments on a 662-long circle), whose
                // opposite traversal cancelled the vote (turbine.stp,
                // PK-648).
                let on_edge = exact_tol.max(0.25 * p.distance_to(q));
                if dp > on_edge || dq > on_edge {
                    continue;
                }
                let mut d = sq - sp;
                if closed {
                    d -= length * (d / length).round();
                }
                // A genuine on-curve chord advances by roughly its own
                // chord length (arc/chord ≤ π/2 even for a semicircle
                // chord). A boundary segment that merely CONNECTS the
                // curve's two endpoints — the closing chord of a
                // crescent face, whose vertices are shared with the
                // edge — projects successfully at both ends yet fakes
                // a near-full-length advance, cancelling the true
                // traversal and flipping the vote (measured: grind-cut
                // cap d = −0.016 instead of +0.854 on a 0.87 rim).
                let chord = p.distance_to(q);
                if d.abs() > chord * 4.0 {
                    continue;
                }
                advance += d;
            }
        }
        advance
    };

    // Face-pair parity votes over shared BRep edges, accumulated by
    // weight (summed matched boundary length). Same-signed advance on
    // both sides means one of the pair is inside-out.
    let mut votes: PairWindingVotes = HashMap::new(); // (same, opposite)
    for (_, edge) in brep.edges.iter() {
        let Some(curve) = edge.curve.as_ref() else {
            continue;
        };
        if !(curve.t_start.is_finite() && curve.t_end.is_finite()) {
            continue;
        }
        let he = edge.half_edge;
        let (Some(h1), Some(h2)) = (
            brep.half_edges.get(he),
            brep.half_edges
                .get(he)
                .and_then(|h| brep.half_edges.get(h.twin)),
        ) else {
            continue;
        };
        let (Some(l1), Some(l2)) = (brep.loops.get(h1.loop_), brep.loops.get(h2.loop_)) else {
            continue;
        };
        if l1.face == l2.face {
            continue; // seam edge: no cross-face constraint
        }
        let (Some(&f1), Some(&f2)) = (mesh_index_of.get(&l1.face), mesh_index_of.get(&l2.face))
        else {
            continue;
        };
        const ORIENT_SAMPLES: usize = 33;
        let eval: &dyn popkern_core::geometry::curve::CurveEval = &curve.geometry;
        let polyline: Vec<(f64, Point3)> = (0..ORIENT_SAMPLES)
            .map(|i| {
                let t = curve.t_start
                    + (curve.t_end - curve.t_start) * i as f64 / (ORIENT_SAMPLES - 1) as f64;
                (t, eval.point_at(t))
            })
            .collect();
        let d1 = advance_along(&boundary_segments[f1], &segment_starts[f1], &polyline, eval);
        let d2 = advance_along(&boundary_segments[f2], &segment_starts[f2], &polyline, eval);
        // Require a non-trivial matched length on BOTH sides before
        // trusting the vote (a lone straggler segment is sew noise).
        let mut edge_len = 0.0;
        for w in polyline.windows(2) {
            edge_len += w[0].1.distance_to(&w[1].1);
        }
        let significant = edge_len * 1.0e-3;
        if d1.abs() < significant || d2.abs() < significant {
            if popkern_core::trace::trace_enabled("mesh-orientation") {
                eprintln!(
                    "  [orient] pair ({f1},{f2}) advance below floor: d1={d1:.6} d2={d2:.6} \
                     edge_len={edge_len:.4}"
                );
            }
            continue;
        }
        let weight = d1.abs().min(d2.abs());
        if popkern_core::trace::trace_enabled("mesh-orientation") {
            eprintln!(
                "  [orient] vote ({f1},{f2}) {} w={weight:.5} d1={d1:.5} d2={d2:.5} \
                 edge_len={edge_len:.4}",
                if d1 * d2 > 0.0 { "same" } else { "opposite" }
            );
        }
        let entry = votes.entry((f1.min(f2), f1.max(f2))).or_default();
        if d1 * d2 > 0.0 {
            entry.0 += weight;
        } else {
            entry.1 += weight;
        }
    }

    // Majority-resolved adjacency, sorted for deterministic traversal.
    let mut adj: Vec<Vec<(usize, bool)>> = vec![Vec::new(); n]; // (other, same_direction)
    let mut pairs: Vec<_> = votes.iter().collect();
    pairs.sort_by_key(|(k, _)| **k);
    for (&(f1, f2), &(same, opposite)) in pairs {
        if same == opposite {
            // Tied evidence — no trustworthy constraint. Leave the
            // pair unlinked rather than guess.
            if popkern_core::trace::trace_enabled("mesh-orientation") {
                eprintln!("  [orient] pair ({f1},{f2}) tied votes {same}:{opposite} — unlinked");
            }
            continue;
        }
        let same_dir = same > opposite;
        if same.min(opposite) > 0.0 && popkern_core::trace::trace_enabled("mesh-orientation") {
            eprintln!(
                "  [orient] pair ({f1},{f2}) conflicted votes same={same:.3} opposite={opposite:.3} — majority {}",
                if same_dir { "same" } else { "opposite" }
            );
        }
        adj[f1].push((f2, same_dir));
        adj[f2].push((f1, same_dir));
    }
    // BFS parity propagation per component (deterministic: seeds in
    // index order, adjacency sorted above).
    let mut flip: Vec<Option<bool>> = vec![None; n];
    for seed in 0..n {
        if flip[seed].is_some() {
            continue;
        }
        flip[seed] = Some(false);
        let mut queue = std::collections::VecDeque::from([seed]);
        let mut component = vec![seed];
        while let Some(f) = queue.pop_front() {
            let f_flip = flip[f].unwrap_or(false);
            for &(g, same_dir) in &adj[f] {
                // Opposite traversal (same_dir=false) is CONSISTENT:
                // g keeps f's parity. Same-direction traversal means
                // one of the two is inside-out relative to the other.
                let want = if same_dir { !f_flip } else { f_flip };
                match flip[g] {
                    None => {
                        flip[g] = Some(want);
                        component.push(g);
                        queue.push_back(g);
                    }
                    Some(assigned) => {
                        if assigned != want
                            && popkern_core::trace::trace_enabled("mesh-orientation")
                        {
                            eprintln!(
                                "  [orient] parity cycle conflict at faces ({f},{g}) — keeping first assignment"
                            );
                        }
                    }
                }
            }
        }
        // Component-global side: an OUTER boundary component encloses
        // positive volume (outward facing); a VOID boundary — a
        // component enclosed inside another component of the same body
        // (ray-crossing parity, position-only, so it needs no facing
        // decisions) — encloses negative volume, facing into the void.
        // The former unconditional-positive rule meshed every hollow
        // body's inner skin inside-out (see winding_side).
        let mut vol6 = 0.0;
        for &f in &component {
            let s = if flip[f] == Some(true) { -1.0 } else { 1.0 };
            let fm = &face_meshes[f];
            for t in &fm.triangles {
                let (a, b, c) = (
                    fm.vertices[t[0]].position,
                    fm.vertices[t[1]].position,
                    fm.vertices[t[2]].position,
                );
                vol6 += s
                    * (a.x * (b.y * c.z - b.z * c.y) - a.y * (b.x * c.z - b.z * c.x)
                        + a.z * (b.x * c.y - b.y * c.x));
            }
        }
        let enclosed = crate::winding_side::component_is_enclosed(&component, face_meshes);
        // Only a CLOSED component's signed volume says which side it
        // faces. An open one — faces cut off from the rest of the
        // body by missing votes — keeps the winding its faces' flags
        // gave it (measured: turbine.stp's fillet strip, unlinked
        // from its planes, was flipped on an open strip's volume and
        // the solid meshed 66 % over; PK-648).
        let closed = component_is_closed(brep, &component, face_meshes);
        let wrong_side = closed && if enclosed { vol6 > 0.0 } else { vol6 < 0.0 };
        if popkern_core::trace::trace_enabled("mesh-orientation") {
            eprintln!(
                "  [orient] component {component:?} vol6 {vol6:.4} enclosed={enclosed} \
                 closed={closed} wrong_side={wrong_side}"
            );
        }
        if wrong_side {
            for &f in &component {
                flip[f] = flip[f].map(|x| !x);
            }
        }
    }
    for (fi, fm) in face_meshes.iter_mut().enumerate() {
        if flip[fi] == Some(true) {
            for t in &mut fm.triangles {
                t.swap(1, 2);
            }
            for v in &mut fm.vertices {
                v.normal = -v.normal;
            }
        }
    }
}

/// The face mesh's boundary as directed segments: triangle edges used
/// once, welded by position so duplicated boundary vertices at one
/// position still match. Sorted by position bits: the counting map
/// iterates in its per-instance random order, and the votes sum these
/// segments' advances, so an unsorted list made the vote weights vary
/// in the last bits from run to run (PK-796).
fn boundary_segments(fm: &crate::mesh::FaceMesh) -> Vec<(Point3, Point3)> {
    use std::collections::HashMap;
    let key = |p: &Point3| -> PositionBitKey { (p.x.to_bits(), p.y.to_bits(), p.z.to_bits()) };
    let mut counts: HashMap<(PositionBitKey, PositionBitKey), i64> = HashMap::new();
    for t in &fm.triangles {
        for (a, b) in [(t[0], t[1]), (t[1], t[2]), (t[2], t[0])] {
            let (ka, kb) = (key(&fm.vertices[a].position), key(&fm.vertices[b].position));
            if ka == kb {
                continue;
            }
            if ka < kb {
                *counts.entry((ka, kb)).or_default() += 1;
            } else {
                *counts.entry((kb, ka)).or_default() -= 1;
            }
        }
    }
    let from_bits = |k: (u64, u64, u64)| -> Point3 {
        Point3::new(
            f64::from_bits(k.0),
            f64::from_bits(k.1),
            f64::from_bits(k.2),
        )
    };
    let mut segments: Vec<(Point3, Point3)> = counts
        .into_iter()
        .filter(|&(_, c)| c != 0)
        .map(|((ka, kb), c)| {
            if c > 0 {
                (from_bits(ka), from_bits(kb))
            } else {
                (from_bits(kb), from_bits(ka))
            }
        })
        .collect();
    segments.sort_unstable_by_key(|(p, q)| (key(p), key(q)));
    segments
}

/// Whether a component is topologically closed: every edge of every
/// member face has its other face in the component too. (Balancing
/// the meshes' boundary segments read closed bodies as open — face
/// meshes still carry T-junctions here, before the crack zipper — and
/// left them inside-out: jack_in_the_box's 3-face bodies and io1-tc
/// meshed at negative volume, caught by the corpus OCC volume gate.)
fn component_is_closed(
    brep: &BRep,
    component: &[usize],
    face_meshes: &[crate::mesh::FaceMesh],
) -> bool {
    let faces: std::collections::HashSet<_> =
        component.iter().map(|&f| face_meshes[f].face).collect();
    faces.iter().all(|&face| {
        let Some(data) = brep.faces.get(face) else {
            return false;
        };
        std::iter::once(data.outer_loop)
            .chain(data.inner_loops.iter().copied())
            .all(|lh| {
                let Some(start) = brep.loops.get(lh).map(|l| l.half_edge) else {
                    return false;
                };
                let mut he = start;
                for _ in 0..=brep.half_edges.len() {
                    let Some(h) = brep.half_edges.get(he) else {
                        return false;
                    };
                    let other = brep
                        .half_edges
                        .get(h.twin)
                        .and_then(|t| brep.loops.get(t.loop_))
                        .map(|l| l.face);
                    if !other.is_some_and(|o| faces.contains(&o)) {
                        return false;
                    }
                    he = h.next;
                    if he == start {
                        break;
                    }
                }
                true
            })
    })
}

#[cfg(test)]
#[path = "tessellate_orientation_tests.rs"]
mod tests;
