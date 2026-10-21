//! Slit and point-loop retirement passes for boolean outputs (split
//! from `spur_edge.rs` under the 500-line file budget): same-face
//! slit-pair excision, pure-slit face retirement, and point-loop
//! retirement (PK-465).

use popkern_core::topology::BRep;

use crate::serve_region::ServeRegion;

/// Excise interwoven same-face slit PAIRS: an edge both of whose
/// half-edges lie in the SAME loop of a NON-closed face is an
/// out-and-back spur woven into a real boundary (the walk-level
/// InteriorGraze class surviving to the brep — measured: the angled
/// tap's crest face carried its discarded above-plate weave inside
/// its kept outer loop, and the next boolean rejected the body as an
/// open shell, PK-465). Cycle surgery removes the pair: the span
/// between the two occurrences becomes its own loop (a genuine
/// enclosed region reached through the slit), empty spans drop.
/// Closed/periodic faces keep their seams (a seam has exactly this
/// half-edge shape and is load-bearing).
pub fn excise_same_face_slit_pairs(brep: &mut BRep) -> usize {
    excise_same_face_slit_pairs_in(brep, &ServeRegion::Whole)
}

/// [`excise_same_face_slit_pairs`] on the faces `region` admits (the
/// surgery edits one face; a loop it splits off stays on that face).
pub(crate) fn excise_same_face_slit_pairs_in(brep: &mut BRep, region: &ServeRegion) -> usize {
    use popkern_core::geometry::surface::SurfaceGeometry;
    let face_is_closed = |brep: &BRep, fh: popkern_core::topology::FaceHandle| -> bool {
        let Some(face) = brep.faces.get(fh) else {
            return true; // unknown: leave alone
        };
        let Some(surf) = face.surface.as_ref() else {
            return true;
        };
        match &surf.geometry {
            SurfaceGeometry::Cylinder { .. }
            | SurfaceGeometry::Cone { .. }
            | SurfaceGeometry::Sphere { .. }
            | SurfaceGeometry::Torus { .. } => true,
            SurfaceGeometry::Nurbs(_) => {
                use popkern_core::geometry::surface::SurfaceEval;
                let sfc = &surf.geometry;
                let ((u0, u1), (v0, v1)) = sfc.parameter_range();
                let gate = brep.tolerance.linear().max(1.0e-9) * 10.0;
                let closed_u = (0..=4).all(|i| {
                    let v = v0 + (v1 - v0) * i as f64 / 4.0;
                    sfc.point_at(u0, v).distance_to(&sfc.point_at(u1, v)) < gate
                });
                let closed_v = (0..=4).all(|i| {
                    let u = u0 + (u1 - u0) * i as f64 / 4.0;
                    sfc.point_at(u, v0).distance_to(&sfc.point_at(u, v1)) < gate
                });
                closed_u || closed_v
            }
            _ => false,
        }
    };
    let mut removed = 0usize;
    let mut skip: std::collections::BTreeSet<popkern_core::topology::EdgeHandle> =
        std::collections::BTreeSet::new();
    for _ in 0..brep.edges.len().max(1) {
        // Find one slit pair: both half-edges on the same face
        // (same loop, or two loops of one face), skipping candidates
        // whose surgery already failed.
        let mut target = None;
        for (eh, edge) in &brep.edges {
            if skip.contains(&eh) {
                continue;
            }
            let Some(he) = brep.half_edges.get(edge.half_edge) else {
                continue;
            };
            let Some(tw) = brep.half_edges.get(he.twin) else {
                continue;
            };
            if he.twin == edge.half_edge {
                continue;
            }
            let (Some(l1), Some(l2)) = (brep.loops.get(he.loop_), brep.loops.get(tw.loop_)) else {
                continue;
            };
            if l1.face != l2.face || !region.admits_face(l1.face) {
                continue;
            }
            if face_is_closed(brep, l1.face) {
                if popkern_core::trace::trace_enabled("spur") {
                    eprintln!(
                        "[spur] slit candidate {eh:?} on {:?} skipped: face closed",
                        l1.face
                    );
                }
                continue;
            }
            target = Some((eh, edge.half_edge, he.twin, he.loop_, tw.loop_, l1.face));
            break;
        }
        let Some((eh, h, t, lh1, lh2, fh)) = target else {
            break;
        };
        skip.insert(eh);
        let snapshot = brep.clone();
        let collect_cycle =
            |brep: &BRep, lh: popkern_core::topology::LoopHandle| -> Option<Vec<_>> {
                let l = brep.loops.get(lh)?;
                let mut cycle = Vec::new();
                let mut cur = l.half_edge;
                for _ in 0..brep.half_edges.len().max(1) {
                    cycle.push(cur);
                    cur = brep.half_edges.get(cur)?.next;
                    if cur == l.half_edge {
                        return Some(cycle);
                    }
                }
                None
            };
        if lh1 == lh2 {
            // Same loop: shed the slit; the span between the two
            // occurrences becomes its own loop.
            let Some(cycle) = collect_cycle(brep, lh1) else {
                continue;
            };
            let (Some(i), Some(j)) = (
                cycle.iter().position(|&x| x == h),
                cycle.iter().position(|&x| x == t),
            ) else {
                continue;
            };
            let (i, j) = if i < j { (i, j) } else { (j, i) };
            let main: Vec<_> = cycle[..i]
                .iter()
                .chain(cycle[j + 1..].iter())
                .copied()
                .collect();
            let side: Vec<_> = cycle[i + 1..j].to_vec();
            if main.is_empty() {
                continue; // whole loop is the slit: the pure-slit passes own it
            }
            for w in 0..main.len() {
                let a = main[w];
                let b = main[(w + 1) % main.len()];
                brep.half_edges[a].next = b;
                brep.half_edges[b].prev = a;
            }
            if let Some(l) = brep.loops.get_mut(lh1) {
                l.half_edge = main[0];
            }
            if !side.is_empty() {
                let new_loop = brep.loops.insert(popkern_core::topology::Loop {
                    half_edge: side[0],
                    face: fh,
                });
                for w in 0..side.len() {
                    let a = side[w];
                    let b = side[(w + 1) % side.len()];
                    brep.half_edges[a].next = b;
                    brep.half_edges[b].prev = a;
                    brep.half_edges[a].loop_ = new_loop;
                }
                if let Some(face) = brep.faces.get_mut(fh) {
                    face.inner_loops.push(new_loop);
                }
            }
        } else {
            // Two loops of one face: removing the pair merges them.
            let (Some(ca), Some(cb)) = (collect_cycle(brep, lh1), collect_cycle(brep, lh2)) else {
                continue;
            };
            let (Some(i), Some(j)) = (
                ca.iter().position(|&x| x == h),
                cb.iter().position(|&x| x == t),
            ) else {
                continue;
            };
            let merged: Vec<_> = ca[..i]
                .iter()
                .chain(cb[j + 1..].iter())
                .chain(cb[..j].iter())
                .chain(ca[i + 1..].iter())
                .copied()
                .collect();
            if merged.is_empty() {
                continue;
            }
            for w in 0..merged.len() {
                let a = merged[w];
                let b = merged[(w + 1) % merged.len()];
                brep.half_edges[a].next = b;
                brep.half_edges[b].prev = a;
                brep.half_edges[a].loop_ = lh1;
            }
            if let Some(l) = brep.loops.get_mut(lh1) {
                l.half_edge = merged[0];
            }
            // Retire lh2 (whichever role it had).
            if let Some(face) = brep.faces.get_mut(fh) {
                if face.outer_loop == lh2 {
                    // The outer merged into lh1: promote lh1.
                    face.inner_loops.retain(|&x| x != lh1);
                    face.outer_loop = lh1;
                } else {
                    face.inner_loops.retain(|&x| x != lh2);
                }
            }
            brep.loops.remove(lh2);
        }
        // Drop the pair + edge; re-anchor vertices.
        let verts = [brep.half_edges[h].origin, brep.half_edges[t].origin];
        brep.remove_half_edge_reanchored(h);
        brep.remove_half_edge_reanchored(t);
        brep.edges.remove(eh);
        for v in verts {
            match brep
                .half_edges
                .iter()
                .find(|(_, he)| he.origin == v)
                .map(|(hh, _)| hh)
            {
                Some(out_he) => {
                    if let Some(vd) = brep.vertices.get_mut(v) {
                        vd.outgoing = out_he;
                    }
                }
                None => {
                    brep.vertices.remove(v);
                }
            }
        }
        let before = snapshot.validate_collect_structural().len();
        let after = brep.validate_collect_structural().len();
        if after > before {
            if popkern_core::trace::trace_enabled("spur") {
                eprintln!(
                    "[spur] slit pair {h:?}/{t:?} ROLLBACK: {:?}",
                    brep.validate_collect_structural().first()
                );
            }
            *brep = snapshot;
            continue;
        }
        if popkern_core::trace::trace_enabled("spur") {
            eprintln!("[spur] same-face slit pair {h:?}/{t:?} excised from loop {lh1:?}");
        }
        removed += 1;
    }
    removed
}

/// Retire faces whose EVERY edge borders the face on BOTH sides — a
/// pure slit-weave enclosing zero area (a discarded region's
/// boundary walked out-and-back and mis-kept; measured: the angled
/// tap's above-plate crest weave served as a 12-edge zero-area face
/// and the next boolean rejected the body as an open shell,
/// PK-465). Removing such a face removes only its own loops and
/// edges: by construction nothing else references them.
pub fn retire_pure_slit_faces(brep: &mut BRep) -> usize {
    retire_pure_slit_faces_in(brep, &ServeRegion::Whole)
}

/// [`retire_pure_slit_faces`] on the faces `region` admits (a slit
/// loop or face retires alone).
pub(crate) fn retire_pure_slit_faces_in(brep: &mut BRep, region: &ServeRegion) -> usize {
    let mut removed = 0usize;
    loop {
        // A SLIT LOOP: every half-edge's twin lies in the SAME loop —
        // an out-and-back weave enclosing zero area. It can be a
        // whole face (every loop slit) or a single bogus loop
        // attached to a real face; either way the slit loop bounds
        // nothing and retires whole.
        let mut target: Option<(
            popkern_core::topology::FaceHandle,
            Vec<popkern_core::topology::LoopHandle>,
        )> = None;
        'faces: for (fh, face) in brep.faces.iter().filter(|&(f, _)| region.admits_face(f)) {
            let mut all_loops = vec![face.outer_loop];
            all_loops.extend(face.inner_loops.iter().copied());
            for &lh in &all_loops {
                let Some(l) = brep.loops.get(lh) else {
                    continue;
                };
                let mut heh = l.half_edge;
                let mut is_slit = true;
                let mut any = false;
                for _ in 0..brep.half_edges.len().max(1) {
                    let Some(he) = brep.half_edges.get(heh) else {
                        is_slit = false;
                        break;
                    };
                    let Some(tw) = brep.half_edges.get(he.twin) else {
                        is_slit = false;
                        break;
                    };
                    if tw.loop_ != lh || he.twin == heh {
                        is_slit = false;
                        break;
                    }
                    any = true;
                    heh = he.next;
                    if heh == l.half_edge {
                        break;
                    }
                }
                if is_slit && any {
                    if lh == face.outer_loop {
                        // A slit OUTER loop retires the whole face
                        // only when every other loop is also slit —
                        // otherwise leave it visible (louder bug).
                        let all_slit = all_loops.len() == 1;
                        if all_slit {
                            target = Some((fh, all_loops.clone()));
                            break 'faces;
                        }
                        continue;
                    }
                    target = Some((fh, vec![lh]));
                    break 'faces;
                }
            }
        }
        let Some((fh, loops)) = target else {
            break;
        };
        let whole_face = brep
            .faces
            .get(fh)
            .map(|f| loops.contains(&f.outer_loop))
            .unwrap_or(false);
        if popkern_core::trace::trace_enabled("spur") {
            eprintln!(
                "[spur] pure-slit {} {fh:?} retired ({} loop(s))",
                if whole_face { "face" } else { "loop of" },
                loops.len()
            );
        }
        // Collect all hes/edges/vertices of the face.
        let mut hes = Vec::new();
        for &lh in &loops {
            if let Some(l) = brep.loops.get(lh) {
                let mut heh = l.half_edge;
                for _ in 0..brep.half_edges.len().max(1) {
                    let Some(he) = brep.half_edges.get(heh) else {
                        break;
                    };
                    hes.push(heh);
                    heh = he.next;
                    if heh == l.half_edge {
                        break;
                    }
                }
            }
        }
        let edges: std::collections::BTreeSet<_> = hes
            .iter()
            .filter_map(|&h| brep.half_edges.get(h).map(|he| he.edge))
            .collect();
        let verts: std::collections::BTreeSet<_> = hes
            .iter()
            .filter_map(|&h| brep.half_edges.get(h).map(|he| he.origin))
            .collect();
        for h in hes {
            brep.remove_half_edge_reanchored(h);
        }
        for e in edges {
            brep.edges.remove(e);
        }
        for &lh in &loops {
            if let Some(face) = brep.faces.get_mut(fh) {
                face.inner_loops.retain(|&x| x != lh);
            }
            brep.loops.remove(lh);
        }
        if whole_face {
            brep.faces.remove(fh);
        }
        for v in verts {
            match brep
                .half_edges
                .iter()
                .find(|(_, he)| he.origin == v)
                .map(|(hh, _)| hh)
            {
                Some(out_he) => {
                    if let Some(vd) = brep.vertices.get_mut(v) {
                        vd.outgoing = out_he;
                    }
                }
                None => {
                    brep.vertices.remove(v);
                }
            }
        }
        removed += 1;
    }
    removed
}
