//! Assembly cleanup and provenance: unwired-half-edge removal and
//! output-face provenance construction.
//!
//! Split from `cellular/assemble.rs` to keep it within the
//! structural line ratchet.

use std::collections::HashMap;

use popkern_core::constructors::builder::BRepBuilder;
use popkern_core::topology::{BRep, HalfEdgeHandle, VertexHandle};

use crate::merge::fragment_classify::FragmentSelection;
use crate::merge::graph::Operand;
use crate::merge::walk::Fragment;

/// Remove half-edges not wired into any loop, and their orphaned edges/vertices.
pub(super) fn cleanup_unwired(builder: &mut BRepBuilder) {
    let brep = builder.brep_mut();
    let default_loop = popkern_core::topology::LoopHandle::default();

    let dangling: Vec<HalfEdgeHandle> = brep
        .half_edges
        .iter()
        .filter(|(_, he)| he.loop_ == default_loop)
        .map(|(h, _)| h)
        .collect();

    let twins: Vec<(HalfEdgeHandle, HalfEdgeHandle)> = dangling
        .iter()
        .filter_map(|&h| brep.half_edges.get(h).map(|he| (h, he.twin)))
        .collect();

    for (h, twin_h) in &twins {
        if let Some(twin) = brep.half_edges.get_mut(*twin_h) {
            twin.twin = *twin_h;
        }
        // The surviving edge KEEPS its curve: edge geometry is a
        // property of the edge, not of its twin wiring — clearing it
        // here silently degraded seam-style survivors (a boss/plate
        // intersection circle lost its geometry entirely).
        brep.remove_half_edge_reanchored(*h);
    }

    // Remove fully-orphaned edges
    let dead_edges: Vec<_> = brep
        .edges
        .iter()
        .filter(|(_, e)| !brep.half_edges.contains_key(e.half_edge))
        .filter(|(eh, _)| !brep.half_edges.iter().any(|(_, he)| he.edge == *eh))
        .map(|(eh, _)| eh)
        .collect();
    for eh in dead_edges {
        brep.edges.remove(eh);
    }

    // Fix vertex outgoing handles
    let outgoing: HashMap<VertexHandle, HalfEdgeHandle> = brep
        .half_edges
        .iter()
        .map(|(h, he)| (he.origin, h))
        .collect();
    for (vh, v) in &mut brep.vertices {
        if !brep.half_edges.contains_key(v.outgoing) {
            v.outgoing = outgoing
                .get(&vh)
                .copied()
                .unwrap_or(HalfEdgeHandle::default());
        }
    }

    // Remove orphan vertices
    let orphans: Vec<VertexHandle> = brep
        .vertices
        .iter()
        .filter(|(_, v)| v.outgoing == HalfEdgeHandle::default())
        .map(|(h, _)| h)
        .collect();
    for vh in orphans {
        brep.vertices.remove(vh);
    }
}

/// Partition the assembled faces into edge-connected components and
/// register one shell (and solid) per component.
///
/// A boolean can DISCONNECT the workpiece — a flush through-slot
/// subtract leaves two closed boxes — and an assembly that files
/// every face under one shell fails Euler-Poincaré by construction
/// (χ = 2 per closed component, the formula expects 2 per SHELL), so
/// the whole cellular result fell back to graph-merge (measured on
/// the slotted-block corpus part: correct 12-face fragment population
/// discarded, graph-merge served 48 open edges).
///
/// Mutually exterior components each become their own solid. A
/// component whose AABB nests inside another's is classified by ray
/// parity against the CONTAINING component's face subset (PK-217):
/// an interior point of the nested component's enclosed region that
/// reads Inside the container is a CAVITY (an inner shell of the
/// container's solid — a moat's void, a closed hollow); one that
/// reads Outside is a FLOATING PIECE in a through-opening (its own
/// solid). Parity runs on `point_in_solid_final` with the face-AABB
/// slice filtered to the container — the slice IS the face filter.
/// Any ambiguity (no interior point found, OnBoundary parity, ray
/// failure, cavity whose encloser is itself a cavity) DECLINES the
/// partition (returns false, single-shell result stands, the Euler
/// gate routes it to graph-merge exactly as before) with the reason
/// observable under the boolean trace topic.
pub(super) fn partition_shells(
    builder: &mut BRepBuilder,
    tolerance: &popkern_core::tolerance::Tolerance,
) -> bool {
    let tolerance_linear = tolerance.linear();
    let brep = builder.brep_mut();
    let faces: Vec<_> = brep.faces.keys().collect();
    if faces.len() < 2 {
        return true;
    }
    let index: HashMap<_, usize> = faces.iter().enumerate().map(|(i, &f)| (f, i)).collect();

    // Union-find over faces joined by a shared (twinned) edge.
    let mut parent: Vec<usize> = (0..faces.len()).collect();
    fn find(parent: &mut [usize], i: usize) -> usize {
        let mut r = i;
        while parent[r] != r {
            r = parent[r];
        }
        let mut c = i;
        while parent[c] != r {
            let next = parent[c];
            parent[c] = r;
            c = next;
        }
        r
    }
    for (h, he) in brep.half_edges.iter() {
        let twin = he.twin;
        if twin == h || !brep.half_edges.contains_key(twin) {
            continue;
        }
        let (Some(l1), Some(l2)) = (
            brep.loops.get(he.loop_),
            brep.loops.get(brep.half_edges[twin].loop_),
        ) else {
            continue;
        };
        if let (Some(&i1), Some(&i2)) = (index.get(&l1.face), index.get(&l2.face)) {
            let (r1, r2) = (find(&mut parent, i1), find(&mut parent, i2));
            if r1 != r2 {
                parent[r1] = r2;
            }
        }
    }
    let mut component_of: HashMap<usize, usize> = HashMap::new();
    let mut roots: Vec<usize> = Vec::new();
    for i in 0..faces.len() {
        let r = find(&mut parent, i);
        let c = *component_of.entry(r).or_insert_with(|| {
            roots.push(r);
            roots.len() - 1
        });
        let _ = c;
    }
    if roots.len() < 2 {
        return true;
    }

    // Component AABBs from wired half-edge origins.
    let mut aabbs: Vec<([f64; 3], [f64; 3])> = vec![([f64::MAX; 3], [f64::MIN; 3]); roots.len()];
    // Vertex sums per component (each wired origin once per outgoing
    // half-edge — a weighting that keeps the mean inside the hull).
    let mut vertex_sums: Vec<([f64; 3], f64)> = vec![([0.0; 3], 0.0); roots.len()];
    for (_, he) in brep.half_edges.iter() {
        let Some(lp) = brep.loops.get(he.loop_) else {
            continue;
        };
        let Some(&fi) = index.get(&lp.face) else {
            continue;
        };
        let r = find(&mut parent, fi);
        let c = component_of[&r];
        let p = brep.vertices[he.origin].position;
        for (k, v) in [p.x, p.y, p.z].into_iter().enumerate() {
            aabbs[c].0[k] = aabbs[c].0[k].min(v);
            aabbs[c].1[k] = aabbs[c].1[k].max(v);
            vertex_sums[c].0[k] += v;
        }
        vertex_sums[c].1 += 1.0;
    }
    let tol = tolerance_linear.max(1e-12);
    let nested = |a: &([f64; 3], [f64; 3]), b: &([f64; 3], [f64; 3])| -> bool {
        (0..3).all(|k| a.0[k] > b.0[k] - tol && a.1[k] < b.1[k] + tol)
    };
    let face_component: Vec<(popkern_core::topology::FaceHandle, usize)> = faces
        .iter()
        .enumerate()
        .map(|(i, &f)| (f, component_of[&find(&mut parent, i)]))
        .collect();

    // Enclosure classification (PK-217), needed only when AABBs nest.
    // `encloser[c]` = the component whose solid a cavity shell joins.
    let any_nesting =
        (0..roots.len()).any(|i| (0..roots.len()).any(|j| i != j && nested(&aabbs[i], &aabbs[j])));
    let mut encloser: Vec<Option<usize>> = vec![None; roots.len()];
    if any_nesting {
        use crate::classify::point_in_solid_final;
        use crate::error::PointClassification;
        let all_aabbs = crate::intersect::compute_face_aabbs(brep);
        let comp_of_face: HashMap<popkern_core::topology::FaceHandle, usize> =
            face_component.iter().copied().collect();
        let comp_aabbs = |c: usize| -> Vec<_> {
            all_aabbs
                .iter()
                .filter(|(f, _)| comp_of_face.get(f) == Some(&c))
                .cloned()
                .collect()
        };
        // Interior point of each nested component's ENCLOSED region:
        // the AABB center, else the vertex centroid, whichever first
        // reads Inside against the component's own faces. The centroid
        // is a convex combination of the vertices, so it is interior
        // to every convex piece — a corner tetrahedron's box center is
        // not (PK-705: a yawed slab cut left a corner lump of a plate
        // nested in the plate's box; the center missed it and the
        // partition declined). No candidate → decline (a nested
        // ring-shaped component is out of scope, recorded).
        let interior_point = |c: usize| -> Option<popkern_core::geometry::Point3> {
            let (lo, hi) = &aabbs[c];
            let center = popkern_core::geometry::Point3::new(
                0.5 * (lo[0] + hi[0]),
                0.5 * (lo[1] + hi[1]),
                0.5 * (lo[2] + hi[2]),
            );
            let (sum, n) = vertex_sums[c];
            let centroid = (n > 0.0)
                .then(|| popkern_core::geometry::Point3::new(sum[0] / n, sum[1] / n, sum[2] / n));
            let own = comp_aabbs(c);
            // The centre of the component's FACE boxes: a cavity bounded
            // by a cylinder and two discs has only seam vertices, so its
            // vertex box collapses onto the seam line and both vertex
            // candidates sit on its boundary (PK-745: a plug sealing a
            // cylindrical cup declined as an Euler violation).
            let face_center = own
                .iter()
                .map(|(_, bb)| (bb.min, bb.max))
                .reduce(|(a0, a1), (b0, b1)| {
                    (
                        popkern_core::geometry::Point3::new(
                            a0.x.min(b0.x),
                            a0.y.min(b0.y),
                            a0.z.min(b0.z),
                        ),
                        popkern_core::geometry::Point3::new(
                            a1.x.max(b1.x),
                            a1.y.max(b1.y),
                            a1.z.max(b1.z),
                        ),
                    )
                })
                .map(|(lo, hi)| {
                    popkern_core::geometry::Point3::new(
                        0.5 * (lo.x + hi.x),
                        0.5 * (lo.y + hi.y),
                        0.5 * (lo.z + hi.z),
                    )
                });
            std::iter::once(center)
                .chain(centroid)
                .chain(face_center)
                .find(|q| {
                    matches!(
                        point_in_solid_final(q, brep, &own, tolerance),
                        Ok(PointClassification::Inside)
                    )
                })
        };
        // Enclosure depth by parity against each AABB-containing
        // component's face subset; immediate encloser = the smallest
        // containing component by AABB volume.
        let volume = |c: usize| -> f64 {
            let (lo, hi) = &aabbs[c];
            (hi[0] - lo[0]) * (hi[1] - lo[1]) * (hi[2] - lo[2])
        };
        let mut depth = vec![0usize; roots.len()];
        for i in 0..roots.len() {
            let containers: Vec<usize> = (0..roots.len())
                .filter(|&j| j != i && nested(&aabbs[i], &aabbs[j]))
                .collect();
            if containers.is_empty() {
                continue;
            }
            let Some(q) = interior_point(i) else {
                if popkern_core::trace::trace_enabled("boolean") {
                    eprintln!(
                        "  [assemble] shell partition DECLINED: no interior point for nested component {i}"
                    );
                }
                return false;
            };
            let mut inside: Vec<usize> = Vec::new();
            for j in containers {
                match point_in_solid_final(&q, brep, &comp_aabbs(j), tolerance) {
                    Ok(PointClassification::Inside) => inside.push(j),
                    Ok(PointClassification::Outside) => {}
                    _ => {
                        if popkern_core::trace::trace_enabled("boolean") {
                            eprintln!(
                                "  [assemble] shell partition DECLINED: ambiguous parity of component {i} vs {j}"
                            );
                        }
                        return false;
                    }
                }
            }
            depth[i] = inside.len();
            if depth[i] % 2 == 1 {
                encloser[i] = inside.into_iter().min_by(|&a, &b| {
                    volume(a)
                        .partial_cmp(&volume(b))
                        .unwrap_or(std::cmp::Ordering::Equal)
                });
            }
        }
        // A cavity's encloser must itself be a solid-owning component
        // (even depth) — alternation anomalies decline.
        for (i, enc) in encloser.iter().enumerate() {
            if let Some(e) = *enc
                && depth[e] % 2 == 1
            {
                if popkern_core::trace::trace_enabled("boolean") {
                    eprintln!(
                        "  [assemble] shell partition DECLINED: cavity {i} enclosed by cavity {e}"
                    );
                }
                return false;
            }
        }
    }

    // Wire shells: every component gets a shell; solid-owning
    // components (even enclosure depth) get their own solid, cavity
    // components join their encloser's solid. The first solid-owning
    // component reuses the existing shell/solid.
    let first_own = (0..roots.len()).find(|&c| encloser[c].is_none());
    let Some(first_own) = first_own else {
        return false; // every component claims to be a cavity — lie
    };
    let existing_shell = brep.faces[faces[0]].shell;
    let existing_solid = brep.shells[existing_shell].solid;
    let mut shells: Vec<Option<popkern_core::topology::ShellHandle>> = vec![None; roots.len()];
    let mut solids: Vec<Option<popkern_core::topology::SolidHandle>> = vec![None; roots.len()];
    shells[first_own] = Some(existing_shell);
    solids[first_own] = Some(existing_solid);
    for c in 0..roots.len() {
        if c == first_own || encloser[c].is_some() {
            continue;
        }
        let solid = builder.add_solid();
        solids[c] = Some(solid);
        shells[c] = Some(builder.add_shell(solid));
    }
    let mut cavity_count = 0usize;
    for c in 0..roots.len() {
        if let Some(e) = encloser[c] {
            let Some(owner_solid) = solids[e] else {
                return false;
            };
            shells[c] = Some(builder.add_shell(owner_solid));
            cavity_count += 1;
        }
    }
    // `add_shell` re-points its solid's primary shell at the newest
    // one; cavity shells must not usurp the outer shell's slot.
    let brep = builder.brep_mut();
    for c in 0..roots.len() {
        if encloser[c].is_none()
            && let (Some(sh), Some(so)) = (shells[c], solids[c])
        {
            brep.solids[so].shell = sh;
        }
    }
    // Re-anchor every shell's representative face (component 0's may
    // have pointed at a face now moved to another shell).
    let mut anchored = vec![false; roots.len()];
    for (f, c) in face_component {
        let Some(sh) = shells[c] else { return false };
        brep.faces[f].shell = sh;
        if !anchored[c] {
            brep.shells[sh].face = f;
            anchored[c] = true;
        }
    }
    if popkern_core::trace::trace_enabled("boolean") {
        eprintln!(
            "  [assemble] shell partition: {} component(s) → {} solid(s), {} cavity shell(s)",
            roots.len(),
            roots.len() - cavity_count,
            cavity_count
        );
    }
    true
}

/// Build entity provenance map from assembly information.
///
/// Records which input face each output face came from, and which
/// graph edge source each output vertex came from.
/// Build entity provenance from assembly data.
///
/// Populates face provenance by matching output faces to their source
/// fragments. Called after assembly when provenance tracking is needed.
pub(crate) fn build_face_provenance(
    result: &BRep,
    kept: &[usize],
    fragments: &[Fragment],
    _selections: &[FragmentSelection],
) -> crate::provenance::EntityMap {
    use crate::provenance::types::{EntityMap, FaceProvenance, SourceBody};

    let mut map = EntityMap::default();

    let result_faces: Vec<_> = result.faces.keys().collect();

    // Match output faces to input faces via fragment ordering.
    // The assembly creates faces in the same order as `kept` iteration.
    for (i, &fi) in kept.iter().enumerate() {
        if i >= result_faces.len() {
            break;
        }
        let frag = &fragments[fi];
        let source = match frag.face.operand {
            Operand::A => SourceBody::A,
            Operand::B => SourceBody::B,
        };
        map.faces.push((
            result_faces[i],
            FaceProvenance::Kept {
                source,
                original_face: frag.face.face,
            },
        ));
    }

    map
}
