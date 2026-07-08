//! Order-independent BRep assembly from selected graph fragments.
//!
//! Unlike `merge/assemble.rs`, this assembly pre-creates all edges from
//! a canonical set (sorted by node indices) before wiring loops. This
//! makes the output independent of fragment iteration order, enabling
//! BVH-based face-pair filtering.

use std::collections::HashMap;

use popkern_core::{
    constructors::builder::BRepBuilder,
    geometry::curve::{CurveGeometry, CurveOnEdge},
    tolerance::Tolerance,
    topology::{BRep, HalfEdgeHandle, VertexHandle},
};

use crate::error::{BoolDiagnostics, BoolError, Stage};
use crate::merge::fragment_classify::FragmentSelection;
use crate::merge::graph::{IntersectionGraph, Operand};
use crate::merge::walk::Fragment;

/// Assemble a BRep from selected fragments using order-independent edge creation.
///
/// 1. Collect all vertices needed by kept fragments
/// 2. Collect all unique edges (canonical: min_node, max_node)
/// 3. Create all vertices and edges upfront
/// 4. Wire fragment loops using the pre-created edges
/// 5. Attach surface geometry
pub(crate) fn cellular_assemble(
    graph: &IntersectionGraph,
    fragments: &[Fragment],
    selections: &[FragmentSelection],
    brep_a: &BRep,
    brep_b: &BRep,
    tolerance: Tolerance,
) -> Result<BRep, BoolDiagnostics> {
    let mut builder = BRepBuilder::new(tolerance);
    let solid = builder.add_solid();
    let shell = builder.add_shell(solid);

    // Step 1: Collect all nodes needed by kept fragments.
    let mut needed_nodes = std::collections::BTreeSet::new();
    let kept: Vec<usize> = (0..fragments.len())
        .filter(|&i| selections[i].keep)
        .collect();

    for &fi in &kept {
        let sel = &selections[fi];
        let frag = &fragments[fi];
        for &ghe in frag
            .outer_loop
            .iter()
            .chain(frag.inner_loops.iter().flat_map(|l| l.iter()))
        {
            let ghe = if sel.flip {
                graph.half_edges[ghe].twin
            } else {
                ghe
            };
            needed_nodes.insert(graph.half_edges[ghe].origin);
            let twin = graph.half_edges[ghe].twin;
            needed_nodes.insert(graph.half_edges[twin].origin);
        }
    }

    // Step 2: Create vertices (sorted by node index for determinism).
    let mut node_to_vertex: HashMap<usize, VertexHandle> = HashMap::new();
    for &node in &needed_nodes {
        let vh = builder.add_vertex(graph.nodes[node].position);
        node_to_vertex.insert(node, vh);
    }

    // Step 3: Collect all unique edges. CANONICAL KEY = the graph
    // edge's own identity, min(he, twin) — NOT the node pair: two
    // distinct edges can connect the same nodes (the two arcs of a
    // split intersection circle) and a closed intersection curve is
    // a self-loop whose node pair (n, n) cannot distinguish it from
    // anything else at that node. Node-pair keying merged those,
    // dropping curves and mispairing faces.
    let mut canonical_edges: std::collections::BTreeMap<usize, usize> =
        std::collections::BTreeMap::new();

    for &fi in &kept {
        let sel = &selections[fi];
        let frag = &fragments[fi];
        for &raw_ghe in frag
            .outer_loop
            .iter()
            .chain(frag.inner_loops.iter().flat_map(|l| l.iter()))
        {
            let ghe = if sel.flip {
                graph.half_edges[raw_ghe].twin
            } else {
                raw_ghe
            };
            let twin = graph.half_edges[ghe].twin;
            let key = ghe.min(twin);
            canonical_edges.entry(key).or_insert(key);
        }
    }

    // Step 4: Create all edges upfront (deterministic order from
    // BTreeMap). Map: min-half-edge → (he_fwd, he_twin) where fwd
    // follows the MIN graph half-edge's direction.
    let mut edge_he_map: HashMap<usize, (HalfEdgeHandle, HalfEdgeHandle)> = HashMap::new();

    // DEGENERATE edges — coincident endpoint nodes with no curve on
    // either half (a genuine closed rim carries its circle) — are
    // boolean residue at discarded tool corners. Creating them left
    // zero-length curveless edges in the output (χ' drifted to 1 and
    // the mesher choked); skip creation and splice them out of the
    // fragment loops below.
    let is_degenerate_key = |key: usize| -> bool {
        let na = graph.half_edges[key].origin;
        let nb = graph.half_edges[graph.half_edges[key].twin].origin;
        if graph.nodes[na]
            .position
            .distance_to(&graph.nodes[nb].position)
            >= 1.0e-9
        {
            return false;
        }
        // Coincident endpoints: a genuine closed rim (circle) has a
        // curve with real extent — sample its midpoint. A degenerate
        // stub's curve stays at the point.
        let curve = graph.half_edges[key]
            .curve
            .as_ref()
            .or_else(|| graph.half_edges[graph.half_edges[key].twin].curve.as_ref());
        match curve {
            None => true,
            Some(c) => {
                use popkern_core::geometry::curve::CurveEval;
                let pm = c.geometry.point_at(0.5 * (c.t_start + c.t_end));
                pm.distance_to(&graph.nodes[na].position) < 1.0e-9
            }
        }
    };
    for (&key, &sample_ghe) in &canonical_edges {
        if is_degenerate_key(key) {
            continue;
        }
        let na = graph.half_edges[key].origin;
        let nb = graph.half_edges[graph.half_edges[key].twin].origin;
        let vh_a = node_to_vertex[&na];
        let vh_b = node_to_vertex[&nb];
        let (edge_h, he_fwd, he_twin) = builder.add_edge(vh_a, vh_b);

        edge_he_map.insert(key, (he_fwd, he_twin));
        if std::env::var_os("PK_BOOL_TRACE").is_some() {
            eprintln!(
                "  [assemble] edge key {key}: nodes ({na},{nb}) shared={} curve={}",
                matches!(
                    graph.half_edges[sample_ghe].source,
                    crate::merge::graph::EdgeSource::Shared { .. }
                ),
                graph.half_edges[sample_ghe].curve.is_some()
                    || graph.half_edges[graph.half_edges[sample_ghe].twin]
                        .curve
                        .is_some()
            );
        }

        // Attach curve geometry.
        // For shared (intersection) edges: use the graph's actual intersection
        // curve. Both incident faces share the exact same CurveOnEdge, enforcing
        // the shared geometry principle (spec §4.5, §9.2) for watertightness.
        // For original/split edges and self-loops: use the graph's stored curve
        // if available, otherwise synthesize a line from endpoints.
        let is_shared = matches!(
            graph.half_edges[sample_ghe].source,
            crate::merge::graph::EdgeSource::Shared { .. }
        );
        if is_shared || na == nb {
            // Shared edge or self-loop: use the graph's actual intersection
            // curve for watertightness (spec §4.5, §9.2).
            //
            // Determine which graph half-edge owns the curve and align the
            // parameterization to match the canonical edge direction (na→nb).
            let (curve_he, curve_ref) = if graph.half_edges[sample_ghe].curve.is_some() {
                (sample_ghe, graph.half_edges[sample_ghe].curve.as_ref())
            } else {
                let twin = graph.half_edges[sample_ghe].twin;
                (twin, graph.half_edges[twin].curve.as_ref())
            };

            if let Some(curve) = curve_ref {
                let curve_origin = graph.half_edges[curve_he].origin;
                let mut aligned = curve.clone();
                if curve_origin != na && na != nb {
                    // Curve goes nb→na but we need na→nb: swap parameters
                    std::mem::swap(&mut aligned.t_start, &mut aligned.t_end);
                }
                if std::env::var_os("PK_BOOL_TRACE").is_some() {
                    eprintln!(
                        "  [assemble] set curve on key {sample_ghe}: t [{:.3},{:.3}]",
                        aligned.t_start, aligned.t_end
                    );
                }
                builder.set_edge_curve(edge_h, aligned);
            } else if std::env::var_os("PK_BOOL_TRACE").is_some() {
                eprintln!(
                    "  [assemble] shared/self-loop edge key {sample_ghe} ({na},{nb}): NO CURVE on either half"
                );
            }
        } else {
            // Original/split edge: ADOPT the stored curve when its
            // endpoints coincide with this edge's nodes — an
            // untouched curved boundary edge synthesized as a line
            // degraded arcs to their chords (a cylinder wall's
            // bottom arc and the wedge ellipse both flattened, the
            // chord then generated phantom edge-face hits and the
            // mesher's boundary-first pass bailed on the off-surface
            // curve). Fragments of a SPLIT curved edge whose stored
            // curve spans the whole original fail the endpoint match
            // and keep the line fallback.
            let p_start = graph.nodes[na].position;
            let p_end = graph.nodes[nb].position;
            let stored = graph.half_edges[sample_ghe].curve.as_ref().or_else(|| {
                graph.half_edges[graph.half_edges[sample_ghe].twin]
                    .curve
                    .as_ref()
            });
            let adopted = stored.and_then(|c| {
                use popkern_core::geometry::curve::CurveEval;
                let c0 = c.geometry.point_at(c.t_start);
                let c1 = c.geometry.point_at(c.t_end);
                let tol = 1.0e-6;
                if c0.distance_to(&p_start) < tol && c1.distance_to(&p_end) < tol {
                    Some(c.clone())
                } else if c0.distance_to(&p_end) < tol && c1.distance_to(&p_start) < tol {
                    let mut rev = c.clone();
                    std::mem::swap(&mut rev.t_start, &mut rev.t_end);
                    Some(rev)
                } else {
                    None
                }
            });
            if let Some(curve) = adopted {
                builder.set_edge_curve(edge_h, curve);
            } else {
                let dir = p_end - p_start;
                let len = dir.length();
                if len > 1e-12 {
                    let line = CurveOnEdge {
                        geometry: CurveGeometry::Line {
                            origin: p_start,
                            direction: dir * (1.0 / len),
                        },
                        t_start: 0.0,
                        t_end: len,
                    };
                    builder.set_edge_curve(edge_h, line);
                }
            }
        }
    }

    // Step 5: Wire fragment loops.
    // Track which half-edges have been used (to handle shared edges).
    // When both fwd and twin are used, create a duplicate edge.
    let mut used_hes: std::collections::HashSet<HalfEdgeHandle> = std::collections::HashSet::new();

    for &fi in &kept {
        let sel = &selections[fi];
        let frag = &fragments[fi];
        if std::env::var_os("PK_BOOL_TRACE").is_some() {
            eprintln!(
                "  [assemble] kept fragment {:?} outer {} inner {} flip={}",
                frag.face,
                frag.outer_loop.len(),
                frag.inner_loops.len(),
                sel.flip
            );
        }
        let face_h = builder.add_face(shell);

        // Resolve all graph half-edges to BRep half-edges upfront,
        // before calling builder methods that borrow mutably.
        let mut resolve_loop = |raw_hes: &[usize],
                                flip: bool,
                                builder: &mut BRepBuilder|
         -> Vec<HalfEdgeHandle> {
            // Splice out degenerate zero-length hes (their edges were
            // never created).
            let filtered: Vec<usize> = raw_hes
                .iter()
                .copied()
                .filter(|&ghe| {
                    let twin = graph.half_edges[ghe].twin;
                    !is_degenerate_key(ghe.min(twin))
                })
                .collect();
            let raw_hes: &[usize] = &filtered;
            // Dry-run: would the straight mapping conflict with
            // already-wired half-edges? A loop whose directions were
            // consumed by other fragments (the walker may hand a face
            // the twin-direction ring wholesale) is re-wired REVERSED
            // via the twin directions — the same ring with flipped
            // winding, which any loop tolerates. Substituting twins
            // per-half-edge instead breaks the chain (reversed
            // endpoints), and minting duplicates detaches the ring
            // (measured: a star socket in a counterbore floor went
            // non-manifold on all 12 minted edges).
            let straight_conflicts = raw_hes.iter().any(|&raw_ghe| {
                let ghe = if flip {
                    graph.half_edges[raw_ghe].twin
                } else {
                    raw_ghe
                };
                let twin_idx = graph.half_edges[ghe].twin;
                let key = ghe.min(twin_idx);
                match edge_he_map.get(&key) {
                    Some(&(he_fwd, he_twin)) => {
                        let brep_he = if ghe == key { he_fwd } else { he_twin };
                        used_hes.contains(&brep_he)
                    }
                    None => false,
                }
            });
            let reversed_free = raw_hes.iter().all(|&raw_ghe| {
                let ghe = if flip {
                    graph.half_edges[raw_ghe].twin
                } else {
                    raw_ghe
                };
                // Twin direction of this he.
                let ghe = graph.half_edges[ghe].twin;
                let twin_idx = graph.half_edges[ghe].twin;
                let key = ghe.min(twin_idx);
                match edge_he_map.get(&key) {
                    Some(&(he_fwd, he_twin)) => {
                        let brep_he = if ghe == key { he_fwd } else { he_twin };
                        !used_hes.contains(&brep_he)
                    }
                    None => false,
                }
            });
            let (order, extra_twin): (Vec<usize>, bool) = if straight_conflicts && reversed_free {
                if std::env::var_os("PK_BOOL_TRACE").is_some() {
                    eprintln!(
                        "  [resolve] loop of {} re-wired reversed (direction conflict)",
                        raw_hes.len()
                    );
                }
                (raw_hes.iter().rev().copied().collect(), true)
            } else {
                (raw_hes.to_vec(), false)
            };
            let mut result = Vec::new();
            for &raw_ghe in &order {
                let raw_ghe = if extra_twin {
                    graph.half_edges[raw_ghe].twin
                } else {
                    raw_ghe
                };
                let ghe = if flip {
                    graph.half_edges[raw_ghe].twin
                } else {
                    raw_ghe
                };
                let twin_idx = graph.half_edges[ghe].twin;
                let key = ghe.min(twin_idx);

                let Some(&(he_fwd, he_twin)) = edge_he_map.get(&key) else {
                    if std::env::var_os("PK_BOOL_TRACE").is_some() {
                        eprintln!(
                            "  [resolve] ghe {ghe} (key {key}): edge NOT IN MAP — dropped from loop"
                        );
                    }
                    continue;
                };

                // fwd follows the MIN half-edge's direction; identity
                // (not node comparison) resolves self-loops and
                // double edges correctly.
                let brep_he = if ghe == key { he_fwd } else { he_twin };

                if !used_hes.contains(&brep_he) {
                    used_hes.insert(brep_he);
                    result.push(brep_he);
                } else {
                    // The preferred direction is taken by another
                    // fragment's loop. The twin CANNOT substitute for
                    // an open-chain half-edge — its endpoints are
                    // reversed, so the wired loop breaks exactly there
                    // (measured: a star socket cut into a counterbore
                    // floor lost one ring segment to the swap and the
                    // inner loop failed closure by one star chord).
                    // Twin-substitution is only sound for self-loops
                    // (start == end); everything else mints a fresh
                    // edge pair with the correct orientation.
                    let twin_brep = if brep_he == he_fwd { he_twin } else { he_fwd };
                    let ghe_data = &graph.half_edges[ghe];
                    let is_self_loop = ghe_data.origin == graph.half_edges[ghe_data.twin].origin;
                    if is_self_loop && !used_hes.contains(&twin_brep) {
                        used_hes.insert(twin_brep);
                        result.push(twin_brep);
                    } else {
                        let origin_node = ghe_data.origin;
                        let end_node = graph.half_edges[ghe_data.twin].origin;
                        if let (Some(&vh_start), Some(&vh_end)) = (
                            node_to_vertex.get(&origin_node),
                            node_to_vertex.get(&end_node),
                        ) {
                            let (_dup_edge, dup_fwd, dup_twin) = builder.add_edge(vh_start, vh_end);
                            let _ = dup_twin;
                            used_hes.insert(dup_fwd);
                            result.push(dup_fwd);
                            if std::env::var_os("PK_BOOL_TRACE").is_some() {
                                eprintln!(
                                    "  [resolve] ghe {ghe} (key {key}): direction conflict — minted duplicate edge"
                                );
                            }
                        } else if std::env::var_os("PK_BOOL_TRACE").is_some() {
                            eprintln!(
                                "  [resolve] ghe {ghe} (key {key}): conflict and nodes unmapped — dropped"
                            );
                        }
                    }
                }
            }
            result
        };

        let raw_outer: Vec<usize> = if sel.flip {
            frag.outer_loop.iter().rev().copied().collect()
        } else {
            frag.outer_loop.clone()
        };
        let outer_hes = resolve_loop(&raw_outer, sel.flip, &mut builder);

        if !outer_hes.is_empty() {
            builder.add_loop(&outer_hes, face_h).map_err(|e| {
                BoolDiagnostics::from_error(
                    BoolError::MergeFailed {
                        stage: format!("{:?}", Stage::Merge),
                        reason: format!("cellular outer loop: {e}"),
                    },
                    Stage::Merge,
                )
            })?;
        }

        for inner in &frag.inner_loops {
            let raw_inner: Vec<usize> = if sel.flip {
                inner.iter().rev().copied().collect()
            } else {
                inner.clone()
            };
            let inner_hes = resolve_loop(&raw_inner, sel.flip, &mut builder);
            if !inner_hes.is_empty() {
                if let Err(e) = builder.add_inner_loop(&inner_hes, face_h)
                    && std::env::var_os("PK_BOOL_TRACE").is_some()
                {
                    eprintln!(
                        "  [assemble] add_inner_loop FAILED on {face_h:?} ({} hes): {e}",
                        inner_hes.len()
                    );
                }
            } else if std::env::var_os("PK_BOOL_TRACE").is_some() {
                eprintln!("  [assemble] inner loop resolved EMPTY on {face_h:?}");
            }
        }

        // Attach surface geometry
        let source_brep = match frag.face.operand {
            Operand::A => brep_a,
            Operand::B => brep_b,
        };
        if let Some(face) = source_brep.faces.get(frag.face.face)
            && let Some(ref surface) = face.surface
        {
            let mut surf = surface.clone();
            if sel.flip {
                surf.orientation = !surf.orientation;
            }
            builder.set_face_surface(face_h, surf);
        }
    }

    // Cleanup: remove unwired half-edges, dangling edges, orphan vertices.
    cleanup_unwired(&mut builder);

    Ok(builder.build_unchecked())
}

/// Remove half-edges not wired into any loop, and their orphaned edges/vertices.
fn cleanup_unwired(builder: &mut BRepBuilder) {
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
        brep.half_edges.remove(*h);
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
