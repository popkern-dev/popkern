//! Boundary-driven blend-ring carve for saddle (`OppositeSide`)
//! corners (PK-234).
//!
//! The mitred corner construction leaves the saddle's blend ring on
//! the SHARP model edges: rim chords claim fillet cylinders they sit
//! up to `√2·r − r` away from, while the three-quad decomposition's
//! Gregory patches interpolate the mitred [`BlendBoundary`] 1.3e-1
//! away (measured, notched-cube saddle). The corner-point setback
//! surgery cannot serve this class either — at the reflex junction
//! between the two convex strips the rims never meet, so its
//! 3-corner ring model has no valid corner (measured: the
//! tangent-line crossing sits outside both strips' physical
//! extents).
//!
//! This carve builds the ring the boundary actually describes:
//!
//! * per `FilletEnd` segment — split the strip face's boundary at
//!   the arc's two junction endpoints (they land on the strip's
//!   tangency rails by construction) and `mef` the exact arc curve
//!   across the strip;
//! * per `FaceEdge` segment — `mef` the bridge curve across its
//!   support face between the two junction vertices the arcs
//!   created;
//! * merge the surrendered region (apex slivers + wedges inside the
//!   ring) into a single `2n`-gon blend face via the setback
//!   surgery's region merge.
//!
//! The blend face then flows through the UNCHANGED three-quad
//! dispatch: the poke sees a `2n` ring, the seam re-anchor measures
//! already-anchored (the ring IS the boundary), and the decomposition
//! bonds onto edges that lie on their faces.
//!
//! Transactional: any step error restores the entry snapshot and the
//! caller falls through to the mitred fallback construction with the
//! decline recorded (today's behavior, no regression).

use super::boundary_ring_geom::{
    junction_vertex, select_bridge_support, worst_ring_on_surface_deviation, write_ring_curve,
};
use super::*;

/// Junction acceptance: a boundary junction must land on a strip
/// boundary edge within this band before we split there. Junctions
/// are constructed ON the tangency rails, and the BRep rail curves
/// are the same carve output the boundary read, so genuine landings
/// measure ~1e-9; a 1e-4 gate rejects a boundary/topology mismatch
/// class without admitting mis-splits.
pub(super) const JUNCTION_ON_EDGE_BAND: f64 = 1.0e-4;

/// Ring postcondition: every ring edge curve sample must sit on both
/// adjacent surfaces within this band (the arcs are exact chart isos
/// and the bridges live on their support faces, so genuine rings
/// measure ≤ 1e-9; check 9's default band is 1e-5).
pub(super) const RING_ON_SURFACE_BAND: f64 = 1.0e-5;

/// Carve the boundary ring at a saddle vertex and merge the interior
/// into one `2n`-gon blend face. Returns the blend face; on any
/// decline the BRep is exactly the entry state.
pub(crate) fn try_boundary_ring_carve(
    brep: &mut BRep,
    vertex: VertexHandle,
    boundary: &crate::blend::boundary::BlendBoundary,
    fillet_surfaces: &[FilletSurfaceResult],
    tolerance: &Tolerance,
) -> Result<FaceHandle, String> {
    use crate::blend::boundary::BoundarySegment;

    let n = fillet_surfaces.len();
    if boundary.segments.len() != 2 * n {
        return Err(format!(
            "boundary has {} segments, expected {}",
            boundary.segments.len(),
            2 * n
        ));
    }
    let vertex_position = brep
        .vertices
        .get(vertex)
        .map(|v| v.position)
        .ok_or_else(|| "saddle vertex not in BRep".to_string())?;

    // Segment pattern: alternating FilletEnd / FaceEdge starting at a
    // FilletEnd (the boundary builder's contract).
    for (i, seg) in boundary.segments.iter().enumerate() {
        let ok = if i % 2 == 0 {
            matches!(seg, BoundarySegment::FilletEnd { .. })
        } else {
            matches!(seg, BoundarySegment::FaceEdge { .. })
        };
        if !ok {
            return Err(format!("segment {i} breaks the arc/bridge alternation"));
        }
    }

    let sides: Vec<RingSide<'_>> = boundary
        .segments
        .iter()
        .map(|seg| match seg {
            BoundarySegment::FilletEnd {
                fillet_index,
                curve,
            } => RingSide::Arc {
                fillet_index: *fillet_index,
                curve,
            },
            BoundarySegment::FaceEdge { curve, .. } => RingSide::Bridge { curve },
        })
        .collect();
    let snapshot = brep.clone();
    let carved = carve_ring_sides(
        brep,
        vertex,
        vertex_position,
        &sides,
        fillet_surfaces,
        tolerance,
    )
    .and_then(|ring| {
        super::boundary_ring_geom::normalize_alternating_ring(brep, &ring, n)?;
        Ok(ring.blend_face)
    });
    if carved.is_err() {
        *brep = snapshot;
    }
    carved
}

/// One side of a carved ring.
pub(crate) enum RingSide<'a> {
    /// A strip-end arc, carved across its strip face.
    Arc {
        /// Index into the fillet surfaces.
        fillet_index: usize,
        /// The exact arc.
        curve: &'a NurbsCurve,
    },
    /// A bridge, carved across the support face its junctions share.
    Bridge {
        /// The exact bridge curve.
        curve: &'a NurbsCurve,
    },
}

/// A carved ring: the merged blend face and each side's edge, in side
/// order.
pub(crate) struct CarvedRing {
    /// The merged blend face bounded by the ring.
    pub blend_face: FaceHandle,
    /// Per side, its ring edge.
    pub side_edges: Vec<EdgeHandle>,
    /// Per side, whether it is an arc.
    pub side_is_arc: Vec<bool>,
}

/// Carve a ring given as sides in cyclic order and merge its interior
/// into one blend face. Consecutive arcs meet at a shared junction
/// with no bridge between them: the saddle setback's corner junctions
/// (PK-527). The blend face's ring has one edge per side. On error
/// the caller restores its snapshot.
pub(crate) fn carve_ring_sides(
    brep: &mut BRep,
    vertex: VertexHandle,
    vertex_position: Point3,
    sides: &[RingSide<'_>],
    fillet_surfaces: &[FilletSurfaceResult],
    tolerance: &Tolerance,
) -> Result<CarvedRing, String> {
    carve_ring(
        brep,
        vertex,
        vertex_position,
        sides,
        fillet_surfaces,
        tolerance,
    )
}

fn carve_ring(
    brep: &mut BRep,
    vertex: VertexHandle,
    vertex_position: Point3,
    sides: &[RingSide<'_>],
    fillet_surfaces: &[FilletSurfaceResult],
    tolerance: &Tolerance,
) -> Result<CarvedRing, String> {
    let side_count = sides.len();
    let trace = popkern_core::trace::trace_enabled("corner-carve");

    // ---- arcs: one exact strip-end curve per fillet -----------------
    let mut arc_vertices: Vec<Option<[VertexHandle; 2]>> = vec![None; side_count];
    let mut strip_faces: Vec<FaceHandle> = Vec::with_capacity(side_count);
    let mut side_edges: Vec<Option<EdgeHandle>> = vec![None; side_count];
    let mut ring_edges: std::collections::BTreeSet<EdgeHandle> = Default::default();
    let mut slivers: Vec<FaceHandle> = Vec::new();
    for (k, side) in sides.iter().enumerate() {
        let RingSide::Arc {
            fillet_index,
            curve,
        } = side
        else {
            continue;
        };
        let fillet = fillet_surfaces
            .get(*fillet_index)
            .ok_or_else(|| format!("fillet index {fillet_index} out of range"))?;
        let (t0, t1) = curve.parameter_range();
        let p_start = curve.point_at(t0);
        let p_end = curve.point_at(t1);
        // Bind the strip face by the arc's interior samples (the
        // binding cubic is geometric — interior points only).
        let pseudo_cubic = [
            p_start,
            curve.point_at(t0 + (t1 - t0) / 3.0),
            curve.point_at(t0 + 2.0 * (t1 - t0) / 3.0),
            p_end,
        ];
        let strip_face = setback_arc_binding::find_fillet_face(
            brep,
            fillet,
            &bezier_cubic_curve(&pseudo_cubic),
            vertex_position,
        )?;
        // Per endpoint: split the rail where it exists; where the
        // mitre complex replaced the rail with an off-surface chord
        // (junction interior to the chord-bounded trim), restore the
        // rail via the spur recut. The chord `kef` runs AFTER the arc
        // `mef` (the corner sliver must exist as a face first).
        let mut pending_chord_kefs: Vec<EdgeHandle> = Vec::new();
        let resolve = |brep: &mut BRep,
                       point: Point3,
                       pending: &mut Vec<EdgeHandle>|
         -> Result<VertexHandle, String> {
            match junction_vertex(brep, strip_face, point, tolerance)? {
                Some(v) => Ok(v),
                None => {
                    let recut = super::boundary_ring_recut::rail_spur_recut(
                        brep,
                        strip_face,
                        &fillet.surface,
                        point,
                    )
                    .map_err(|e| format!("rail-spur recut: {e}"))?;
                    pending.push(recut.chord_edge);
                    Ok(recut.vertex)
                }
            }
        };
        let v_start = resolve(brep, p_start, &mut pending_chord_kefs)
            .map_err(|e| format!("arc {k} start junction: {e}"))?;
        let v_end = resolve(brep, p_end, &mut pending_chord_kefs)
            .map_err(|e| format!("arc {k} end junction: {e}"))?;
        if v_start == v_end {
            return Err(format!("arc {k} junctions collapse to one vertex"));
        }
        let he_start = half_edge_at_vertex_on_face(brep, strip_face, v_start)?;
        let he_end = half_edge_at_vertex_on_face(brep, strip_face, v_end)?;
        let mef = brep
            .mef(he_start, he_end)
            .map_err(|e| format!("mef arc {k} across strip: {e}"))?;
        write_ring_curve(brep, mef.edge, curve, v_start, v_end);
        if pending_chord_kefs.is_empty() {
            let d_old = face_centroid(brep, strip_face)?.distance_to(&vertex_position);
            let d_new = face_centroid(brep, mef.face)?.distance_to(&vertex_position);
            slivers.push(if d_new < d_old { mef.face } else { strip_face });
        } else {
            // Recut class: return the corner sliver (overreach band +
            // between-station piece) to the support face by killing
            // the off-surface chord. The region merge's vertex-seeded
            // flood covers the interior afterwards.
            for chord in pending_chord_kefs.drain(..) {
                let he = brep
                    .edges
                    .get(chord)
                    .map(|e| e.half_edge)
                    .ok_or_else(|| format!("chord {chord:?} vanished before kef"))?;
                // The merged face's geometry must be the SUPPORT's
                // (the far side): the merge returns the overreach
                // sliver to the support plane. Keeping the sliver's
                // strip surface instead poisoned the whole support
                // face (measured: the entire cube face x=1 inherited
                // the sliver surface and every far cube edge failed
                // check 9 at 8.8e-1). The chord borders sliver (the
                // strip-split corner side, which owns the chord after
                // the arc mef) and support; identify the support as
                // the side that is NOT the arc-mef sliver — the face
                // NOT equal to either arc-split product.
                let f_a = brep.loops[brep.half_edges[he].loop_].face;
                let f_b = brep.loops[brep.half_edges[brep.half_edges[he].twin].loop_].face;
                let support_side = if f_a == strip_face || f_a == mef.face {
                    f_b
                } else {
                    f_a
                };
                let support_surface = brep.faces.get(support_side).and_then(|f| f.surface.clone());
                let kef = brep
                    .kef(he)
                    .map_err(|e| format!("kef mitre chord on arc {k}: {e}"))?;
                let survivor = brep.loops[kef.surviving_loop].face;
                if let Some(face_data) = brep.faces.get_mut(survivor) {
                    face_data.surface = support_surface;
                }
            }
        }
        ring_edges.insert(mef.edge);
        side_edges[k] = Some(mef.edge);
        arc_vertices[k] = Some([v_start, v_end]);
        strip_faces.push(strip_face);
        if trace {
            eprintln!(
                "  [boundary-ring] arc {k} strip {strip_face:?} \
                 ({:.4},{:.4},{:.4})->({:.4},{:.4},{:.4})",
                p_start.x, p_start.y, p_start.z, p_end.x, p_end.y, p_end.z
            );
        }
    }

    // ---- junctions: bridges, or arcs meeting at a shared vertex ------
    // A bridge connects the previous arc's END to the next arc's START
    // (the ring is in cyclic order). Consecutive arcs share their
    // junction vertex (the saddle setback's corners).
    let arc_end = |i: usize, end: usize| -> Result<VertexHandle, String> {
        arc_vertices[i]
            .map(|v| v[end])
            .ok_or_else(|| format!("side {i} is not an arc next to a junction"))
    };
    for k in 0..side_count {
        let next = (k + 1) % side_count;
        if matches!(sides[k], RingSide::Arc { .. }) && matches!(sides[next], RingSide::Arc { .. }) {
            if arc_end(k, 1)? != arc_end(next, 0)? {
                return Err(format!("arcs {k} and {next} do not meet at their junction"));
            }
            continue;
        }
        let RingSide::Bridge { curve } = sides[k] else {
            continue;
        };
        let v_a = arc_end((k + side_count - 1) % side_count, 1)?;
        let v_b = arc_end(next, 0)?;
        if v_a == v_b {
            // Degenerate bridge: the ring closes through the shared
            // junction. The blend face then has fewer than 2n edges
            // and the poke contract breaks — decline until a fixture
            // exercises this class.
            return Err(format!("bridge {k} is degenerate (shared junction vertex)"));
        }
        // Support face: the boundary's handle cannot discriminate
        // coplanar splits of one support plane (measured: it named a
        // mitre-wedge piece), and chord kefs move junction ownership.
        // Select topologically: the non-strip face carrying BOTH
        // junction vertices whose surface contains the bridge.
        let support = select_bridge_support(brep, v_a, v_b, &strip_faces, curve)
            .map_err(|e| format!("bridge {k}: {e}"))?;
        let (t0, t1) = curve.parameter_range();
        let b_start = curve.point_at(t0);
        let b_end = curve.point_at(t1);
        let pa = brep.vertices[v_a].position;
        // Bridge orientation in the boundary may run either way.
        let (v_from, v_to) = if b_start.distance_to(&pa) <= b_end.distance_to(&pa) {
            (v_a, v_b)
        } else {
            (v_b, v_a)
        };
        let d_from = brep.vertices[v_from].position.distance_to(&b_start);
        let d_to = brep.vertices[v_to].position.distance_to(&b_end);
        if d_from > JUNCTION_ON_EDGE_BAND || d_to > JUNCTION_ON_EDGE_BAND {
            return Err(format!(
                "bridge {k} endpoints miss the junction vertices \
                 (start {d_from:.3e}, end {d_to:.3e})"
            ));
        }
        let he_from = half_edge_at_vertex_on_face(brep, support, v_from)
            .map_err(|e| format!("bridge {k} start not on support {support:?}: {e}"))?;
        let he_to = half_edge_at_vertex_on_face(brep, support, v_to)
            .map_err(|e| format!("bridge {k} end not on support {support:?}: {e}"))?;
        let mef = brep
            .mef(he_from, he_to)
            .map_err(|e| format!("mef bridge {k} across support: {e}"))?;
        write_ring_curve(brep, mef.edge, curve, v_from, v_to);
        // A support face with holes: each hole belongs to the side
        // whose boundary contains it, never by default to the split
        // face (PK-389: the boss hole rode the corner sliver into the
        // cap).
        let rehomed = brep.rehome_inner_loops_after_split(support, mef.face);
        if trace && !rehomed.is_empty() {
            eprintln!(
                "  [boundary-ring] bridge {k}: re-homed {} inner loop(s) across the split",
                rehomed.len()
            );
        }
        let d_old = face_centroid(brep, support)?.distance_to(&vertex_position);
        let d_new = face_centroid(brep, mef.face)?.distance_to(&vertex_position);
        slivers.push(if d_new < d_old { mef.face } else { support });
        ring_edges.insert(mef.edge);
        side_edges[k] = Some(mef.edge);
        if trace {
            eprintln!(
                "  [boundary-ring] bridge {k} support {support:?} \
                 ({:.4},{:.4},{:.4})->({:.4},{:.4},{:.4})",
                b_start.x, b_start.y, b_start.z, b_end.x, b_end.y, b_end.z
            );
        }
    }

    // ---- merge the surrendered interior ------------------------------
    let blend_face = merge_faces(brep, vertex, &slivers, &ring_edges)
        .map_err(|e| format!("saddle ring region merge: {e}"))?;

    // ---- postconditions ----------------------------------------------
    let ring_len = face_loop_half_edges(brep, blend_face)?.len();
    if ring_len != side_count {
        return Err(format!(
            "blend face ring has {ring_len} edges, expected {side_count}"
        ));
    }
    let worst = worst_ring_on_surface_deviation(brep, &ring_edges, blend_face)?;
    if worst > RING_ON_SURFACE_BAND {
        return Err(format!(
            "ring edge sits {worst:.3e} off an adjacent surface (band {RING_ON_SURFACE_BAND:.0e})"
        ));
    }
    if trace {
        eprintln!(
            "  [boundary-ring] blend face {blend_face:?} ring {ring_len} edges, \
             worst on-surface deviation {worst:.3e}"
        );
    }
    let _ = tolerance;
    let side_edges: Vec<EdgeHandle> = side_edges
        .into_iter()
        .collect::<Option<_>>()
        .ok_or("a ring side was not carved")?;
    Ok(CarvedRing {
        blend_face,
        side_is_arc: sides
            .iter()
            .map(|side| matches!(side, RingSide::Arc { .. }))
            .collect(),
        side_edges,
    })
}
