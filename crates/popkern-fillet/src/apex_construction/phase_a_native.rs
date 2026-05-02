//! Phase 6c-A2 §28.H Session 3 — interior-anchor topology surgery.
//!
//! Phase A's perp-walk-math places anchors on cube edges at distance
//! `r_k` from V along each cube edge.  For uniform-radius corners,
//! these anchor positions COINCIDE with the planner's apex_vertex
//! targets — Phase A's existing path produces the right topology.
//!
//! For variable-radius corners, the planner's targets sit at
//! sphere-fit positions INTERIOR to host planes, e.g.,
//! `(r_min, r_min, 0)` on the V-side of a host plane sub-face,
//! NOT on any cube edge.  Phase A's perp-walks can't place
//! interior vertices.
//!
//! This module supplies the missing capability via the `mev` +
//! `mev` + `mef` sequence per V-incident plane sub-face.
//!
//! ## Session 3 status (2026-05-02): R1 risk realized
//!
//! The mev+mev+mef sequence creates the planner-position vertices
//! correctly (3 chords on 3 V-incident planes for cube_mixed),
//! but the resulting topology leaves V with **4 same-loop bridges**
//! after stage 3's `kef_cross_face_v_edges` runs.  Each per-plane
//! mef chord splits the V-side sub-face into a V-corner triangle
//! and a rest piece, with 2 cross-face spurs between them.  When
//! stage 3 kef's the first spur on a plane, the V-corner triangle
//! merges into the rest — turning the OTHER spur on that plane
//! into a same-loop bridge.  Same effect cascades for the cube-
//! edge stubs.  Net: V has 4 same-loop bridges that
//! `kev_v_bridge` (which expects valence ≤ 1) cannot remove.
//!
//! `kev` itself can't help: the bridge endpoints (T_a, T_b, or
//! V14/V11/V9 perp-anchors) all have valence ≥ 2 due to the chord
//! or cube-edge connections, so kev's "destination valence 1"
//! precondition fails.
//!
//! See ADR `docs/decisions/2026-04-30-phase-6c-a2-direct-apex-
//! construction.md` §36.5 and §36.9 R1.
//!
//! ## What this lands as-is
//!
//! * Correct topology surgery for the variable-radius case is
//!   committed (the mev+mev+mef sequence runs, creates 6 new
//!   BRep vertices at planner sphere-fit positions, and 3 chord
//!   edges).  When the apex blend loop later runs A2, the
//!   `resolve_apex_vertices` snap finds zero-distance matches at
//!   these new vertices.
//! * Stage 3's `kev_v_bridge` then errors with `V has valence 4`
//!   on cube_mixed.  `try_a2_apex_construction` catches this as
//!   `ApexBuildError::PlanMismatch` and falls back to legacy via
//!   the existing `BlendWarning::A2FellBackToLegacy` mechanism.
//! * Matrix unchanged at fallbacks=4 (no regression).  The
//!   foundation for §28.H proper is in place; Session 3.5 will
//!   resolve the multi-bridge V topology before stage 3 runs.
//!
//! ## Session 3.5 (next) — multi-bridge V dissolve
//!
//! Add a custom Euler operation that, before stage 3, walks V's
//! same-loop bridges and removes them by:
//!
//! 1. For each chord T_a-T_b: find its two endpoints' valences.
//!    If both are valence 2 (spur + chord), kef the chord to
//!    drop one endpoint to valence 1.  No, kef requires cross-
//!    face — chord is same-loop after spur kef.  Need a different
//!    primitive.
//! 2. Custom: "merge same-loop bridge into face's outer loop" —
//!    removes the bridge's 2 half-edges and patches loop's
//!    next/prev pointers to skip them.  Reduces bridge endpoints'
//!    valence by 1.  Repeat until V has valence ≤ 1, then
//!    `kev_v_bridge` runs cleanly.
//!
//! This requires a new `BRep::dissolve_same_loop_bridge` method
//! in `popkern-core/src/euler/`.  Tracked as Session 3.5 in the
//! pickup checklist.

use std::collections::HashMap;

use popkern_core::geometry::Point3;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::topology::{BRep, FaceHandle, HalfEdgeHandle, VertexHandle};

/// Per V-incident plane sub-face, place 2 planner targets via
/// `mev` and connect them with a `mef` chord.  Skips uniform-
/// radius corners (where Phase A already produces the correct
/// topology) and any V where the planner emitted no targets.
///
/// Returns the count of chord edges created.  Caller can use
/// this for telemetry / regression detection.
///
/// `apex_v` is the apex vertex (cube/pyramid corner).  Must
/// still be a valid BRep vertex at call time (Phase A doesn't
/// remove V — that happens in stage 3).
///
/// `planner_targets` is the `apex_vertex` slice from the planner
/// pre-pass (`plan_apex_targets_pre_phase_a` output).  Each
/// entry is one Point3 on a host plane.  Targets are routed to
/// the correct V-incident plane by perpendicular distance to
/// each plane (those with distance < `1e-6` are "on" that plane).
///
/// **Returns `Err` if** any mev/mef step fails (e.g., the chord
/// endpoints aren't on the same loop).  Caller rolls back via
/// the existing `apply_blend_geometry` failure path; the corner
/// then uses the legacy fallback at the apex blend loop.
pub(crate) fn carve_apex_region_native(
    brep: &mut BRep,
    apex_v: VertexHandle,
    planner_targets: &[Point3],
    radii: &[f64],
) -> Result<usize, String> {
    let trace = std::env::var_os("POPKERN_A2_NATIVE_TRACE").is_some();

    // Skip uniform-radius corners: Phase A's perp-walk-math
    // already produces planner-compatible anchors there.  The
    // mev+mev+mef sequence would just create duplicate vertices
    // at perp-anchor positions, which validate_all rejects.
    if is_uniform_radius(radii) {
        if trace {
            eprintln!("[a2-native] V={apex_v:?}: uniform radii — skip");
        }
        return Ok(0);
    }

    // Group planner targets by V-incident plane sub-face.  Each
    // V-incident plane has exactly 2 targets (T_a from one fillet,
    // T_b from the next, in BlendBoundary order).
    let v_incident_planes = collect_v_incident_planes(brep, apex_v);
    if v_incident_planes.is_empty() {
        if trace {
            eprintln!("[a2-native] V={apex_v:?}: no V-incident planes found");
        }
        return Ok(0);
    }

    let mut targets_per_plane: HashMap<FaceHandle, Vec<Point3>> = HashMap::new();
    for face in &v_incident_planes {
        let on_plane: Vec<Point3> = planner_targets
            .iter()
            .filter(|t| point_on_plane_face(brep, *face, **t).unwrap_or(false))
            .copied()
            .collect();
        if on_plane.len() == 2 {
            targets_per_plane.insert(*face, on_plane);
        } else if trace {
            eprintln!(
                "[a2-native] V={apex_v:?} face={face:?}: expected 2 targets on plane, found {} — skip face",
                on_plane.len()
            );
        }
    }

    if targets_per_plane.is_empty() {
        if trace {
            eprintln!("[a2-native] V={apex_v:?}: no plane has the expected 2 targets");
        }
        return Ok(0);
    }

    // For each V-incident plane with 2 targets: mev/mev/mef.
    let mut chords_created = 0usize;
    for (face, targets) in &targets_per_plane {
        let t_a = targets[0];
        let t_b = targets[1];
        if trace {
            eprintln!(
                "[a2-native] V={apex_v:?} face={face:?}: mev/mev/mef T_a={t_a:?} T_b={t_b:?}"
            );
        }
        carve_face(brep, apex_v, *face, t_a, t_b)?;
        chords_created += 1;
    }

    if trace {
        eprintln!(
            "[a2-native] V={apex_v:?}: {chords_created} chords created across {} V-incident planes",
            targets_per_plane.len()
        );
    }
    Ok(chords_created)
}

/// Carve one V-incident plane sub-face: mev V→T_a, mev V→T_b,
/// mef T_a→T_b.  Topology becomes a V-corner triangle (V, T_a,
/// T_b connected by spurs + chord) plus the rest of the sub-face.
fn carve_face(
    brep: &mut BRep,
    apex_v: VertexHandle,
    face: FaceHandle,
    t_a: Point3,
    t_b: Point3,
) -> Result<(), String> {
    // Find V's outgoing half-edge on `face`'s outer loop.
    let he_v_initial = find_v_outgoing_in_face(brep, apex_v, face)
        .ok_or_else(|| format!("V={apex_v:?} not on face={face:?}'s outer loop"))?;

    // mev V→T_a, splicing the spur before V's outgoing HE.
    let mev_a = brep
        .mev(he_v_initial, t_a)
        .map_err(|e| format!("mev V→T_a failed: {e:?}"))?;

    // After mev_a, V has a new outgoing HE (mev_a.half_edge =
    // V→T_a).  V's "old" outgoing (he_v_initial) is now spliced
    // after the spur (it's still V's outgoing on the same face,
    // just appears later in the loop).  Use he_v_initial again
    // for the second mev — same anchor point in V's incidence
    // structure.
    let mev_b = brep
        .mev(he_v_initial, t_b)
        .map_err(|e| format!("mev V→T_b failed: {e:?}"))?;

    // mef T_a→T_b uses the half-edges whose origins are T_a and
    // T_b respectively.  mev_a.half_edge_twin is T_a→V (origin
    // T_a); same for mev_b.
    let mef_chord = brep
        .mef(mev_a.half_edge_twin, mev_b.half_edge_twin)
        .map_err(|e| format!("mef T_a→T_b chord failed: {e:?}"))?;

    // §28.H Session 3.6 — stitch T_a and T_b to the perp-anchor
    // mesh.  Without these chords, T_a/T_b are dangling interior
    // vertices that Session 3.5's dissolve cannot leave in a valid
    // cap face — the loop traversal hits invalid handles.  Each
    // chord connects one planner-position vertex to its nearest
    // perp-anchor on the same V-incident plane (via the rest-of-F
    // sub-face on the non-V-corner side of the just-mef'd chord).
    //
    // mef_chord identifies the new face (V-corner triangle).  The
    // OTHER side of the chord is the "rest of F" — its outer loop
    // contains both T_a, T_b, and the V-incident plane's perp-
    // anchors.  We mef from T to the nearest perp-anchor across
    // rest_of_F.
    // §28.H Session 3.7 — pass mef_chord's half-edges directly
    // to stitch.  This is robust across cube/pyramid topologies:
    // mef_chord.half_edge (T_a→T_b) is on rest-of-F's loop by
    // mef's contract; mev_b.half_edge_twin (T_b→V) is on rest-of-
    // F's loop too (it stays with old_loop after mef).  No
    // heuristic needed — we know exactly which half-edges go with
    // rest-of-F.
    let _ = face; // unused under Session 3.7 — kept for future API
    let _ = apex_v;
    let stitch_a = stitch_t_via_he(brep, mev_a.vertex, mef_chord.half_edge);
    let stitch_b = stitch_t_via_he(brep, mev_b.vertex, mev_b.half_edge_twin);
    if std::env::var_os("POPKERN_A2_NATIVE_TRACE").is_some() {
        eprintln!(
            "[a2-native]   stitch T_a → perp-anchor: {:?}; T_b → perp-anchor: {:?}",
            stitch_a, stitch_b
        );
    }

    Ok(())
}

/// Find an outgoing half-edge of `apex_v` on `face`'s outer loop.
fn find_v_outgoing_in_face(
    brep: &BRep,
    apex_v: VertexHandle,
    face: FaceHandle,
) -> Option<HalfEdgeHandle> {
    let face_data = brep.faces.get(face)?;
    let outer_loop = face_data.outer_loop;
    let loop_data = brep.loops.get(outer_loop)?;
    let start = loop_data.half_edge;
    let mut cursor = start;
    let max_iters = brep.half_edges.len() + 1;
    for _ in 0..max_iters {
        let he = brep.half_edges.get(cursor)?;
        if he.origin == apex_v {
            return Some(cursor);
        }
        cursor = he.next;
        if cursor == start {
            break;
        }
    }
    None
}

/// Collect V-incident planar sub-faces (sub-faces whose outer
/// loop contains V AND whose surface is a Plane).
fn collect_v_incident_planes(brep: &BRep, apex_v: VertexHandle) -> Vec<FaceHandle> {
    let mut out: Vec<FaceHandle> = Vec::new();
    for (fh, face) in &brep.faces {
        // Only planar host faces — cyl/cone/etc. have different
        // boundary semantics and shouldn't get the mev/mef
        // treatment.
        let is_plane = matches!(
            face.surface.as_ref().map(|s| &s.geometry),
            Some(SurfaceGeometry::Plane { .. })
        );
        if !is_plane {
            continue;
        }
        if find_v_outgoing_in_face(brep, apex_v, fh).is_some() {
            out.push(fh);
        }
    }
    out
}

/// Test whether `pt` lies on `face`'s plane (perpendicular
/// distance < `1e-6` m).  Returns `None` if `face` isn't a plane.
fn point_on_plane_face(brep: &BRep, face: FaceHandle, pt: Point3) -> Option<bool> {
    let face_data = brep.faces.get(face)?;
    let SurfaceGeometry::Plane { origin, normal, .. } = face_data.surface.as_ref()?.geometry else {
        return None;
    };
    let d = popkern_core::geometry::Vec3::new(pt.x - origin.x, pt.y - origin.y, pt.z - origin.z);
    Some(d.dot(&normal).abs() < 1e-6)
}

/// §28.H Session 3.7 — stitch a planner-position vertex to the
/// nearest perp-anchor via `mef`, using a known-correct outgoing
/// half-edge on the rest-of-F sub-face.
///
/// `t_v` is one of T_a/T_b created by `carve_face`'s mev (origin
/// of `t_he_on_rest`).  The caller passes the exact half-edge on
/// rest-of-F (the chord's half-edge for T_a, the spur's twin for
/// T_b — both definitionally on rest-of-F per `mef`'s splice
/// contract).  This avoids the heuristic "max loop size" that
/// failed on pyramid topology in Session 3.6.
///
/// Walks the rest-of-F loop, picks the nearest non-T, non-V
/// vertex to `t_v` in 3D, and `mef`s a chord from `t_v` to that
/// vertex.  Returns the perp-anchor handle on success, or a
/// reason string on failure (failures are non-fatal — caller
/// continues).
fn stitch_t_via_he(
    brep: &mut BRep,
    t_v: VertexHandle,
    t_he_on_rest: HalfEdgeHandle,
) -> Result<VertexHandle, String> {
    if !brep.half_edges.contains_key(t_he_on_rest) {
        return Err(format!("t_he_on_rest={t_he_on_rest:?} no longer valid"));
    }
    if brep.half_edges[t_he_on_rest].origin != t_v {
        return Err(format!(
            "t_he_on_rest={t_he_on_rest:?} doesn't originate at t_v={t_v:?}"
        ));
    }

    // Walk the loop containing t_he_on_rest, collect ALL non-T_v
    // candidates with their distance from T_v.  Sort by distance,
    // then pick the closest non-adjacent one.  Adjacency is via
    // next/prev pointers — mef requires he1 and he2 to be non-
    // adjacent in the loop (otherwise it creates a degenerate
    // duplicate-edge or zero-length triangle).
    let loop_h = brep.half_edges[t_he_on_rest].loop_;
    let start = brep.loops[loop_h].half_edge;
    let t_pos = brep.vertices[t_v].position;
    let next_of_t = brep.half_edges[t_he_on_rest].next;
    let prev_of_t = brep.half_edges[t_he_on_rest].prev;

    let mut candidates: Vec<(VertexHandle, HalfEdgeHandle, f64)> = Vec::new();
    let mut cursor = start;
    let mut iters = 0;
    let max_iters = 64;
    while let Some(h) = brep.half_edges.get(cursor) {
        let v = h.origin;
        if v != t_v && cursor != next_of_t && cursor != prev_of_t {
            let p = brep.vertices[v].position;
            let d = t_pos.distance_to(&p);
            candidates.push((v, cursor, d));
        }
        cursor = h.next;
        iters += 1;
        if cursor == start || iters >= max_iters {
            break;
        }
    }
    candidates.sort_by(|a, b| a.2.partial_cmp(&b.2).unwrap_or(std::cmp::Ordering::Equal));
    let (perp_vh, perp_he, _dist) = candidates
        .into_iter()
        .next()
        .ok_or_else(|| format!("T={t_v:?}: no non-adjacent candidates on rest-of-F loop"))?;

    brep.mef(t_he_on_rest, perp_he)
        .map_err(|e| format!("mef T→perp-anchor failed: {e:?}"))?;
    Ok(perp_vh)
}

/// True when all radii are within 1e-6 of each other.
fn is_uniform_radius(radii: &[f64]) -> bool {
    if radii.is_empty() {
        return true;
    }
    let r0 = radii[0];
    radii.iter().all(|&r| (r - r0).abs() < 1e-6)
}

/// §28.H Session 3.6 — connect a planner-position vertex `t_v` to
/// the nearest perp-anchor on the same V-incident plane via a
/// `mef` chord across the "rest of F" sub-face.
///
/// `t_v` is one of T_a/T_b created by `carve_face`'s mev.  After
/// the mef chord T_a-T_b in `carve_face`, `t_v` is on TWO loops:
/// the V-corner triangle and "rest of F" (the other half of the
/// original V-incident plane sub-face).
///
/// We pick a half-edge of `t_v` whose loop is NOT the V-corner
/// triangle (i.e., the "rest of F" loop), then find the nearest
/// perp-anchor on that loop, and `mef` between them.  Splits
/// "rest of F" into 2 sub-faces.
///
/// Returns the perp-anchor that was connected on success, or a
/// reason string on failure.  Failures are non-fatal: the caller
/// continues without this stitch (the topology may still be
/// usable, or the apex blend loop will fall back to legacy).
///
/// Superseded by `stitch_t_via_he` in Session 3.7 (which takes
/// known-correct half-edge handles instead of heuristic-picking).
/// Retained `#[allow(dead_code)]` for diagnostic re-enable.
#[allow(dead_code)]
fn stitch_t_to_perp_anchor(
    brep: &mut BRep,
    t_v: VertexHandle,
    _original_face: FaceHandle,
    apex_v: VertexHandle,
) -> Result<VertexHandle, String> {
    // Find a half-edge originating at t_v.  After carve_face,
    // t_v has multiple outgoings: the spur to V (in V-corner
    // triangle) and another HE on "rest of F" (the loop reached
    // by walking around T's incidences).
    //
    // Identify which loop is "rest of F": the one that does NOT
    // contain V's spur half-edges.  Specifically, V-corner
    // triangle has V's outgoing half-edges (V→T_b for one of
    // them); rest-of-F doesn't.
    let outgoings: Vec<HalfEdgeHandle> = brep
        .half_edges
        .iter()
        .filter(|(_, h)| h.origin == t_v)
        .map(|(h, _)| h)
        .collect();

    // Both loops containing T_v have V on them (V-corner has
    // V via spurs; rest-of-F has V via spurs + cube-edge stubs).
    // Distinguish by loop size: V-corner triangle has 3 half-
    // edges; rest-of-F has many more.  Pick T_v's outgoing whose
    // loop has the most half-edges.
    let _ = apex_v; // unused now but kept for future heuristics
    let mut t_he_on_rest: Option<HalfEdgeHandle> = None;
    let mut max_loop_size = 0usize;
    for &he in &outgoings {
        let Some(he_data) = brep.half_edges.get(he) else {
            continue;
        };
        let loop_h = he_data.loop_;
        let Some(loop_data) = brep.loops.get(loop_h) else {
            continue;
        };
        let start = loop_data.half_edge;
        let mut cursor = start;
        let mut iters = 0;
        let max_iters = 64;
        while brep.half_edges.contains_key(cursor) {
            cursor = brep.half_edges[cursor].next;
            iters += 1;
            if cursor == start || iters >= max_iters {
                break;
            }
        }
        if iters > max_loop_size {
            max_loop_size = iters;
            t_he_on_rest = Some(he);
        }
    }
    let t_he_on_rest = t_he_on_rest.ok_or_else(|| format!("T={t_v:?}: no outgoing half-edges"))?;

    // Walk the loop containing t_he_on_rest, collect candidate
    // perp-anchors (existing vertices that are NOT V, T_a, T_b).
    let loop_h = brep.half_edges[t_he_on_rest].loop_;
    let start = brep.loops[loop_h].half_edge;
    let mut cursor = start;
    let mut iters = 0;
    let max_iters = 64;
    let mut candidates: Vec<(VertexHandle, HalfEdgeHandle)> = Vec::new();
    let t_pos = brep.vertices[t_v].position;
    while let Some(h) = brep.half_edges.get(cursor) {
        let v = h.origin;
        if v != t_v && v != apex_v {
            candidates.push((v, cursor));
        }
        cursor = h.next;
        iters += 1;
        if cursor == start || iters >= max_iters {
            break;
        }
    }

    // Pick the candidate closest to t_v in 3D.
    let mut best: Option<(VertexHandle, HalfEdgeHandle, f64)> = None;
    for (vh, he) in &candidates {
        let p = brep.vertices[*vh].position;
        let d = t_pos.distance_to(&p);
        match best {
            Some((_, _, b_d)) if b_d <= d => {}
            _ => best = Some((*vh, *he, d)),
        }
    }
    let (perp_vh, perp_he, _dist) =
        best.ok_or_else(|| format!("T={t_v:?}: no perp-anchor candidates on rest-of-F loop"))?;

    // mef from T's HE to perp-anchor's HE on the same loop.
    brep.mef(t_he_on_rest, perp_he)
        .map_err(|e| format!("mef T→perp-anchor failed: {e:?}"))?;
    Ok(perp_vh)
}
