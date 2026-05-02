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
    brep.mef(mev_a.half_edge_twin, mev_b.half_edge_twin)
        .map_err(|e| format!("mef T_a→T_b chord failed: {e:?}"))?;

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

/// True when all radii are within 1e-6 of each other.
fn is_uniform_radius(radii: &[f64]) -> bool {
    if radii.is_empty() {
        return true;
    }
    let r0 = radii[0];
    radii.iter().all(|&r| (r - r0).abs() < 1e-6)
}
