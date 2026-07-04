//! Phase A perpendicular edge splitting for fillet topology construction.

use super::*;
use std::collections::HashMap;

pub(super) struct PerpSplitData {
    pub(super) face_1_anchors: (VertexHandle, VertexHandle),
    pub(super) face_2_anchors: (VertexHandle, VertexHandle),
}

pub(super) struct PerpendicularSplitState {
    pub(super) perp_data: Vec<Option<PerpSplitData>>,
    pub(super) descendants: HashMap<EdgeHandle, Vec<EdgeHandle>>,
}

pub(super) fn run_perpendicular_splits(
    brep: &BRep,
    result: &mut BRep,
    edge_geometries: &[EdgeGeometry],
    tolerance: &Tolerance,
    use_direct_apex_construction: bool,
    failures: &mut Vec<(EdgeHandle, BlendFailure)>,
) -> PerpendicularSplitState {
    // ── Phase A: ALL perpendicular edge splits ─ only modifies edges (inserts vertices), not faces.
    let mut perp_data: Vec<Option<PerpSplitData>> = Vec::new();
    let mut closed_handled: std::collections::HashSet<usize> = std::collections::HashSet::new();
    // Round 5r+: snapshot before Phase A mutates result; anchor positions come from here.
    let phase_a_snapshot = result.clone();

    // §28.H Session 1 — planner pre-pass.  Run the DirectApex planner on
    // the pre-Phase-A BRep state to compute apex_vertex targets per
    // V.  Subsequent §28.H sessions (2-3) will thread these targets
    // into Phase A's anchor placement.  Currently the result is
    // computed and discarded — proves the pre-pass runs without
    // regression on every fixture.  Gated by `use_direct_apex_
    // construction` so legacy is unaffected.
    let _planner_targets_per_apex: std::collections::HashMap<
        VertexHandle,
        Vec<popkern_core::geometry::Point3>,
    > = if use_direct_apex_construction {
        crate::apex_construction::planner_prepass::plan_apex_targets_pre_phase_a(
            &phase_a_snapshot,
            edge_geometries,
            tolerance,
        )
    } else {
        std::collections::HashMap::new()
    };
    // Per-original-edge descendants after Phase A, for Phase C consumed-face lookup.
    let mut descendants: HashMap<EdgeHandle, Vec<EdgeHandle>> = HashMap::new();
    for eg in edge_geometries {
        descendants.insert(eg.analysis.edge, vec![eg.analysis.edge]);
    }

    for (idx, eg) in edge_geometries.iter().enumerate() {
        // Closed fillet edges (e.g. cylindrical hole rims) don't
        // have perpendicular boundary edges to split. Dispatch to
        // the closed-fillet surgery path and mark the edge as
        // already handled — skip Phase B/C for it.
        if closed_fillet::is_closed_fillet_edge(result, eg.analysis.edge) {
            match closed_fillet::apply_closed_fillet(result, eg) {
                Ok(()) => {
                    closed_handled.insert(idx);
                    perp_data.push(None);
                    continue;
                }
                Err(failure) => {
                    failures.push((eg.analysis.edge, failure));
                    perp_data.push(None);
                    continue;
                }
            }
        }

        // Only enable the on-cyl anchor placement (#136 commercial-
        // grade fix) when the fillet surface is a Cylinder.  Torus,
        // cone, and other fillet kinds produce different analytic
        // trim curves that would mis-place anchors via the cyl×plane
        // crossing helper.
        // Anchor placement must follow the TRIM CURVES whenever the
        // geometry pipeline produced them — radius-based perp-walk
        // math places anchors where an arc cross-section would land,
        // which disagrees with the G2 sweep's contact lines (PK-095:
        // uniform 0.06 seam between every strip and its faces on the
        // pentagon; the Cylinder-only gate dates from when only the
        // analytic path had trustworthy trims).
        // §28.H direct-apex construction builds its Gregory corner
        // boundary from perp-walk anchor positions; feeding it
        // trim-crossing anchors regresses the pyramid3 native
        // fixtures. Everywhere else the trims are authoritative.
        // Trim curves feed anchor placement for the analytic cylinder
        // path (legacy, diagnostic crossing only) and the G2 NURBS
        // sweeps (PK-095: crossing-derived anchors — perp-walk math
        // disagrees with the sweep's contact lines by ~0.06 on the
        // pentagon). The direct-apex native construction (§28.H)
        // builds its Gregory boundary from perp-walk positions, so
        // crossing anchors stay off there.
        let trim_for_anchors = matches!(
            eg.surface,
            SurfaceGeometry::Cylinder { .. } | SurfaceGeometry::Nurbs(_)
        );
        // Cylinder sweeps join the trim-crossing convention for the
        // setback path (PK-097: the G1 analytic strips carve at
        // perp-walk radius, leaving the setback arcs beyond the main
        // face's boundary). The §28.H direct-apex machinery keeps
        // perp-walk positions on all surface kinds.
        let anchor_on_trim = matches!(
            eg.surface,
            SurfaceGeometry::Nurbs(_) | SurfaceGeometry::Cylinder { .. }
        ) && !use_direct_apex_construction;

        // §28.H Session 2 — build per-anchor planner override.  For
        // each anchor slot (V_start vs V_end on face_1 vs face_2),
        // pick the planner target closest to the perp-walk-math
        // position that would otherwise be used.  When the closest
        // match is far (variable-radius, planner emits sphere-fit
        // positions interior to faces), the override falls back to
        // perp-walk-math at insertion time.
        let override_ = if use_direct_apex_construction && !_planner_targets_per_apex.is_empty() {
            Some(build_perp_anchor_override(
                &phase_a_snapshot,
                eg.analysis.edge,
                eg.analysis.face_1,
                eg.analysis.face_2,
                eg.radius_at_start,
                eg.radius_at_end,
                &_planner_targets_per_apex,
            ))
        } else {
            None
        };

        match helpers::split_perpendicular_edges_by_radius_with_overrides(
            &phase_a_snapshot,
            result,
            eg.analysis.face_1,
            eg.analysis.face_2,
            eg.analysis.edge,
            trim_for_anchors.then_some(&eg.trim_3d_1),
            trim_for_anchors.then_some(&eg.trim_3d_2),
            anchor_on_trim,
            eg.radius_at_start,
            eg.radius_at_end,
            tolerance,
            override_.as_ref(),
        ) {
            Ok(info) => {
                // Propagate every split that happened in this
                // iteration to the descendants map: any fillet
                // edge whose current descendants list contained
                // `old_edge` now tracks `edge_a` and `edge_b`
                // instead.
                for (old_edge, edge_a, edge_b) in &info.edge_splits {
                    for desc_list in descendants.values_mut() {
                        if let Some(pos) = desc_list.iter().position(|e| e == old_edge) {
                            desc_list.splice(pos..=pos, [*edge_a, *edge_b]);
                        }
                    }
                }
                perp_data.push(Some(PerpSplitData {
                    face_1_anchors: info.face_1_anchors,
                    face_2_anchors: info.face_2_anchors,
                }))
            }
            Err(e) => {
                let loc = crate::geometry_compute::edge_midpoint_by_handle(brep, eg.analysis.edge);
                failures.push((
                    eg.analysis.edge,
                    BlendFailure::InvalidTopology {
                        reason: e,
                        location: loc,
                    },
                ));
                perp_data.push(None);
            }
        }
    }

    PerpendicularSplitState {
        perp_data,
        descendants,
    }
}

/// §28.H Session 2 — build per-anchor planner override for one
/// fillet edge.  For each anchor slot (V_start vs V_end on face_1
/// vs face_2), pick the planner target nearest the perp-walk-math
/// reference position by greedy distance.  When a vertex isn't an
/// apex (no planner targets emitted), the corresponding override
/// fields stay `None` and Phase A uses perp-walk math.
///
/// Tolerance for "near": `radius * 1.5` — generous enough that
/// uniform-radius targets (which coincide with perp-walk
/// positions) match cleanly, but tight enough that variable-
/// radius targets (sphere-fit, far from perp-walk) don't get
/// forced onto the wrong slot.  When `None` is returned for a
/// slot due to no near match, Phase A falls back to perp-walk
/// math.
fn build_perp_anchor_override(
    brep: &popkern_core::topology::BRep,
    fillet_edge: EdgeHandle,
    face_1: FaceHandle,
    face_2: FaceHandle,
    radius_at_start: f64,
    radius_at_end: f64,
    planner_targets: &HashMap<VertexHandle, Vec<popkern_core::geometry::Point3>>,
) -> helpers::PerpAnchorOverride {
    use popkern_core::geometry::Point3;
    use popkern_core::geometry::surface::SurfaceGeometry;

    // Identify V_start and V_end of the fillet edge.
    let Some(edge_data) = brep.edges.get(fillet_edge) else {
        return helpers::PerpAnchorOverride::default();
    };
    let primary = edge_data.half_edge;
    let Some(he) = brep.half_edges.get(primary) else {
        return helpers::PerpAnchorOverride::default();
    };
    let Some(twin_he) = brep.half_edges.get(he.twin) else {
        return helpers::PerpAnchorOverride::default();
    };
    let v_start = he.origin;
    let v_end = twin_he.origin;

    // Helper: signed perpendicular distance from `pt` to a Plane
    // face's plane.  For non-plane host faces, returns None — the
    // override only fires for plane faces today (the §28.H scope).
    fn plane_distance(
        brep: &popkern_core::topology::BRep,
        fh: FaceHandle,
        pt: Point3,
    ) -> Option<f64> {
        let face = brep.faces.get(fh)?;
        let SurfaceGeometry::Plane { origin, normal, .. } = face.surface.as_ref()?.geometry else {
            return None;
        };
        let d =
            popkern_core::geometry::Vec3::new(pt.x - origin.x, pt.y - origin.y, pt.z - origin.z);
        Some(d.dot(&normal).abs())
    }

    // For the targets associated with V (start or end of fillet
    // edge), pick the one closest to face_1's plane and the one
    // closest to face_2's plane.  Constrain to targets within
    // `radius * 1.5` of V (so we don't pull from a different
    // apex's targets).
    fn split_targets_by_face(
        brep: &popkern_core::topology::BRep,
        targets: Option<&Vec<Point3>>,
        v_pos: Option<Point3>,
        face_1: FaceHandle,
        face_2: FaceHandle,
        radius_band: f64,
    ) -> (Option<Point3>, Option<Point3>) {
        let Some(targets) = targets else {
            return (None, None);
        };
        let Some(v_pos) = v_pos else {
            return (None, None);
        };
        // Filter to targets within radius_band of V.
        let candidates: Vec<Point3> = targets
            .iter()
            .copied()
            .filter(|t| v_pos.distance_to(t) <= radius_band)
            .collect();

        // For each face, find the candidate with smallest
        // perpendicular distance to that face's plane.
        let on_face_1 = candidates
            .iter()
            .filter_map(|t| plane_distance(brep, face_1, *t).map(|d| (*t, d)))
            .min_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
        let on_face_2 = candidates
            .iter()
            .filter_map(|t| plane_distance(brep, face_2, *t).map(|d| (*t, d)))
            .min_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
        // Only accept when the closest is actually ON the plane
        // (within 1e-6 m).  Targets that don't lie on either face's
        // plane shouldn't be assigned to a face.
        let on_face_1 = on_face_1.filter(|&(_, d)| d < 1e-6).map(|(t, _)| t);
        let on_face_2 = on_face_2.filter(|&(_, d)| d < 1e-6).map(|(t, _)| t);
        // Avoid double-claiming: if the same target lands on
        // both faces (shared-edge case), null face_2.
        let on_face_2 = match (on_face_1, on_face_2) {
            (Some(t1), Some(t2)) if t1.distance_to(&t2) < 1e-9 => None,
            _ => on_face_2,
        };
        (on_face_1, on_face_2)
    }

    let v_start_pos = brep.vertices.get(v_start).map(|v| v.position);
    let v_end_pos = brep.vertices.get(v_end).map(|v| v.position);
    let r_start_band = radius_at_start * 1.5;
    let r_end_band = radius_at_end * 1.5;
    let (v_start_face_1, v_start_face_2) = split_targets_by_face(
        brep,
        planner_targets.get(&v_start),
        v_start_pos,
        face_1,
        face_2,
        r_start_band,
    );
    let (v_end_face_1, v_end_face_2) = split_targets_by_face(
        brep,
        planner_targets.get(&v_end),
        v_end_pos,
        face_1,
        face_2,
        r_end_band,
    );
    helpers::PerpAnchorOverride {
        v_start_face_1,
        v_start_face_2,
        v_end_face_1,
        v_end_face_2,
    }
}
