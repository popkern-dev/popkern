//! In-place routing for parameter-only feature edits (PK-280 stage 3).
//!
//! `modify_feature` normally invalidates the cache from the edited
//! feature onward, so the next evaluation re-executes the feature in
//! full (a fillet re-runs intersection, contact walking, corner
//! construction …). When the edit is a pure PARAMETER tweak of an
//! analytic dress-up feature — fillet radius, chamfer distance, draft
//! angle — and the cache holds the states just before and at the
//! feature, the crate-level in-place modifiers re-solve the geometry
//! on the cached post-feature state directly (measured 35–2693× over
//! remove+recreate for the fillet path). The refreshed snapshot
//! replaces the cache entry, so only DOWNSTREAM features replay.
//!
//! Honesty rules: the route only fires for the default configuration
//! (config overrides re-resolve parameters), only when both cache
//! entries exist, and any decline anywhere falls back to today's
//! full-invalidation replay — the in-place path can refuse, never
//! misbuild (each crate modifier validates before serving).

use popkern_core::topology::{BRep, FaceHandle};

use crate::FeatureTree;
use crate::feature::Feature;
use crate::node::FeatureId;
use crate::parameter::{ParameterId, ParameterTable};

impl FeatureTree {
    /// Attempt the in-place route for an `old → new` feature edit.
    /// Returns the refreshed post-feature snapshot on success.
    pub(crate) fn try_parameter_tweak_in_place(
        &self,
        id: FeatureId,
        old: &Feature,
        new: &Feature,
    ) -> Option<BRep> {
        let trace = popkern_core::trace::trace_enabled("history");
        if self.active_config.is_some() {
            return None;
        }
        let pos = self.timeline.iter().position(|&fid| fid == id)?;
        if pos == 0 {
            return None;
        }
        let state_before = match self.cache.get(None, self.timeline[pos - 1]) {
            Some(s) => s,
            None => {
                if trace {
                    eprintln!(
                        "  [tweak-route] no cache at pos-1 ({:?})",
                        self.timeline[pos - 1]
                    );
                }
                return None;
            }
        };
        let state_at = match self.cache.get(None, id) {
            Some(s) => s,
            None => {
                if trace {
                    eprintln!("  [tweak-route] no cache at feature {id:?}");
                }
                return None;
            }
        };
        let table = &self.parameter_table;

        match (old, new) {
            (
                Feature::Fillet {
                    edges: e_old,
                    radius: r_old,
                    recipe: rc_old,
                },
                Feature::Fillet {
                    edges: e_new,
                    radius: r_new,
                    recipe: rc_new,
                },
            ) if e_old == e_new
                && rc_old == rc_new
                && rc_old.radius_law == crate::feature_blend::FilletRadiusLaw::Constant =>
            {
                let (old_v, new_v) = resolved_pair(table, *r_old, *r_new)?;
                if (old_v - new_v).abs() < 1.0e-15 {
                    return None;
                }
                let mut state = state_at.clone();
                for name in e_old {
                    let edge = match crate::resolve::resolve_edge(name, state_before) {
                        Ok(e) => e,
                        Err(e) => {
                            if trace {
                                eprintln!("  [tweak-route] edge resolve failed: {e}");
                            }
                            return None;
                        }
                    };
                    let supports = translate_supports(state_before, &state, edge)?;
                    let Some(face) = blend_face_between(&state, supports, StripKind::Fillet) else {
                        if trace {
                            eprintln!(
                                "  [tweak-route] no unique fillet strip between {supports:?}"
                            );
                        }
                        return None;
                    };
                    state = match popkern_blend::modify_fillet_radius(
                        &state,
                        face,
                        popkern_blend::FilletRadius::Constant(new_v),
                        state.tolerance,
                    ) {
                        Ok(s) => s,
                        Err(e) => {
                            if trace {
                                eprintln!("  [tweak-route] modify_fillet_radius declined: {e}");
                            }
                            return None;
                        }
                    };
                }
                Some(state)
            }
            (
                Feature::Chamfer {
                    edges: e_old,
                    distance: d_old,
                    kind: k_old,
                },
                Feature::Chamfer {
                    edges: e_new,
                    distance: d_new,
                    kind: k_new,
                },
            ) if e_old == e_new
                && k_old == k_new
                && *k_old == crate::feature_blend::ChamferKindDesc::Equal =>
            {
                let (old_v, new_v) = resolved_pair(table, *d_old, *d_new)?;
                if (old_v - new_v).abs() < 1.0e-15 {
                    return None;
                }
                let mut state = state_at.clone();
                for name in e_old {
                    let edge = crate::resolve::resolve_edge(name, state_before).ok()?;
                    let supports = translate_supports(state_before, &state, edge)?;
                    let face = blend_face_between(&state, supports, StripKind::Chamfer)?;
                    state = popkern_blend::modify_chamfer_distance(
                        &state,
                        face,
                        popkern_blend::ChamferKind::Equal(new_v),
                        state.tolerance,
                    )
                    .ok()?;
                }
                Some(state)
            }
            (
                Feature::Draft {
                    faces: f_old,
                    angle: a_old,
                    neutral_plane: n_old,
                },
                Feature::Draft {
                    faces: f_new,
                    angle: a_new,
                    neutral_plane: n_new,
                },
            ) if f_old == f_new && n_old == n_new => {
                let (old_v, new_v) = resolved_pair(table, *a_old, *a_new)?;
                if (old_v - new_v).abs() < 1.0e-15 {
                    return None;
                }
                let neutral_face = crate::resolve::resolve_face(n_old, state_before).ok()?;
                let (origin, normal) = plane_of(state_before, neutral_face)?;
                let mut state = state_at.clone();
                for name in f_old {
                    // Draft preserves face handles; resolve against the
                    // pre-state and reuse the handle in the post-state.
                    let face = crate::resolve::resolve_face(name, state_before).ok()?;
                    state = popkern_draft::modify_draft_angle(
                        &state,
                        face,
                        old_v,
                        new_v,
                        normal,
                        popkern_draft::NeutralReference::Plane { origin, normal },
                        state.tolerance,
                    )
                    .ok()?
                    .brep;
                }
                Some(state)
            }
            _ => None,
        }
    }
}

/// Resolve both parameter ids; `None` when either is missing.
fn resolved_pair(table: &ParameterTable, old: ParameterId, new: ParameterId) -> Option<(f64, f64)> {
    Some((table.get(old)?.value, table.get(new)?.value))
}

/// The two faces flanking an edge.
fn edge_faces(
    brep: &BRep,
    edge: popkern_core::topology::EdgeHandle,
) -> Option<(FaceHandle, FaceHandle)> {
    let he = brep.edges.get(edge)?.half_edge;
    let hed = brep.half_edges.get(he)?;
    let twin = brep.half_edges.get(hed.twin)?;
    let f1 = brep.loops.get(hed.loop_)?.face;
    let f2 = brep.loops.get(twin.loop_)?.face;
    Some((f1, f2))
}

/// The edge's two support faces from the PRE-state, re-located in the
/// POST-state by surface geometry. Planar supports only (the analytic
/// in-place modifiers require planar supports anyway).
fn translate_supports(
    before: &BRep,
    after: &BRep,
    edge: popkern_core::topology::EdgeHandle,
) -> Option<(FaceHandle, FaceHandle)> {
    let (b1, b2) = edge_faces(before, edge)?;
    let plane =
        |f: FaceHandle| -> Option<(popkern_core::geometry::Point3, popkern_core::geometry::Vec3)> {
            match &before.faces.get(f)?.surface.as_ref()?.geometry {
                popkern_core::geometry::surface::SurfaceGeometry::Plane {
                    origin, normal, ..
                } => Some((*origin, *normal)),
                _ => None,
            }
        };
    let (o1, n1) = plane(b1)?;
    let (o2, n2) = plane(b2)?;
    let tol = after.tolerance.linear();
    Some((
        face_with_plane(after, o1, n1, tol)?,
        face_with_plane(after, o2, n2, tol)?,
    ))
}

/// Re-locate a support face across the feature execution by SURFACE
/// geometry: the operations rebuild face handles (measured: a cube
/// fillet re-homed both supports), but the support PLANES are
/// preserved. `None` when the plane doesn't match a unique face.
fn face_with_plane(
    brep: &BRep,
    origin: popkern_core::geometry::Point3,
    normal: popkern_core::geometry::Vec3,
    tol: f64,
) -> Option<FaceHandle> {
    let mut matches = Vec::new();
    for (f, fd) in &brep.faces {
        let Some(sof) = fd.surface.as_ref() else {
            continue;
        };
        let popkern_core::geometry::surface::SurfaceGeometry::Plane {
            origin: o,
            normal: n,
            ..
        } = &sof.geometry
        else {
            continue;
        };
        let parallel = n.dot(&normal).abs() > 1.0 - 1.0e-9;
        let offset = (origin - *o).dot(n).abs();
        if parallel && offset < tol * 10.0 {
            matches.push(f);
        }
    }
    if matches.len() == 1 {
        Some(matches[0])
    } else {
        None
    }
}

/// What kind of strip the route is looking for between two supports.
#[derive(Clone, Copy)]
enum StripKind {
    Fillet,
    Chamfer,
}

/// The unique blend face in `brep` adjacent to BOTH supports — the
/// strip the feature authored between them. Adjacency alone is
/// ambiguous (a box's END faces neighbour both supports too,
/// measured 3 candidates on the filleted box edge), so candidates
/// filter by kind: a fillet strip identifies as `BlendKind::Fillet`;
/// a chamfer strip is the plane meeting BOTH supports at an oblique
/// dihedral (end faces sit perpendicular). `None` (→ full-replay
/// fallback) when zero or several candidates survive.
fn blend_face_between(
    brep: &BRep,
    supports: (FaceHandle, FaceHandle),
    kind: StripKind,
) -> Option<FaceHandle> {
    let support_normal = |f: FaceHandle| -> Option<popkern_core::geometry::Vec3> {
        match &brep.faces.get(f)?.surface.as_ref()?.geometry {
            popkern_core::geometry::surface::SurfaceGeometry::Plane { normal, .. } => Some(*normal),
            _ => None,
        }
    };
    let mut candidates = Vec::new();
    for f in brep.faces.keys() {
        if f == supports.0 || f == supports.1 {
            continue;
        }
        let neighbours = popkern_blend::identify::blend_neighbours(brep, f);
        if !(neighbours.contains(&supports.0) && neighbours.contains(&supports.1)) {
            continue;
        }
        if popkern_core::trace::trace_enabled("history") {
            eprintln!(
                "  [tweak-route] adj candidate {f:?} identify {:?}",
                popkern_blend::identify::identify_blend(brep, f).map(|b| b.kind)
            );
        }
        let kind_matches = match kind {
            StripKind::Fillet => popkern_blend::identify::identify_blend(brep, f)
                .is_some_and(|b| b.kind == popkern_blend::identify::BlendKind::Fillet),
            StripKind::Chamfer => {
                let oblique = |sf: FaceHandle| -> bool {
                    match (support_normal(sf), support_normal(f)) {
                        (Some(ns), Some(nc)) => {
                            let c = ns.dot(&nc).abs();
                            (0.05..=0.95).contains(&c)
                        }
                        _ => false,
                    }
                };
                oblique(supports.0) && oblique(supports.1)
            }
        };
        if kind_matches {
            candidates.push(f);
        }
    }
    if candidates.len() == 1 {
        Some(candidates[0])
    } else {
        None
    }
}

/// Origin + normal of a planar face.
fn plane_of(
    brep: &BRep,
    face: FaceHandle,
) -> Option<(popkern_core::geometry::Point3, popkern_core::geometry::Vec3)> {
    match &brep.faces.get(face)?.surface.as_ref()?.geometry {
        popkern_core::geometry::surface::SurfaceGeometry::Plane { origin, normal, .. } => {
            Some((*origin, *normal))
        }
        _ => None,
    }
}

#[cfg(test)]
#[path = "modify_inplace_route_tests.rs"]
mod tests;
