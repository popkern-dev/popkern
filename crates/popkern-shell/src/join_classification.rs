//! Gap / smooth / resolvable edge classification for join construction
//! (design §5 Stage 3a, PK-325).
//!
//! Classification is MEASURED from the offset geometry, not derived from
//! orientation-flag algebra (the PK-330 lesson): per station on each
//! interior non-seam edge, the displacement to each adjacent offset
//! surface is an exact projection, and the local offset-plane
//! intersection is tested against the loop-derived into-face directions.
//! An edge is a **gap** exactly when the intersection lies behind both
//! natural footprints — the offset surfaces separate and a join face is
//! required. This one rule covers outward (material-convex), inward
//! (material-concave / erosion) and variable-distance ladders.

use std::collections::{BTreeMap, BTreeSet, HashMap};

use crate::error::{JoinDecline, ShellError, ShellWarning};
use crate::skin_map::SkinMap;
use crate::types::JoinStyle;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, VertexHandle};

/// Fractions of the edge parameter range sampled per classification.
/// Interior stations only: endpoint stations can sit on tangency
/// degeneracies (the wavy plate's top edges flatten toward corners).
const STATION_FRACTIONS: [f64; 5] = [0.1, 0.3, 0.5, 0.7, 0.9];

/// Per-edge join classification.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum EdgeJoinClass {
    /// Both offset surfaces displace the edge identically (G1 pair or
    /// seam) — the shared edge continues, no join.
    Smooth,
    /// The offset surfaces still intersect within their natural
    /// footprints — the Stage-2 recompute ladder serves the edge.
    Resolvable,
    /// The offset surfaces separate over their natural footprints — a
    /// join face is required.
    Gap,
    /// The class varies along the edge or the displacement magnitudes
    /// disagree: no uniform join exists (declines per-edge).
    Irregular,
}

/// The geometric constraint a corner sector's target position solves
/// (PK-326 E1): the closest point to the ORIGINAL vertex on the
/// constraining recomputed entity.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum CornerTarget {
    /// Single-face sector: exact projection onto this OUTER face's
    /// offset surface.
    Face(FaceHandle),
    /// Sector spanning two faces across one resolvable/smooth edge:
    /// nearest point on the two offset surfaces' intersection.
    EdgePair(FaceHandle, FaceHandle),
}

/// Which corner mechanism serves a gap-corner vertex.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum CornerKind {
    /// Every non-seam incident edge is a gap: sphere-window fan peel
    /// (PK-325 Stage 3a).
    AllGap,
    /// Gap and resolvable/smooth edges meet: the polar cone has empty
    /// interior, so the corner resolves by transversal trimming — a
    /// 2-sector peel whose connector edge is the intersection curve of
    /// the two cut-wedge surfaces (PK-326 E1).
    MixedTwoSector,
    /// A free-rim (sheet boundary) vertex with gap edges: the strip
    /// dead-ends at the rim, closed by its canal end cross-section
    /// circle — a 2-sector peel whose connector twin joins the free
    /// rim ring (PK-326, open-sheet thicken).
    BoundaryTwoSector,
}

/// Per-corner plan entry: mechanism plus the retained sector's target
/// constraint (consumed by vertex recompute BEFORE surgery so
/// resolvable-edge recomputes trim against final positions).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct CornerPlanEntry {
    pub kind: CornerKind,
    pub retained: CornerTarget,
}

/// A gap edge scheduled for join construction, in SKIN handles.
#[derive(Debug, Clone)]
pub(crate) struct GapEdge {
    /// The corresponding edge in the input body.
    pub outer_edge: EdgeHandle,
    /// Outer face on the side of `edge.half_edge`.
    pub outer_face1: FaceHandle,
    /// Outer face on the twin side.
    pub outer_face2: FaceHandle,
    /// Measured rolling-ball radius (uniform displacement magnitude).
    pub ball_radius: f64,
}

/// The per-skin join plan: what to construct, what was contained.
#[derive(Debug, Clone)]
pub(crate) struct JoinPlan {
    /// The style in effect after containment downgrades.
    pub effective: JoinStyle,
    /// Gap edges by SKIN edge handle (sorted iteration).
    pub gap_edges: BTreeMap<EdgeHandle, GapEdge>,
    /// Every interior edge's measured class, by SKIN edge handle —
    /// corner construction partitions vertex fans with it.
    pub classes: BTreeMap<EdgeHandle, EdgeJoinClass>,
    /// Gap corners: SKIN vertex handle → corner mechanism + retained
    /// sector target (seams live in the retained sector).
    pub gap_corners: BTreeMap<VertexHandle, CornerPlanEntry>,
    /// Containment / containment-adjacent warnings.
    pub warnings: Vec<ShellWarning>,
}

impl JoinPlan {
    /// Skin edges the Arc path removes from the Stage-2 recompute ladder.
    pub(crate) fn arc_skip_edges(&self) -> BTreeSet<EdgeHandle> {
        if self.effective == JoinStyle::Arc {
            self.gap_edges.keys().copied().collect()
        } else {
            BTreeSet::new()
        }
    }
}

/// Classify every interior edge of the skin and plan gap joins.
///
/// Runs on the freshly cloned skin (original geometry still in place)
/// against the computed offset surfaces. Under `JoinStyle::Arc`, skins
/// whose gap corners reach an unserved vertex class (mixed corners, rim
/// vertices with gap edges — PK-326) downgrade to
/// `JoinStyle::Intersection` with warnings naming each blocking vertex.
///
/// # Errors
///
/// Returns [`ShellError::JoinFailed`] when the effective style is `Arc`
/// and an edge classifies `Irregular` (variable convexity or
/// displacement step along one edge): no uniform join exists, and
/// serving without it would hand back stale geometry.
pub(crate) fn classify_joins(
    skin: &BRep,
    skin_map: &SkinMap,
    offset_surfaces: &HashMap<FaceHandle, SurfaceGeometry>,
    tolerance: &Tolerance,
    requested: JoinStyle,
    serve_rim_corners: bool,
) -> Result<JoinPlan, ShellError> {
    let reverse_face = skin_map.reverse_face_map();
    let reverse_edge = skin_map.reverse_edge_map();
    let reverse_vertex = skin_map.reverse_vertex_map();

    let mut classes: BTreeMap<EdgeHandle, EdgeJoinClass> = BTreeMap::new();
    let mut gap_edges: BTreeMap<EdgeHandle, GapEdge> = BTreeMap::new();
    let mut irregular: Vec<JoinDecline> = Vec::new();

    let edge_handles: Vec<EdgeHandle> = {
        let mut v: Vec<_> = skin.edges.keys().collect();
        v.sort();
        v
    };

    for eh in edge_handles {
        let Some(sides) = interior_edge_sides(skin, eh, &reverse_face) else {
            continue;
        };
        let (outer_f1, outer_f2) = sides;
        if outer_f1 == outer_f2 {
            continue; // seam: single face, never a join
        }
        let (Some(off1), Some(off2)) = (
            offset_surfaces.get(&outer_f1),
            offset_surfaces.get(&outer_f2),
        ) else {
            continue; // rim / unclassified side: rim machinery owns it
        };
        match classify_edge(skin, eh, off1, off2, tolerance) {
            Some((EdgeJoinClass::Gap, radius)) => {
                let Some(&outer_eh) = reverse_edge.get(&eh) else {
                    continue;
                };
                classes.insert(eh, EdgeJoinClass::Gap);
                gap_edges.insert(
                    eh,
                    GapEdge {
                        outer_edge: outer_eh,
                        outer_face1: outer_f1,
                        outer_face2: outer_f2,
                        ball_radius: radius,
                    },
                );
            }
            Some((EdgeJoinClass::Irregular, _)) => {
                classes.insert(eh, EdgeJoinClass::Irregular);
                if let Some(&outer_eh) = reverse_edge.get(&eh) {
                    irregular.push(JoinDecline {
                        edge: outer_eh,
                        mechanism: "edge class varies along the edge (convexity or \
                                    displacement step changes sign); no uniform join \
                                    window exists"
                            .into(),
                        arc_alternative: false,
                    });
                }
            }
            Some((class, _)) => {
                classes.insert(eh, class);
            }
            None => {
                classes.insert(eh, EdgeJoinClass::Resolvable);
            }
        }
    }

    let mut plan = JoinPlan {
        effective: requested,
        gap_edges,
        classes: classes.clone(),
        gap_corners: BTreeMap::new(),
        warnings: Vec::new(),
    };

    if plan.gap_edges.is_empty() && irregular.is_empty() {
        return Ok(plan);
    }

    // Irregular edges (variable convexity or variable per-face
    // distance — no uniform rolling ball) contain the whole skin to
    // the sharp join rather than declining: the extension/intersection
    // path serves exactly what it served before this stage existed,
    // and the warning records the boundary (design §10 E4).
    if requested == JoinStyle::Arc && !irregular.is_empty() {
        plan.effective = JoinStyle::Intersection;
        for decline in irregular {
            plan.warnings.push(ShellWarning::JoinStyleDowngradedAtEdge {
                edge: decline.edge,
                mechanism: decline.mechanism,
            });
        }
        return Ok(plan);
    }

    // Cone-point honesty (Arc only): a dilated cone apex stays SHARP —
    // the Minkowski-faithful result would cap it with a sphere window
    // (PK-326). The apex has no incident gap edges (its edges are
    // seams), so no corner machinery reaches it; the warning records
    // the deviation instead of leaving it silent. Dilation is measured
    // from the exact offset arm: the offset apex displaces AGAINST the
    // axis (which points from the apex into the body).
    if requested == JoinStyle::Arc {
        let mut cone_faces: Vec<(FaceHandle, &SurfaceGeometry)> = offset_surfaces
            .iter()
            .filter(|(_, off)| matches!(off, SurfaceGeometry::Cone { .. }))
            .map(|(fh, off)| (*fh, off))
            .collect();
        cone_faces.sort_by_key(|(fh, _)| *fh);
        for (outer_fh, off) in cone_faces {
            let Some(&skin_fh) = skin_map.faces.get(&outer_fh) else {
                continue;
            };
            let Some(SurfaceGeometry::Cone {
                apex: orig_apex,
                axis,
                ..
            }) = skin
                .faces
                .get(skin_fh)
                .and_then(|f| f.surface.as_ref())
                .map(|s| &s.geometry)
            else {
                continue;
            };
            let SurfaceGeometry::Cone {
                apex: offset_apex, ..
            } = off
            else {
                continue;
            };
            let shift = *offset_apex - *orig_apex;
            if shift.dot(axis) >= 0.0 {
                continue; // erosion or unmoved: the sharp apex is exact
            }
            let mut apex_vertices: Vec<VertexHandle> = skin
                .vertices
                .iter()
                .filter(|(_, v)| v.position.distance_to(orig_apex) <= tolerance.linear() * 10.0)
                .map(|(vh, _)| vh)
                .collect();
            apex_vertices.sort();
            for vh in apex_vertices {
                let outer_vh = reverse_vertex.get(&vh).copied().unwrap_or(vh);
                plan.warnings
                    .push(ShellWarning::ConeApexSharp { vertex: outer_vh });
            }
        }
    }

    // Corner analysis (Arc only): every vertex touching a gap edge must
    // be a corner class the arc machinery serves, else the whole skin
    // downgrades (temporary containment, PK-326 is the replacement).
    if requested == JoinStyle::Arc {
        let mut downgrade: Vec<(VertexHandle, String)> = Vec::new();
        let mut corners: BTreeMap<VertexHandle, CornerPlanEntry> = BTreeMap::new();

        let mut gap_vertices: BTreeSet<VertexHandle> = BTreeSet::new();
        for eh in plan.gap_edges.keys() {
            let edge = &skin.edges[*eh];
            let he = &skin.half_edges[edge.half_edge];
            let twin = &skin.half_edges[he.twin];
            gap_vertices.insert(he.origin);
            gap_vertices.insert(twin.origin);
        }

        let classified: BTreeSet<FaceHandle> = offset_surfaces.keys().copied().collect();
        for vh in gap_vertices {
            match plan_corner(
                skin,
                vh,
                &classes,
                &reverse_face,
                &classified,
                serve_rim_corners,
            ) {
                Ok(entry) => {
                    corners.insert(vh, entry);
                }
                Err(mechanism) => {
                    let outer_vh = reverse_vertex.get(&vh).copied().unwrap_or(vh);
                    downgrade.push((outer_vh, mechanism));
                }
            }
        }

        if downgrade.is_empty() {
            plan.gap_corners = corners;
        } else {
            plan.effective = JoinStyle::Intersection;
            for (vertex, mechanism) in downgrade {
                plan.warnings
                    .push(ShellWarning::JoinStyleDowngraded { vertex, mechanism });
            }
        }
    }

    Ok(plan)
}

/// After the Stage-2 recompute ladder ran under an effective
/// `Intersection` style, convert its silent per-edge declines on GAP
/// edges into the typed transactional decline (design §5: per-edge
/// `JoinFailed` naming the arc alternative).
pub(crate) fn intersection_gap_postcheck(
    plan: &JoinPlan,
    undone_edges: &[EdgeHandle],
) -> Result<(), ShellError> {
    if plan.effective != JoinStyle::Intersection {
        return Ok(());
    }
    let undone: BTreeSet<EdgeHandle> = undone_edges.iter().copied().collect();
    let declines: Vec<JoinDecline> = plan
        .gap_edges
        .iter()
        .filter(|(eh, _)| undone.contains(eh))
        .map(|(_, gap)| JoinDecline {
            edge: gap.outer_edge,
            mechanism: "offset surfaces separate at this gap edge and the \
                        window extensions do not meet (domain-bounded or \
                        non-intersecting offsets)"
                .into(),
            arc_alternative: true,
        })
        .collect();
    if declines.is_empty() {
        Ok(())
    } else {
        Err(ShellError::JoinFailed {
            style: JoinStyle::Intersection,
            declines,
        })
    }
}

/// Resolve the two OUTER faces adjacent to a skin edge, or `None` for
/// boundary/loopless/rim edges.
fn interior_edge_sides(
    skin: &BRep,
    eh: EdgeHandle,
    reverse_face: &HashMap<FaceHandle, FaceHandle>,
) -> Option<(FaceHandle, FaceHandle)> {
    let edge = skin.edges.get(eh)?;
    let he = skin.half_edges.get(edge.half_edge)?;
    let twin = skin.half_edges.get(he.twin)?;
    let f1 = skin.loops.get(he.loop_)?.face;
    let f2 = skin.loops.get(twin.loop_)?.face;
    let of1 = reverse_face.get(&f1).copied()?;
    let of2 = reverse_face.get(&f2).copied()?;
    Some((of1, of2))
}

#[path = "join_classification_corner_plan.rs"]
mod join_classification_corner_plan;
#[path = "join_classification_edges.rs"]
mod join_classification_edges;
use join_classification_corner_plan::plan_corner;
use join_classification_edges::classify_edge;

#[cfg(test)]
#[path = "join_classification_tests.rs"]
mod tests;
