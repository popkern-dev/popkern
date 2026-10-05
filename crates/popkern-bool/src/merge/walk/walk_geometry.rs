//! Geometry helpers for loop walking.
//!
//! Angular comparison, tangent computation, and surface normal lookup
//! used by the walk module to select the correct next half-edge at
//! intersection nodes.

use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::geometry::{CurveEval, Vec3};
use popkern_core::topology::BRep;
use popkern_math::PortableMath;

use crate::merge::graph::{FaceRef, IntersectionGraph, Operand};

/// Which walk axis the angular pick measures about on CURVED faces.
///
/// Planar faces always use [`planar_walk_normal`] (material-left);
/// this only selects the curved-face axis.
#[derive(Clone, Copy, PartialEq, Eq)]
pub(super) enum WalkAxis {
    /// The inward face normal — material-on-left, matching the planar
    /// convention. Used for junction TURNS, where candidates bound
    /// DIFFERENT regions and the hand decides which region the loop
    /// closes (PK-439: the right-hand walk closed a merged region and
    /// stranded the kept sector's crossings).
    MaterialLeft,
    /// The raw outward surface normal — the historical hand. Used for
    /// self-loop direction CLAIMS, where both candidates close the
    /// SAME region and only cross-walker agreement matters: the chart
    /// walker and the assembly reversal machinery are tuned around
    /// this claim hand.
    OutwardClaim,
}

/// Pick the next edge from candidates using per-face CW angle.
///
/// Projects edge tangent vectors onto the surface tangent plane at the
/// node, then picks the candidate with the smallest clockwise angle from
/// the incoming direction.  Using the surface normal (instead of the
/// cross product of edge directions) is essential for curved surfaces
/// where edges are not coplanar.
pub(super) fn pick_by_angle(
    graph: &IntersectionGraph,
    he_idx: usize,
    candidates: &[usize],
    face: &FaceRef,
    brep_a: &BRep,
    brep_b: &BRep,
) -> Option<usize> {
    pick_by_angle_ext(
        graph,
        he_idx,
        candidates,
        face,
        brep_a,
        brep_b,
        WalkAxis::MaterialLeft,
    )
}

/// [`pick_by_angle`] with an explicit curved-face [`WalkAxis`].
#[allow(clippy::too_many_arguments)]
pub(super) fn pick_by_angle_ext(
    graph: &IntersectionGraph,
    he_idx: usize,
    candidates: &[usize],
    face: &FaceRef,
    brep_a: &BRep,
    brep_b: &BRep,
    axis: WalkAxis,
) -> Option<usize> {
    let twin_idx = graph.half_edges[he_idx].twin;
    let end_node = graph.half_edges[twin_idx].origin;
    let node_pos = graph.nodes[end_node].position;

    let origin_node = graph.half_edges[he_idx].origin;
    // Incoming direction: the TRUE curve tangent at the arrival node
    // when the edge carries a curve (chords of curved edges point
    // across the geometry — two arcs of a split intersection circle
    // share identical chords and become indistinguishable). The curve
    // may live on the TWIN (pairs store it on the forward half only);
    // the nearest-end check below is direction-agnostic, so a borrowed
    // curve needs no extra flip. Self-loop handling falls out of the
    // same rule.
    let incoming_dir = if origin_node == end_node {
        incoming_tangent_for_self_loop(graph, he_idx)
    } else if let Some(curve) = graph.half_edges[he_idx]
        .curve
        .as_ref()
        .or_else(|| graph.half_edges[twin_idx].curve.as_ref())
    {
        let p_start = curve.geometry.point_at(curve.t_start);
        let p_end = curve.geometry.point_at(curve.t_end);
        let at_end = p_end.distance_to(&node_pos) <= p_start.distance_to(&node_pos);
        if at_end {
            curve.geometry.tangent_at(curve.t_end) * -1.0
        } else {
            curve.geometry.tangent_at(curve.t_start)
        }
    } else {
        graph.nodes[origin_node].position - node_pos
    };

    // Try to get the surface normal from the face's geometry.
    // For curved surfaces this is critical for correct angle computation.
    // PLANAR faces use the face's own boundary winding (Newell of the
    // stored outer loop, which points INTO the material in this
    // kernel's convention) instead: the walk consumes boundary
    // half-edges in their stored directions, so fragment loops wind
    // CCW about exactly this vector. The previous candidate-derived
    // guess was self-fulfilling — `incoming × first_candidate` always
    // ranks the first candidate at +90°, so whichever edge happened
    // to be registered first at the node won every straddling pick
    // (measured: a flush coplanar union's kept slope fragment doubled
    // back along a collinear seam record instead of continuing down
    // its own boundary, dead-ending the walk).
    let curved_axis = || match axis {
        WalkAxis::MaterialLeft => curved_walk_normal(graph, end_node, face, brep_a, brep_b),
        WalkAxis::OutwardClaim => surface_normal_at_node(graph, end_node, face, brep_a, brep_b),
    };
    let normal = planar_walk_normal(face, brep_a, brep_b)
        .or_else(curved_axis)
        .unwrap_or_else(|| {
            // Fallback: compute from incoming direction and candidates.
            // For self-loop candidates, use curve tangent instead of position.
            let mut n = Vec3::new(0.0, 0.0, 1.0);
            for &c in candidates {
                let twin_c = graph.half_edges[c].twin;
                let is_self_loop_c = graph.half_edges[c].origin == graph.half_edges[twin_c].origin;
                let out = if is_self_loop_c {
                    outgoing_tangent_for_self_loop(graph, c)
                } else {
                    graph.nodes[graph.half_edges[twin_c].origin].position - node_pos
                };
                if let Some(computed) = incoming_dir.cross(&out).normalized() {
                    n = computed;
                    break;
                }
            }
            n
        });

    let in_proj = incoming_dir - normal * incoming_dir.dot(&normal);
    let ref_dir = match in_proj.normalized() {
        Some(d) => d,
        None => return Some(candidates[0]),
    };

    let mut best = None;
    let mut best_angle = f64::MAX;

    for &c in candidates {
        let twin_c = graph.half_edges[c].twin;
        let end_pos = graph.nodes[graph.half_edges[twin_c].origin].position;
        // For self-loop candidates (origin == end), use the curve tangent
        // at the start parameter as the outgoing direction.
        let out_dir = if graph.half_edges[c].origin == graph.half_edges[twin_c].origin {
            outgoing_tangent_for_self_loop(graph, c)
        } else if let Some(curve) = graph.half_edges[c]
            .curve
            .as_ref()
            .or_else(|| graph.half_edges[twin_c].curve.as_ref())
        {
            // True departure tangent for curved candidates. The curve
            // may live on the TWIN (pairs store it on the forward half
            // only); a chord fallback there mis-headed every borrowed
            // major arc — its chord points the way the SHORT way round
            // goes, tying with the short piece at the same node
            // (PK-439: the flush-rim junction measured both onward
            // arcs at +90° and the wall walk merged two regions).
            // The nearest-end check is direction-agnostic, so a
            // borrowed curve needs no extra flip.
            let p_start = curve.geometry.point_at(curve.t_start);
            let p_end = curve.geometry.point_at(curve.t_end);
            let at_start = p_start.distance_to(&node_pos) <= p_end.distance_to(&node_pos);
            if at_start {
                curve.geometry.tangent_at(curve.t_start)
            } else {
                curve.geometry.tangent_at(curve.t_end) * -1.0
            }
        } else {
            end_pos - node_pos
        };
        let out_proj = out_dir - normal * out_dir.dot(&normal);
        // A degenerate candidate (zero-length tangency self-loop, or a
        // direction parallel to the normal) has no in-plane heading —
        // rank it last, matching the historical atan2(0,0)→wrap
        // behavior so touch-point self-loops still get entered when
        // they are the only remaining candidates.
        let mut angle = match out_proj.normalized() {
            Some(out_unit) => {
                let cos_val = out_unit.dot(&ref_dir);
                let sin_val = ref_dir.cross(&out_unit).dot(&normal);
                let a = sin_val.patan2(cos_val);
                // A candidate nearly COLLINEAR with the arrival edge
                // (angle at the 0/2π boundary, cos ≈ +1 against the
                // backward reference) is a double-back along a
                // geometrically overlapping edge. Construction noise
                // between two representations of the same 3D line (an
                // intersection record vs. the original boundary it
                // overlaps) is far above the exact-zero wrap, so the
                // pick would flip between "wins everything" and
                // "loses everything" on the noise sign. Rank it last
                // deterministically.
                const COLLINEAR_EPS: f64 = 1e-7;
                if cos_val > 0.0 && sin_val.abs() < COLLINEAR_EPS {
                    std::f64::consts::TAU
                } else {
                    a
                }
            }
            None => 0.0,
        };

        if angle <= 1e-10 {
            angle += std::f64::consts::TAU;
        }

        if popkern_core::trace::trace_enabled("fragment-walk") {
            eprintln!("      [angle] cand {c}: angle={angle:.3e}");
        }
        if angle < best_angle {
            best_angle = angle;
            best = Some(c);
        }
    }

    best
}

/// Walk normal for a PLANAR face: the NEGATED Newell of the stored
/// outer loop.
///
/// The stored loop's Newell IS the outward normal (kernel authoring
/// convention — constructors and sweep, confirmed by the PK-104
/// winding census), so this vector points INTO the material. That is
/// not a mistake: `pick_by_angle` measures candidates
/// counter-clockwise from the BACKWARD direction of arrival, and
/// material-on-left traversal needs the smallest CLOCKWISE angle
/// about the OUTWARD normal — equivalently the smallest
/// counter-clockwise angle about this inward vector. The pair
/// (backward reference, inward axis) is exactly the classical DCEL
/// rotational successor for face-on-left walks; see the chart walker
/// (`chart_walk::next_in_chart`), which implements the same
/// comparison in oriented UV. (Verified against four seam junctions
/// of a flush coplanar union plus the convex-corner base case.)
/// Returns `None` for non-planar faces — curved charts keep the
/// tuned UV-normal path.
fn planar_walk_normal(face: &FaceRef, brep_a: &BRep, brep_b: &BRep) -> Option<Vec3> {
    use popkern_core::geometry::surface::SurfaceGeometry;

    let brep = match face.operand {
        Operand::A => brep_a,
        Operand::B => brep_b,
    };
    let face_data = brep.faces.get(face.face)?;
    let surface = face_data.surface.as_ref()?;
    if !matches!(surface.geometry, SurfaceGeometry::Plane { .. }) {
        return None;
    }
    let pts = brep.face_outer_loop_positions(face.face);
    if pts.len() < 3 {
        return None;
    }
    let mut n = Vec3::new(0.0, 0.0, 0.0);
    for i in 0..pts.len() {
        let p = pts[i];
        let q = pts[(i + 1) % pts.len()];
        n.x += (p.y - q.y) * (p.z + q.z);
        n.y += (p.z - q.z) * (p.x + q.x);
        n.z += (p.x - q.x) * (p.y + q.y);
    }
    (n * -1.0).normalized()
}

/// Walk axis for a CURVED face: the INWARD face normal at the node.
///
/// Mirrors [`planar_walk_normal`]'s material-on-left convention:
/// `pick_by_angle` takes the smallest counter-clockwise angle about
/// this axis from the backward arrival direction, and material-on-left
/// traversal needs that angle measured about the INWARD normal. The
/// raw surface normal is orientation-adjusted to the face's outward
/// normal first, then negated.
///
/// The historical axis was the raw outward surface normal —
/// material-on-RIGHT, the mirror of the planar convention. Interior
/// islands never noticed: both directions of an island ring walk as
/// separate loops and winding assigns hole/piece roles, so the turn
/// hand cancels out. A BOUNDARY-TO-BOUNDARY trim on a curved face is
/// where the hand is load-bearing: the face outline's edges carry only
/// ONE direction on this face (their twins live on the adjacent
/// faces), so the kept region can only close material-left, and the
/// right-handed walk spent the trims' crossing directions on the
/// doomed cut-off regions instead — the kept sector stranded (PK-439,
/// flush-cap prism ∪ hub: the wall's rim arc is single-direction on
/// the wall, and the walk merged the back sector with a dropped
/// strip).
fn curved_walk_normal(
    graph: &IntersectionGraph,
    node_idx: usize,
    face: &FaceRef,
    brep_a: &BRep,
    brep_b: &BRep,
) -> Option<Vec3> {
    let n = surface_normal_at_node(graph, node_idx, face, brep_a, brep_b)?;
    let brep = match face.operand {
        Operand::A => brep_a,
        Operand::B => brep_b,
    };
    let outward_aligned = brep
        .faces
        .get(face.face)
        .and_then(|f| f.surface.as_ref())
        .map(|s| s.orientation)
        .unwrap_or(true);
    Some(if outward_aligned { -n } else { n })
}

/// Compute the incoming tangent direction for a self-loop edge.
///
/// For a self-loop (origin == end node), the position-based direction is zero.
/// Instead, evaluate the curve geometry at its end parameter to get the tangent.
/// The incoming direction is negated (pointing back towards where we came from).
pub(super) fn incoming_tangent_for_self_loop(graph: &IntersectionGraph, he_idx: usize) -> Vec3 {
    // Try the forward half-edge's curve first, then the twin's.
    if let Some(ref curve) = graph.half_edges[he_idx].curve {
        let tangent = curve.geometry.tangent_at(curve.t_end);
        // Negate: "incoming" direction at the end node points backwards.
        return tangent * -1.0;
    }
    let twin_idx = graph.half_edges[he_idx].twin;
    if let Some(ref curve) = graph.half_edges[twin_idx].curve {
        // Twin runs in the opposite direction; its tangent at t_start
        // corresponds to the same endpoint. Negate for incoming.
        let tangent = curve.geometry.tangent_at(curve.t_start);
        return tangent;
    }
    // No curve geometry — fall back to zero (pick_by_angle will use candidates[0]).
    Vec3::new(0.0, 0.0, 0.0)
}

/// Compute the outgoing tangent direction for a self-loop candidate edge.
///
/// Similar to [`incoming_tangent_for_self_loop`] but returns the tangent at
/// the start of the curve (the direction the edge departs from the node).
pub(super) fn outgoing_tangent_for_self_loop(graph: &IntersectionGraph, he_idx: usize) -> Vec3 {
    if let Some(ref curve) = graph.half_edges[he_idx].curve {
        return curve.geometry.tangent_at(curve.t_start);
    }
    let twin_idx = graph.half_edges[he_idx].twin;
    if let Some(ref curve) = graph.half_edges[twin_idx].curve {
        // Twin's end tangent (negated) gives the forward direction.
        return curve.geometry.tangent_at(curve.t_end) * -1.0;
    }
    Vec3::new(0.0, 0.0, 0.0)
}

/// Look up the surface normal at a graph node's position on the given face.
///
/// Uses the node's stored UV coordinates (if available) to evaluate the
/// surface normal from the face's geometry.  Returns `None` if the face
/// has no surface geometry or the node has no UV coordinates for this face.
pub(super) fn surface_normal_at_node(
    graph: &IntersectionGraph,
    node_idx: usize,
    face: &FaceRef,
    brep_a: &BRep,
    brep_b: &BRep,
) -> Option<Vec3> {
    use popkern_core::geometry::surface::SurfaceGeometry;

    let brep = match face.operand {
        Operand::A => brep_a,
        Operand::B => brep_b,
    };
    let face_data = brep.faces.get(face.face)?;
    let surface = face_data.surface.as_ref()?;
    let geom = &surface.geometry;
    // Try UV coords first.
    if let Some(&(u, v)) = graph.nodes[node_idx].uv_coords.get(face) {
        let n = geom.normal_at(u, v);
        return n.normalized();
    }

    // Fallback: compute surface normal directly from the 3D position for
    // curved analytic surfaces where the cross-product fallback may fail
    // (e.g. self-loop edges produce zero incoming direction).
    // Planes are excluded: the cross-product fallback works reliably for
    // planar faces and changing the normal source could alter edge selection
    // (measured: both the oriented normal and the raw chart normal broke
    // seven boolean suites — the tuned walk behavior depends on the
    // candidate-derived fallback for planar faces).
    let pos = graph.nodes[node_idx].position;
    match geom {
        SurfaceGeometry::Plane { .. } => None,
        SurfaceGeometry::Cylinder { origin, axis, .. } => {
            let to_pt = pos - *origin;
            let along_axis = *axis * to_pt.dot(axis);
            let radial = to_pt - along_axis;
            radial.normalized()
        }
        SurfaceGeometry::Sphere { center, .. } => {
            let radial = pos - *center;
            radial.normalized()
        }
        _ => None,
    }
}

// Loop measures (Newell, signed area, outer containment) live in
// `walk_loop_measures`; re-exported here so existing
// `walk_geometry::` paths keep working.
pub(crate) use super::walk_loop_measures::loop_newell;
pub(super) use super::walk_loop_measures::{loop_outside_outer, signed_area};
