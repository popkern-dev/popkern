//! Fillet end runoff onto planar end faces.
//!
//! A G² sweep is built perpendicular to its spine, so at an edge end
//! where the filleted edge meets an end face at a non-right angle the
//! sweep overshoots past that face's plane (measured on the pentagon:
//! the strip's end border dips 0.012–0.028 below the base). The strip
//! face's boundary meanwhile still runs through the sharp end-corner
//! vertex along the end face's edges, so the rendered strip curls
//! under the floor and the floor's corner pokes out past the fillet.
//!
//! The classical treatment is END RUNOFF: trim the fillet at its
//! intersection curve with the end face, and cut the end face's sharp
//! corner back to the same curve. The runoff curve's endpoints are the
//! trim-crossing anchors (where the contact lines meet the end face),
//! so the surgery is local: one mef across the strip along the runoff
//! curve, then the corner sliver (whose surface region lies in the end
//! face's plane) is merged into the end face.
//!
//! Invariants: Euler-operator-only surgery, snapshot rollback per
//! corner, `validate()`-checked by the caller's pipeline gates. The
//! pass declines (with a recorded reason) for non-planar end faces or
//! when no overshoot is detected — perpendicular ends (cube fixtures)
//! are untouched by construction.

use popkern_core::geometry::curve::{CurveGeometry, CurveOnEdge};
use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::geometry::{NurbsCurve, Point3, SurfaceGeometry, Vec3};
use popkern_core::topology::{BRep, FaceHandle, HalfEdgeHandle, VertexHandle};

/// One applied runoff, for diagnostics.
#[derive(Debug)]
pub(super) struct EndRunoffOutcome {
    /// Strip face that was trimmed.
    pub strip_face: FaceHandle,
    /// End face whose corner was absorbed.
    pub end_face: FaceHandle,
    /// Worst overshoot depth that triggered the runoff.
    pub overshoot: f64,
    /// Why the surgery was rolled back, when it failed.
    pub error: Option<String>,
}

/// Apply end runoff to every fillet strip corner that overshoots a
/// planar end face. Returns the applied outcomes; declines are
/// silent per-corner (no overshoot means nothing to do).
pub(super) fn apply_fillet_end_runoff(brep: &mut BRep) -> Vec<EndRunoffOutcome> {
    let mut outcomes = Vec::new();
    // Corners whose surgery failed: skipped on rescans so one bad
    // corner can't wedge the pass or hide the others. Each failure
    // rolls back to its own snapshot and is reported as an outcome
    // with the error recorded.
    let mut failed: Vec<VertexHandle> = Vec::new();
    loop {
        let Some(candidate) = find_runoff_candidate(brep, &failed) else {
            break;
        };
        let snapshot = brep.clone();
        match runoff_one_corner(brep, &candidate) {
            Ok(()) => outcomes.push(EndRunoffOutcome {
                strip_face: candidate.strip_face,
                end_face: candidate.end_face,
                overshoot: candidate.overshoot,
                error: None,
            }),
            Err(e) => {
                *brep = snapshot;
                failed.push(candidate.corner);
                outcomes.push(EndRunoffOutcome {
                    strip_face: candidate.strip_face,
                    end_face: candidate.end_face,
                    overshoot: candidate.overshoot,
                    error: Some(e),
                });
            }
        }
        if outcomes.len() > 64 {
            break;
        }
    }
    outcomes
}

struct RunoffCandidate {
    strip_face: FaceHandle,
    end_face: FaceHandle,
    corner: VertexHandle,
    /// Loop neighbors of `corner` within the strip loop (anchor
    /// vertices on the contact lines).
    anchor_prev: VertexHandle,
    anchor_next: VertexHandle,
    plane_origin: Point3,
    plane_normal: Vec3,
    overshoot: f64,
}

/// Find one strip corner whose sweep overshoots a planar end face.
fn find_runoff_candidate(brep: &BRep, skip: &[VertexHandle]) -> Option<RunoffCandidate> {
    for (fh, face) in brep.faces.iter() {
        let Some(SurfaceGeometry::Nurbs(_)) = face.surface.as_ref().map(|s| &s.geometry) else {
            continue;
        };
        let surface = &face.surface.as_ref().unwrap().geometry;
        let loop_ = face.outer_loop;
        let start = brep.loops.get(loop_)?.half_edge;
        let mut he = start;
        for _ in 0..10_000 {
            let data = &brep.half_edges[he];
            let next = data.next;
            let he_out = next; // corner → anchor_next
            let corner = brep.half_edges[next].origin;
            let anchor_prev = data.origin;
            let anchor_next = brep.half_edges[brep.half_edges[next].next].origin;
            // Both edges around `corner` must be shared with the SAME
            // plane face (the end face).
            let f_in = face_across(brep, he);
            let f_out = face_across(brep, he_out);
            he = next;
            let (Some(f_end), Some(f_end2)) = (f_in, f_out) else {
                continue;
            };
            if f_end != f_end2 || f_end == fh {
                continue;
            }
            let Some(SurfaceGeometry::Plane { origin, normal, .. }) = brep
                .faces
                .get(f_end)
                .and_then(|f| f.surface.as_ref())
                .map(|s| &s.geometry)
            else {
                continue;
            };
            let (origin, normal) = (*origin, *normal);
            // Overshoot test: the sweep's near-corner border must
            // cross the end plane. Sample the border iso nearest the
            // corner (whichever of u=0 / u=1 is closer).
            let pc = brep.vertices[corner].position;
            let u_iso = if surface.point_at(0.0, 0.5).distance_to(&pc)
                < surface.point_at(1.0, 0.5).distance_to(&pc)
            {
                0.0
            } else {
                1.0
            };
            let signed = |p: Point3| (p - origin).dot(&normal) / normal.length();
            // Sign-convention-free overshoot test: the strip's
            // mid-surface is on the material side of the end plane;
            // any end-border point on the OPPOSITE side is overshoot
            // (the whole border may lie beyond the plane, so a
            // sign-change-on-the-border test misses it).
            let mid_sign = signed(surface.point_at(0.5, 0.5)).signum();
            let overshoot = (0..=16)
                .map(|i| signed(surface.point_at(u_iso, i as f64 / 16.0)))
                .filter(|d| d.signum() != mid_sign)
                .map(f64::abs)
                .fold(0.0_f64, f64::max);
            if overshoot > 1.0e-6 && !skip.contains(&corner) {
                return Some(RunoffCandidate {
                    strip_face: fh,
                    end_face: f_end,
                    corner,
                    anchor_prev,
                    anchor_next,
                    plane_origin: origin,
                    plane_normal: normal,
                    overshoot,
                });
            }
            if he == start {
                break;
            }
        }
    }
    None
}

/// Face on the other side of a half-edge.
fn face_across(brep: &BRep, he: HalfEdgeHandle) -> Option<FaceHandle> {
    let twin = brep.half_edges.get(he)?.twin;
    let loop_ = brep.half_edges.get(twin)?.loop_;
    Some(brep.loops.get(loop_)?.face)
}

/// Trim one strip corner at the end face's plane and absorb the
/// corner sliver into the end face.
fn runoff_one_corner(brep: &mut BRep, cand: &RunoffCandidate) -> Result<(), String> {
    let surface = brep
        .faces
        .get(cand.strip_face)
        .and_then(|f| f.surface.as_ref())
        .map(|s| s.geometry.clone())
        .ok_or("strip lost surface")?;
    let p_prev = brep.vertices[cand.anchor_prev].position;
    let p_next = brep.vertices[cand.anchor_next].position;
    let signed =
        |p: Point3| (p - cand.plane_origin).dot(&cand.plane_normal) / cand.plane_normal.length();

    // Runoff curve: march across the sweep between the two anchors,
    // finding the plane crossing per station by bisection along the
    // sweep's u direction. Anchor UVs come from point inversion.
    use popkern_core::geometry::uv_projection::project_point_to_uv;
    let (ua, va) = project_point_to_uv(&surface, &p_prev).ok_or("anchor_prev inversion")?;
    let (ub, vb) = project_point_to_uv(&surface, &p_next).ok_or("anchor_next inversion")?;
    const N: usize = 33;
    let mut pts: Vec<Point3> = Vec::with_capacity(N);
    for i in 0..N {
        let f = i as f64 / (N - 1) as f64;
        let v = va + (vb - va) * f;
        // Bracket the crossing along u: from the strip's mid-surface
        // (material side) toward the overshooting border.
        let u_border = if (ua + ub) * 0.5 > 0.5 { 1.0 } else { 0.0 };
        let (mut lo, mut hi) = (0.5, u_border);
        let d_lo = signed(surface.point_at(lo, v));
        let d_hi = signed(surface.point_at(hi, v));
        if d_lo.abs() < 1.0e-12 {
            pts.push(surface.point_at(lo, v));
            continue;
        }
        if d_lo.signum() == d_hi.signum() {
            // No crossing at this station (ends exactly on the
            // plane): take the border point projected to the plane.
            let p = surface.point_at(hi, v);
            let d = signed(p);
            let n = cand.plane_normal;
            let n2 = n.dot(&n);
            pts.push(Point3::new(
                p.x - n.x * d / n2.sqrt(),
                p.y - n.y * d / n2.sqrt(),
                p.z - n.z * d / n2.sqrt(),
            ));
            continue;
        }
        for _ in 0..60 {
            let mid = 0.5 * (lo + hi);
            let dm = signed(surface.point_at(mid, v));
            if dm.signum() == d_lo.signum() {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        pts.push(surface.point_at(0.5 * (lo + hi), v));
    }
    // Endpoints snapped to the topology vertices.
    pts[0] = p_prev;
    let last = pts.len() - 1;
    pts[last] = p_next;

    // Degree-1 interpolating NURBS through the runoff points.
    let mut knots = Vec::with_capacity(N + 2);
    knots.push(0.0);
    for i in 0..N {
        knots.push(i as f64 / (N - 1) as f64);
    }
    knots.push(1.0);
    let runoff = NurbsCurve {
        degree: 1,
        knots,
        control_points: pts,
        weights: vec![1.0; N],
    };

    // Surgery: mef anchor_prev → anchor_next across the strip face.
    let (he_a, he_b) =
        half_edge_pair_in_loop(brep, cand.strip_face, cand.anchor_prev, cand.anchor_next)?;
    let mef = brep
        .mef(he_a, he_b)
        .map_err(|e| format!("runoff mef: {e}"))?;
    if let Some(edge) = brep.edges.get_mut(mef.edge) {
        edge.curve = Some(CurveOnEdge {
            geometry: CurveGeometry::Nurbs(runoff),
            t_start: 0.0,
            t_end: 1.0,
        });
    }
    // The sliver = whichever face contains the corner vertex. Its
    // loop must be exactly the triangle [runoff edge, corner edge,
    // corner edge] — anything else means the mef cut somewhere
    // unexpected.
    let sliver = if loop_contains_vertex(brep, mef.face, cand.corner) {
        mef.face
    } else {
        cand.strip_face
    };
    let sliver_len = loop_len(brep, sliver)?;
    if sliver_len != 3 {
        return Err(format!("runoff sliver has {sliver_len} edges, expected 3"));
    }
    // Merge the sliver into the end face. kef keeps the face of the
    // half-edge it is given and kills the twin's — pass the END
    // FACE's half-edge so the sliver (which carries the sweep
    // surface) is the one that dies. (Passing the sliver side gave
    // the merged floor the SWEEP surface — the floor rendered as an
    // extrapolated fillet diving into the solid.)
    let shared = shared_half_edge(brep, sliver, cand.end_face)
        .ok_or("no shared edge between sliver and end face")?;
    let end_side = brep.half_edges[shared].twin;
    brep.kef(end_side).map_err(|e| format!("runoff kef: {e}"))?;
    // Residual bridge edges at the corner: dissolve to fixpoint.
    loop {
        let mut bridge = None;
        let Some(face) = brep.faces.get(cand.end_face) else {
            break;
        };
        let start = brep.loops[face.outer_loop].half_edge;
        let mut he = start;
        for _ in 0..10_000 {
            let data = &brep.half_edges[he];
            let twin = data.twin;
            if brep.half_edges[twin].loop_ == data.loop_ {
                bridge = Some(he);
                break;
            }
            he = data.next;
            if he == start {
                break;
            }
        }
        let Some(he) = bridge else { break };
        brep.dissolve_same_loop_bridge(he)
            .map_err(|e| format!("runoff bridge dissolve: {e}"))?;
    }
    // The corner vertex may be orphaned or 2-valent now; drop it if
    // orphaned (2-valent corners stay — they are legal).
    let orphan = brep
        .half_edges
        .iter()
        .all(|(_, he)| he.origin != cand.corner);
    if orphan {
        brep.vertices.remove(cand.corner);
    }
    Ok(())
}

/// Half-edge pair (at `a`, at `b`) in the same loop of `face`.
fn half_edge_pair_in_loop(
    brep: &BRep,
    face: FaceHandle,
    a: VertexHandle,
    b: VertexHandle,
) -> Result<(HalfEdgeHandle, HalfEdgeHandle), String> {
    let face_data = brep.faces.get(face).ok_or("face missing")?;
    let mut loops = vec![face_data.outer_loop];
    loops.extend(face_data.inner_loops.iter().copied());
    for lh in loops {
        let start = brep.loops.get(lh).ok_or("loop missing")?.half_edge;
        let (mut he_a, mut he_b) = (None, None);
        let mut he = start;
        for _ in 0..10_000 {
            let data = &brep.half_edges[he];
            if data.origin == a {
                he_a = Some(he);
            }
            if data.origin == b {
                he_b = Some(he);
            }
            he = data.next;
            if he == start {
                break;
            }
        }
        if let (Some(x), Some(y)) = (he_a, he_b) {
            return Ok((x, y));
        }
    }
    Err("anchors not in a common loop".into())
}

fn loop_len(brep: &BRep, face: FaceHandle) -> Result<usize, String> {
    let face_data = brep.faces.get(face).ok_or("face missing")?;
    let start = brep
        .loops
        .get(face_data.outer_loop)
        .ok_or("loop missing")?
        .half_edge;
    let mut n = 0;
    let mut he = start;
    for _ in 0..10_000 {
        n += 1;
        he = brep.half_edges[he].next;
        if he == start {
            break;
        }
    }
    Ok(n)
}

fn loop_contains_vertex(brep: &BRep, face: FaceHandle, v: VertexHandle) -> bool {
    let Some(face_data) = brep.faces.get(face) else {
        return false;
    };
    let Some(loop_) = brep.loops.get(face_data.outer_loop) else {
        return false;
    };
    let start = loop_.half_edge;
    let mut he = start;
    for _ in 0..10_000 {
        let data = &brep.half_edges[he];
        if data.origin == v {
            return true;
        }
        he = data.next;
        if he == start {
            break;
        }
    }
    false
}

/// Any half-edge of `face_a` whose twin lies in `face_b`.
fn shared_half_edge(brep: &BRep, face_a: FaceHandle, face_b: FaceHandle) -> Option<HalfEdgeHandle> {
    let face_data = brep.faces.get(face_a)?;
    let start = brep.loops.get(face_data.outer_loop)?.half_edge;
    let mut he = start;
    for _ in 0..10_000 {
        let data = &brep.half_edges[he];
        let twin_loop = brep.half_edges[data.twin].loop_;
        if brep.loops[twin_loop].face == face_b {
            return Some(he);
        }
        he = data.next;
        if he == start {
            break;
        }
    }
    None
}
