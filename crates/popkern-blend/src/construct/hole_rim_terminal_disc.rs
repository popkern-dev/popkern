//! Hole-rim reflex terminal discs (PK-185 / PK-186).
//!
//! **The class.** A fillet on an inner-loop (hole-rim) edge segment
//! whose end vertex is a REFLEX face corner — a rectangular hole
//! corner — with an un-blended transversal hole edge. Every boundary
//! edge incident to the corner on the hole-side support walks along
//! the hole, away from the trim rail, so the anchor-split stage landed the rail end
//! as an interior-rail mev spur ([`super::helpers::ReflexEndAnchor`])
//! and the face-trim stage cut the band between spur tips. After the consumed-merge stage
//! merge, the strip's end boundary runs spur → corner → stub through
//! the ORIGINAL sharp corner — a straight polyline that disagrees
//! with the strip cylinder by exactly `(√2−1)·r` at the corner (the
//! PK-148 stale-corner disease, but here the corner must STAY sharp:
//! it belongs to the un-blended transversal edge).
//!
//! **Repair, by construction.** The strip's true terminal boundary
//! is the exact cross-section arc of its cylinder at the end
//! station, and the region between that arc and the spur/stub
//! polyline is a flat face in the transversal plane (perpendicular
//! to the strip axis through the corner): the rounded-over material
//! ends on a quarter-disc that extends the transversal wall's plane
//! beyond its boundary. One `mef` from the stub anchor to the spur
//! tip across the strip face carves exactly that disc:
//!
//! - the NEW face is the curvilinear triangle (arc, spur, stub) —
//!   it keeps the sharp corner and receives the transversal Plane;
//! - the strip keeps its handle, its cylinder, and gains the arc as
//!   its terminal boundary;
//! - the new edge is authored as the exact minor arc
//!   (`CurveGeometry::Ellipse` with explicit basis, PK-148
//!   convention).
//!
//! **Material side.** A CONVEX rim (pocket opening) removes material
//! within the strip's axial span, so the disc's outward normal
//! points INTO the span; a CONCAVE junction (rib base) adds material
//! within the span, so it points OUT of it. `Mixed` convexity
//! declines.
//!
//! **Invariants:** one mef per end (χ-neutral: E+1, F+1, L+1);
//! semantic ownership — the disc is a new face with the transversal
//! plane, the strip keeps its handle and surface; the corner vertex,
//! spur, and stub edges keep their handles. **Fallback:** any
//! structural or geometric doubt pushes a per-edge failure, which
//! declines the whole operation transactionally (serving the
//! spur/stub polyline as the strip boundary would ship the stale
//! corner). **Tolerance:** the caller's blend tolerance; arc
//! postconditions gate at `linear × 10` (the PK-148 gate).

use popkern_math::PortableMath;
use std::collections::HashMap;

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::{SurfaceGeometry, SurfaceOnFace};
use popkern_core::geometry::{CurveGeometry, CurveOnEdge, Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle, VertexHandle};

use super::EdgeGeometry;
use super::perpendicular_splits::PerpSplitData;
use crate::error::BlendFailure;

/// Author the terminal disc for every reflex hole-corner end
/// recorded in the anchor-split stage. Returns the authored disc faces (protected
/// from sliver absorption downstream). Any decline pushes a per-edge
/// failure — the caller's transactional gate rolls the operation
/// back.
pub(super) fn author_reflex_end_discs(
    result: &mut BRep,
    edge_geometries: &[EdgeGeometry],
    perp_data: &[Option<PerpSplitData>],
    fillet_faces: &HashMap<usize, FaceHandle>,
    tolerance: &Tolerance,
    failures: &mut Vec<(EdgeHandle, BlendFailure)>,
) -> Vec<FaceHandle> {
    let trace = popkern_core::trace::trace_enabled("rim-disc");
    let mut discs = Vec::new();
    for (idx, pd) in perp_data.iter().enumerate() {
        let Some(pd) = pd else { continue };
        if pd.reflex_ends.is_empty() {
            continue;
        }
        let eg = &edge_geometries[idx];
        let fillet_edge = eg.analysis.edge;
        let fail = |reason: String, location: Option<Point3>| {
            (
                fillet_edge,
                BlendFailure::InvalidTopology { reason, location },
            )
        };
        let Some(&strip) = fillet_faces.get(&idx) else {
            failures.push(fail(
                "hole-rim terminal disc: fillet strip face not recorded by the consumed-merge stage".into(),
                None,
            ));
            continue;
        };
        let SurfaceGeometry::Cylinder {
            origin,
            axis,
            radius,
        } = &eg.surface
        else {
            failures.push(fail(
                "hole-rim terminal disc requires a cylindrical strip (non-cylinder \
                 reflex ends are a recorded decline)"
                    .into(),
                None,
            ));
            continue;
        };
        let Some(axis_hat) = axis.normalized() else {
            failures.push(fail("hole-rim terminal disc: degenerate axis".into(), None));
            continue;
        };
        let span_sign = match eg.analysis.convexity {
            crate::analysis::Convexity::Convex => -1.0,
            crate::analysis::Convexity::Concave => 1.0,
            crate::analysis::Convexity::Mixed => {
                failures.push(fail(
                    "hole-rim terminal disc: Mixed convexity has no determined \
                     material side (recorded decline)"
                        .into(),
                    None,
                ));
                continue;
            }
        };
        for re in &pd.reflex_ends {
            let other_anchors = if re.face == eg.analysis.face_1 {
                pd.face_2_anchors
            } else {
                pd.face_1_anchors
            };
            let (a_vertex, same_face_other_end) = if re.at_start {
                let sf = if re.face == eg.analysis.face_1 {
                    pd.face_1_anchors.1
                } else {
                    pd.face_2_anchors.1
                };
                (other_anchors.0, sf)
            } else {
                let sf = if re.face == eg.analysis.face_1 {
                    pd.face_1_anchors.0
                } else {
                    pd.face_2_anchors.0
                };
                (other_anchors.1, sf)
            };
            let pos = |v: VertexHandle| result.vertices.get(v).map(|d| d.position);
            let (Some(p_v), Some(p_b), Some(p_a), Some(p_far)) = (
                pos(re.corner),
                pos(re.spur_tip),
                pos(a_vertex),
                pos(same_face_other_end),
            ) else {
                failures.push(fail(
                    "hole-rim terminal disc: an anchor vertex vanished before the \
                     disc stage"
                        .into(),
                    None,
                ));
                break;
            };
            // Exact station: the arc centre is the axis foot of the
            // rail end.
            let center = *origin + axis_hat * (p_b - *origin).dot(&axis_hat);
            let gate = tolerance.linear().max(1e-9) * 10.0;
            let (ra, rb) = ((p_a - center).length(), (p_b - center).length());
            if (ra - radius).abs() > gate || (rb - radius).abs() > gate {
                failures.push(fail(
                    format!(
                        "hole-rim terminal disc: end anchors off the strip cylinder \
                         (|C-A| = {ra:.6}, |C-B| = {rb:.6}, r = {radius:.6})"
                    ),
                    Some(p_v),
                ));
                break;
            }
            match carve_disc(
                result, strip, re.corner, p_v, p_b, p_a, center, axis_hat, span_sign, p_far,
                tolerance, trace,
            ) {
                Ok(disc) => discs.push(disc),
                Err(reason) => {
                    failures.push(fail(reason, Some(p_v)));
                    break;
                }
            }
        }
    }
    discs
}

/// One end's mef + arc + plane authoring. `p_b` is the spur tip
/// (rail end on the reflex face), `p_a` the stub anchor on the other
/// support, `p_far` the same-reflex-face anchor at the OPPOSITE end
/// (fixes the axial direction of the strip span).
#[allow(clippy::too_many_arguments)]
fn carve_disc(
    result: &mut BRep,
    strip: FaceHandle,
    corner: VertexHandle,
    p_v: Point3,
    p_b: Point3,
    p_a: Point3,
    center: Point3,
    axis_hat: Vec3,
    span_sign: f64,
    p_far: Point3,
    tolerance: &Tolerance,
    trace: bool,
) -> Result<FaceHandle, String> {
    // Locate the spur → corner → stub path in the strip's outer
    // loop: the half-edge whose origin is the corner, flanked by the
    // end anchors.
    let loop_h = result
        .faces
        .get(strip)
        .map(|f| f.outer_loop)
        .ok_or_else(|| "hole-rim terminal disc: strip face vanished".to_string())?;
    let start = result
        .loops
        .get(loop_h)
        .map(|l| l.half_edge)
        .ok_or_else(|| "hole-rim terminal disc: strip loop vanished".to_string())?;
    let mut he_v: Option<HalfEdgeHandle> = None;
    let mut cur = start;
    for _ in 0..result.half_edges.len() + 1 {
        let Some(he) = result.half_edges.get(cur) else {
            break;
        };
        if he.origin == corner {
            he_v = Some(cur);
            break;
        }
        cur = he.next;
        if cur == start {
            break;
        }
    }
    let he_v =
        he_v.ok_or_else(|| "hole-rim terminal disc: corner not on the strip loop".to_string())?;
    let he1 = result.half_edges[he_v].prev;
    let he2 = result.half_edges[he_v].next;
    let o1 = result.half_edges[he1].origin;
    let o2 = result.half_edges[he2].origin;
    let (q1, q2) = (result.vertices[o1].position, result.vertices[o2].position);
    let near = |p: Point3, q: Point3| p.distance_to(&q) < tolerance.linear().max(1e-9) * 100.0;
    if !((near(q1, p_a) && near(q2, p_b)) || (near(q1, p_b) && near(q2, p_a))) {
        return Err(format!(
            "hole-rim terminal disc: strip loop neighbours of the corner are not the \
             end anchors (found ({:.4},{:.4},{:.4}) / ({:.4},{:.4},{:.4}))",
            q1.x, q1.y, q1.z, q2.x, q2.y, q2.z
        ));
    }
    // mef(he1, he2): new edge o1 → o2; the NEW face keeps the
    // o1 → corner → o2 path — the disc triangle with the sharp
    // corner. The strip keeps its handle and its cylinder.
    let mef = result
        .mef(he1, he2)
        .map_err(|e| format!("hole-rim terminal disc mef failed: {e}"))?;
    // Exact minor arc o1 → o2 (the mef edge's primary direction).
    let u_basis = q1 - center;
    let mut v_basis = axis_hat.cross(&u_basis);
    let d = q2 - center;
    let mut sweep = d.dot(&v_basis).patan2(d.dot(&u_basis));
    if sweep < 0.0 {
        v_basis = -v_basis;
        sweep = -sweep;
    }
    if !(1e-9..=std::f64::consts::PI - 1e-9).contains(&sweep) {
        return Err(format!(
            "hole-rim terminal disc: arc sweep {sweep:.4} outside the minor-arc range"
        ));
    }
    let arc = CurveOnEdge {
        geometry: CurveGeometry::Ellipse {
            center,
            u_basis,
            v_basis,
            start_angle: 0.0,
            end_angle: sweep,
        },
        t_start: 0.0,
        t_end: sweep,
    };
    let pos_gate = (tolerance.linear() * 10.0).max(1e-9);
    let hit_start = arc.geometry.point_at(arc.t_start).distance_to(&q1);
    let hit_end = arc.geometry.point_at(arc.t_end).distance_to(&q2);
    if hit_start > pos_gate || hit_end > pos_gate {
        return Err(format!(
            "hole-rim terminal disc: authored arc misses its endpoints \
             (start {hit_start:.3e}, end {hit_end:.3e})"
        ));
    }
    if let Some(e) = result.edges.get_mut(mef.edge) {
        e.curve = Some(arc);
    }
    // Disc plane: through the corner, normal along the axis with the
    // material side fixed by convexity (convex removes material
    // within the span → normal INTO the span; concave adds it →
    // normal OUT of the span).
    let axial = |p: Point3| (p - center).dot(&axis_hat);
    let toward_end = if axial(p_far) <= 0.0 { 1.0 } else { -1.0 };
    let normal = axis_hat * (toward_end * span_sign);
    let u_axis = (p_b - center)
        .normalized()
        .ok_or_else(|| "hole-rim terminal disc: degenerate radial direction".to_string())?;
    let v_axis = normal.cross(&u_axis);
    if let Some(f) = result.faces.get_mut(mef.face) {
        f.surface = Some(SurfaceOnFace {
            geometry: SurfaceGeometry::Plane {
                origin: p_v,
                normal,
                u_axis,
                v_axis,
            },
            orientation: true,
        });
    }
    let _ = result.recompute_face_pcurves(mef.face);
    let _ = result.recompute_face_pcurves(strip);
    if trace {
        eprintln!(
            "[rim-disc] strip {strip:?}: disc {:?} at corner ({:.3},{:.3},{:.3}) \
             sweep {sweep:.4} normal ({:.2},{:.2},{:.2})",
            mef.face, p_v.x, p_v.y, p_v.z, normal.x, normal.y, normal.z
        );
    }
    Ok(mef.face)
}
