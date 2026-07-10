//! Stage 3: Classify face fragments as inside/outside/on-boundary.
//!
//! Uses multi-ray point-in-solid voting to determine the spatial relationship
//! of each face fragment to the opposite solid. Each ray is cast from a sample
//! point on the fragment, and transversal crossings with the solid's faces are
//! counted. Odd crossings indicate inside, even indicate outside.

use std::collections::{HashMap, HashSet};

use popkern_core::{
    geometry::{
        Point3, Vec3,
        curve::CurveGeometry,
        intersection::{
            SurfaceCurveIntersection, intersect_surface_curve, surface_subdivision::Aabb,
        },
        surface::{SurfaceEval, SurfaceGeometry},
    },
    tolerance::Tolerance,
    topology::{BRep, EdgeHandle, FaceHandle},
};

use crate::error::{BoolError, Classification, PointClassification};

// ---------------------------------------------------------------------------
// 1. Ray-AABB intersection (slab test)
// ---------------------------------------------------------------------------

/// Test whether a ray intersects an axis-aligned bounding box using the slab method.
///
/// Returns `true` if the ray `origin + t * direction` (for t >= 0) passes through
/// or touches `aabb`. Handles axis-aligned rays (direction component = 0) correctly.
fn ray_intersects_aabb(origin: &Point3, direction: &Vec3, aabb: &Aabb) -> bool {
    let mut t_min = f64::NEG_INFINITY;
    let mut t_max = f64::INFINITY;

    let axes = [
        (origin.x, direction.x, aabb.min.x, aabb.max.x),
        (origin.y, direction.y, aabb.min.y, aabb.max.y),
        (origin.z, direction.z, aabb.min.z, aabb.max.z),
    ];

    for (o, d, bmin, bmax) in axes {
        if d.abs() < 1e-15 {
            // Ray is parallel to this slab — miss if origin outside.
            // (1e-15: numerical floor, no tolerance in scope)
            if o < bmin || o > bmax {
                return false;
            }
        } else {
            let inv_d = 1.0 / d;
            let mut t1 = (bmin - o) * inv_d;
            let mut t2 = (bmax - o) * inv_d;
            if t1 > t2 {
                std::mem::swap(&mut t1, &mut t2);
            }
            t_min = t_min.max(t1);
            t_max = t_max.min(t2);
            if t_min > t_max {
                return false;
            }
        }
    }

    // Ray hits if the exit is at t >= 0.
    t_max >= 0.0_f64.max(t_min)
}

// ---------------------------------------------------------------------------
// 2. Face vertex bounding range helper
// ---------------------------------------------------------------------------

/// Compute a finite UV bounding range for a face by examining its vertex
/// positions projected onto the surface parameter space.
///
/// For infinite-range surfaces (planes), projects face loop vertices onto
/// the surface's u/v axes to determine a finite bounding rectangle with
/// generous padding. Returns the surface's own parameter range if it is
/// already finite.
fn face_parameter_range(brep: &BRep, face_handle: FaceHandle) -> ((f64, f64), (f64, f64)) {
    let face = &brep.faces[face_handle];
    let surface_on_face = match face.surface.as_ref() {
        Some(s) => s,
        None => return ((-10.0, 10.0), (-10.0, 10.0)),
    };
    let ((u_min, u_max), (v_min, v_max)) = surface_on_face.geometry.parameter_range();

    // If both ranges are finite, use them directly.
    if u_min.is_finite() && u_max.is_finite() && v_min.is_finite() && v_max.is_finite() {
        return ((u_min, u_max), (v_min, v_max));
    }

    // For infinite surfaces (planes), project boundary SAMPLE
    // positions (vertices + edge-curve samples) onto the surface's
    // u/v axes. Vertex-only projection collapses to a point for a
    // face bounded by a closed rim — a cylinder cap has ONE vertex —
    // so ray classification never counted cap crossings and every
    // point below such a body classified Inside via the lone wall
    // crossing (measured: (12, 0, −0.5) under a Ø70 cylinder).
    // Vertex positions bound the box for ordinary faces — for the
    // ray-crossing counter the box acts as the face-trim proxy, so a
    // WIDER box creates phantom crossings (measured: curve-sampling
    // every face broke disjoint-cube CSG). Only faces whose vertex
    // polygon is degenerate (closed rims: a cylinder cap has ONE
    // vertex) need curve samples — without them their box collapses
    // to the seam point and cap crossings never count (a point below
    // such a body classified Inside via the lone wall crossing).
    let mut positions = face_vertex_positions(brep, face_handle);
    if positions.len() < 3 {
        positions = brep.face_outer_loop_sample_positions(face_handle, 8);
    }
    if positions.is_empty() {
        return ((-100.0, 100.0), (-100.0, 100.0));
    }

    // Extract plane origin and axes for projection.
    let (plane_origin, u_axis, v_axis) = match &surface_on_face.geometry {
        SurfaceGeometry::Plane {
            origin,
            u_axis,
            v_axis,
            ..
        } => (*origin, *u_axis, *v_axis),
        _ => {
            // Non-plane infinite surface — use large default range.
            return ((-100.0, 100.0), (-100.0, 100.0));
        }
    };

    // Project each vertex onto the plane's u/v axes.
    let mut u_lo = f64::INFINITY;
    let mut u_hi = f64::NEG_INFINITY;
    let mut v_lo = f64::INFINITY;
    let mut v_hi = f64::NEG_INFINITY;

    for p in &positions {
        let d = *p - plane_origin;
        let u_val = d.dot(&u_axis);
        let v_val = d.dot(&v_axis);
        u_lo = u_lo.min(u_val);
        u_hi = u_hi.max(u_val);
        v_lo = v_lo.min(v_val);
        v_hi = v_hi.max(v_val);
    }

    // Add generous padding.
    let u_extent = (u_hi - u_lo).max(1.0);
    let v_extent = (v_hi - v_lo).max(1.0);
    u_lo -= u_extent;
    u_hi += u_extent;
    v_lo -= v_extent;
    v_hi += v_extent;

    // Override only the infinite components.
    let u_range = (
        if u_min.is_finite() { u_min } else { u_lo },
        if u_max.is_finite() { u_max } else { u_hi },
    );
    let v_range = (
        if v_min.is_finite() { v_min } else { v_lo },
        if v_max.is_finite() { v_max } else { v_hi },
    );

    (u_range, v_range)
}

/// Collect vertex positions from the outer loop of a face by walking half-edges.
///
/// Uses a bounded iteration limit to prevent infinite loops on broken topology.
fn face_vertex_positions(brep: &BRep, face_handle: FaceHandle) -> Vec<Point3> {
    let face = &brep.faces[face_handle];
    let loop_handle = face.outer_loop;
    let loop_data = match brep.loops.get(loop_handle) {
        Some(l) => l,
        None => return Vec::new(),
    };
    let start_he = loop_data.half_edge;
    if !brep.half_edges.contains_key(start_he) {
        return Vec::new();
    }

    let mut positions = Vec::new();
    let mut current = start_he;
    let limit = brep.half_edges.len() + 1;
    for _ in 0..limit {
        let he = &brep.half_edges[current];
        positions.push(brep.vertices[he.origin].position);
        current = he.next;
        if current == start_he {
            break;
        }
    }
    positions
}

// ---------------------------------------------------------------------------
// 3. Single-ray crossing count
// ---------------------------------------------------------------------------

/// Cast a single ray and count transversal crossings with the faces of a solid.
///
/// Returns `Ok(count)` with the number of transversal face crossings, or
/// `Err(())` if any intersection is ambiguous (near-tangent or near an edge).
fn ray_crossing_count(
    origin: &Point3,
    direction: &Vec3,
    brep: &BRep,
    face_aabbs: &[(FaceHandle, Aabb)],
    tolerance: &Tolerance,
    skip_hole_hits: bool,
) -> Result<usize, ()> {
    let ray_curve = CurveGeometry::Line {
        origin: *origin,
        direction: *direction,
    };
    let t_range = (0.0, 1e6);
    let mut count = 0usize;

    for (face_handle, aabb) in face_aabbs {
        if !ray_intersects_aabb(origin, direction, aabb) {
            continue;
        }

        let face = &brep.faces[*face_handle];
        let surface_on_face = match face.surface.as_ref() {
            Some(s) => s,
            None => continue,
        };
        let surface = &surface_on_face.geometry;

        let ((u_min, u_max), (v_min, v_max)) = face_parameter_range(brep, *face_handle);

        let result = intersect_surface_curve(
            surface,
            (u_min, u_max),
            (v_min, v_max),
            &ray_curve,
            t_range,
            tolerance,
        )
        .map_err(|_| ())?;

        for isect in &result.intersections {
            match isect {
                SurfaceCurveIntersection::Point { t, u, v, .. } => {
                    // Only forward intersections (t > 0).
                    if *t < tolerance.linear() {
                        continue;
                    }

                    let normal = surface.normal_at(*u, *v);
                    let n_len = normal.length();
                    if n_len < tolerance.numerical_floor() {
                        // Degenerate normal — skip.
                        continue;
                    }
                    let n_unit = normal * (1.0 / n_len);
                    let d_len = direction.length();
                    if d_len < tolerance.numerical_floor() {
                        return Err(());
                    }
                    let d_unit = *direction * (1.0 / d_len);

                    let dot = n_unit.dot(&d_unit).abs();

                    // If the hit is nearly tangent, the ray is inconclusive.
                    if dot < tolerance.angular() {
                        return Err(());
                    }

                    // A hit inside one of the face's INNER LOOPS is
                    // a hole, not a crossing: the parameter box is the
                    // trim proxy and covers holes (measured: a point
                    // in a counterbore void classified Inside because
                    // the ray "hit" the holed top face through its
                    // own hole; the tool's proud band was kept and the
                    // result went non-manifold). Scoped to PLANE
                    // faces with inner loops — the general trim-aware
                    // gate was measured to break simple cube CSG, and
                    // holes on periodic charts have unresolved seam
                    // branches.
                    let is_plane = matches!(
                        surface,
                        popkern_core::geometry::surface::SurfaceGeometry::Plane { .. }
                    );
                    // CYLINDER hits: accept only within the face's
                    // EXACT angular span (vertex-derived, +tol). The
                    // padded parameter window is solver headroom, and
                    // rays passing BEHIND a partial arc wall "hit"
                    // the phantom half of the cylinder (measured: an
                    // arm-socket tool's discarded fragments sampled
                    // near a 180° arc front classified Inside and the
                    // assembly went non-manifold). Full closed rims
                    // have degenerate vertex spans and keep the
                    // padded behavior; FINAL bodies only, like the
                    // plane gate below.
                    if skip_hole_hits
                        && let popkern_core::geometry::surface::SurfaceGeometry::Cylinder {
                            origin: cyl_origin,
                            axis: cyl_axis,
                            ..
                        } = surface
                    {
                        // Parameter conventions differ between the
                        // solver window and the stored chart, so gate
                        // GEOMETRICALLY: the hit's 3D point (from the
                        // ray, not the u/v) must fall within the
                        // face's boundary-derived angular span in the
                        // surface's own basis. The span comes from
                        // CURVE SAMPLES, not vertices alone: a face
                        // spanning exactly half the cylinder has all
                        // its vertices ON the atan2 ±π cut, where the
                        // covered half is decided by the sign of a
                        // zero (measured: a boolean emitted −0.0 in a
                        // vertex coordinate and the span flipped to
                        // the phantom half, rejecting every interior
                        // hit — a prow-rear block read Outside at an
                        // interior point and a seat-cutter wall was
                        // dropped). The span is the minimal circular
                        // arc covering the samples (largest-gap), so
                        // cut placement cannot bias it; a full rim's
                        // gap is one sample spacing, which yields
                        // span ≈ τ and keeps the padded behavior.
                        let samples = brep.face_boundary_sample_points(*face_handle, 8);
                        if samples.len() >= 3 {
                            let (perp1, perp2) =
                                popkern_core::geometry::analytic_surfaces::orthonormal_basis(
                                    *cyl_axis,
                                );
                            let ang = |p: &popkern_core::geometry::Point3| -> f64 {
                                let d = *p - *cyl_origin;
                                d.dot(&perp2).atan2(d.dot(&perp1))
                            };
                            let mut angs: Vec<f64> = samples.iter().map(&ang).collect();
                            angs.sort_by(|a, b| a.total_cmp(b));
                            let tau = std::f64::consts::TAU;
                            let mut max_gap = 0.0;
                            let mut span_start = angs[0];
                            for i in 0..angs.len() {
                                let a = angs[i];
                                let b = if i + 1 < angs.len() {
                                    angs[i + 1]
                                } else {
                                    angs[0] + tau
                                };
                                if b - a > max_gap {
                                    max_gap = b - a;
                                    span_start = b;
                                }
                            }
                            let span = tau - max_gap;
                            if span > 1.0e-6 && span <= std::f64::consts::PI + 1.0e-6 {
                                let hit3d = *origin + *direction * *t;
                                let rel = (ang(&hit3d) - span_start).rem_euclid(tau);
                                let margin = 1.0e-6;
                                if rel > span + margin && rel < tau - margin {
                                    if std::env::var_os("PK_RAY_TRACE").is_some() {
                                        eprintln!(
                                            "  [ray-rej-cyl] {face_handle:?} rel={rel:.3} span={span:.3} from {span_start:.3}"
                                        );
                                    }
                                    continue;
                                }
                            }
                        }
                    }
                    // AXIAL gate for cylinder and cone hits: the hit's
                    // coordinate along the surface axis must fall
                    // within the face's boundary-derived axial span
                    // (+one relative margin). The padded parameter box
                    // is solver headroom, and near a chamfered rim a
                    // near-horizontal ray legitimately crosses the cap
                    // and cone, then ALSO "hits" the wall's infinite
                    // cylinder BELOW the weld (measured: rays from
                    // under a rim-chamfered column counted cap + cone
                    // + phantom wall at v=0.29 against a wall trimmed
                    // to [1,59] — every below-body tool fragment
                    // classified Inside and the subtract kept a skirt;
                    // the mirror case accepted an extended-cone hit at
                    // z=2.3 on a cone banded [0,1]). Geometric, like
                    // the angular gate: boundary curve samples, not
                    // pcurves, so seam branches cannot bias it. FINAL
                    // bodies only.
                    if skip_hole_hits {
                        let axis_frame = match surface {
                            popkern_core::geometry::surface::SurfaceGeometry::Cylinder {
                                origin: o,
                                axis,
                                ..
                            } => Some((*o, *axis)),
                            popkern_core::geometry::surface::SurfaceGeometry::Cone {
                                apex,
                                axis,
                                ..
                            } => Some((*apex, *axis)),
                            _ => None,
                        };
                        if let Some((axis_origin, axis_dir)) = axis_frame {
                            let samples = brep.face_boundary_sample_points(*face_handle, 8);
                            if samples.len() >= 3 {
                                let axial = |p: &popkern_core::geometry::Point3| -> f64 {
                                    (*p - axis_origin).dot(&axis_dir)
                                };
                                let (mut lo, mut hi) = (f64::MAX, f64::MIN);
                                for sp in &samples {
                                    let a = axial(sp);
                                    lo = lo.min(a);
                                    hi = hi.max(a);
                                }
                                let hit3d = *origin + *direction * *t;
                                let a = axial(&hit3d);
                                // Margin: 5% of the span (+tolerance).
                                // The coincident-overlap machinery
                                // probes at ~1e-3 of flush interfaces
                                // that sit exactly ON a face's axial
                                // end — a tight margin rejected those
                                // knife-edge hits and re-tuned flush
                                // verdicts (measured: the flush gusset
                                // union lost a cooriented OnBoundary
                                // and kept a duplicate skin). The rim
                                // phantoms this gate exists for sit
                                // 20–130% outside their spans.
                                let margin = tolerance.linear().max(0.05 * (hi - lo));
                                if a < lo - margin || a > hi + margin {
                                    if std::env::var_os("PK_RAY_TRACE").is_some() {
                                        eprintln!(
                                            "  [ray-rej-axial] {face_handle:?} a={a:.3} span=[{lo:.3},{hi:.3}]"
                                        );
                                    }
                                    continue;
                                }
                            }
                        }
                    }
                    // A hit OUTSIDE the face's outer pcurve polygon is
                    // not on the face at all — the padded parameter box
                    // is a solver window, not a trim (measured: a ray
                    // passing 6 mm UNDER a pocket wall crossed its
                    // padded box and a buried tool fragment classified
                    // Outside, dropping two pocket walls). Scoped to
                    // PLANE charts (affine pcurves, no seam branches)
                    // and to FINAL bodies (`skip_hole_hits`): mid-split
                    // states carry STALE pcurves — recompute runs at
                    // the end of the operation — and a stale polygon
                    // rejects hits that ARE on the face (measured:
                    // cube ∩ cube misclassified during the stage
                    // classification).
                    if skip_hole_hits
                        && is_plane
                        && !hit_inside_outer_loop(brep, *face_handle, *u, *v)
                    {
                        if std::env::var_os("PK_RAY_TRACE").is_some() {
                            eprintln!("  [ray-rej] face {face_handle:?} uv=({u:.2},{v:.2})");
                        }
                        continue;
                    }
                    if skip_hole_hits
                        && is_plane
                        && hit_inside_inner_loop(brep, *face_handle, *u, *v)
                    {
                        continue;
                    }
                    if std::env::var_os("PK_RAY_ACC").is_some() {
                        let h = *origin + *direction * *t;
                        eprintln!(
                            "  [ray-acc] {face_handle:?} t={t:.2} hit=({:.2},{:.2},{:.2})",
                            h.x, h.y, h.z
                        );
                    }
                    // Transversal crossing.
                    count += 1;
                }
                SurfaceCurveIntersection::Overlap { .. } => {
                    // Ray lies on the surface — inconclusive.
                    return Err(());
                }
            }
        }
    }

    Ok(count)
}

/// Sample a loop's pcurves into a UV polygon and even-odd test
/// `(u, v)` against it. Returns `None` when the loop has no usable
/// pcurve polygon (fewer than 3 samples) — the caller must treat
/// that as "no verdict", not "outside".
fn uv_point_in_loop(
    brep: &BRep,
    lh: popkern_core::topology::LoopHandle,
    u: f64,
    v: f64,
) -> Option<bool> {
    let loop_data = brep.loops.get(lh)?;
    let start = loop_data.half_edge;
    if !brep.half_edges.contains_key(start) {
        return None;
    }
    let mut ring: Vec<(f64, f64)> = Vec::new();
    let mut he = start;
    let cap = brep.half_edges.len() + 1;
    for _ in 0..cap {
        let d = &brep.half_edges[he];
        if let Some(pc) = d.pcurve.as_ref() {
            let (t0, t1) = pc.parameter_range();
            for i in 0..8 {
                let t = t0 + (t1 - t0) * i as f64 / 8.0;
                let p = pc.point_at(t);
                ring.push((p.x, p.y));
            }
        }
        he = d.next;
        if he == start {
            break;
        }
    }
    if ring.len() < 3 {
        return None;
    }
    let mut crossings = 0;
    let n = ring.len();
    for i in 0..n {
        let (x0, y0) = ring[i];
        let (x1, y1) = ring[(i + 1) % n];
        if (y0 > v) != (y1 > v) {
            let xc = x0 + (v - y0) / (y1 - y0) * (x1 - x0);
            if xc > u {
                crossings += 1;
            }
        }
    }
    Some(crossings % 2 == 1)
}

/// True when the surface hit at `(u, v)` lands inside the face's
/// OUTER pcurve polygon. The verdict is trusted when the outer ring
/// is COMPLETE (every half-edge carries a pcurve). Linear pcurves
/// contribute their exact endpoints; CURVED pcurves are sampled
/// densely (64 stations — sagitta ≲ 0.2% of the span), and hits
/// within a dead-band of the polygon boundary abstain (keep the
/// hit) so sampling error can never reject a rim-grazing crossing
/// (measured: an 8-station circle misjudged a 2.7 mm band at r=35
/// and broke star-socket caps). A hit clearly outside — like the
/// 0.67-off phantom on an ellipse-fronted top face that flipped a
/// tool fragment to Inside — is rejected. Partial rings
/// mid-pipeline still abstain entirely.
fn hit_inside_outer_loop(brep: &BRep, face: FaceHandle, u: f64, v: f64) -> bool {
    let Some(face_data) = brep.faces.get(face) else {
        return true;
    };
    let Some(loop_data) = brep.loops.get(face_data.outer_loop) else {
        return true;
    };
    let start = loop_data.half_edge;
    if !brep.half_edges.contains_key(start) {
        return true;
    }
    // Collect each half-edge's pcurve segment, then chain them
    // orientation-aware: a stored pcurve may run against the loop
    // traversal (twin-shared spans), so each segment is flipped to
    // continue from the previous tail.
    let mut has_curved = false;
    let mut segs: Vec<((f64, f64), (f64, f64))> = Vec::new();
    let mut he = start;
    let cap = brep.half_edges.len() + 1;
    for _ in 0..cap {
        let d = &brep.half_edges[he];
        let Some(pc) = d.pcurve.as_ref() else {
            return true; // incomplete ring — no verdict
        };
        let linear_content = match &pc.geometry {
            popkern_core::geometry::pcurve::PcurveGeometry::Linear { .. } => true,
            popkern_core::geometry::pcurve::PcurveGeometry::Nurbs(c) => c.degree == 1,
            _ => false,
        };
        let (t0, t1) = pc.parameter_range();
        if linear_content {
            let p0 = pc.point_at(t0);
            let p1 = pc.point_at(t1);
            segs.push(((p0.x, p0.y), (p1.x, p1.y)));
        } else {
            has_curved = true;
            // Curved pcurve: densify into chained sub-segments.
            let n = 64;
            let mut prev = pc.point_at(t0);
            for i in 1..=n {
                let t = t0 + (t1 - t0) * i as f64 / n as f64;
                let q = pc.point_at(t);
                segs.push(((prev.x, prev.y), (q.x, q.y)));
                prev = q;
            }
        }
        he = d.next;
        if he == start {
            break;
        }
    }
    if segs.len() < 3 {
        return true;
    }
    let mut ring: Vec<(f64, f64)> = Vec::with_capacity(segs.len());
    let d2 = |a: (f64, f64), b: (f64, f64)| (a.0 - b.0).powi(2) + (a.1 - b.1).powi(2);
    let (mut tail, first_tail) = {
        // Orient the first segment by whichever end the SECOND
        // segment continues from.
        let (a0, a1) = segs[0];
        let (b0, b1) = segs[1];
        let fwd = d2(a1, b0).min(d2(a1, b1));
        let rev = d2(a0, b0).min(d2(a0, b1));
        if fwd <= rev {
            ring.push(a0);
            (a1, a0)
        } else {
            ring.push(a1);
            (a0, a1)
        }
    };
    for &(p0, p1) in segs.iter().skip(1) {
        if d2(tail, p0) <= d2(tail, p1) {
            ring.push(p0);
            tail = p1;
        } else {
            ring.push(p1);
            tail = p0;
        }
    }
    let _ = first_tail;
    // Dead-band: within sampling error of the boundary, abstain
    // (keep the hit) — the polygon is approximate for curved rims.
    let n = ring.len();
    let (mut lo_u, mut hi_u, mut lo_v, mut hi_v) = (f64::MAX, f64::MIN, f64::MAX, f64::MIN);
    for &(x, y) in &ring {
        lo_u = lo_u.min(x);
        hi_u = hi_u.max(x);
        lo_v = lo_v.min(y);
        hi_v = hi_v.max(y);
    }
    // The dead-band exists for SAMPLING error only: linear rings are
    // exact and get none (a 1 mm band on an exact ring re-admitted
    // the very phantoms the gate exists to reject).
    let band = if has_curved {
        0.002 * ((hi_u - lo_u) + (hi_v - lo_v)).max(1e-9)
    } else {
        0.0
    };
    if band > 0.0 {
        let mut min_d = f64::MAX;
        for i in 0..n {
            let (x0, y0) = ring[i];
            let (x1, y1) = ring[(i + 1) % n];
            let (dx, dy) = (x1 - x0, y1 - y0);
            let l2 = (dx * dx + dy * dy).max(1e-300);
            let t = (((u - x0) * dx + (v - y0) * dy) / l2).clamp(0.0, 1.0);
            let d = ((u - x0 - dx * t).powi(2) + (v - y0 - dy * t).powi(2)).sqrt();
            min_d = min_d.min(d);
        }
        if min_d < band {
            return true;
        }
    }
    let mut crossings = 0;
    for i in 0..n {
        let (x0, y0) = ring[i];
        let (x1, y1) = ring[(i + 1) % n];
        if (y0 > v) != (y1 > v) {
            let xc = x0 + (v - y0) / (y1 - y0) * (x1 - x0);
            if xc > u {
                crossings += 1;
            }
        }
    }
    crossings % 2 == 1
}

/// True when the surface hit at `(u, v)` lands inside one of the
/// face's INNER-loop polygons (pcurve-sampled, even-odd).
fn hit_inside_inner_loop(brep: &BRep, face: FaceHandle, u: f64, v: f64) -> bool {
    let Some(face_data) = brep.faces.get(face) else {
        return false;
    };
    face_data
        .inner_loops
        .iter()
        .any(|&lh| uv_point_in_loop(brep, lh, u, v) == Some(true))
}

// ---------------------------------------------------------------------------
// 4. Point-in-solid with multi-ray voting
// ---------------------------------------------------------------------------

/// Directions for the multi-ray voting scheme.
///
/// These are intentionally non-axis-aligned to reduce the chance of hitting
/// edges, vertices, or coincident surfaces.
const RAY_DIRECTIONS: [Vec3; 6] = [
    Vec3 {
        x: 1.0,
        y: 0.1,
        z: 0.05,
    },
    Vec3 {
        x: 0.05,
        y: 1.0,
        z: 0.1,
    },
    Vec3 {
        x: 0.1,
        y: 0.05,
        z: 1.0,
    },
    Vec3 {
        x: -0.7,
        y: 0.5,
        z: 0.3,
    },
    Vec3 {
        x: 0.3,
        y: -0.7,
        z: 0.5,
    },
    Vec3 {
        x: 0.5,
        y: 0.3,
        z: -0.7,
    },
];

/// Classify a 3D point as inside, outside, or on the boundary of a solid.
///
/// Uses multi-ray voting: casts several rays from `point` and counts
/// transversal crossings with the solid's faces. An odd count means inside,
/// even means outside. A majority vote from conclusive rays determines the
/// result.
///
/// # Errors
///
/// Returns [`BoolError::ClassificationFailed`] if no conclusive ray can be
/// found after exhausting all candidate directions.
pub fn point_in_solid(
    point: &Point3,
    solid: &BRep,
    face_aabbs: &[(FaceHandle, Aabb)],
    tolerance: &Tolerance,
) -> Result<PointClassification, BoolError> {
    point_in_solid_impl(point, solid, face_aabbs, tolerance, false)
}

/// [`point_in_solid`] for FINAL bodies (boolean inputs/outputs, not
/// mid-split states): a hit landing inside a face's inner loop is a
/// genuine VOID and is skipped. Mid-split "inner loops" are covered
/// by sibling island faces, where the double-hit self-corrects by
/// parity — skipping there flips classifications (measured: a plain
/// wedge subtract classified everything Inside and returned an empty
/// body).
pub fn point_in_solid_final(
    point: &Point3,
    solid: &BRep,
    face_aabbs: &[(FaceHandle, Aabb)],
    tolerance: &Tolerance,
) -> Result<PointClassification, BoolError> {
    point_in_solid_impl(point, solid, face_aabbs, tolerance, true)
}

fn point_in_solid_impl(
    point: &Point3,
    solid: &BRep,
    face_aabbs: &[(FaceHandle, Aabb)],
    tolerance: &Tolerance,
    skip_hole_hits: bool,
) -> Result<PointClassification, BoolError> {
    let mut inside_votes = 0u32;
    let mut outside_votes = 0u32;

    let ray_trace = std::env::var_os("PK_RAY_TRACE").is_some();
    for direction in &RAY_DIRECTIONS {
        if ray_trace {
            let r = ray_crossing_count(
                point,
                direction,
                solid,
                face_aabbs,
                tolerance,
                skip_hole_hits,
            );
            eprintln!(
                "  [ray] p=({:.2},{:.2},{:.2}) d=({:.2},{:.2},{:.2}) -> {r:?}",
                point.x, point.y, point.z, direction.x, direction.y, direction.z
            );
        }
        match ray_crossing_count(
            point,
            direction,
            solid,
            face_aabbs,
            tolerance,
            skip_hole_hits,
        ) {
            Ok(count) => {
                if count % 2 == 1 {
                    inside_votes += 1;
                } else {
                    outside_votes += 1;
                }

                // Early exit: if we have a clear majority from 3+ rays, stop.
                if inside_votes >= 2 || outside_votes >= 2 {
                    break;
                }
            }
            Err(()) => {
                // Inconclusive ray — try next direction.
                continue;
            }
        }
    }

    if inside_votes > outside_votes {
        Ok(PointClassification::Inside)
    } else if outside_votes > inside_votes {
        Ok(PointClassification::Outside)
    } else {
        // Tied or no conclusive rays at all — treat as on-boundary.
        // This is a heuristic; a more robust approach would add more rays.
        Ok(PointClassification::OnBoundary)
    }
}

// ---------------------------------------------------------------------------
// 5. Classify all fragments
// ---------------------------------------------------------------------------

/// Classify every face in `fragments_brep` relative to `other_brep`.
///
/// For each face, picks a sample point (surface midpoint or first vertex) and
/// calls [`point_in_solid`] against the other solid. The resulting
/// [`PointClassification`] is mapped to [`Classification`]:
///
/// - `Inside` -> `Classification::Inside`
/// - `Outside` -> `Classification::Outside`
/// - `OnBoundary` -> `Classification::OnBoundaryCooriented` (conservative default)
///
/// # Arguments
///
/// * `fragments_brep` — the BRep whose faces are being classified.
/// * `other_brep` — the solid to test containment against.
/// * `other_face_aabbs` — precomputed face AABBs for `other_brep`.
/// * `intersection_edges` — edges that lie on the intersection (used for
///   future propagation optimization; currently unused).
/// * `tolerance` — geometric tolerance for all distance and coincidence tests.
pub(crate) fn classify_fragments(
    fragments_brep: &BRep,
    other_brep: &BRep,
    other_face_aabbs: &[(FaceHandle, Aabb)],
    intersection_edges: &HashSet<EdgeHandle>,
    tolerance: &Tolerance,
) -> Result<HashMap<FaceHandle, Classification>, BoolError> {
    // Use topological propagation: classify one seed per connected component,
    // then propagate across non-intersection edges. Falls back to per-face
    // ray casting if propagation fails.
    match crate::classify_propagation::classify_with_propagation(
        fragments_brep,
        other_brep,
        other_face_aabbs,
        intersection_edges,
        tolerance,
    ) {
        Ok((result, _ray_count)) => Ok(result),
        Err(_) => {
            // Propagation failed — fall back to per-face classification
            let mut result = HashMap::new();
            for face_handle in fragments_brep.faces.keys() {
                let sample = pick_sample_point(fragments_brep, face_handle);
                let pc = point_in_solid(&sample, other_brep, other_face_aabbs, tolerance)?;
                let class = match pc {
                    PointClassification::Inside => Classification::Inside,
                    PointClassification::Outside => Classification::Outside,
                    PointClassification::OnBoundary => Classification::OnBoundaryCooriented,
                };
                result.insert(face_handle, class);
            }
            Ok(result)
        }
    }
}

/// Pick a representative sample point on a face for classification.
///
/// Tries UV-centroid sampling first (using pcurves to find a point guaranteed
/// inside the face boundary), then falls back to edge midpoint (for faces with
/// inner loops) or vertex centroid (for simple faces). Falls back to the surface
/// parameter midpoint for faces without vertices.
pub(crate) fn pick_sample_point(brep: &BRep, face_handle: FaceHandle) -> Point3 {
    let face = &brep.faces[face_handle];

    // Try UV-centroid sampling for faces with surface geometry and pcurves.
    if let Some(pt) = uv_centroid_sample(brep, face_handle) {
        return pt;
    }

    // For faces with inner loops, use the midpoint of an outer-loop edge.
    // The centroid can fall inside a hole, giving wrong classification.
    if !face.inner_loops.is_empty() {
        let loop_data = match brep.loops.get(face.outer_loop) {
            Some(l) => l,
            None => return Point3::new(0.0, 0.0, 0.0),
        };
        let start_he = loop_data.half_edge;
        if let Some(he) = brep.half_edges.get(start_he) {
            let origin = brep.vertices[he.origin].position;
            if let Some(next_he) = brep.half_edges.get(he.next) {
                let dest = brep.vertices[next_he.origin].position;
                return Point3::new(
                    (origin.x + dest.x) * 0.5,
                    (origin.y + dest.y) * 0.5,
                    (origin.z + dest.z) * 0.5,
                );
            }
        }
    }

    // Use vertex centroid — reliable for convex polygonal faces without holes.
    let positions = face_vertex_positions(brep, face_handle);
    if !positions.is_empty() {
        let n = positions.len() as f64;
        let sum = positions.iter().fold(Vec3::new(0.0, 0.0, 0.0), |acc, p| {
            Vec3::new(acc.x + p.x, acc.y + p.y, acc.z + p.z)
        });
        return Point3::new(sum.x / n, sum.y / n, sum.z / n);
    }

    // Fall back to surface parameter midpoint.
    let face = &brep.faces[face_handle];
    if let Some(surface_on_face) = face.surface.as_ref() {
        let ((u_min, u_max), (v_min, v_max)) = face_parameter_range(brep, face_handle);
        let u_mid = (u_min + u_max) * 0.5;
        let v_mid = (v_min + v_max) * 0.5;
        return surface_on_face.geometry.point_at(u_mid, v_mid);
    }

    // Last resort — origin.
    Point3::new(0.0, 0.0, 0.0)
}

#[path = "classify_uv_sample.rs"]
mod uv_sample;
use uv_sample::uv_centroid_sample;

#[cfg(test)]
#[path = "classify_tests.rs"]
mod tests;

#[cfg(test)]
#[path = "classify_edge_tests.rs"]
mod edge_tests;
