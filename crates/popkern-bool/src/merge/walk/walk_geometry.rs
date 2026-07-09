//! Geometry helpers for loop walking.
//!
//! Angular comparison, tangent computation, and surface normal lookup
//! used by the walk module to select the correct next half-edge at
//! intersection nodes.

use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::geometry::{CurveEval, Vec3};
use popkern_core::topology::BRep;

use crate::merge::graph::{FaceRef, IntersectionGraph, Operand};

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
    let twin_idx = graph.half_edges[he_idx].twin;
    let end_node = graph.half_edges[twin_idx].origin;
    let node_pos = graph.nodes[end_node].position;

    let origin_node = graph.half_edges[he_idx].origin;
    // Incoming direction: the TRUE curve tangent at the arrival node
    // when the edge carries a curve (chords of curved edges point
    // across the geometry — two arcs of a split intersection circle
    // share identical chords and become indistinguishable). Self-loop
    // handling falls out of the same rule.
    let incoming_dir = if origin_node == end_node {
        incoming_tangent_for_self_loop(graph, he_idx)
    } else if let Some(curve) = &graph.half_edges[he_idx].curve {
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
    let normal = planar_walk_normal(face, brep_a, brep_b)
        .or_else(|| surface_normal_at_node(graph, end_node, face, brep_a, brep_b))
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
        } else if let Some(curve) = &graph.half_edges[c].curve {
            // True departure tangent for curved candidates.
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
                let a = sin_val.atan2(cos_val);
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

        if std::env::var_os("PK_WALK_TRACE").is_some() {
            eprintln!("      [angle] cand {c}: angle={angle:.3e}");
        }
        if angle < best_angle {
            best_angle = angle;
            best = Some(c);
        }
    }

    best
}

/// Walk normal for a PLANAR face, derived from the face's stored
/// outer-loop winding via Newell's method — NEGATED, so it points OUT
/// of the material (this kernel's loop Newell points inward).
///
/// This is the walk-consistent handedness: fragment loops reuse the
/// face's boundary half-edges in their registered directions, so every
/// valid loop winds counter-clockwise about the inward Newell vector.
/// `pick_by_angle` measures candidates counter-clockwise from the
/// BACKWARD direction of arrival, and the tightest fragment on the
/// walk's material side is the candidate with the smallest CLOCKWISE
/// angle about the inward vector — which is exactly the smallest
/// counter-clockwise angle about this outward vector, matching the
/// existing comparison (verified against four seam junctions of a
/// flush coplanar union plus the convex-corner base case). Returns
/// `None` for non-planar faces — curved charts keep the tuned
/// UV-normal path.
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

/// Compute a signed area proxy for a loop using Newell's method.
///
/// Positive = CCW (outer), negative = CW (inner/hole).
pub(super) fn signed_area(graph: &IntersectionGraph, loop_hes: &[usize]) -> f64 {
    let normal = loop_newell(graph, loop_hes);
    normal.length() * normal.z.signum()
}

/// The (borrow-signed) Newell vector of a loop: |N| = 2·area, and
/// the DIRECTION encodes winding — two coplanar loops wind opposite
/// exactly when their Newell vectors point opposite ways. This is
/// the orientation-safe comparison (the historical `sign(N.z)` trick
/// degenerates on vertical faces).
pub(crate) fn loop_newell(graph: &IntersectionGraph, loop_hes: &[usize]) -> Vec3 {
    // Sum cross products of consecutive positions projected onto an
    // approximate normal. CURVED half-edges are SAMPLED along their
    // span — origin-only polygons measure a curved loop as its
    // chord polygon, and an ellipse-fronted face lost the outer-loop
    // election to its own seat ring because the chord quad
    // under-measured the true area (the bulge is ~30% of the face).
    let mut positions: Vec<popkern_core::geometry::Point3> = Vec::new();
    for &he_idx in loop_hes {
        let he = &graph.half_edges[he_idx];
        let origin_pos = graph.nodes[he.origin].position;
        let own = he.curve.is_some();
        let curve = he
            .curve
            .as_ref()
            .or_else(|| graph.half_edges[he.twin].curve.as_ref());
        let curved = curve
            .map(|c| match &c.geometry {
                popkern_core::geometry::curve::CurveGeometry::Line { .. } => false,
                popkern_core::geometry::curve::CurveGeometry::Nurbs(n) => n.degree > 1,
                _ => true,
            })
            .unwrap_or(false);
        if curved {
            use popkern_core::geometry::curve::CurveEval;
            let c = curve.expect("curved implies curve");
            let (t0, t1) = (c.t_start, c.t_end);
            // Walk the samples in THIS half-edge's traversal
            // direction: a borrowed (twin) curve runs opposite.
            let n_s = 8;
            for i in 0..n_s {
                let f = i as f64 / n_s as f64;
                let t = if own {
                    t0 + (t1 - t0) * f
                } else {
                    t1 - (t1 - t0) * f
                };
                positions.push(c.geometry.point_at(t));
            }
        } else {
            positions.push(origin_pos);
        }
    }

    // Check if all positions are coincident (self-loop structure).
    // If so, sample curve geometry to compute area.
    let all_coincident = positions.len() >= 2
        && positions
            .iter()
            .all(|p| p.distance_to(&positions[0]) < 1e-10);

    if positions.len() < 3 || all_coincident {
        // Self-loop edge (closed curve like a circle): estimate area from
        // the curve geometry by sampling points along the curve.
        // Try each half-edge in the loop until one has curve geometry.
        for &he_idx in loop_hes {
            let he = &graph.half_edges[he_idx];
            // Curve geometry may live on either member of the pair —
            // a rim half-edge without its own curve measured area 0,
            // lost the outer-loop election to its 12-segment socket
            // hole, and the cap fragment came out inside-out. When
            // the curve is BORROWED from the twin, the samples run
            // opposite this loop's traversal, so the SIGN of the
            // area must flip (hole/piece role assignment reads the
            // sign; an unflipped borrow made cap inner loops wind
            // with their outer and reopened the socket seams).
            let own = he.curve.is_some();
            let curve = he
                .curve
                .as_ref()
                .or_else(|| graph.half_edges[he.twin].curve.as_ref());
            let borrow_sign = if own { 1.0 } else { -1.0 };
            if let Some(curve) = curve {
                use popkern_core::geometry::curve::CurveEval;
                let n_samples = 16;
                let (t0, t1) = (curve.t_start, curve.t_end);
                let mut pts = Vec::with_capacity(n_samples);
                for i in 0..n_samples {
                    let t = t0 + (t1 - t0) * (i as f64) / (n_samples as f64);
                    pts.push(curve.geometry.point_at(t));
                }
                // Newell's method on sampled points. Scale matches the
                // polygon branch below (|N| = 2·area) so closed rims
                // and polygonal loops compete on one measure.
                let mut nx = 0.0;
                let mut ny = 0.0;
                let mut nz = 0.0;
                for i in 0..pts.len() {
                    let cur = &pts[i];
                    let next = &pts[(i + 1) % pts.len()];
                    nx += (cur.y - next.y) * (cur.z + next.z);
                    ny += (cur.z - next.z) * (cur.x + next.x);
                    nz += (cur.x - next.x) * (cur.y + next.y);
                }
                let normal = Vec3::new(nx, ny, nz);
                return normal * borrow_sign;
            }
        }
        return Vec3::new(0.0, 0.0, 0.0);
    }

    // Newell's method for 3D polygon normal (proportional to area).
    let mut nx = 0.0;
    let mut ny = 0.0;
    let mut nz = 0.0;

    let n = positions.len();
    for i in 0..n {
        let cur = &positions[i];
        let next = &positions[(i + 1) % n];
        nx += (cur.y - next.y) * (cur.z + next.z);
        ny += (cur.z - next.z) * (cur.x + next.x);
        nz += (cur.x - next.x) * (cur.y + next.y);
    }

    // The length of the vector is proportional to twice the area;
    // the direction encodes winding.
    Vec3::new(nx, ny, nz)
}

/// True when `lp`'s nodes are (by strict majority) OUTSIDE the outer
/// loop's sampled polygon — such a loop is a disjoint region, never
/// a hole of this fragment. Points within a small band of the
/// boundary are skipped (loops may share junction nodes with the
/// outer ring). Projection: the face plane's two dominant axes from
/// the outer Newell vector.
pub(super) fn loop_outside_outer(
    graph: &IntersectionGraph,
    lp: &[usize],
    outer: &[usize],
    outer_newell: &Vec3,
) -> bool {
    // 2D chart: drop the dominant normal component.
    let n = *outer_newell;
    let (ax, ay) = if n.x.abs() >= n.y.abs() && n.x.abs() >= n.z.abs() {
        (1, 2)
    } else if n.y.abs() >= n.z.abs() {
        (0, 2)
    } else {
        (0, 1)
    };
    let coord = |p: &popkern_core::geometry::Point3, a: usize| match a {
        0 => p.x,
        1 => p.y,
        _ => p.z,
    };
    // Sampled outer polygon (curved half-edges sampled like
    // loop_newell does, via node origins + curve stations).
    let mut ring: Vec<(f64, f64)> = Vec::new();
    for &he_idx in outer {
        let he = &graph.half_edges[he_idx];
        let own = he.curve.is_some();
        let curve = he
            .curve
            .as_ref()
            .or_else(|| graph.half_edges[he.twin].curve.as_ref());
        let curved = curve
            .map(|c| match &c.geometry {
                popkern_core::geometry::curve::CurveGeometry::Line { .. } => false,
                popkern_core::geometry::curve::CurveGeometry::Nurbs(nc) => nc.degree > 1,
                _ => true,
            })
            .unwrap_or(false);
        if curved {
            use popkern_core::geometry::curve::CurveEval;
            let c = curve.expect("curved implies curve");
            let (t0, t1) = (c.t_start, c.t_end);
            for i in 0..8 {
                let f = i as f64 / 8.0;
                let t = if own {
                    t0 + (t1 - t0) * f
                } else {
                    t1 - (t1 - t0) * f
                };
                let q = c.geometry.point_at(t);
                ring.push((coord(&q, ax), coord(&q, ay)));
            }
        } else {
            let p = graph.nodes[he.origin].position;
            ring.push((coord(&p, ax), coord(&p, ay)));
        }
    }
    if ring.len() < 3 {
        return false;
    }
    // Scale-aware boundary band.
    let (mut lo, mut hi) = ((f64::MAX, f64::MAX), (f64::MIN, f64::MIN));
    for &(x, y) in &ring {
        lo = (lo.0.min(x), lo.1.min(y));
        hi = (hi.0.max(x), hi.1.max(y));
    }
    let band = 1.0e-3 * ((hi.0 - lo.0) + (hi.1 - lo.1)).max(1e-9);
    let mut inside_n = 0usize;
    let mut outside_n = 0usize;
    for &he_idx in lp {
        let p = graph.nodes[graph.half_edges[he_idx].origin].position;
        let (px, py) = (coord(&p, ax), coord(&p, ay));
        // Distance to the ring: skip near-boundary points.
        let mut min_d = f64::MAX;
        for i in 0..ring.len() {
            let (x0, y0) = ring[i];
            let (x1, y1) = ring[(i + 1) % ring.len()];
            let (dx, dy) = (x1 - x0, y1 - y0);
            let l2 = (dx * dx + dy * dy).max(1e-300);
            let t = (((px - x0) * dx + (py - y0) * dy) / l2).clamp(0.0, 1.0);
            let d = ((px - x0 - dx * t).powi(2) + (py - y0 - dy * t).powi(2)).sqrt();
            min_d = min_d.min(d);
        }
        if min_d < band {
            continue;
        }
        let mut crossings = 0;
        for i in 0..ring.len() {
            let (x0, y0) = ring[i];
            let (x1, y1) = ring[(i + 1) % ring.len()];
            if (y0 > py) != (y1 > py) {
                let xc = x0 + (py - y0) / (y1 - y0) * (x1 - x0);
                if xc > px {
                    crossings += 1;
                }
            }
        }
        if crossings % 2 == 1 {
            inside_n += 1;
        } else {
            outside_n += 1;
        }
    }
    // Require an UNAMBIGUOUS verdict: every tested node outside and
    // enough of them to trust. A mixed or thin verdict keeps the
    // existing size-rule behavior (measured: a through-slot's ring
    // on a long tapered side face produced noisy outside votes and
    // the demotion blocked the tunnel — χ' 0 → 2).
    inside_n == 0 && outside_n >= 3
}
