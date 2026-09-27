//! Seam bookkeeping for the closed-fillet surgery: seam-edge lookup,
//! seam-consistent contact anchors, and seam-curve rebuild. Split out
//! of `closed_fillet.rs` (500-line ratchet).

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::topology::{BRep, EdgeHandle, HalfEdgeHandle, LoopHandle, VertexHandle};

/// Seam-consistent anchor on a contact circle (PK-285): the point on
/// the exact contact circle radially nearest `v_old_pos` — NOT the
/// trim curve's own start point. For a closed rim the trim NURBS
/// inherits its parameterization origin from the offset-surface
/// intersection marcher's seed, which is an arbitrary, ULP-sensitive
/// artifact (measured: a 7e-15 extrude-height difference moved the
/// anchor from the seam angle to the diametrically opposite side, and
/// the rebuilt seam line crossed the cylinder as a chord — check-9
/// caught it at deviation = radius). Anchoring at the old rim-seam
/// vertex keeps the shortened seam edge a generator of its surface.
pub(super) fn seam_consistent_anchor(
    circle: &Option<(Point3, Vec3, f64)>,
    trim: &popkern_core::geometry::NurbsCurve,
    v_old_pos: Point3,
) -> Point3 {
    if let Some((center, normal, radius)) = circle {
        let d = v_old_pos - *center;
        let radial = d - *normal * d.dot(normal);
        if let Some(dir) = radial.normalized() {
            return *center + dir * *radius;
        }
    }
    // No exact circle (non-circular closed trim): closest sampled
    // point on the trim curve to the old vertex — same
    // seam-consistency contract, sampled.
    let (t0, t1) = trim.parameter_range();
    let mut best = trim.point_at(t0);
    let mut best_dist = f64::INFINITY;
    for i in 0..=64 {
        let t = t0 + (t1 - t0) * (i as f64) / 64.0;
        let p = trim.point_at(t);
        let dist = p.distance_to(&v_old_pos);
        if dist < best_dist {
            best_dist = dist;
            best = p;
        }
    }
    best
}

/// Find a half-edge in `loop_handle` other than `skip` whose
/// `origin` is `vertex` — i.e. an edge incident to the rim vertex
/// from a non-rim direction. For the cylinder case this returns
/// the seam edge half-edge.
pub(super) fn find_seam_edge_at_vertex(
    brep: &BRep,
    loop_handle: LoopHandle,
    skip: HalfEdgeHandle,
    vertex: VertexHandle,
) -> Option<EdgeHandle> {
    let start = brep.loops.get(loop_handle)?.half_edge;
    let mut cursor = start;
    let mut iters = 0;
    loop {
        iters += 1;
        if iters > 10_000 {
            return None;
        }
        if cursor != skip {
            let h = brep.half_edges.get(cursor)?;
            if h.origin == vertex {
                return Some(h.edge);
            }
        }
        let next = brep.half_edges.get(cursor)?.next;
        cursor = next;
        if cursor == start {
            break;
        }
    }
    None
}

/// Re-window a seam edge's curve to its vertices' current positions
/// after one endpoint moved to the new contact-circle height. The
/// canonical "forward" direction is preserved: the curve starts at
/// `edge.half_edge.origin` and ends at the twin's `origin`.
///
/// A straight seam (cylinder, cone) is rebuilt as the line between
/// the vertices. A curved seam — a sphere's meridian, a torus's tube
/// circle — keeps its geometry and takes the parameter window whose
/// ends project to the moved vertices: replacing it with the chord
/// (the old behaviour for every form) left a capped ball's seam 0.29
/// off its sphere after the cap-rim fillet (PK-418). When either
/// vertex does not lie on the curve within the body's tolerance the
/// chord is authored as before, and the validator reports it.
pub(super) fn rebuild_seam_curve(brep: &mut BRep, edge: EdgeHandle) {
    use popkern_core::geometry::curve::{CurveGeometry, CurveOnEdge};

    let Some(edata) = brep.edges.get(edge) else {
        return;
    };
    let he = edata.half_edge;
    let Some(he_data) = brep.half_edges.get(he).cloned() else {
        return;
    };
    let Some(twin_data) = brep.half_edges.get(he_data.twin).cloned() else {
        return;
    };

    let v_start_pos = brep.vertices[he_data.origin].position;
    let v_end_pos = brep.vertices[twin_data.origin].position;
    let direction = Vec3::new(
        v_end_pos.x - v_start_pos.x,
        v_end_pos.y - v_start_pos.y,
        v_end_pos.z - v_start_pos.z,
    );
    if direction.length() < 1e-12 {
        return;
    }

    if let Some(curve) = edata.curve.clone()
        && !matches!(curve.geometry, CurveGeometry::Line { .. })
        && let Some(windowed) = rewindowed(&curve, v_start_pos, v_end_pos, brep.tolerance.linear())
    {
        brep.edges[edge].curve = Some(windowed);
        return;
    }

    brep.edges[edge].curve = Some(CurveOnEdge {
        geometry: CurveGeometry::Line {
            origin: v_start_pos,
            direction,
        },
        t_start: 0.0,
        t_end: 1.0,
    });
}

/// `curve` over the parameter window running from the point nearest
/// `start` to the point nearest `end`, when both lie on the curve
/// within `gate` and the window keeps the curve's direction.
fn rewindowed(
    curve: &popkern_core::geometry::curve::CurveOnEdge,
    start: Point3,
    end: Point3,
    gate: f64,
) -> Option<popkern_core::geometry::curve::CurveOnEdge> {
    let eval: &dyn CurveEval = &curve.geometry;
    let t_a = eval.closest_point(&start)?;
    let t_b = eval.closest_point(&end)?;
    if eval.point_at(t_a).distance_to(&start) > gate || eval.point_at(t_b).distance_to(&end) > gate
    {
        return None;
    }
    // The old window ran start → end; the new one must too.
    let forward = (curve.t_end - curve.t_start).signum();
    if (t_b - t_a).abs() < 1e-12 || (t_b - t_a).signum() != forward {
        return None;
    }
    Some(popkern_core::geometry::curve::CurveOnEdge {
        geometry: curve.geometry.clone(),
        t_start: t_a,
        t_end: t_b,
    })
}

/// Author the strip's transverse seam curve between the two contact
/// anchors `p_out` (outer-loop vertex) and `p_in` (inner-ring
/// vertex).
///
/// Exact-first (PK-295): on a torus strip the seam is a tube
/// cross-section — an exact circular arc in the plane spanned by the
/// surface axis and the station's radial direction. The 17-point
/// polyline fallback carries ~2.4e-4 chord sag at unit scales, which
/// is 100× the boolean's 2e-6 junction snap: an edge-face hit
/// against the polyline lands off every true intersection curve, the
/// pierce node never joins the record graph, and the fragments
/// flanking the seam cannot close (measured: the θ=0 torus bulge of
/// a rim-filleted cylinder minus a cuboid dropped out of the
/// subtraction, stranding three self-twinned half-edges).
///
/// Non-toroidal strips keep the UV-interpolated polyline (branch
/// alignment per PK-285), authored by the caller.
pub(super) fn exact_torus_seam_arc(
    surface: &popkern_core::geometry::surface::SurfaceGeometry,
    p_out: Point3,
    p_in: Point3,
    tolerance: &popkern_core::tolerance::Tolerance,
) -> Option<popkern_core::geometry::CurveOnEdge> {
    use popkern_core::geometry::curve::CurveGeometry;
    let popkern_core::geometry::surface::SurfaceGeometry::Torus {
        center,
        axis,
        major_radius,
        minor_radius,
    } = surface
    else {
        return None;
    };
    let d = p_out - *center;
    let h = d.dot(axis);
    let radial = d - *axis * h;
    let radial_len = radial.length();
    if radial_len < 1e-12 {
        return None;
    }
    let radial_dir = radial * (1.0 / radial_len);
    let c = *center + radial_dir * *major_radius;
    // Tube-circle normal: the center circle's tangent at this station.
    let normal = axis.cross(&radial_dir);
    // The Arc evaluator walks ASCENDING angles in its own
    // `circle_basis` frame; pick whichever normal orientation makes
    // the p_out → p_in sweep positive, recomputing angles in that
    // frame (a mirrored frame is NOT a sign flip of the original
    // basis — circle_basis chooses its own e1).
    let mut chosen: Option<(Vec3, f64, f64)> = None;
    for n in [normal, -normal] {
        let (e1, e2) = popkern_core::geometry::analytic_curves::circle_basis(n);
        let ang = |p: Point3| -> f64 {
            let q = p - c;
            q.dot(&e2).atan2(q.dot(&e1))
        };
        let a0 = ang(p_out);
        let mut sweep = ang(p_in) - a0;
        while sweep > std::f64::consts::PI {
            sweep -= std::f64::consts::TAU;
        }
        while sweep < -std::f64::consts::PI {
            sweep += std::f64::consts::TAU;
        }
        if sweep > 1e-9 {
            chosen = Some((n, a0, sweep));
            break;
        }
    }
    let (n, a0, sweep) = chosen?;
    let arc = popkern_core::geometry::CurveOnEdge {
        geometry: CurveGeometry::Arc {
            center: c,
            normal: n,
            radius: *minor_radius,
            start_angle: a0,
            end_angle: a0 + sweep,
        },
        t_start: a0,
        t_end: a0 + sweep,
    };
    // Endpoint honesty gate: the arc must START at p_out and END at
    // p_in within the authoring band, else decline to the fallback.
    let pos_tol = (tolerance.linear() * 10.0).max(1e-9);
    let s = (&arc.geometry as &dyn CurveEval).point_at(arc.t_start);
    let e = (&arc.geometry as &dyn CurveEval).point_at(arc.t_end);
    if s.distance_to(&p_out) > pos_tol || e.distance_to(&p_in) > pos_tol {
        return None;
    }
    Some(arc)
}

/// Author the strip's transverse seam curve: exact torus tube arc
/// first ([`exact_torus_seam_arc`]), else the 17-point
/// UV-interpolated polyline with PK-285 period-branch alignment
/// (moved verbatim from `closed_fillet.rs` step 14).
pub(super) fn author_strip_seam_curve(
    surface: &popkern_core::geometry::surface::SurfaceGeometry,
    p_out: Point3,
    p_in: Point3,
    tolerance: &popkern_core::tolerance::Tolerance,
) -> Option<popkern_core::geometry::CurveOnEdge> {
    if let Some(arc) = exact_torus_seam_arc(surface, p_out, p_in, tolerance) {
        return Some(arc);
    }
    use popkern_core::geometry::surface::SurfaceEval;
    let uv_of = |p: &Point3| popkern_core::geometry::uv_projection::project_point_to_uv(surface, p);
    let (uv_a, mut uv_b) = (uv_of(&p_out)?, uv_of(&p_in)?);
    // The two projections can land on DIFFERENT period branches (one
    // endpoint at u = 2π, the other at u = 0): interpolating across
    // that gap sweeps the "cross-section" a full lap around the
    // chart, and the seam's pcurves become a diagonal through every
    // u (measured: OCC flagged the strip wire self-intersecting).
    // Align u_b to u_a's branch — the seam is a single station.
    if let Some(period) = surface.parameter_periodicity().0 {
        uv_b.0 -= period * ((uv_b.0 - uv_a.0) / period).round();
    }
    // The v direction needs the same branch discipline (PK-285): on
    // a torus strip the outer contact point sits at tube angle
    // v = 0 ≡ 2π, and atan2(±ε, ·) puts the raw projection on either
    // side of the wrap by ULP luck. Unaligned, the seam interpolates
    // the LONG way around the tube (measured: 270° instead of 90°;
    // its 17-point polyline chords then sag 2.2e-2 into the tube,
    // the boundary-first tessellator rejects the ring and the strip
    // cracks open along both contact circles).
    if let Some(period) = surface.parameter_periodicity().1 {
        uv_b.1 -= period * ((uv_b.1 - uv_a.1) / period).round();
    }
    let n = 17usize;
    let mut pts: Vec<Point3> = (0..n)
        .map(|i| {
            let f = i as f64 / (n - 1) as f64;
            let u = uv_a.0 + (uv_b.0 - uv_a.0) * f;
            let v = uv_a.1 + (uv_b.1 - uv_a.1) * f;
            surface.point_at(u, v)
        })
        .collect();
    // Endpoints exactly at the vertices.
    pts[0] = p_out;
    pts[n - 1] = p_in;
    let mut knots = Vec::with_capacity(n + 2);
    knots.push(0.0);
    for i in 0..n {
        knots.push(i as f64 / (n - 1) as f64);
    }
    knots.push(1.0);
    Some(popkern_core::geometry::CurveOnEdge {
        geometry: popkern_core::geometry::CurveGeometry::Nurbs(
            popkern_core::geometry::NurbsCurve {
                degree: 1,
                knots,
                control_points: pts,
                weights: vec![1.0; n],
            },
        ),
        t_start: 0.0,
        t_end: 1.0,
    })
}

/// Orient a contact circle's axis so its angle parameterization runs
/// the way the replaced rim half-edge `old_he` walked its support
/// face.
///
/// The new contact-circle edge's canonical half-edge sits on the
/// support face in place of `old_he` (the rim shrinks onto a coaxial
/// circle; the loop keeps its direction), and it traverses the Arc in
/// increasing angle about the stored axis. The axis the trim fit
/// hands over has an arbitrary sign, so without this the support loop
/// could wind backwards; the serving-boundary flag reconcile then
/// "honestly" flipped the face — the chamfered shank tip and its wall
/// served inside-out whenever the input rim was outward (measured on
/// the M6 bolt tip once the head ∪ shank union stopped inverting the
/// shank). The direction authority is the old half-edge's own pcurve
/// trace on its face (what the canonical-forward audit reads); with
/// no pcurve the axis is kept.
pub(super) fn orient_circle_to_rim(
    brep: &BRep,
    old_he: HalfEdgeHandle,
    circle: Option<(Point3, Vec3, f64)>,
) -> Option<(Point3, Vec3, f64)> {
    let (center, normal, radius) = circle?;
    let Some(axis) = half_edge_winding_axis(brep, old_he) else {
        return circle;
    };
    if axis.dot(&normal) < 0.0 {
        Some((center, normal * -1.0, radius))
    } else {
        circle
    }
}

/// Newell axis of a closed half-edge's traversal, from its pcurve
/// mapped through its face's surface (right-hand rule: the loop runs
/// counter-clockwise about the returned axis). `None` without a
/// pcurve or surface, or when the trace is degenerate.
fn half_edge_winding_axis(brep: &BRep, he: HalfEdgeHandle) -> Option<Vec3> {
    use popkern_core::geometry::SurfaceEval as _;
    const SAMPLES: usize = 16;
    let h = brep.half_edges.get(he)?;
    let pc = h.pcurve.as_ref()?;
    let face = brep.loops.get(h.loop_)?.face;
    let surface = &brep.faces.get(face)?.surface.as_ref()?.geometry;
    let (t0, t1) = pc.parameter_range();
    let pts: Vec<Point3> = (0..SAMPLES)
        .map(|k| {
            let uv = pc.point_at(t0 + (t1 - t0) * k as f64 / SAMPLES as f64);
            surface.point_at(uv.x, uv.y)
        })
        .collect();
    let mut n = Vec3::new(0.0, 0.0, 0.0);
    for (i, a) in pts.iter().enumerate() {
        let b = pts[(i + 1) % SAMPLES];
        n = n + Vec3::new(
            (a.y - b.y) * (a.z + b.z),
            (a.z - b.z) * (a.x + b.x),
            (a.x - b.x) * (a.y + b.y),
        );
    }
    n.normalized()
}
