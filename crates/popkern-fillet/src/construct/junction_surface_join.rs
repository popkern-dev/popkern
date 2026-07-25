//! Joint strip surface across smooth-through junctions (PK-177).
//!
//! At a smooth-through junction the chain continues over ONE
//! envelope, but the two pieces' degree-1 lofts interpolate chords
//! between DIFFERENT station sets — their mutual gap along the merge
//! boundary is the chord-vs-chord sagitta (~s²/8R, measured 6.6e-5
//! on the boss cross-seam boundary), irreducible by any curve fit
//! against two disagreeing surfaces. The honest object is a SINGLE
//! loft over the union of both pieces' stations.
//!
//! FULL-PASS contract (the surgical variant was falsified: a local
//! junction-gap acceptance hid served curves at 2.09e-2 off the
//! replacement): after the joint surface installs on both faces,
//! every edge curve bounding either face REFITS against it
//! (reconcile + same-face projection), and acceptance is GLOBAL —
//! the worst witness-measured curve-to-surface deviation over ALL
//! affected edges must strictly improve. Transactional.

use popkern_core::geometry::surface::{SurfaceGeometry, SurfaceOnFace};
use popkern_core::geometry::{NurbsSurface, Point3, SurfaceEval, Vec3};
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle};

/// Circumcenter of three points (the station-arc centre when the
/// points sample a circle). `None` for collinear samples.
fn circumcenter(p0: &Point3, pm: &Point3, p2: &Point3) -> Option<Point3> {
    let a = *pm - *p0;
    let b = *p2 - *p0;
    let n = a.cross(&b);
    let n2 = n.dot(&n);
    if n2 < 1e-24 {
        return None;
    }
    let num = b.cross(&n) * a.dot(&a) + n.cross(&a) * b.dot(&b);
    Some(*p0 + num * (0.5 / n2))
}

/// Per-station (rail0, arc-mid, rail1) triples from a strip loft,
/// layout-independently.
fn extract_station_triples(surface: &NurbsSurface) -> Vec<(Point3, Point3, Point3)> {
    let geometry = SurfaceGeometry::Nurbs(surface.clone());
    let n = surface.num_u.max(1);
    let ((u0, u1), (v0, v1)) = geometry.parameter_range();
    (0..n)
        .map(|i| {
            let u = u0 + (u1 - u0) * (i as f64) / ((n - 1).max(1) as f64);
            (
                geometry.point_at(u, v0),
                geometry.point_at(u, 0.5 * (v0 + v1)),
                geometry.point_at(u, v1),
            )
        })
        .collect()
}

/// Edges bounding a face (outer loop).
fn face_edges(brep: &BRep, face: FaceHandle) -> Vec<EdgeHandle> {
    let mut edges = Vec::new();
    let Some(fd) = brep.faces.get(face) else {
        return edges;
    };
    let Some(l) = brep.loops.get(fd.outer_loop) else {
        return edges;
    };
    let start = l.half_edge;
    let mut cur = start;
    for _ in 0..brep.half_edges.len() {
        let Some(he) = brep.half_edges.get(cur) else {
            break;
        };
        if !edges.contains(&he.edge) {
            edges.push(he.edge);
        }
        cur = he.next;
        if cur == start {
            break;
        }
    }
    edges
}

/// Worst witness-measured curve-to-adjacent-surface deviation over
/// the given edges. `label` names the phase in the gated trace.
fn global_deviation(brep: &BRep, edges: &[EdgeHandle], label: &str) -> f64 {
    let trace = std::env::var_os("PK_RECONCILE_TRACE").is_some();
    use popkern_core::geometry::curve::CurveEval;
    let mut worst = 0.0_f64;
    for &edge in edges {
        let Some(ed) = brep.edges.get(edge) else {
            continue;
        };
        let Some(curve) = ed.curve.as_ref() else {
            continue;
        };
        let Some(he) = brep.half_edges.get(ed.half_edge) else {
            continue;
        };
        let Some(tw) = brep.half_edges.get(he.twin) else {
            continue;
        };
        for face in [
            brep.loops.get(he.loop_).map(|l| l.face),
            brep.loops.get(tw.loop_).map(|l| l.face),
        ]
        .into_iter()
        .flatten()
        {
            let Some(sof) = brep.faces.get(face).and_then(|f| f.surface.as_ref()) else {
                continue;
            };
            let mut prev: Option<(f64, f64)> = None;
            for i in 0..=8 {
                let t = curve.t_start + (curve.t_end - curve.t_start) * (i as f64) / 8.0;
                let pt = curve.geometry.point_at(t);
                let cold =
                    popkern_core::geometry::uv_projection::project_point_to_uv(&sof.geometry, &pt);
                let seeded = prev.and_then(|(pu, pv)| {
                    popkern_core::geometry::uv_projection::project_point_to_uv_near(
                        &sof.geometry,
                        &pt,
                        pu,
                        pv,
                    )
                });
                if let Some((d, uv)) = [cold, seeded]
                    .into_iter()
                    .flatten()
                    .map(|(u, v)| (pt.distance_to(&sof.geometry.point_at(u, v)), (u, v)))
                    .min_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal))
                {
                    prev = Some(uv);
                    if d > worst {
                        worst = d;
                    }
                    if trace && d > 1.0e-4 {
                        eprintln!(
                            "[junction-join]   {label} {edge:?} on {face:?}: d={d:.3e} at i={i}"
                        );
                    }
                }
            }
        }
    }
    worst
}

pub(crate) struct SmoothJunctionJoin<'a> {
    pub surface_a: &'a SurfaceGeometry,
    pub surface_b: &'a SurfaceGeometry,
    pub v_pos: Point3,
    pub chain_direction: Vec3,
    pub radius: f64,
    pub trace: bool,
}

/// Join the two strip faces meeting at a smooth junction onto one
/// loft, refit their boundary curves, and accept on GLOBAL
/// improvement. Returns true when the joint committed.
pub(crate) fn join_smooth_junction_strips(brep: &mut BRep, join: SmoothJunctionJoin<'_>) -> bool {
    let SmoothJunctionJoin {
        surface_a,
        surface_b,
        v_pos,
        chain_direction,
        radius,
        trace,
    } = join;
    let (SurfaceGeometry::Nurbs(loft_a), SurfaceGeometry::Nurbs(loft_b)) = (surface_a, surface_b)
    else {
        return false;
    };
    let sig_a = format!("{surface_a:?}");
    let sig_b = format!("{surface_b:?}");
    let mut face_a: Option<FaceHandle> = None;
    let mut face_b: Option<FaceHandle> = None;
    for (fh, fd) in brep.faces.iter() {
        let Some(sof) = fd.surface.as_ref() else {
            continue;
        };
        let sig = format!("{:?}", sof.geometry);
        if sig == sig_a {
            face_a = Some(fh);
        } else if sig == sig_b {
            face_b = Some(fh);
        }
    }
    let (Some(face_a), Some(face_b)) = (face_a, face_b) else {
        return false;
    };
    // Union of station rows, ordered along the chain, deduplicated.
    let mut rows: Vec<(f64, [Point3; 5], f64)> = Vec::new();
    for triples in [
        extract_station_triples(loft_a),
        extract_station_triples(loft_b),
    ] {
        for (p0, pm, p2) in triples {
            let Some(centre) = circumcenter(&p0, &pm, &p2) else {
                continue;
            };
            if (centre.distance_to(&p0) - radius).abs() > radius * 0.05 {
                continue;
            }
            let Some((cps, w)) = crate::station_geometry::bi_rational_arc_row(&centre, &p0, &p2)
            else {
                continue;
            };
            if cps[2].distance_to(&pm) > radius * 0.2 {
                continue;
            }
            rows.push(((p0 - v_pos).dot(&chain_direction), cps, w));
        }
    }
    if rows.len() < 4 {
        return false;
    }
    rows.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
    let min_spacing = radius * 1e-3;
    let mut kept: Vec<([Point3; 5], f64)> = Vec::new();
    let mut last_s = f64::NEG_INFINITY;
    for (s, mut cps, w) in rows {
        if s - last_s < min_spacing {
            continue;
        }
        last_s = s;
        // Rail-orientation continuity: the two pieces' v-directions
        // can oppose across the junction — a mixed loft twists
        // (measured: global deviation 2.07e-2 -> 2.4e-2 without the
        // flip). Each row's rail 0 must continue the previous row's
        // rail 0.
        if let Some((prev, _)) = kept.last()
            && cps[0].distance_to(&prev[0]) > cps[4].distance_to(&prev[0])
        {
            cps.reverse();
        }
        kept.push((cps, w));
    }
    if kept.len() < 4 {
        return false;
    }
    let n = kept.len();
    let mut control_points = Vec::with_capacity(n * 5);
    let mut weights = Vec::with_capacity(n * 5);
    for (cps, w) in &kept {
        control_points.extend_from_slice(cps);
        weights.extend_from_slice(&[1.0, *w, 1.0, *w, 1.0]);
    }
    let joint = SurfaceGeometry::Nurbs(NurbsSurface {
        degree_u: 1,
        degree_v: 2,
        knots_u: crate::fillet_surface::build_clamped_uniform_knots(n, 1),
        knots_v: vec![0.0, 0.0, 0.0, 0.5, 0.5, 1.0, 1.0, 1.0],
        control_points,
        num_u: n,
        num_v: 5,
        weights,
    });
    // GLOBAL acceptance scope: every edge bounding either face,
    // EXCLUDING edges that touch a corner-class vertex (3+ distinct
    // incident surfaces): their endpoint deviations are the
    // recorded carve-pending debt no junction join can improve —
    // with them in scope the before-measure saturates at 2.07e-2
    // and drowns the join's signal (measured).
    let corner_class = |v: popkern_core::topology::VertexHandle| -> bool {
        let mut sigs: Vec<String> = Vec::new();
        for (_, he) in brep.half_edges.iter() {
            if he.origin != v {
                continue;
            }
            if let Some(l) = brep.loops.get(he.loop_)
                && let Some(fd) = brep.faces.get(l.face)
                && let Some(sof) = fd.surface.as_ref()
            {
                let sig = format!("{:?}", sof.geometry);
                if !sigs.contains(&sig) {
                    sigs.push(sig);
                }
            }
        }
        sigs.len() >= 3
    };
    let mut scope = face_edges(brep, face_a);
    for e in face_edges(brep, face_b) {
        if !scope.contains(&e) {
            scope.push(e);
        }
    }
    scope.retain(|&e| {
        brep.edges.get(e).is_some_and(|ed| {
            brep.half_edges.get(ed.half_edge).is_some_and(|he| {
                brep.half_edges
                    .get(he.twin)
                    .is_some_and(|tw| !corner_class(he.origin) && !corner_class(tw.origin))
            })
        })
    });
    if scope.is_empty() {
        return false;
    }
    let before = global_deviation(brep, &scope, "before");
    let mut cand = brep.clone();
    for face in [face_a, face_b] {
        if let Some(fd) = cand.faces.get_mut(face) {
            let orientation = fd.surface.as_ref().map(|s| s.orientation).unwrap_or(true);
            fd.surface = Some(SurfaceOnFace {
                geometry: joint.clone(),
                orientation,
            });
        }
    }
    // Refit every scoped curve against the joint surface, then heal
    // pcurves.
    super::edge_reconcile::reconcile_edge_curves(&mut cand);
    for &edge in &scope {
        let _ = super::seam_self_fit::try_same_face_seam_fit(&mut cand, edge);
    }
    cand.fix_stale_pcurves();
    let after = global_deviation(&cand, &scope, "after");
    if trace {
        eprintln!(
            "[junction-join] {n} joint station(s), {} scoped edge(s); global deviation \
             {before:.3e} -> {after:.3e}",
            scope.len()
        );
    }
    if !after.is_finite() || after >= before || cand.validate().is_err() {
        return false;
    }
    *brep = cand;
    true
}

pub(crate) struct SmoothJunctionContext {
    pub vertex: popkern_core::topology::VertexHandle,
    pub a: usize,
    pub b: usize,
    pub v_pos: Point3,
    pub ta: Vec3,
    pub tb: Vec3,
    pub radius: f64,
    pub trace: bool,
}

/// Smooth-through junction tail: settle the junction vertex onto
/// its incident-surface intersection (corner snap fires only for
/// true corners; a stale vertex measured a 2.07e-2 join-chord
/// bulge), then optionally attempt the joint loft (opt-in
/// PK177_JOIN=1 — the global-deviation acceptance currently
/// declines because its witness measure stalls at a pre-existing
/// 2.07e-2 artifact; measurement redesign recorded in PK-177).
pub(crate) fn handle_smooth_junction(
    brep: &mut BRep,
    edge_geometries: &[crate::construct::EdgeGeometry],
    ctx: SmoothJunctionContext,
    tolerance: &popkern_core::tolerance::Tolerance,
) -> usize {
    let SmoothJunctionContext {
        vertex,
        a,
        b,
        v_pos,
        ta,
        tb,
        radius,
        trace,
    } = ctx;
    if ta.dot(&tb).abs() < 0.99 {
        return 0;
    }
    let mut absorbed = 0usize;
    if let Some(n) =
        super::corner_chord_settle::settle_vertex_on_incident_surfaces(brep, vertex, tolerance)
        && n > 0
    {
        absorbed += n;
        if trace {
            eprintln!("[junction] {vertex:?} smooth junction vertex settled");
        }
    }
    // Neighbor stale-vertex settle (mirrors the corner pass):
    // two-surface vertices within a few loop hops of the junction
    // can sit at stale original-boundary positions — measured: the
    // cross-seam chord's far endpoint at the original wall position,
    // 2.07e-2 off BOTH mega strips, at the boss top junction.
    let mut neighbors: Vec<popkern_core::topology::VertexHandle> = Vec::new();
    for (_, he) in brep.half_edges.iter() {
        if he.origin != vertex {
            continue;
        }
        let mut cur = he.next;
        for _ in 0..3 {
            let Some(d) = brep.half_edges.get(cur) else {
                break;
            };
            if d.origin != vertex && !neighbors.contains(&d.origin) {
                neighbors.push(d.origin);
            }
            cur = d.next;
        }
    }
    for v in neighbors {
        let mut sigs: Vec<String> = Vec::new();
        let mut nurbs = 0usize;
        let mut cap_like = false;
        for (_, he) in brep.half_edges.iter() {
            if he.origin != v {
                continue;
            }
            // Both sides of every incident edge: a cap on the TWIN
            // side evades an origin-only scan, and moving a
            // cap-boundary vertex broke the cube-uniform cap G1
            // gate to 1.0 (measured).
            for l in [
                brep.loops.get(he.loop_),
                brep.half_edges
                    .get(he.twin)
                    .and_then(|tw| brep.loops.get(tw.loop_)),
            ]
            .into_iter()
            .flatten()
            {
                let Some(sof) = brep.faces.get(l.face).and_then(|fd| fd.surface.as_ref()) else {
                    continue;
                };
                let sig = format!("{:?}", sof.geometry);
                if !sigs.contains(&sig) {
                    match &sof.geometry {
                        SurfaceGeometry::Nurbs(_) => nurbs += 1,
                        SurfaceGeometry::Sphere { .. } | SurfaceGeometry::Gregory(_) => {
                            cap_like = true;
                        }
                        _ => {}
                    }
                    sigs.push(sig);
                }
            }
        }
        // 2 = stale boundary vertex; 3 = the junction TRIPLE with
        // the exact signature 2 strips + 1 analytic support
        // (measured: a 2-surface settle put the shared rail
        // endpoint on ONE strip's intersection, 1.9e-2 off the
        // other; the unconstrained 3-surface settle regressed the
        // cube-uniform cap G1 gate). Cap-adjacent vertices and
        // 4+-surface hole corners never settle here.
        let triple = sigs.len() == 3 && nurbs == 2 && !cap_like;
        if !(sigs.len() == 2 && !cap_like || triple) {
            continue;
        }
        if let Some(n) =
            super::corner_chord_settle::settle_vertex_on_incident_surfaces(brep, v, tolerance)
            && n > 0
        {
            absorbed += n;
            if trace {
                eprintln!("[junction] {vertex:?}: stale neighbor {v:?} settled");
            }
        }
    }
    if std::env::var_os("PK177_JOIN").is_some()
        && join_smooth_junction_strips(
            brep,
            SmoothJunctionJoin {
                surface_a: &edge_geometries[a].surface,
                surface_b: &edge_geometries[b].surface,
                v_pos,
                chain_direction: tb,
                radius,
                trace,
            },
        )
    {
        absorbed += 1;
    }
    absorbed
}
