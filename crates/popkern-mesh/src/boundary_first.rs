//! Boundary-first face tessellation (PK-095 closure).
//!
//! Face meshes sampled independently leave T-junction cracks along
//! every shared edge (the cap-rim "flap ring"). This module makes
//! adjacent faces watertight BY CONSTRUCTION: every edge is sampled
//! ONCE into a canonical polyline (adaptive to the 3D curve's
//! flatness), and each face's boundary vertices ARE those canonical
//! points — both sides of an edge use bit-identical positions. The
//! face interior is filled with grid points and triangulated against
//! the fixed boundary by a Bowyer–Watson Delaunay in UV, filtered to
//! the trim polygon.
//!
//! Faces without complete pcurve loops (or without edge curves) fall
//! back to the legacy grid+clip path in `tessellate_face`.
//!
//! Observability: set `PK_BF_TRACE=1` to log, per face falling back
//! to the legacy tessellator, WHY it fell back (missing pcurve, no
//! canonical samples, off-surface residual, uncovered ring segment,
//! periodic full-wrap ring).

use std::collections::HashMap;

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle};

use crate::TessellationOptions;
use crate::mesh::{FaceMesh, MeshVertex};

/// Canonical per-edge samples: parameter fractions (0..=1 of the
/// edge's curve window) and their 3D points. Both faces adjacent to
/// the edge consume the same set.
pub(crate) struct EdgeSamples {
    /// Curve-window fractions, ascending, including 0 and 1.
    pub fractions: Vec<f64>,
    /// 3D points at those fractions.
    pub points: Vec<Point3>,
}

/// Sample every edge with a finite curve window once, adaptively by
/// 3D chord height.
pub(crate) fn sample_all_edges(brep: &BRep, chord_height: f64) -> HashMap<EdgeHandle, EdgeSamples> {
    let mut map = HashMap::new();
    for (eh, edge) in brep.edges.iter() {
        let Some(curve) = edge.curve.as_ref() else {
            continue;
        };
        if !(curve.t_start.is_finite() && curve.t_end.is_finite()) {
            continue;
        }
        let eval: &dyn CurveEval = &curve.geometry;
        let at = |f: f64| eval.point_at(curve.t_start + (curve.t_end - curve.t_start) * f);
        // Adaptive subdivision by chord height, depth-limited.
        let mut fractions = vec![0.0, 1.0];
        let mut i = 0;
        while i + 1 < fractions.len() {
            let (a, b) = (fractions[i], fractions[i + 1]);
            let (pa, pm, pb) = (at(a), at(0.5 * (a + b)), at(b));
            let chord_mid = Point3::new(
                0.5 * (pa.x + pb.x),
                0.5 * (pa.y + pb.y),
                0.5 * (pa.z + pb.z),
            );
            if pm.distance_to(&chord_mid) > chord_height && fractions.len() < 256 && (b - a) > 1e-3
            {
                fractions.insert(i + 1, 0.5 * (a + b));
            } else {
                i += 1;
            }
        }
        let mut points: Vec<Point3> = fractions.iter().map(|&f| at(f)).collect();
        // Snap the endpoint samples to the shared VERTEX positions:
        // curve-endpoint evaluations can differ from the topological
        // vertex by ~1e-8 (refit/snap noise), and two edges meeting
        // at a vertex would then produce near-miss boundary points
        // that straddle the weld rounding.
        let he = edge.half_edge;
        if let (Some(first), Some(last)) = (points.first().copied(), points.last().copied()) {
            let va = brep.vertices[brep.half_edges[he].origin].position;
            let vb = brep.vertices[brep.half_edges[brep.half_edges[he].twin].origin].position;
            let (pa, pb) = if first.distance_to(&va) + last.distance_to(&vb)
                <= first.distance_to(&vb) + last.distance_to(&va)
            {
                (va, vb)
            } else {
                (vb, va)
            };
            let snap_gate = 1.0e-6;
            if first.distance_to(&pa) < snap_gate {
                points[0] = pa;
            }
            let n = points.len() - 1;
            if points[n].distance_to(&pb) < snap_gate {
                points[n] = pb;
            }
        }
        map.insert(eh, EdgeSamples { fractions, points });
    }
    map
}

/// One boundary ring: ordered UV polygon with matched 3D positions.
struct BoundaryRing {
    uv: Vec<(f64, f64)>,
    pos: Vec<Point3>,
}

/// Tessellate a face against canonical edge samples. Returns `None`
/// when the face lacks the data (pcurves / curves) — caller falls
/// back to the legacy path.
pub(crate) fn tessellate_face_boundary_first(
    brep: &BRep,
    face_handle: FaceHandle,
    edge_samples: &HashMap<EdgeHandle, EdgeSamples>,
    options: &TessellationOptions,
) -> Option<FaceMesh> {
    let face = brep.faces.get(face_handle)?;
    let surface_on_face = face.surface.as_ref()?;
    let surface = &surface_on_face.geometry;
    let orientation = surface_on_face.orientation;

    // Build boundary rings from canonical samples.
    let mut rings: Vec<BoundaryRing> = Vec::new();
    let mut loops = vec![face.outer_loop];
    loops.extend(face.inner_loops.iter().copied());
    for lh in loops {
        let mut ring = build_ring(brep, lh, edge_samples)?;
        // Periodic surfaces (analytic cylinders: u = angle) can have
        // pcurves authored in different period branches; a ring that
        // crosses the seam jumps ~2π in the UV chart and Delaunay
        // cannot cover the jump segment. Unwrap u to the branch
        // nearest each previous sample (positions are canonical 3D
        // points and unaffected; cos/sin accept any branch).
        {
            use popkern_core::geometry::surface::SurfaceEval as _;
            let (period_u, period_v) = surface.parameter_periodicity();
            let mut unwrap_dir = |period: f64, pick: fn(&mut (f64, f64)) -> &mut f64| -> bool {
                for i in 1..ring.uv.len() {
                    let prev = *pick(&mut ring.uv[i - 1].clone());
                    let cur = pick(&mut ring.uv[i]);
                    while *cur - prev > 0.5 * period {
                        *cur -= period;
                    }
                    while *cur - prev < -0.5 * period {
                        *cur += period;
                    }
                }
                // Net winding = a genuine full wrap; boundary-first
                // cannot chart it.
                match (ring.uv.first(), ring.uv.last()) {
                    (Some(f), Some(l)) => {
                        let (mut fc, mut lc) = (*f, *l);
                        (*pick(&mut lc) - *pick(&mut fc)).abs() <= 0.5 * period
                    }
                    _ => true,
                }
            };
            if let Some(p) = period_u {
                if !unwrap_dir(p, |uv| &mut uv.0) {
                    if std::env::var_os("PK_BF_TRACE").is_some() {
                        eprintln!("    [bf] full u-wrap ring — cannot chart");
                    }
                    return None;
                }
            }
            if let Some(p) = period_v {
                if !unwrap_dir(p, |uv| &mut uv.1) {
                    if std::env::var_os("PK_BF_TRACE").is_some() {
                        eprintln!("    [bf] full v-wrap ring — cannot chart");
                    }
                    return None;
                }
            }
        }
        if ring.uv.len() >= 3 {
            rings.push(ring);
        }
    }
    if rings.is_empty() {
        return None;
    }
    // Collinear boundary runs (straight edges sampled at several
    // stations) sit exactly on one line in UV; Bowyer–Watson then
    // connects the run's endpoints directly and the interior
    // stations never join the hull — the coverage check below would
    // reject the whole face. Nudge each strictly-collinear boundary
    // point INWARD by an invisible epsilon in UV (3D positions stay
    // canonical, so shared-edge closure is unaffected).
    for ring in &mut rings {
        let n = ring.uv.len();
        if n < 3 {
            continue;
        }
        let mut scale = 0.0_f64;
        for i in 0..n {
            let a = ring.uv[i];
            let b = ring.uv[(i + 1) % n];
            scale = scale.max((a.0 - b.0).abs().max((a.1 - b.1).abs()));
        }
        let eps = scale.max(1e-9) * 1.0e-7;
        let mut area2 = 0.0;
        for i in 0..n {
            let (x1, y1) = ring.uv[i];
            let (x2, y2) = ring.uv[(i + 1) % n];
            area2 += x1 * y2 - x2 * y1;
        }
        let inward = if area2 > 0.0 { 1.0 } else { -1.0 };
        let orig: Vec<(f64, f64)> = ring.uv.clone();
        for i in 0..n {
            let p = orig[i];
            let q = orig[(i + n - 1) % n];
            let r = orig[(i + 1) % n];
            let cross = (p.0 - q.0) * (r.1 - q.1) - (p.1 - q.1) * (r.0 - q.0);
            let span = ((r.0 - q.0).powi(2) + (r.1 - q.1).powi(2))
                .sqrt()
                .max(1e-300);
            if (cross / span).abs() < eps {
                // Left normal of q→r is inward for CCW rings.
                let (dx, dy) = ((r.0 - q.0) / span, (r.1 - q.1) / span);
                ring.uv[i] = (
                    p.0 - dy * eps * inward * -1.0,
                    p.1 + dx * eps * inward * -1.0,
                );
            }
        }
    }

    // Point set: boundary points (fixed) + interior grid points.
    let mut uv_pts: Vec<(f64, f64)> = Vec::new();
    let mut pos_pts: Vec<Option<Point3>> = Vec::new(); // None → eval surface
    let mut ring_ranges: Vec<(usize, usize)> = Vec::new();
    for ring in &rings {
        let start = uv_pts.len();
        uv_pts.extend(ring.uv.iter().copied());
        pos_pts.extend(ring.pos.iter().copied().map(Some));
        ring_ranges.push((start, uv_pts.len()));
    }

    // Interior grid from the outer ring's bbox; keep points strictly
    // inside the trim region and away from the boundary (a third of
    // a cell) so Delaunay doesn't create slivers against the rings.
    let (mut u0, mut u1, mut v0, mut v1) = (
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::INFINITY,
        f64::NEG_INFINITY,
    );
    for &(u, v) in &rings[0].uv {
        u0 = u0.min(u);
        u1 = u1.max(u);
        v0 = v0.min(v);
        v1 = v1.max(v);
    }
    let n = options.min_subdivisions.max(1);
    // Refine like the legacy path would for curved surfaces: sample
    // the surface diagonal and bump density if it deviates from the
    // chord more than chord_height.
    let mut n = n;
    {
        let p00 = surface.point_at(u0, v0);
        let p11 = surface.point_at(u1, v1);
        let pm = surface.point_at(0.5 * (u0 + u1), 0.5 * (v0 + v1));
        let chord_mid = Point3::new(
            0.5 * (p00.x + p11.x),
            0.5 * (p00.y + p11.y),
            0.5 * (p00.z + p11.z),
        );
        let dev = pm.distance_to(&chord_mid);
        while (dev / (n as f64 * n as f64)) > options.chord_height && n < options.max_subdivisions {
            n *= 2;
        }
    }
    let cell = ((u1 - u0) / n as f64).max((v1 - v0) / n as f64).max(1e-12);
    let clearance = cell / 3.0;
    for i in 1..n {
        for j in 1..n {
            let u = u0 + (u1 - u0) * i as f64 / n as f64;
            let v = v0 + (v1 - v0) * j as f64 / n as f64;
            if point_in_region((u, v), &rings) && dist_to_rings((u, v), &rings) > clearance {
                uv_pts.push((u, v));
                pos_pts.push(None);
            }
        }
    }

    // Delaunay in UV.
    let tris = bowyer_watson(&uv_pts)?;

    // Filter to the trim region and enforce winding.
    let mut triangles: Vec<[usize; 3]> = Vec::new();
    for t in tris {
        let (a, b, c) = (uv_pts[t[0]], uv_pts[t[1]], uv_pts[t[2]]);
        let cx = (a.0 + b.0 + c.0) / 3.0;
        let cy = (a.1 + b.1 + c.1) / 3.0;
        if !point_in_region((cx, cy), &rings) {
            continue;
        }
        let signed = (b.0 - a.0) * (c.1 - a.1) - (c.0 - a.0) * (b.1 - a.1);
        if signed.abs() < 1e-16 {
            continue;
        }
        if signed > 0.0 {
            triangles.push(t);
        } else {
            triangles.push([t[0], t[2], t[1]]);
        }
    }
    // Boundary-edge coverage check: every ring segment must appear as
    // a triangle edge, else the Delaunay cut a concave corner — fall
    // back to the legacy path rather than ship a gap.
    let mut edge_set = std::collections::HashSet::new();
    for t in &triangles {
        for (x, y) in [(t[0], t[1]), (t[1], t[2]), (t[2], t[0])] {
            edge_set.insert((x.min(y), x.max(y)));
        }
    }
    for &(start, end) in &ring_ranges {
        let m = end - start;
        for k in 0..m {
            let a = start + k;
            let b = start + (k + 1) % m;
            if !edge_set.contains(&(a.min(b), a.max(b))) {
                if std::env::var_os("PK_BF_TRACE").is_some() {
                    eprintln!(
                        "    [bf] ring segment {a}-{b} uncovered at uv ({:.4},{:.4})-({:.4},{:.4})",
                        uv_pts[a].0, uv_pts[a].1, uv_pts[b].0, uv_pts[b].1
                    );
                }
                return None;
            }
        }
    }

    // Materialize vertices.
    let mut vertices = Vec::with_capacity(uv_pts.len());
    for (i, &(u, v)) in uv_pts.iter().enumerate() {
        let position = pos_pts[i].unwrap_or_else(|| surface.point_at(u, v));
        let mut normal = surface.normal_at(u, v);
        let len = normal.length();
        if len > 1e-300 {
            normal = Vec3::new(normal.x / len, normal.y / len, normal.z / len);
        }
        if !orientation {
            normal = -normal;
        }
        vertices.push(MeshVertex {
            position,
            normal,
            uv: (u, v),
        });
    }
    if !orientation {
        for t in &mut triangles {
            t.swap(1, 2);
        }
    }
    Some(FaceMesh {
        face: face_handle,
        vertices,
        triangles,
    })
}

/// Walk one loop, concatenating canonical edge samples in loop order.
/// UV comes from each half-edge's pcurve at the matching fraction
/// (reversed when the half-edge runs against the edge curve); 3D is
/// the shared canonical point.
fn build_ring(
    brep: &BRep,
    loop_handle: popkern_core::topology::LoopHandle,
    edge_samples: &HashMap<EdgeHandle, EdgeSamples>,
) -> Option<BoundaryRing> {
    let start = brep.loops.get(loop_handle)?.half_edge;
    let mut uv = Vec::new();
    let mut pos = Vec::new();
    let mut he = start;
    for _ in 0..10_000 {
        let data = &brep.half_edges[he];
        let eh = data.edge;
        let Some(pc) = data.pcurve.as_ref() else {
            if std::env::var_os("PK_BF_TRACE").is_some() {
                eprintln!("    [bf] edge {eh:?}: missing pcurve");
            }
            return None;
        };
        let Some(samples) = edge_samples.get(&eh) else {
            if std::env::var_os("PK_BF_TRACE").is_some() {
                eprintln!("    [bf] edge {eh:?}: no canonical samples");
            }
            return None;
        };
        let origin = brep.vertices[data.origin].position;

        // Direction: does the canonical polyline start at this
        // half-edge's origin?
        let forward = samples.points.first()?.distance_to(&origin)
            <= samples.points.last()?.distance_to(&origin);
        let m = samples.points.len();
        // Emit all but the last point (the next half-edge provides
        // its own start, which is the shared vertex).
        for k in 0..(m - 1) {
            let idx = if forward { k } else { m - 1 - k };
            let f = samples.fractions[idx];
            // Pcurve param at the same fraction of ITS window,
            // mirrored when the pcurve runs against the walk
            // direction. The pcurve's own direction is detected from
            // its endpoint images.
            let p_start = pc.point_at(pc.t_start);
            let p_end = pc.point_at(pc.t_end);
            let sof = brep.faces[brep.loops[data.loop_].face].surface.as_ref()?;
            let s_start = sof.geometry.point_at(p_start.x, p_start.y);
            let s_end = sof.geometry.point_at(p_end.x, p_end.y);
            let pcurve_forward_along_walk =
                s_start.distance_to(&origin) <= s_end.distance_to(&origin);
            // Fraction along the WALK is k/(m-1); convert to the
            // curve fraction (f) orientation… the pcurve maps the
            // same geometric edge, so use the walk fraction with the
            // pcurve's own direction.
            let walk_f = k as f64 / (m - 1) as f64;
            let _ = f;
            let tp = if pcurve_forward_along_walk {
                pc.t_start + (pc.t_end - pc.t_start) * walk_f
            } else {
                pc.t_end - (pc.t_end - pc.t_start) * walk_f
            };
            // The pcurve is only a SEED: refine each UV by point
            // inversion of the canonical 3D sample so imprecise
            // pcurves (straight chords for bulging curves) still
            // yield exact boundary UVs. If the refined image stays
            // far from the canonical point, the curve genuinely does
            // not lie on this surface (sphere seam pole chord) —
            // fall back to the legacy path.
            let uvp = pc.point_at(tp);
            let sof = brep.faces[brep.loops[data.loop_].face].surface.as_ref()?;
            let target = samples.points[idx];
            let (mut uu, mut vv) = (uvp.x, uvp.y);
            // Analytic cylinders: invert analytically (atan2 in the
            // surface's own basis) — the generic Newton projection
            // wanders across period branches and reversed-branch
            // pcurve seeds walked rings the long way around the seam
            // (span 2π). Canonical branch here; the ring-level
            // unwrap restores continuity.
            if let popkern_core::geometry::SurfaceGeometry::Cylinder { origin, axis, .. } =
                &sof.geometry
            {
                let (perp1, perp2) =
                    popkern_core::geometry::analytic_surfaces::orthonormal_basis(*axis);
                let d = target - *origin;
                (uu, vv) = (d.dot(&perp2).atan2(d.dot(&perp1)), d.dot(axis));
            } else if sof.geometry.point_at(uu, vv).distance_to(&target) > 1.0e-9 {
                if let Some((ru, rv)) =
                    popkern_core::geometry::uv_projection::project_point_to_uv_near(
                        &sof.geometry,
                        &target,
                        uu,
                        vv,
                    )
                {
                    (uu, vv) = (ru, rv);
                }
            }
            let residual = sof.geometry.point_at(uu, vv).distance_to(&target);
            let scale = samples
                .points
                .first()?
                .distance_to(samples.points.last()?)
                .max(1e-9);
            if residual > 1.0e-3 * scale.max(1.0) {
                if std::env::var_os("PK_BF_TRACE").is_some() {
                    eprintln!("    [bf] edge {eh:?}: off-surface residual {residual:.3e}");
                }
                return None;
            }
            uv.push((uu, vv));
            pos.push(target);
        }
        he = data.next;
        if he == start {
            break;
        }
    }
    Some(BoundaryRing { uv, pos })
}

/// Even-odd point-in-region across all rings (outer + holes).
fn point_in_region(p: (f64, f64), rings: &[BoundaryRing]) -> bool {
    let mut inside = false;
    for ring in rings {
        let n = ring.uv.len();
        for i in 0..n {
            let (x1, y1) = ring.uv[i];
            let (x2, y2) = ring.uv[(i + 1) % n];
            if (y1 > p.1) != (y2 > p.1) {
                let x_int = x1 + (p.1 - y1) / (y2 - y1) * (x2 - x1);
                if p.0 < x_int {
                    inside = !inside;
                }
            }
        }
    }
    inside
}

fn dist_to_rings(p: (f64, f64), rings: &[BoundaryRing]) -> f64 {
    let mut best = f64::INFINITY;
    for ring in rings {
        let n = ring.uv.len();
        for i in 0..n {
            let a = ring.uv[i];
            let b = ring.uv[(i + 1) % n];
            let (dx, dy) = (b.0 - a.0, b.1 - a.1);
            let l2 = dx * dx + dy * dy;
            let t = if l2 < 1e-300 {
                0.0
            } else {
                (((p.0 - a.0) * dx + (p.1 - a.1) * dy) / l2).clamp(0.0, 1.0)
            };
            let (qx, qy) = (a.0 + dx * t, a.1 + dy * t);
            best = best.min(((p.0 - qx).powi(2) + (p.1 - qy).powi(2)).sqrt());
        }
    }
    best
}

/// Bowyer–Watson Delaunay triangulation over a point set. Returns
/// `None` on degenerate input (all collinear etc.).
fn bowyer_watson(points: &[(f64, f64)]) -> Option<Vec<[usize; 3]>> {
    if points.len() < 3 {
        return None;
    }
    // Super-triangle enclosing everything.
    let (mut x0, mut x1, mut y0, mut y1) = (
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::INFINITY,
        f64::NEG_INFINITY,
    );
    for &(x, y) in points {
        x0 = x0.min(x);
        x1 = x1.max(x);
        y0 = y0.min(y);
        y1 = y1.max(y);
    }
    let d = ((x1 - x0).max(y1 - y0)).max(1e-9) * 16.0;
    let cx = 0.5 * (x0 + x1);
    let cy = 0.5 * (y0 + y1);
    let mut pts: Vec<(f64, f64)> = points.to_vec();
    let s0 = pts.len();
    pts.push((cx - d, cy - d * 0.6));
    pts.push((cx + d, cy - d * 0.6));
    pts.push((cx, cy + d));

    let mut tris: Vec<[usize; 3]> = vec![[s0, s0 + 1, s0 + 2]];
    for (pi, &p) in points.iter().enumerate() {
        // Find triangles whose circumcircle contains p.
        let mut bad: Vec<usize> = Vec::new();
        for (ti, t) in tris.iter().enumerate() {
            if in_circumcircle(p, pts[t[0]], pts[t[1]], pts[t[2]]) {
                bad.push(ti);
            }
        }
        // Boundary of the cavity: edges appearing exactly once.
        let mut edge_count: HashMap<(usize, usize), usize> = HashMap::new();
        let mut edge_dir: HashMap<(usize, usize), (usize, usize)> = HashMap::new();
        for &ti in &bad {
            let t = tris[ti];
            for (a, b) in [(t[0], t[1]), (t[1], t[2]), (t[2], t[0])] {
                let key = (a.min(b), a.max(b));
                *edge_count.entry(key).or_insert(0) += 1;
                edge_dir.insert(key, (a, b));
            }
        }
        // Remove bad triangles (descending index).
        bad.sort_unstable_by(|a, b| b.cmp(a));
        for ti in bad {
            tris.swap_remove(ti);
        }
        for (key, count) in edge_count {
            if count == 1 {
                let (a, b) = edge_dir[&key];
                tris.push([a, b, pi]);
            }
        }
    }
    tris.retain(|t| t.iter().all(|&i| i < s0));
    if tris.is_empty() { None } else { Some(tris) }
}

fn in_circumcircle(p: (f64, f64), a: (f64, f64), b: (f64, f64), c: (f64, f64)) -> bool {
    // Standard incircle determinant with orientation normalization.
    let orient = (b.0 - a.0) * (c.1 - a.1) - (b.1 - a.1) * (c.0 - a.0);
    if orient.abs() < 1e-300 {
        return false;
    }
    let (ax, ay) = (a.0 - p.0, a.1 - p.1);
    let (bx, by) = (b.0 - p.0, b.1 - p.1);
    let (cx, cy) = (c.0 - p.0, c.1 - p.1);
    let det = (ax * ax + ay * ay) * (bx * cy - by * cx) - (bx * bx + by * by) * (ax * cy - ay * cx)
        + (cx * cx + cy * cy) * (ax * by - ay * bx);
    if orient > 0.0 { det > 0.0 } else { det < 0.0 }
}
