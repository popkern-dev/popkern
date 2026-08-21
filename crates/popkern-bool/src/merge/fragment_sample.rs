//! Fragment sample point computation for classification.
//!
//! Provides functions to compute representative 3D points for face fragments,
//! used by the point-in-solid classifier to determine fragment membership.

use popkern_core::{
    geometry::{Point3, surface::SurfaceEval},
    topology::BRep,
};

use super::graph::{EdgeSource, FaceRef, IntersectionGraph, Operand};

/// Pick a representative sample point for a fragment.
///
/// Prefers the midpoint of a non-shared (original or split) edge, which is
/// guaranteed to lie well inside the fragment and away from the intersection
/// curve.  Falls back to the loop centroid if every edge is shared.
pub(crate) fn sample_point_for_fragment(graph: &IntersectionGraph, loop_hes: &[usize]) -> Point3 {
    let midpoint_of = |he_idx: usize| -> Point3 {
        let he = &graph.half_edges[he_idx];
        let twin_idx = he.twin;
        // Curve-true midpoint: a curved edge's CHORD midpoint is off
        // the curve, and it can land exactly ON the opposite operand's
        // boundary (PK-295, measured: the rim arc through θ=π between
        // wall crossings at (−0.9,±0.436,1.8) chord-mids to
        // (−0.9,0,1.8) — precisely on the cutter's x=−0.9 wall — and
        // the knife-edge in/out verdict discarded three of the four
        // torus bulge fragments). Graph splits rewindow the stored
        // curve per piece, so the mid-parameter evaluation is honest.
        let curve = he
            .curve
            .as_ref()
            .or(graph.half_edges[twin_idx].curve.as_ref());
        if let Some(c) = curve {
            use popkern_core::geometry::curve::CurveEval;
            return (&c.geometry as &dyn CurveEval).point_at(0.5 * (c.t_start + c.t_end));
        }
        let start = graph.nodes[he.origin].position;
        let end = graph.nodes[graph.half_edges[twin_idx].origin].position;
        Point3::new(
            (start.x + end.x) * 0.5,
            (start.y + end.y) * 0.5,
            (start.z + end.z) * 0.5,
        )
    };
    let picked = |he_idx: usize, p: Point3| -> Point3 {
        if std::env::var_os("PK_SAMPLE_TRACE").is_some() {
            let c = graph.half_edges[he_idx]
                .curve
                .as_ref()
                .or(graph.half_edges[graph.half_edges[he_idx].twin]
                    .curve
                    .as_ref())
                .map(|c| (c.t_start, c.t_end));
            eprintln!(
                "  [sample] he {he_idx} src {:?} w {c:?} -> ({:.3},{:.3},{:.3})",
                graph.half_edges[he_idx].source, p.x, p.y, p.z
            );
        }
        p
    };
    // Try split edges first.
    for &he_idx in loop_hes {
        if matches!(
            graph.half_edges[he_idx].source,
            EdgeSource::SplitA { .. } | EdgeSource::SplitB { .. }
        ) {
            return picked(he_idx, midpoint_of(he_idx));
        }
    }
    // Then try original edges.
    for &he_idx in loop_hes {
        if matches!(
            graph.half_edges[he_idx].source,
            EdgeSource::OriginalA(_) | EdgeSource::OriginalB(_)
        ) {
            return picked(he_idx, midpoint_of(he_idx));
        }
    }
    if std::env::var_os("PK_SAMPLE_TRACE").is_some() {
        eprintln!("  [sample] centroid fallback over {} hes", loop_hes.len());
    }
    compute_centroid(graph, loop_hes)
}

/// Compute the centroid of a loop's vertex positions.
pub(crate) fn compute_centroid(graph: &IntersectionGraph, loop_hes: &[usize]) -> Point3 {
    if loop_hes.is_empty() {
        return Point3::new(0.0, 0.0, 0.0);
    }

    let mut sum_x = 0.0;
    let mut sum_y = 0.0;
    let mut sum_z = 0.0;

    for &he_idx in loop_hes {
        let origin = graph.half_edges[he_idx].origin;
        let pos = graph.nodes[origin].position;
        sum_x += pos.x;
        sum_y += pos.y;
        sum_z += pos.z;
    }

    let n = loop_hes.len() as f64;
    Point3::new(sum_x / n, sum_y / n, sum_z / n)
}

/// Sample a point in the fragment's UV-space interior, then evaluate
/// the face surface to get a 3D point.
///
/// Computes the UV centroid from pcurve vertices (area-weighted polygon
/// centroid in UV space). For planar faces, this reliably maps to the
/// face interior. Returns `None` if pcurves or surface geometry are
/// unavailable.
pub(crate) fn sample_point_uv(
    graph: &IntersectionGraph,
    loop_hes: &[usize],
    face: &FaceRef,
    brep_a: &BRep,
    brep_b: &BRep,
) -> Option<Point3> {
    let brep = match face.operand {
        Operand::A => brep_a,
        Operand::B => brep_b,
    };
    let face_data = brep.faces.get(face.face)?;
    let surface = face_data.surface.as_ref()?;

    let mut uv_verts: Vec<(f64, f64)> = Vec::new();
    for &he_idx in loop_hes {
        let he = &graph.half_edges[he_idx];
        let (pcurve, use_t_end) = if he.face.as_ref() == Some(face) {
            (he.pcurve.as_ref(), false)
        } else {
            let twin = &graph.half_edges[he.twin];
            if twin.face.as_ref() == Some(face) {
                (twin.pcurve.as_ref(), true)
            } else {
                (he.pcurve.as_ref(), false)
            }
        };
        if let Some(pc) = pcurve {
            // Endpoint parameters live at the pcurve's intrinsic
            // range; uniform across all variants (Linear/Nurbs).
            let (t_start, t_end) = pc.parameter_range();
            let t = if use_t_end { t_end } else { t_start };
            let uv_pt = pc.point_at(t);
            uv_verts.push((uv_pt.x, uv_pt.y));
        }
    }

    if uv_verts.len() < 3 {
        let n = loop_hes.len() as f64;
        if n < 3.0 {
            return None;
        }
        let mut cx = 0.0;
        let mut cy = 0.0;
        let mut cz = 0.0;
        for &he_idx in loop_hes {
            let pos = graph.nodes[graph.half_edges[he_idx].origin].position;
            cx += pos.x;
            cy += pos.y;
            cz += pos.z;
        }
        return Some(Point3::new(cx / n, cy / n, cz / n));
    }

    let n = uv_verts.len();
    let mut area_2 = 0.0_f64;
    let mut cu = 0.0_f64;
    let mut cv = 0.0_f64;
    for i in 0..n {
        let (u0, v0) = uv_verts[i];
        let (u1, v1) = uv_verts[(i + 1) % n];
        let cross = u0 * v1 - u1 * v0;
        area_2 += cross;
        cu += (u0 + u1) * cross;
        cv += (v0 + v1) * cross;
    }

    if area_2.abs() < 1e-20 {
        let u = uv_verts.iter().map(|(u, _)| u).sum::<f64>() / n as f64;
        let v = uv_verts.iter().map(|(_, v)| v).sum::<f64>() / n as f64;
        return Some(surface.geometry.point_at(u, v));
    }

    let u = cu / (3.0 * area_2);
    let v = cv / (3.0 * area_2);
    // The polygon centroid of a NON-CONVEX loop can lie outside the
    // loop (PK-295, measured: a cutter wall's notched frame fragment —
    // a square outline with the pocket strip notched out — centroids
    // into the notch, maps into the pocket interior, reads Inside the
    // other operand, and the keep-biased retry ships a frame that is
    // almost entirely outside it). A sample that does not lie in the
    // fragment's own UV polygon classifies a DIFFERENT region — only
    // in-polygon samples are honest. Fall back to inward-nudged edge
    // midpoints when the centroid misses.
    if point_in_uv_polygon(u, v, &uv_verts) {
        return Some(surface.geometry.point_at(u, v));
    }
    let orient = area_2.signum();
    for i in 0..n {
        let (u0, v0) = uv_verts[i];
        let (u1, v1) = uv_verts[(i + 1) % n];
        let (mu, mv) = (0.5 * (u0 + u1), 0.5 * (v0 + v1));
        let (du, dv) = (u1 - u0, v1 - v0);
        let len = (du * du + dv * dv).sqrt();
        if len < 1e-12 {
            continue;
        }
        // Inward normal for a CCW (positive-area) polygon is the left
        // normal; flipped for CW.
        let (nu, nv) = (-dv / len * orient, du / len * orient);
        for scale in [1e-3, 1e-2, 1e-1] {
            let step = scale * len;
            let (su, sv) = (mu + nu * step, mv + nv * step);
            if point_in_uv_polygon(su, sv, &uv_verts) {
                return Some(surface.geometry.point_at(su, sv));
            }
        }
    }
    None
}

/// Even-odd containment of a UV point in a polygon.
fn point_in_uv_polygon(u: f64, v: f64, poly: &[(f64, f64)]) -> bool {
    let n = poly.len();
    let mut inside = false;
    for i in 0..n {
        let (u0, v0) = poly[i];
        let (u1, v1) = poly[(i + 1) % n];
        if (v0 > v) != (v1 > v) {
            let cross_u = u0 + (v - v0) / (v1 - v0) * (u1 - u0);
            if u < cross_u {
                inside = !inside;
            }
        }
    }
    inside
}
