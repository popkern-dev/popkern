//! Face lineage for boolean results (PK-633).
//!
//! The stable names `entity_naming` assigns are kept unchanged (a
//! `Generated` name stays `Generated`: saved references to it keep
//! resolving). Beside each output face this records WHICH input faces
//! it is a piece of, so a label or spelled name bound to an input face
//! resolves to its pieces after the boolean
//! (`popkern_core::naming_resolve::resolve_lineage`).
//!
//! The relation is decided geometrically and independently of how the
//! pipeline produced the face: an interior point of the output face
//! (inside its own trim, `BRep::face_uv_trim_containment`) lies ON an
//! input face's surface (within the boolean's tolerance, normals
//! parallel) and INSIDE that input face's trim. A fragment of a face
//! always satisfies this however much the boolean cut it — the case the
//! centroid match misses (a pin's cylinder trimmed by the plate it was
//! united with). Coincident faces of both operands (coplanar caps)
//! record both. A face that is a piece of no input face (an
//! intersection-born face) records nothing.

use popkern_core::geometry::intersection::surface_subdivision::Aabb;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::uv_projection::project_point_to_uv;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};
use popkern_core::validation::FaceUvContainment;

/// A point strictly inside `face` (inside its own trim) and the
/// surface normal there, with the face's orientation applied.
pub(crate) fn interior_point(brep: &BRep, face: FaceHandle) -> Option<(Point3, Vec3)> {
    interior_point_where(brep, face, |_| true)
}

/// [`interior_point`]'s first candidate inside the trim that `keep`
/// accepts (a second point away from the first).
pub(crate) fn interior_point_where(
    brep: &BRep,
    face: FaceHandle,
    keep: impl Fn(&Point3) -> bool,
) -> Option<(Point3, Vec3)> {
    let surface = brep.faces.get(face)?.surface.as_ref()?;
    let geometry = &surface.geometry;
    let samples = brep.face_outer_loop_sample_positions(face, 4);
    let uvs: Vec<(f64, f64)> = samples
        .iter()
        .filter_map(|p| project_point_to_uv(geometry, p))
        .collect();
    if uvs.is_empty() {
        return None;
    }
    let (mut u0, mut u1, mut v0, mut v1) = (
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::INFINITY,
        f64::NEG_INFINITY,
    );
    for (u, v) in &uvs {
        u0 = u0.min(*u);
        u1 = u1.max(*u);
        v0 = v0.min(*v);
        v1 = v1.max(*v);
    }
    let n = uvs.len() as f64;
    let centre = (
        uvs.iter().map(|p| p.0).sum::<f64>() / n,
        uvs.iter().map(|p| p.1).sum::<f64>() / n,
    );
    // The loop's uv centroid first, then a grid over its uv bounds.
    let grid = (1..8).flat_map(|i| (1..8).map(move |j| (i, j)));
    let candidates = std::iter::once(centre).chain(grid.map(|(i, j)| {
        (
            u0 + (u1 - u0) * f64::from(i) / 8.0,
            v0 + (v1 - v0) * f64::from(j) / 8.0,
        )
    }));
    for (u, v) in candidates {
        if keep(&geometry.point_at(u, v))
            && brep.face_uv_trim_containment(face, geometry, u, v)
                == Some(FaceUvContainment::Interior)
        {
            // normal_at is not unit on every surface (a cylinder's is
            // radius-scaled): normalise before comparing directions.
            let raw = geometry.normal_at(u, v);
            let normal = raw.normalized().unwrap_or(raw);
            let sign = if surface.orientation { 1.0 } else { -1.0 };
            return Some((geometry.point_at(u, v), normal * sign));
        }
    }
    None
}

pub(crate) fn padded(aabb: &Aabb, pad: f64) -> Aabb {
    Aabb {
        min: Point3::new(aabb.min.x - pad, aabb.min.y - pad, aabb.min.z - pad),
        max: Point3::new(aabb.max.x + pad, aabb.max.y + pad, aabb.max.z + pad),
    }
}

pub(crate) fn contains(aabb: &Aabb, p: &Point3) -> bool {
    (aabb.min.x..=aabb.max.x).contains(&p.x)
        && (aabb.min.y..=aabb.max.y).contains(&p.y)
        && (aabb.min.z..=aabb.max.z).contains(&p.z)
}

/// Whether the surface point `p` (normal `n`) lies on `face` of `input`:
/// on its surface, normals parallel, inside its trim.
pub(crate) fn lies_on(input: &BRep, face: FaceHandle, p: &Point3, n: &Vec3, tol: f64) -> bool {
    let Some(surface) = input.faces.get(face).and_then(|f| f.surface.as_ref()) else {
        return false;
    };
    let geometry: &SurfaceGeometry = &surface.geometry;
    let Some((u, v)) = project_point_to_uv(geometry, p) else {
        return false;
    };
    if (geometry.point_at(u, v) - *p).length() > tol {
        return false;
    }
    let raw = geometry.normal_at(u, v);
    let m = raw.normalized().unwrap_or(raw);
    if n.dot(&m).abs() < 1.0 - 1e-6 {
        return false;
    }
    matches!(
        input.face_uv_trim_containment(face, geometry, u, v),
        Some(FaceUvContainment::Interior | FaceUvContainment::Boundary)
    )
}

/// Record, beside every face of `result`, the input faces of `a` and
/// the `tools` it is a piece of (their names and ancestry, as values).
pub(crate) fn record_face_lineage(
    result: &mut BRep,
    a: &BRep,
    tools: &[&BRep],
    tol: &Tolerance,
    carry: Option<&crate::entity_lineage_carry::CarryIndex>,
) {
    let sources = crate::entity_sources::FaceSources::new(a, tools, tol, carry);
    let faces: Vec<FaceHandle> = result.faces.keys().collect();
    for face in faces {
        let Some(name_id) = result.faces.get(face).and_then(|f| f.name) else {
            continue;
        };
        // A carried face keeps `a`'s name and its lineage (scheme 4).
        if carry.is_some()
            && crate::entity_lineage_carry::CarryIndex::keeps_a_name(result, a, name_id)
        {
            continue;
        }
        let lineage = crate::entity_sources::lineage_of(&sources.of(result, face));
        result.name_pool.add_lineage(name_id, lineage);
    }
}

/// Whether parameter `t` falls inside the edge range `[t0, t1]` of
/// `curve` (period-aware for closed conics).
pub(crate) fn in_range(
    curve: &popkern_core::geometry::curve::CurveGeometry,
    t: f64,
    t0: f64,
    t1: f64,
) -> bool {
    use popkern_core::geometry::curve::CurveGeometry;
    let (lo, hi) = (t0.min(t1), t0.max(t1));
    let slack = 1e-9 * (hi - lo).abs().max(1.0);
    match curve {
        CurveGeometry::Circle { .. } | CurveGeometry::Ellipse { .. } => {
            let period = std::f64::consts::TAU;
            if hi - lo >= period - slack {
                return true;
            }
            let shifted = lo + (t - lo).rem_euclid(period);
            shifted <= hi + slack || shifted >= lo + period - slack
        }
        _ => t >= lo - slack && t <= hi + slack,
    }
}

/// The parameter of `point`'s foot on `curve`, as the lineage test reads
/// it (distance, range and tangent are checked by the caller). Lines,
/// straight two-point NURBS and circles/arcs are projected exactly —
/// the multi-start Newton search ran 5 starts per box-accepted input
/// edge for every result edge, and a slotted disc's rim circles accept
/// every edge (PK-803). Newton on a line lands on the same unclamped
/// projection in one step; a circle's angle is measured in the curve's
/// own frame (its points at 0 and pi/2), and the range test unwraps it.
pub(crate) fn foot(
    curve: &popkern_core::geometry::curve::CurveGeometry,
    point: &Point3,
) -> Option<f64> {
    use popkern_core::geometry::PortableMath;
    use popkern_core::geometry::curve::{CurveEval, CurveGeometry, closest_point_on_curve};
    match curve {
        CurveGeometry::Line { .. } => curve.closest_point(point),
        CurveGeometry::Circle { center, radius, .. }
        | CurveGeometry::Arc { center, radius, .. }
            if *radius > 0.0 =>
        {
            let u = curve.point_at(0.0) - *center;
            let v = curve.point_at(std::f64::consts::FRAC_PI_2) - *center;
            let w = *point - *center;
            let (x, y) = (w.dot(&u), w.dot(&v));
            let floor = f64::EPSILON * radius * radius;
            if x * x + y * y <= floor * floor {
                return closest_point_on_curve(curve, point);
            }
            Some(y.patan2(x))
        }
        CurveGeometry::Nurbs(n)
            if n.degree == 1
                && n.control_points.len() == 2
                && n.knots.len() == 4
                && n.weights.iter().all(|w| *w == n.weights[0]) =>
        {
            let (p0, p1) = (n.control_points[0], n.control_points[1]);
            let d = p1 - p0;
            let len2 = d.dot(&d);
            if len2 <= 0.0 {
                return closest_point_on_curve(curve, point);
            }
            let s = (*point - p0).dot(&d) / len2;
            Some(n.knots[1] + s * (n.knots[2] - n.knots[1]))
        }
        _ => closest_point_on_curve(curve, point),
    }
}

/// Record, beside every edge of `result`, the input edges of `a` and
/// the `tools` it is a piece of (their names and ancestry, as values).
pub(crate) fn record_edge_lineage(
    result: &mut BRep,
    a: &BRep,
    tools: &[&BRep],
    tol: &Tolerance,
    carry: Option<&crate::entity_lineage_carry::CarryIndex>,
) {
    let sources = crate::entity_sources::EdgeSources::new(a, tools, tol);
    let edges: Vec<_> = result.edges.keys().collect();
    for edge in edges {
        let Some(name_id) = result.edges.get(edge).and_then(|e| e.name) else {
            continue;
        };
        // A carried edge keeps `a`'s name and its lineage (scheme 4).
        if carry.is_some()
            && crate::entity_lineage_carry::CarryIndex::keeps_a_name(result, a, name_id)
        {
            continue;
        }
        let lineage = crate::entity_sources::lineage_of(&sources.of(result, edge));
        result.name_pool.add_lineage(name_id, lineage);
    }
}
