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
use popkern_core::naming::EntityName;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};
use popkern_core::validation::FaceUvContainment;

/// A point strictly inside `face` (inside its own trim) and the
/// surface normal there, with the face's orientation applied.
fn interior_point(brep: &BRep, face: FaceHandle) -> Option<(Point3, Vec3)> {
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
        if brep.face_uv_trim_containment(face, geometry, u, v) == Some(FaceUvContainment::Interior)
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

fn padded(aabb: &Aabb, pad: f64) -> Aabb {
    Aabb {
        min: Point3::new(aabb.min.x - pad, aabb.min.y - pad, aabb.min.z - pad),
        max: Point3::new(aabb.max.x + pad, aabb.max.y + pad, aabb.max.z + pad),
    }
}

fn contains(aabb: &Aabb, p: &Point3) -> bool {
    (aabb.min.x..=aabb.max.x).contains(&p.x)
        && (aabb.min.y..=aabb.max.y).contains(&p.y)
        && (aabb.min.z..=aabb.max.z).contains(&p.z)
}

/// Whether the surface point `p` (normal `n`) lies on `face` of `input`:
/// on its surface, normals parallel, inside its trim.
fn lies_on(input: &BRep, face: FaceHandle, p: &Point3, n: &Vec3, tol: f64) -> bool {
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
/// `b` it is a piece of (their names and ancestry, as values).
pub(crate) fn record_face_lineage(result: &mut BRep, a: &BRep, b: &BRep, tol: &Tolerance) {
    // The names' own match tolerance: the boolean's, widened for the
    // projection round-off on curved inputs.
    let on = 10.0 * tol.linear();
    let inputs: Vec<(&BRep, Vec<(FaceHandle, Aabb)>)> = [a, b]
        .into_iter()
        .map(|brep| {
            (
                brep,
                crate::intersect::compute_face_aabbs(brep)
                    .into_iter()
                    .map(|(f, bb)| (f, padded(&bb, on)))
                    .collect(),
            )
        })
        .collect();
    let faces: Vec<FaceHandle> = result.faces.keys().collect();
    for face in faces {
        let Some(name_id) = result.faces.get(face).and_then(|f| f.name) else {
            continue;
        };
        let Some((point, normal)) = interior_point(result, face) else {
            continue;
        };

        let mut lineage: Vec<EntityName> = Vec::new();
        for (input, boxes) in &inputs {
            for (candidate, aabb) in boxes {
                if !contains(aabb, &point) || !lies_on(input, *candidate, &point, &normal, on) {
                    continue;
                }
                let Some(id) = input.faces.get(*candidate).and_then(|f| f.name) else {
                    continue;
                };
                if let Some(value) = input.name_pool.get(id) {
                    lineage.push(value.clone());
                }
                lineage.extend(input.name_pool.ancestry(id));
            }
        }
        result.name_pool.add_lineage(name_id, lineage);
    }
}

/// Whether parameter `t` falls inside the edge range `[t0, t1]` of
/// `curve` (period-aware for closed conics).
fn in_range(
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

/// Record, beside every edge of `result`, the input edges of `a` and
/// `b` it is a piece of: its curve midpoint lies on the input edge's
/// curve, inside its range, tangents parallel.
pub(crate) fn record_edge_lineage(result: &mut BRep, a: &BRep, b: &BRep, tol: &Tolerance) {
    use popkern_core::geometry::curve::{CurveEval, closest_point_on_curve};
    let on = 10.0 * tol.linear();
    let midpoint = |brep: &BRep, edge: popkern_core::topology::EdgeHandle| {
        let curve = brep.edges.get(edge)?.curve.as_ref()?;
        let t = 0.5 * (curve.t_start + curve.t_end);
        let tangent = curve.geometry.tangent_at(t);
        Some((
            curve.geometry.point_at(t),
            tangent.normalized().unwrap_or(tangent),
        ))
    };
    // Each input edge's curve bounds over its range, padded: the cheap
    // reject before a projection.
    let inputs: Vec<(&BRep, Vec<(&popkern_core::topology::Edge, Aabb)>)> = [a, b]
        .into_iter()
        .map(|brep| {
            let boxes = brep
                .edges
                .values()
                .filter_map(|e| {
                    let c = e.curve.as_ref()?;
                    let bb = c
                        .geometry
                        .curve_bounding_box((c.t_start.min(c.t_end), c.t_start.max(c.t_end)));
                    Some((e, padded(&bb, on)))
                })
                .collect();
            (brep, boxes)
        })
        .collect();
    let edges: Vec<_> = result.edges.keys().collect();
    for edge in edges {
        let Some(name_id) = result.edges.get(edge).and_then(|e| e.name) else {
            continue;
        };
        let Some((point, tangent)) = midpoint(result, edge) else {
            continue;
        };
        let mut lineage: Vec<EntityName> = Vec::new();
        for (input, boxes) in &inputs {
            for (data, aabb) in boxes {
                if !contains(aabb, &point) {
                    continue;
                }
                let Some(curve) = data.curve.as_ref() else {
                    continue;
                };
                let Some(id) = data.name else { continue };
                let Some(t) = closest_point_on_curve(&curve.geometry, &point) else {
                    continue;
                };
                if (curve.geometry.point_at(t) - point).length() > on
                    || !in_range(&curve.geometry, t, curve.t_start, curve.t_end)
                {
                    continue;
                }
                let along = curve.geometry.tangent_at(t);
                let along = along.normalized().unwrap_or(along);
                if along.dot(&tangent).abs() < 1.0 - 1e-6 {
                    continue;
                }
                if let Some(value) = input.name_pool.get(id) {
                    lineage.push(value.clone());
                }
                lineage.extend(input.name_pool.ancestry(id));
            }
        }
        result.name_pool.add_lineage(name_id, lineage);
    }
}
