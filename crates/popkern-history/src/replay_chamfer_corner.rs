//! A chamfer's corners closed flat (PK-1217).
//!
//! Where three or more bevels meet at a corner, the blend runs their
//! planes on until they meet: a point at a cube's corner, a ridge where
//! four edges meet (the convex body cut by the bevel half-spaces alone).
//! The mainstream kernels close such a corner with one flat face through
//! the points where the bevels' setback lines meet on each original
//! face: a cube of side 10 chamfered 2 on every edge has the equilateral
//! triangles through its inset squares' corners, x + y + z ≤ 26 at
//! (10, 10, 10), and measures 2408/3 where the meeting planes leave 808.
//!
//! The corner tip is found from the chamfered body: the vertices all of
//! whose faces are bevel faces (planes no input face carries), joined
//! along their edges into clusters (the apex, or the ridge), and the
//! vertices one edge away from a cluster, the corner points on the
//! original faces. Those must be coplanar within the linear tolerance.
//! The tip is then cut away by a prism standing on that plane over the
//! corner polygon, widened half again about its centroid and reaching
//! well past the cluster: the corner is convex (its faces are bevels of
//! a convex corner), so the body beyond the plane there is the tip
//! alone, and the prism takes exactly it. A corner whose points are not
//! coplanar, whose bevels are not planar, or whose cut declines keeps
//! its apex, and the caller records why.

use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::history::TrackingContext;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle, VertexHandle};

/// Close every bevel corner of `body` (chamfered from `before`) flat;
/// returns why each corner left at its apex stayed so.
pub(crate) fn flatten_corners(
    before: &BRep,
    body: &mut BRep,
    tracking: &mut TrackingContext,
    tolerance: Tolerance,
) -> Vec<String> {
    let lin = tolerance.linear();
    let mut kept = Vec::new();
    // Clusters are found once, on the body as the blend left it; each
    // cut leaves the others' vertices where they were, found again by
    // position.
    let corners = corner_clusters(before, body, lin);
    for (cluster, points) in corners {
        let Some(tool) = corner_tool(&cluster, &points, lin, tolerance) else {
            kept.push(format!(
                "corner at {:?}: its {} corner points are not one plane",
                cluster.first(),
                points.len()
            ));
            continue;
        };
        match popkern_bool::boolean(
            popkern_bool::error::BoolOperation::Subtraction,
            body.clone(),
            tool,
            tolerance,
            Some(&mut *tracking),
        ) {
            Ok(cut) => *body = cut,
            Err(e) => kept.push(format!(
                "corner at {:?}: cut declined: {e}",
                cluster.first()
            )),
        }
    }
    kept
}

/// Each corner's cluster positions and its corner points' positions.
fn corner_clusters(before: &BRep, body: &BRep, lin: f64) -> Vec<(Vec<Point3>, Vec<Point3>)> {
    let original = |plane: (&Point3, &Vec3)| {
        before.faces.values().any(|f| {
            matches!(
                f.surface.as_ref().map(|s| &s.geometry),
                Some(SurfaceGeometry::Plane { origin, normal, .. })
                    if normal.cross(plane.1).length() <= 1e-9
                        && (*plane.0 - *origin).dot(normal).abs() <= lin
            )
        })
    };
    let bevel = |f: FaceHandle| match body
        .faces
        .get(f)
        .and_then(|d| d.surface.as_ref())
        .map(|s| &s.geometry)
    {
        Some(SurfaceGeometry::Plane { origin, normal, .. }) => !original((origin, normal)),
        _ => false,
    };
    let faces_at = |v: VertexHandle| -> Vec<FaceHandle> {
        let mut out = Vec::new();
        for he in body.half_edges.values() {
            if he.origin == v
                && let Some(l) = body.loops.get(he.loop_)
                && !out.contains(&l.face)
            {
                out.push(l.face);
            }
        }
        out
    };
    let neighbours = |v: VertexHandle| -> Vec<VertexHandle> {
        body.half_edges
            .values()
            .filter(|he| he.origin == v)
            .filter_map(|he| body.half_edges.get(he.next).map(|n| n.origin))
            .collect()
    };
    let interior: Vec<VertexHandle> = body
        .vertices
        .keys()
        .filter(|&v| {
            let fs = faces_at(v);
            fs.len() >= 3 && fs.iter().all(|&f| bevel(f))
        })
        .collect();
    let mut seen: Vec<VertexHandle> = Vec::new();
    let mut out = Vec::new();
    for &v in &interior {
        if seen.contains(&v) {
            continue;
        }
        let mut cluster = vec![v];
        let mut i = 0;
        while i < cluster.len() {
            for n in neighbours(cluster[i]) {
                if interior.contains(&n) && !cluster.contains(&n) {
                    cluster.push(n);
                }
            }
            i += 1;
        }
        seen.extend(cluster.iter().copied());
        let mut points: Vec<VertexHandle> = Vec::new();
        for &c in &cluster {
            for n in neighbours(c) {
                if !cluster.contains(&n) && !points.contains(&n) {
                    points.push(n);
                }
            }
        }
        if points.len() >= 3 {
            let at = |v: &VertexHandle| body.vertices[*v].position;
            out.push((
                cluster.iter().map(at).collect(),
                points.iter().map(at).collect(),
            ));
        }
    }
    out
}

/// The prism that cuts a corner's tip: on the plane of `points`, over
/// their polygon widened about its centroid, reaching past `cluster`.
/// `None` when the points are not coplanar.
fn corner_tool(
    cluster: &[Point3],
    points: &[Point3],
    lin: f64,
    tolerance: Tolerance,
) -> Option<BRep> {
    let n = points.len() as f64;
    let centroid = points.iter().fold(Point3::new(0.0, 0.0, 0.0), |acc, p| {
        Point3::new(acc.x + p.x / n, acc.y + p.y / n, acc.z + p.z / n)
    });
    // The plane through the first three non-collinear points, checked
    // against the rest.
    let normal = points
        .iter()
        .flat_map(|a| points.iter().map(move |b| (*a, *b)))
        .find_map(|(a, b)| (a - centroid).cross(&(b - centroid)).normalized())?;
    if points
        .iter()
        .any(|p| (*p - centroid).dot(&normal).abs() > lin)
    {
        return None;
    }
    // Facing the tip.
    let rise = |q: &Point3| (*q - centroid).dot(&normal);
    let tip = cluster
        .iter()
        .map(rise)
        .fold(0.0_f64, |m, r| if r.abs() > m.abs() { r } else { m });
    if tip.abs() <= lin {
        return None;
    }
    let normal = if tip > 0.0 { normal } else { -normal };
    let height = 3.0 * tip.abs();
    // The polygon in angular order about its centroid, widened by half.
    let (e1, e2) = popkern_core::geometry::analytic_curves::circle_basis(normal);
    let mut polygon: Vec<Point3> = points
        .iter()
        .map(|p| centroid + (*p - centroid) * 1.5)
        .collect();
    polygon.sort_by(|a, b| {
        let angle = |p: &Point3| {
            let d = *p - centroid;
            popkern_core::geometry::PortableMath::patan2(d.dot(&e2), d.dot(&e1))
        };
        angle(a).total_cmp(&angle(b))
    });
    BRep::prism(&polygon, normal, height, tolerance).ok()
}
