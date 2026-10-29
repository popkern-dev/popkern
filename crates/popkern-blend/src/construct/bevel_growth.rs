//! Growing planar bevels to their requested size (PK-1001).
//!
//! **The class.** Planar chamfers so wide that a face between bevels is
//! consumed whole: on a 6-valent cone, two adjacent slant-edge bevels
//! at d = 0.3 meet beyond the narrow face between them, so the exact
//! result has no piece of it. The per-edge pipeline trims every face
//! by its rails and cannot delete one (it declines with "no sub-face
//! containing the fillet edge"); at a smaller size the same request
//! serves exactly.
//!
//! **Growth.** The chamfer driver serves the request at a fraction of
//! its size and hands this stage the served body plus each bevel
//! face's FINAL plane. A bevel's normal does not depend on its size
//! (equal, two-distance and distance-angle chamfers all scale their
//! rails together), so each bevel plane translates along its normal by
//! an offset affine in the growth parameter `τ ∈ [0, 1]`; every other
//! face stays put. A vertex on three planes moves affinely with them,
//! and so does every edge's signed length along its fixed direction
//! (`n_F × n_G`), so the next change of topology is the first root of
//! an affine function: exact, not sampled.
//!
//! **Events,** made by Euler operators at the instant an edge's length
//! reaches zero (`bevel_growth_events`). The vanishing piece merges into
//! one vertex: a face of that edge with no area left is consumed (`kef`
//! into the neighbour across its longest side, then `kve` of its
//! vanished sides; a triangle shrinks to a point, a thin wall's side
//! face bevelled along both long sides to a segment that becomes the
//! ridge between its bevels), otherwise the lone edge is contracted
//! (`kve`). A merged vertex on more than three planes then splits into
//! the ridge tree of its planes' velocities (`split_vertex` per ridge):
//! the four-plane case is the edge flip. Served corners where more than
//! three planes meet (equal bevels meeting sharp faces) split the same
//! way at the start, or stay one vertex while their planes keep
//! concurring, which is checked over every interval between events.
//! Anything else (more vanishing edges than one, a face that would drop
//! below three sides, a ridge that does not grow) declines.
//!
//! **Dissolved vertices.** A served vertex in the middle of a straight
//! edge has no third plane to move with: where both its faces take part
//! in the growth, one `kve` merges it into an end of its edge first.
//!
//! **Invariants and postconditions.** Every moving vertex lies on its
//! three planes at every event and at the end; every edge at a moving
//! vertex is convex (one orientation sign across the region) and keeps
//! a positive length; no vertex off a bevel lies in the region its face
//! swept (a far piece of a non-convex body the growth would cut
//! through). Faces keep their handles; a consumed face is reported.
//! **Fallback:** any failure returns `Err` with the reason and the
//! caller keeps its original decline. **Tolerance:** positions and
//! lengths gate at the caller's `linear × 10`.

#[path = "bevel_growth_events.rs"]
mod events;
#[path = "bevel_growth_sweep.rs"]
mod sweep;

use super::mitre_apex_corner::{Plane, destination, outgoing};
use super::mitre_apex_scaffold::three_plane_point;
use popkern_core::geometry::curve::CurveOnEdge;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, VertexHandle};
use std::collections::{BTreeMap, BTreeSet};

/// What the growth did.
#[derive(Debug, Default, Clone, PartialEq, Eq)]
pub(crate) struct BevelGrowthReport {
    /// Edges that shrank away and flipped.
    pub(crate) edge_flips: usize,
    /// Faces the bevels consumed whole.
    pub(crate) consumed_faces: Vec<FaceHandle>,
}

/// Every face's plane as `n · x = c` with `n` agreeing with the face's
/// loop winding, and `c` at the start and the end of the growth.
pub(super) struct Motion {
    planes: BTreeMap<FaceHandle, (Vec3, f64, f64)>,
}

impl Motion {
    fn plane(&self, face: FaceHandle, tau: f64) -> Option<Plane> {
        let &(n, c0, c1) = self.planes.get(&face)?;
        let c = c0 + (c1 - c0) * tau;
        Some((Point3::new(0.0, 0.0, 0.0) + n * c, n))
    }

    fn normal(&self, face: FaceHandle) -> Option<Vec3> {
        self.planes.get(&face).map(|p| p.0)
    }

    /// How fast `face`'s offset changes per unit `τ` (0 when fixed).
    fn rate(&self, face: FaceHandle) -> Option<f64> {
        self.planes.get(&face).map(|p| p.2 - p.1)
    }

    fn grows(&self, face: FaceHandle) -> bool {
        self.planes.get(&face).is_some_and(|p| p.1 != p.2)
    }
}

/// The plane of a planar `face`, its normal oriented by its outer
/// loop's winding (Newell).
fn wound_plane(out: &BRep, face: FaceHandle) -> Option<(Vec3, f64)> {
    let surface = &out.faces.get(face)?.surface.as_ref()?.geometry;
    let SurfaceGeometry::Plane { origin, normal, .. } = surface else {
        return None;
    };
    let normal = normal.normalized()?;
    let ring = out.face_outer_loop_positions(face);
    let mut newell = Vec3::new(0.0, 0.0, 0.0);
    for (k, a) in ring.iter().enumerate() {
        let b = ring[(k + 1) % ring.len()];
        newell =
            newell + (*a - Point3::new(0.0, 0.0, 0.0)).cross(&(b - Point3::new(0.0, 0.0, 0.0)));
    }
    let n = if newell.dot(&normal) < 0.0 {
        -normal
    } else {
        normal
    };
    Some((n, n.dot(&(*origin - Point3::new(0.0, 0.0, 0.0)))))
}

/// The distinct faces around `vertex`.
pub(super) fn faces_at(out: &BRep, vertex: VertexHandle) -> BTreeSet<FaceHandle> {
    outgoing(out, vertex)
        .into_iter()
        .filter_map(|he| out.loops.get(out.half_edges.get(he)?.loop_).map(|l| l.face))
        .collect()
}

/// The two faces of `edge`, left of its primary half-edge first, and
/// its origin and destination.
pub(super) fn edge_ends(
    out: &BRep,
    edge: EdgeHandle,
) -> Option<(FaceHandle, FaceHandle, VertexHandle, VertexHandle)> {
    let he = out.edges.get(edge)?.half_edge;
    let h = out.half_edges.get(he)?;
    let twin = out.half_edges.get(h.twin)?;
    Some((
        out.loops.get(h.loop_)?.face,
        out.loops.get(twin.loop_)?.face,
        h.origin,
        destination(out, he)?,
    ))
}

/// Where `vertex` is at `tau`: fixed unless a growing face is around
/// it, then the point of its three planes.
fn position(out: &BRep, motion: &Motion, vertex: VertexHandle, tau: f64) -> Result<Point3, String> {
    let faces = faces_at(out, vertex);
    if !faces.iter().any(|&f| motion.grows(f)) {
        return out
            .vertices
            .get(vertex)
            .map(|v| v.position)
            .ok_or_else(|| "stale vertex".into());
    }
    let planes: Vec<Plane> = faces
        .iter()
        .map(|&f| {
            motion
                .plane(f, tau)
                .ok_or("a face at a moving vertex is not planar")
        })
        .collect::<Result<_, _>>()?;
    if planes.len() < 3 {
        return Err(format!("a moving vertex lies on {} planes", planes.len()));
    }
    // More than three planes concur only where `concurrent_corners`
    // proved it: any well-conditioned three of them give the point.
    let triple =
        |a: usize, b: usize, c: usize| planes[a].1.dot(&planes[b].1.cross(&planes[c].1)).abs();
    let mut best = (0, 1, 2);
    for a in 0..planes.len() {
        for b in a + 1..planes.len() {
            for c in b + 1..planes.len() {
                if triple(a, b, c) > triple(best.0, best.1, best.2) {
                    best = (a, b, c);
                }
            }
        }
    }
    three_plane_point(&planes[best.0], &planes[best.1], &planes[best.2])
        .ok_or_else(|| "a moving vertex's planes share a line".into())
}

/// Every moving vertex on more than three planes must stay on all of
/// them until the next event: its point and the planes' offsets are
/// affine in `τ`, so meeting them at `from` and `to` means meeting them
/// throughout (the symmetric corner of two equal bevels and two sharp
/// faces).
fn concurrent_corners(
    out: &BRep,
    motion: &Motion,
    from: f64,
    to: f64,
    gate: f64,
) -> Result<(), String> {
    for v in out.vertices.keys().filter(|&v| moves(out, motion, v)) {
        let faces = faces_at(out, v);
        if faces.len() <= 3 {
            continue;
        }
        for t in [from, to] {
            let p = position(out, motion, v, t)?;
            for &f in &faces {
                let (o, n) = motion
                    .plane(f, t)
                    .ok_or("a face at a moving vertex is not planar")?;
                if (p - o).dot(&n).abs() > gate {
                    return Err(format!(
                        "a moving vertex on {} planes splits apart (a ridge tree, not a growth)",
                        faces.len()
                    ));
                }
            }
        }
    }
    Ok(())
}

/// Dissolve the vertices in the middle of a straight edge (two edges,
/// collinear, the same two faces) whose faces both carry a moving
/// vertex: such a vertex has no third plane to move with. One `kve`
/// merges it into an end of its edge.
fn dissolve_mid_edge_vertices(out: &mut BRep, motion: &Motion, gate: f64) -> Result<(), String> {
    loop {
        let busy: BTreeSet<FaceHandle> = out
            .vertices
            .keys()
            .filter(|&w| moves(out, motion, w))
            .flat_map(|w| faces_at(out, w))
            .collect();
        let found = out.vertices.keys().find_map(|v| {
            let fan = outgoing(out, v);
            let [a, b] = fan[..] else {
                return None;
            };
            let (pa, pb) = (destination(out, a)?, destination(out, b)?);
            let p = out.vertices.get(v)?.position;
            let (qa, qb) = (
                out.vertices.get(pa)?.position,
                out.vertices.get(pb)?.position,
            );
            let line = (qb - qa).normalized()?;
            let off = (p - qa) - line * (p - qa).dot(&line);
            let faces = faces_at(out, v);
            (off.length() <= gate && faces.len() == 2 && faces.iter().all(|f| busy.contains(f)))
                .then(|| out.half_edges.get(a).map(|h| h.twin))?
        });
        let Some(toward) = found else {
            return Ok(());
        };
        out.kve(toward).map_err(|e| format!("dissolve kve: {e}"))?;
    }
}

/// Whether `vertex` moves (a growing face is around it).
fn moves(out: &BRep, motion: &Motion, vertex: VertexHandle) -> bool {
    faces_at(out, vertex).iter().any(|&f| motion.grows(f))
}

/// The edges with a moving end.
fn tracked_edges(out: &BRep, motion: &Motion) -> Vec<EdgeHandle> {
    out.edges
        .keys()
        .filter(|&e| {
            edge_ends(out, e)
                .is_some_and(|(_, _, o, d)| moves(out, motion, o) || moves(out, motion, d))
        })
        .collect()
}

/// `edge`'s length along its convex direction `n_F × n_G` at `tau`,
/// times the region's orientation sign `sigma`.
pub(super) fn signed_length(
    out: &BRep,
    motion: &Motion,
    edge: EdgeHandle,
    tau: f64,
    sigma: f64,
) -> Result<f64, String> {
    let (f, g, o, d) = edge_ends(out, edge).ok_or("stale edge")?;
    let (nf, ng) = (
        motion
            .normal(f)
            .ok_or("a face at a moving edge is not planar")?,
        motion
            .normal(g)
            .ok_or("a face at a moving edge is not planar")?,
    );
    let dir = nf
        .cross(&ng)
        .normalized()
        .ok_or("an edge between parallel faces")?;
    let span = position(out, motion, d, tau)? - position(out, motion, o, tau)?;
    Ok(sigma * span.dot(&dir))
}

/// Move every moving vertex to its position at `tau`.
fn place(out: &mut BRep, motion: &Motion, tau: f64) -> Result<(), String> {
    let moving: Vec<VertexHandle> = out
        .vertices
        .keys()
        .filter(|&v| moves(out, motion, v))
        .collect();
    for v in moving {
        let p = position(out, motion, v, tau)?;
        out.vertices.get_mut(v).ok_or("stale vertex")?.position = p;
    }
    Ok(())
}

/// Grow the bevel faces of `out` to their `end` planes (see the module
/// docs). `out` is left unspecified on `Err`.
pub(crate) fn grow_bevels(
    out: &mut BRep,
    growing: &[(FaceHandle, Plane)],
    gate: f64,
) -> Result<BevelGrowthReport, String> {
    let mut planes = BTreeMap::new();
    for face in out.faces.keys() {
        if let Some((n, c)) = wound_plane(out, face) {
            planes.insert(face, (n, c, c));
        }
    }
    for &(face, (o, n)) in growing {
        let entry = planes.get_mut(&face).ok_or("a bevel face is not planar")?;
        if (entry.0.dot(&n).abs() - 1.0).abs() > 1e-9 {
            return Err("a bevel's end plane is not parallel to its face".into());
        }
        entry.2 = entry.0.dot(&(o - Point3::new(0.0, 0.0, 0.0)));
    }
    let motion = Motion { planes };
    dissolve_mid_edge_vertices(out, &motion, gate)?;
    // One orientation sign across the region: every edge at a moving
    // vertex convex.
    let mut sigma = 0.0f64;
    for e in tracked_edges(out, &motion) {
        let raw = signed_length(out, &motion, e, 0.0, 1.0)?;
        if raw.abs() <= gate || (sigma != 0.0 && raw.signum() != sigma) {
            return Err("the growth region has a concave or vanishing edge".into());
        }
        sigma = raw.signum();
    }
    // A served corner on more than three planes (equal bevels meeting
    // sharp faces) splits now if its planes part as they grow.
    let crowded: Vec<VertexHandle> = out
        .vertices
        .keys()
        .filter(|&v| moves(out, &motion, v) && faces_at(out, v).len() > 3)
        .collect();
    for v in crowded {
        events::expand(out, &motion, v, 0.0, sigma, gate)?;
    }
    let mut report = BevelGrowthReport::default();
    let mut swept: BTreeMap<FaceHandle, Vec<Point3>> = BTreeMap::new();
    let record = |out: &BRep, swept: &mut BTreeMap<FaceHandle, Vec<Point3>>| {
        for &(face, _) in growing {
            if out.faces.contains_key(face) {
                swept
                    .entry(face)
                    .or_default()
                    .extend(out.face_outer_loop_positions(face));
            }
        }
    };
    record(out, &mut swept);
    let start = out.clone();
    let mut tau = 0.0f64;
    let budget = 4 * out.edges.len() + 8;
    for _ in 0..budget {
        let mut next: Option<(f64, EdgeHandle)> = None;
        for e in tracked_edges(out, &motion) {
            let (l0, l1) = (
                signed_length(out, &motion, e, tau, sigma)?,
                signed_length(out, &motion, e, 1.0, sigma)?,
            );
            if l1.abs() <= gate {
                return Err("an edge vanishes exactly at the requested size".into());
            }
            if l1 < 0.0 {
                let t = tau + (1.0 - tau) * (l0 / (l0 - l1)).max(0.0);
                if next.is_none_or(|(best, _)| t < best) {
                    next = Some((t, e));
                }
            }
        }
        concurrent_corners(out, &motion, tau, next.map_or(1.0, |n| n.0), gate)?;
        let Some((t, edge)) = next else {
            place(out, &motion, 1.0)?;
            record(out, &mut swept);
            finish(out, &motion, growing, gate)?;
            sweep::sweep_guard(&start, out, &motion, &swept, gate)?;
            return Ok(report);
        };
        place(out, &motion, t)?;
        record(out, &mut swept);
        events::resolve(out, &motion, edge, t, sigma, gate, &mut report)?;
        tau = t;
    }
    Err("the growth did not settle within its event budget".into())
}

/// Set the bevel surfaces to their end planes and re-author every edge
/// at a moving vertex as its exact segment; then check every moving
/// vertex on its planes.
fn finish(
    out: &mut BRep,
    motion: &Motion,
    growing: &[(FaceHandle, Plane)],
    gate: f64,
) -> Result<(), String> {
    for &(face, (o, _)) in growing {
        let Some(f) = out.faces.get_mut(face) else {
            continue;
        };
        if let Some(surface) = f.surface.as_mut()
            && let SurfaceGeometry::Plane { origin, normal, .. } = &mut surface.geometry
        {
            let shift = (o - *origin).dot(normal) / normal.dot(normal);
            *origin = *origin + *normal * shift;
        }
    }
    for e in tracked_edges(out, motion) {
        let (_, _, o, d) = edge_ends(out, e).ok_or("stale edge")?;
        let (po, pd) = (
            out.vertices.get(o).ok_or("stale vertex")?.position,
            out.vertices.get(d).ok_or("stale vertex")?.position,
        );
        let segment = CurveOnEdge::line_segment(po, pd, gate / 10.0).ok_or("degenerate segment")?;
        out.edges.get_mut(e).ok_or("stale edge")?.curve = Some(segment);
    }
    for v in out.vertices.keys().filter(|&v| moves(out, motion, v)) {
        let p = out.vertices[v].position;
        for f in faces_at(out, v) {
            let (o, n) = motion
                .plane(f, 1.0)
                .ok_or("a face at a moving vertex is not planar")?;
            if (p - o).dot(&n).abs() > gate {
                return Err("a moving vertex is off a face after the growth".into());
            }
        }
    }
    Ok(())
}
