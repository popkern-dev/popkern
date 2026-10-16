//! The exact rolling-ball corner (`CornerConstruction::RollingBall`,
//! PK-861).
//!
//! Where three rounds of one radius `r` meet on planar supports, the
//! ball tangent to all three supports has its centre `C` on all three
//! round axes, and the corner is exactly the spherical triangle of that
//! ball between the rounds' end cross-sections: each round ends on the
//! plane through `C` perpendicular to its edge (its cross-section there
//! is a great circle of the ball), and the triangle's vertices are the
//! ball's feet on the supports. G1 to every round by construction.
//!
//! The setback construction already builds this topology at such a
//! corner: its setback is the farther rail crossing, and both rails of
//! a round cross their neighbours at the ball's feet, which lie on the
//! round's own cross-section through `C` (`P_F = C - r n_F` and the
//! edge lies in `F`), so the bridges vanish and the cap boundary is the
//! three end arcs. Measured on a tetrahedron's corners: every cap
//! boundary vertex within 4e-4·r of the ball; the KP20 cap's own centre
//! 6e-3·r inside it. This pass replaces that cap by the ball:
//! eligibility is checked on the served body (every face across the
//! cap boundary a round cylinder of one radius, the axes concurrent at
//! `C`, every boundary vertex and edge sample on the ball), the cap
//! patches merge into one face (Euler `kef`, then `kev` for the spokes'
//! centre), and that face takes the sphere and fresh pcurves. A corner
//! that is not eligible keeps its setback cap, or fails the fillet when
//! the preference is strict.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry, SurfaceOnFace};
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle, HalfEdgeHandle, VertexHandle};
use std::collections::{BTreeMap, BTreeSet};

use crate::types::{BlendResult, ContinuityLevel, CornerMethod, CornerPatchInfo};

/// One corner's ball: centre, radius, the round faces around it.
struct Ball {
    center: Point3,
    radius: f64,
    rounds: usize,
}

fn loop_half_edges(brep: &BRep, face: FaceHandle) -> Vec<HalfEdgeHandle> {
    let Some(f) = brep.faces.get(face) else {
        return Vec::new();
    };
    let mut out = Vec::new();
    for lh in std::iter::once(f.outer_loop).chain(f.inner_loops.iter().copied()) {
        let Some(start) = brep.loops.get(lh).map(|l| l.half_edge) else {
            continue;
        };
        let mut he = start;
        for _ in 0..=brep.half_edges.len() {
            out.push(he);
            let Some(next) = brep.half_edges.get(he).map(|h| h.next) else {
                break;
            };
            he = next;
            if he == start {
                break;
            }
        }
    }
    out
}

fn face_of(brep: &BRep, he: HalfEdgeHandle) -> Option<FaceHandle> {
    let h = brep.half_edges.get(he)?;
    Some(brep.loops.get(h.loop_)?.face)
}

/// The point nearest all `axes` (least squares), solving
/// `Σ (I − a aᵀ) x = Σ (I − a aᵀ) o`.
fn nearest_to_axes(axes: &[(Point3, Vec3)]) -> Option<Point3> {
    let mut m = [[0.0f64; 3]; 3];
    let mut b = [0.0f64; 3];
    for (o, a) in axes {
        let a = a.normalized()?;
        let av = [a.x, a.y, a.z];
        let ov = [o.x, o.y, o.z];
        for r in 0..3 {
            for c in 0..3 {
                let p = f64::from(u8::from(r == c)) - av[r] * av[c];
                m[r][c] += p;
                b[r] += p * ov[c];
            }
        }
    }
    let det = |m: [[f64; 3]; 3]| {
        m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
            - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
    };
    let d = det(m);
    if d.abs() < 1e-9 {
        return None;
    }
    let mut x = [0.0; 3];
    for (col, xc) in x.iter_mut().enumerate() {
        let mut mc = m;
        for row in 0..3 {
            mc[row][col] = b[row];
        }
        *xc = det(mc) / d;
    }
    Some(Point3::new(x[0], x[1], x[2]))
}

fn axis_distance(p: &Point3, (o, a): &(Point3, Vec3)) -> f64 {
    let a = a.normalized().unwrap_or(*a);
    let d = *p - *o;
    (d - a * d.dot(&a)).length()
}

/// The ball of a corner whose cap is `cap`, or why there is none.
fn corner_ball(brep: &BRep, cap: &BTreeSet<FaceHandle>, tol: f64) -> Result<Ball, String> {
    let mut rounds: BTreeMap<FaceHandle, (Point3, Vec3, f64)> = BTreeMap::new();
    let mut boundary = Vec::new();
    for &face in cap {
        for he in loop_half_edges(brep, face) {
            let twin = brep
                .half_edges
                .get(he)
                .map(|h| h.twin)
                .ok_or("stale half-edge")?;
            let other = face_of(brep, twin).ok_or("stale twin")?;
            if cap.contains(&other) {
                continue;
            }
            boundary.push(he);
            let geometry = brep
                .faces
                .get(other)
                .and_then(|f| f.surface.as_ref())
                .map(|s| &s.geometry)
                .ok_or("a cap neighbour has no surface")?;
            let SurfaceGeometry::Cylinder {
                origin,
                axis,
                radius,
            } = geometry
            else {
                return Err(format!(
                    "the cap borders a {} face, not only round cylinders",
                    geometry.kind_name()
                ));
            };
            rounds.insert(other, (*origin, *axis, *radius));
        }
    }
    if rounds.len() != 3 {
        return Err(format!("the cap borders {} rounds, not 3", rounds.len()));
    }
    let radii: Vec<f64> = rounds.values().map(|r| r.2).collect();
    let radius = radii.iter().sum::<f64>() / radii.len() as f64;
    let gate = tol.max(1e-7 * radius);
    if radii.iter().any(|r| (r - radius).abs() > gate) {
        return Err(format!("the rounds' radii differ: {radii:?}"));
    }
    let axes: Vec<(Point3, Vec3)> = rounds.values().map(|r| (r.0, r.1)).collect();
    let center = nearest_to_axes(&axes).ok_or("the round axes are parallel")?;
    let off_axis = axes
        .iter()
        .map(|a| axis_distance(&center, a))
        .fold(0.0f64, f64::max);
    if off_axis > gate {
        return Err(format!("the round axes miss one centre by {off_axis:.3e}"));
    }
    // Every cap boundary vertex and edge sample on the ball.
    let mut worst = 0.0f64;
    for he in boundary {
        let h = brep.half_edges.get(he).ok_or("stale half-edge")?;
        let p = brep.vertices.get(h.origin).ok_or("stale vertex")?.position;
        worst = worst.max((p.distance_to(&center) - radius).abs());
        if let Some(curve) = brep.edges.get(h.edge).and_then(|e| e.curve.as_ref()) {
            for k in 1..8 {
                let t = curve.t_start + (curve.t_end - curve.t_start) * k as f64 / 8.0;
                let q = curve.geometry.point_at(t);
                worst = worst.max((q.distance_to(&center) - radius).abs());
            }
        }
    }
    // The setback's end arcs are fitted to the ball's cross-sections:
    // they hold it to their own fit accuracy.
    let on_ball = (100.0 * tol).max(1e-5 * radius);
    if worst > on_ball {
        return Err(format!(
            "the cap boundary leaves the ball by {worst:.3e} (band {on_ball:.1e})"
        ));
    }
    Ok(Ball {
        center,
        radius,
        rounds: rounds.len(),
    })
}

/// Merge the cap faces into one: `kef` across every inner edge, then
/// `kev` the spurs (the KP20 spokes' centre). Returns the survivor.
fn merge_cap(brep: &mut BRep, cap: &BTreeSet<FaceHandle>) -> Result<FaceHandle, String> {
    let mut faces: BTreeSet<FaceHandle> = cap.clone();
    loop {
        let inner = faces.iter().find_map(|&f| {
            loop_half_edges(brep, f).into_iter().find(|&he| {
                let twin = brep.half_edges.get(he).map(|h| h.twin);
                let other = twin.and_then(|t| face_of(brep, t));
                other.is_some_and(|o| o != f && faces.contains(&o))
            })
        });
        let Some(he) = inner else {
            break;
        };
        let killed = brep.kef(he).map_err(|e| format!("kef cap seam: {e}"))?;
        faces.remove(&killed.removed_face);
    }
    let [survivor] = faces.iter().copied().collect::<Vec<_>>()[..] else {
        return Err(format!("the cap merged into {} faces", faces.len()));
    };
    // Spurs: edges whose both sides now lie in the survivor's loop.
    loop {
        let spur = loop_half_edges(brep, survivor).into_iter().find(|&he| {
            let Some(h) = brep.half_edges.get(he) else {
                return false;
            };
            let Some(t) = brep.half_edges.get(h.twin) else {
                return false;
            };
            // The spur's tip: the destination with valence 1.
            t.loop_ == h.loop_
                && brep
                    .half_edges
                    .values()
                    .filter(|x| x.origin == t.origin)
                    .count()
                    == 1
        });
        let Some(he) = spur else {
            break;
        };
        brep.kev(he).map_err(|e| format!("kev cap centre: {e}"))?;
    }
    Ok(survivor)
}

/// Replace each eligible corner's setback cap by its rolling ball.
/// `Err((vertex, reason))` names a corner that is not eligible when
/// `strict`; otherwise such corners keep their caps. Returns how many
/// corners became balls.
pub(crate) fn settle_rolling_ball_corners(
    result: &mut BlendResult,
    tolerance: &Tolerance,
    strict: bool,
) -> Result<usize, (VertexHandle, String)> {
    let tol = tolerance.linear().max(1e-9);
    let mut by_corner: BTreeMap<VertexHandle, Vec<usize>> = BTreeMap::new();
    for (k, patch) in result.corner_patches.iter().enumerate() {
        by_corner.entry(patch.vertex).or_default().push(k);
    }
    let mut settled = 0;
    let mut drop: BTreeSet<usize> = BTreeSet::new();
    for (vertex, patches) in by_corner {
        let cap: BTreeSet<FaceHandle> = patches
            .iter()
            .map(|&k| result.corner_patches[k].face)
            .collect();
        // A corner already served as the exact ball (cube-class corners:
        // the vertex blend's surface fix) stays as it is.
        let already_ball = cap.iter().all(|&f| {
            matches!(
                result
                    .brep
                    .faces
                    .get(f)
                    .and_then(|x| x.surface.as_ref())
                    .map(|s| &s.geometry),
                Some(SurfaceGeometry::Sphere { .. })
            )
        });
        if already_ball {
            continue;
        }
        let ball = match corner_ball(&result.brep, &cap, tol) {
            Ok(ball) => ball,
            Err(why) if strict => return Err((vertex, why)),
            Err(_) => continue,
        };
        // Outward normal sense from the cap as served.
        let first = *cap
            .iter()
            .next()
            .ok_or((vertex, "an empty cap".to_string()))?;
        let probe = loop_half_edges(&result.brep, first)
            .first()
            .and_then(|&he| result.brep.half_edges.get(he))
            .and_then(|h| result.brep.vertices.get(h.origin))
            .map(|v| v.position)
            .ok_or((vertex, "a cap without vertices".to_string()))?;
        let outward_from_center = {
            let sof = result.brep.faces[first].surface.as_ref();
            sof.and_then(|s| {
                let (u, v) = s.geometry.closest_point(&probe)?;
                let n = s.geometry.normal_at(u, v) * if s.orientation { 1.0 } else { -1.0 };
                Some(n.dot(&(probe - ball.center)) > 0.0)
            })
            .unwrap_or(true)
        };
        let mut staged = result.brep.clone();
        let face = merge_cap(&mut staged, &cap).map_err(|why| (vertex, why))?;
        let sphere = SurfaceGeometry::Sphere {
            center: ball.center,
            radius: ball.radius,
            axis: Vec3::Z,
        };
        let chart_outward = sphere
            .closest_point(&probe)
            .map(|(u, v)| sphere.normal_at(u, v).dot(&(probe - ball.center)) > 0.0)
            .unwrap_or(true);
        staged.faces[face].surface = Some(SurfaceOnFace {
            geometry: sphere,
            orientation: chart_outward == outward_from_center,
        });
        for he in loop_half_edges(&staged, face) {
            if let Some(h) = staged.half_edges.get_mut(he) {
                h.pcurve = None;
            }
        }
        staged.face_use_pcurves.retain(|key, _| key.face != face);
        staged.recompute_pcurves_for_lenient(&[face]);
        if let Err(e) = staged.validate() {
            if strict {
                return Err((vertex, format!("the ball corner fails validation: {e}")));
            }
            continue;
        }
        result.brep = staged;
        let keep = patches[0];
        drop.extend(patches.iter().skip(1).copied());
        let info: &mut CornerPatchInfo = &mut result.corner_patches[keep];
        info.face = face;
        info.n_sides = ball.rounds;
        info.method = CornerMethod::RollingBallSphere;
        info.achieved_continuity = ContinuityLevel::G1;
        info.cap_sphere = Some((ball.center, ball.radius));
        info.setback_cap = false;
        info.exact_boundary_cap = false;
        info.sub_patch_provenance = None;
        settled += 1;
    }
    let mut k = 0;
    result.corner_patches.retain(|_| {
        let kept = !drop.contains(&k);
        k += 1;
        kept
    });
    Ok(settled)
}
