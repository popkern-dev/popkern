//! Tolerance-scoped T-junction zippering (PK-208, completed by
//! PK-334 as curve-parameter station merging).
//!
//! Faces sample a shared BRep edge from several legitimate sources —
//! canonical per-edge stations, encroachment enrichment, clip-cut
//! T-vertices, legacy grid stations — and the sewing pass lands every
//! chain ON the shared edge curve. Two chains with different station
//! sets therefore sit on the SAME curve at STAGGERED positions:
//! position-welding can never pair them, and chord-distance matching
//! is the wrong geometry outright (an on-curve station sits a full
//! chord-sagitta off the opposing chain's chord — measured 1.5e-2 on
//! a 0.36 rad rim chord, far beyond any honest matching band).
//!
//! This pass merges the station sets by CURVE PARAMETER: an open mesh
//! edge whose endpoints lie on a two-sided BRep edge's curve is split
//! at every opposing open-boundary vertex that also lies on that
//! curve, strictly inside the open edge's parameter interval. Both
//! chains end up sharing every station and the weld closes.
//!
//! # Why the scoping is the whole design
//!
//! Zippering unconditionally would blind the watertightness gates: a
//! REAL geometric hole (a flush-boolean gap, a dropped face) would be
//! silently stitched shut, and "watertight" would stop being
//! evidence. The zipper therefore only acts where the BRep TOPOLOGY
//! says the mesh must be closed: on open mesh edges lying along a
//! TWO-SIDED BRep edge (distinct faces on both half-edges). There,
//! openness can only be a sampling artifact. Membership on the curve
//! is held to a band covering the discretization's own accepted sag
//! (the tessellation chord height) or 1.5x a RECORDED entity
//! tolerance — the check-9 contract. A hole whose opposing chain is
//! absent or off the curve stays open — honest.
//!
//! The pass reports how many junctions it closed via
//! [`crate::mesh::TriMesh::zipped_t_junctions`]; a consumer that
//! wants to gate on "watertight without repair" can assert that
//! count is zero.

use std::collections::HashMap;

use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::CurveEval;
use popkern_core::topology::BRep;

use crate::mesh::TriMesh;

/// Weld-grid resolution — must match the welded-open-edge census the
/// example gates use.
const WELD: f64 = 1.0e7;
/// Maximum zip rounds; each round can expose junctions the previous
/// round's splits created (a triangle carries at most one fan per
/// round, and a corner triangle can own two open edges).
const MAX_ROUNDS: usize = 8;
/// Interior window as a fraction of the open edge's parameter span: a
/// candidate at (or numerically at) an endpoint is already welded.
const END_FRACTION: f64 = 1.0e-6;
/// Fan cap per open edge per round; leftover stations close in later
/// rounds.
const MAX_FAN: usize = 32;

/// Default matching band for a chain's distance OFF the shared curve
/// with no recorded tolerance: covers sampling/sewing offsets while
/// staying far under any real geometric gap worth catching.
const DEFAULT_BAND: f64 = 2.0e-4;

/// Sampled polyline of one two-sided BRep edge, with cumulative
/// arc-length stations for parameter location.
struct EdgePolyline {
    pts: Vec<Point3>,
    /// `cum[i]` = arc length from `pts[0]` to `pts[i]`.
    cum: Vec<f64>,
    band: f64,
    /// Curve endpoints coincide (full rim): parameter intervals may
    /// wrap through the polyline's start.
    closed: bool,
}

/// Close T-junction cracks along two-sided BRep edges. Returns the
/// number of junction splits performed.
pub(crate) fn zip_recorded_tolerance_cracks(
    brep: &BRep,
    mesh: &mut TriMesh,
    chord_height: f64,
) -> usize {
    let mut closed_edges: Vec<EdgePolyline> = Vec::new();
    for (_eh, edge) in brep.edges.iter() {
        let two_sided = (|| {
            let he = brep.half_edges.get(edge.half_edge)?;
            let twin = brep.half_edges.get(he.twin)?;
            let fa = brep.loops.get(he.loop_)?.face;
            let fb = brep.loops.get(twin.loop_)?.face;
            Some(fa != fb)
        })()
        .unwrap_or(false);
        if !two_sided {
            continue;
        }
        let Some(curve) = edge.curve.as_ref() else {
            continue;
        };
        if !(curve.t_start.is_finite() && curve.t_end.is_finite()) {
            continue;
        }
        let ce: &dyn CurveEval = &curve.geometry;
        let pts: Vec<Point3> = (0..=64)
            .map(|i| ce.point_at(curve.t_start + (curve.t_end - curve.t_start) * i as f64 / 64.0))
            .collect();
        let mut cum = Vec::with_capacity(pts.len());
        let mut total = 0.0;
        cum.push(0.0);
        for w in pts.windows(2) {
            total += w[0].distance_to(&w[1]);
            cum.push(total);
        }
        if total < 1e-12 {
            continue;
        }
        // The band must cover the DISCRETIZATION's own accepted sag:
        // chain stations sit up to the tessellation chord height off
        // the exact curve (PK-299, measured 4.2e-3 on a rim sampled
        // at ~10° spans). Geometry finer than the chord height is
        // below the mesh's resolution contract, so openness within it
        // along a two-sided edge is a sampling artifact by
        // construction.
        let band = DEFAULT_BAND
            .max(chord_height)
            .max(edge.tolerance.map(|t| 1.5 * t).unwrap_or(0.0));
        let closed = pts
            .first()
            .zip(pts.last())
            .map(|(a, b)| a.distance_to(b) < 1e-9)
            .unwrap_or(false);
        closed_edges.push(EdgePolyline {
            pts,
            cum,
            band,
            closed,
        });
    }
    if closed_edges.is_empty() {
        return 0;
    }

    let mut zipped = 0usize;
    for round in 0..MAX_ROUNDS {
        let splits = zip_round(mesh, &closed_edges);
        if std::env::var_os("PK_ZIP_TRACE").is_some() {
            eprintln!("  [zip] round {round}: {splits} splits");
        }
        if splits == 0 {
            break;
        }
        zipped += splits;
    }
    zipped
}

/// Locate `p` on `edge`: `(arc position, distance off the polyline)`.
fn locate_on(edge: &EdgePolyline, p: &Point3) -> (f64, f64) {
    let mut best = (0.0, f64::MAX);
    for (i, w) in edge.pts.windows(2).enumerate() {
        let (f, foot) = project_point_segment(p, &w[0], &w[1]);
        let d = foot.distance_to(p);
        if d < best.1 {
            best = (edge.cum[i] + (edge.cum[i + 1] - edge.cum[i]) * f, d);
        }
    }
    best
}

/// Best on-curve location of `p` across all polylines within band:
/// `(edge index, arc position)`.
fn locate(closed_edges: &[EdgePolyline], p: &Point3) -> Option<(usize, f64)> {
    let mut best: Option<(usize, f64, f64)> = None;
    for (ei, edge) in closed_edges.iter().enumerate() {
        let (s, d) = locate_on(edge, p);
        if d <= edge.band && best.map(|(_, _, bd)| d < bd).unwrap_or(true) {
            best = Some((ei, s, d));
        }
    }
    best.map(|(ei, s, _)| (ei, s))
}

fn zip_round(mesh: &mut TriMesh, closed_edges: &[EdgePolyline]) -> usize {
    // Weld census: open edges are used by exactly one triangle.
    let mut weld: HashMap<(i64, i64, i64), usize> = HashMap::new();
    let mut weld_id = |p: Point3| -> usize {
        let k = (
            (p.x * WELD).round() as i64,
            (p.y * WELD).round() as i64,
            (p.z * WELD).round() as i64,
        );
        let next = weld.len();
        *weld.entry(k).or_insert(next)
    };
    let mut ids: Vec<[usize; 3]> = Vec::with_capacity(mesh.triangles.len());
    for t in &mesh.triangles {
        ids.push([
            weld_id(mesh.vertices[t[0]].position),
            weld_id(mesh.vertices[t[1]].position),
            weld_id(mesh.vertices[t[2]].position),
        ]);
    }
    let mut edge_use: HashMap<(usize, usize), usize> = HashMap::new();
    for w in &ids {
        if w[0] == w[1] || w[1] == w[2] || w[0] == w[2] {
            continue;
        }
        for (a, b) in [(w[0], w[1]), (w[1], w[2]), (w[2], w[0])] {
            *edge_use.entry((a.min(b), a.max(b))).or_insert(0) += 1;
        }
    }
    // Open boundary vertices (weld-id -> a representative position).
    let mut open_vertex_pos: HashMap<usize, Point3> = HashMap::new();
    let mut open_edges: Vec<(usize, usize, usize)> = Vec::new(); // (tri, corner, corner+1)
    for (ti, t) in mesh.triangles.iter().enumerate() {
        let w = ids[ti];
        if w[0] == w[1] || w[1] == w[2] || w[0] == w[2] {
            continue;
        }
        for k in 0..3 {
            let (a, b) = (w[k], w[(k + 1) % 3]);
            if edge_use.get(&(a.min(b), a.max(b))) == Some(&1) {
                open_edges.push((ti, k, (k + 1) % 3));
                open_vertex_pos
                    .entry(a)
                    .or_insert(mesh.vertices[t[k]].position);
                open_vertex_pos
                    .entry(b)
                    .or_insert(mesh.vertices[t[(k + 1) % 3]].position);
            }
        }
    }
    if open_edges.is_empty() {
        return 0;
    }
    // Canonical candidate order: `open_vertex_pos` is a HashMap, and
    // draining it in map order made split choices hash-ordered
    // (measured call-to-call triangle-count drift, PK-277). Weld ids
    // are insertion-ordered, so sorting by id is deterministic; the
    // per-edge fan then orders by arc position with the id as the
    // tie-break.
    let mut candidates: Vec<(usize, Point3)> = open_vertex_pos.into_iter().collect();
    candidates.sort_unstable_by_key(|&(wv, _)| wv);
    // On-curve location per candidate, computed once.
    let located: Vec<Option<(usize, f64)>> = candidates
        .iter()
        .map(|(_, pv)| locate(closed_edges, pv))
        .collect();

    // One fan per triangle per round keeps the bookkeeping simple;
    // subsequent rounds finish the job (a corner triangle can own two
    // open edges).
    let mut split_tris: std::collections::BTreeSet<usize> = std::collections::BTreeSet::new();
    let mut new_triangles: Vec<([usize; 3], popkern_core::topology::FaceHandle)> = Vec::new();
    let mut splits = 0usize;
    let mut stats = (0usize, 0usize, 0usize);
    for (ti, ka, kb) in open_edges {
        stats.0 += 1;
        if split_tris.contains(&ti) {
            continue;
        }
        let tri = mesh.triangles[ti];
        let (ia, ib) = (tri[ka], tri[kb]);
        let (pa, pb) = (mesh.vertices[ia].position, mesh.vertices[ib].position);
        let (wa, wb) = (ids[ti][ka], ids[ti][kb]);
        if pa.distance_to(&pb) < 1e-12 {
            continue;
        }
        // Both endpoints must lie ON the same two-sided edge's curve
        // (within its band) — the topology then says the mesh must
        // close here, and the shared curve parameterizes the merge.
        let (Some((ea, s_a)), Some((eb, s_b))) =
            (locate(closed_edges, &pa), locate(closed_edges, &pb))
        else {
            stats.1 += 1;
            trace_skip("off-curve endpoint", &pa, &pb);
            continue;
        };
        if ea != eb {
            stats.1 += 1;
            trace_skip("endpoints on different edges", &pa, &pb);
            continue;
        }
        let edge = &closed_edges[ea];
        let total = edge.cum.last().copied().unwrap_or(0.0);
        // Parameter interval from a to b, taking the SHORTER arc on
        // closed curves (open mesh edges are single chords, always
        // short relative to a full rim).
        let direct = (s_b - s_a).abs();
        let wrapping = edge.closed && direct > 0.5 * total;
        let span = if wrapping { total - direct } else { direct };
        if span < 1e-12 {
            continue;
        }
        let increasing = if wrapping { s_b < s_a } else { s_b > s_a };
        let key_of = |s: f64| -> f64 {
            let raw = if increasing { s - s_a } else { s_a - s };
            if edge.closed {
                raw.rem_euclid(total)
            } else {
                raw
            }
        };
        // Every opposing open vertex on the same curve strictly inside
        // the interval splits the edge — merged in one fan, ordered by
        // arc position.
        let mut fan: Vec<(f64, usize, Point3)> = Vec::new();
        for (ci, (wv, pv)) in candidates.iter().enumerate() {
            if *wv == wa || *wv == wb {
                continue;
            }
            let Some((ec, s_c)) = located[ci] else {
                continue;
            };
            if ec != ea {
                continue;
            }
            let key = key_of(s_c);
            if key <= span * END_FRACTION || key >= span * (1.0 - END_FRACTION) {
                continue;
            }
            if pv.distance_to(&pa) < 1e-9 || pv.distance_to(&pb) < 1e-9 {
                continue;
            }
            fan.push((key, *wv, *pv));
        }
        if fan.is_empty() {
            stats.2 += 1;
            trace_skip("no on-curve candidate", &pa, &pb);
            continue;
        }
        fan.sort_by(|x, y| {
            x.0.partial_cmp(&y.0)
                .unwrap_or(std::cmp::Ordering::Equal)
                .then(x.1.cmp(&y.1))
        });
        fan.dedup_by_key(|e| e.1);
        fan.truncate(MAX_FAN);
        // Split: replace (a, b, c) with the fan (a, v1, c),
        // (v1, v2, c), …, (vk, b, c). Each split vertex takes the
        // OPPOSING side's exact position, so the weld pairs them; the
        // displacement off the original chord is bounded by the
        // matching band's sag contract. Winding follows the original
        // a->b->c.
        let kc = 3 - ka - kb;
        let ic = tri[kc];
        let face = mesh.face_map[ti];
        let nb_normal = mesh.vertices[ib].normal;
        let mut prev = ia;
        for (_, _, pv) in &fan {
            let vi = mesh.vertices.len();
            let mut nv = mesh.vertices[ia].clone();
            nv.position = *pv;
            nv.normal = popkern_core::geometry::Vec3::new(
                0.5 * (nv.normal.x + nb_normal.x),
                0.5 * (nv.normal.y + nb_normal.y),
                0.5 * (nv.normal.z + nb_normal.z),
            );
            mesh.vertices.push(nv);
            if prev == ia {
                mesh.triangles[ti] = [ia, vi, ic];
            } else {
                new_triangles.push(([prev, vi, ic], face));
            }
            prev = vi;
        }
        new_triangles.push(([prev, ib, ic], face));
        split_tris.insert(ti);
        splits += fan.len();
    }
    for (tri, face) in new_triangles {
        mesh.triangles.push(tri);
        mesh.face_map.push(face);
    }
    if std::env::var_os("PK_ZIP_TRACE").is_some() {
        eprintln!(
            "  [zip] open={} off_curve={} no_candidate={} splits={splits}",
            stats.0, stats.1, stats.2
        );
    }
    splits
}

fn trace_skip(reason: &str, pa: &Point3, pb: &Point3) {
    if std::env::var_os("PK_ZIP_TRACE").is_some() {
        eprintln!(
            "  [zip]   {reason}: ({:.4},{:.4},{:.4})-({:.4},{:.4},{:.4})",
            pa.x, pa.y, pa.z, pb.x, pb.y, pb.z
        );
    }
}

fn project_point_segment(p: &Point3, a: &Point3, b: &Point3) -> (f64, Point3) {
    let ab = *b - *a;
    let len2 = ab.dot(&ab);
    if len2 < 1e-300 {
        return (0.0, *a);
    }
    let f = ((*p - *a).dot(&ab) / len2).clamp(0.0, 1.0);
    (
        f,
        Point3::new(a.x + ab.x * f, a.y + ab.y * f, a.z + ab.z * f),
    )
}
