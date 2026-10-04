//! Facet-deviation refinement for boundary-first face meshes.
//!
//! The interior lattice is sized by iso-line chord sag, which bounds
//! neither a lattice cell's diagonal facet nor the triangles that
//! bridge a region the lattice cannot reach (a 3.5 mm tab on an
//! r 37.5 cylinder wall spanned side to side, facets 0.036 off at a
//! 0.005 chord height; PK-651). The chord height bounds the FACET:
//! every triangle whose centroid lies farther than the chord height
//! from the surface gets a Steiner point (the midpoint of its most
//! sagging interior edge, else its chart centroid), inserted into the
//! constrained Delaunay triangulation incrementally
//! (`boundary_first_cdt::ConstrainedMesh`). Ring segments are
//! constraints, so a Steiner point beside a boundary segment never
//! forces it to split: the boundary samples, shared with the
//! neighbouring faces, are never touched and the mesh stays
//! watertight by construction.
//!
//! Transactional per round: a round that folds the face onto itself is
//! discarded and the previous round's mesh is kept.

use popkern_core::geometry::Point3;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};

use super::boundary_first_cdt::ConstrainedMesh;
use super::boundary_first_ring_index::RingIndex;

/// Refinement rounds; each inserts at most one Steiner point per
/// over-deviation triangle.
const ROUNDS: usize = 8;

/// Interior points one face may gain from refinement (a guard against
/// a surface the chord height cannot resolve at any density).
const POINT_BUDGET: usize = 200_000;

/// The constrained, region-filtered triangulation of one face chart.
pub(super) struct Chart<'a> {
    /// The face (diagnostics only).
    pub face: popkern_core::topology::FaceHandle,
    /// The face's surface.
    pub surface: &'a SurfaceGeometry,
    /// Metric scale of the chart.
    pub scale: (f64, f64),
    /// The chart's Delaunay scale (the metric, anisotropically
    /// compressed on ruled lattices).
    pub tri_scale: (f64, f64),
    /// Region, clearance and containment queries.
    pub ring_index: &'a RingIndex,
    /// Each ring's index range in the point list.
    pub ring_ranges: &'a [(usize, usize)],
}

/// Refine `triangles` over `uv_pts` (with exact positions `pos_pts`
/// where known) until every facet centroid lies within `chord_height`
/// of the surface, the rounds run out, or the point budget is spent.
pub(super) fn refine_by_deviation(
    chart: &Chart<'_>,
    uv_pts: &mut Vec<(f64, f64)>,
    pos_pts: &mut Vec<Option<Point3>>,
    triangles: &mut Vec<[usize; 3]>,
    chord_height: f64,
) {
    if matches!(chart.surface, SurfaceGeometry::Plane { .. }) || chord_height <= 0.0 {
        return;
    }
    // Refinement inserts into the constrained triangulation, so every
    // ring segment must already be an edge (a face that lost one goes
    // to the caller's coverage split instead).
    if !covers_rings(triangles, chart.ring_ranges) {
        return;
    }
    let segments: Vec<(usize, usize)> = chart
        .ring_ranges
        .iter()
        .flat_map(|&(start, end)| {
            (start..end).map(move |a| (a, start + (a + 1 - start) % (end - start)))
        })
        .collect();
    let (ku, kv) = chart.tri_scale;
    let scaled: Vec<(f64, f64)> = uv_pts.iter().map(|&(u, v)| (u * ku, v * kv)).collect();
    let Some(mut cdt) = ConstrainedMesh::new(&scaled, triangles, &segments) else {
        return;
    };
    let entry = uv_pts.len();
    // Only a fold refinement INTRODUCES discards a round.
    let folded_before = folds(triangles, uv_pts, pos_pts, chart);
    for round in 0..ROUNDS {
        let added = steiner_candidates(chart, uv_pts, pos_pts, triangles, chord_height);
        if added.is_empty() || uv_pts.len() + added.len() > entry + POINT_BUDGET {
            break;
        }
        let kept = uv_pts.len();
        for q in &added {
            if cdt.insert((q.0 * ku, q.1 * kv)).is_some() {
                uv_pts.push(*q);
                pos_pts.push(None);
            }
        }
        if uv_pts.len() == kept {
            break;
        }
        let refined = cdt.triangles();
        if !folded_before && folds(&refined, uv_pts, pos_pts, chart) {
            uv_pts.truncate(kept);
            pos_pts.truncate(kept);
            if popkern_core::trace::trace_enabled("boundary-first") {
                eprintln!(
                    "    [bf] face {:?} deviation refinement round {round} discarded (fold)",
                    chart.face
                );
            }
            break;
        }
        *triangles = refined;
        if popkern_core::trace::trace_enabled("boundary-first") {
            eprintln!(
                "    [bf] face {:?} deviation refinement round {round}: {} Steiner point(s)",
                chart.face,
                uv_pts.len() - kept
            );
        }
    }
}

/// Chart centroids of the triangles whose facet strays more than
/// `chord_height` from the surface, one per neighbourhood.
fn steiner_candidates(
    chart: &Chart<'_>,
    uv_pts: &[(f64, f64)],
    pos_pts: &[Option<Point3>],
    triangles: &[[usize; 3]],
    chord_height: f64,
) -> Vec<(f64, f64)> {
    let (su, sv) = chart.scale;
    let point =
        |i: usize| pos_pts[i].unwrap_or_else(|| chart.surface.point_at(uv_pts[i].0, uv_pts[i].1));
    let mut added: Vec<(f64, f64)> = Vec::new();
    // Candidates of one round keep a quarter of their triangle's
    // shortest edge apart (bucketed by that spacing).
    let mut buckets: std::collections::HashMap<(i64, i64), Vec<usize>> =
        std::collections::HashMap::new();
    for t in triangles {
        let (a, b, c) = (point(t[0]), point(t[1]), point(t[2]));
        let centroid = Point3::new(
            (a.x + b.x + c.x) / 3.0,
            (a.y + b.y + c.y) / 3.0,
            (a.z + b.z + c.z) / 3.0,
        );
        let (ua, ub, uc) = (uv_pts[t[0]], uv_pts[t[1]], uv_pts[t[2]]);
        let cuv = ((ua.0 + ub.0 + uc.0) / 3.0, (ua.1 + ub.1 + uc.1) / 3.0);
        if normal_gap(chart.surface, cuv, centroid) <= chord_height {
            continue;
        }
        // The facet's worst chord: the interior edge whose midpoint
        // sags most. Bisecting it halves the curved extent per round
        // (a centroid shrinks it by a third and leaves the long edge
        // in place, slivering a narrow tab instead of resolving it).
        // Ring segments are shared with the neighbouring faces and
        // never split here.
        let mut worst = (chord_height, None);
        for (i, j) in [(t[0], t[1]), (t[1], t[2]), (t[2], t[0])] {
            if is_ring_segment(chart.ring_ranges, i, j) {
                continue;
            }
            let (p, q) = (uv_pts[i], uv_pts[j]);
            let m = (0.5 * (p.0 + q.0), 0.5 * (p.1 + q.1));
            let (pi, pj) = (point(i), point(j));
            let chord_mid = Point3::new(
                0.5 * (pi.x + pj.x),
                0.5 * (pi.y + pj.y),
                0.5 * (pi.z + pj.z),
            );
            let sag = normal_gap(chart.surface, m, chord_mid);
            if sag > worst.0 {
                worst = (sag, Some(m));
            }
        }
        let usable =
            |p: (f64, f64)| chart.ring_index.inside(p) && !chart.ring_index.within_clearance(p);
        let cuv = match worst.1 {
            Some(m) if usable(m) => m,
            _ if usable(cuv) => cuv,
            _ => continue,
        };
        let spacing = 0.25
            * [(ua, ub), (ub, uc), (uc, ua)]
                .iter()
                .map(|(p, q)| (su * (p.0 - q.0)).hypot(sv * (p.1 - q.1)))
                .fold(f64::INFINITY, f64::min);
        if spacing.is_nan() || spacing <= 0.0 {
            continue;
        }
        let m = (cuv.0 * su, cuv.1 * sv);
        let cell = |x: f64| (x / spacing).floor() as i64;
        let (cx, cy) = (cell(m.0), cell(m.1));
        let crowded = (cx - 1..=cx + 1).any(|i| {
            (cy - 1..=cy + 1).any(|j| {
                buckets.get(&(i, j)).is_some_and(|ks| {
                    ks.iter().any(|&k| {
                        let q = added[k];
                        (su * q.0 - m.0).hypot(sv * q.1 - m.1) < spacing
                    })
                })
            })
        });
        if crowded {
            continue;
        }
        buckets.entry((cx, cy)).or_default().push(added.len());
        added.push(cuv);
    }
    // Insert in Morton order over the metric chart: Bowyer–Watson walks
    // from the last insertion, and candidates in old-triangle order
    // scattered those walks across the chart (23k points: 1.4 s per
    // re-triangulation instead of tens of milliseconds; PK-651).
    let (mut x0, mut x1, mut y0, mut y1) = (
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::INFINITY,
        f64::NEG_INFINITY,
    );
    for &(u, v) in &added {
        x0 = x0.min(u * su);
        x1 = x1.max(u * su);
        y0 = y0.min(v * sv);
        y1 = y1.max(v * sv);
    }
    let span = (x1 - x0).max(y1 - y0).max(1e-300);
    let morton = |&(u, v): &(f64, f64)| -> u64 {
        let q = |t: f64| ((t / span).clamp(0.0, 1.0) * 65535.0) as u64;
        let (mut x, mut y) = (q(u * su - x0), q(v * sv - y0));
        let mut code = 0u64;
        for bit in 0..16 {
            code |= ((x & 1) << (2 * bit)) | ((y & 1) << (2 * bit + 1));
            x >>= 1;
            y >>= 1;
        }
        code
    };
    added.sort_by_key(morton);
    added
}

/// Whether every ring segment is a triangle edge.
fn covers_rings(triangles: &[[usize; 3]], ring_ranges: &[(usize, usize)]) -> bool {
    let edges: std::collections::HashSet<(usize, usize)> = triangles
        .iter()
        .flat_map(|t| (0..3).map(move |k| (t[k].min(t[(k + 1) % 3]), t[k].max(t[(k + 1) % 3]))))
        .collect();
    ring_ranges.iter().all(|&(start, end)| {
        (start..end).all(|a| {
            let b = start + (a + 1 - start) % (end - start);
            edges.contains(&(a.min(b), a.max(b)))
        })
    })
}

/// Whether some 3-D edge of the triangulation is walked the same way
/// by two of its triangles (positions welded at 1e-9 of the extent):
/// the face covers part of itself twice.
fn folds(
    triangles: &[[usize; 3]],
    uv_pts: &[(f64, f64)],
    pos_pts: &[Option<Point3>],
    chart: &Chart<'_>,
) -> bool {
    let positions: Vec<Point3> = (0..uv_pts.len())
        .map(|i| pos_pts[i].unwrap_or_else(|| chart.surface.point_at(uv_pts[i].0, uv_pts[i].1)))
        .collect();
    let extent = positions
        .iter()
        .fold(0.0_f64, |m, p| {
            m.max(p.x.abs()).max(p.y.abs()).max(p.z.abs())
        })
        .max(1e-12);
    let s = 1e9 / extent;
    let q = |p: &Point3| {
        (
            (p.x * s).round() as i64,
            (p.y * s).round() as i64,
            (p.z * s).round() as i64,
        )
    };
    let mut seen = std::collections::HashSet::new();
    for t in triangles {
        for (i, j) in [(t[0], t[1]), (t[1], t[2]), (t[2], t[0])] {
            let (a, b) = (q(&positions[i]), q(&positions[j]));
            if a != b && !seen.insert((a, b)) {
                return true;
            }
        }
    }
    false
}

/// Whether `(i, j)` are consecutive points of one ring.
fn is_ring_segment(ring_ranges: &[(usize, usize)], i: usize, j: usize) -> bool {
    ring_ranges.iter().any(|&(start, end)| {
        let inside = |k: usize| k >= start && k < end;
        inside(i) && inside(j) && {
            let next = |k: usize| start + (k + 1 - start) % (end - start);
            next(i) == j || next(j) == i
        }
    })
}

/// The distance of `p` from the surface along the normal at `uv`: the
/// facet error at a facet point whose chart image is `uv`. The plain
/// distance to `S(uv)` adds the chart's tangential distortion (a
/// non-affine helical NURBS parameterization), which is no facet error
/// and over-refined a 48-turn thread flank about tenfold (PK-651).
pub(super) fn normal_gap(surface: &SurfaceGeometry, uv: (f64, f64), p: Point3) -> f64 {
    let s = surface.point_at(uv.0, uv.1);
    let n = surface.normal_at(uv.0, uv.1);
    let len = n.length();
    if len > 1e-300 {
        ((p - s).dot(&n) / len).abs()
    } else {
        p.distance_to(&s)
    }
}
