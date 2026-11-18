//! Facet-deviation refinement for boundary-first face meshes.
//!
//! The interior lattice is sized by iso-line chord sag, which bounds
//! neither a lattice cell's diagonal facet nor the triangles that
//! bridge a region the lattice cannot reach (a 3.5 mm tab on an
//! r 37.5 cylinder wall spanned side to side, facets 0.036 off at a
//! 0.005 chord height; PK-651). The chord height bounds the FACET:
//! every triangle whose centroid or interior chord midpoint lies
//! farther than the chord height from the surface gets a Steiner point (the midpoint of its most
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
use popkern_core::geometry::surface::SurfaceEval;
use popkern_math::PortableMath;

use super::boundary_first_cdt::ConstrainedMesh;
use super::boundary_first_ring_index::RingIndex;

/// Refinement rounds; each inserts at most one Steiner point per
/// over-deviation triangle.
const ROUNDS: usize = 8;

/// Share of the chord height under which a facet's largest edge sag
/// bounds its interior within the height (a quadratic surface's facet
/// strays at most 4/3 of its largest edge sag).
const INTERIOR_BOUND: f64 = 0.75;

/// Each edge's midpoint sag and chart midpoint, by its ends' indices.
type EdgeSags = popkern_math::hash::FastHashMap<(usize, usize), (f64, (f64, f64))>;

/// Interior points one face may gain from refinement (a guard against
/// a surface the chord height cannot resolve at any density).
const POINT_BUDGET: usize = 200_000;

/// The constrained, region-filtered triangulation of one face chart.
pub(super) struct Chart<'a> {
    /// The face (diagnostics only).
    pub face: popkern_core::topology::FaceHandle,
    /// The face's surface (through a regular chart for an atlas cap).
    pub surface: &'a dyn SurfaceEval,
    /// Whether the surface is a plane (no deviation to refine).
    pub planar: bool,
    /// Metric scale of the chart.
    pub scale: (f64, f64),
    /// The chart's Delaunay scale (the metric, anisotropically
    /// compressed on ruled lattices).
    pub tri_scale: (f64, f64, f64),
    /// Region, clearance and containment queries.
    pub ring_index: &'a RingIndex,
    /// Each ring's index range in the point list.
    pub ring_ranges: &'a [(usize, usize)],
    /// The point each point stands for in the triangulation (coincident
    /// ring points share one; `boundary_first_triangulate`).
    pub alias: &'a [usize],
}

/// Refine `triangles` over `uv_pts` (with exact positions `pos_pts`
/// where known) until every facet centroid and interior chord midpoint
/// lies within `chord_height` of the surface, the rounds run out, or the point budget is spent.
pub(super) fn refine_by_deviation(
    chart: &Chart<'_>,
    uv_pts: &mut Vec<(f64, f64)>,
    pos_pts: &mut Vec<Option<Point3>>,
    triangles: &mut Vec<[usize; 3]>,
    chord_height: f64,
) {
    if chart.planar || chord_height <= 0.0 {
        return;
    }
    // Refinement inserts into the constrained triangulation, so every
    // ring segment must already be an edge (a face that lost one goes
    // to the caller's coverage split instead).
    // Ring segments between the points the triangulation kept: a point
    // coincident with another (a seam's out-and-back, a ring sample
    // repeated) stands for its alias, and a segment between two of one
    // group is no segment. Read by index alone, a cap ring holding a
    // coincident pair never counted as covered and its triangles went
    // unrefined (a radius-40 ball's caps, 8.5e-3 off at a 5e-3 chord
    // height; PK-1274).
    let segments = ring_segments(chart);
    if !covers_rings(triangles, &segments) {
        return;
    }
    let segment_set: popkern_math::hash::FastHashSet<(usize, usize)> = segments
        .iter()
        .map(|&(a, b)| (a.min(b), a.max(b)))
        .collect();
    let (ku, kuv, kv) = chart.tri_scale;
    let scaled: Vec<(f64, f64)> = uv_pts
        .iter()
        .map(|&(u, v)| (u * ku + v * kuv, v * kv))
        .collect();
    let Some(mut cdt) = ConstrainedMesh::new(&scaled, triangles, &segments) else {
        return;
    };
    let entry = uv_pts.len();
    let mut sags = EdgeSags::default();
    // Only a fold refinement INTRODUCES discards a round.
    let folded_before = folds(triangles, uv_pts, pos_pts, chart);
    for round in 0..ROUNDS {
        let added = steiner_candidates(
            chart,
            (uv_pts, pos_pts),
            triangles,
            (&segment_set, &mut sags),
            chord_height,
        );
        if added.is_empty() || uv_pts.len() + added.len() > entry + POINT_BUDGET {
            break;
        }
        let kept = uv_pts.len();
        for q in &added {
            if cdt.insert((q.0 * ku + q.1 * kuv, q.1 * kv)).is_some() {
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

/// Steiner points of the triangles whose facet strays more than
/// `chord_height` from the surface (at the centroid or an interior
/// chord's midpoint), one per neighbourhood.
fn steiner_candidates(
    chart: &Chart<'_>,
    (uv_pts, pos_pts): (&[(f64, f64)], &[Option<Point3>]),
    triangles: &[[usize; 3]],
    (segments, sags): (
        &popkern_math::hash::FastHashSet<(usize, usize)>,
        &mut EdgeSags,
    ),
    chord_height: f64,
) -> Vec<(f64, f64)> {
    let (su, sv) = chart.scale;
    let point =
        |i: usize| pos_pts[i].unwrap_or_else(|| chart.surface.point_at(uv_pts[i].0, uv_pts[i].1));
    let mut added: Vec<(f64, f64)> = Vec::new();
    // Candidates of one round keep a quarter of their triangle's
    // shortest edge apart (bucketed by that spacing).
    let mut buckets: popkern_math::hash::FastHashMap<(i64, i64), Vec<usize>> =
        popkern_math::hash::FastHashMap::default();
    for t in triangles {
        let (a, b, c) = (point(t[0]), point(t[1]), point(t[2]));
        let (ua, ub, uc) = (uv_pts[t[0]], uv_pts[t[1]], uv_pts[t[2]]);
        // Each edge's sag at its midpoint, measured once (two triangles
        // share an interior edge, and edges outlive a round).
        let mut sag_of = |i: usize, j: usize| -> (f64, (f64, f64)) {
            *sags.entry((i.min(j), i.max(j))).or_insert_with(|| {
                let (p, q) = (uv_pts[i], uv_pts[j]);
                let m = (0.5 * (p.0 + q.0), 0.5 * (p.1 + q.1));
                let (pi, pj) = (point(i), point(j));
                let chord_mid = Point3::new(
                    0.5 * (pi.x + pj.x),
                    0.5 * (pi.y + pj.y),
                    0.5 * (pi.z + pj.z),
                );
                (normal_gap(chart.surface, m, chord_mid), m)
            })
        };
        // The facet's worst chord: the interior edge whose midpoint
        // sags most. Bisecting it halves the curved extent per round
        // (a centroid shrinks it by a third and leaves the long edge
        // in place, slivering a narrow tab instead of resolving it).
        // Ring segments are shared with the neighbouring faces and
        // never split here. A facet strays past the chord height where
        // an interior chord sags past it (a right triangle's farthest
        // point from the surface is its long edge's midpoint: a ball's
        // band lattice diagonals sagged 5.6e-3 at a radius of 10 and
        // 8.5e-3 at 40 against a 5e-3 chord height while the centroids
        // read within it; PK-1274) or where its interior does.
        let mut worst = (chord_height, None);
        let mut largest = 0.0_f64;
        for (i, j) in [(t[0], t[1]), (t[1], t[2]), (t[2], t[0])] {
            let (sag, m) = sag_of(i, j);
            largest = largest.max(sag);
            if sag > worst.0 && !segments.contains(&(i.min(j), i.max(j))) {
                worst = (sag, Some(m));
            }
        }
        // Over a locally quadratic surface a facet strays at most 4/3 of
        // its largest edge sag (an equilateral facet on a sphere, its
        // circumcentre): within 3/4 of the height there, its interior is
        // within the height, and only the rest are evaluated inside (the
        // centroid and an acute facet's circumcentre cost a chart
        // inversion each on an atlas cap; read for every facet, the
        // paraboloid cup's measure ran 1.7x slower).
        if worst.1.is_none() && largest <= INTERIOR_BOUND * chord_height {
            continue;
        }
        let centroid = Point3::new(
            (a.x + b.x + c.x) / 3.0,
            (a.y + b.y + c.y) / 3.0,
            (a.z + b.z + c.z) / 3.0,
        );
        let cuv = ((ua.0 + ub.0 + uc.0) / 3.0, (ua.1 + ub.1 + uc.1) / 3.0);
        // The facet's interior stray: at its centroid, and at its
        // circumcentre when that lies inside (an acute facet's farthest
        // point from a sphere: a radius-40 ball's band read 5.1e-3 there
        // with its centroid and chords within 5e-3; PK-1274).
        let interior_gap = match circumcentre_weights(a, b, c) {
            Some(w) => {
                let at = |p: f64, q: f64, r: f64| w[0] * p + w[1] * q + w[2] * r;
                let p = Point3::new(at(a.x, b.x, c.x), at(a.y, b.y, c.y), at(a.z, b.z, c.z));
                let quv = (at(ua.0, ub.0, uc.0), at(ua.1, ub.1, uc.1));
                normal_gap(chart.surface, cuv, centroid).max(normal_gap(chart.surface, quv, p))
            }
            None => normal_gap(chart.surface, cuv, centroid),
        };
        if interior_gap <= chord_height && worst.1.is_none() {
            continue;
        }
        let usable =
            |p: (f64, f64)| chart.ring_index.inside(p) && !chart.ring_index.within_clearance(p);
        // The clearance band keeps lattice points from slivering against
        // the boundary, but a facet past the chord height whose chords all
        // sag within it bulges in its interior: its centroid may sit in
        // the band, never outside the trim. PK-569: two coil facets at the
        // seam by the ground cut sagged 5.06e-3 and 5.29e-3 against the
        // 5e-3 chord with their centroids in the band.
        let cuv = match worst.1 {
            Some(m) if usable(m) => m,
            _ if usable(cuv) => cuv,
            // A chord past the height with its midpoint and the centroid
            // both in the band: the chord's midpoint, inside the trim (a
            // radius-40 ball's cap triangle beside its row sagged 8.5e-3
            // there; PK-1274).
            Some(m) if chart.ring_index.inside(m) => m,
            // No chord sags past the height: the bulge is inside the
            // facet, and only its centroid cures it.
            None if chart.ring_index.inside(cuv) => cuv,
            _ => continue,
        };
        let spacing = 0.25
            * [(ua, ub), (ub, uc), (uc, ua)]
                .iter()
                .map(|(p, q)| (su * (p.0 - q.0)).phypot(sv * (p.1 - q.1)))
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
                        (su * q.0 - m.0).phypot(sv * q.1 - m.1) < spacing
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
fn covers_rings(triangles: &[[usize; 3]], segments: &[(usize, usize)]) -> bool {
    let edges: popkern_math::hash::FastHashSet<(usize, usize)> = triangles
        .iter()
        .flat_map(|t| (0..3).map(move |k| (t[k].min(t[(k + 1) % 3]), t[k].max(t[(k + 1) % 3]))))
        .collect();
    segments
        .iter()
        .all(|&(a, b)| edges.contains(&(a.min(b), a.max(b))))
}

/// The barycentric weights of the circumcentre of the triangle `a b c`
/// when it lies strictly inside (an acute triangle), else `None`.
fn circumcentre_weights(a: Point3, b: Point3, c: Point3) -> Option<[f64; 3]> {
    let (ab, ac) = (b - a, c - a);
    let n = ab.cross(&ac);
    let nn = n.dot(&n);
    if nn <= 0.0 || !nn.is_finite() {
        return None;
    }
    // a + s·ab + t·ac, from (p − a)·ab = |ab|²/2 and (p − a)·ac = |ac|²/2;
    // the system's determinant is |ab × ac|².
    let (bb, cc, d) = (ab.dot(&ab), ac.dot(&ac), ab.dot(&ac));
    let s = cc * (bb - d) / (2.0 * nn);
    let t = bb * (cc - d) / (2.0 * nn);
    let w = [1.0 - s - t, s, t];
    w.iter().all(|&x| x > 0.0).then_some(w)
}

/// Each ring segment between the points standing for its ends (their
/// aliases), segments within one coincident group left out.
fn ring_segments(chart: &Chart<'_>) -> Vec<(usize, usize)> {
    let alias = |k: usize| chart.alias.get(k).copied().unwrap_or(k);
    chart
        .ring_ranges
        .iter()
        .flat_map(|&(start, end)| {
            (start..end).map(move |a| (alias(a), alias(start + (a + 1 - start) % (end - start))))
        })
        .filter(|&(a, b)| a != b)
        .collect()
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
    let mut seen = popkern_math::hash::FastHashSet::default();
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

/// The distance of `p` from the surface along the normal at `uv`: the
/// facet error at a facet point whose chart image is `uv`. The plain
/// distance to `S(uv)` adds the chart's tangential distortion (a
/// non-affine helical NURBS parameterization), which is no facet error
/// and over-refined a 48-turn thread flank about tenfold (PK-651).
pub(super) fn normal_gap(surface: &dyn SurfaceEval, uv: (f64, f64), p: Point3) -> f64 {
    let s = surface.point_at(uv.0, uv.1);
    let n = surface.normal_at(uv.0, uv.1);
    let len = n.length();
    if len > 1e-300 {
        ((p - s).dot(&n) / len).abs()
    } else {
        p.distance_to(&s)
    }
}
