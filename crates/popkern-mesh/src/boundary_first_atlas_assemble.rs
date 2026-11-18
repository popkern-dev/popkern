//! Assembling the atlas's pieces from a face's rings (chart-service.md
//! §6a, step 4): arcs between row samples, joined along each cap row in
//! the sense the face flag gives.

use popkern_core::geometry::Point3;
use popkern_core::geometry::chart::atlas::{along_offset, canonical_along};
use popkern_core::geometry::chart::ops::evaluate;
use popkern_core::geometry::chart::{Axis, ChartDescriptor};
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::topology::EdgeHandle;
use popkern_math::PortableMath;

use super::super::BoundaryRing;
use super::{Cut, Piece, Side};

/// The pieces' rings, in the face chart (the cap pieces' converted to
/// their projection by the mesher): `(piece, rings)`.
pub(super) type Pieces = Vec<(Piece, Vec<BoundaryRing>)>;

/// Assemble the pieces from the rings (see the module docs). `None` when a
/// walk along a row finds no arc entering its piece (a configuration the
/// face flag contradicts).
pub(super) fn assemble(
    cuts: &[Cut],
    (built, ring_sides): (&[BoundaryRing], &[Vec<Side>]),
    (orientation, chord, angular, steps): (bool, f64, Option<f64>, &[f64]),
    surface: (&SurfaceGeometry, &ChartDescriptor),
    row_in_face: &dyn Fn(usize) -> bool,
) -> Option<Pieces> {
    let mut arcs: Vec<Arc> = Vec::new();
    let mut whole: Vec<(Piece, BoundaryRing)> = Vec::new();
    for (r, ring) in built.iter().enumerate() {
        let sides = &ring_sides[r];
        let rows: Vec<usize> = (0..sides.len())
            .filter(|&i| matches!(sides[i], Side::Row(_)))
            .collect();
        if rows.is_empty() {
            let piece = match sides.first() {
                Some(Side::Cap(k)) => Piece::Cap(*k),
                _ => Piece::Band,
            };
            whole.push((piece, ring.clone()));
            continue;
        }
        for (m, &start) in rows.iter().enumerate() {
            let end = rows[(m + 1) % rows.len()];
            let len = (end + sides.len() - start) % sides.len();
            let len = if len == 0 { sides.len() } else { len };
            let idx: Vec<usize> = (0..=len).map(|q| (start + q) % sides.len()).collect();
            let interior = idx[1..idx.len() - 1].iter().map(|&i| sides[i]).next();
            let row_k = |i: usize| match sides[i] {
                Side::Row(k) => k,
                _ => usize::MAX,
            };
            let piece = match interior {
                Some(Side::Cap(k)) => Piece::Cap(k),
                Some(_) => Piece::Band,
                // A segment lying along the row: the cap's side.
                None => Piece::Cap(row_k(start)),
            };
            arcs.push(Arc {
                ring: r,
                idx,
                piece,
                from_row: row_k(start),
                to_row: row_k(end),
            });
        }
    }
    let mut out: Pieces = Vec::new();
    let mut add = |piece: Piece, ring: BoundaryRing| match out.iter_mut().find(|(p, _)| *p == piece)
    {
        Some((_, rings)) => rings.push(ring),
        None => out.push((piece, vec![ring])),
    };
    for (piece, ring) in whole {
        add(piece, ring);
    }
    // Walk each piece's arcs into rings.
    let pieces: Vec<Piece> = (0..cuts.len())
        .map(Piece::Cap)
        .chain(std::iter::once(Piece::Band))
        .collect();
    for &piece in &pieces {
        let mut used = vec![false; arcs.len()];
        for first in 0..arcs.len() {
            if used[first] || arcs[first].piece != piece {
                continue;
            }
            let mut ring = BoundaryRing {
                uv: Vec::new(),
                pos: Vec::new(),
                sources: Vec::new(),
            };
            let mut cur = first;
            for _ in 0..=arcs.len() {
                used[cur] = true;
                let arc = &arcs[cur];
                let src = &built[arc.ring];
                for &i in &arc.idx[..arc.idx.len() - 1] {
                    ring.uv.push(src.uv[i]);
                    ring.pos.push(src.pos[i]);
                    ring.sources.push(src.sources[i]);
                }
                let end_sample = *arc.idx.last()?;
                let end = *arc.idx.last()?;
                let k = arc.to_row;
                let cut = cuts.get(k)?;
                // The arc after it on its own ring.
                let next_on_ring = arcs
                    .iter()
                    .position(|a| a.ring == arc.ring && a.idx[0] == end)?;
                let next = if arcs[next_on_ring].piece == piece {
                    next_on_ring
                } else {
                    // The crossing itself, then the row to the next.
                    ring.uv.push(src.uv[end_sample]);
                    ring.pos.push(src.pos[end_sample]);
                    ring.sources.push(src.sources[end_sample]);
                    let walk = walk_sense(cut, piece, orientation);
                    let at = cut.along(src.uv[end]);
                    let target = next_entry(&arcs, built, (k, cut), piece, at, walk)?;
                    let t_uv = built[arcs[target].ring].uv[arcs[target].idx[0]];
                    // A segment two pieces share is generated from the same
                    // crossing (the smaller ring/sample key) in its own
                    // sense, so both read the same points to the bit; the
                    // other piece takes them reversed.
                    let (here, there) = ((arc.ring, end), (arcs[target].ring, arcs[target].idx[0]));
                    let samples = if here <= there {
                        row_samples(
                            cut,
                            surface,
                            at,
                            cut.along(t_uv),
                            walk,
                            (chord, angular, steps[k]),
                        )
                    } else {
                        let mut s = row_samples(
                            cut,
                            surface,
                            cut.along(t_uv),
                            at,
                            -walk,
                            (chord, angular, steps[k]),
                        );
                        s.reverse();
                        s
                    };
                    for (uv, p) in samples {
                        ring.uv.push(uv);
                        ring.pos.push(p);
                        ring.sources.push(row_source());
                    }
                    target
                };
                if next == first {
                    break;
                }
                if used[next] {
                    return None;
                }
                cur = next;
            }
            add(piece, ring);
        }
    }
    // Rows no ring crosses: wholly in the face (both pieces bounded by
    // them) or out of it.
    for (k, cut) in cuts.iter().enumerate() {
        let crossed = arcs.iter().any(|a| a.from_row == k || a.to_row == k);
        if crossed || !row_in_face(k) {
            continue;
        }
        for piece in [Piece::Cap(k), Piece::Band] {
            let walk = walk_sense(cut, piece, orientation);
            let period = cut.period?;
            let lo = cut.start;
            let (start, end) = if walk > 0.0 {
                (lo, lo + period)
            } else {
                (lo + period, lo)
            };
            // The start's position read at the row's own start for both
            // pieces (the band walks from a period on): the pieces weld
            // there by position bits.
            let first_uv = cut.uv(start);
            let mut ring = BoundaryRing {
                uv: vec![first_uv],
                pos: vec![evaluate(surface.1, surface.0, cut.uv(lo))],
                sources: vec![row_source()],
            };
            for (uv, p) in row_samples(cut, surface, start, end, walk, (chord, angular, steps[k])) {
                ring.uv.push(uv);
                ring.pos.push(p);
                ring.sources.push(row_source());
            }
            add(piece, ring);
        }
    }
    Some(out)
}

/// The source recorded for a row sample: no edge (rows are the atlas's,
/// sampled once per face; a split request on one is dropped).
pub(super) fn row_source() -> (EdgeHandle, f64, f64) {
    (EdgeHandle::default(), 0.0, 0.0)
}

/// One arc of a ring between two row samples (inclusive).
struct Arc {
    ring: usize,
    idx: Vec<usize>,
    piece: Piece,
    from_row: usize,
    to_row: usize,
}

/// The direction along the row (`±1` in the along coordinate) in which
/// `piece`'s boundary continues after an arc leaves it: the face lies
/// left of its loops when the flag is set (right otherwise), so the walk
/// keeps the piece on that side.
fn walk_sense(cut: &Cut, piece: Piece, orientation: bool) -> f64 {
    // Leaving the cap moves toward the band (+toward_band in the row
    // coordinate); leaving the band moves away from it.
    let dc = match piece {
        Piece::Cap(_) => cut.toward_band,
        Piece::Band => -cut.toward_band,
    };
    // Left of a step (0, dc) in (along, coordinate) when the coordinate is
    // v is (−dc, 0); when it is u, the step (dc, 0) has (0, dc) on its left.
    let left = match cut.cap.singular.axis {
        Axis::V => -dc,
        Axis::U => dc,
    };
    if orientation { left } else { -left }
}

/// The next arc entering `piece` from row `k` in the walk direction from
/// `at`, by the along coordinate (circular when the row closes).
fn next_entry(
    arcs: &[Arc],
    built: &[BoundaryRing],
    (k, cut): (usize, &Cut),
    piece: Piece,
    at: f64,
    walk: f64,
) -> Option<usize> {
    let surface_period = cut.period;
    arcs.iter()
        .enumerate()
        .filter(|(_, a)| a.piece == piece && a.from_row == k)
        .map(|(i, a)| {
            let along = cut.along(built[a.ring].uv[a.idx[0]]);
            (i, walk_offset(surface_period, at, along, walk))
        })
        .filter(|&(_, d)| d > 0.0)
        .min_by(|x, y| x.1.total_cmp(&y.1))
        .map(|(i, _)| i)
}

/// Samples along row `cut` strictly between `from` and `to` (along
/// coordinate, `walk` its sign; the ends are the crossings' own samples),
/// split by halving until each chord's midpoint lies within `chord` of the
/// row in space.
fn row_samples(
    cut: &Cut,
    (surface, descriptor): (&SurfaceGeometry, &ChartDescriptor),
    from: f64,
    to: f64,
    walk: f64,
    (chord, angular, max_step): (f64, Option<f64>, f64),
) -> Vec<((f64, f64), Point3)> {
    let end = match cut.period {
        Some(_) => from + walk * walk_offset(cut.period, from, to, walk),
        None => to,
    };
    // Generated in the +along sense so both pieces sharing a segment get
    // the same points (one reversed).
    let (lo, hi) = if from <= end {
        (from, end)
    } else {
        (end, from)
    };
    // On the row's canonical branch: the two pieces walk a segment from
    // crossings whose parameters can lie a period apart, and must read
    // the same points to the bit to weld; the parameters go back to the
    // walk's branch for the piece's chart.
    let shift = lo - canonical_along(cut.period, cut.start, lo);
    let (lo, hi) = (lo - shift, hi - shift);
    // Canonical along values with their points; parameters on the walk's
    // branch are made from them only at the end (recovering the canonical
    // value from a shifted one is not exact, and the pieces would split
    // the row at points a rounding apart).
    let at = |a: f64| (a, evaluate(descriptor, surface, cut.uv(a)));
    // The row is sampled as edges are: by chord height and by the
    // surface normal's turn under the angular facet bound.
    let normal_at = |a: f64| {
        let uv = cut.uv(a);
        let (su, sv) = surface.first_partials(uv.0, uv.1);
        su.cross(&sv)
            .normalized()
            .unwrap_or(popkern_core::geometry::Vec3::new(0.0, 0.0, 1.0))
    };
    let mut canon = vec![at(lo), at(hi)];
    let mut i = 0;
    while i + 1 < canon.len() && canon.len() < 4096 {
        let (a, b) = (canon[i].0, canon[i + 1].0);
        let m = at(0.5 * (a + b));
        let chord_mid = Point3::new(
            0.5 * (canon[i].1.x + canon[i + 1].1.x),
            0.5 * (canon[i].1.y + canon[i + 1].1.y),
            0.5 * (canon[i].1.z + canon[i + 1].1.z),
        );
        let turn = |x: f64, y: f64| {
            let (n1, n2) = (normal_at(x), normal_at(y));
            n1.dot(&n2).clamp(-1.0, 1.0).pacos()
        };
        let over_angle = angular.is_some_and(|bound| turn(a, b) > bound);
        let too_long = canon[i].1.distance_to(&canon[i + 1].1) > max_step;
        if m.1.distance_to(&chord_mid) > chord || over_angle || too_long {
            canon.insert(i + 1, m);
        } else {
            i += 1;
        }
    }
    let mut pts: Vec<((f64, f64), Point3)> = canon
        .into_iter()
        .map(|(a, p)| (cut.uv(a + shift), p))
        .collect();
    if from > end {
        pts.reverse();
    }
    pts.remove(0);
    pts.pop();
    pts
}

/// How far the walk along a row runs from `from` to the crossing at `to`:
/// the offset along the row in the walk's sense, and a crossing at the
/// same place along a closed row a whole turn away (a seam crossing the
/// row: down one use to the pole and back up the other), never no
/// distance, which would leave the piece between with no width.
fn walk_offset(period: Option<f64>, from: f64, to: f64, walk: f64) -> f64 {
    let d = along_offset(period, from, to, walk);
    match period {
        Some(p) if d <= 1e-12 * p => d + p,
        _ => d,
    }
}
