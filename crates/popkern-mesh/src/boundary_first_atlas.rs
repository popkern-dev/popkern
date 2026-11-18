//! The atlas (chart-service.md §6a): a face reaching a chart singularity
//! is meshed in pieces, each in a chart with no singular row in it. Every
//! piece's rings are cut where they cross a cap row (the row at which the
//! normal has turned a quarter right angle from the singular point's,
//! `chart::atlas::cap_row`); a cap piece is meshed in the pole projection,
//! the band in the face chart. No ring ever runs along a pole row, so the
//! lat-long chart's pole machinery (which way a ring runs along a row
//! that is one point in space) has nothing to decide.
//!
//! Pieces are assembled chart-free of the pole: each ring is split into
//! arcs at its samples on a cap row; at a crossing, the piece the arc
//! leaves continues along the row to the next crossing whose arc enters
//! it, the side read from the face flag (the face lies left of its loops
//! in the face chart, and the pole projection keeps that orientation). A
//! row no ring crosses is wholly in the face or out of it, by the
//! kernel's point-in-face authority at a point of it.

use popkern_core::geometry::Point3;
use popkern_core::geometry::chart::atlas::{CAP_TURNS, CapRow, cap_row_at};
use popkern_core::geometry::chart::{Axis, ChartDescriptor, SingularKind};
use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle};
use popkern_math::PortableMath;

use super::*;

#[path = "boundary_first_atlas_assemble.rs"]
mod assemble;
#[path = "boundary_first_atlas_cap.rs"]
mod cap;
use assemble::{Pieces, assemble};

/// Where a ring sample lies.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Side {
    /// Inside cap `k`.
    Cap(usize),
    /// On cap `k`'s row.
    Row(usize),
    /// In the band (no cap).
    Band,
}

/// One piece's side, for assembly.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Piece {
    Cap(usize),
    Band,
}

/// A cap row with its geometry for side tests and row walks.
struct Cut {
    cap: CapRow,
    /// Sample within this of the row's value is on it.
    eps: f64,
    /// `+1` when the band lies at larger values of the row coordinate.
    toward_band: f64,
    /// The period along the row (when the chart closes along it) and
    /// the along coordinate's domain start.
    period: Option<f64>,
    start: f64,
}

impl Cut {
    fn coordinate(&self, uv: (f64, f64)) -> f64 {
        match self.cap.singular.axis {
            Axis::V => uv.1,
            Axis::U => uv.0,
        }
    }
    /// The coordinate along the row (the collapsed one at the pole).
    fn along(&self, uv: (f64, f64)) -> f64 {
        match self.cap.singular.axis {
            Axis::V => uv.0,
            Axis::U => uv.1,
        }
    }
    fn uv(&self, along: f64) -> (f64, f64) {
        match self.cap.singular.axis {
            Axis::V => (along, self.cap.at),
            Axis::U => (self.cap.at, along),
        }
    }
    fn in_cap(&self, uv: (f64, f64)) -> bool {
        (self.coordinate(uv) - self.cap.at) * self.toward_band < -self.eps
    }
    fn on_row(&self, uv: (f64, f64)) -> bool {
        (self.coordinate(uv) - self.cap.at).abs() <= self.eps
    }
}

/// The side of a ring sample, read from its position in space (its
/// surface parameters by projection), as the crossing bisection reads
/// the edge: a pcurve's parameters can sit off the row by their fit
/// where the canonical sample lies on it.
fn side_at(cuts: &[Cut], surface: &SurfaceGeometry, uv: (f64, f64), p: &Point3) -> Side {
    let uv = surface.closest_point(p).map_or(uv, |q| {
        // The projection's along coordinate is arbitrary; only the row
        // coordinate is read.
        match cuts.first().map(|c| c.cap.singular.axis) {
            Some(Axis::U) => (q.0, uv.1),
            _ => (uv.0, q.1),
        }
    });
    side(cuts, uv)
}

fn side(cuts: &[Cut], uv: (f64, f64)) -> Side {
    for (k, cut) in cuts.iter().enumerate() {
        if cut.on_row(uv) {
            return Side::Row(k);
        }
        if cut.in_cap(uv) {
            return Side::Cap(k);
        }
    }
    Side::Band
}

/// Tessellate `face` on the atlas, or `None` when it reaches no chart
/// singularity (the caller's path meshes it). `Some(None)`: the atlas
/// owns the face but declined this round (it requested canonical samples
/// at its row crossings, or a piece declined).
pub(super) fn tessellate_on_atlas(
    brep: &BRep,
    face_handle: FaceHandle,
    (surface, orientation): (&SurfaceGeometry, bool),
    built: &[BoundaryRing],
    (descriptor, options): (&ChartDescriptor, &TessellationOptions),
    split_requests: &mut Vec<SplitRequest>,
) -> Option<Option<FaceMesh>> {
    let cut_at = |cap: CapRow| {
        let ((u0, _), (v0, _)) = descriptor.domain;
        let (along_axis, lo) = match cap.singular.axis {
            Axis::V => (Axis::U, u0),
            Axis::U => (Axis::V, v0),
        };
        Cut {
            eps: 1e-9 * (cap.at - cap.singular.at).abs(),
            toward_band: (cap.at - cap.singular.at).signum(),
            period: descriptor.period(along_axis).map(|p| p.period),
            start: if lo.is_finite() { lo } else { 0.0 },
            cap,
        }
    };
    // Each singular row's cap row at the candidate turn lying farthest
    // from every ring's turning points in the row coordinate (a loop
    // running along the row flips its crossings as its samples refine:
    // the dome's cut circle grazed the 45° row and crossed it twice in one
    // round, not at all in the next).
    let cuts: Vec<Cut> = descriptor
        .singular_rows
        .iter()
        .filter(|row| row.kind == SingularKind::Chart)
        .filter_map(|row| {
            let candidates: Vec<(Cut, f64)> = CAP_TURNS
                .iter()
                .filter_map(|&turn| cap_row_at(descriptor, surface, row, turn).ok())
                .map(|cap| {
                    let cut = cut_at(cap);
                    let clearance = turning_clearance(&cut, surface, built);
                    (cut, clearance)
                })
                .collect();
            // The first, in preference order, clear of every turning point
            // by a twentieth of its cap's span (a margin far above the
            // samples' resolution, so the choice holds as they refine);
            // the clearest when none is.
            let clear =
                |(cut, c): &(Cut, f64)| *c > 0.05 * (cut.cap.at - cut.cap.singular.at).abs();
            let pick = candidates.iter().position(clear).or_else(|| {
                (0..candidates.len()).max_by(|&a, &b| candidates[a].1.total_cmp(&candidates[b].1))
            })?;
            candidates.into_iter().nth(pick).map(|(cut, _)| cut)
        })
        .collect();
    if cuts.is_empty() {
        return None;
    }
    // The whole surface within one cap (its normal never turns the cap's
    // quarter right angle): the face is one cap piece, no row to cut.
    if let Some(k) = cuts.iter().position(|c| c.cap.whole) {
        let pieces = vec![(Piece::Cap(k), built.to_vec())];
        return Some(cap::mesh_pieces(
            brep,
            face_handle,
            (surface, orientation),
            (descriptor, options),
            (&cuts, pieces),
            split_requests,
        ));
    }
    // A cap the face reaches: a ring sample in it, or its row in the face.
    // A row no ring crosses is wholly in the face or out of it, read by
    // parity in the band chart (no pole in it): from a point just on the
    // face side of a ring sample straight to the row, an even number of
    // ring crossings puts the row in the face. Not the trim classifier,
    // which answers Outside on a row inside a dome holding both poles
    // (PK-1273, until W2c moves it onto the chart service).
    let row_in_face = |cut: &Cut| row_in_face_by_parity(cut, surface, built, orientation);
    let reached: Vec<bool> = cuts
        .iter()
        .enumerate()
        .map(|(k, cut)| {
            built.iter().any(|r| {
                r.uv.iter()
                    .zip(&r.pos)
                    .any(|(&uv, p)| side_at(&cuts, surface, uv, p) == Side::Cap(k))
            }) || row_in_face(cut)
        })
        .collect();
    if !reached.iter().any(|&r| r) {
        return None;
    }
    let trace = popkern_core::trace::trace_enabled("atlas");
    // Crossings without a sample on the row: request one, decline.
    let Ok(sides) = ring_sides(brep, (surface, &cuts), built, split_requests) else {
        return Some(None);
    };
    let Some(sides) = sides else {
        if trace {
            eprintln!("    [atlas] face {face_handle:?}: crossing samples requested");
        }
        return Some(None);
    };
    if trace {
        for (r, sd) in sides.iter().enumerate() {
            let count = |f: &dyn Fn(&Side) -> bool| sd.iter().filter(|x| f(x)).count();
            eprintln!(
                "    [atlas] face {face_handle:?}: ring {r} of {}: cap {} row {} band {}",
                sd.len(),
                count(&|x| matches!(x, Side::Cap(_))),
                count(&|x| matches!(x, Side::Row(_))),
                count(&|x| matches!(x, Side::Band))
            );
        }
    }
    // Rows are sampled to their cap's lattice (a row segment has no edge
    // for the cap to split): at most one and a half estimated cells.
    let steps: Vec<f64> = cuts
        .iter()
        .map(|c| {
            cap::estimated_cap_cell(descriptor, surface, c, options)
                .map_or(f64::INFINITY, |c| 1.5 * c)
        })
        .collect();
    let pieces = assemble(
        &cuts,
        (built, &sides),
        (
            orientation,
            options.chord_height,
            options.angular_bound(),
            &steps,
        ),
        (surface, descriptor),
        &|k| row_in_face(&cuts[k]),
    );
    if trace {
        eprintln!(
            "    [atlas] face {face_handle:?}: caps {:?} reached {reached:?}, pieces {:?}",
            cuts.iter().map(|c| c.cap.at).collect::<Vec<_>>(),
            pieces.as_ref().map(|p| p
                .iter()
                .map(|(piece, rings)| (
                    *piece,
                    rings.iter().map(|r| r.uv.len()).collect::<Vec<_>>()
                ))
                .collect::<Vec<_>>())
        );
    }
    let pieces = pieces?;
    Some(cap::mesh_pieces(
        brep,
        face_handle,
        (surface, orientation),
        (descriptor, options),
        (&cuts, pieces),
        split_requests,
    ))
}

/// Every ring sample's side, with each crossing of a cap row resolved:
/// a crossing between two samples is requested as a canonical sample at
/// its curve fraction (bisected on the edge in space), or, within the
/// encroachment dedup of an existing sample (which the driver would not
/// insert), that sample snapped onto the row. `None` when a sample was
/// requested (the face declines this round).
fn ring_sides(
    brep: &BRep,
    (surface, cuts): (&SurfaceGeometry, &[Cut]),
    built: &[BoundaryRing],
    split_requests: &mut Vec<SplitRequest>,
) -> Result<Option<Vec<Vec<Side>>>, ()> {
    let mut requested = false;
    let mut all = Vec::with_capacity(built.len());
    for ring in built {
        let n = ring.uv.len();
        let mut sides: Vec<Side> = (0..n)
            .map(|i| side_at(cuts, surface, ring.uv[i], &ring.pos[i]))
            .collect();
        for i in 0..n {
            let j = (i + 1) % n;
            let crossing = match (sides[i], sides[j]) {
                (Side::Cap(k), x) | (x, Side::Cap(k)) if x != Side::Cap(k) && x != Side::Row(k) => {
                    Some(k)
                }
                _ => None,
            };
            let Some(k) = crossing else {
                continue;
            };
            let (edge, fa, fb) = ring.sources[i];
            let Some(f) = crossing_fraction(brep, surface, &cuts[k], edge, (fa, fb)) else {
                // A crossing the edge cannot locate: the face cannot be cut
                // there, so the atlas declines it (never a ring left uncut
                // across a row, which would mesh the wrong region).
                if popkern_core::trace::trace_enabled("atlas") {
                    eprintln!(
                        "    [atlas] crossing on {edge:?} in ({fa}, {fb}) not located: {:?} -> {:?}",
                        sides[i], sides[j]
                    );
                }
                return Err(());
            };
            let dedup = crate::encroachment::SAMPLE_DEDUP;
            if (f - fa).abs() < dedup {
                sides[i] = Side::Row(k);
            } else if (f - fb).abs() < dedup {
                sides[j] = Side::Row(k);
            } else {
                split_requests.push((edge, f));
                requested = true;
            }
        }
        all.push(sides);
    }
    Ok((!requested).then_some(all))
}

/// The curve fraction in `(fa, fb)` where `edge` crosses `cut`'s row.
fn crossing_fraction(
    brep: &BRep,
    surface: &SurfaceGeometry,
    cut: &Cut,
    edge: EdgeHandle,
    (fa, fb): (f64, f64),
) -> Option<f64> {
    let curve = brep.edges.get(edge)?.curve.as_ref()?;
    let gap = |f: f64| {
        let p = curve
            .geometry
            .point_at(curve.t_start + (curve.t_end - curve.t_start) * f);
        surface
            .closest_point(&p)
            .map(|uv| (cut.coordinate(uv) - cut.cap.at) * cut.toward_band)
    };
    let (mut lo, mut hi) = (fa, fb);
    let (mut glo, ghi) = (gap(lo)?, gap(hi)?);
    if glo.signum() == ghi.signum() {
        return None;
    }
    for _ in 0..80 {
        let mid = 0.5 * (lo + hi);
        let g = gap(mid)?;
        if g.signum() == glo.signum() {
            lo = mid;
            glo = g;
        } else {
            hi = mid;
        }
        if (hi - lo).abs() <= 1e-15 * (fa - fb).abs().max(1e-300) {
            break;
        }
    }
    Some(0.5 * (lo + hi))
}

/// How far the nearest turning point of any ring (a sample where its row
/// coordinate stops rising or falling, read from its position) lies from
/// `cut`'s row, in the row coordinate (infinite with none).
fn turning_clearance(cut: &Cut, surface: &SurfaceGeometry, built: &[BoundaryRing]) -> f64 {
    let mut clearance = f64::INFINITY;
    for ring in built {
        let c: Vec<f64> = ring
            .uv
            .iter()
            .zip(&ring.pos)
            .map(|(&uv, p)| cut.coordinate(surface.closest_point(p).unwrap_or(uv)))
            .collect();
        let n = c.len();
        for i in 0..n {
            let (a, b, d) = (c[(i + n - 1) % n], c[i], c[(i + 1) % n]);
            if (b - a) * (d - b) <= 0.0 {
                clearance = clearance.min((b - cut.cap.at).abs());
            }
        }
    }
    clearance
}

/// Whether `cut`'s row lies in the face, by parity in the band chart (see
/// the call site): a reference point a hair to the face side of a ring
/// sample lying in the band (the face lies left of its loops when the
/// flag is set), the vertical segment from it to the row at its own
/// along coordinate, and the rings' crossings of that segment. `false`
/// when no ring sample lies in the band.
fn row_in_face_by_parity(
    cut: &Cut,
    surface: &SurfaceGeometry,
    built: &[BoundaryRing],
    orientation: bool,
) -> bool {
    let coordinate_of =
        |uv: (f64, f64), p: &Point3| cut.coordinate(surface.closest_point(p).unwrap_or(uv));
    let in_band = |c: f64| (c - cut.cap.at) * cut.toward_band > cut.eps;
    // The reference: the band sample of the longest ring segment (a
    // well-conditioned left side), offset a thousandth of its length.
    let mut reference: Option<((f64, f64), f64)> = None;
    for ring in built {
        let n = ring.uv.len();
        for i in 0..n {
            let (a, b) = (ring.uv[i], ring.uv[(i + 1) % n]);
            let len = (b.0 - a.0).phypot(b.1 - a.1);
            if len > 0.0
                && in_band(coordinate_of(a, &ring.pos[i]))
                && in_band(coordinate_of(b, &ring.pos[(i + 1) % n]))
                && reference.is_none_or(|(_, l)| len > l)
            {
                let m = (0.5 * (a.0 + b.0), 0.5 * (a.1 + b.1));
                let side = if orientation { 1.0 } else { -1.0 };
                let left = (-(b.1 - a.1) / len * side, (b.0 - a.0) / len * side);
                reference = Some(((m.0 + 1e-3 * len * left.0, m.1 + 1e-3 * len * left.1), len));
            }
        }
    }
    let Some((r, _)) = reference else {
        return false;
    };
    // The vertical segment from r to the row, in (along, coordinate).
    let (ra, rc) = (cut.along(r), cut.coordinate(r));
    let (lo, hi) = (rc.min(cut.cap.at), rc.max(cut.cap.at));
    let period = cut.period;
    let mut crossings = 0usize;
    for ring in built {
        let n = ring.uv.len();
        for i in 0..n {
            let (a, b) = (ring.uv[i], ring.uv[(i + 1) % n]);
            let (aa, ac, mut ba, bc) = (
                cut.along(a),
                cut.coordinate(a),
                cut.along(b),
                cut.coordinate(b),
            );
            // The segment on one branch, then every branch reaching ra.
            if let Some(p) = period {
                ba = aa + popkern_core::geometry::chart::atlas::along_offset(Some(p), aa, ba, 1.0);
                if ba - aa > 0.5 * p {
                    ba -= p;
                }
            }
            let shifts: Vec<f64> = match period {
                Some(p) => {
                    let k = ((ra - aa.min(ba)) / p).floor();
                    vec![k * p, (k - 1.0) * p, (k + 1.0) * p]
                }
                None => vec![0.0],
            };
            for s in shifts {
                let (x0, x1) = (aa + s, ba + s);
                if (x0 - ra) * (x1 - ra) >= 0.0 || x0 == x1 {
                    continue;
                }
                let t = (ra - x0) / (x1 - x0);
                let c = ac + (bc - ac) * t;
                if c > lo && c < hi {
                    crossings += 1;
                }
            }
        }
    }
    crossings.is_multiple_of(2)
}
