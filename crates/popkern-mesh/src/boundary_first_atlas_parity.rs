//! Whether a cap row no ring crosses lies in the face, read by parity in
//! the band chart (chart-service.md §6a, step 2, until W2c moves the
//! point-in-face authority onto the chart service, PK-1273).

use popkern_core::geometry::Point3;
use popkern_core::geometry::chart::Axis;
use popkern_core::geometry::chart::atlas::along_offset;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_math::PortableMath;

use super::super::BoundaryRing;
use super::Cut;

/// `b` moved by whole periods along the row to the branch nearest `a`.
fn along_step(cut: &Cut, a: (f64, f64), b: (f64, f64)) -> (f64, f64) {
    let Some(p) = cut.period else { return b };
    let d = along_offset(Some(p), cut.along(a), cut.along(b), 1.0);
    let d = if d > 0.5 * p { d - p } else { d };
    let along = cut.along(a) + d;
    match cut.cap.singular.axis {
        Axis::V => (along, b.1),
        Axis::U => (b.0, along),
    }
}

/// Whether `cut`'s row lies in the face, by parity in the band chart (see
/// the call site): a reference point a hair to the face side of a ring
/// sample lying in the band (the face lies left of its loops when the
/// flag is set), the vertical segment from it to the row at its own
/// along coordinate, and the rings' crossings of that segment. `false`
/// when no ring sample lies in the band.
pub(super) fn row_in_face_by_parity(
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
            // The step on one branch: across the seam a raw step runs back
            // nearly a period (a dome's rim: it was the longest, its left
            // the wrong pole).
            let (a, b) = (
                ring.uv[i],
                along_step(cut, ring.uv[i], ring.uv[(i + 1) % n]),
            );
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
        // Each sample's along coordinate unwrapped once round the ring
        // (each step to the branch nearest the last): a sample two
        // segments share is one number, so the half-open test counts a
        // ray through it once. Unwrapped per segment, its two copies
        // differed by a rounding and the ray through it counted both (a
        // dome's inner face read the lower row as in the face).
        let mut along = Vec::with_capacity(n + 1);
        if let Some(first) = ring.uv.first() {
            along.push(cut.along(*first));
        }
        for i in 0..n {
            let (prev, next) = (along[i], cut.along(ring.uv[(i + 1) % n]));
            along.push(match period {
                Some(p) => {
                    let d = popkern_core::geometry::chart::atlas::along_offset(
                        Some(p),
                        prev,
                        next,
                        1.0,
                    );
                    prev + if d > 0.5 * p { d - p } else { d }
                }
                None => next,
            });
        }
        for i in 0..n {
            let (ac, bc) = (
                cut.coordinate(ring.uv[i]),
                cut.coordinate(ring.uv[(i + 1) % n]),
            );
            let (aa, ba) = (along[i], along[i + 1]);
            let shifts: Vec<f64> = match period {
                Some(p) => {
                    let k = ((ra - aa.min(ba)) / p).floor();
                    vec![k * p, (k - 1.0) * p, (k + 1.0) * p]
                }
                None => vec![0.0],
            };
            for s in shifts {
                let (x0, x1) = (aa + s, ba + s);
                // Half-open in along (a segment holds its low end, not its
                // high one): a ray through a ring sample counts it once.
                // Strict on both ends, a ray through a sample counted
                // neither segment (a sphere zone's top and bottom rims
                // share their azimuths, and the reference sat on one: the
                // lower row read as in the face, PK-1225).
                if (x0 <= ra) == (x1 <= ra) {
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
