//! Fillet strips around closed rims (PK-890,
//! `docs/architecture/closed-nurbs-rings.md` §5): a cubic loft that
//! interpolates a cyclic station sequence as a periodic ring.
//!
//! A clamped loft through a closed sequence meets itself at the seam
//! with a tangent kink, and the PK-141 chart overlap past its ends
//! makes it hold its seam twice. This loft is the periodic uniform
//! cubic through the stations, stored in periodic form
//! ([`NurbsSurface::u_period`]): C2 across its seam like everywhere
//! else, and periodic for every pcurve and boolean consumer. Cubic, not
//! quadratic: between stations a quadratic strip parted from the
//! cubic rails through the same contacts by 6e-5 on a 25 deg angled
//! rim, past the 1e-5 rail gate.

use popkern_core::geometry::{NurbsSurface, Point3, SurfaceEval};

use crate::contact::ContactSolution;
use crate::g2_cross_section::{
    G2Error, G2Request, g2_cross_section_bezier, g2_sweep_from_contact_solutions,
};

/// Gauss–Seidel sweeps for the cyclic tridiagonal interpolation
/// system. The loft requires its diagonal to exceed 1.5 times the
/// off-diagonal sum, so each sweep contracts by at least 2/3 and 128
/// sweeps reach roundoff for any station count.
const CYCLIC_SWEEPS: usize = 128;

/// The periodic cubic loft through `rows`, a cyclic station sequence
/// (the closing station NOT repeated). `u` runs along the stations
/// with period 1 over the domain `[0, 1]` by chord length (the mean of
/// the first and last columns' chords), station 0 on the seam at
/// `u = 0`; `v` keeps the rows' arc layout (`degree_v`,
/// `knots_v`, m points and weights per row). `None` for fewer than
/// three stations, ragged rows, a repeated station, spacing too uneven
/// to interpolate, or a degenerate weight.
pub(crate) fn closed_station_loft(
    rows: &[Vec<Point3>],
    row_weights: &[Vec<f64>],
    degree_v: usize,
    knots_v: &[f64],
) -> Option<NurbsSurface> {
    let n = rows.len();
    if n < 3 || row_weights.len() != n {
        return None;
    }
    let m = rows[0].len();
    if m < 2
        || knots_v.len() != m + degree_v + 1
        || rows.iter().any(|r| r.len() != m)
        || row_weights.iter().any(|w| w.len() != m)
    {
        return None;
    }
    let data: Vec<Vec<[f64; 4]>> = rows
        .iter()
        .zip(row_weights)
        .map(|(row, ws)| {
            row.iter()
                .zip(ws)
                .map(|(p, &w)| [p.x * w, p.y * w, p.z * w, w])
                .collect()
        })
        .collect();
    // Chord-length station parameters on the period [0, 1): the mean
    // of the two rails' chords. A uniform parameter on an angled rim's
    // unevenly spaced stations left the strip 3e-5 off rails that the
    // walker's own chord-length check had at 2.5e-6.
    let chord = |i: usize| {
        let (a, b) = (&rows[i], &rows[(i + 1) % n]);
        0.5 * (a[0].distance_to(&b[0]) + a[m - 1].distance_to(&b[m - 1]))
    };
    let steps: Vec<f64> = (0..n).map(chord).collect();
    let total: f64 = steps.iter().sum();
    if !(total > 0.0 && steps.iter().all(|&h| h > 0.0)) {
        return None;
    }
    let h: Vec<f64> = steps.iter().map(|&s| s / total).collect();
    let mut tau = Vec::with_capacity(n);
    let mut acc = 0.0;
    for &step in &h {
        tau.push(acc);
        acc += step;
    }
    // Periodic cubic with simple knots at the stations: station i sits
    // on the middle knot of its control row Q_i, so
    // D_i = a_i Q_{i-1} + b_i Q_i + c_i Q_{i+1} with the at-knot values
    // of the cubic B-splines (spacings h around station i).
    let hm = |k: isize| h[k.rem_euclid(n as isize) as usize];
    let coeffs: Vec<(f64, f64, f64)> = (0..n as isize)
        .map(|i| {
            let (h2, h1, h0, h_1) = (hm(i - 2), hm(i - 1), hm(i), hm(i + 1));
            let a = h0 * h0 / ((h1 + h0) * (h2 + h1 + h0));
            let c = h1 * h1 / ((h1 + h0) * (h1 + h0 + h_1));
            (a, 1.0 - a - c, c)
        })
        .collect();
    // Gauss–Seidel needs the diagonal to dominate; spacing that uneven
    // is not a station walk.
    if coeffs.iter().any(|&(a, b, c)| b <= 1.5 * (a + c)) {
        return None;
    }
    let mut q = data.clone();
    for _ in 0..CYCLIC_SWEEPS {
        for (i, &(a, b, c)) in coeffs.iter().enumerate() {
            let (prev, next) = ((i + n - 1) % n, (i + 1) % n);
            for j in 0..m {
                for k in 0..4 {
                    q[i][j][k] = (data[i][j][k] - a * q[prev][j][k] - c * q[next][j][k]) / b;
                }
            }
        }
    }
    // Periodic form: row k is Q_{k-1} (cyclic) for k = 0..n+2, so the
    // last three rows repeat the first three, over the station knots
    // extended one period each way: knot j is tau_{j-3}. Row k's middle
    // knot is tau_{k-1}.
    let knot = |j: usize| {
        let i = j as isize - 3;
        let turns = i.div_euclid(n as isize) as f64;
        tau[i.rem_euclid(n as isize) as usize] + turns
    };
    let knots_u: Vec<f64> = (0..n + 7).map(knot).collect();
    let mut control_points = Vec::with_capacity((n + 3) * m);
    let mut weights = Vec::with_capacity((n + 3) * m);
    for k in 0..n + 3 {
        for h in &q[(k + n - 1) % n] {
            let w = h[3];
            if !(w.is_finite() && w.abs() > f64::MIN_POSITIVE) {
                return None;
            }
            control_points.push(Point3::new(h[0] / w, h[1] / w, h[2] / w));
            weights.push(w);
        }
    }
    Some(NurbsSurface {
        degree_u: 3,
        degree_v,
        knots_u,
        knots_v: knots_v.to_vec(),
        control_points,
        num_u: n + 3,
        num_v: m,
        weights,
    })
}

/// Stations closer than this fraction of the median station spacing
/// are one station: a closed walk ends on a repeat of its first
/// contact (measured ~1e-7 apart on an r 5 rim at 0.4 spacing).
const SAME_STATION: f64 = 1e-3;

/// The periodic ring through a station-row loft (`u` along the
/// stations, one control row per station, as `build_nurbs_sweep`
/// builds it) when both rails close: their last footprint repeats the
/// first. The repeat is dropped and the remaining stations lofted
/// with [`closed_station_loft`]. `None` for an open sequence or a
/// loft that is not one row per station.
pub(crate) fn closed_ring_through_stations(
    loft: &NurbsSurface,
    footprint_1: &[Point3],
    footprint_2: &[Point3],
) -> Option<NurbsSurface> {
    let n = loft.num_u;
    let m = loft.num_v;
    if n < 4
        || footprint_1.len() != n
        || footprint_2.len() != n
        || loft.control_points.len() != n * m
    {
        return None;
    }
    if !(stations_close(footprint_1) && stations_close(footprint_2)) {
        return None;
    }
    let rows: Vec<Vec<Point3>> = (0..n - 1)
        .map(|i| loft.control_points[i * m..(i + 1) * m].to_vec())
        .collect();
    let row_weights: Vec<Vec<f64>> = (0..n - 1)
        .map(|i| loft.weights[i * m..(i + 1) * m].to_vec())
        .collect();
    closed_station_loft(&rows, &row_weights, loft.degree_v, &loft.knots_v)
}

/// True when a station sequence ends on a repeat of its first station
/// ([`SAME_STATION`] of its median spacing).
fn stations_close(pts: &[Point3]) -> bool {
    let (Some(first), Some(last)) = (pts.first(), pts.last()) else {
        return false;
    };
    let mut steps: Vec<f64> = pts.windows(2).map(|w| w[0].distance_to(&w[1])).collect();
    steps.sort_by(f64::total_cmp);
    let median = steps.get(steps.len() / 2).copied().unwrap_or(0.0);
    pts.len() >= 4 && median > 0.0 && first.distance_to(last) <= SAME_STATION * median
}

/// The G2 strip through `solutions`: the periodic ring through their
/// G2 cross-sections when the stations close (PK-884), else the
/// clamped cubic loft [`g2_sweep_from_contact_solutions`]. A clamped
/// loft through a closed rim held its seam twice and closed only to
/// its fit (1.9e-5); its seam ends projected to opposite chart ends.
pub(crate) fn g2_sweep_or_ring(
    solutions: &[ContactSolution],
    surface_0: &dyn SurfaceEval,
    surface_1: &dyn SurfaceEval,
) -> Result<NurbsSurface, G2Error> {
    let contacts = |pick: fn(&ContactSolution) -> Point3| -> Vec<Point3> {
        solutions.iter().map(pick).collect()
    };
    let closed =
        stations_close(&contacts(|s| s.contact0)) && stations_close(&contacts(|s| s.contact1));
    if closed {
        let mut sections = Vec::with_capacity(solutions.len() - 1);
        for solution in &solutions[..solutions.len() - 1] {
            let request = G2Request {
                surface_0,
                surface_1,
                solution,
            };
            sections.push(g2_cross_section_bezier(&request)?.curve);
        }
        let rows: Vec<Vec<Point3>> = sections.iter().map(|c| c.control_points.clone()).collect();
        let weights: Vec<Vec<f64>> = sections.iter().map(|c| c.weights.clone()).collect();
        if let Some(first) = sections.first()
            && let Some(ring) = closed_station_loft(&rows, &weights, first.degree, &first.knots)
        {
            return Ok(ring);
        }
    }
    g2_sweep_from_contact_solutions(solutions, surface_0, surface_1)
}

#[cfg(test)]
#[path = "closed_station_loft_tests.rs"]
mod tests;
