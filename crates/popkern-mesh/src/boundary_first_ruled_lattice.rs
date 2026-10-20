//! Interior lattice of RULED boundary-first charts (cylinder and
//! cone walls, extruded NURBS), split from `boundary_first.rs` under
//! its frozen line budget (PK-651).

use popkern_core::geometry::Point3;

use super::boundary_first_ring_index::RingIndex;
use super::{BoundaryRing, dist_to_rings};
use crate::TessellationOptions;

/// Interior points of a ruled chart over the bounding box
/// `(u0, u1, v0, v1)`: columns at the boundary's own u-stations, rows
/// at its v-cadence, filtered by `ring_index` (region and clearance),
/// pushed onto `uv_pts`/`pos_pts`.
///
/// Returns the factor by which the chart's Delaunay compresses v: the
/// lattice's median column spacing over its widest row gap (at most 1). Rows at
/// the boundary's v-cadence can sit 14x farther apart than the
/// columns (0.88 x 12.5 mm on a slotted r 37.5 wall), and an isotropic
/// Delaunay fans any off-row point across many columns (facets 0.035
/// off the surface; PK-651). Along a ruling the facet error does not
/// grow with length, so a lattice-square metric is exact there.
#[allow(clippy::too_many_arguments)]
pub(super) fn ruled_lattice(
    rings: &[BoundaryRing],
    ring_index: &RingIndex,
    (u0, u1, v0, v1): (f64, f64, f64, f64),
    (su, sv): (f64, f64),
    cells: (f64, f64),
    options: &TessellationOptions,
    uv_pts: &mut Vec<(f64, f64)>,
    pos_pts: &mut Vec<Option<Point3>>,
) -> f64 {
    // RULED charts (cylinder/cone walls, extruded NURBS): the
    // interior lattice reuses the BOUNDARY's u-stations — the
    // adaptive edge sampling already made those as dense as the
    // curvature demands, while the bbox grid is capped at
    // max_subdivisions and its coarser u-spacing carries 16× the
    // chord error (measured: 4.6e-3 vs a 2e-3 volume gate on the
    // chamfered column). Along the ruling v is STRAIGHT, so a few
    // rows cost no accuracy and give the Delaunay interior
    // support against trim bites (measured: the breached-bore
    // wall membraned with boundary-only triangulation even after
    // the trim splits).
    let mut u_stations: Vec<f64> = rings
        .iter()
        .flat_map(|r| r.uv.iter().map(|&(u, _)| u))
        .collect();
    u_stations.sort_by(f64::total_cmp);
    // Dedup PROPORTIONALLY to the chart span: the two rims invert
    // to near-identical u (differing ~1e-8 from independent
    // atan2 of separately-sampled 3D points), and near-duplicate
    // interior points corrupt the Bowyer-Watson cavity
    // (measured: non-Delaunay 0.63-rad chords beside the seam
    // with a dense kept lattice — the PK-146 seam band).
    // Near-twin rim stations (0.03 mm apart on socket) merge at a quarter u-cell, keeping a tab's columns under the cap (PK-651).
    let station_merge = ((u1 - u0).abs() * 1.0e-7)
        .max(1e-9)
        .max(0.25 * cells.0 / su);
    u_stations.dedup_by(|a, b| (*a - *b).abs() < station_merge);
    if u_stations.len() > 512 {
        let step = u_stations.len() as f64 / 512.0;
        u_stations = (0..512)
            .map(|i| u_stations[(i as f64 * step) as usize])
            .collect();
    }
    let n_v = options.min_subdivisions.max(4);
    let seam_trace =
        popkern_core::trace::trace_enabled("boundary-first-seam") && u_stations.len() > 100;
    let mut kill = [0usize; 3]; // region, clearance, kept — seam band only
    // Interior v-stations CONFORM to the boundary's v-cadence:
    // with a fixed sparse row count, ruling-side boundary samples
    // (a refined seam at half the row spacing) sit exactly midway
    // between rows and the whole seam neighborhood becomes one
    // near-cocircular family — every triangulation is legal and
    // Bowyer-Watson picks long-u fans through the empty inter-row
    // band (measured: 0.63-rad chords beside the bore seam,
    // volume +1.3e-3 at any chord tolerance, ZERO circumcircle
    // violations — the fans are legal). Taking the rows from the
    // ring's own interior v-values tiles the lattice conformally
    // with the boundary; the fixed grid remains the fallback for
    // boundaries with no interior v-samples.
    let mut v_stations: Vec<f64> = rings
        .iter()
        .flat_map(|r| r.uv.iter().map(|&(_, v)| v))
        .filter(|&v| v > v0 + 1e-9 && v < v1 - 1e-9)
        .collect();
    v_stations.sort_by(f64::total_cmp);
    // Rows merge at a thousandth of the lattice row cell: the two
    // rails' paired stations invert from independently sampled 3D
    // points and differ by inversion noise, which a 1e-7-of-span merge
    // did not always absorb — a 0.32 fillet wall's mid-rail stations
    // read v 1.15999998 / 1.16000002 (3.4e-8 apart) and the lattice
    // emitted two coincident rows, 30 zero-area slivers (PK-626).
    // Genuine row clusters (a skew rim's 0.009 cadence) sit far above.
    let v_merge = ((v1 - v0).abs() * 1.0e-7)
        .max(1e-9)
        .max(1.0e-3 * cells.1 / sv);
    v_stations.dedup_by(|a, b| (*a - *b).abs() < v_merge);
    if v_stations.len() < n_v - 1 {
        v_stations = (1..n_v)
            .map(|j| v0 + (v1 - v0) * j as f64 / n_v as f64)
            .collect();
    } else if v_stations.len() > 64 {
        let step = v_stations.len() as f64 / 64.0;
        v_stations = (0..64)
            .map(|i| v_stations[(i as f64 * step) as usize])
            .collect();
    }
    for &u in &u_stations {
        if u <= u0 || u >= u1 {
            continue;
        }
        for &v in &v_stations {
            let in_band = u - u0 < 0.8 || u1 - u < 0.8;
            let detail = seam_trace && u1 - u < 0.9 && (v - 1.0).abs() < 0.01;
            if !ring_index.inside((u, v)) {
                if seam_trace && in_band {
                    kill[0] += 1;
                }
                if detail {
                    eprintln!("    [bf-seam] ({u:.3},{v:.3}) KILLED region");
                }
            } else if ring_index.within_clearance((u, v)) {
                if seam_trace && in_band {
                    kill[1] += 1;
                }
                if detail {
                    eprintln!(
                        "    [bf-seam] ({u:.3},{v:.3}) KILLED clearance d={:.4}",
                        dist_to_rings((u, v), rings, su, sv)
                    );
                }
            } else {
                if seam_trace && in_band {
                    kill[2] += 1;
                }
                if detail {
                    eprintln!("    [bf-seam] ({u:.3},{v:.3}) KEPT");
                }
                uv_pts.push((u, v));
                pos_pts.push(None);
            }
        }
    }
    if seam_trace {
        eprintln!(
            "  [bf-seam] stations={} n_v={n_v} cells={cells:?} band kills: region={} clearance={} kept={}",
            u_stations.len(),
            kill[0],
            kill[1],
            kill[2]
        );
    }
    // Median metric gap of the kept columns.
    let median_gap = |stations: &[f64], scale: f64| -> Option<f64> {
        let mut gaps: Vec<f64> = stations
            .windows(2)
            .map(|w| scale * (w[1] - w[0]).abs())
            .filter(|&g| g > 0.0)
            .collect();
        gaps.sort_by(f64::total_cmp);
        gaps.get(gaps.len() / 2).copied()
    };
    // Rows are as sparse as the boundary's v-cadence, and unevenly so
    // (3.45 apart beside a 0.009-spaced cluster where a skew rim runs):
    // the WIDEST row gap must look no taller than a column gap.
    let mut rows: Vec<f64> = Vec::with_capacity(v_stations.len() + 2);
    rows.push(v0);
    rows.extend(v_stations.iter().copied());
    rows.push(v1);
    let max_row_gap = rows
        .windows(2)
        .map(|w| sv * (w[1] - w[0]).abs())
        .fold(0.0_f64, f64::max);
    match median_gap(&u_stations, su) {
        Some(cu) if max_row_gap > cu => cu / max_row_gap,
        _ => 1.0,
    }
}

/// Whether the chart's v iso-lines are straight across the whole box
/// (five u positions, each within `chord_height` of its end-to-end
/// chord at every [`v_line_fractions`] sample): the condition under
/// which compressing v in the Delaunay metric cannot lengthen a
/// facet's error. The ruled-path test samples the u-midline only.
pub(super) fn ruled_everywhere(
    surface: &popkern_core::geometry::surface::SurfaceGeometry,
    (u0, u1, v0, v1): (f64, f64, f64, f64),
    chord_height: f64,
) -> bool {
    let fractions = v_line_fractions(surface, v0, v1);
    [0.1, 0.3, 0.5, 0.7, 0.9].iter().all(|&fu| {
        v_line_straight(
            surface,
            u0 + (u1 - u0) * fu,
            (v0, v1),
            &fractions,
            chord_height,
        )
    })
}

/// Whether the v iso-line at `u` lies within `chord_height` of its
/// end-to-end chord at each of `fractions` of `(v0, v1)`.
pub(super) fn v_line_straight(
    surface: &popkern_core::geometry::surface::SurfaceGeometry,
    u: f64,
    (v0, v1): (f64, f64),
    fractions: &[f64],
    chord_height: f64,
) -> bool {
    use popkern_core::geometry::surface::SurfaceEval;
    let (pa, pb) = (surface.point_at(u, v0), surface.point_at(u, v1));
    fractions.iter().all(|&fv| {
        let p = surface.point_at(u, v0 + (v1 - v0) * fv);
        let on_chord = Point3::new(
            pa.x + (pb.x - pa.x) * fv,
            pa.y + (pb.y - pa.y) * fv,
            pa.z + (pb.z - pa.z) * fv,
        );
        p.distance_to(&on_chord) <= chord_height
    })
}

/// Fractions of `(v0, v1)` at which a v iso-line's straightness is
/// sampled: its quarters, and on a NURBS chart the quarters of every v
/// knot span the window meets. A line wound through many knot spans
/// (a tap's thread flank, a 142 mm helix six turns round) can return
/// to its chord at fixed fractions: the quarters alone read 5.7e-4 on
/// a line 6.2 mm off it, the chart took the ruled Delaunay metric, and
/// its facets sagged 0.23 at a 0.011 chord (PK-900).
pub(super) fn v_line_fractions(
    surface: &popkern_core::geometry::surface::SurfaceGeometry,
    v0: f64,
    v1: f64,
) -> Vec<f64> {
    let mut fractions = vec![0.25, 0.5, 0.75];
    if let popkern_core::geometry::surface::SurfaceGeometry::Nurbs(n) = surface {
        fractions.extend(super::boundary_first_lattice::span_fractions(
            &n.knots_v, v0, v1,
        ));
    }
    fractions
}

#[cfg(test)]
mod tests {
    use super::*;
    use popkern_core::topology::EdgeHandle;

    /// A quarter-turn r 0.08 wall strip, u in [0, pi/2], v in [1, 1.32],
    /// whose two rails carry their mid stations at v 1.16 -+ 1.7e-8 —
    /// the inversion noise of independently sampled rails (PK-626).
    fn noisy_rail_strip() -> BoundaryRing {
        let (u1, r) = (std::f64::consts::FRAC_PI_2, 0.08);
        let mut uv = Vec::new();
        for k in 0..=8 {
            uv.push((u1 * k as f64 / 8.0, 1.0));
        }
        uv.extend([(u1, 1.08), (u1, 1.16 - 1.7e-8), (u1, 1.24)]);
        for k in (0..=8).rev() {
            uv.push((u1 * k as f64 / 8.0, 1.32));
        }
        uv.extend([(0.0, 1.24), (0.0, 1.16 + 1.7e-8), (0.0, 1.08)]);
        let pos = uv
            .iter()
            .map(|&(u, v)| Point3::new(r * u.cos(), r * u.sin(), v))
            .collect();
        let sources = vec![(EdgeHandle::default(), 0.0, 1.0); uv.len()];
        BoundaryRing { uv, pos, sources }
    }

    #[test]
    fn paired_rail_stations_within_inversion_noise_make_one_row() {
        let rings = vec![noisy_rail_strip()];
        let (su, sv) = (0.08, 1.0);
        let cells = (su * std::f64::consts::FRAC_PI_2 / 8.0, 0.08);
        let index = RingIndex::new(&rings, su, sv, cells);
        let (mut uv_pts, mut pos_pts) = (Vec::new(), Vec::new());
        ruled_lattice(
            &rings,
            &index,
            (0.0, std::f64::consts::FRAC_PI_2, 1.0, 1.32),
            (su, sv),
            cells,
            &TessellationOptions::default(),
            &mut uv_pts,
            &mut pos_pts,
        );
        let mut rows: Vec<f64> = uv_pts.iter().map(|&(_, v)| v).collect();
        rows.sort_by(f64::total_cmp);
        rows.dedup();
        assert_eq!(rows.len(), 3, "one interior row per rail station: {rows:?}");
    }
}
