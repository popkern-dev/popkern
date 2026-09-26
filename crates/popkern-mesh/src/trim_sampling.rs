//! Trim-loop boundary sampling (split out of `trim.rs` for the
//! 500-line ratchet): the canonical-first loop sampler and the
//! adaptive pcurve subdivision fallback.

use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::uv_projection::project_point_to_uv;
use popkern_core::topology::{BRep, LoopHandle};

/// Sample a loop; canonical stations first (PK-299), adaptive
/// fallback on decline.
pub(crate) fn sample_loop_with_samples(
    brep: &BRep,
    loop_handle: LoopHandle,
    surface: &SurfaceGeometry,
    chord_height: f64,
    edge_samples: Option<
        &std::collections::HashMap<
            popkern_core::topology::EdgeHandle,
            crate::boundary_first::EdgeSamples,
        >,
    >,
) -> (Vec<(f64, f64)>, Vec<Point3>, bool) {
    if let Some(samples) = edge_samples
        && let Some((uv, p3)) = crate::trim_canonical::sample_loop_canonical(
            brep,
            loop_handle,
            surface,
            samples,
            chord_height,
        )
    {
        return (uv, p3, true);
    }
    let Some(loop_) = brep.loops.get(loop_handle) else {
        return (Vec::new(), Vec::new(), false);
    };
    let start_he = loop_.half_edge;
    if !brep.half_edges.contains_key(start_he) {
        return (Vec::new(), Vec::new(), false);
    }

    // Check if ALL half-edges in this loop have pcurves.  If any is
    // missing, use the 3D-to-UV projection fallback for the ENTIRE
    // loop to avoid mixing two UV sources (pcurve vs projection) which
    // can produce inconsistent coordinates and broken polygons.
    let all_have_pcurves = {
        let mut check = start_he;
        let mut all = true;
        loop {
            if brep.half_edges[check].pcurve.is_none() {
                all = false;
                break;
            }
            check = brep.half_edges[check].next;
            if check == start_he {
                break;
            }
        }
        all
    };

    let mut uv_points = Vec::new();
    let mut points_3d = Vec::new();
    let mut current = start_he;

    loop {
        let he = &brep.half_edges[current];

        if all_have_pcurves {
            if let Some(ref pcurve) = he.pcurve {
                let prev_len = uv_points.len();
                // A member spanning a large fraction of a periodic axis
                // must sample FINER than half the period: its UV image
                // may be perfectly straight (a full rim is one chord at
                // constant v), the adaptive subdivision then emits only
                // the two endpoints, and `unwrap_periodic` reads the
                // full-period step as a seam wrap and FOLDS the member
                // onto its start — the trim polygon collapses and the
                // face meshes untrimmed (PK-332: the seam-slice crack
                // on every revolution primitive traces here).
                let (t0p, t1p) = pcurve.parameter_range();
                let a = pcurve.point_at(t0p);
                let b = pcurve.point_at(t1p);
                let (period_u, period_v) = surface.periodicity();
                let wide =
                    |d: f64, p: Option<f64>| p.is_some_and(|p| p > 0.0 && d.abs() >= 0.4 * p);
                if wide(b.x - a.x, period_u) || wide(b.y - a.y, period_v) {
                    let n = 32;
                    for i in 0..=n {
                        let t = t0p + (t1p - t0p) * (i as f64) / (n as f64);
                        let q = pcurve.point_at(t);
                        uv_points.push((q.x, q.y));
                    }
                } else {
                    adaptive_sample_pcurve(pcurve, surface, chord_height, &mut uv_points);
                }
                // Evaluate surface at each new UV point to get 3D positions.
                for uv in &uv_points[prev_len..] {
                    points_3d.push(surface.point_at(uv.0, uv.1));
                }
            }
        } else {
            // Fallback: project 3D edge curve (or vertex) to UV.
            let edge = &brep.edges[he.edge];
            if let Some(ref curve_on_edge) = edge.curve {
                let samples = 32_usize;
                let curve_start = curve_on_edge.geometry.point_at(curve_on_edge.t_start);
                let he_origin = brep.vertices[he.origin].position;
                let forward = brep.tolerance.points_coincident(&curve_start, &he_origin);
                let (t_lo, t_hi) = if forward {
                    (curve_on_edge.t_start, curve_on_edge.t_end)
                } else {
                    (curve_on_edge.t_end, curve_on_edge.t_start)
                };

                for i in 0..samples {
                    let t = t_lo + (t_hi - t_lo) * (i as f64 / samples as f64);
                    let p = curve_on_edge.geometry.point_at(t);
                    if let Some(uv) = project_point_to_uv(surface, &p) {
                        uv_points.push(uv);
                        points_3d.push(p);
                    }
                }
            } else {
                let p = brep.vertices[he.origin].position;
                if let Some(uv) = project_point_to_uv(surface, &p) {
                    uv_points.push(uv);
                    points_3d.push(p);
                }
            }
        }

        current = he.next;
        if current == start_he {
            break;
        }
    }

    (uv_points, points_3d, all_have_pcurves)
}

/// Adaptively sample a pcurve using chord-height subdivision in both
/// its uv chart and its 3D image on `surface`.
///
/// Pushes the start point, then recursively subdivides spans whose uv
/// chord-height error, or whose surface image's 3D sag, exceeds
/// `chord_height`, up to [`MAX_TRIM_SUBDIVISION_DEPTH`]. Finally pushes
/// the end point.
///
/// The 3D criterion is what the trim clip needs: vertices created where
/// the trim cuts grid triangles take their 3D position by linear
/// interpolation along the trim polyline, so the polyline's 3D chords
/// must follow the surface (PK-560: a helical thread flank's trim is a
/// straight uv line over several turns; sampled by uv straightness it
/// got two points, and the 3D chord between them cut across the hole —
/// clip vertices landed on the hole axis, 2.57 off the flank). The span
/// is probed at its quarter points as well as its midpoint: after an
/// even number of helix turns the midpoint alone lies back on the chord.
pub(crate) fn adaptive_sample_pcurve(
    pcurve: &popkern_core::geometry::pcurve::Pcurve,
    surface: &SurfaceGeometry,
    chord_height: f64,
    output: &mut Vec<(f64, f64)>,
) {
    let (t_start, t_end) = pcurve.parameter_range();
    let p_start = pcurve.point_at(t_start);
    output.push((p_start.x, p_start.y));

    subdivide_pcurve(pcurve, surface, t_start, t_end, chord_height, 0, output);

    let p_end = pcurve.point_at(t_end);
    output.push((p_end.x, p_end.y));
}

/// Deepest trim subdivision (4096 spans per pcurve).
const MAX_TRIM_SUBDIVISION_DEPTH: usize = 12;

/// Depth to which the historical uv chord-height criterion applies.
const UV_SUBDIVISION_DEPTH: usize = 6;

/// Recursive subdivision helper for [`adaptive_sample_pcurve`].
fn subdivide_pcurve(
    pcurve: &popkern_core::geometry::pcurve::Pcurve,
    surface: &SurfaceGeometry,
    t0: f64,
    t1: f64,
    chord_height: f64,
    depth: usize,
    output: &mut Vec<(f64, f64)>,
) {
    if depth >= MAX_TRIM_SUBDIVISION_DEPTH {
        return;
    }

    let t_mid = 0.5 * (t0 + t1);
    let p0 = pcurve.point_at(t0);
    let p1 = pcurve.point_at(t1);
    let p_mid = pcurve.point_at(t_mid);

    // uv chord-height error of the midpoint.
    let uv_error = if depth < UV_SUBDIVISION_DEPTH {
        let dx_chord = p1.x - p0.x;
        let dy_chord = p1.y - p0.y;
        let chord_len = (dx_chord * dx_chord + dy_chord * dy_chord).sqrt();
        let dx_mid = p_mid.x - p0.x;
        let dy_mid = p_mid.y - p0.y;
        if chord_len > 1e-30 {
            (dx_chord * dy_mid - dy_chord * dx_mid).abs() / chord_len
        } else {
            (dx_mid * dx_mid + dy_mid * dy_mid).sqrt()
        }
    } else {
        0.0
    };

    // 3D sag of the surface image at the quarter and mid points.
    let image = |t: f64| {
        let q = pcurve.point_at(t);
        surface.point_at(q.x, q.y)
    };
    let (s0, s1) = (image(t0), image(t1));
    let sag = [0.25, 0.5, 0.75]
        .iter()
        .map(|&f| {
            let on_chord = Point3::new(
                s0.x + f * (s1.x - s0.x),
                s0.y + f * (s1.y - s0.y),
                s0.z + f * (s1.z - s0.z),
            );
            image(t0 + f * (t1 - t0)).distance_to(&on_chord)
        })
        .fold(0.0_f64, f64::max);

    if uv_error > chord_height || sag > chord_height {
        subdivide_pcurve(pcurve, surface, t0, t_mid, chord_height, depth + 1, output);
        output.push((p_mid.x, p_mid.y));
        subdivide_pcurve(pcurve, surface, t_mid, t1, chord_height, depth + 1, output);
    }
}
