//! FAR-tier exact-junction search along a record's source curve
//! (split from `clip_endpoint_far.rs` under the 500-line structural
//! budget; the re-trim that consumes it stays there).

use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::curve_on_surface::CurveOnSurface;
use popkern_core::geometry::intersection::{SurfaceCurveIntersection, intersect_surface_curve};
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};

use crate::clip_endpoint_refine::{clone_curve_geometry, trim_boundary_edges};
use crate::closed_source_window::SourceWindow;
use crate::curve_projection::project_onto_curve_dense;

/// FAR-TIER junction search (see the call site): the nearest exact
/// boundary-edge × other-surface junction that lies ON the record's
/// source curve, unambiguous by 2×. Returns the junction point and
/// its parameter on the source curve.
#[allow(clippy::too_many_arguments)]
pub(crate) fn far_junction(
    brep_a: &BRep,
    face_a: FaceHandle,
    surf_a: &SurfaceGeometry,
    windows_a: ((f64, f64), (f64, f64)),
    brep_b: &BRep,
    face_b: FaceHandle,
    surf_b: &SurfaceGeometry,
    windows_b: ((f64, f64), (f64, f64)),
    p: &Point3,
    source: &CurveOnSurface,
    t_end: f64,
    tolerance: &Tolerance,
) -> Option<Vec<(Point3, f64)>> {
    let (s0, s1) = source.curve_3d.parameter_range();
    let chord = source
        .curve_3d
        .point_at(s0)
        .distance_to(&source.curve_3d.point_at(s1))
        .max(
            source
                .curve_3d
                .point_at(s0)
                .distance_to(&source.curve_3d.point_at(0.5 * (s0 + s1))),
        );
    let far_cap = 0.25 * chord;
    let fit_gate = tolerance.linear() * 2000.0;
    let window = SourceWindow::of(source);
    let mut cands: Vec<(f64, Point3, f64)> = Vec::new();
    let mut exact_here = false;
    for (brep, face, other_surface, other_windows) in [
        (brep_a, face_a, surf_b, windows_b),
        (brep_b, face_b, surf_a, windows_a),
    ] {
        for eh in trim_boundary_edges(brep, face) {
            let Some(edge) = brep.edges.get(eh) else {
                continue;
            };
            let Some(c) = edge.curve.as_ref() else {
                continue;
            };
            let (t0, t1) = (c.t_start, c.t_end);
            // Prune edges whose curve never comes near the endpoint
            // at far-cap scale.
            let near = (&c.geometry as &dyn CurveEval)
                .closest_point(p)
                .map(|t| {
                    let tc = t.clamp(t0.min(t1), t0.max(t1));
                    c.geometry.point_at(tc).distance_to(p) < far_cap
                })
                .unwrap_or(false);
            if !near {
                continue;
            }
            let range = (t0.min(t1), t0.max(t1));
            let Ok(result) = intersect_surface_curve(
                other_surface,
                other_windows.0,
                other_windows.1,
                &clone_curve_geometry(&c.geometry),
                range,
                tolerance,
            ) else {
                continue;
            };
            let far_trace = popkern_core::trace::trace_enabled("clip-far");
            for isect in &result.intersections {
                let SurfaceCurveIntersection::Point { point, .. } = isect else {
                    continue;
                };
                // TRIPLE-POINT REFINEMENT (PK-480): remove the
                // along-edge 1/sin(crossing angle) slide of the 1D
                // solve (see `clip_triple_junction`).
                let point = &crate::clip_triple_junction::refine_to_triple_point(
                    brep,
                    eh,
                    other_surface,
                    *point,
                    tolerance,
                );
                let d = point.distance_to(p);
                if d <= tolerance.linear() {
                    // The endpoint IS an exact junction: nothing was
                    // mislocated, so no far correction applies (the
                    // next junction along the curve is a different
                    // crossing — PK-632).
                    exact_here = true;
                    continue;
                }
                if d > far_cap {
                    continue;
                }
                // The junction must lie ON the source curve. Dense
                // projection — the trait `closest_point` clamps to
                // the seam basin on closed sources (measured,
                // PK-472: the mirrored mouth chamfer's crest circle
                // junction projected to t=0 with dj=2.73 and the
                // stop-short retrim never armed).
                let Some(t_star) = project_onto_curve_dense(&source.curve_3d, point) else {
                    if far_trace {
                        eprintln!(
                            "    [clip-far] p ({:.4},{:.4},{:.4}) junction ({:.4},{:.4},{:.4}) d {d:.3e}: no curve projection",
                            p.x, p.y, p.z, point.x, point.y, point.z
                        );
                    }
                    continue;
                };
                let t_star = t_star.clamp(s0.min(s1), s0.max(s1));
                let dj = source.curve_3d.point_at(t_star).distance_to(point);
                if far_trace {
                    eprintln!(
                        "    [clip-far] p ({:.4},{:.4},{:.4}) junction ({:.4},{:.4},{:.4}) d {d:.3e} t* {t_star:.4} dj {dj:.3e} gate {fit_gate:.3e}",
                        p.x, p.y, p.z, point.x, point.y, point.z
                    );
                }
                if dj > fit_gate {
                    continue;
                }
                // Closed source: the junction's lift next to the end
                // it corrects (PK-673 — a window through the seam).
                cands.push((d, *point, window.lift_near(t_star, t_end)));
            }
        }
    }
    if exact_here {
        return Some(Vec::new());
    }
    // Among boundary junctions ON the source curve, the clip's sloppy
    // crossing belongs to the one nearest in PARAMETER — the crossing
    // it mislocated. 3D-nearest is the wrong discriminator in a
    // junction cluster (the angled tap's runout wedge has six
    // junctions within 0.15, and a 2× spatial-ambiguity gate refused
    // every correction there).
    cands.sort_by(|x, y| {
        (x.2 - t_end)
            .abs()
            .total_cmp(&(y.2 - t_end).abs())
            .then(x.0.total_cmp(&y.0))
    });
    // Arc-consistency guard instead of a parameter-fraction cap: the
    // correction's midpoint must stay within the same spatial cap as
    // the junction itself — a fold-back or a junction far along the
    // curve in arc terms is never the mislocated-crossing class.
    let n_before = cands.len();
    cands.retain(|c| {
        let tm = 0.5 * (c.2 + t_end);
        window.point_at(&source.curve_3d, tm).distance_to(p) <= far_cap
    });
    if popkern_core::trace::trace_enabled("clip-far") && n_before > 0 {
        eprintln!(
            "    [clip-far] p ({:.4},{:.4},{:.4}) t_end {t_end:.4} src [{s0:.4},{s1:.4}] far_cap {far_cap:.3}: {} of {n_before} kept",
            p.x,
            p.y,
            p.z,
            cands.len()
        );
    }
    Some(cands.into_iter().map(|c| (c.1, c.2)).collect())
}
