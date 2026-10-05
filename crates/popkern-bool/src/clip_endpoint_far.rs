//! FAR-tier clip endpoint machinery: exact-junction search along a
//! record's source curve, dense curve projection, and the
//! shortening-only re-trim. Split from `clip_endpoint_refine.rs`
//! under the 500-line structural budget; the rationale comments live
//! with the code.

use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::curve_on_surface::CurveOnSurface;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::uv_projection::project_point_to_uv_near;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};

use super::clip_endpoint_refine::{endpoint, trim_boundary_edges};
use crate::closed_source_window::{SourceWindow, fold_toward};
use popkern_core::geometry::surface::SurfaceEval as _;

pub(crate) use crate::clip_far_junction::far_junction;
pub(crate) use crate::curve_projection::project_onto_curve_dense;

/// Re-trim `piece` along its SOURCE curve so the given end lands at
/// the exact junction: all three curves (3D + both uv charts) are
/// re-sampled from the source over the corrected parameter window and
/// re-fit, with the end sample pinned to the junction. Returns false
/// (piece untouched) when the source parameterization does not match
/// the piece or the correction would degenerate the span.
#[allow(clippy::too_many_arguments)]
pub(crate) fn retrim_piece_to_junction(
    piece: &mut CurveOnSurface,
    source: &CurveOnSurface,
    end: bool,
    t_star: f64,
    j: Point3,
    brep_a: &BRep,
    face_a: FaceHandle,
    surf_a: &SurfaceGeometry,
    brep_b: &BRep,
    face_b: FaceHandle,
    surf_b: &SurfaceGeometry,
    tolerance: &Tolerance,
) -> bool {
    let fit_gate = tolerance.linear() * 2000.0;
    let far_trace = popkern_core::trace::trace_enabled("clip-far");
    // Work entirely in SOURCE parameters. The piece's own parameter
    // domain is not trustworthy here: an earlier endpoint refinement
    // re-fits the piece in its own coordinates (measured on the
    // angled tap: a seam-split piece carried domain [0,1] against a
    // source spanning [-4.49,-1.89], and every window computed from
    // piece parameters was garbage). Both piece endpoints are instead
    // PROJECTED onto the source; the projections certify that the
    // piece lies on the source curve.
    let this_p = endpoint(piece, end);
    let other_p = endpoint(piece, !end);
    let (Some(t_this), Some(t_other)) = (
        project_onto_curve_dense(&source.curve_3d, &this_p),
        project_onto_curve_dense(&source.curve_3d, &other_p),
    ) else {
        return false;
    };
    let on_source = |t: f64, q: &Point3| source.curve_3d.point_at(t).distance_to(q) < fit_gate;
    let mid_p = {
        let (p0, p1) = piece.curve_3d.parameter_range();
        piece.curve_3d.point_at(0.5 * (p0 + p1))
    };
    let t_mid = project_onto_curve_dense(&source.curve_3d, &mid_p);
    let Some(t_mid) = t_mid.filter(|&t| on_source(t, &mid_p)) else {
        if far_trace {
            eprintln!("    [clip-far] retrim refused: piece does not lie on the source curve");
        }
        return false;
    };
    if !on_source(t_this, &this_p) || !on_source(t_other, &other_p) {
        if far_trace {
            eprintln!("    [clip-far] retrim refused: piece does not lie on the source curve");
        }
        return false;
    }
    // A CLOSED source's window may run through its parameter seam
    // (see `closed_source_window`): unwrap the ends along the piece's
    // own traversal and lift the junction next to this end. Open
    // sources and seam-free pieces keep their raw parameters.
    let window = SourceWindow::of(source);
    let (t_this, t_other) = window.unwrap_piece(t_this, t_other, t_mid);
    let t_star = window.lift_near(t_star, t_this);
    // Corrected window in source parameters, oriented the way the
    // piece traverses: the far end stays, this end moves to the
    // junction's parameter. SHORTENING ONLY — pulling an overshot end
    // back onto a junction the clip ran past can never invent
    // geometry, while EXTENDING past the clipped end can run over a
    // legitimate mid-face termination (missing in-face stretches are
    // the pierce-authority windows' job).
    let (n0, n1) = if end {
        (t_other, t_star)
    } else {
        (t_star, t_other)
    };
    let dir_old = if end {
        t_this - t_other
    } else {
        t_other - t_this
    };
    let span_old = dir_old.abs();
    let extending = (t_star - t_other).abs() >= (t_this - t_other).abs();
    if (n1 - n0).abs() < span_old * 0.05
        || !(n1 - n0).is_finite()
        || (n1 - n0).signum() != dir_old.signum()
    {
        if far_trace {
            eprintln!(
                "    [clip-far] retrim refused: window [{n0:.4},{n1:.4}] vs [{t_other:.4},{t_this:.4}] degenerate/reversed"
            );
        }
        return false;
    }
    // EXTENSION tier (PK-472): a graze crossing can also be
    // mislocated SHORT of the true junction — the m6 nut's mirrored
    // mouth chamfer's crest×cone circle stopped 3.7e-2 before the
    // helix junction; the shortening-only rule refused, the far-tier
    // weld hauled the NODE while the curve stayed short, and the
    // 3.7e-2 curve-vs-node lie stranded the walk (5 sides).
    // Extension INVENTS geometry, so the bar is strict: the exact
    // junction already lies ON the source curve (caller-certified),
    // and EVERY added-span sample must RIDE a face boundary within
    // the graze slope gate — the added run is the mislocated
    // boundary hug itself, never an in-face stretch (a genuine
    // mid-face termination has no boundary alongside).
    if extending && !source.is_closed {
        // Extension on an OPEN source is scoped to the stop-short
        // MICRO class (PK-480): the junction must lie ON the source
        // curve at machine precision (a genuine pass-through — the
        // deep-plunged tap's cone×flank arcs stop 2.4e-3 short of
        // their crest triple point with dj ≈ 4e-8, while a fit-slop
        // coincidence reads 1e-4+), and the added arc stays under the
        // near-tangent junction-slide cap. The chamfered rod's cap
        // records that regressed the unscoped extension needed ~3e-1
        // runs alongside boundaries and stay refused by the cap; the
        // boundary-ride certificate below still applies to whatever
        // passes.
        let source_chord = {
            let (s0, s1) = source.curve_3d.parameter_range();
            source
                .curve_3d
                .point_at(s0)
                .distance_to(&source.curve_3d.point_at(s1))
                .max(
                    source
                        .curve_3d
                        .point_at(s0)
                        .distance_to(&source.curve_3d.point_at(0.5 * (s0 + s1))),
                )
        };
        let dj = window.point_at(&source.curve_3d, t_star).distance_to(&j);
        let added = j.distance_to(&this_p);
        // Both gates are bounded by the source's own scale so a
        // coarse-tolerance corpus (multicurve cubes at linear = 1e-3)
        // cannot blow them into unit-sized licenses.
        let dj_gate = (tolerance.linear() * 100.0).min(source_chord * 1.0e-3);
        let added_gate = (tolerance.linear() * 10_000.0).min(source_chord * 0.02);
        if dj > dj_gate || added > added_gate {
            if far_trace {
                eprintln!(
                    "    [clip-far] retrim refused: extension on open-source record \
                     ([{n0:.4},{n1:.4}] vs [{t_other:.4},{t_this:.4}], dj {dj:.3e}, added {added:.3e})"
                );
            }
            return false;
        }
    }
    if extending {
        let added_lo = t_this.min(t_star);
        let added_hi = t_this.max(t_star);
        let boundary_curves: Vec<_> = [(brep_a, face_a), (brep_b, face_b)]
            .into_iter()
            .flat_map(|(brep, face)| {
                trim_boundary_edges(brep, face)
                    .into_iter()
                    .filter_map(move |eh| brep.edges.get(eh)?.curve.clone())
            })
            .collect();
        let mut all_ride = true;
        for k in 1..=3 {
            let t = added_lo + (added_hi - added_lo) * k as f64 / 4.0;
            let q = window.point_at(&source.curve_3d, t);
            // The graze deviation grows FROM the junction (zero
            // there, sin(theta) x arc elsewhere) — the gate scales
            // with distance-to-junction alone (measured on the
            // tapped hole's crest circle: slope 0.065 at every
            // sample; a min() with distance-to-old-end strangled
            // the sample nearest the old end, whose deviation is
            // largest by construction).
            let clearance = (tolerance.linear() * 1000.0).max(0.15 * q.distance_to(&j));
            let mut best = f64::MAX;
            for c in &boundary_curves {
                let Some(tc) = (&c.geometry as &dyn CurveEval).closest_point(&q) else {
                    continue;
                };
                let tcc = tc.clamp(c.t_start.min(c.t_end), c.t_start.max(c.t_end));
                best = best.min(c.geometry.point_at(tcc).distance_to(&q));
            }
            if far_trace {
                eprintln!(
                    "    [clip-far] extension sample ({:.4},{:.4},{:.4}) boundary d {best:.3e} clearance {clearance:.3e}",
                    q.x, q.y, q.z
                );
            }
            if best > clearance {
                all_ride = false;
                break;
            }
        }
        if !all_ride {
            if far_trace {
                eprintln!(
                    "    [clip-far] retrim refused: extension span leaves the boundary ([{n0:.4},{n1:.4}] vs [{t_other:.4},{t_this:.4}])"
                );
            }
            return false;
        }
    }
    // The stretch being CUT must not read plainly in-face on both
    // charts — a shortening that removes genuinely-kept geometry is a
    // mis-identified junction, not an overshoot (measured: the
    // notched-cube saddle's fillet records shortened onto junctions
    // their curves pass through, and the subtraction went
    // non-manifold). A sample RIDING a face boundary ABSTAINS: UV
    // containment within boundary slop is a coin flip, and the
    // retrim exists exactly to fix mislocated boundary runs
    // (measured, PK-472: the m6 nut's crest×cone graze circle
    // overshot its chain junction by 4.2e-2; every removed-span
    // sample hugged the crest face's boundary helix, the coin flips
    // read in-face, and the refused retrim left a backward-facing
    // tip the gap closure could only mis-bridge). Distance to the
    // faces' boundary edges in 3D is the reliability meter — the
    // same instrument `intervals_clear_of_boundary` uses.
    {
        let boundary_curves: Vec<_> = [(brep_a, face_a), (brep_b, face_b)]
            .into_iter()
            .flat_map(|(brep, face)| {
                trim_boundary_edges(brep, face)
                    .into_iter()
                    .filter_map(move |eh| brep.edges.get(eh)?.curve.clone())
            })
            .collect();
        // A graze mislocation's samples diverge from the boundary
        // LINEARLY at slope sin(crossing angle) ≪ 1 (measured on the
        // m6 nut's crest×cone circle: 6.2e-4 / 1.2e-3 / 1.9e-3 over
        // a 4.2e-2 correction — slope 0.06, a 3.4° graze). The gate
        // therefore scales with the sample's own distance from the
        // junction: shallower than ~8.5° abstains, a transversal
        // crossing's samples (slope ~1) keep full veto power.
        let rides_boundary = |q: &Point3| -> bool {
            let clearance = (tolerance.linear() * 1000.0).max(0.15 * q.distance_to(&j));
            let mut best = f64::MAX;
            for c in &boundary_curves {
                let Some(t) = (&c.geometry as &dyn CurveEval).closest_point(q) else {
                    continue;
                };
                let tc = t.clamp(c.t_start.min(c.t_end), c.t_start.max(c.t_end));
                best = best.min(c.geometry.point_at(tc).distance_to(q));
            }
            if far_trace {
                eprintln!(
                    "    [clip-far] removed-span sample ({:.4},{:.4},{:.4}) boundary d {best:.3e} clearance {clearance:.3e}",
                    q.x, q.y, q.z
                );
            }
            best < clearance
        };
        let mut inside_votes = 0;
        for k in 1..=3 {
            let t = t_star + (t_this - t_star) * k as f64 / 4.0;
            if rides_boundary(&window.point_at(&source.curve_3d, t)) {
                continue;
            }
            let in_a = crate::split_inner_loop::point_inside_face_region_uv(
                brep_a,
                face_a,
                window.point_at(&source.curve_uv1, t),
                &brep_a.tolerance,
            );
            let in_b = crate::split_inner_loop::point_inside_face_region_uv(
                brep_b,
                face_b,
                window.point_at(&source.curve_uv2, t),
                &brep_b.tolerance,
            );
            if in_a && in_b {
                inside_votes += 1;
            }
        }
        if inside_votes >= 2 {
            if far_trace {
                eprintln!("    [clip-far] retrim refused: removed span reads in-face");
            }
            return false;
        }
    }
    // A window through a CLOSED source's seam re-trims exactly.
    if window.exact_seam_retrim(
        piece,
        source,
        (n0, n1),
        end,
        j,
        other_p,
        (surf_a, surf_b),
        tolerance,
    ) {
        return true;
    }
    let (idx_this, idx_other) = if end {
        (16usize, 0usize)
    } else {
        (0usize, 16usize)
    };
    let n = 17usize;
    let params: Vec<f64> = (0..n)
        .map(|k| n0 + (n1 - n0) * k as f64 / (n - 1) as f64)
        .collect();
    // The refit keeps the PIECE's endpoints authoritative: the far
    // end stays exactly where it is (it may itself sit on a refined
    // junction), the near end lands exactly on the junction.
    let mut pts: Vec<Point3> = params
        .iter()
        .map(|&t| window.point_at(&source.curve_3d, t))
        .collect();
    pts[idx_this] = j;
    pts[idx_other] = other_p;
    let Ok(fitted_3d) = popkern_core::geometry::fitting::fit_curve(&pts, &params, 3, n) else {
        return false;
    };
    let mut fitted_uv = Vec::with_capacity(2);
    for (src_uv, surf) in [(&source.curve_uv1, surf_a), (&source.curve_uv2, surf_b)] {
        let mut uv_pts: Vec<Point3> = params.iter().map(|&t| window.point_at(src_uv, t)).collect();
        window.unwrap_chart_run(&mut uv_pts, surf);
        // Analytic projections return the canonical period branch;
        // pin each end on the branch of the sample it replaces.
        let (pu, pv) = surf.parameter_periodicity();
        for (idx, q) in [(idx_this, &j), (idx_other, &other_p)] {
            let near = uv_pts[idx];
            if let Some((nu, nv)) = project_point_to_uv_near(surf, q, near.x, near.y) {
                uv_pts[idx] = Point3::new(
                    fold_toward(nu, near.x, pu),
                    fold_toward(nv, near.y, pv),
                    0.0,
                );
            }
        }
        let Ok(f) = popkern_core::geometry::fitting::fit_curve(&uv_pts, &params, 3, n) else {
            return false;
        };
        fitted_uv.push(f);
    }
    piece.curve_3d = fitted_3d;
    piece.curve_uv2 = fitted_uv.pop().unwrap_or_else(|| piece.curve_uv2.clone());
    piece.curve_uv1 = fitted_uv.pop().unwrap_or_else(|| piece.curve_uv1.clone());
    true
}
