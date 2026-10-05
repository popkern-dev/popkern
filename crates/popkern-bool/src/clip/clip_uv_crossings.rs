//! Boundary-crossing detection for the UV clip — split from
//! `clip_uv.rs` (500-line cap ratchet). The seam-dance, period-branch
//! normalization, and torn-boundary guards live with the finder.

use popkern_core::geometry::curve::{CurveEval, CurveGeometry};
use popkern_core::geometry::intersection::intersect_curves;
use popkern_core::geometry::nurbs_curve::NurbsCurve;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};

use super::clip_branch_normalize;
use super::clip_crossing_reports::collapse_crossing_reports;
use super::clip_crossings_rescue::sampled_uv_crossings;
use super::clip_lift_crossings::extend_with_lifted_loop_crossings;
use super::clip_seam_shift::{is_near_seam, shift_curve_u};
use super::clip_uv::{TWO_PI, UvCrossing};
use super::face_uv_boundary;

/// Fraction of the linear tolerance an analytic boundary pcurve's
/// NURBS form may deviate from it, in model space.
const BOUNDARY_CONVERSION_FRACTION: f64 = 1e-3;

pub(crate) fn find_uv_crossings(
    brep: &BRep,
    face: FaceHandle,
    curve_uv: &NurbsCurve,
    tolerance: &Tolerance,
) -> Option<Vec<UvCrossing>> {
    // Chart-cell route (PK-304): doubly-periodic faces compute
    // crossings from folded-containment window ends — the ±π dance
    // below cannot represent their two-axis lattice.
    if let Some(cell_crossings) = super::chart_cell::cell_crossings(brep, face, curve_uv, tolerance)
    {
        return Some(cell_crossings);
    }
    let boundary = face_uv_boundary(brep, face)?;
    // Deliberately OUTER-LOOP-ONLY: the split stage's tuned
    // hole-interaction machinery (the PK-185/186 inner-loop band
    // splits, match_inner_loop_cut) expects records that run THROUGH
    // holes at full span and partitions them itself — clipping
    // records AT hole rims here stranded three tuned
    // countersink/socket fixtures in merge (measured 2026-09-03).
    // Hole REGION semantics live in containment
    // (`point_inside_face_uv`), which is hole-aware.
    let mut crossings = Vec::new();

    let curve_range = curve_uv.parameter_range();
    if (curve_range.1 - curve_range.0).abs() < tolerance.numerical_floor() {
        return Some(crossings);
    }

    let in_periodic_range = |curve: &NurbsCurve| -> bool {
        curve
            .control_points
            .iter()
            .any(|p| p.x > std::f64::consts::PI)
    };
    // Boundary pcurves converted within BOUNDARY_CONVERSION_FRACTION of
    // the linear tolerance in model space (PK-553): the former coarse
    // 17-sample cubic displaced analytic trims by up to 1.3e-2 on
    // corpus conic trims, so crossings landed off the boundary edges
    // the split graph must meet.
    let face_surface = brep
        .faces
        .get(face)
        .and_then(|f| f.surface.as_ref())
        .map(|s| &s.geometry);
    let linear = BOUNDARY_CONVERSION_FRACTION * tolerance.linear();
    let boundary_nurbs: Vec<(popkern_core::topology::HalfEdgeHandle, NurbsCurve)> = boundary
        .iter()
        .map(|(h, pc)| {
            let nurbs = match face_surface {
                Some(surface) => pc.to_nurbs_on_surface(surface, linear).0,
                None => pc.to_nurbs(),
            };
            (*h, nurbs)
        })
        .collect();
    // The seam-shift dance is ONLY meaningful on periodic surfaces:
    // its "near seam" / "in periodic range" tests are raw u-value
    // thresholds (π/4, π) and fire spuriously on PLANAR faces whose
    // UV coordinates are plain lengths (measured: a 20mm-wide tool
    // wall got its intersection segment shifted by π and clipped to
    // ±5.78 instead of ±10, so the pocket never met the wall
    // corners and the boolean output carried degenerate stubs).
    let surface_is_periodic = brep
        .faces
        .get(face)
        .and_then(|f| f.surface.as_ref())
        .map(|s| {
            matches!(
                s.geometry,
                popkern_core::geometry::surface::SurfaceGeometry::Cylinder { .. }
                    | popkern_core::geometry::surface::SurfaceGeometry::Sphere { .. }
                    | popkern_core::geometry::surface::SurfaceGeometry::Torus { .. }
                    | popkern_core::geometry::surface::SurfaceGeometry::Cone { .. }
            )
        })
        .unwrap_or(false);
    let curve_uv = &clip_branch_normalize::normalize_query_branch(
        brep,
        face,
        &boundary_nurbs,
        curve_uv,
        surface_is_periodic,
    );
    let needs_shift = surface_is_periodic
        && (is_near_seam(curve_uv) || boundary_nurbs.iter().any(|(_, pc)| is_near_seam(pc)))
        && (in_periodic_range(curve_uv)
            || boundary_nurbs.iter().any(|(_, pc)| in_periodic_range(pc)));
    let shift = std::f64::consts::PI;

    let collect_crossings = |shifted: bool| -> Vec<UvCrossing> {
        let working_curve = if shifted {
            shift_curve_u(curve_uv, shift)
        } else {
            curve_uv.clone()
        };
        let mut found = Vec::new();
        for (he_handle, boundary_pcurve) in &boundary_nurbs {
            let boundary_range = boundary_pcurve.parameter_range();
            if (boundary_range.1 - boundary_range.0).abs() < tolerance.numerical_floor() {
                continue;
            }

            let twin_he = brep.half_edges[*he_handle].twin;
            if brep.half_edges.contains_key(twin_he)
                && brep.half_edges[twin_he].loop_ == brep.half_edges[*he_handle].loop_
            {
                continue;
            }

            let working_boundary = if shifted {
                shift_curve_u(boundary_pcurve, shift)
            } else {
                boundary_pcurve.clone()
            };
            // TORN-BOUNDARY guard (PK-304): the ±π dance's per-point
            // wrap legitimately HEALS a boundary stored the long way
            // (span shrinks), but on a short-literal boundary it
            // fabricates a bridge sweeping most of the chart (span
            // GROWS — measured: a torus quadrant pcurve u 3.142→1.571
            // became 0.000→4.712 and its phantom "crossings" claimed
            // intersection arcs for patches on the far side of the
            // surface). A shifted image that widened its u-span is a
            // tear, never a heal — skip it; the unshifted pass judges
            // that boundary honestly.
            if shifted {
                let uspan = |c: &NurbsCurve| {
                    let (lo, hi) = c
                        .control_points
                        .iter()
                        .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), p| {
                            (lo.min(p.x), hi.max(p.x))
                        });
                    hi - lo
                };
                if uspan(&working_boundary) > uspan(boundary_pcurve) + 1e-9 {
                    if popkern_core::trace::trace_enabled("ring-crossing") {
                        eprintln!("      [xing-torn-bnd] he {he_handle:?} skipped (span grew)");
                    }
                    continue;
                }
            }

            if popkern_core::trace::trace_enabled("ring-crossing") {
                use popkern_core::geometry::curve::CurveEval;
                let b0 = working_boundary.point_at(boundary_range.0);
                let b1 = working_boundary.point_at(boundary_range.1);
                let mids: Vec<String> = (1..5)
                    .map(|k| {
                        let t = boundary_range.0
                            + (boundary_range.1 - boundary_range.0) * k as f64 / 5.0;
                        let p = working_boundary.point_at(t);
                        format!("({:.3},{:.3})", p.x, p.y)
                    })
                    .collect();
                eprintln!(
                    "      [xing-try] he {he_handle:?} bnd({:.6},{:.6})->({:.6},{:.6}) mids {} brange({:.6e},{:.6e}) bcps={} qcps={}",
                    b0.x,
                    b0.y,
                    b1.x,
                    b1.y,
                    mids.join(" "),
                    boundary_range.0,
                    boundary_range.1,
                    working_boundary.control_points.len(),
                    working_curve.control_points.len()
                );
            }
            let result = intersect_curves(
                &CurveGeometry::Nurbs(working_boundary),
                boundary_range,
                &CurveGeometry::Nurbs(working_curve.clone()),
                curve_range,
                tolerance,
            );
            if popkern_core::trace::trace_enabled("ring-crossing") {
                match &result {
                    Ok(r) => {
                        eprintln!(
                            "      [xing-raw] he {he_handle:?} -> {} hit(s) qrange {curve_range:?} qcp0 {:?} qcpN {:?}",
                            r.intersections.len(),
                            working_curve.control_points.first(),
                            working_curve.control_points.last()
                        );
                    }
                    Err(e) => eprintln!("      [xing-err] he {he_handle:?} -> {e}"),
                }
            }

            if let Ok(result) = result {
                use popkern_core::geometry::intersection::CurveIntersection;
                for hit in &result.intersections {
                    match hit {
                        CurveIntersection::Point { t1, t2, .. } => {
                            // TEAR CONSISTENCY: the fold wraps each
                            // control point by rem_euclid, which
                            // TEARS a curve spanning the wrap — its
                            // NURBS evaluation then sweeps spurious
                            // bridge paths, and a crossing found on a
                            // bridge has a t-parameter pointing at
                            // unrelated geometry (measured: the
                            // gusset-web ellipse's dance crossings
                            // landed on the extended cylinder's
                            // phantom re-crossing half a period from
                            // the boundary they claimed to cross). A
                            // fold-found crossing is honest only
                            // where the folded QUERY evaluation
                            // equals the per-point fold of the
                            // original evaluation. The check is
                            // query-only: BOUNDARY pcurves may
                            // legitimately encode a seam-spanning
                            // edge "the long way" in the chart, and
                            // tearing them is precisely the healing
                            // the dance exists for; query records
                            // come from the marcher continuously
                            // unwrapped, so a torn query is always a
                            // bridge artifact.
                            if shifted {
                                use popkern_core::geometry::curve::CurveEval;
                                let untorn = |orig: &NurbsCurve, work: &NurbsCurve, t: f64| {
                                    let o = orig.point_at(t);
                                    let w = work.point_at(t);
                                    let fold = (o.x + shift).rem_euclid(TWO_PI);
                                    let du = w.x - fold;
                                    (du - TWO_PI * (du / TWO_PI).round()).abs() < 1e-6
                                        && (w.y - o.y).abs() < 1e-6 * (1.0 + o.y.abs())
                                };
                                if !untorn(curve_uv, &working_curve, *t2) {
                                    if popkern_core::trace::trace_enabled("ring-crossing") {
                                        eprintln!(
                                            "      [xing-torn] he {he_handle:?} t2={t2:.4} rejected"
                                        );
                                    }
                                    continue;
                                }
                            }
                            if popkern_core::trace::trace_enabled("ring-crossing") {
                                use popkern_core::geometry::curve::CurveEval;
                                let q = working_curve.point_at(*t2);
                                let qb = boundary_pcurve.point_at(*t1);
                                eprintln!(
                                    "      [xing] he {he_handle:?} t2={t2:.4} uv=({:.4},{:.4}) bnd-uv=({:.4},{:.4})",
                                    q.x, q.y, qb.x, qb.y
                                );
                            }
                            found.push(UvCrossing {
                                t_curve: *t2,
                                boundary_he: *he_handle,
                                t_boundary: *t1,
                            });
                        }
                        CurveIntersection::Overlap { t1_range, t2_range } => {
                            // A COLLINEAR overlap with a boundary
                            // edge means the curve runs ALONG the
                            // face's rim for that stretch. The
                            // honest entry/exit candidates are the
                            // overlap ENDS — the old single midpoint
                            // "crossing" clipped a flush gusset
                            // cap's plate record to a fake corner
                            // 1.2 mm inside the true one and the
                            // union ring never closed.
                            found.push(UvCrossing {
                                t_curve: t2_range.0,
                                boundary_he: *he_handle,
                                t_boundary: t1_range.0,
                            });
                            found.push(UvCrossing {
                                t_curve: t2_range.1,
                                boundary_he: *he_handle,
                                t_boundary: t1_range.1,
                            });
                        }
                    }
                }
            }
        }
        found
    };

    let period_of = || -> f64 {
        use popkern_core::geometry::SurfaceEval as _;
        brep.faces
            .get(face)
            .and_then(|f| f.surface.as_ref())
            .and_then(|s| s.geometry.parameter_periodicity().0)
            .unwrap_or(std::f64::consts::TAU)
    };

    // The plain (unshifted) pass ALWAYS runs: its crossings are honest
    // by construction (branch normalization above put the query on the
    // boundary's period branch). The π-shifted seam dance AUGMENTS it
    // for crossings the plain pass loses near the seam — it must never
    // REPLACE it, because the shift's per-control-point rem_euclid
    // fold TEARS any single-span boundary pcurve spanning the full
    // period (a cylinder rim stored as one degree-1 span [2π → 0]
    // folds to the degenerate [π → π] and its crossing vanishes).
    // Measured (PK-295, rim-filleted cylinder minus cuboid): every
    // wall×cylinder line record lost its bottom-rim crossing in the
    // shifted-only pass, the odd crossing count read as Outside, all
    // eight records dropped, no face was ever split, and the
    // subtraction served a five-face open box as a clean result.
    crossings.extend(collect_crossings(false));
    if needs_shift {
        let span = (curve_range.1 - curve_range.0).abs().max(1e-12);
        let dedup_eps = span * 1e-6;
        let extra: Vec<UvCrossing> = collect_crossings(true)
            .into_iter()
            .filter(|c| {
                !crossings.iter().any(|p| {
                    p.boundary_he == c.boundary_he && (p.t_curve - c.t_curve).abs() < dedup_eps
                })
            })
            .collect();
        crossings.extend(extra);
    }
    // Loops on multi-turn charts: crossings on every lift (PK-673).
    extend_with_lifted_loop_crossings(
        &mut crossings,
        brep,
        &boundary_nurbs,
        curve_uv,
        period_of(),
        tolerance,
    );
    if crossings.is_empty() && surface_is_periodic {
        // The precise curve×curve intersector can miss real
        // crossings when the query's UV extent dwarfs the trim (a
        // near-axis-parallel plane×cylinder ellipse spans ±1000 in v
        // while the face is 40 wide — measured: the gusset-web weld
        // line found no crossings here, fell through to the seam
        // dance, and the torn-polygon crossings clipped the record
        // to the phantom re-crossing). Rescue with an in/out
        // TRANSITION SCAN against the face containment predicate —
        // honest regardless of pcurve period branches. Scoped to
        // queries whose UV extent dwarfs the trim (the starvation
        // signature): firing it on comparable-scale pairs changed
        // zero-crossing containment verdicts the cone suites are
        // tuned to.
        let bbox_diag = |pts: &mut dyn Iterator<Item = (f64, f64)>| -> f64 {
            let (mut lu, mut hu, mut lv, mut hv) = (f64::MAX, f64::MIN, f64::MAX, f64::MIN);
            for (u, v) in pts {
                lu = lu.min(u);
                hu = hu.max(u);
                lv = lv.min(v);
                hv = hv.max(v);
            }
            if lu > hu {
                return 0.0;
            }
            ((hu - lu).powi(2) + (hv - lv).powi(2)).sqrt()
        };
        let q_diag = bbox_diag(&mut curve_uv.control_points.iter().map(|p| (p.x, p.y)));
        let b_diag = bbox_diag(
            &mut boundary_nurbs
                .iter()
                .flat_map(|(_, pc)| pc.control_points.iter().map(|p| (p.x, p.y))),
        );
        if q_diag.is_finite() && b_diag > 0.0 && q_diag > 8.0 * b_diag {
            crossings.extend(sampled_uv_crossings(
                brep,
                face,
                curve_uv,
                &boundary_nurbs,
                period_of(),
                tolerance,
            ));
        }
    }
    if crossings.is_empty() && surface_is_periodic {
        // The branch heuristic is a guess; when the chosen branch
        // finds NO crossings, a curve genuinely crossing the boundary
        // in the other branch would be misjudged by the containment
        // probe and clipped away whole (measured: two of the four
        // fillet-cylinder ellipses of a rotated lip union sat one
        // period from their quarter-arc charts and were discarded, so
        // the island ring never closed). Retry the other branch and
        // trust whichever finds real crossings — verified like the
        // first pass.
        crossings.extend(collect_crossings(!needs_shift));
    }

    crossings.sort_by(|a, b| {
        a.t_curve
            .partial_cmp(&b.t_curve)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    Some(collapse_crossing_reports(
        brep, face, curve_uv, crossings, tolerance,
    ))
}
