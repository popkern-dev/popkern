//! Boundary-crossing detection for the UV clip — split from
//! `clip_uv.rs` (500-line cap ratchet). The seam-dance, period-branch
//! normalization, and torn-boundary guards live with the finder.

use popkern_core::geometry::chart::ChartDescriptor;
use popkern_core::geometry::curve::{CurveEval, CurveGeometry};
use popkern_core::geometry::intersection::intersect_curves;
use popkern_core::geometry::nurbs_curve::NurbsCurve;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};

use super::clip_crossing_graze::drop_grazes;
use super::clip_crossing_reports::{collapse_crossing_reports, resolve_contact_runs};
use super::clip_lift_crossings::extend_with_lifted_loop_crossings;
use super::clip_uv::UvCrossing;
use super::face_uv_boundary;

/// Fraction of the linear tolerance an analytic boundary pcurve's
/// NURBS form may deviate from it, in model space.
const BOUNDARY_CONVERSION_FRACTION: f64 = 1e-3;

/// Multiple of the linear tolerance within which a crossing's query
/// and boundary chart points must agree. Generous — real crossings
/// agree to the intersector's tolerance; the gate only rejects
/// grossly separated (collinear-artifact) reports.
const COINCIDENCE_BAND: f64 = 100.0;

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
    // A plane's chart has no period (a lift by one finds nothing true
    // there), so a record no boundary edge can reach crosses none:
    // decided before any edge is converted (PK-1168).
    let on_plane = brep
        .faces
        .get(face)
        .and_then(|f| f.surface.as_ref())
        .is_some_and(|s| {
            matches!(
                s.geometry,
                popkern_core::geometry::surface::SurfaceGeometry::Plane { .. }
            )
        });
    let pad = super::clip_rim_reach::RIM_REACH_TOLERANCES * tolerance.linear();
    if on_plane
        && boundary
            .iter()
            .all(|(_, pc)| super::clip_rim_reach::rim_clear_of_query(pc, curve_uv, pad))
    {
        return Some(crossings);
    }

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
    // On a rotational chart the window's crossings decide (W2c,
    // chart-service.md §6b): every whole-period lift of the record
    // against the boundary as its edges lie, one crossing per place and
    // per pass. The finder below serves the other charts.
    if super::clip_window_crossings::rotational_chart(brep, face) {
        let period = face_surface
            .and_then(|s| ChartDescriptor::evaluation_periods_of(s).0)
            .unwrap_or(std::f64::consts::TAU);
        return Some(super::clip_window_crossings::window_crossings(
            brep,
            face,
            (face_surface, &boundary_nurbs),
            curve_uv,
            period,
            tolerance,
        ));
    }
    let collect_crossings = || -> Vec<UvCrossing> {
        let working_curve = curve_uv.clone();
        let mut found = Vec::new();
        for (he_handle, boundary_pcurve) in &boundary_nurbs {
            let boundary_range = boundary_pcurve.parameter_range();
            if (boundary_range.1 - boundary_range.0).abs() < tolerance.numerical_floor() {
                continue;
            }

            let seam_use = brep.half_edges[*he_handle].mate().is_some_and(|twin_he| {
                brep.half_edges.contains_key(twin_he)
                    && brep.half_edges[twin_he].loop_ == brep.half_edges[*he_handle].loop_
            });
            if seam_use {
                continue;
            }

            let working_boundary = boundary_pcurve.clone();

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
                &CurveGeometry::Nurbs(working_boundary.clone()),
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

            // COINCIDENCE: a hit puts query and boundary at ONE chart
            // point, in the frame it was found in — collinear pieces on
            // one iso-line made the intersector report points 2 rad apart
            // (a B-side cap circle claimed a strip's rim edge at u = 3.11
            // against the edge's u = 1.07).
            let dishonest = |t1: f64, t2: f64| -> Option<&'static str> {
                use popkern_core::geometry::curve::CurveEval;
                let w = working_curve.point_at(t2);
                let b = working_boundary.point_at(t1);
                let gap = ((w.x - b.x).powi(2) + (w.y - b.y).powi(2)).sqrt();
                if gap > COINCIDENCE_BAND * tolerance.linear() {
                    return Some("apart");
                }
                // Within the linear reach in space.
                let curves = (&working_boundary, &working_curve);
                super::clip_crossing_gap::afar(face_surface, curves, (t1, t2), 0.0, tolerance)
                    .then_some("afar")
            };
            if let Ok(result) = result {
                use popkern_core::geometry::intersection::CurveIntersection;
                for hit in &result.intersections {
                    match hit {
                        CurveIntersection::Point { t1, t2, .. } => {
                            let curves = (&working_boundary, &working_curve);
                            if !super::clip_overlap_ends::crosses_on(face_surface, curves, *t2) {
                                continue;
                            }
                            if let Some(why) = dishonest(*t1, *t2) {
                                if popkern_core::trace::trace_enabled("ring-crossing") {
                                    eprintln!(
                                        "      [xing-{why}] he {he_handle:?} t2={t2:.4} rejected"
                                    );
                                }
                                continue;
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
                            if dishonest(t1_range.0, t2_range.0)
                                .or(dishonest(t1_range.1, t2_range.1))
                                .is_some()
                            {
                                continue;
                            }
                            // A COLLINEAR overlap with a boundary
                            // edge means the curve runs ALONG the
                            // face's rim for that stretch. The
                            // honest entry/exit candidates are the
                            // overlap ENDS — the old single midpoint
                            // "crossing" clipped a flush gusset
                            // cap's plate record to a fake corner
                            // 1.2 mm inside the true one and the
                            // union ring never closed.
                            found.extend(super::clip_overlap_ends::overlap_exits(
                                super::clip_overlap_ends::chart_periods(face_surface),
                                (&working_boundary, &working_curve),
                                *he_handle,
                                *t1_range,
                                *t2_range,
                            ));
                        }
                    }
                }
            }
        }
        found
    };

    let period_of = || -> f64 {
        brep.faces
            .get(face)
            .and_then(|f| f.surface.as_ref())
            .and_then(|s| ChartDescriptor::evaluation_periods_of(&s.geometry).0)
            .unwrap_or(std::f64::consts::TAU)
    };

    crossings.extend(collect_crossings());
    // Loops on multi-turn charts: crossings on every lift (PK-673),
    // each one place in space (the passes above judge their own).
    let direct = crossings.len();
    extend_with_lifted_loop_crossings(
        &mut crossings,
        brep,
        &boundary_nurbs,
        curve_uv,
        period_of(),
        tolerance,
    );
    let (lifted, surface) = ((&boundary_nurbs, curve_uv), face_surface);
    super::clip_crossing_gap::keep_one_place(
        surface,
        lifted.0,
        lifted.1,
        &mut crossings,
        direct,
        tolerance,
    );
    crossings.sort_by(|a, b| {
        a.t_curve
            .partial_cmp(&b.t_curve)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    let crossings = collapse_crossing_reports(brep, face, curve_uv, crossings, tolerance);
    let crossings = resolve_contact_runs(brep, face, curve_uv, crossings, tolerance);
    let on = (face_surface, &boundary_nurbs[..]);
    let mut crossings = drop_grazes(on, curve_uv, crossings, tolerance);
    // A closed query on the boundary at its closure point (PK-775).
    let period = None;
    let period_v = {
        brep.faces
            .get(face)
            .and_then(|f| f.surface.as_ref())
            .and_then(|s| ChartDescriptor::evaluation_periods_of(&s.geometry).1)
    };
    if let Some(c) = super::clip_lift_crossings::closure_crossing(
        brep,
        &boundary_nurbs,
        curve_uv,
        &crossings,
        (period, period_v),
    ) {
        crossings.insert(0, c);
    }
    // On a rotational chart the window's crossings decide (W2c,
    // chart-service.md §6b step 3): every whole-period lift of the record
    // against the boundary as its edges lie, one crossing per place and
    // pass. The finder above runs in shadow until it is deleted; with the
    // shadow log set its disagreements are logged.
    Some(crossings)
}
