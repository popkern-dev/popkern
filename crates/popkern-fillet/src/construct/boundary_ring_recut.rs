//! Rail-spur recut for interior boundary-ring junctions (PK-234,
//! generalized to swept NURBS strips in PK-237).
//!
//! Strip insertion leaves the corner as a mitre complex: on the
//! reflex side a strip's support boundary is a straight CHORD from
//! the far rail point to the mitre tip (measured on the notched-cube
//! saddle: f10's cube-face boundary runs (1,1,0.52)→(1,0.675,0.65),
//! 2.3e-2 off the strip's own cylinder mid-span), so the pushed
//! boundary junction lands INTERIOR to the chord-bounded trim and
//! there is no rail edge to split.
//!
//! The recut restores the true rail as topology:
//!
//! 1. `mev` a rail spur from the chord's on-rail far endpoint to the
//!    junction foot (exact rail segment: the analytic axis line for
//!    cylinder/cone strips, the strip chart's exact border iso for
//!    swept NURBS strips — variable-radius saddles reach this class,
//!    PK-237);
//! 2. the caller `mef`s the arc from the foot — the strip splits
//!    into the main strip and a corner sliver (overreach band +
//!    between-station piece);
//! 3. `kef` the chord — the sliver returns to the SUPPORT face, the
//!    spur becomes the restored strip↔support rail, and the junction
//!    sits on the support face's loop for the bridge `mef`.
//!
//! Everything inside the future ring (the sliver's corner portion,
//! the mitre wedges) is absorbed into the blend face by the region
//! merge afterwards; the transient sliver-in-support ownership never
//! survives to a served result.

use super::*;
use popkern_core::geometry::Vec3;
use popkern_core::geometry::curve::{CurveEval, CurveGeometry, CurveOnEdge};
use popkern_core::geometry::nurbs_curve::NurbsCurve;

/// How far a candidate chord endpoint may sit off the rail and still
/// anchor the spur (the far endpoint is where the carve was still
/// clean, so genuine anchors measure ~1e-12; the mitre-tip endpoint
/// sits ~1e-1 off).
const SPUR_ANCHOR_ON_RAIL_BAND: f64 = 1.0e-6;

/// How far the pushed junction may sit off the selected rail model
/// (the junction carries the trim fit's wobble; measured 1e-3 class
/// on the notched-cube saddle, exact on interior-iso stations).
const JUNCTION_OFF_RAIL_BAND: f64 = 5.0e-3;

/// Result of a rail-spur recut: the junction vertex (spur tip on the
/// restored rail) plus the chord edge to `kef` after the arc `mef`
/// (the bridge phase re-derives the support face topologically).
pub(super) struct RailSpurRecut {
    pub vertex: VertexHandle,
    pub chord_edge: EdgeHandle,
}

/// The strip's rail near the junction: a straight axis line for the
/// analytic strip classes, or the chart's exact border iso curve for
/// swept NURBS strips.
enum RailModel {
    Line { dir: Vec3 },
    BorderIso { rail: NurbsCurve },
}

impl RailModel {
    fn distance_to(&self, p: &Point3, through: &Point3) -> f64 {
        match self {
            RailModel::Line { dir } => distance_to_line(p, through, dir),
            RailModel::BorderIso { rail } => {
                let t = closest_param_on_curve(rail, p);
                rail.point_at(t).distance_to(p)
            }
        }
    }
}

/// Attempt the recut on `strip_face` for a junction `point` whose
/// rail edge is missing. The chord is identified geometrically: the
/// unique strip-loop edge departing the rail (one endpoint on it,
/// one off) whose far-side face's surface contains the junction.
pub(super) fn rail_spur_recut(
    brep: &mut BRep,
    strip_face: FaceHandle,
    strip_surface: &popkern_core::geometry::SurfaceGeometry,
    point: Point3,
) -> Result<RailSpurRecut, String> {
    use popkern_core::geometry::surface::SurfaceEval;
    use popkern_core::geometry::uv_projection::project_point_to_uv;

    let rail = rail_model(strip_surface, &point)
        .ok_or_else(|| "strip surface has no rail model".to_string())?;

    // The chord: departs the rail, and its far-side face is the
    // support the junction/bridge live on.
    let mut chord: Option<(EdgeHandle, VertexHandle, f64, FaceHandle)> = None;
    for he in face_loop_half_edges(brep, strip_face)? {
        let hd = &brep.half_edges[he];
        let edge = hd.edge;
        let he0 = brep.edges[edge].half_edge;
        let v_a = brep.half_edges[he0].origin;
        let v_b = brep.half_edges[brep.half_edges[he0].twin].origin;
        let d_a = rail.distance_to(&brep.vertices[v_a].position, &point);
        let d_b = rail.distance_to(&brep.vertices[v_b].position, &point);
        let (v_on, d_on, d_off) = if d_a <= d_b {
            (v_a, d_a, d_b)
        } else {
            (v_b, d_b, d_a)
        };
        if d_on > SPUR_ANCHOR_ON_RAIL_BAND || d_off <= SPUR_ANCHOR_ON_RAIL_BAND {
            continue;
        }
        let other = brep.loops[brep.half_edges[hd.twin].loop_].face;
        let Some(surface) = brep.faces.get(other).and_then(|f| f.surface.as_ref()) else {
            continue;
        };
        let contains = project_point_to_uv(&surface.geometry, &point)
            .map(|(u, v)| {
                (&surface.geometry as &dyn SurfaceEval)
                    .point_at(u, v)
                    .distance_to(&point)
            })
            .is_some_and(|d| d <= 1.0e-6);
        if !contains {
            continue;
        }
        if let Some((existing, _, _, _)) = chord {
            return Err(format!(
                "ambiguous rail-departing chords on {strip_face:?} ({existing:?} vs {edge:?})"
            ));
        }
        chord = Some((edge, v_on, d_on, other));
    }
    let (chord_edge, v_far, _d_far, _support_face) = chord.ok_or_else(|| {
        format!(
            "no rail-departing chord with a junction-bearing far side on {strip_face:?} \
             — the mitre complex here is not the chord class"
        )
    })?;

    // Junction foot: exact projection onto the rail model (the
    // pushed junction carries the trim fit's wobble; the restored
    // rail is the clean geometry through the far anchor).
    let p_far = brep.vertices[v_far].position;
    let (foot, spur_curve) = match &rail {
        RailModel::Line { dir } => {
            let foot = project_to_line(&point, &p_far, dir);
            (foot, None)
        }
        RailModel::BorderIso { rail } => {
            let t_foot = closest_param_on_curve(rail, &point);
            let foot = rail.point_at(t_foot);
            if foot.distance_to(&point) > JUNCTION_OFF_RAIL_BAND {
                return Err(format!(
                    "junction sits {:.3e} off the strip border iso (band {JUNCTION_OFF_RAIL_BAND:.0e})",
                    foot.distance_to(&point)
                ));
            }
            let t_far = closest_param_on_curve(rail, &p_far);
            (foot, Some((rail.clone(), t_far, t_foot)))
        }
    };

    let he_far = half_edge_at_vertex_on_face(brep, strip_face, v_far)?;
    let mev = brep
        .mev(he_far, foot)
        .map_err(|e| format!("mev rail spur: {e}"))?;
    if let Some(edge) = brep.edges.get_mut(mev.edge) {
        edge.curve = match spur_curve {
            None => CurveOnEdge::line_segment(p_far, foot, 1.0e-12),
            Some((rail_curve, t_far, t_foot)) => Some(CurveOnEdge {
                geometry: CurveGeometry::Nurbs(rail_curve),
                t_start: t_far,
                t_end: t_foot,
            }),
        };
    }
    Ok(RailSpurRecut {
        vertex: mev.vertex,
        chord_edge,
    })
}

/// The strip's rail model near `point`: the analytic axis line for
/// cylinder/cone strips; for swept NURBS strips, the chart border
/// iso curve nearest the junction (exact on the surface by
/// identity). `None` only for surface families with no rail notion.
fn rail_model(
    surface: &popkern_core::geometry::SurfaceGeometry,
    point: &Point3,
) -> Option<RailModel> {
    use popkern_core::geometry::SurfaceGeometry;
    use popkern_core::geometry::surface::SurfaceEval;
    match surface {
        SurfaceGeometry::Cylinder { axis, .. } | SurfaceGeometry::Cone { axis, .. } => {
            axis.normalized().map(|dir| RailModel::Line { dir })
        }
        SurfaceGeometry::Nurbs(n) => {
            let ((u_lo, u_hi), (v_lo, v_hi)) = n.parameter_range();
            let borders = [
                n.iso_curve_at_u(u_lo),
                n.iso_curve_at_u(u_hi),
                n.iso_curve_at_v(v_lo),
                n.iso_curve_at_v(v_hi),
            ];
            borders
                .into_iter()
                .map(|rail| {
                    let t = closest_param_on_curve(&rail, point);
                    let d = rail.point_at(t).distance_to(point);
                    (d, rail)
                })
                .min_by(|(a, _), (b, _)| a.total_cmp(b))
                .map(|(_, rail)| RailModel::BorderIso { rail })
        }
        _ => None,
    }
}

/// Grid + window-refined closest parameter on a curve (65-sample
/// coarse pass, 24 shrinking-window rounds — the boundary-ring
/// convention for point-to-curve queries).
fn closest_param_on_curve(curve: &NurbsCurve, p: &Point3) -> f64 {
    const S: usize = 64;
    let (t0, t1) = curve.parameter_range();
    let mut best = (f64::INFINITY, t0);
    for i in 0..=S {
        let t = t0 + (t1 - t0) * i as f64 / S as f64;
        let d = curve.point_at(t).distance_to(p);
        if d < best.0 {
            best = (d, t);
        }
    }
    let mut w = (t1 - t0).abs() / S as f64;
    for _ in 0..24 {
        let lo = (best.1 - w).max(t0.min(t1));
        let hi = (best.1 + w).min(t0.max(t1));
        for i in 0..=16 {
            let t = lo + (hi - lo) * i as f64 / 16.0;
            let d = curve.point_at(t).distance_to(p);
            if d < best.0 {
                best = (d, t);
            }
        }
        w *= 0.35;
    }
    best.1
}

fn distance_to_line(p: &Point3, on_line: &Point3, dir: &Vec3) -> f64 {
    let w = *p - *on_line;
    let along = w.dot(dir);
    (w - *dir * along).length()
}

fn project_to_line(p: &Point3, on_line: &Point3, dir: &Vec3) -> Point3 {
    let w = *p - *on_line;
    let along = w.dot(dir);
    *on_line + *dir * along
}
