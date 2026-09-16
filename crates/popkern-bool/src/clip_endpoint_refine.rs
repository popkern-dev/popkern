//! Clipped-record endpoint refinement onto face-boundary junctions.
//!
//! A clipped record piece's endpoints are evaluated on the FITTED
//! record curve at the crossing parameter and carry its lateral slop
//! (~1e-4..1e-3 measured). The boundary edge the record crosses is
//! split at the edge-face hit's exact Newton solve — so the merge
//! graph minted TWO nodes for one junction (1.79e-4 apart on the
//! chamfered-rod fixture), the record dead-ended at its own node,
//! the chart walk's graze filter deleted it, and the face never
//! partitioned (PK-460: every headed fastener declined; a node-level
//! weld alone then failed shell-gap validation because the CURVES
//! still ended at the unrefined position).
//!
//! The fix at the source: when a clipped endpoint lies within the
//! record-fit slop band of a boundary EDGE of either face, re-solve
//! the junction exactly — the boundary edge's 3D curve against the
//! OTHER face's surface — and move the endpoint (3D and both UV
//! charts) onto that solution. Both provenances then mint the same
//! graph node within snap, and the curves agree with it.

use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::{CurveEval, CurveGeometry};
use popkern_core::geometry::curve_on_surface::CurveOnSurface;
use popkern_core::geometry::intersection::{SurfaceCurveIntersection, intersect_surface_curve};
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::uv_projection::project_point_to_uv_near;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};

use crate::clip_endpoint_far::{far_junction, project_onto_curve_dense, retrim_piece_to_junction};

/// Refine both endpoints of `piece` onto exact boundary junctions
/// where they lie within the record-fit slop band of a face-boundary
/// edge. No-op for endpoints away from every boundary.
#[allow(clippy::too_many_arguments)]
pub(crate) fn refine_endpoints_onto_boundaries(
    brep_a: &BRep,
    face_a: FaceHandle,
    surf_a: &SurfaceGeometry,
    windows_a: ((f64, f64), (f64, f64)),
    brep_b: &BRep,
    face_b: FaceHandle,
    surf_b: &SurfaceGeometry,
    windows_b: ((f64, f64), (f64, f64)),
    piece: &mut CurveOnSurface,
    source: &CurveOnSurface,
    tolerance: &Tolerance,
) {
    refine_selected_endpoints_onto_boundaries(
        brep_a,
        face_a,
        surf_a,
        windows_a,
        brep_b,
        face_b,
        surf_b,
        windows_b,
        piece,
        source,
        tolerance,
        (true, true),
    );
}

/// [`refine_endpoints_onto_boundaries`] restricted to selected ends —
/// `(start, end)`. A trim-residue piece refines ONLY the end it
/// shares with its kept sibling's cut bound (PK-472): its free end is
/// a mid-face slit tip that legitimately stands OFF the boundaries,
/// and the near-tier drag glued such tips onto junctions (measured:
/// the m6 mouth chamfer re-stranded 4 sides).
#[allow(clippy::too_many_arguments)]
pub(crate) fn refine_selected_endpoints_onto_boundaries(
    brep_a: &BRep,
    face_a: FaceHandle,
    surf_a: &SurfaceGeometry,
    windows_a: ((f64, f64), (f64, f64)),
    brep_b: &BRep,
    face_b: FaceHandle,
    surf_b: &SurfaceGeometry,
    windows_b: ((f64, f64), (f64, f64)),
    piece: &mut CurveOnSurface,
    source: &CurveOnSurface,
    tolerance: &Tolerance,
    ends: (bool, bool),
) {
    // The record-fit slop class measures 2e-5..1.8e-4; the band stays
    // just above it. A generous band regressed cube corpora by
    // dragging record ends that legitimately terminate NEAR (not on)
    // a boundary — a mid-face T-junction 1e-3 from an edge — onto
    // the boundary junction (measured: two_cubes_subtracted strands).
    // Tolerance-relative alone is not enough: coarse-tolerance
    // corpora (the multicurve cubes run at linear = 1e-3) turned
    // 250× linear into 0.25 and every mid-face junction within a
    // quarter-unit of an edge was dragged onto it. Each candidate
    // edge caps the band at 2% of its own chord — a refinement that
    // moves an endpoint a large fraction of the edge it lands on is
    // never the fit-slop class.
    let band = tolerance.linear() * 250.0;
    for end in [false, true] {
        if !(if end { ends.1 } else { ends.0 }) {
            continue;
        }
        let p = endpoint(piece, end);
        // Outward unit tangent of the record at this endpoint: an
        // endpoint move ALONG it is a parameter re-trim of the same
        // curve, not a lateral drag — a shallow record × boundary
        // crossing amplifies the lateral fit slop into an along-curve
        // junction miss by 1/sin(crossing angle) (measured, PK-464:
        // the mouth-chamfer circle met the nut's crest-band helix at
        // slope 0.16 and the window end sat 6.2e-4 along the record
        // from the true junction — beyond the lateral band). The
        // lateral component keeps the tuned band; the total move gets
        // 8× headroom.
        let tan_out = endpoint_outward_tangent(piece, end);
        let refined = refine_against_face(
            brep_a, face_a, surf_b, windows_b, &p, band, &tan_out, tolerance,
        )
        .or_else(|| {
            refine_against_face(
                brep_b, face_b, surf_a, windows_a, &p, band, &tan_out, tolerance,
            )
        });
        if let Some(r) = refined
            && r.distance_to(&p) > tolerance.linear()
            && r.distance_to(&p) < band * 64.0
        {
            if popkern_core::trace::trace_enabled("clip") {
                eprintln!(
                    "    [clip-refine] end {} ({:.6},{:.6},{:.6}) -> ({:.6},{:.6},{:.6}) d={:.2e}",
                    end as u8,
                    p.x,
                    p.y,
                    p.z,
                    r.x,
                    r.y,
                    r.z,
                    r.distance_to(&p)
                );
            }
            set_endpoint_3d_uv(piece, end, r, surf_a, surf_b, tolerance);
            continue;
        }
        // FAR RE-TRIM: a record crossing a boundary NEAR-TANGENT puts
        // the clip's crossing parameter far along the curve from the
        // true junction — beyond any lateral-slop band (measured,
        // PK-465: the angled tap's thread exit crosses the mouth rim
        // at helix pitch angle; endpoints landed 0.04–0.54 from
        // existing exact junctions, in BOTH directions — stopped
        // short and overshot past the rim). The tangent-ray tests are
        // meaningless at this scale, so the certificate changes: the
        // exact edge×surface junction must lie ON the record's SOURCE
        // (pre-clip) curve within the fit-slop gate, and the endpoint
        // is then RE-TRIMMED along the source to the junction's own
        // parameter — a parameter correction of the same clip, never
        // a lateral drag, so a mid-face T-junction (off the curve)
        // can never be captured.
        // This endpoint's parameter ON THE SOURCE (the piece's own
        // parameter domain may have been re-fit by an earlier
        // refinement and no longer align with the source).
        let Some(t_on_source) = project_onto_curve_dense(&source.curve_3d, &p) else {
            continue;
        };
        if source.curve_3d.point_at(t_on_source).distance_to(&p) > tolerance.linear() * 2000.0 {
            continue;
        }
        if let Some(cands) = far_junction(
            brep_a,
            face_a,
            surf_a,
            windows_a,
            brep_b,
            face_b,
            surf_b,
            windows_b,
            &p,
            source,
            t_on_source,
            tolerance,
        ) {
            for (j, t_star) in cands {
                if retrim_piece_to_junction(
                    piece, source, end, t_star, j, brep_a, face_a, surf_a, brep_b, face_b, surf_b,
                    tolerance,
                ) {
                    if popkern_core::trace::trace_enabled("clip") {
                        eprintln!(
                            "    [clip-retrim] end {} ({:.4},{:.4},{:.4}) -> junction ({:.4},{:.4},{:.4}) t* {t_star:.4}",
                            end as u8, p.x, p.y, p.z, j.x, j.y, j.z
                        );
                    }
                    break;
                }
            }
        }
    }
}

pub(crate) fn endpoint(piece: &CurveOnSurface, end: bool) -> Point3 {
    let (t0, t1) = piece.curve_3d.parameter_range();
    piece.curve_3d.point_at(if end { t1 } else { t0 })
}

/// Unit tangent at the piece's endpoint pointing OUTWARD (past the
/// end), by an inside finite difference. Zero vector on degenerate
/// spans — the caller's lateral test then reduces to the plain
/// distance gate.
fn endpoint_outward_tangent(piece: &CurveOnSurface, end: bool) -> popkern_core::geometry::Vec3 {
    let (t0, t1) = piece.curve_3d.parameter_range();
    let span = t1 - t0;
    let (t_at, t_in) = if end {
        (t1, t1 - span * 1e-3)
    } else {
        (t0, t0 + span * 1e-3)
    };
    let d = piece.curve_3d.point_at(t_at) - piece.curve_3d.point_at(t_in);
    d.normalized()
        .unwrap_or(popkern_core::geometry::Vec3::zero())
}

/// Exact junction near `p`: the nearest boundary edge of `face`
/// (within the slop band) intersected with `other_surface`.
pub(crate) use crate::clip_boundary_edges::{
    boundary_edges, boundary_edges_outer, trim_boundary_edges,
};

#[allow(clippy::too_many_arguments)]
fn refine_against_face(
    brep: &BRep,
    face: FaceHandle,
    other_surface: &SurfaceGeometry,
    other_windows: ((f64, f64), (f64, f64)),
    p: &Point3,
    band: f64,
    tan_out: &popkern_core::geometry::Vec3,
    tolerance: &Tolerance,
) -> Option<Point3> {
    let mut best: Option<Point3> = None;
    for eh in boundary_edges_outer(brep, face) {
        if let Some(edge) = brep.edges.get(eh)
            && let Some(c) = edge.curve.as_ref()
        {
            // Nearness gate by TRUE distance to the edge curve —
            // sparse samples of a rim circle sit up to a quarter
            // circumference apart and a point ON the rim between
            // them reads far.
            let (t0, t1) = (c.t_start, c.t_end);
            let chord = c
                .geometry
                .point_at(t0)
                .distance_to(&c.geometry.point_at(t1))
                .max(
                    c.geometry
                        .point_at(t0)
                        .distance_to(&c.geometry.point_at(0.5 * (t0 + t1))),
                );
            let band_eff = band.min(0.02 * chord);
            // The near gate matches the widened LATERAL acceptance
            // (4x): the fit slop that displaces the endpoint from the
            // junction displaces it from the boundary edge by the
            // same lateral amount.
            let near = (&c.geometry as &dyn popkern_core::geometry::curve::CurveEval)
                .closest_point(p)
                .map(|t| {
                    let tc = t.clamp(t0.min(t1), t0.max(t1));
                    c.geometry.point_at(tc).distance_to(p) < band_eff * 8.0
                })
                .unwrap_or(false);
            if near {
                let range = (t0.min(t1), t0.max(t1));
                if let Ok(result) = intersect_surface_curve(
                    other_surface,
                    other_windows.0,
                    other_windows.1,
                    &clone_curve_geometry(&c.geometry),
                    range,
                    tolerance,
                ) {
                    for isect in &result.intersections {
                        if let SurfaceCurveIntersection::Point { point, t, .. } = isect {
                            let d = point.distance_to(p);
                            // Tangent-decomposed acceptance: within
                            // the tuned band as before, or laterally
                            // within it with the along-record move
                            // allowed the measured amplification —
                            // the lateral fit slop maps to an
                            // along-curve junction miss by
                            // 1/sin(crossing angle) between the
                            // record and the boundary edge (capped at
                            // 16×; measured 2.9e-3 at a 0.16-slope
                            // helix × circle crossing, PK-464).
                            let mv = *point - *p;
                            let along = mv.dot(tan_out);
                            let lateral = (mv - *tan_out * along).length();
                            let span = c.t_end - c.t_start;
                            let edge_tan = (c.geometry.point_at(t + span.abs() * 1e-4)
                                - c.geometry.point_at(t - span.abs() * 1e-4))
                            .normalized();
                            let sin_angle =
                                edge_tan.map(|et| et.cross(tan_out).length()).unwrap_or(1.0);
                            let allowed = band_eff / sin_angle.max(1.0 / 16.0);
                            let accept = d < band_eff || (lateral < band_eff && d < allowed);
                            if accept
                                && best.is_none_or(|b| point.distance_to(p) < b.distance_to(p))
                            {
                                best = Some(*point);
                            }
                        }
                    }
                }
            }
        }
    }
    best
}

pub(crate) fn clone_curve_geometry(g: &CurveGeometry) -> CurveGeometry {
    g.clone()
}

/// Move one endpoint of the piece: the 3D end control point, and each
/// UV chart's end control point re-projected near its current value
/// (clamped NURBS — the endpoint IS a control point).
///
/// A LARGE move (beyond 50× linear) re-authors the whole piece by
/// interpolation through on-curve samples with the corrected end
/// sample instead: the bare end-control-point move shifts every
/// point of the final span by the endpoint's delta VECTOR, whose
/// direction rotates off the local tangent along the span — a
/// 6.2e-4 tangential end correction on a 0.44-rad circle arc bent
/// the span 4.0e-5 off both carrier surfaces and the serve failed
/// the curve-surface gate (PK-464). Interpolated samples stay on
/// the fitted curve, so the correction is confined to the last
/// segment, whose ends are both on-surface.
pub(crate) fn set_endpoint_3d_uv(
    piece: &mut CurveOnSurface,
    end: bool,
    r: Point3,
    surf_a: &SurfaceGeometry,
    surf_b: &SurfaceGeometry,
    tolerance: &Tolerance,
) {
    let set_cp = |c: &mut popkern_core::geometry::nurbs_curve::NurbsCurve, p: Point3| {
        if end {
            if let Some(cp) = c.control_points.last_mut() {
                *cp = p;
            }
        } else if let Some(cp) = c.control_points.first_mut() {
            *cp = p;
        }
    };
    // Branch-honest projection: `project_point_to_uv_near` ignores
    // the hint on analytic surfaces and returns the CANONICAL period
    // branch, while the piece's uv pcurve may legitimately live on an
    // unwrapped branch (a seam-split segment of a closed torus lens
    // loop carries v in [2π, 4π] on its second-lap stretch). Pinning
    // a canonical-branch endpoint onto such a curve makes a
    // mixed-branch pcurve whose first segment jumps a whole period —
    // the chart walk's junction direction then points a period away
    // and the wedge matching U-turns the chain, double-claiming the
    // twin sides (measured, PK-484: the rotated torus × torus union's
    // B-flank walked both sides of two lens-chain edges and the
    // assembly stranded 4 record sides). Unwrap each projected
    // coordinate onto the branch of the uv sample it replaces.
    let unwrap_toward = |surf: &SurfaceGeometry, p: (f64, f64), near: (f64, f64)| -> (f64, f64) {
        use popkern_core::geometry::surface::SurfaceEval as _;
        let (pu, pv) = surf.parameter_periodicity();
        let fold = |x: f64, anchor: f64, period: Option<f64>| match period {
            Some(pp) if pp.is_finite() && pp > 0.0 => x + ((anchor - x) / pp).round() * pp,
            _ => x,
        };
        (fold(p.0, near.0, pu), fold(p.1, near.1, pv))
    };
    let move_d = endpoint(piece, end).distance_to(&r);
    if move_d > tolerance.linear() * 50.0 {
        let n = 17usize;
        let (t0, t1) = piece.curve_3d.parameter_range();
        let mut pts: Vec<Point3> = (0..n)
            .map(|k| {
                piece
                    .curve_3d
                    .point_at(t0 + (t1 - t0) * k as f64 / (n - 1) as f64)
            })
            .collect();
        let idx = if end { n - 1 } else { 0 };
        pts[idx] = r;
        // Keep the piece's parameter domain: downstream bookkeeping
        // (window t ranges, node parameters, closure flags) reads it.
        let params: Vec<f64> = (0..n)
            .map(|k| t0 + (t1 - t0) * k as f64 / (n - 1) as f64)
            .collect();
        if let Ok(fitted) = popkern_core::geometry::fitting::fit_curve(&pts, &params, 3, n) {
            piece.curve_3d = fitted;
        } else {
            set_cp(&mut piece.curve_3d, r);
        }
        for (uv_curve, surf) in [
            (&mut piece.curve_uv1, surf_a),
            (&mut piece.curve_uv2, surf_b),
        ] {
            let (u0, u1) = uv_curve.parameter_range();
            let mut uv_pts: Vec<Point3> = (0..n)
                .map(|k| uv_curve.point_at(u0 + (u1 - u0) * k as f64 / (n - 1) as f64))
                .collect();
            let uv_params: Vec<f64> = (0..n)
                .map(|k| u0 + (u1 - u0) * k as f64 / (n - 1) as f64)
                .collect();
            let near = uv_pts[idx];
            if let Some((nu, nv)) = project_point_to_uv_near(surf, &r, near.x, near.y) {
                let (nu, nv) = unwrap_toward(surf, (nu, nv), (near.x, near.y));
                uv_pts[idx] = Point3::new(nu, nv, 0.0);
                if let Ok(fitted) =
                    popkern_core::geometry::fitting::fit_curve(&uv_pts, &uv_params, 3, n)
                {
                    *uv_curve = fitted;
                } else {
                    set_cp(uv_curve, Point3::new(nu, nv, 0.0));
                }
            }
        }
        return;
    }
    set_cp(&mut piece.curve_3d, r);
    for (uv_curve, surf) in [
        (&mut piece.curve_uv1, surf_a),
        (&mut piece.curve_uv2, surf_b),
    ] {
        let (t0, t1) = uv_curve.parameter_range();
        let cur = uv_curve.point_at(if end { t1 } else { t0 });
        if let Some((nu, nv)) = project_point_to_uv_near(surf, &r, cur.x, cur.y) {
            let (nu, nv) = unwrap_toward(surf, (nu, nv), (cur.x, cur.y));
            set_cp(uv_curve, Point3::new(nu, nv, 0.0));
        }
    }
}
