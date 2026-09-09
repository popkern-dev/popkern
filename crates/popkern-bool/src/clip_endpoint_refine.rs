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
    tolerance: &Tolerance,
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
        let p = endpoint(piece, end);
        // The record lives on BOTH faces; a boundary junction on
        // either face's rim is solved against the OTHER face's
        // surface (the record's second carrier).
        let refined = refine_against_face(brep_a, face_a, surf_b, windows_b, &p, band, tolerance)
            .or_else(|| {
                refine_against_face(brep_b, face_b, surf_a, windows_a, &p, band, tolerance)
            });
        if let Some(r) = refined
            && r.distance_to(&p) > tolerance.linear()
            && r.distance_to(&p) < band
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
            set_endpoint_3d_uv(piece, end, r, surf_a, surf_b);
        }
    }
}

fn endpoint(piece: &CurveOnSurface, end: bool) -> Point3 {
    let (t0, t1) = piece.curve_3d.parameter_range();
    piece.curve_3d.point_at(if end { t1 } else { t0 })
}

/// Exact junction near `p`: the nearest boundary edge of `face`
/// (within the slop band) intersected with `other_surface`.
fn refine_against_face(
    brep: &BRep,
    face: FaceHandle,
    other_surface: &SurfaceGeometry,
    other_windows: ((f64, f64), (f64, f64)),
    p: &Point3,
    band: f64,
    tolerance: &Tolerance,
) -> Option<Point3> {
    let face_data = brep.faces.get(face)?;
    let lp = brep.loops.get(face_data.outer_loop)?;
    let start = lp.half_edge;
    let mut he = start;
    let mut best: Option<Point3> = None;
    for _ in 0..brep.half_edges.len().max(1) {
        let h = brep.half_edges.get(he)?;
        if let Some(edge) = brep.edges.get(h.edge)
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
            let near = (&c.geometry as &dyn popkern_core::geometry::curve::CurveEval)
                .closest_point(p)
                .map(|t| {
                    let tc = t.clamp(t0.min(t1), t0.max(t1));
                    c.geometry.point_at(tc).distance_to(p) < band_eff
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
                        if let SurfaceCurveIntersection::Point { point, .. } = isect {
                            let d = point.distance_to(p);
                            if d < band_eff
                                && best.is_none_or(|b| point.distance_to(p) < b.distance_to(p))
                            {
                                best = Some(*point);
                            }
                            let _ = d;
                        }
                    }
                }
            }
        }
        he = h.next;
        if he == start {
            break;
        }
    }
    best
}

fn clone_curve_geometry(g: &CurveGeometry) -> CurveGeometry {
    g.clone()
}

/// Move one endpoint of the piece: the 3D end control point, and each
/// UV chart's end control point re-projected near its current value
/// (clamped NURBS — the endpoint IS a control point).
fn set_endpoint_3d_uv(
    piece: &mut CurveOnSurface,
    end: bool,
    r: Point3,
    surf_a: &SurfaceGeometry,
    surf_b: &SurfaceGeometry,
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
    set_cp(&mut piece.curve_3d, r);
    for (uv_curve, surf) in [
        (&mut piece.curve_uv1, surf_a),
        (&mut piece.curve_uv2, surf_b),
    ] {
        let (t0, t1) = uv_curve.parameter_range();
        let cur = uv_curve.point_at(if end { t1 } else { t0 });
        if let Some((nu, nv)) = project_point_to_uv_near(surf, &r, cur.x, cur.y) {
            set_cp(uv_curve, Point3::new(nu, nv, 0.0));
        }
    }
}
