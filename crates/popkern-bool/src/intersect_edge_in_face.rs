//! Operand edges lying in the other operand's face (edge-in-face
//! interference, PK-971 / plan §4.3a).
//!
//! When an edge of one operand lies, within tolerance, inside a face of
//! the other, the face must be imprinted along it: the edge bounds the
//! two regions of that face. The face-face solver does not reliably
//! deliver that curve, because it runs exactly along the boundary of
//! the edge's own face, where trim clipping keeps or drops it on
//! rounding. Measured: a hinge leaf's strip edge lying in the other
//! leaf's coplanar-opposite bottom face had no record, so the bottom
//! face's fragment spanned the 0.04-wide sliver and the cell complex
//! declined an ambiguous radial order on the knuckle rim (PK-1017).
//!
//! This pass recovers the curve from the interference itself, as
//! OCCT's edge/face interference records an edge "in" a face: for
//! every edge of X that lies on a planar face f of Y and inside f's
//! region along its whole length, each face g of X at the edge (other
//! than one coincident with f, whose overlap the coincidence engine
//! owns) gets a record (g, f) along the edge — unless the solver
//! already delivered one covering it. The record's uv curves are exact
//! where the geometry allows (affine projection onto planes) and
//! verified against the 3D curve otherwise; an unverifiable image adds
//! nothing, leaving today's behaviour.
//!
//! Scope: planar receiving faces, the measured cases (the region test
//! is exact on non-periodic charts only), and edges lying on the plane
//! along their whole length; each span of such an edge inside the
//! face's region is imprinted. Curved receiving faces, and edges that
//! leave the plane part-way, are follow-ups.

use popkern_core::geometry::curve::{CurveEval, CurveGeometry, closest_point_in_window};
use popkern_core::geometry::curve_on_surface::CurveOnSurface;
use popkern_core::geometry::intersection::coincident_analytic::are_surfaces_coincident_analytic;
use popkern_core::geometry::intersection::surface_subdivision::Aabb;
use popkern_core::geometry::nurbs_curve::NurbsCurve;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle};

use crate::intersect::{EndpointBoundaryType, IntersectionRecord};
use crate::intersect_edge_in_face_spans::{PlaneChart, inside_spans};

/// Samples along an edge for the on-face and coverage tests.
const EDGE_SAMPLES: usize = 16;
/// Samples for a uv image that is not an exact projection.
const UV_SAMPLES: usize = 32;

/// Append the records of every operand edge lying in the other
/// operand's planar face that the solver did not deliver.
pub(crate) fn recover_edge_in_face_records(
    a: &BRep,
    b: &BRep,
    aabbs_a: &[(FaceHandle, Aabb)],
    aabbs_b: &[(FaceHandle, Aabb)],
    records: &mut Vec<IntersectionRecord>,
    tolerance: &Tolerance,
) {
    let mut found = Vec::new();
    for (x_is_a, x, y, aabbs_y) in [(true, a, b, aabbs_b), (false, b, a, aabbs_a)] {
        // The receiving faces' hull, padded as the per-face test pads
        // the sample box: an edge whose conservative box stays clear of
        // it meets none of them (PK-1032: rejected before sampling).
        let corners: Vec<Point3> = aabbs_y.iter().flat_map(|(_, f)| [f.min, f.max]).collect();
        let region = sample_box(&corners, 2.0 * tolerance.linear());
        for (eh, edge) in &x.edges {
            let Some(curve) = edge.curve.as_ref() else {
                continue;
            };
            let mut hull = vec![
                curve.geometry.point_at(curve.t_start),
                curve.geometry.point_at(curve.t_end),
            ];
            match &curve.geometry {
                // The whole control polygon holds the span, unclipped.
                CurveGeometry::Nurbs(n) if n.weights.iter().all(|w| *w > 0.0) => {
                    hull.extend_from_slice(&n.control_points)
                }
                _ => crate::intersect::conservative_curve_extent(curve, &mut hull),
            }
            if !sample_box(&hull, 0.0).overlaps(&region, 0.0) {
                continue;
            }
            if !curve_meets_vertices(x, edge.half_edge, curve, tolerance) {
                continue;
            }
            let samples: Vec<Point3> = (0..=EDGE_SAMPLES)
                .map(|k| {
                    let s = k as f64 / EDGE_SAMPLES as f64;
                    curve
                        .geometry
                        .point_at(curve.t_start + s * (curve.t_end - curve.t_start))
                })
                .collect();
            let ebox = sample_box(&samples, tolerance.linear());
            for (fh, fbox) in aabbs_y {
                if !ebox.overlaps(fbox, tolerance.linear()) {
                    continue;
                }
                let Some(chart) = plane_chart(y, *fh) else {
                    continue;
                };
                // The whole edge on the plane; its spans in the region.
                if !samples.iter().all(|p| {
                    (*p - chart.chart.origin).dot(&chart.normal).abs() <= tolerance.linear()
                }) {
                    continue;
                }
                let window = (
                    curve.t_start.min(curve.t_end),
                    curve.t_start.max(curve.t_end),
                );
                let spans = inside_spans(&curve.geometry, window, y, *fh, &chart.chart, tolerance);
                let faces = edge_faces(x, eh);
                let pair = |g: FaceHandle| if x_is_a { (g, *fh) } else { (*fh, g) };
                for (span, &(g, he)) in spans
                    .iter()
                    .flat_map(|s| faces.iter().map(move |gh| (*s, gh)))
                {
                    let (face_a, face_b) = pair(g);
                    let span_samples: Vec<Point3> = (0..=EDGE_SAMPLES)
                        .map(|k| {
                            let s = k as f64 / EDGE_SAMPLES as f64;
                            curve.geometry.point_at(span.0 + s * (span.1 - span.0))
                        })
                        .collect();
                    // A record of either face along the edge already
                    // imprints `fh` there: the clip gives a cut riding
                    // the edge to one of its faces (half-open), so a
                    // second copy would have both claim it.
                    let delivered = faces.iter().any(|&(g2, _)| {
                        let (fa, fb) = pair(g2);
                        covered(records, fa, fb, &span_samples, tolerance)
                            || covered(&found, fa, fb, &span_samples, tolerance)
                    });
                    if coincident_faces(x, g, y, *fh, tolerance) || delivered {
                        continue;
                    }
                    let Some(record) =
                        edge_record(x, g, he, y, *fh, curve, span, x_is_a, tolerance)
                    else {
                        if popkern_core::trace::trace_enabled("boolean") {
                            eprintln!(
                                "  [edge-in-face] {eh:?} in {fh:?}: uv image unverified, skipped"
                            );
                        }
                        continue;
                    };
                    if popkern_core::trace::trace_enabled("boolean") {
                        eprintln!(
                            "  [edge-in-face] {} edge {eh:?} lies in {fh:?}: record ({face_a:?},{face_b:?})",
                            if x_is_a { "A" } else { "B" }
                        );
                    }
                    found.push(record);
                }
            }
        }
    }
    records.extend(found);
}

/// Whether the edge's curve ends at its vertices (either way round).
/// A body whose vertices were moved without its curves has edges that
/// lie nowhere near its faces' boundaries; reading those curves would
/// imprint where the body is not.
fn curve_meets_vertices(
    brep: &BRep,
    he: HalfEdgeHandle,
    curve: &popkern_core::geometry::CurveOnEdge,
    tolerance: &Tolerance,
) -> bool {
    let vertex = |h: HalfEdgeHandle| {
        brep.half_edges
            .get(h)
            .and_then(|d| brep.vertices.get(d.origin))
            .map(|v| v.position)
    };
    let Some(twin) = brep.half_edges.get(he).and_then(|h| h.mate()) else {
        return false;
    };
    let (Some(p), Some(q)) = (vertex(he), vertex(twin)) else {
        return false;
    };
    let (s, e) = (
        curve.geometry.point_at(curve.t_start),
        curve.geometry.point_at(curve.t_end),
    );
    let reach = tolerance.linear() * 10.0;
    (s.distance_to(&p) <= reach && e.distance_to(&q) <= reach)
        || (s.distance_to(&q) <= reach && e.distance_to(&p) <= reach)
}

/// Box of `pts`, padded by `pad`.
fn sample_box(pts: &[Point3], pad: f64) -> Aabb {
    let mut lo = Point3::new(f64::MAX, f64::MAX, f64::MAX);
    let mut hi = Point3::new(f64::MIN, f64::MIN, f64::MIN);
    for p in pts {
        lo = Point3::new(lo.x.min(p.x), lo.y.min(p.y), lo.z.min(p.z));
        hi = Point3::new(hi.x.max(p.x), hi.y.max(p.y), hi.z.max(p.z));
    }
    let d = Vec3::new(pad, pad, pad);
    Aabb {
        min: lo - d,
        max: hi + d,
    }
}

/// Face `fh`'s plane chart and normal, when its surface is a plane.
struct PlanarFace {
    chart: PlaneChart,
    normal: Vec3,
}

/// Face `fh`'s plane, if it is planar.
fn plane_chart(brep: &BRep, fh: FaceHandle) -> Option<PlanarFace> {
    let surf = brep.faces.get(fh)?.surface.as_ref()?;
    let SurfaceGeometry::Plane {
        origin,
        normal,
        u_axis,
        v_axis,
    } = &surf.geometry
    else {
        return None;
    };
    Some(PlanarFace {
        chart: PlaneChart {
            origin: *origin,
            u_axis: *u_axis,
            v_axis: *v_axis,
        },
        normal: *normal,
    })
}

/// The faces at an edge, each with its half-edge there (a seam's face
/// once).
fn edge_faces(brep: &BRep, eh: EdgeHandle) -> Vec<(FaceHandle, HalfEdgeHandle)> {
    let mut out: Vec<(FaceHandle, HalfEdgeHandle)> = Vec::new();
    let Some(he0) = brep.edges.get(eh).map(|e| e.half_edge) else {
        return out;
    };
    let twin = brep.half_edges.get(he0).and_then(|h| h.mate());
    for he in std::iter::once(he0).chain(twin) {
        let Some(face) = brep
            .half_edges
            .get(he)
            .and_then(|h| brep.loops.get(h.loop_))
            .map(|l| l.face)
        else {
            continue;
        };
        if !out.iter().any(|(f, _)| *f == face) {
            out.push((face, he));
        }
    }
    out
}

/// Whether faces `g` of `x` and `f` of `y` lie on one surface.
fn coincident_faces(x: &BRep, g: FaceHandle, y: &BRep, f: FaceHandle, tol: &Tolerance) -> bool {
    match (
        x.faces.get(g).and_then(|f| f.surface.as_ref()),
        y.faces.get(f).and_then(|f| f.surface.as_ref()),
    ) {
        (Some(sg), Some(sf)) => are_surfaces_coincident_analytic(&sg.geometry, &sf.geometry, tol),
        _ => false,
    }
}

/// Whether a record of the pair already runs along every sample. Read
/// on the record's curve, not a polyline of it: a chord of a circle's
/// 64-gon sags far past tolerance (1.6e-3 at radius 1.3).
fn covered(
    records: &[IntersectionRecord],
    face_a: FaceHandle,
    face_b: FaceHandle,
    samples: &[Point3],
    tolerance: &Tolerance,
) -> bool {
    let reach = tolerance.linear() * 10.0;
    let pair: Vec<&IntersectionRecord> = records
        .iter()
        .filter(|r| r.face_a == face_a && r.face_b == face_b)
        .collect();
    !pair.is_empty()
        && samples.iter().all(|p| {
            pair.iter().any(|r| {
                let c = &r.curve.curve_3d;
                let (t0, t1) = c.parameter_range();
                c.point_at(closest_point_in_window(c, p, t0, t1))
                    .distance_to(p)
                    <= reach
            })
        })
}

/// The record (g, f) along the edge, oriented as the edge's curve.
#[allow(clippy::too_many_arguments)]
fn edge_record(
    x: &BRep,
    g: FaceHandle,
    he: HalfEdgeHandle,
    y: &BRep,
    f: FaceHandle,
    curve: &popkern_core::geometry::CurveOnEdge,
    (t0, t1): (f64, f64),
    x_is_a: bool,
    tolerance: &Tolerance,
) -> Option<IntersectionRecord> {
    let window = (
        curve.t_start.min(curve.t_end),
        curve.t_start.max(curve.t_end),
    );
    let curve_3d = curve.geometry.to_nurbs_in_range((t0, t1));
    let surf_g = &x.faces.get(g)?.surface.as_ref()?.geometry;
    let surf_f = &y.faces.get(f)?.surface.as_ref()?.geometry;
    let uv_f = plane_image(surf_f, &curve_3d)?;
    let uv_g = match plane_image(surf_g, &curve_3d) {
        Some(uv) => uv,
        None => pcurve_image(x, he, surf_g, &curve.geometry, window, (t0, t1), tolerance)?,
    };
    let mid = curve_3d.point_at(0.5 * (t0 + t1));
    let (ug, vg) = surf_g.closest_point(&mid)?;
    let (uf, vf) = surf_f.closest_point(&mid)?;
    let (ng, nf) = (surf_g.normal_at(ug, vg), surf_f.normal_at(uf, vf));
    let is_tangent = popkern_core::adaptive_precision::adaptive_cross_ratio_sq_lt(
        (ng.x, ng.y, ng.z),
        (nf.x, nf.y, nf.z),
        tolerance.angular(),
    );
    let (face_a, face_b, curve_uv1, curve_uv2) = if x_is_a {
        (g, f, uv_g, uv_f)
    } else {
        (f, g, uv_f, uv_g)
    };
    Some(IntersectionRecord {
        face_a,
        face_b,
        curve: CurveOnSurface {
            curve_3d,
            curve_uv1,
            curve_uv2,
            is_closed: false,
        },
        is_tangent,
        is_coincident: false,
        start_boundary: EndpointBoundaryType::OnVertex,
        end_boundary: EndpointBoundaryType::OnVertex,
        residue_owner_a: None,
    })
}

/// A 3D NURBS curve's uv image on a plane: the plane chart is affine,
/// so projecting the control points (same knots, same weights) is the
/// exact image.
fn plane_image(surf: &SurfaceGeometry, c: &NurbsCurve) -> Option<NurbsCurve> {
    let SurfaceGeometry::Plane {
        origin,
        u_axis,
        v_axis,
        ..
    } = surf
    else {
        return None;
    };
    Some(NurbsCurve {
        degree: c.degree,
        knots: c.knots.clone(),
        control_points: c
            .control_points
            .iter()
            .map(|p| {
                let d = *p - *origin;
                Point3::new(d.dot(u_axis), d.dot(v_axis), 0.0)
            })
            .collect(),
        weights: c.weights.clone(),
    })
}

/// The edge's own pcurve on `g` (half-edge `he`), resampled on the 3D
/// curve's parameter and verified to image the curve within
/// tolerance; `None` when it does not.
fn pcurve_image(
    x: &BRep,
    he: HalfEdgeHandle,
    surf: &SurfaceGeometry,
    c: &popkern_core::geometry::CurveGeometry,
    (w0, w1): (f64, f64),
    (t0, t1): (f64, f64),
    tolerance: &Tolerance,
) -> Option<NurbsCurve> {
    let pc = x.half_edges.get(he)?.pcurve.as_ref()?;
    let start = c.point_at(w0);
    // The pcurve spans the whole edge window, either way along it.
    let at = |s: f64| pc.point_at(pc.t_start + s * (pc.t_end - pc.t_start));
    let image = |q: Point3| surf.point_at(q.x, q.y);
    let forward = image(at(0.0)).distance_to(&start) <= image(at(1.0)).distance_to(&start);
    let uv_at = |t: f64| {
        let s = (t - w0) / (w1 - w0);
        if forward { at(s) } else { at(1.0 - s) }
    };
    let n = if matches!(c, popkern_core::geometry::CurveGeometry::Line { .. }) {
        1
    } else {
        UV_SAMPLES
    };
    let params: Vec<f64> = (0..=n)
        .map(|k| t0 + (t1 - t0) * k as f64 / n as f64)
        .collect();
    let pts: Vec<Point3> = params
        .iter()
        .map(|&t| {
            let q = uv_at(t);
            Point3::new(q.x, q.y, 0.0)
        })
        .collect();
    // The polyline's image must follow the 3D curve, at its samples and
    // between them.
    for k in 0..n {
        for s in [0.0, 0.5] {
            let t = params[k] + s * (params[k + 1] - params[k]);
            let q = pts[k] + (pts[k + 1] - pts[k]) * s;
            if image(q).distance_to(&c.point_at(t)) > tolerance.linear() {
                return None;
            }
        }
    }
    if image(pts[n]).distance_to(&c.point_at(t1)) > tolerance.linear() {
        return None;
    }
    let mut knots = vec![params[0]];
    knots.extend(params.iter().copied());
    knots.push(params[n]);
    Some(NurbsCurve {
        degree: 1,
        knots,
        control_points: pts,
        weights: vec![1.0; n + 1],
    })
}

#[cfg(test)]
#[path = "intersect_edge_in_face_tests.rs"]
mod tests;
