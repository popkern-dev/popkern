use crate::geometry::Point3;
use crate::geometry::curve::CurveEval;
use crate::geometry::intersection::surface_surface::SurfaceSurfaceIntersection;
use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};
use crate::topology::{BRep, FaceHandle, LoopHandle};
use std::f64::consts::TAU;

impl BRep {
    pub(super) fn intersection_overlaps_face_trims(
        &self,
        face_a: FaceHandle,
        surface_a: &SurfaceGeometry,
        face_b: FaceHandle,
        surface_b: &SurfaceGeometry,
        intersection: &SurfaceSurfaceIntersection,
        linear_tol: f64,
    ) -> bool {
        if let SurfaceSurfaceIntersection::Curve(curve) = intersection {
            if !surface_uses_iterative_projection(surface_a)
                && !surface_uses_iterative_projection(surface_b)
            {
                let points = intersection_sample_points(surface_a, intersection);
                return self.sampled_points_overlap_faces(
                    face_a, surface_a, face_b, surface_b, &points, linear_tol,
                );
            }

            let (t0, t1) = curve.curve_3d.parameter_range();
            let (uv1_t0, uv1_t1) = curve.curve_uv1.parameter_range();
            let (uv2_t0, uv2_t1) = curve.curve_uv2.parameter_range();
            if !t0.is_finite()
                || !t1.is_finite()
                || !uv1_t0.is_finite()
                || !uv1_t1.is_finite()
                || !uv2_t0.is_finite()
                || !uv2_t1.is_finite()
                || (t1 - t0).abs() < f64::EPSILON
            {
                return true;
            }

            const SAMPLE_COUNT: usize = 9;
            return (0..SAMPLE_COUNT).any(|i| {
                let alpha = i as f64 / (SAMPLE_COUNT - 1) as f64;
                let uv_a = curve.curve_uv1.point_at(uv1_t0 + (uv1_t1 - uv1_t0) * alpha);
                let uv_b = curve.curve_uv2.point_at(uv2_t0 + (uv2_t1 - uv2_t0) * alpha);
                let trim_a_contains = self.face_trim_classify_uv(
                    face_a,
                    surface_a,
                    UvPoint {
                        u: uv_a.x,
                        v: uv_a.y,
                    },
                );
                let trim_b_contains = self.face_trim_classify_uv(
                    face_b,
                    surface_b,
                    UvPoint {
                        u: uv_b.x,
                        v: uv_b.y,
                    },
                );
                trace_trim_overlap_sample(
                    face_a,
                    surface_a,
                    trim_a_contains,
                    UvPoint {
                        u: uv_a.x,
                        v: uv_a.y,
                    },
                    face_b,
                    surface_b,
                    trim_b_contains,
                    UvPoint {
                        u: uv_b.x,
                        v: uv_b.y,
                    },
                );
                trim_classes_overlap(trim_a_contains, trim_b_contains)
            });
        }

        let points = intersection_sample_points(surface_a, intersection);
        if points.is_empty() {
            return true;
        }

        self.sampled_points_overlap_faces(face_a, surface_a, face_b, surface_b, &points, linear_tol)
    }

    pub(super) fn sampled_points_overlap_faces(
        &self,
        face_a: FaceHandle,
        surface_a: &SurfaceGeometry,
        face_b: FaceHandle,
        surface_b: &SurfaceGeometry,
        points: &[Point3],
        linear_tol: f64,
    ) -> bool {
        let mut boundary_contacts = Vec::new();
        for point in points {
            let trim_a_contains =
                self.face_trim_classify_point(face_a, surface_a, point, linear_tol);
            let trim_b_contains =
                self.face_trim_classify_point(face_b, surface_b, point, linear_tol);
            trace_trim_overlap_point(
                face_a,
                surface_a,
                trim_a_contains,
                face_b,
                surface_b,
                trim_b_contains,
                point,
            );
            match (trim_a_contains, trim_b_contains) {
                (Some(TrimContainment::Outside), _) | (_, Some(TrimContainment::Outside)) => {}
                (Some(TrimContainment::Boundary), Some(TrimContainment::Boundary)) => {
                    if !self.faces_share_vertex_at_point(face_a, face_b, point, linear_tol) {
                        boundary_contacts.push(*point);
                    }
                }
                _ => return true,
            }
        }

        boundary_contacts.iter().enumerate().any(|(i, a)| {
            boundary_contacts
                .iter()
                .skip(i + 1)
                .any(|b| a.distance_to(b) > linear_tol * 10.0)
        })
    }

    #[allow(dead_code)]
    fn trim_classes_overlap_at_point(
        &self,
        face_a: FaceHandle,
        face_b: FaceHandle,
        a: Option<TrimContainment>,
        b: Option<TrimContainment>,
        point: &Point3,
        linear_tol: f64,
    ) -> bool {
        match (a, b) {
            (Some(TrimContainment::Outside), _) | (_, Some(TrimContainment::Outside)) => false,
            (Some(TrimContainment::Boundary), Some(TrimContainment::Boundary)) => {
                !self.faces_share_vertex_at_point(face_a, face_b, point, linear_tol)
            }
            _ => true,
        }
    }

    fn faces_share_vertex_at_point(
        &self,
        face_a: FaceHandle,
        face_b: FaceHandle,
        point: &Point3,
        linear_tol: f64,
    ) -> bool {
        let Some(vertices_a) = self.face_outer_loop_vertices(face_a) else {
            return false;
        };
        let Some(vertices_b) = self.face_outer_loop_vertices(face_b) else {
            return false;
        };
        vertices_a.iter().any(|vertex| {
            vertices_b.contains(vertex)
                && self.vertices[*vertex].position.distance_to(point) <= linear_tol * 10.0
        })
    }

    fn face_outer_loop_vertices(
        &self,
        face: FaceHandle,
    ) -> Option<Vec<crate::topology::VertexHandle>> {
        let face_data = self.faces.get(face)?;
        let loop_data = self.loops.get(face_data.outer_loop)?;
        let start = loop_data.half_edge;
        if !self.half_edges.contains_key(start) {
            return None;
        }

        let mut out = Vec::new();
        let mut cursor = start;
        for _ in 0..=self.half_edges.len() {
            let half_edge = self.half_edges.get(cursor)?;
            out.push(half_edge.origin);
            cursor = half_edge.next;
            if cursor == start {
                return Some(out);
            }
        }
        None
    }

    fn face_trim_classify_uv(
        &self,
        face_handle: FaceHandle,
        surface: &SurfaceGeometry,
        query: UvPoint,
    ) -> Option<TrimContainment> {
        let face = self.faces.get(face_handle)?;
        let param_tol = self
            .tolerance
            .parametric_base()
            .max(self.tolerance.linear())
            .max(1.0e-9)
            * 10.0;

        let outer = self.loop_uv_polygon(face.outer_loop, surface, query)?;
        let mut containment = point_in_uv_loop(query, &outer, param_tol);
        if containment == TrimContainment::Outside {
            return Some(TrimContainment::Outside);
        }

        for &inner_loop in &face.inner_loops {
            if let Some(inner) = self.loop_uv_polygon(inner_loop, surface, query) {
                match point_in_uv_loop(query, &inner, param_tol) {
                    TrimContainment::Interior => return Some(TrimContainment::Outside),
                    TrimContainment::Boundary => containment = TrimContainment::Boundary,
                    TrimContainment::Outside => {}
                }
            }
        }

        Some(containment)
    }

    fn face_trim_classify_point(
        &self,
        face_handle: FaceHandle,
        surface: &SurfaceGeometry,
        point: &Point3,
        linear_tol: f64,
    ) -> Option<TrimContainment> {
        let (u, v) = crate::geometry::uv_projection::project_point_to_uv(surface, point)?;
        let projected = surface.point_at(u, v);
        if projected.distance_to(point) > linear_tol * 10.0 {
            return Some(TrimContainment::Outside);
        }

        self.face_trim_classify_uv(face_handle, surface, UvPoint { u, v })
    }

    fn loop_uv_polygon(
        &self,
        loop_handle: LoopHandle,
        surface: &SurfaceGeometry,
        query: UvPoint,
    ) -> Option<Vec<UvPoint>> {
        let loop_ = self.loops.get(loop_handle)?;
        let start_he = loop_.half_edge;
        if !self.half_edges.contains_key(start_he) {
            return None;
        }

        const SAMPLES_PER_PCURVE: usize = 6;
        let mut polygon = Vec::new();
        let mut current = start_he;
        let walk_bound = self.half_edges.len().max(1);
        for _ in 0..=walk_bound {
            let he = self.half_edges.get(current)?;
            if let Some(pcurve) = &he.pcurve {
                for sample in 0..SAMPLES_PER_PCURVE {
                    let t = pcurve.t_start
                        + (pcurve.t_end - pcurve.t_start)
                            * (sample as f64 / SAMPLES_PER_PCURVE as f64);
                    let uv = pcurve.point_at(t);
                    push_distinct_uv(
                        &mut polygon,
                        unwrap_uv_for_query(
                            UvPoint { u: uv.x, v: uv.y },
                            query,
                            surface_u_is_periodic(surface),
                        ),
                    );
                }
            } else {
                let vertex = self.vertices.get(he.origin)?;
                let (u, v) =
                    crate::geometry::uv_projection::project_point_to_uv(surface, &vertex.position)?;
                push_distinct_uv(
                    &mut polygon,
                    unwrap_uv_for_query(UvPoint { u, v }, query, surface_u_is_periodic(surface)),
                );
            }

            current = he.next;
            if current == start_he {
                break;
            }
        }

        if polygon.len() >= 3 {
            Some(polygon)
        } else {
            None
        }
    }
}

#[derive(Debug, Clone, Copy)]
struct UvPoint {
    u: f64,
    v: f64,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum TrimContainment {
    Outside,
    Boundary,
    Interior,
}

fn trim_classes_overlap(a: Option<TrimContainment>, b: Option<TrimContainment>) -> bool {
    match (a, b) {
        (Some(TrimContainment::Outside), _) | (_, Some(TrimContainment::Outside)) => false,
        _ => true,
    }
}

fn intersection_sample_points(
    surface_a: &SurfaceGeometry,
    intersection: &SurfaceSurfaceIntersection,
) -> Vec<Point3> {
    match intersection {
        SurfaceSurfaceIntersection::Curve(curve) => {
            let (t0, t1) = curve.curve_3d.parameter_range();
            if !t0.is_finite() || !t1.is_finite() || (t1 - t0).abs() < f64::EPSILON {
                return Vec::new();
            }

            const SAMPLE_COUNT: usize = 9;
            (0..SAMPLE_COUNT)
                .map(|i| {
                    let alpha = i as f64 / (SAMPLE_COUNT - 1) as f64;
                    curve.curve_3d.point_at(t0 + (t1 - t0) * alpha)
                })
                .collect()
        }
        SurfaceSurfaceIntersection::CoincidentRegion {
            u1_range, v1_range, ..
        } => {
            let mut points = Vec::with_capacity(9);
            for u_alpha in [0.0, 0.5, 1.0] {
                for v_alpha in [0.0, 0.5, 1.0] {
                    let u = u1_range.0 + (u1_range.1 - u1_range.0) * u_alpha;
                    let v = v1_range.0 + (v1_range.1 - v1_range.0) * v_alpha;
                    points.push(surface_a.point_at(u, v));
                }
            }
            points
        }
    }
}

fn surface_u_is_periodic(surface: &SurfaceGeometry) -> bool {
    matches!(
        surface,
        SurfaceGeometry::Cylinder { .. }
            | SurfaceGeometry::Cone { .. }
            | SurfaceGeometry::Sphere { .. }
            | SurfaceGeometry::Torus { .. }
    )
}

fn surface_uses_iterative_projection(surface: &SurfaceGeometry) -> bool {
    matches!(
        surface,
        SurfaceGeometry::Nurbs(_) | SurfaceGeometry::Gregory(_)
    )
}

fn trace_trim_overlap_sample(
    face_a: FaceHandle,
    surface_a: &SurfaceGeometry,
    trim_a_contains: Option<TrimContainment>,
    uv_a: UvPoint,
    face_b: FaceHandle,
    surface_b: &SurfaceGeometry,
    trim_b_contains: Option<TrimContainment>,
    uv_b: UvPoint,
) {
    if std::env::var_os("POPKERN_SELF_INTERSECTION_TRACE").is_none() {
        return;
    }
    eprintln!(
        "[self-intersection trim] {face_a:?}/{} uv=({:.6},{:.6}) inside={trim_a_contains:?} \
         {face_b:?}/{} uv=({:.6},{:.6}) inside={trim_b_contains:?}",
        surface_kind(surface_a),
        uv_a.u,
        uv_a.v,
        surface_kind(surface_b),
        uv_b.u,
        uv_b.v,
    );
}

fn trace_trim_overlap_point(
    face_a: FaceHandle,
    surface_a: &SurfaceGeometry,
    trim_a_contains: Option<TrimContainment>,
    face_b: FaceHandle,
    surface_b: &SurfaceGeometry,
    trim_b_contains: Option<TrimContainment>,
    point: &Point3,
) {
    if std::env::var_os("POPKERN_SELF_INTERSECTION_TRACE").is_none() {
        return;
    }
    eprintln!(
        "[self-intersection trim] {face_a:?}/{} inside={trim_a_contains:?} \
         {face_b:?}/{} inside={trim_b_contains:?} point=({:.6},{:.6},{:.6})",
        surface_kind(surface_a),
        surface_kind(surface_b),
        point.x,
        point.y,
        point.z,
    );
}

fn surface_kind(surface: &SurfaceGeometry) -> &'static str {
    match surface {
        SurfaceGeometry::Plane { .. } => "Plane",
        SurfaceGeometry::Cylinder { .. } => "Cylinder",
        SurfaceGeometry::Cone { .. } => "Cone",
        SurfaceGeometry::Sphere { .. } => "Sphere",
        SurfaceGeometry::Torus { .. } => "Torus",
        SurfaceGeometry::Nurbs(_) => "Nurbs",
        SurfaceGeometry::Gregory(_) => "Gregory",
    }
}

fn unwrap_uv_for_query(uv: UvPoint, query: UvPoint, periodic_u: bool) -> UvPoint {
    if !periodic_u {
        return uv;
    }

    let mut u = uv.u;
    while u - query.u > TAU * 0.5 {
        u -= TAU;
    }
    while query.u - u > TAU * 0.5 {
        u += TAU;
    }
    UvPoint { u, v: uv.v }
}

fn push_distinct_uv(polygon: &mut Vec<UvPoint>, uv: UvPoint) {
    let duplicate = polygon
        .last()
        .is_some_and(|last| (last.u - uv.u).hypot(last.v - uv.v) < 1.0e-12);
    if !duplicate {
        polygon.push(uv);
    }
}

fn point_in_uv_loop(point: UvPoint, loop_polygon: &[UvPoint], tolerance: f64) -> TrimContainment {
    if loop_polygon.len() < 3 {
        return TrimContainment::Outside;
    }

    let mut previous = *loop_polygon.last().expect("length checked above");
    for &current in loop_polygon {
        if distance_to_uv_segment(point, previous, current) <= tolerance {
            return TrimContainment::Boundary;
        }
        previous = current;
    }

    let mut inside = false;
    previous = *loop_polygon.last().expect("length checked above");
    for &current in loop_polygon {
        let crosses_v = (current.v > point.v) != (previous.v > point.v);
        if crosses_v {
            let u_at_v = (previous.u - current.u) * (point.v - current.v)
                / (previous.v - current.v)
                + current.u;
            if point.u < u_at_v {
                inside = !inside;
            }
        }
        previous = current;
    }
    if inside {
        TrimContainment::Interior
    } else {
        TrimContainment::Outside
    }
}

fn distance_to_uv_segment(point: UvPoint, a: UvPoint, b: UvPoint) -> f64 {
    let ab_u = b.u - a.u;
    let ab_v = b.v - a.v;
    let ap_u = point.u - a.u;
    let ap_v = point.v - a.v;
    let len_sq = ab_u * ab_u + ab_v * ab_v;
    if len_sq <= f64::EPSILON {
        return ap_u.hypot(ap_v);
    }

    let t = ((ap_u * ab_u + ap_v * ab_v) / len_sq).clamp(0.0, 1.0);
    let closest = UvPoint {
        u: a.u + ab_u * t,
        v: a.v + ab_v * t,
    };
    (point.u - closest.u).hypot(point.v - closest.v)
}
