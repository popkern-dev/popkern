use crate::geometry::Point3;
use crate::geometry::curve::CurveEval;
use crate::geometry::intersection::surface_surface::SurfaceSurfaceIntersection;
use crate::geometry::surface::{SurfaceEval, SurfaceGeometry};
use crate::topology::{BRep, FaceHandle, LoopHandle};
use self_intersection_trim_trace::{
    TrimOverlapTraceSample, trace_trim_overlap_point, trace_trim_overlap_sample,
};
use self_intersection_uv_containment::{
    TrimContainment, UvPoint, point_in_uv_loop_periodic, push_distinct_uv, unwrap_uv_for_query,
};

#[path = "self_intersection_trim_trace.rs"]
mod self_intersection_trim_trace;
#[path = "self_intersection_uv_containment.rs"]
mod self_intersection_uv_containment;

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

            // Same two-tier semantics as `sampled_points_overlap_faces`
            // below: any interior overlap flags immediately; contact
            // where BOTH trims classify Boundary is served-boundary
            // contact, legitimate at a shared topological vertex
            // (measured: variable-radius N=3 corner — the two strips'
            // internally tangent sections touch at the shared setback
            // corner; 1 terminal sample Boundary∧Boundary, all others
            // Outside) and a defect only with spatial extent.
            const SAMPLE_COUNT: usize = 9;
            let mut boundary_contacts: Vec<Point3> = Vec::new();
            for i in 0..SAMPLE_COUNT {
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
                    TrimOverlapTraceSample {
                        face: face_a,
                        surface: surface_a,
                        trim_contains: trim_a_contains,
                        uv: UvPoint {
                            u: uv_a.x,
                            v: uv_a.y,
                        },
                    },
                    TrimOverlapTraceSample {
                        face: face_b,
                        surface: surface_b,
                        trim_contains: trim_b_contains,
                        uv: UvPoint {
                            u: uv_b.x,
                            v: uv_b.y,
                        },
                    },
                );
                match (trim_a_contains, trim_b_contains) {
                    (Some(TrimContainment::Outside), _) | (_, Some(TrimContainment::Outside)) => {}
                    (Some(TrimContainment::Boundary), Some(TrimContainment::Boundary)) => {
                        let point = surface_a.point_at(uv_a.x, uv_a.y);
                        if !self.faces_share_vertex_at_point(face_a, face_b, &point, linear_tol) {
                            boundary_contacts.push(point);
                        }
                    }
                    _ => return true,
                }
            }
            return boundary_contacts.iter().enumerate().any(|(i, a)| {
                boundary_contacts
                    .iter()
                    .skip(i + 1)
                    .any(|b| a.distance_to(b) > linear_tol * 10.0)
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

        let periodic_u = surface_u_is_periodic(surface);
        let outer = self.loop_uv_polygon(face.outer_loop, surface, query)?;
        let mut containment = point_in_uv_loop_periodic(query, &outer, param_tol, periodic_u);
        if containment == TrimContainment::Outside {
            return Some(TrimContainment::Outside);
        }

        for &inner_loop in &face.inner_loops {
            if let Some(inner) = self.loop_uv_polygon(inner_loop, surface, query) {
                match point_in_uv_loop_periodic(query, &inner, param_tol, periodic_u) {
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
            // Chain-continuous unwrapping: each sample resolves its
            // period branch against the PREVIOUS polygon point (the
            // first against the query, only to anchor the overall
            // branch). Unwrapping every sample independently against
            // the query folds any loop spanning more than half a
            // period into a self-overlapping polygon, fabricating
            // containment for query points the trim never covers
            // (measured: half-period cylinder wall vs rear facet on
            // the laptop-stand block — phantom self-intersection).
            let periodic_u = surface_u_is_periodic(surface);
            if let Some(pcurve) = &he.pcurve {
                for sample in 0..SAMPLES_PER_PCURVE {
                    let t = pcurve.t_start
                        + (pcurve.t_end - pcurve.t_start)
                            * (sample as f64 / SAMPLES_PER_PCURVE as f64);
                    let uv = pcurve.point_at(t);
                    let reference = polygon.last().copied().unwrap_or(query);
                    push_distinct_uv(
                        &mut polygon,
                        unwrap_uv_for_query(UvPoint { u: uv.x, v: uv.y }, reference, periodic_u),
                    );
                }
            } else {
                let vertex = self.vertices.get(he.origin)?;
                let (u, v) =
                    crate::geometry::uv_projection::project_point_to_uv(surface, &vertex.position)?;
                let reference = polygon.last().copied().unwrap_or(query);
                push_distinct_uv(
                    &mut polygon,
                    unwrap_uv_for_query(UvPoint { u, v }, reference, periodic_u),
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
