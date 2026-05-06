//! Self-intersection validation check (#11).
//!
//! Detects when non-adjacent faces in a solid intersect each other.

use crate::error::ValidationError;
use crate::geometry::Point3;
use crate::geometry::intersection::degeneracy::DefaultDegeneracyHandler;
use crate::geometry::intersection::intersect_surfaces_with_handler;
use crate::geometry::intersection::surface_subdivision::Aabb;
use crate::geometry::surface::SurfaceEval;
use crate::tolerance::Tolerance;
use crate::topology::BRep;
use crate::validation::bvh::Bvh;
use std::collections::HashSet;

impl BRep {
    /// Check #11: no face intersects another non-adjacent face.
    ///
    /// The check proceeds in three stages:
    ///
    /// 1. Build a set of adjacent face-index pairs by walking shared half-edges.
    /// 2. Sample each face surface on a 5×5 grid to compute an AABB; build a
    ///    BVH from those AABBs for efficient spatial queries.
    /// 3. For each non-adjacent pair whose AABBs overlap (from BVH query), call
    ///    `intersect_surfaces_with_handler` with `DefaultDegeneracyHandler` and
    ///    a coarse tolerance (2× linear) to test for actual surface intersection.
    ///    Any intersections found are reported as a `SelfIntersection` error.
    pub fn validate_self_intersection(&self) -> Result<(), ValidationError> {
        let face_keys: Vec<_> = self.faces.keys().collect();
        let n = face_keys.len();
        if n < 2 {
            return Ok(());
        }

        // ------------------------------------------------------------------ //
        // 1. Build adjacency set (by face index pairs, not handles)           //
        // ------------------------------------------------------------------ //
        let mut adjacent: HashSet<(usize, usize)> = HashSet::new();

        for (i, &fk) in face_keys.iter().enumerate() {
            let face = &self.faces[fk];
            let start_he_handle = self.loops[face.outer_loop].half_edge;
            // Empty loop — skip
            if !self.half_edges.contains_key(start_he_handle) {
                continue;
            }
            let mut he_handle = start_he_handle;
            // Walk the outer loop; guard against cycles larger than the mesh.
            let bound = self.half_edges.len().max(1);
            for _ in 0..=bound {
                let he = match self.half_edges.get(he_handle) {
                    Some(h) => h,
                    None => break,
                };
                let twin_handle = he.twin;
                if let Some(twin) = self.half_edges.get(twin_handle) {
                    let twin_loop_handle = twin.loop_;
                    if let Some(twin_loop) = self.loops.get(twin_loop_handle) {
                        let twin_face = twin_loop.face;
                        if let Some(j) = face_keys.iter().position(|&k| k == twin_face)
                            && i != j
                        {
                            let pair = if i < j { (i, j) } else { (j, i) };
                            adjacent.insert(pair);
                        }
                    }
                }
                he_handle = he.next;
                if he_handle == start_he_handle {
                    break;
                }
            }
        }

        // ------------------------------------------------------------------ //
        // 2. Compute face AABBs via 5×5 sampling; build BVH                  //
        // ------------------------------------------------------------------ //
        let sample_params = [0.0, 0.25, 0.5, 0.75, 1.0];

        let aabbs: Vec<Option<Aabb>> = face_keys
            .iter()
            .map(|&fk| {
                let face = &self.faces[fk];
                face.surface.as_ref().and_then(|s| {
                    let positions = self.face_outer_loop_positions(fk);
                    let ((u0, u1), (v0, v1)) =
                        crate::geometry::uv_projection::tight_bounded_parameter_range(
                            &positions,
                            &s.geometry,
                        );
                    let mut pts = Vec::with_capacity(25);
                    for &tp in &sample_params {
                        let u = u0 + tp * (u1 - u0);
                        for &tq in &sample_params {
                            let v = v0 + tq * (v1 - v0);
                            pts.push(s.geometry.point_at(u, v));
                        }
                    }
                    Aabb::from_points(&pts)
                })
            })
            .collect();

        // Build BVH items: only include faces that have valid AABBs.
        let bvh_items: Vec<(usize, Aabb)> = aabbs
            .iter()
            .enumerate()
            .filter_map(|(i, aabb)| aabb.as_ref().map(|a| (i, a.clone())))
            .collect();

        if bvh_items.len() < 2 {
            return Ok(());
        }

        let bvh = Bvh::build(&bvh_items);

        // Use 2× linear tolerance for intersection detection: tight enough to catch
        // genuine self-intersections, loose enough for fast intersection testing.
        // This is conservative: may report near-misses as intersections (false positives)
        // but won't miss actual intersections within tolerance (no false negatives).
        let linear_tol = self.tolerance.linear();
        let coarse_tol = Tolerance::default()
            .with_linear(linear_tol * 2.0)
            .with_angular(self.tolerance.angular())
            .with_parametric(self.tolerance.parametric_base());

        // ------------------------------------------------------------------ //
        // 3. Test non-adjacent, AABB-overlapping pairs with intersect_surfaces //
        // ------------------------------------------------------------------ //
        // Validation budget: cap the number of pairs tested to prevent
        // validation from being more expensive than the operation that
        // produced the BRep. 200 pairs covers typical solids; beyond that,
        // diminishing returns on defect detection.
        const MAX_PAIRS_CHECKED: usize = 200;
        let overlapping = bvh.find_overlapping_pairs(linear_tol);
        let mut pairs_checked = 0usize;

        for (i, j) in overlapping {
            if pairs_checked >= MAX_PAIRS_CHECKED {
                break;
            }
            let pair = if i < j { (i, j) } else { (j, i) };
            if adjacent.contains(&pair) {
                continue;
            }

            let surf_i = match self.faces[face_keys[i]].surface.as_ref() {
                Some(s) => s,
                None => continue,
            };
            let surf_j = match self.faces[face_keys[j]].surface.as_ref() {
                Some(s) => s,
                None => continue,
            };

            // Use exact (zero-padding) ranges for intersection testing to
            // avoid false positives at shared edges where padding from two
            // perpendicular faces would overlap.
            let positions_i = self.face_outer_loop_positions(face_keys[i]);
            let ((u0i, u1i), (v0i, v1i)) =
                crate::geometry::uv_projection::exact_bounded_parameter_range(
                    &positions_i,
                    &surf_i.geometry,
                );
            let positions_j = self.face_outer_loop_positions(face_keys[j]);
            let ((u0j, u1j), (v0j, v1j)) =
                crate::geometry::uv_projection::exact_bounded_parameter_range(
                    &positions_j,
                    &surf_j.geometry,
                );

            // Skip faces with degenerate (zero-area) parameter ranges.
            if (u1i - u0i).abs() < f64::EPSILON || (v1i - v0i).abs() < f64::EPSILON {
                continue;
            }
            if (u1j - u0j).abs() < f64::EPSILON || (v1j - v0j).abs() < f64::EPSILON {
                continue;
            }

            pairs_checked += 1;
            trace_self_intersection_candidate(
                self,
                face_keys[i],
                face_keys[j],
                pairs_checked,
                (u0i, u1i),
                (v0i, v1i),
                (u0j, u1j),
                (v0j, v1j),
            );
            if let Some(overlaps) = self.bounded_cylinder_cylinder_trim_overlap(
                face_keys[i],
                &surf_i.geometry,
                (u0i, u1i),
                (v0i, v1i),
                face_keys[j],
                &surf_j.geometry,
                (u0j, u1j),
                (v0j, v1j),
                linear_tol,
            ) {
                if overlaps {
                    return Err(ValidationError::SelfIntersection {
                        face_a: face_keys[i],
                        face_b: face_keys[j],
                    });
                }
                continue;
            }
            let result = intersect_surfaces_with_handler(
                &surf_i.geometry,
                (u0i, u1i),
                (v0i, v1i),
                &surf_j.geometry,
                (u0j, u1j),
                (v0j, v1j),
                &coarse_tol,
                &DefaultDegeneracyHandler,
            );

            match result {
                Ok(ref r)
                    if r.intersections.iter().any(|intersection| {
                        self.intersection_overlaps_face_trims(
                            face_keys[i],
                            &surf_i.geometry,
                            face_keys[j],
                            &surf_j.geometry,
                            intersection,
                            linear_tol,
                        )
                    }) =>
                {
                    if std::env::var_os("POPKERN_SELF_INTERSECTION_TRACE").is_some() {
                        self.trace_self_intersection_pair(face_keys[i], face_keys[j]);
                    }
                    return Err(ValidationError::SelfIntersection {
                        face_a: face_keys[i],
                        face_b: face_keys[j],
                    });
                }
                _ => {}
            }
        }

        Ok(())
    }

    #[allow(clippy::too_many_arguments)]
    fn bounded_cylinder_cylinder_trim_overlap(
        &self,
        face_a: crate::topology::FaceHandle,
        surface_a: &crate::geometry::surface::SurfaceGeometry,
        u_range_a: (f64, f64),
        v_range_a: (f64, f64),
        face_b: crate::topology::FaceHandle,
        surface_b: &crate::geometry::surface::SurfaceGeometry,
        u_range_b: (f64, f64),
        v_range_b: (f64, f64),
        linear_tol: f64,
    ) -> Option<bool> {
        let samples = cylinder_cylinder_near_samples(
            surface_a, u_range_a, v_range_a, surface_b, u_range_b, v_range_b, linear_tol,
        )?;
        if samples.is_empty() {
            return Some(false);
        }
        Some(self.sampled_points_overlap_faces(
            face_a, surface_a, face_b, surface_b, &samples, linear_tol,
        ))
    }

    fn trace_self_intersection_pair(
        &self,
        face_a: crate::topology::FaceHandle,
        face_b: crate::topology::FaceHandle,
    ) {
        eprintln!(
            "[self-intersection pair] {face_a:?}/{} {face_b:?}/{}",
            trace_surface_kind(self, face_a),
            trace_surface_kind(self, face_b),
        );
        for face in [face_a, face_b] {
            let Some(face_data) = self.faces.get(face) else {
                continue;
            };
            let Some(loop_data) = self.loops.get(face_data.outer_loop) else {
                continue;
            };
            let mut entries = Vec::new();
            let mut cursor = loop_data.half_edge;
            for _ in 0..=self.half_edges.len() {
                let Some(he) = self.half_edges.get(cursor) else {
                    break;
                };
                let p = self.vertices[he.origin].position;
                entries.push(format!(
                    "{cursor:?}:e{:?}:v{:?}=({:.4},{:.4},{:.4})",
                    he.edge, he.origin, p.x, p.y, p.z,
                ));
                cursor = he.next;
                if cursor == loop_data.half_edge {
                    break;
                }
            }
            eprintln!(
                "[self-intersection loop] {face:?}: {}",
                entries.join(" -> ")
            );
        }
    }

    /// Check for twisted surfaces: faces where the normal reverses direction.
    ///
    /// Samples the surface normal on a grid and checks for sign changes in
    /// the determinant of the Jacobian (partial_u × partial_v). A sign
    /// change indicates the surface folds over itself within the face.
    pub fn validate_twisted_surfaces(&self) -> Result<(), ValidationError> {
        use crate::geometry::surface::SurfaceEval;

        for (fh, face) in &self.faces {
            let surface = match face.surface.as_ref() {
                Some(s) => &s.geometry,
                None => continue,
            };

            // Get face parameter range
            let positions = self.face_outer_loop_positions(fh);
            if positions.is_empty() {
                continue;
            }
            let ((u_min, u_max), (v_min, v_max)) =
                crate::geometry::uv_projection::tight_bounded_parameter_range(&positions, surface);

            // Sample normals on a 4x4 grid
            let n = 4;
            let mut first_sign = None;
            let mut twisted = false;

            for i in 0..n {
                let u = u_min + (u_max - u_min) * (i as f64 + 0.5) / (n as f64);
                for j in 0..n {
                    let v = v_min + (v_max - v_min) * (j as f64 + 0.5) / (n as f64);
                    let normal = surface.normal_at(u, v);
                    let len = normal.length();
                    if len < 1e-12 {
                        continue; // degenerate normal, skip
                    }

                    // Use z-component sign as a proxy for orientation
                    // (for general surfaces, compare against first sample's direction)
                    let sign = match &first_sign {
                        None => {
                            first_sign = Some(normal);
                            true
                        }
                        Some(first) => {
                            let dot = normal.dot(first);
                            dot > 0.0
                        }
                    };

                    if first_sign.is_some() && !sign {
                        twisted = true;
                        break;
                    }
                }
                if twisted {
                    break;
                }
            }

            if twisted {
                return Err(ValidationError::TwistedSurface { face: fh });
            }
        }

        Ok(())
    }

    /// §5.3 item #10 — check every Gregory-surface face for
    /// within-face self-intersection (fold-over or curl-around).
    /// Delegates per-face to
    /// [`GregoryPatch::self_intersects_with_tolerance`].
    pub fn validate_gregory_self_intersection(&self) -> Result<(), ValidationError> {
        use crate::geometry::surface::SurfaceGeometry;

        let tol_linear = self.tolerance.linear().max(1e-6);
        for (fh, face) in &self.faces {
            let Some(sof) = face.surface.as_ref() else {
                continue;
            };
            let SurfaceGeometry::Gregory(ref patch) = sof.geometry else {
                continue;
            };
            if patch.self_intersects_with_tolerance(tol_linear) {
                return Err(ValidationError::TwistedSurface { face: fh });
            }
        }
        Ok(())
    }
}

#[allow(clippy::too_many_arguments)]
fn cylinder_cylinder_near_samples(
    surface_a: &crate::geometry::surface::SurfaceGeometry,
    u_range_a: (f64, f64),
    v_range_a: (f64, f64),
    surface_b: &crate::geometry::surface::SurfaceGeometry,
    u_range_b: (f64, f64),
    v_range_b: (f64, f64),
    linear_tol: f64,
) -> Option<Vec<Point3>> {
    use crate::geometry::surface::SurfaceGeometry;

    let SurfaceGeometry::Cylinder {
        origin: origin_a,
        axis: axis_a,
        radius: radius_a,
    } = surface_a
    else {
        return None;
    };
    let SurfaceGeometry::Cylinder {
        origin: origin_b,
        axis: axis_b,
        radius: radius_b,
    } = surface_b
    else {
        return None;
    };

    const GRID: usize = 17;
    let threshold = (linear_tol * 10.0).max(1.0e-8);
    let mut samples = Vec::new();
    // Validation containment for unsupported skew cylinder-cylinder pairs:
    // use bounded implicit-distance sampling plus trim classification instead
    // of the generic marcher, whose 4-D subdivision can be disproportionate
    // for small fillet cylinders that only meet at trimmed boundaries.
    sample_cylinder_against_cylinder(
        surface_a,
        u_range_a,
        v_range_a,
        *origin_b,
        *axis_b,
        *radius_b,
        threshold,
        &mut samples,
    );
    sample_cylinder_against_cylinder(
        surface_b,
        u_range_b,
        v_range_b,
        *origin_a,
        *axis_a,
        *radius_a,
        threshold,
        &mut samples,
    );

    if samples.is_empty() {
        for i in 0..GRID {
            let alpha = i as f64 / (GRID - 1) as f64;
            for j in 0..GRID {
                let beta = j as f64 / (GRID - 1) as f64;
                let u = u_range_a.0 + (u_range_a.1 - u_range_a.0) * alpha;
                let v = v_range_a.0 + (v_range_a.1 - v_range_a.0) * beta;
                let p = surface_a.point_at(u, v);
                let distance = distance_to_cylinder(p, *origin_b, *axis_b, *radius_b);
                if distance <= threshold {
                    samples.push(p);
                }
            }
        }
    }

    Some(samples)
}

fn sample_cylinder_against_cylinder(
    sampled_surface: &crate::geometry::surface::SurfaceGeometry,
    u_range: (f64, f64),
    v_range: (f64, f64),
    other_origin: Point3,
    other_axis: crate::geometry::Vec3,
    other_radius: f64,
    threshold: f64,
    samples: &mut Vec<Point3>,
) {
    const GRID: usize = 17;
    for i in 0..GRID {
        let alpha = i as f64 / (GRID - 1) as f64;
        for j in 0..GRID {
            let beta = j as f64 / (GRID - 1) as f64;
            let u = u_range.0 + (u_range.1 - u_range.0) * alpha;
            let v = v_range.0 + (v_range.1 - v_range.0) * beta;
            let p = sampled_surface.point_at(u, v);
            let distance = distance_to_cylinder(p, other_origin, other_axis, other_radius);
            if distance <= threshold {
                samples.push(p);
            }
        }
    }
}

fn distance_to_cylinder(
    point: Point3,
    origin: Point3,
    axis: crate::geometry::Vec3,
    radius: f64,
) -> f64 {
    let offset = point - origin;
    let axial = offset.dot(&axis);
    let radial = offset - axis * axial;
    (radial.length() - radius).abs()
}

fn trace_surface_kind(brep: &BRep, face: crate::topology::FaceHandle) -> &'static str {
    use crate::geometry::surface::SurfaceGeometry;

    match brep
        .faces
        .get(face)
        .and_then(|f| f.surface.as_ref())
        .map(|s| &s.geometry)
    {
        Some(SurfaceGeometry::Plane { .. }) => "Plane",
        Some(SurfaceGeometry::Cylinder { .. }) => "Cylinder",
        Some(SurfaceGeometry::Cone { .. }) => "Cone",
        Some(SurfaceGeometry::Sphere { .. }) => "Sphere",
        Some(SurfaceGeometry::Torus { .. }) => "Torus",
        Some(SurfaceGeometry::Nurbs(_)) => "Nurbs",
        Some(SurfaceGeometry::Gregory(_)) => "Gregory",
        None => "None",
    }
}

fn trace_self_intersection_candidate(
    brep: &BRep,
    face_a: crate::topology::FaceHandle,
    face_b: crate::topology::FaceHandle,
    index: usize,
    u_range_a: (f64, f64),
    v_range_a: (f64, f64),
    u_range_b: (f64, f64),
    v_range_b: (f64, f64),
) {
    if std::env::var_os("POPKERN_SELF_INTERSECTION_TRACE").is_none() {
        return;
    }
    eprintln!(
        "[self-intersection candidate] #{index} {face_a:?}/{} u=({:.6},{:.6}) v=({:.6},{:.6}) \
         {face_b:?}/{} u=({:.6},{:.6}) v=({:.6},{:.6})",
        trace_surface_kind(brep, face_a),
        u_range_a.0,
        u_range_a.1,
        v_range_a.0,
        v_range_a.1,
        trace_surface_kind(brep, face_b),
        u_range_b.0,
        u_range_b.1,
        v_range_b.0,
        v_range_b.1,
    );
    if let Some(detail) = trace_cylinder_detail(brep, face_a) {
        eprintln!("[self-intersection candidate]   {face_a:?} {detail}");
    }
    if let Some(detail) = trace_cylinder_detail(brep, face_b) {
        eprintln!("[self-intersection candidate]   {face_b:?} {detail}");
    }
}

fn trace_cylinder_detail(brep: &BRep, face: crate::topology::FaceHandle) -> Option<String> {
    use crate::geometry::surface::SurfaceGeometry;

    match &brep.faces.get(face)?.surface.as_ref()?.geometry {
        SurfaceGeometry::Cylinder {
            origin,
            axis,
            radius,
        } => Some(format!(
            "Cylinder origin=({:.6},{:.6},{:.6}) axis=({:.6},{:.6},{:.6}) r={:.6}",
            origin.x, origin.y, origin.z, axis.x, axis.y, axis.z, radius
        )),
        _ => None,
    }
}

#[cfg(test)]
#[path = "self_intersection_tests.rs"]
mod tests;
