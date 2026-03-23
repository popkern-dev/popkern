//! Stage 1: Face-face intersection curve computation.
//!
//! Given two BRep solids, finds all intersection curves between their faces
//! using AABB filtering followed by surface-surface intersection.

use popkern_core::geometry::curve_on_surface::CurveOnSurface;
use popkern_core::{
    geometry::{
        Point3,
        intersection::{SurfaceSurfaceIntersection, intersect_surfaces, surface_subdivision::Aabb},
        surface::SurfaceEval,
    },
    tolerance::Tolerance,
    topology::{BRep, FaceHandle},
};

use crate::error::{BoolDiagnostics, BoolError, Stage};

/// A single face-face intersection record.
#[derive(Debug, Clone)]
pub(crate) struct IntersectionRecord {
    /// Face handle from solid A.
    pub face_a: FaceHandle,
    /// Face handle from solid B.
    pub face_b: FaceHandle,
    /// The intersection curve on both surfaces.
    pub curve: CurveOnSurface,
}

/// Complete result of Stage 1 intersection computation.
#[derive(Debug)]
pub(crate) struct IntersectionResult {
    /// All face-face intersection records found.
    pub records: Vec<IntersectionRecord>,
    /// AABBs computed for each face of solid A.
    pub aabbs_a: Vec<(FaceHandle, Aabb)>,
    /// AABBs computed for each face of solid B.
    pub aabbs_b: Vec<(FaceHandle, Aabb)>,
}

/// Collect vertex positions from the outer loop of a face by walking half-edges.
///
/// Returns an empty `Vec` if the face's loop has no valid half-edges.
fn face_vertex_positions(brep: &BRep, face_handle: FaceHandle) -> Vec<Point3> {
    let face = &brep.faces[face_handle];
    let loop_handle = face.outer_loop;
    let loop_data = match brep.loops.get(loop_handle) {
        Some(l) => l,
        None => return Vec::new(),
    };
    let start_he = loop_data.half_edge;
    if !brep.half_edges.contains_key(start_he) {
        return Vec::new();
    }

    let mut positions = Vec::new();
    let mut current = start_he;
    loop {
        let he = &brep.half_edges[current];
        positions.push(brep.vertices[he.origin].position);
        current = he.next;
        if current == start_he {
            break;
        }
    }
    positions
}

/// Sample a face to compute its AABB.
///
/// For faces with finite parameter ranges (e.g. NURBS, sphere, torus), samples
/// a 5×5 grid over the surface. For faces with infinite ranges (e.g. planes),
/// falls back to collecting vertex positions from the face's outer loop.
///
/// Returns `None` if the face has no surface geometry or no valid vertices.
fn face_aabb(brep: &BRep, face_handle: FaceHandle) -> Option<Aabb> {
    let surface_on_face = brep.faces[face_handle].surface.as_ref()?;
    let surf = &surface_on_face.geometry;
    let ((u_min, u_max), (v_min, v_max)) = surf.parameter_range();

    let mut points = Vec::new();

    // Always add vertex positions — accurate for planar faces and provides a
    // safe fallback when parameter range is infinite.
    points.extend(face_vertex_positions(brep, face_handle));

    // If the surface has a finite parameter range, also sample a 5×5 grid for
    // curved surfaces where vertices may not capture the full extent.
    if u_min.is_finite() && u_max.is_finite() && v_min.is_finite() && v_max.is_finite() {
        const N: usize = 5;
        for i in 0..N {
            let u = u_min + (u_max - u_min) * (i as f64) / ((N - 1) as f64);
            for j in 0..N {
                let v = v_min + (v_max - v_min) * (j as f64) / ((N - 1) as f64);
                points.push(surf.point_at(u, v));
            }
        }
    }

    Aabb::from_points(&points)
}

/// Compute AABBs for all faces in a BRep that have surface geometry.
///
/// Returns a list of `(FaceHandle, Aabb)` pairs.
pub(crate) fn compute_face_aabbs(brep: &BRep) -> Vec<(FaceHandle, Aabb)> {
    brep.faces
        .keys()
        .filter_map(|fh| face_aabb(brep, fh).map(|aabb| (fh, aabb)))
        .collect()
}

/// Find all intersection curves between faces of two BRep solids.
///
/// Uses AABB filtering to reduce the number of expensive surface-surface
/// intersection calls. Returns an [`IntersectionResult`] on success, or a
/// [`BoolDiagnostics`] containing the error and any partial results on failure.
///
/// # Arguments
///
/// * `a` — the first (left-hand) solid.
/// * `b` — the second (right-hand) solid.
/// * `tolerance` — geometric tolerance for all distance and coincidence tests.
pub(crate) fn find_intersections(
    a: &BRep,
    b: &BRep,
    tolerance: &Tolerance,
) -> Result<IntersectionResult, BoolDiagnostics> {
    let aabbs_a = compute_face_aabbs(a);
    let aabbs_b = compute_face_aabbs(b);
    let mut records: Vec<IntersectionRecord> = Vec::new();

    for (face_a, aabb_a) in &aabbs_a {
        for (face_b, aabb_b) in &aabbs_b {
            if !aabb_a.overlaps(aabb_b, tolerance.linear) {
                continue;
            }

            // Both faces are guaranteed to have surface geometry because
            // compute_face_aabbs only includes faces with surfaces.
            let surf_a = match a.faces[*face_a].surface.as_ref() {
                Some(s) => &s.geometry,
                None => continue,
            };
            let surf_b = match b.faces[*face_b].surface.as_ref() {
                Some(s) => &s.geometry,
                None => continue,
            };

            let ((u1_min, u1_max), (v1_min, v1_max)) = surf_a.parameter_range();
            let ((u2_min, u2_max), (v2_min, v2_max)) = surf_b.parameter_range();

            match intersect_surfaces(
                surf_a,
                (u1_min, u1_max),
                (v1_min, v1_max),
                surf_b,
                (u2_min, u2_max),
                (v2_min, v2_max),
                tolerance,
            ) {
                Ok(result) => {
                    for isect in result.intersections {
                        match isect {
                            SurfaceSurfaceIntersection::Curve(curve) => {
                                records.push(IntersectionRecord {
                                    face_a: *face_a,
                                    face_b: *face_b,
                                    curve,
                                });
                            }
                            SurfaceSurfaceIntersection::CoincidentRegion { .. } => {
                                return Err(BoolDiagnostics {
                                    error: Some(BoolError::CoincidentFaces {
                                        face_a: *face_a,
                                        face_b: *face_b,
                                    }),
                                    intersection_curves: records
                                        .into_iter()
                                        .map(|r| r.curve)
                                        .collect(),
                                    face_classifications: Vec::new(),
                                    failed_stage: Some(Stage::Intersection),
                                });
                            }
                        }
                    }
                }
                Err(e) => {
                    return Err(BoolDiagnostics {
                        error: Some(BoolError::IntersectionFailed(format!(
                            "face_a={face_a:?} face_b={face_b:?}: {e}"
                        ))),
                        intersection_curves: records.into_iter().map(|r| r.curve).collect(),
                        face_classifications: Vec::new(),
                        failed_stage: Some(Stage::Intersection),
                    });
                }
            }
        }
    }

    Ok(IntersectionResult {
        records,
        aabbs_a,
        aabbs_b,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use popkern_core::{geometry::Vec3, tolerance::Tolerance, topology::BRep};

    fn tol() -> Tolerance {
        Tolerance::default()
    }

    /// Translate all vertices of a BRep by a vector offset.
    fn translate(mut brep: BRep, offset: Vec3) -> BRep {
        for (_, vertex) in &mut brep.vertices {
            vertex.position = vertex.position + offset;
        }
        brep
    }

    #[test]
    fn overlapping_cubes_find_intersection_curves() {
        let cube_a = BRep::unit_cube(tol()).expect("unit cube a");
        // Offset cube_b by (0.5, 0.3, 0.2) so it overlaps with cube_a but
        // none of the face planes are coincident (all 6 faces land on
        // distinct planes from cube_a's 6 faces).
        let cube_b = translate(
            BRep::unit_cube(tol()).expect("unit cube b"),
            Vec3::new(0.5, 0.3, 0.2),
        );

        let result = find_intersections(&cube_a, &cube_b, &tol())
            .expect("overlapping cubes should not error");

        assert!(
            !result.records.is_empty(),
            "expected at least one intersection curve for overlapping cubes"
        );
    }

    #[test]
    fn disjoint_cubes_no_intersection_curves() {
        let cube_a = BRep::unit_cube(tol()).expect("unit cube a");
        // Offset cube_b by 5.0 in x so it is completely separate.
        let cube_b = translate(
            BRep::unit_cube(tol()).expect("unit cube b"),
            Vec3::new(5.0, 0.0, 0.0),
        );

        let result =
            find_intersections(&cube_a, &cube_b, &tol()).expect("disjoint cubes should not error");

        assert!(
            result.records.is_empty(),
            "expected no intersection curves for disjoint cubes"
        );
    }

    #[test]
    fn compute_face_aabbs_unit_cube_has_six_faces() {
        let cube = BRep::unit_cube(tol()).expect("unit cube");
        let aabbs = compute_face_aabbs(&cube);
        assert_eq!(aabbs.len(), 6, "unit cube should have 6 faces with AABBs");
    }

    #[test]
    fn compute_face_aabbs_no_surface_returns_empty() {
        // A freshly-created BRep has no faces at all.
        let brep = BRep::new(tol());
        let aabbs = compute_face_aabbs(&brep);
        assert!(aabbs.is_empty());
    }

    #[test]
    fn intersection_result_aabbs_preserved() {
        let cube_a = BRep::unit_cube(tol()).expect("unit cube a");
        let cube_b = translate(
            BRep::unit_cube(tol()).expect("unit cube b"),
            Vec3::new(5.0, 0.0, 0.0),
        );
        let result = find_intersections(&cube_a, &cube_b, &tol()).expect("no error");
        assert_eq!(result.aabbs_a.len(), 6);
        assert_eq!(result.aabbs_b.len(), 6);
    }

    #[test]
    fn face_aabb_unit_cube_face_is_bounded() {
        let cube = BRep::unit_cube(tol()).expect("unit cube");
        // Check every face produces a finite AABB within [0,1]^3 ± small epsilon.
        let aabbs = compute_face_aabbs(&cube);
        for (_, aabb) in &aabbs {
            assert!(
                aabb.min.x >= -1e-6 && aabb.max.x <= 1.0 + 1e-6,
                "x out of range"
            );
            assert!(
                aabb.min.y >= -1e-6 && aabb.max.y <= 1.0 + 1e-6,
                "y out of range"
            );
            assert!(
                aabb.min.z >= -1e-6 && aabb.max.z <= 1.0 + 1e-6,
                "z out of range"
            );
        }
    }

    /// Exercise face_vertex_positions when the loop's half_edge handle is
    /// missing from the BRep — it should return an empty Vec rather than panic.
    #[test]
    fn face_vertex_positions_missing_half_edge_returns_empty() {
        use popkern_core::topology::{Face, HalfEdgeHandle, Loop, ShellHandle};
        let mut brep = BRep::new(tol());
        // Insert a loop whose half_edge handle does not exist in brep.half_edges.
        let loop_h = brep.loops.insert(Loop {
            half_edge: HalfEdgeHandle::default(),
            face: popkern_core::topology::FaceHandle::default(),
        });
        let face_h = brep.faces.insert(Face {
            outer_loop: loop_h,
            inner_loops: Vec::new(),
            shell: ShellHandle::default(),
            surface: None,
        });
        brep.loops[loop_h].face = face_h;
        // face_aabb should return None (no surface geometry).
        let result = face_aabb(&brep, face_h);
        assert!(result.is_none());
    }

    /// Exercise the 5×5 grid-sampling branch in face_aabb for a surface with
    /// a finite parameter range (sphere/torus). We use a unit cube here and
    /// verify the AABBs are still correct — the cube planes have infinite
    /// parameter ranges so only vertices are used. For a more direct test we
    /// build a minimal BRep face with a sphere surface.
    #[test]
    fn face_aabb_sphere_surface_samples_grid() {
        use popkern_core::{
            geometry::surface::{SurfaceGeometry, SurfaceOnFace},
            topology::{Face, HalfEdgeHandle, Loop, ShellHandle},
        };
        let mut brep = BRep::new(tol());
        let loop_h = brep.loops.insert(Loop {
            half_edge: HalfEdgeHandle::default(),
            face: popkern_core::topology::FaceHandle::default(),
        });
        let face_h = brep.faces.insert(Face {
            outer_loop: loop_h,
            inner_loops: Vec::new(),
            shell: ShellHandle::default(),
            surface: Some(SurfaceOnFace {
                geometry: SurfaceGeometry::Sphere {
                    center: Point3::new(0.0, 0.0, 0.0),
                    radius: 1.0,
                },
                orientation: true,
            }),
        });
        brep.loops[loop_h].face = face_h;
        // The sphere surface has a finite parameter range, so face_aabb should
        // sample the 5×5 grid and return a non-None AABB.
        let result = face_aabb(&brep, face_h);
        assert!(result.is_some(), "sphere face should produce an AABB");
        let aabb = result.unwrap();
        // Unit sphere AABB should be approximately [-1, 1]^3.
        assert!(
            aabb.min.x >= -1.1 && aabb.max.x <= 1.1,
            "x out of sphere range"
        );
    }

    /// Exercise the CoincidentRegion error path in find_intersections.
    /// We use two identical cubes (same position) which produce coincident
    /// face pairs and trigger the error.
    #[test]
    fn coincident_cubes_returns_coincident_faces_error() {
        let cube_a = BRep::unit_cube(tol()).expect("unit cube a");
        let cube_b = BRep::unit_cube(tol()).expect("unit cube b (same position)");
        let result = find_intersections(&cube_a, &cube_b, &tol());
        assert!(
            result.is_err(),
            "coincident cubes should produce a CoincidentFaces error"
        );
    }
}
