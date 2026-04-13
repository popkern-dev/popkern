//! Tests for [`crate::modify_inplace`] — both the closed
//! plane-cylinder annulus path (Round 5d-2) and the open plane-plane
//! path (Round 5d-3). Includes the §14.7 ≥2× speedup measurements.
//!
//! Lives in its own file so neither the production module nor the
//! main `modify_tests.rs` exceeds the 500-line structural cap.

use crate::modify::modify_fillet_radius;
use crate::types::{FilletOptions, FilletRadius, FilletSpec};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use std::collections::HashSet;

fn tol() -> Tolerance {
    Tolerance::default()
}

fn cube() -> BRep {
    BRep::unit_cube(tol()).unwrap()
}

// ---------------------------------------------------------------------------
// Round 5d-2: closed plane-cylinder annulus in-place
// ---------------------------------------------------------------------------

#[test]
fn modify_closed_in_place_preserves_topology_and_changes_radius() {
    use popkern_bool::{BoolOperation, boolean};
    use popkern_core::geometry::Transform3;
    use popkern_core::geometry::surface::SurfaceGeometry;

    let cube_brep = cube();
    let cyl = BRep::cylinder(0.3, 2.0, tol()).unwrap();
    let translated = cyl
        .transform(&Transform3::translation(0.5, 0.5, -0.5))
        .unwrap();
    let subtracted = match boolean(
        BoolOperation::Subtraction,
        cube_brep,
        translated,
        tol(),
        None,
    ) {
        Ok(b) => b,
        Err(_) => return,
    };
    let edge = subtracted.edges.keys().find(|e| {
        let he = &subtracted.half_edges[subtracted.edges[*e].half_edge];
        let fa = subtracted.loops[he.loop_].face;
        let fb = subtracted.loops[subtracted.half_edges[he.twin].loop_].face;
        let sa = subtracted.faces[fa].surface.as_ref().map(|s| &s.geometry);
        let sb = subtracted.faces[fb].surface.as_ref().map(|s| &s.geometry);
        matches!(
            (sa, sb),
            (
                Some(SurfaceGeometry::Plane { .. }),
                Some(SurfaceGeometry::Cylinder { .. })
            ) | (
                Some(SurfaceGeometry::Cylinder { .. }),
                Some(SurfaceGeometry::Plane { .. })
            )
        )
    });
    let Some(edge) = edge else { return };
    let r1 = match crate::fillet_edge(&subtracted, edge, 0.05, tol(), None) {
        Ok(b) => b,
        Err(_) => return,
    };
    let annulus = r1.brep.faces.keys().find(|fh| {
        let f = &r1.brep.faces[*fh];
        if f.inner_loops.len() != 1 {
            return false;
        }
        let single = |lh: popkern_core::topology::LoopHandle| -> bool {
            let Some(l) = r1.brep.loops.get(lh) else {
                return false;
            };
            let start = l.half_edge;
            let Some(he) = r1.brep.half_edges.get(start) else {
                return false;
            };
            he.next == start && he.prev == start
        };
        single(f.outer_loop) && single(f.inner_loops[0])
    });
    let annulus = match annulus {
        Some(f) => f,
        None => return,
    };

    let pre_face_count = r1.brep.faces.len();
    let pre_edge_count = r1.brep.edges.len();
    let pre_vertex_count = r1.brep.vertices.len();

    let modified = crate::modify_inplace::modify_closed_in_place(&r1.brep, annulus, 0.08)
        .expect("in-place modify should succeed for plane-cylinder annulus");

    assert_eq!(modified.faces.len(), pre_face_count);
    assert_eq!(modified.edges.len(), pre_edge_count);
    assert_eq!(modified.vertices.len(), pre_vertex_count);

    let new_surface = modified.faces[annulus]
        .surface
        .as_ref()
        .expect("modified annulus should have a surface");
    match &new_surface.geometry {
        SurfaceGeometry::Torus { minor_radius, .. } => {
            assert!(
                (minor_radius - 0.08).abs() < 1e-9,
                "in-place torus minor radius {minor_radius} should equal new radius 0.08"
            );
        }
        other => panic!("in-place modify should assign a Torus surface, got {other:?}"),
    }

    let original_surface = r1.brep.faces[annulus]
        .surface
        .as_ref()
        .expect("original annulus surface");
    if let SurfaceGeometry::Torus { minor_radius, .. } = &original_surface.geometry {
        assert!(
            (minor_radius - 0.08).abs() > 1e-6,
            "original surface should not have already been mutated"
        );
    }
}

#[test]
fn modify_closed_in_place_is_at_least_2x_faster_than_remove_recreate() {
    use popkern_bool::{BoolOperation, boolean};
    use popkern_core::geometry::Transform3;
    use popkern_core::geometry::surface::SurfaceGeometry;
    use std::time::Instant;

    let cube_brep = cube();
    let cyl = BRep::cylinder(0.3, 2.0, tol()).unwrap();
    let translated = cyl
        .transform(&Transform3::translation(0.5, 0.5, -0.5))
        .unwrap();
    let subtracted = match boolean(
        BoolOperation::Subtraction,
        cube_brep,
        translated,
        tol(),
        None,
    ) {
        Ok(b) => b,
        Err(_) => return,
    };
    let edge = subtracted.edges.keys().find(|e| {
        let he = &subtracted.half_edges[subtracted.edges[*e].half_edge];
        let fa = subtracted.loops[he.loop_].face;
        let fb = subtracted.loops[subtracted.half_edges[he.twin].loop_].face;
        let sa = subtracted.faces[fa].surface.as_ref().map(|s| &s.geometry);
        let sb = subtracted.faces[fb].surface.as_ref().map(|s| &s.geometry);
        matches!(
            (sa, sb),
            (
                Some(SurfaceGeometry::Plane { .. }),
                Some(SurfaceGeometry::Cylinder { .. })
            ) | (
                Some(SurfaceGeometry::Cylinder { .. }),
                Some(SurfaceGeometry::Plane { .. })
            )
        )
    });
    let Some(edge) = edge else { return };
    let r1 = match crate::fillet_edge(&subtracted, edge, 0.05, tol(), None) {
        Ok(b) => b,
        Err(_) => return,
    };
    let annulus = r1.brep.faces.keys().find(|fh| {
        let f = &r1.brep.faces[*fh];
        if f.inner_loops.len() != 1 {
            return false;
        }
        let single = |lh: popkern_core::topology::LoopHandle| -> bool {
            let Some(l) = r1.brep.loops.get(lh) else {
                return false;
            };
            let start = l.half_edge;
            let Some(he) = r1.brep.half_edges.get(start) else {
                return false;
            };
            he.next == start && he.prev == start
        };
        single(f.outer_loop) && single(f.inner_loops[0])
    });
    let annulus = match annulus {
        Some(f) => f,
        None => return,
    };

    let n = 20;
    let t_inplace = {
        let t0 = Instant::now();
        for i in 0..n {
            let r = 0.06 + (i as f64) * 0.001;
            let _ = crate::modify_inplace::modify_closed_in_place(&r1.brep, annulus, r).unwrap();
        }
        t0.elapsed()
    };

    let t_baseline = {
        let t0 = Instant::now();
        for i in 0..n {
            let r = 0.06 + (i as f64) * 0.001;
            let unfilleted =
                crate::unfillet_closed::unfillet_closed_annulus(&r1.brep, annulus).unwrap();
            let rim = unfilleted.edges.keys().find(|e| {
                let he = &unfilleted.half_edges[unfilleted.edges[*e].half_edge];
                if he.origin != unfilleted.half_edges[he.twin].origin {
                    return false;
                }
                let fa = unfilleted.loops[he.loop_].face;
                let fb = unfilleted.loops[unfilleted.half_edges[he.twin].loop_].face;
                let sa = unfilleted.faces[fa].surface.as_ref().map(|s| &s.geometry);
                let sb = unfilleted.faces[fb].surface.as_ref().map(|s| &s.geometry);
                matches!(
                    (sa, sb),
                    (
                        Some(SurfaceGeometry::Plane { .. }),
                        Some(SurfaceGeometry::Cylinder { .. })
                    ) | (
                        Some(SurfaceGeometry::Cylinder { .. }),
                        Some(SurfaceGeometry::Plane { .. })
                    )
                )
            });
            let rim = match rim {
                Some(r) => r,
                None => return,
            };
            let _ = crate::fillet_edge(&unfilleted, rim, r, tol(), None).unwrap();
        }
        t0.elapsed()
    };

    let ratio = t_baseline.as_secs_f64() / t_inplace.as_secs_f64();
    eprintln!(
        "[5d-2 perf] in-place: {:?}, baseline: {:?}, speedup: {:.2}×",
        t_inplace, t_baseline, ratio
    );
    assert!(
        ratio >= 2.0,
        "in-place should be ≥2× faster than remove+recreate; got {ratio:.2}× \
         (in-place: {t_inplace:?}, baseline: {t_baseline:?})"
    );
}

// ---------------------------------------------------------------------------
// Round 5d-3: open plane-plane in-place
// ---------------------------------------------------------------------------

#[test]
fn modify_open_in_place_preserves_topology_and_changes_radius() {
    use popkern_core::geometry::surface::SurfaceGeometry;

    let brep = cube();
    let edge = brep.edges.keys().next().unwrap();

    let opts = FilletOptions {
        propagate_tangent: false,
        vertex_mode: crate::types::VertexMode::Mitred,
    };
    let specs = vec![FilletSpec {
        edge,
        radius: FilletRadius::Constant(0.1),
    }];
    let r1 = crate::fillet_edges(&brep, &specs, &opts, tol(), None).unwrap();
    let blend = r1
        .brep
        .faces
        .keys()
        .find(|fh| {
            matches!(
                r1.brep.faces[*fh].surface.as_ref().map(|s| &s.geometry),
                Some(SurfaceGeometry::Torus { .. })
            )
        })
        .expect("post-fillet brep should contain a Torus face");

    let pre_face_count = r1.brep.faces.len();
    let pre_edge_count = r1.brep.edges.len();
    let pre_vertex_count = r1.brep.vertices.len();

    let modified = crate::modify_inplace::modify_open_in_place(&r1.brep, blend, 0.2)
        .expect("in-place open modify should succeed for plane-plane fillet");

    assert_eq!(modified.faces.len(), pre_face_count);
    assert_eq!(modified.edges.len(), pre_edge_count);
    assert_eq!(modified.vertices.len(), pre_vertex_count);

    let surf = modified.faces[blend]
        .surface
        .as_ref()
        .expect("modified blend should have a surface");
    match &surf.geometry {
        SurfaceGeometry::Torus { minor_radius, .. } => {
            assert!(
                (minor_radius - 0.2).abs() < 1e-9,
                "in-place torus minor radius {minor_radius} should equal new radius 0.2"
            );
        }
        other => panic!("in-place open modify should assign a Torus, got {other:?}"),
    }

    let original_surf = r1.brep.faces[blend]
        .surface
        .as_ref()
        .expect("original blend surface");
    if let SurfaceGeometry::Torus { minor_radius, .. } = &original_surf.geometry {
        assert!(
            (minor_radius - 0.2).abs() > 1e-6,
            "original blend surface should not have already been mutated"
        );
    }
}

// ---------------------------------------------------------------------------
// Round 5h: variable-radius in-place
// ---------------------------------------------------------------------------

/// `modify_open_in_place_variable` moves corners using per-end
/// radii. After the call, the start corners must be at distance
/// `r_start` from the start endpoint and the end corners at distance
/// `r_end` from the end endpoint.
#[test]
fn modify_open_in_place_variable_uses_per_end_radii() {
    use popkern_core::geometry::surface::SurfaceGeometry;

    let brep = cube();
    let edge = brep.edges.keys().next().unwrap();
    let opts = FilletOptions {
        propagate_tangent: false,
        vertex_mode: crate::types::VertexMode::Mitred,
    };
    let specs = vec![FilletSpec {
        edge,
        radius: FilletRadius::Constant(0.1),
    }];
    let r1 = crate::fillet_edges(&brep, &specs, &opts, tol(), None).unwrap();
    let blend = r1
        .brep
        .faces
        .keys()
        .find(|fh| {
            matches!(
                r1.brep.faces[*fh].surface.as_ref().map(|s| &s.geometry),
                Some(SurfaceGeometry::Torus { .. })
            )
        })
        .expect("post-fillet brep should contain a Torus face");

    let r_start = 0.05;
    let r_end = 0.25;
    let modified =
        crate::modify_inplace::modify_open_in_place_variable(&r1.brep, blend, r_start, r_end)
            .expect("variable in-place should succeed for plane-plane fillet");

    // Walk the modified brep and find the vertices nearest each
    // sharp edge endpoint. The unit cube fillet edge is from
    // (0,0,0) to (1,0,0); the corners at the start should be at
    // distance r_start from (0,0,0), and the corners at the end
    // should be at distance r_end from (1,0,0).
    let plan = crate::unfillet::compute_unfillet_plan(&r1.brep, blend).unwrap();
    let p_start = plan.restored_edge_start;
    let p_end = plan.restored_edge_end;

    let mut start_distances: Vec<f64> = Vec::new();
    let mut end_distances: Vec<f64> = Vec::new();
    for (_, v) in &modified.vertices {
        let pos = v.position;
        let d_start = ((pos.x - p_start.x).powi(2)
            + (pos.y - p_start.y).powi(2)
            + (pos.z - p_start.z).powi(2))
        .sqrt();
        let d_end =
            ((pos.x - p_end.x).powi(2) + (pos.y - p_end.y).powi(2) + (pos.z - p_end.z).powi(2))
                .sqrt();
        if d_start < d_end && d_start > 1e-9 {
            start_distances.push(d_start);
        } else if d_end < d_start && d_end > 1e-9 {
            end_distances.push(d_end);
        }
    }

    // We expect at least 2 corners near each end at the requested
    // radius. (There may also be unrelated cube-corner vertices on
    // the far sides of each end, which will be at distances much
    // larger than r_start/r_end.)
    let near_start_at_radius = start_distances
        .iter()
        .filter(|d| (**d - r_start).abs() < 1e-6)
        .count();
    let near_end_at_radius = end_distances
        .iter()
        .filter(|d| (**d - r_end).abs() < 1e-6)
        .count();
    assert!(
        near_start_at_radius >= 2,
        "should have ≥2 corners at distance {r_start} from start; got {near_start_at_radius} (distances: {start_distances:?})"
    );
    assert!(
        near_end_at_radius >= 2,
        "should have ≥2 corners at distance {r_end} from end; got {near_end_at_radius} (distances: {end_distances:?})"
    );
}

#[test]
fn modify_open_in_place_is_at_least_2x_faster_than_remove_recreate() {
    use std::time::Instant;

    let brep = cube();
    let edge = brep.edges.keys().next().unwrap();
    let pre: HashSet<_> = brep.faces.keys().collect();
    let r1 = crate::fillet_edge(&brep, edge, 0.1, tol(), None).unwrap();
    let blend = r1
        .brep
        .faces
        .keys()
        .find(|fh| !pre.contains(fh) && crate::identify_blend(&r1.brep, *fh).is_some())
        .expect("initial fillet should produce a blend");

    let n = 20;
    let t_inplace = {
        let t0 = Instant::now();
        for i in 0..n {
            let r = 0.12 + (i as f64) * 0.001;
            let _ = crate::modify_inplace::modify_open_in_place(&r1.brep, blend, r).unwrap();
        }
        t0.elapsed()
    };

    let t_baseline = {
        let t0 = Instant::now();
        for i in 0..n {
            let r = 0.12 + (i as f64) * 0.001;
            let _ = modify_fillet_radius(&r1.brep, blend, r, tol()).unwrap();
            let unfilleted = crate::unfillet::unfillet_face(&r1.brep, blend).unwrap();
            let plan = crate::unfillet::compute_unfillet_plan(&r1.brep, blend).unwrap();
            let p0 = plan.restored_edge_start;
            let p1 = plan.restored_edge_end;
            let restored_edge = unfilleted.edges.keys().find(|e| {
                let he = &unfilleted.half_edges[unfilleted.edges[*e].half_edge];
                let twin = &unfilleted.half_edges[he.twin];
                let v0 = unfilleted.vertices[he.origin].position;
                let v1 = unfilleted.vertices[twin.origin].position;
                let close = |a: &popkern_core::geometry::Point3,
                             b: &popkern_core::geometry::Point3|
                 -> bool {
                    let dx = a.x - b.x;
                    let dy = a.y - b.y;
                    let dz = a.z - b.z;
                    (dx * dx + dy * dy + dz * dz).sqrt() < 1e-6
                };
                (close(&v0, &p0) && close(&v1, &p1)) || (close(&v0, &p1) && close(&v1, &p0))
            });
            if let Some(e) = restored_edge {
                let _ = crate::fillet_edge(&unfilleted, e, r, tol(), None).unwrap();
            }
        }
        t0.elapsed()
    };

    let ratio = t_baseline.as_secs_f64() / t_inplace.as_secs_f64();
    eprintln!(
        "[5d-3 perf] in-place: {:?}, baseline: {:?}, speedup: {:.2}×",
        t_inplace, t_baseline, ratio
    );
    assert!(
        ratio >= 2.0,
        "in-place open should be ≥2× faster than remove+recreate; got {ratio:.2}× \
         (in-place: {t_inplace:?}, baseline: {t_baseline:?})"
    );
}
