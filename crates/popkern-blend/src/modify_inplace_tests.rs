//! Tests for [`crate::modify_inplace`] — both the closed
//! plane-cylinder annulus path  and the open plane-plane
//! path . Includes the §14.7 ≥2× speedup measurements.
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
// closed plane-cylinder annulus in-place
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
        SurfaceGeometry::Cylinder { radius, .. } => {
            assert!(
                (radius - 0.08).abs() < 1e-9,
                "in-place cylinder radius {radius} should equal new radius 0.08"
            );
        }
        other => panic!("in-place modify should assign a Cylinder surface, got {other:?}"),
    }

    let original_surface = r1.brep.faces[annulus]
        .surface
        .as_ref()
        .expect("original annulus surface");
    if let SurfaceGeometry::Cylinder { radius, .. } = &original_surface.geometry {
        assert!(
            (radius - 0.08).abs() > 1e-6,
            "original surface should not have already been mutated"
        );
    }

    // PK-241: the closed-annulus surface swap must leave the strip
    // flag shell-honest (its rewrite stamps a provisional flag,
    // then reconciles against the support tangency).
    assert_strip_flag_shell_honest(&modified, annulus);
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
// open plane-plane in-place
// ---------------------------------------------------------------------------

#[test]
fn modify_open_in_place_preserves_topology_and_changes_radius() {
    use popkern_core::geometry::surface::SurfaceGeometry;

    let brep = cube();
    let edge = brep.edges.keys().next().unwrap();

    let opts = FilletOptions {
        propagate_tangent: false,
        vertex_mode: crate::types::VertexMode::Mitred,
        ..FilletOptions::default()
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
                Some(SurfaceGeometry::Cylinder { .. })
            )
        })
        .expect("post-fillet brep should contain a Cylinder face");

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
        SurfaceGeometry::Cylinder { radius, .. } => {
            assert!(
                (radius - 0.2).abs() < 1e-9,
                "in-place cylinder radius {radius} should equal new radius 0.2"
            );
        }
        other => panic!("in-place open modify should assign a Cylinder, got {other:?}"),
    }

    let original_surf = r1.brep.faces[blend]
        .surface
        .as_ref()
        .expect("original blend surface");
    if let SurfaceGeometry::Cylinder { radius, .. } = &original_surf.geometry {
        assert!(
            (radius - 0.2).abs() > 1e-6,
            "original blend surface should not have already been mutated"
        );
    }
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

// ---------------------------------------------------------------------------
// PK-241: shell-honest strip orientation flags survive in-place modify
// ---------------------------------------------------------------------------

/// Assert `strip`'s flag-folded natural normal AGREES with a
/// tangent support face's flag-folded natural at a shared boundary
/// sample — the manifold-shell honesty invariant the strip authors
/// must preserve (PK-238/PK-241). Fails when any strip author
/// stamps a provisional flag without reconciling.
fn assert_strip_flag_shell_honest(brep: &BRep, strip: popkern_core::topology::FaceHandle) {
    use popkern_core::geometry::surface::SurfaceEval;
    let strip_sof = brep.faces[strip].surface.as_ref().expect("strip surface");
    let outer = brep.faces[strip].outer_loop;
    let start = brep.loops[outer].half_edge;
    let mut best: Option<f64> = None;
    let mut cur = start;
    for _ in 0..brep.half_edges.len() + 1 {
        let hd = &brep.half_edges[cur];
        let td = &brep.half_edges[hd.twin];
        let other = brep.loops[td.loop_].face;
        if other != strip
            && let Some(sof) = brep.faces[other].surface.as_ref()
        {
            let p0 = brep.vertices[hd.origin].position;
            let p1 = brep.vertices[td.origin].position;
            let mid = popkern_core::geometry::Point3::new(
                (p0.x + p1.x) * 0.5,
                (p0.y + p1.y) * 0.5,
                (p0.z + p1.z) * 0.5,
            );
            let folded = |s: &popkern_core::geometry::surface::SurfaceOnFace| {
                let eval = &s.geometry as &dyn SurfaceEval;
                let (u, v) = eval.closest_point(&mid)?;
                let n = eval.normal_at(u, v).normalized()?;
                Some(n * if s.orientation { 1.0 } else { -1.0 })
            };
            if let (Some(ns), Some(no)) = (folded(strip_sof), folded(sof)) {
                let d = ns.dot(&no);
                if best.is_none_or(|b: f64| d.abs() > b.abs()) {
                    best = Some(d);
                }
            }
        }
        cur = hd.next;
        if cur == start {
            break;
        }
    }
    let witness = best.expect("strip should have a surfaced neighbour");
    assert!(
        witness.abs() > 0.9,
        "no tangency witness on the strip boundary (best |cos| {:.3})",
        witness.abs()
    );
    assert!(
        witness > 0.0,
        "strip orientation flag disagrees with its tangent support (cos {witness:.3}) — \
         a strip author stamped a provisional flag without reconciling (PK-241)"
    );
}

#[test]
fn modify_open_in_place_keeps_strip_flag_shell_honest() {
    use popkern_core::geometry::surface::SurfaceGeometry;

    let brep = cube();
    let edge = brep.edges.keys().next().unwrap();
    let opts = FilletOptions {
        propagate_tangent: false,
        vertex_mode: crate::types::VertexMode::Mitred,
        ..FilletOptions::default()
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
                Some(SurfaceGeometry::Cylinder { .. })
            )
        })
        .expect("post-fillet brep should contain a Cylinder face");
    // Fresh commit is honest (PK-238's primary-author reconcile)…
    assert_strip_flag_shell_honest(&r1.brep, blend);
    // …and the in-place radius swap must preserve it (its surface
    // rewrite stamps a provisional flag, then reconciles).
    let modified = crate::modify_inplace::modify_open_in_place(&r1.brep, blend, 0.2)
        .expect("in-place open modify should succeed");
    assert_strip_flag_shell_honest(&modified, blend);
}
