//! The generative feature on a small cantilever (PK-956): a 12 × 4 × 4
//! mm block clamped on its x = 0 face and pulled down on its x = 12
//! face, 1 mm cells, a few updates: the body is valid and encloses its
//! mesh, the kept face keeps its name, the cache serves an unchanged
//! design, a deferred host declines an unsolved one, and the observer
//! sees and can stop the solve.

use std::ops::ControlFlow;

use popkern_core::geometry::Point3;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::history::TrackingContext;
use popkern_core::naming::{EntityName, EntityRole, OperationId, RoleCategory};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};
use popkern_generative::{
    GenerativeError, GenerativeRequest, LoadGroup, SupportGroup, SupportKind, defer_with, generate,
    generate_cached, name_design, observe_with,
};

fn tol() -> Tolerance {
    Tolerance::default()
}

/// The block, every face named `Generated/i` under operation 1.
fn block() -> BRep {
    let corners: Vec<Point3> = [0.0, 12.0]
        .iter()
        .flat_map(|&x| {
            [0.0, 4.0]
                .into_iter()
                .flat_map(move |y| [0.0, 4.0].into_iter().map(move |z| Point3::new(x, y, z)))
        })
        .collect();
    let part = popkern_part::Part::polyhedron(&corners, tol()).unwrap_or_else(|e| panic!("{e}"));
    let mut brep = part.brep().clone();
    let faces: Vec<FaceHandle> = brep.faces.keys().collect();
    for (i, f) in faces.into_iter().enumerate() {
        let id = brep.name_pool.intern(
            EntityName {
                operation: OperationId(1),
                role: EntityRole {
                    category: RoleCategory::Generated,
                    index: i as u32,
                },
                ancestors: Vec::new(),
            },
            None,
        );
        brep.faces[f].name = Some(id);
    }
    brep
}

/// The planar face whose outward normal is ±x (`sense`) at x = `x`.
fn end_face(brep: &BRep, sense: f64, x: f64) -> FaceHandle {
    brep.faces
        .iter()
        .find(
            |(_, f)| match f.surface.as_ref().map(|s| (&s.geometry, s.orientation)) {
                Some((SurfaceGeometry::Plane { origin, normal, .. }, o)) => {
                    let n = if o { normal.x } else { -normal.x };
                    n * sense > 0.99 && (origin.x - x).abs() < 1e-9
                }
                _ => false,
            },
        )
        .map(|(h, _)| h)
        .unwrap_or_else(|| panic!("no face at x = {x}"))
}

fn request(brep: &BRep, volume_fraction: f64) -> GenerativeRequest {
    let wall = end_face(brep, -1.0, 0.0);
    let tip = end_face(brep, 1.0, 12.0);
    GenerativeRequest {
        material: "PLA (FDM)".into(),
        supports: vec![SupportGroup {
            faces: vec![wall],
            kind: SupportKind::Fixed,
        }],
        loads: vec![LoadGroup {
            faces: vec![tip],
            force_n: 10.0,
            direction: [0.0, 0.0, -1.0],
        }],
        keep: vec![wall, tip],
        volume_fraction,
        filter: 1.5,
        keep_depth: Some(1.0),
        cell: Some(1.0),
        build_direction: None,
        max_iterations: 150,
        robust: false,
        max_cells: None,
        mm_per_unit: 1.0,
    }
}

fn volume(brep: &BRep) -> f64 {
    popkern_measure::body_volume(brep, &tol()).unwrap_or_else(|e| panic!("{e}"))
}

#[test]
fn the_design_is_a_valid_body_that_encloses_its_mesh() {
    let space = block();
    let (result, design) = generate(&space, &request(&space, 0.4), &mut |_| {
        ControlFlow::Continue(())
    })
    .unwrap_or_else(|e| panic!("{e}"));
    result.brep.validate().unwrap_or_else(|e| panic!("{e}"));
    let v = volume(&result.brep);
    let mesh = result.report.facet.mesh_volume;
    assert!((v - mesh).abs() <= 1e-9 * mesh, "body {v} mesh {mesh}");
    // Thresholded at ½, the body holds about the design's material.
    let space_volume = 12.0 * 4.0 * 4.0;
    assert!(
        (v / space_volume - design.volume_fraction).abs() < 0.15,
        "body {} against design {}",
        v / space_volume,
        design.volume_fraction
    );
    assert_eq!(result.sources.len(), result.brep.faces.len());
    let wall = end_face(&space, -1.0, 0.0);
    assert!(
        result.sources.contains(&Some(wall)),
        "the wall face is kept"
    );
    assert!(
        result.sources.contains(&None),
        "the design has its own surface"
    );
    // Faceted: every face is a plane.
    assert!(result.brep.faces.values().all(|f| matches!(
        f.surface.as_ref().map(|s| &s.geometry),
        Some(SurfaceGeometry::Plane { .. })
    )));
}

#[test]
fn the_kept_face_is_named_as_modified_from_the_design_space() {
    let space = block();
    let (result, _) = generate(&space, &request(&space, 0.4), &mut |_| {
        ControlFlow::Continue(())
    })
    .unwrap_or_else(|e| panic!("{e}"));
    let mut body = result.brep.clone();
    let mut tracking = TrackingContext::new(OperationId(7));
    name_design(&space, &mut body, &result.sources, &mut tracking);
    let wall = end_face(&space, -1.0, 0.0);
    let wall_name = space.faces[wall].name.unwrap_or_else(|| panic!("named"));
    let modified: Vec<&EntityName> = body
        .faces
        .values()
        .filter_map(|f| body.name_pool.get(f.name?))
        .filter(|n| n.role.category == RoleCategory::Modified)
        .collect();
    // Every patch on a design space face (a 4 mm block keeps some on
    // every face), each from its source.
    assert_eq!(
        modified.len(),
        result.sources.iter().flatten().count(),
        "{modified:?}"
    );
    assert!(modified.iter().any(|n| n.ancestors == vec![wall_name]));
    assert!(modified.iter().all(|n| n.operation == OperationId(7)));
    // Every entity is named.
    assert!(body.faces.values().all(|f| f.name.is_some()));
    assert!(body.edges.values().all(|e| e.name.is_some()));
    assert!(body.vertices.values().all(|v| v.name.is_some()));
    let history = tracking.finish();
    let old = space
        .name_pool
        .get(wall_name)
        .cloned()
        .unwrap_or_else(|| panic!("pooled"));
    let patches = result.sources.iter().filter(|s| **s == Some(wall)).count();
    assert_eq!(history.modified.get(&old).map(Vec::len), Some(patches));
    // The largest patch keeps the source's place in the design space.
    let ordinal = space.faces.keys().position(|h| h == wall).unwrap_or(0) as u32;
    assert!(
        modified
            .iter()
            .any(|n| n.ancestors == vec![wall_name] && n.role.index == ordinal)
    );
}

#[test]
fn the_cache_serves_an_unchanged_design_and_a_deferred_host_declines_a_new_one() {
    let space = block();
    let first = request(&space, 0.45);
    let (a, hit_a) = generate_cached(&space, &first).unwrap_or_else(|e| panic!("{e}"));
    let (b, hit_b) = generate_cached(&space, &first).unwrap_or_else(|e| panic!("{e}"));
    assert!(!hit_a && hit_b);
    assert_eq!(a.brep.faces.len(), b.brep.faces.len());
    assert_eq!(a.sources, b.sources);
    assert_eq!(volume(&a.brep).to_bits(), volume(&b.brep).to_bits());
    // Deferred: the cached design still serves, a new one declines.
    let served = defer_with(true, || generate_cached(&space, &first));
    assert!(matches!(served, Ok((_, true))));
    let declined = defer_with(true, || generate_cached(&space, &request(&space, 0.35)));
    match declined {
        Err(e @ GenerativeError::NeedsSolve) => assert_eq!(e.kind(), "NeedsSolve"),
        other => panic!("{other:?}"),
    }
}

#[test]
fn the_observer_sees_each_update_and_can_stop_the_solve() {
    let space = block();
    let counts = std::rc::Rc::new(std::cell::RefCell::new(Vec::new()));
    let seen = std::rc::Rc::clone(&counts);
    let solved = observe_with(
        move |u| {
            seen.borrow_mut().push(u.count);
            ControlFlow::Continue(())
        },
        || generate_cached(&space, &request(&space, 0.5)),
    );
    let (result, _) = solved.unwrap_or_else(|e| panic!("{e}"));
    let counts = counts.borrow().clone();
    assert_eq!(counts.len(), result.report.iterations);
    assert_eq!(counts, (1..=counts.len()).collect::<Vec<_>>());
    let stopped = observe_with(
        |u| {
            if u.count == 2 {
                ControlFlow::Break(())
            } else {
                ControlFlow::Continue(())
            }
        },
        || generate_cached(&space, &request(&space, 0.55)),
    );
    match stopped {
        Err(e) => assert_eq!(e.kind(), "Cancelled"),
        Ok(_) => panic!("served after a stop"),
    }
}

/// A 50 × 16 × 36 mm block with a bore along x on grid lines (y = 0,
/// z = 18): the rotated cylinder leaves the end faces' corners an ulp
/// off the grid planes through them, and the solver's clipping reads the
/// block short on the first grid (PK-963). The design (all of it, one
/// update: the fallback is under test, not the optimiser) serves on the
/// shifted grid and says so.
#[test]
fn a_design_space_that_reads_short_on_the_grid_is_solved_on_a_shifted_one() {
    use popkern_core::geometry::Transform3;
    let corners: Vec<Point3> = [0.0, 50.0]
        .iter()
        .flat_map(|&x| {
            [-8.0, 8.0]
                .into_iter()
                .flat_map(move |y| [0.0, 36.0].into_iter().map(move |z| Point3::new(x, y, z)))
        })
        .collect();
    let block = popkern_part::Part::polyhedron(&corners, tol()).unwrap_or_else(|e| panic!("{e}"));
    let bore = BRep::cylinder(3.0, 52.0, tol())
        .and_then(|c| {
            c.transform(
                &Transform3::translation(-1.0, 0.0, 18.0)
                    .compose(&Transform3::rotation_y(std::f64::consts::FRAC_PI_2)),
            )
        })
        .unwrap_or_else(|e| panic!("{e}"));
    let bore = popkern_part::Part::from_brep(bore, tol(), "bore").unwrap_or_else(|e| panic!("{e}"));
    let space = block
        .subtract(&bore)
        .unwrap_or_else(|e| panic!("{e}"))
        .brep()
        .clone();
    let wall = end_face(&space, -1.0, 0.0);
    let tip = end_face(&space, 1.0, 50.0);
    let request = GenerativeRequest {
        material: "PLA (FDM)".into(),
        supports: vec![SupportGroup {
            faces: vec![wall],
            kind: SupportKind::Fixed,
        }],
        loads: vec![LoadGroup {
            faces: vec![tip],
            force_n: 10.0,
            direction: [0.0, 0.0, -1.0],
        }],
        keep: vec![wall, tip],
        volume_fraction: 1.0,
        filter: 6.0,
        keep_depth: Some(2.0),
        cell: Some(2.0),
        build_direction: None,
        max_iterations: 1,
        robust: false,
        max_cells: None,
        mm_per_unit: 1.0,
    };
    let (result, _) = generate(&space, &request, &mut |_| ControlFlow::Continue(()))
        .unwrap_or_else(|e| panic!("{} {e}", e.kind()));
    result.brep.validate().unwrap_or_else(|e| panic!("{e}"));
    assert!(
        result
            .report
            .warnings
            .iter()
            .any(|w| w.starts_with("AnalysisGridShifted")),
        "{:?}",
        result.report.warnings
    );
    // The wall's main patch (the face that keeps the wall's name) stays
    // on x = 0 within the shift, a thousandth of a 2 mm cell; patches
    // beside the bore follow its cut cells' fill on either grid.
    let area = |f| {
        let p = result.brep.face_outer_loop_positions(f);
        let mut n = [0.0f64; 3];
        for i in 0..p.len() {
            let (a, b) = (p[i], p[(i + 1) % p.len()]);
            n[0] += (a.y - b.y) * (a.z + b.z);
            n[1] += (a.z - b.z) * (a.x + b.x);
            n[2] += (a.x - b.x) * (a.y + b.y);
        }
        (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt()
    };
    let main = result
        .brep
        .faces
        .keys()
        .zip(&result.sources)
        .filter(|(_, s)| **s == Some(wall))
        .map(|(f, _)| f)
        .max_by(|a, b| area(*a).total_cmp(&area(*b)))
        .unwrap_or_else(|| panic!("no wall patch"));
    let worst = result
        .brep
        .face_outer_loop_positions(main)
        .iter()
        .map(|p| p.x.abs())
        .fold(0.0f64, f64::max);
    assert!(worst < 0.005, "the wall moved {worst} mm");
}
