//! Unit harness for the support-cell solve's pinch handling: a spur
//! (out-and-back slit over doubled records) built with public Euler
//! ops on a real cube, exercising the fold-pair split without the
//! 15-second cube_extreme fixture.

use popkern_core::geometry::Point3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};

use super::resolve_pocket_faces;

/// Cube top face pinched with a two-edge spur lying in the x=0 plane.
/// The main loop fits the z=1 carrier; the spur fits only x=0 — the
/// solve must split the pinched slit off as its own cell without
/// Steiner splits (no new vertices) and leave every touched face's
/// records consistent.
fn cube_with_slit_spur() -> (BRep, FaceHandle) {
    let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
    // Find the z=1 top face and its corner at (0,0,1).
    let top = brep
        .faces
        .iter()
        .find(|(_, face)| {
            face.surface.as_ref().is_some_and(|s| {
                use popkern_core::geometry::SurfaceGeometry;
                matches!(
                    &s.geometry,
                    SurfaceGeometry::Plane { origin, normal, .. }
                        if normal.z.abs() > 0.9 && (origin.z - 1.0).abs() < 1e-9
                )
            })
        })
        .map(|(handle, _)| handle)
        .expect("unit cube has a z=1 face");
    // Find a top-loop half-edge whose origin is the corner (0,0,1);
    // `mev` splices the spur before that half-edge.
    let top_loop = brep.faces[top].outer_loop;
    let start = brep.loops[top_loop].half_edge;
    let mut corner_he = None;
    let mut he = start;
    for _ in 0..brep.half_edges.len() {
        let data = &brep.half_edges[he];
        let p = brep.vertices[data.origin].position;
        if p.distance_to(&Point3::new(0.0, 0.0, 1.0)) < 1e-9 {
            corner_he = Some(he);
            break;
        }
        let next = data.next;
        he = next;
        if he == start {
            break;
        }
    }
    let corner_he = corner_he.expect("top loop passes the (0,0,1) corner");
    // Spur into the x=0 plane: corner -> S1 -> S2 (out-and-back).
    let s1 = brep
        .mev(corner_he, Point3::new(0.0, 0.4, 0.7))
        .expect("first spur mev");
    let _s2 = brep
        .mev(s1.half_edge_twin, Point3::new(0.0, 0.7, 0.4))
        .expect("second spur mev");
    // The face under solve is surface-less pocket scaffolding.
    brep.faces[top].surface = None;
    (brep, top)
}

#[test]
fn pinched_slit_spur_splits_without_steiner_vertices() {
    let (mut brep, _top) = cube_with_slit_spur();
    let vertices_before = brep.vertices.len();
    let report = resolve_pocket_faces(&mut brep, &[]);
    assert!(
        report.declined.is_empty(),
        "pinched spur must solve, not decline: {:?}",
        report.declined
    );
    assert!(
        report.cells_created >= 1,
        "the x=0 spur must split off as its own cell (cells={})",
        report.cells_created
    );
    assert_eq!(
        brep.vertices.len(),
        vertices_before,
        "pinch resolution must not mint Steiner vertices"
    );
    // Every face this solve surfaced must contain its own boundary.
    for (handle, face) in brep.faces.iter() {
        let Some(surface) = face.surface.as_ref() else {
            continue;
        };
        use popkern_core::geometry::surface::SurfaceEval;
        let Some(loop_) = brep.loops.get(face.outer_loop) else {
            continue;
        };
        let start = loop_.half_edge;
        let mut he = start;
        for _ in 0..brep.half_edges.len() {
            let Some(data) = brep.half_edges.get(he) else {
                break;
            };
            let p = brep.vertices[data.origin].position;
            let d =
                popkern_core::geometry::uv_projection::project_point_to_uv(&surface.geometry, &p)
                    .map(|(u, v)| surface.geometry.point_at(u, v).distance_to(&p))
                    .unwrap_or(f64::INFINITY);
            assert!(
                d < 1e-6,
                "face {handle:?} boundary vertex {p:?} is {d:.3e} off its surface"
            );
            he = data.next;
            if he == start {
                break;
            }
        }
    }
    assert!(
        brep.validate().is_ok(),
        "solved model must validate: {:?}",
        brep.validate()
    );
}

/// The carved spur cell must carry the x=0 carrier and hold the spur
/// vertices on its loop. (The parent keeps the authored closing
/// appendage and stays in the tolerated bare-wrap class — recovering
/// it is the fold-pair operator's remaining work, see the PK-074
/// ledger.)
#[test]
fn pinched_slit_spur_cell_carries_its_carrier() {
    let (mut brep, top) = cube_with_slit_spur();
    let report = resolve_pocket_faces(&mut brep, &[]);
    assert!(report.declined.is_empty(), "{:?}", report.declined);
    let spur_cell_exists = brep.faces.iter().any(|(handle, face)| {
        handle != top
            && face.surface.as_ref().is_some_and(|s| {
                use popkern_core::geometry::SurfaceGeometry;
                matches!(
                    &s.geometry,
                    SurfaceGeometry::Plane { origin, normal, .. }
                        if normal.x.abs() > 0.9 && origin.x.abs() < 1e-9
                )
            })
            && {
                let Some(loop_) = brep.loops.get(face.outer_loop) else {
                    return false;
                };
                let start = loop_.half_edge;
                let mut he = start;
                let mut touches_spur = false;
                for _ in 0..brep.half_edges.len() {
                    let Some(data) = brep.half_edges.get(he) else {
                        break;
                    };
                    let p = brep.vertices[data.origin].position;
                    if p.distance_to(&Point3::new(0.0, 0.4, 0.7)) < 1e-9
                        || p.distance_to(&Point3::new(0.0, 0.7, 0.4)) < 1e-9
                    {
                        touches_spur = true;
                    }
                    he = data.next;
                    if he == start {
                        break;
                    }
                }
                touches_spur
            }
    });
    assert!(
        spur_cell_exists,
        "a split-off cell must carry the x=0 carrier and the spur vertices"
    );
}
