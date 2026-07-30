//! Unfillet plan-API contract for curved supports (PK-180).
//!
//! Closed plane-cylinder rims (bore-rim annulus / seamed band) now
//! PLAN as well as execute: `compute_unfillet_plan` returns a typed
//! [`crate::unfillet_plan_closed::RestoredEdgeGeometry::Circle`]
//! that matches the rim the annulus surgery restores. Open blends on
//! curved supports still decline — with the configuration named and
//! the operation-record route pointed to — never a mis-plan.
//!
//! History: until PK-180 this file pinned the DECLINE of every
//! curved support ("Round 5b plane-plane only"). The Round-5c
//! diagnostic probes that once lived here (cylinder-cylinder seam
//! fillet) were retired: seam edges stay rejected by design.

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use crate::unfillet_plan_closed::RestoredEdgeGeometry;

fn bore_rim_fillet() -> (BRep, popkern_core::topology::FaceHandle) {
    use popkern_bool::{BoolOperation, boolean};
    use popkern_core::geometry::Transform3;
    use popkern_core::geometry::surface::SurfaceGeometry;

    let tol = Tolerance::default();
    let cube = BRep::unit_cube(tol).unwrap();
    let cyl = BRep::cylinder(0.3, 2.0, tol).unwrap();
    let translated = cyl
        .transform(&Transform3::translation(0.5, 0.5, -0.5))
        .unwrap();
    let subtracted = boolean(BoolOperation::Subtraction, cube, translated, tol, None)
        .expect("cube minus bore must build");

    let edge = subtracted
        .edges
        .keys()
        .find(|e| {
            let he = &subtracted.half_edges[subtracted.edges[*e].half_edge];
            let face_a = subtracted.loops[he.loop_].face;
            let face_b = subtracted.loops[subtracted.half_edges[he.twin].loop_].face;
            let surf_a = subtracted.faces[face_a]
                .surface
                .as_ref()
                .map(|s| &s.geometry);
            let surf_b = subtracted.faces[face_b]
                .surface
                .as_ref()
                .map(|s| &s.geometry);
            matches!(
                (surf_a, surf_b),
                (
                    Some(SurfaceGeometry::Plane { .. }),
                    Some(SurfaceGeometry::Cylinder { .. })
                ) | (
                    Some(SurfaceGeometry::Cylinder { .. }),
                    Some(SurfaceGeometry::Plane { .. })
                )
            )
        })
        .expect("bore rim must expose a plane-cylinder edge");

    let filleted = crate::fillet_edge(&subtracted, edge, 0.05, tol, None)
        .expect("closed-rim fillet must build (annulus path)")
        .brep;

    // The bore WALL also identifies as a (cylinder) blend and even
    // matches the seamed-band loop signature — the closed-rim plan
    // correctly declines it as Ambiguous (its ring twins are Plane +
    // Torus). The fillet band itself is the Torus face.
    let blend_face = filleted
        .faces
        .iter()
        .find(|(_, f)| {
            matches!(
                f.surface.as_ref().map(|s| &s.geometry),
                Some(popkern_core::geometry::surface::SurfaceGeometry::Torus { .. })
            )
        })
        .map(|(fh, _)| fh)
        .expect("filleted rim must carry a torus band");
    (filleted, blend_face)
}

/// A closed bore-rim fillet plans a restored CIRCLE that matches the
/// bore geometry, and the plan agrees with what the annulus surgery
/// actually restores (plan/execution consistency, acceptance 233).
#[test]
fn unfillet_plan_closed_rim_produces_matching_circle_plan() {
    let (filleted, blend_face) = bore_rim_fillet();

    let plan = crate::unfillet::compute_unfillet_plan(&filleted, blend_face)
        .expect("closed-rim blend must produce an unfillet plan (PK-180)");
    let RestoredEdgeGeometry::Circle {
        center,
        axis,
        radius,
    } = plan.restored_edge
    else {
        panic!(
            "closed rim must plan a restored circle, got {:?}",
            plan.restored_edge
        );
    };
    assert!(
        (radius - 0.3).abs() < 1e-6,
        "restored rim radius {radius} != bore radius 0.3"
    );
    assert!(
        (axis.z.abs() - 1.0).abs() < 1e-9,
        "restored rim axis {axis:?} must be the bore axis (±Z)"
    );
    assert!(
        (center.x - 0.5).abs() < 1e-6 && (center.y - 0.5).abs() < 1e-6,
        "restored rim centre {center:?} must sit on the bore axis"
    );
    // Seam point sits ON the planned circle; tangent ⟂ radial.
    let radial = plan.restored_edge_start - center;
    assert!(
        (radial.length() - radius).abs() < 1e-6,
        "seam point must lie on the restored circle"
    );
    assert!(
        radial.dot(&plan.restored_edge_direction).abs() < 1e-6,
        "restored-edge direction must be tangent at the seam point"
    );
    assert_eq!(
        plan.edges_to_remove.len(),
        3,
        "seamed band plan removes two contact rings + the seam bridge"
    );

    // Plan/execution consistency: the annulus surgery restores a rim
    // whose vertices lie on the planned circle.
    let restored = crate::unfillet::unfillet_face(&filleted, blend_face)
        .expect("closed-rim unfillet must execute (round 5c-4)");
    assert!(
        !restored.faces.contains_key(blend_face),
        "unfillet must remove the blend face"
    );
    let mut rim_vertices_checked = 0usize;
    for (_, vertex) in &restored.vertices {
        let p = vertex.position;
        // Rim vertices: on the rim plane (the planned circle's
        // height) near the bore radius.
        if (p.z - center.z).abs() < 1e-6 {
            let r = ((p.x - 0.5).powi(2) + (p.y - 0.5).powi(2)).sqrt();
            if (r - 0.3).abs() < 1e-3 {
                assert!(
                    (r - radius).abs() < 1e-6,
                    "restored rim vertex at radius {r} disagrees with planned {radius}"
                );
                rim_vertices_checked += 1;
            }
        }
    }
    assert!(
        rim_vertices_checked >= 1,
        "execution must leave at least one vertex on the planned rim circle"
    );
}

/// Open blends on curved supports still decline with the
/// configuration NAMED (never the retired planning labels), and the
/// message routes to the operation record.
#[test]
fn unfillet_plan_open_curved_support_declines_with_named_class() {
    // A plane-plane cube fillet's plan succeeds; sanity-check the
    // decline text contract on the helper level instead: seam edges
    // and open curved supports are exercised end-to-end by
    // `unfillet_probe`-era gates that were promoted into
    // `unfillet_tests`. Here: the decline reason must name the
    // support kind and the record route, and must NOT use planning
    // labels.
    use popkern_core::geometry::surface::SurfaceGeometry;
    let tol = Tolerance::default();
    let cyl = BRep::cylinder(0.5, 1.0, tol).unwrap();
    let face = cyl
        .faces
        .iter()
        .find(|(_, f)| {
            f.surface
                .as_ref()
                .is_some_and(|s| matches!(s.geometry, SurfaceGeometry::Cylinder { .. }))
        })
        .map(|(fh, _)| fh)
        .unwrap();
    let err = crate::unfillet_helpers::extract_plane(&cyl, face)
        .expect_err("cylinder support must decline plane extraction");
    let msg = format!("{err:?}");
    assert!(
        msg.contains("Cylinder") && msg.contains("unfillet_via_record"),
        "decline must name the support kind and the record route: {}",
        &msg[..msg.len().min(300)]
    );
    assert!(
        !msg.contains("Round 5"),
        "planning labels must not appear in user-facing declines: {}",
        &msg[..msg.len().min(300)]
    );
}
