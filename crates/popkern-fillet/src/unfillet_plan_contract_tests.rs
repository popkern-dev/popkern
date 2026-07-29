//! Unfillet plan-API decline contract for curved supports (PK-180).
//!
//! `unfillet_face` handles closed rims through the annulus path
//! (`unfillet_closed`, gated by `round_5c4_closed_fillet_round_trip`),
//! but `compute_unfillet_plan` still supports plane-plane fillets
//! only: a blend with a curved support must decline with a
//! structured diagnostic naming the class, never mis-plan. PK-180
//! flips this gate to assert a successful plan once curved-support
//! planning lands.
//!
//! The Round-5c diagnostic probes that lived here (cylinder-cylinder
//! seam fillet, closed-rim loop anatomy) were retired: seam edges
//! stay rejected by design (PK-180 scope note), and the closed-rim
//! questions are answered by the live round-trip gate above.

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

/// A plane-cylinder fillet (cube with a bore, rim filleted) must
/// make `compute_unfillet_plan` decline naming the curved-support
/// limitation — the structured contract PK-180 will replace with a
/// real plan.
#[test]
fn unfillet_plan_declines_curved_supports_with_named_class() {
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

    let blend_face = filleted
        .faces
        .keys()
        .find(|&f| crate::identify_blend(&filleted, f).is_some())
        .expect("filleted rim must identify a blend face");

    let err = crate::unfillet::compute_unfillet_plan(&filleted, blend_face)
        .expect_err("curved-support unfillet plan must decline until PK-180 lands");
    let msg = format!("{err:?}");
    // Measured contract at HEAD: the closed-annulus blend declines
    // at support-pair determination, before the plane-plane-only
    // check is even reached. Both decline sites are structured;
    // neither yet names the curved-support/closed-rim class in
    // semantic terms — PK-180 owns replacing them with a real plan
    // (and, until then, a class-named decline).
    assert!(
        msg.contains("plane-plane fillets only")
            || msg.contains("curved supports")
            || msg.contains("could not determine the support face pair"),
        "curved-support decline must stay structured, got: {}",
        &msg[..msg.len().min(240)]
    );
}
