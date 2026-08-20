//! PK-295 regression pins: the rim-filleted cylinder minus a
//! strip-crossing cuboid — the fixture whose graph-merge subtraction
//! once served a five-face open box with a clean `validate_all`
//! verdict (self-twinned half-edges blind every twin-based check).
//!
//! Promoted from the PK-295 diagnostic census. Pins, in one boolean:
//! - ZERO self-twinned half-edges in the served result (the walk
//!   wires both sides of every kept boundary edge);
//! - the operand's closed-fillet strip seam is an exact analytic Arc
//!   (the 17-point polyline's ~2.4e-4 sag put the seam's pierce node
//!   off every true intersection curve);
//! - the result carries the part's faces (cylinder, torus bulges,
//!   bottom cap), not just the cutter's.

use popkern_core::geometry::Transform3;
use popkern_core::geometry::curve::CurveGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_fillet::{FilletOptions, FilletRadius, FilletSpec, fillet_edges};

fn rim_filleted_cylinder(tol: Tolerance) -> BRep {
    let brep = BRep::cylinder(1.0, 2.0, tol).expect("cylinder");
    let mut rim = Vec::new();
    for (e, d) in brep.edges.iter() {
        let he = &brep.half_edges[d.half_edge];
        let a = brep.vertices[he.origin].position;
        let b = brep.vertices[brep.half_edges[he.twin].origin].position;
        let f1 = brep.loops[he.loop_].face;
        let f2 = brep.loops[brep.half_edges[he.twin].loop_].face;
        if f1 != f2 && (a.z - 2.0).abs() < 1e-9 && (b.z - 2.0).abs() < 1e-9 {
            rim.push(e);
        }
    }
    let specs: Vec<FilletSpec> = rim
        .iter()
        .map(|&edge| FilletSpec {
            edge,
            radius: FilletRadius::Constant(0.2),
        })
        .collect();
    fillet_edges(&brep, &specs, &FilletOptions::default(), tol, None)
        .expect("rim fillet must complete")
        .brep
}

fn strip_crossing_cutter(tol: Tolerance) -> BRep {
    BRep::unit_cube(tol)
        .expect("cube")
        .transform(&Transform3::scale_uniform(1.8))
        .expect("scale")
        .transform(&Transform3::translation(-0.9, -0.9, 0.6))
        .expect("translate")
}

#[test]
fn pk295_subtraction_has_no_self_twinned_half_edges() {
    let tol = Tolerance::default();
    let part = rim_filleted_cylinder(tol);

    // Operand pin: the strip's transverse seam is an exact Arc.
    let strip_seam_arcs = part
        .edges
        .iter()
        .filter(|(_, ed)| {
            let he = &part.half_edges[ed.half_edge];
            let f1 = part.loops[he.loop_].face;
            let f2 = part.loops[part.half_edges[he.twin].loop_].face;
            let on_torus = part
                .faces
                .get(f1)
                .and_then(|f| f.surface.as_ref())
                .is_some_and(|s| {
                    matches!(
                        s.geometry,
                        popkern_core::geometry::SurfaceGeometry::Torus { .. }
                    )
                });
            f1 == f2
                && on_torus
                && matches!(
                    ed.curve.as_ref().map(|c| &c.geometry),
                    Some(CurveGeometry::Arc { .. })
                )
        })
        .count();
    assert_eq!(
        strip_seam_arcs, 1,
        "the closed-fillet strip seam must be an exact Arc (PK-295)"
    );

    let cutter = strip_crossing_cutter(tol);
    let result = popkern_bool::subtraction(part, cutter, tol, None).expect("subtraction");

    let self_twins = result
        .half_edges
        .iter()
        .filter(|(heh, he)| he.twin == *heh)
        .count();
    assert_eq!(
        self_twins, 0,
        "served boolean must wire both sides of every kept boundary edge (PK-295)"
    );

    // The part's surface classes must all survive into the result —
    // the historical failure served ONLY the cutter's five planes.
    use popkern_core::geometry::SurfaceGeometry as S;
    let mut cylinders = 0;
    let mut tori = 0;
    for (_, fd) in result.faces.iter() {
        match fd.surface.as_ref().map(|s| &s.geometry) {
            Some(S::Cylinder { .. }) => cylinders += 1,
            Some(S::Torus { .. }) => tori += 1,
            _ => {}
        }
    }
    assert!(
        cylinders >= 1,
        "result must keep the cylinder wall, got {cylinders}"
    );
    assert_eq!(
        tori, 4,
        "result must keep all four torus bulges outside the cutter"
    );
    result.validate().expect("result validates");
}
