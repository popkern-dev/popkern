//! The flush NESTED union class (PK-427): a body whose top face is
//! coplanar with, and strictly contained in, another body's bottom
//! face — a lid picking up its nesting lip. The containing face's
//! imprint ring closes only if the clip keeps the boundary-grazing
//! quarter-arc records on the rounded corners (the windowed flush
//! span), and a later cut through the joined region must respect the
//! ring face's HOLE (a hole-blind containment imprinted the slot
//! inside the hole — OCC `InvalidImbricationOfWires`, 8 open mesh
//! edges).

use popkern_script::prelude::*;

fn lid_with_lip(m: &mut Model) -> Body {
    let lid = m
        .profile(Plane::xy_at(0.72))
        .rect(1.6, 1.0)
        .sketch("lid_outline")
        .expect("sketch")
        .extrude(m, End::Blind(0.08.into()))
        .expect("extrude");
    let lid = lid
        .fillet(m, edges().straight().direction(Dir::Z), 0.15)
        .expect("lid corners");
    let lip = m
        .tool(|m| {
            let lip = m
                .profile(Plane::xy_at(0.66))
                .rect(1.44, 0.84)
                .sketch("lip_outline")?
                .extrude(m, End::Blind(0.06.into()))?;
            lip.fillet(m, edges().straight().direction(Dir::Z), 0.09)
        })
        .expect("lip tool");
    lid.union(m, lip).expect("flush nested union")
}

/// The union serves watertight with the lip ring imprinted as one
/// face carrying an inner loop.
#[test]
fn a_flush_nested_union_serves_watertight() {
    let mut m = Model::new();
    let lid = lid_with_lip(&mut m);
    m.assert_watertight().expect("watertight");
    let brep = m.tree().body_brep(lid.id).cloned().expect("brep");
    let ring_faces_with_hole = brep
        .faces
        .iter()
        .filter(|(_, f)| !f.inner_loops.is_empty())
        .count();
    assert!(
        ring_faces_with_hole >= 1,
        "the flush junction imprints an inner loop"
    );
}
