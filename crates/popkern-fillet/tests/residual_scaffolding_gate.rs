//! PK-153 gate: blend results must not ship unglued construction
//! scaffolding.
//!
//! After `fillet_edges` returns, every surviving
//! `ConstructionResidual` face must (a) own at least three sides —
//! two-sided holding faces and ghost records (faces whose sides were
//! re-owned by production loops) are retired by
//! `residual_scaffolding::retire_unglued_construction_residual_scaffolding`
//! — and (b) be glued to the production shell (at least one side
//! twinning a production-owned side). The retirement must never roll
//! back on these fixtures.

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, TopologyRole};
use popkern_fillet::{FilletOptions, FilletRadius, FilletSpec, fillet_edges};

fn assert_no_unglued_scaffolding(brep: &BRep, label: &str) {
    let mut owned_sides: std::collections::BTreeMap<_, usize> = std::collections::BTreeMap::new();
    for (_, half_edge) in &brep.half_edges {
        if let Some(face) = brep
            .loops
            .get(half_edge.loop_)
            .map(|loop_data| loop_data.face)
        {
            *owned_sides.entry(face).or_default() += 1;
        }
    }
    for (face, data) in &brep.faces {
        if data.role != TopologyRole::ConstructionResidual {
            continue;
        }
        let owned = owned_sides.get(&face).copied().unwrap_or(0);
        assert!(
            owned >= 3,
            "[{label}] residual face {face:?} owns only {owned} side(s) — ghost or degenerate \
             holding face survived scaffolding retirement"
        );
        let glued = brep.half_edges.iter().any(|(_, half_edge)| {
            brep.loops.get(half_edge.loop_).map(|l| l.face) == Some(face)
                && brep
                    .half_edges
                    .get(half_edge.twin)
                    .and_then(|twin| brep.loops.get(twin.loop_))
                    .and_then(|loop_data| brep.faces.get(loop_data.face))
                    .is_some_and(|owner| owner.role != TopologyRole::ConstructionResidual)
        });
        assert!(
            glued,
            "[{label}] residual face {face:?} has no production-twinned side — unglued \
             scaffolding survived retirement"
        );
    }
}

#[test]
fn mixed_radius_cube_ships_no_unglued_scaffolding() {
    for use_rational_gregory in [false, true] {
        let cube = BRep::unit_cube(Tolerance::default()).unwrap();
        let vertex = cube.vertices.keys().next().unwrap();
        let incident_edges: Vec<_> = cube
            .half_edges
            .values()
            .filter(|he| he.origin == vertex)
            .map(|he| he.edge)
            .collect();
        let radii = [0.08, 0.10, 0.12];
        let specs: Vec<FilletSpec> = incident_edges
            .iter()
            .zip(radii.iter())
            .map(|(&e, &r)| FilletSpec {
                edge: e,
                radius: FilletRadius::Constant(r),
            })
            .collect();
        let opts = FilletOptions {
            use_rational_gregory,
            ..FilletOptions::default()
        };
        let result = fillet_edges(&cube, &specs, &opts, Tolerance::default(), None)
            .expect("mixed-radius cube fillet must build");
        assert!(
            !result.repair_report.residual_scaffolding.rolled_back,
            "scaffolding retirement must not roll back on the mixed cube \
             (use_rational_gregory={use_rational_gregory})"
        );
        assert_no_unglued_scaffolding(&result.brep, &format!("gregory={use_rational_gregory}"));
    }
}
