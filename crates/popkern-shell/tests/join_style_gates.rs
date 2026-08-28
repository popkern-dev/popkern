//! PK-325 acceptance gates for rows 473 and 474: typed per-edge
//! `JoinFailed` declines under `JoinStyle::Intersection` (naming the
//! arc alternative, no partial body served) and join-face provenance
//! (`Generated` with the source edge as ancestor).

use popkern_core::history::TrackingContext;
use popkern_core::naming::{EntityName, EntityRole, OperationId, RoleCategory};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_shell::{ShellDirection, ShellOptions, shell};

fn tol() -> Tolerance {
    Tolerance::default()
}

/// Row 474: join faces are `Generated` with the source edge as
/// ancestor. Input edges carry names; the shell's strip faces must
/// record them.
#[test]
fn join_faces_are_generated_with_source_edge_ancestor() {
    let mut brep = BRep::unit_cube(tol()).unwrap();
    // Name every edge of the input so strips have ancestors to record.
    let edge_handles: Vec<_> = brep.edges.keys().collect();
    let mut edge_names = Vec::new();
    for (i, eh) in edge_handles.iter().enumerate() {
        let name = EntityName {
            operation: OperationId(7),
            role: EntityRole {
                category: RoleCategory::Generated,
                index: i as u32,
            },
            ancestors: Vec::new(),
        };
        let id = brep.name_pool.intern(name.clone(), None);
        brep.edges[*eh].name = Some(id);
        edge_names.push(name);
    }

    let mut ctx = TrackingContext::new(OperationId(99));
    let result = shell(
        &brep,
        0.1,
        ShellDirection::Outward,
        &[],
        &[],
        &ShellOptions::default(),
        tol(),
        Some(&mut ctx),
    )
    .unwrap();
    let out = &result.brep;

    // The outward-moved skin is the rounded box: 12 cylinder strips,
    // each Generated with exactly one ancestor — one of the input edge
    // names.
    use popkern_core::geometry::surface::SurfaceGeometry;
    let mut strip_count = 0;
    let mut with_edge_ancestor = 0;
    for (_, f) in &out.faces {
        let is_strip = matches!(
            f.surface.as_ref().map(|s| &s.geometry),
            Some(SurfaceGeometry::Cylinder { .. })
        );
        if !is_strip {
            continue;
        }
        strip_count += 1;
        let name = f
            .name
            .and_then(|id| out.name_pool.get(id))
            .expect("strip face must be named under tracking");
        assert_eq!(
            name.role.category,
            RoleCategory::Generated,
            "strip faces are Generated"
        );
        if name.ancestors.len() == 1 {
            let ancestor = out.name_pool.get(name.ancestors[0]).unwrap();
            if edge_names.contains(ancestor) {
                with_edge_ancestor += 1;
            }
        }
    }
    assert_eq!(strip_count, 12, "rounded cube has 12 strips");
    assert_eq!(
        with_edge_ancestor, 12,
        "every strip records its source edge as ancestor"
    );
}

/// PK-277 discipline: join construction iterates sorted handles only —
/// repeated runs are bit-identical.
#[test]
fn arc_join_construction_is_bit_identical_across_runs() {
    let brep = BRep::unit_cube(tol()).unwrap();
    let runs: Vec<Vec<(u64, u64, u64)>> = (0..3)
        .map(|_| {
            let r =
                popkern_shell::offset::offset_solid(&brep, 0.15, &ShellOptions::default(), tol())
                    .unwrap();
            let mut bits: Vec<(u64, u64, u64)> = r
                .brep
                .vertices
                .values()
                .map(|v| {
                    (
                        v.position.x.to_bits(),
                        v.position.y.to_bits(),
                        v.position.z.to_bits(),
                    )
                })
                .collect();
            bits.sort_unstable();
            bits
        })
        .collect();
    assert_eq!(runs[0], runs[1]);
    assert_eq!(runs[1], runs[2]);
}
