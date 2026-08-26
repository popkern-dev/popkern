//! PK-328 exact-provenance naming gates: the moved skin's faces carry
//! their SOURCE FACE as ancestor (the centroid heuristic lost this for
//! every real thickness), the un-moved skin stays `Modified`, and the
//! interactive dual-closure path names identically to the plain path
//! (its `None` join sink was the recorded PK-328 gap).

use popkern_core::cancellation::CancellationToken;
use popkern_core::history::TrackingContext;
use popkern_core::naming::{EntityName, EntityRole, OperationId, RoleCategory};
use popkern_core::progress::OperationControl;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_shell::{ShellDirection, ShellOptions, shell, shell_with_control};

fn tol() -> Tolerance {
    Tolerance::default()
}

/// A unit cube with every face carrying a distinct input name.
fn named_cube() -> (BRep, Vec<EntityName>) {
    let mut brep = BRep::unit_cube(tol()).unwrap();
    let face_handles: Vec<_> = brep.faces.keys().collect();
    let mut names = Vec::new();
    for (i, fh) in face_handles.iter().enumerate() {
        let name = EntityName {
            operation: OperationId(7),
            role: EntityRole {
                category: RoleCategory::Generated,
                index: i as u32,
            },
            ancestors: Vec::new(),
        };
        let id = brep.name_pool.intern(name.clone(), None);
        brep.faces[*fh].name = Some(id);
        names.push(name);
    }
    (brep, names)
}

/// Count result faces per (category, has-a-source-face-ancestor) for a
/// named-cube shell result.
fn census(out: &BRep, face_names: &[EntityName]) -> (usize, usize, usize) {
    let mut modified_with_ancestor = 0;
    let mut generated_with_face_ancestor = 0;
    let mut generated_plain = 0;
    for (_, f) in &out.faces {
        let name = f
            .name
            .and_then(|id| out.name_pool.get(id))
            .expect("every face named under tracking");
        let face_ancestor = name.ancestors.len() == 1
            && out
                .name_pool
                .get(name.ancestors[0])
                .is_some_and(|a| face_names.contains(a));
        match (name.role.category, face_ancestor) {
            (RoleCategory::Modified, true) => modified_with_ancestor += 1,
            (RoleCategory::Generated, true) => generated_with_face_ancestor += 1,
            (RoleCategory::Generated, false) => generated_plain += 1,
            other => panic!("unexpected naming shape {other:?}"),
        }
    }
    (
        modified_with_ancestor,
        generated_with_face_ancestor,
        generated_plain,
    )
}

/// Classic inward shell: 6 un-moved faces `Modified` from their source
/// face, 6 moved faces `Generated` WITH their source face as ancestor.
/// The heuristic named the moved skin ancestor-less (its match radius
/// was 100·tol, the displacement is the wall thickness).
#[test]
fn moved_skin_faces_carry_source_face_ancestor() {
    let (brep, face_names) = named_cube();
    let mut ctx = TrackingContext::new(OperationId(99));
    let result = shell(
        &brep,
        0.1,
        ShellDirection::Inward,
        &[],
        &[],
        &ShellOptions::default(),
        tol(),
        Some(&mut ctx),
    )
    .unwrap();
    let (modified, offset, plain) = census(&result.brep, &face_names);
    assert_eq!(modified, 6, "un-moved outer skin is Modified-from-source");
    assert_eq!(offset, 6, "moved skin is Generated with source ancestor");
    assert_eq!(plain, 0, "no ancestor-less faces on a closed cube shell");
}

/// Midplane: BOTH skins moved — the outward skin continues each source
/// face (`Modified`), the inward skin is `Generated` from it. Nothing
/// falls through to ancestor-less names.
#[test]
fn midplane_names_both_skins_from_their_source_faces() {
    let (brep, face_names) = named_cube();
    let mut ctx = TrackingContext::new(OperationId(99));
    let result = shell(
        &brep,
        0.2,
        ShellDirection::Midplane,
        &[],
        &[],
        &ShellOptions::default(),
        tol(),
        Some(&mut ctx),
    )
    .unwrap();
    let (modified, offset, plain) = census(&result.brep, &face_names);
    assert_eq!(modified, 6, "outward midplane skin is Modified-from-source");
    assert_eq!(offset, 6, "inward midplane skin carries source ancestry");
    // The outward skin is a rounded box: 12 edge strips + 8 corner
    // spheres are join faces whose (unnamed) sources are edges and
    // vertices, so they stay plain Generated here.
    assert_eq!(plain, 20);
}

/// The interactive dual-closure path (with a live, never-fired token —
/// the detached worker actually runs) must name identically to the
/// plain path, join-face ancestors included. This closes the recorded
/// `None`-sink gap in `interactive.rs`.
#[test]
fn interactive_path_names_identically_to_plain_path() {
    // Outward shell so join strips exist; name faces AND edges.
    let (mut brep, _) = named_cube();
    let edge_handles: Vec<_> = brep.edges.keys().collect();
    for (i, eh) in edge_handles.iter().enumerate() {
        let name = EntityName {
            operation: OperationId(8),
            role: EntityRole {
                category: RoleCategory::Generated,
                index: 100 + i as u32,
            },
            ancestors: Vec::new(),
        };
        let id = brep.name_pool.intern(name.clone(), None);
        brep.edges[*eh].name = Some(id);
    }

    let name_census = |out: &BRep| -> Vec<(String, Vec<String>)> {
        let mut rows: Vec<_> = out
            .faces
            .iter()
            .map(|(_, f)| {
                let n = f
                    .name
                    .and_then(|id| out.name_pool.get(id))
                    .expect("face named");
                let ancestors: Vec<String> = n
                    .ancestors
                    .iter()
                    .filter_map(|a| out.name_pool.get(*a))
                    .map(|a| format!("{a:?}"))
                    .collect();
                (format!("{:?}/{}", n.role.category, n.role.index), ancestors)
            })
            .collect();
        rows.sort();
        rows
    };

    let mut ctx_plain = TrackingContext::new(OperationId(99));
    let plain = shell(
        &brep,
        0.1,
        ShellDirection::Outward,
        &[],
        &[],
        &ShellOptions::default(),
        tol(),
        Some(&mut ctx_plain),
    )
    .unwrap();

    let token = CancellationToken::new();
    let mut ctx_ctl = TrackingContext::new(OperationId(99));
    let controlled = shell_with_control(
        &brep,
        0.1,
        ShellDirection::Outward,
        &[],
        &[],
        &ShellOptions::default(),
        tol(),
        Some(&mut ctx_ctl),
        OperationControl::with_cancel(&token),
    )
    .unwrap();

    let plain_rows = name_census(&plain.brep);
    let controlled_rows = name_census(&controlled.brep);
    assert_eq!(plain_rows, controlled_rows, "identical name assignment");
    // The strips must actually carry their edge ancestors on BOTH paths
    // (an all-plain-Generated census would compare equal vacuously).
    let with_ancestors = plain_rows.iter().filter(|r| !r.1.is_empty()).count();
    assert!(
        with_ancestors >= 12,
        "expected at least the 12 strip ancestors, got {with_ancestors}"
    );
}
