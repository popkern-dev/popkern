//! Shared builders for the history integration tests: the reference
//! document whose frozen v3 bytes live in `tests/fixtures/`.

#![allow(dead_code)]

use popkern_constraint::{Constraint, Sketch};
use popkern_core::naming::{EntityName, EntityRole, OperationId, RoleCategory};
use popkern_history::*;

/// A solved unit-square constraint sketch as JSON.
pub fn rectangle_sketch_json() -> String {
    let mut sk = Sketch::new();
    let p0 = sk.add_point(0.0, 0.0);
    let p1 = sk.add_point(1.0, 0.0);
    let p2 = sk.add_point(1.0, 1.0);
    let p3 = sk.add_point(0.0, 1.0);
    sk.add_line(p0, p1);
    sk.add_line(p1, p2);
    sk.add_line(p2, p3);
    sk.add_line(p3, p0);
    sk.add_constraint(Constraint::Fixed(p0));
    sk.add_constraint(Constraint::Fixed(p1));
    sk.add_constraint(Constraint::Fixed(p2));
    sk.add_constraint(Constraint::Fixed(p3));
    sk.solve().expect("rectangle should solve");
    serde_json::to_string(&sk).expect("rectangle should serialize")
}

/// The generated entity `index` of operation `op_id`, in `category`.
pub fn entity(op_id: u64, category: RoleCategory, index: u32) -> EntityName {
    EntityName {
        operation: OperationId(op_id),
        role: EntityRole { category, index },
        ancestors: Vec::new(),
    }
}

/// The reference document: a unit block (sketch op 1, extrude op 2)
/// whose height is an expression of a global `size` parameter, its
/// four vertical edges filleted (op 3) with a G2 recipe, a suppressed
/// chamfer (op 4), a `tall` configuration overriding `size`, and an
/// explicit tolerance-independent rollback bar at the end.
///
/// The frozen bytes in `tests/fixtures/reference-document-v3.pkn`
/// were written from exactly this builder by the version-3 writer
/// (commit 94f7a842) before the label table, naming-scheme field,
/// build tolerance, and `RawBody` landed. Keep the builder in step
/// with the fixture's content, never the other way round.
pub fn author_reference_document() -> FeatureTree {
    let mut tree = FeatureTree::new();
    let size = tree
        .add_global_parameter("size", 1.0, Unit::Meters)
        .expect("size");
    tree.add_feature(
        Feature::Sketch {
            plane: SketchPlaneDesc::XY,
            sketch_data: SketchPayload::Constraint(rectangle_sketch_json()),
        },
        "Base Sketch",
    );
    let height = tree
        .add_global_parameter("height", 1.0, Unit::Meters)
        .expect("height");
    tree.set_parameter_expression(
        height,
        Expression::BinOp(
            Box::new(Expression::Param(size)),
            BinOp::Mul,
            Box::new(Expression::Literal(0.75)),
        ),
    )
    .expect("height expression");
    tree.add_feature(
        Feature::Extrude {
            face: entity(1, RoleCategory::Generated, 0),
            direction: Some(Vec3Desc::Literal([0.0, 0.0, 1.0])),
            end_condition: ExtrudeEnd::Blind { distance: height },
            mode: OperationMode::NewBody,
            draft_angle: None,
            second_direction: None,
        },
        "Block",
    );
    let radius = tree
        .add_global_parameter("radius", 0.1, Unit::Meters)
        .expect("radius");
    tree.replay(popkern_core::tolerance::Tolerance::default())
        .expect("block replays");
    let edges = vertical_edge_names(tree.current_brep().expect("block body"));
    tree.add_feature(
        Feature::Fillet {
            edges,
            radius,
            recipe: FilletRecipeDesc {
                continuity: ContinuityKind::G2,
                ..Default::default()
            },
        },
        "Vertical Rounds",
    );
    let distance = tree
        .add_global_parameter("distance", 0.05, Unit::Meters)
        .expect("distance");
    let chamfer = tree.add_feature(
        Feature::Chamfer {
            edges: vec![entity(2, RoleCategory::Generated, 0)],
            distance,
            kind: ChamferKindDesc::Equal,
        },
        "Top Break",
    );
    tree.suppress(chamfer);
    let tall = tree.add_config("tall", "double size").expect("config");
    tree.set_config_override(tall, size, ParamOverride::Literal(2.0))
        .expect("override");
    tree
}

/// Stable names of the block's edges whose tangent is along ±Z.
pub fn vertical_edge_names(brep: &popkern_core::topology::BRep) -> Vec<EntityName> {
    let mut names = Vec::new();
    for (handle, edge) in &brep.edges {
        let tangent = popkern_measure::edge_tangent(brep, handle).expect("edge tangent");
        if tangent.z.abs() > 0.99 {
            let nid = edge.name.expect("block edges are named");
            names.push(brep.name_pool.get(nid).expect("name pool").clone());
        }
    }
    assert_eq!(names.len(), 4, "a block has four vertical edges");
    names
}
