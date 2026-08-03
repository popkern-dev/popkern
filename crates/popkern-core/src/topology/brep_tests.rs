use super::*;
use crate::precision::EntityHandle;
use crate::tolerance::Tolerance;

#[test]
fn new_brep_is_empty() {
    let brep = BRep::new(Tolerance::default());
    assert_eq!(brep.vertices.len(), 0);
    assert_eq!(brep.half_edges.len(), 0);
    assert_eq!(brep.edges.len(), 0);
    assert_eq!(brep.loops.len(), 0);
    assert_eq!(brep.faces.len(), 0);
    assert_eq!(brep.shells.len(), 0);
    assert_eq!(brep.solids.len(), 0);
}

#[test]
fn brep_has_tolerance() {
    let tol = Tolerance::default().with_linear(1e-3);
    let brep = BRep::new(tol);
    assert_eq!(brep.tolerance.linear(), 1e-3);
}

#[test]
fn brep_serializes_round_trip() {
    let brep = BRep::new(Tolerance::default());
    let json = serde_json::to_string(&brep).unwrap();
    let brep2: BRep = serde_json::from_str(&json).unwrap();
    assert_eq!(brep2.vertices.len(), 0);
}

#[test]
fn brep_transform_translates_vertices() {
    use crate::geometry::transform::Transform3;
    let brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let t = Transform3::translation(10.0, 0.0, 0.0);
    let transformed = brep.transform(&t).unwrap();
    for (_, v) in &transformed.vertices {
        assert!(v.position.x >= 10.0 - 1e-10);
    }
}

#[test]
fn brep_transform_uniform_scale_keeps_unit_normals_and_scales_radii() {
    use crate::geometry::surface::SurfaceGeometry;
    use crate::geometry::transform::Transform3;
    let brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let t = Transform3::translation(0.6, 0.6, 0.6).compose(&Transform3::scale_uniform(0.6));
    let scaled = brep.transform(&t).unwrap();
    for (_, face) in &scaled.faces {
        let Some(sof) = face.surface.as_ref() else {
            continue;
        };
        if let SurfaceGeometry::Plane {
            normal,
            u_axis,
            v_axis,
            ..
        } = &sof.geometry
        {
            assert!(
                (normal.length() - 1.0).abs() < 1e-12,
                "scaled plane normal must stay unit, got |n| = {} (PK-231: the 0.6 \
                 notch normal displaced every fillet offset)",
                normal.length()
            );
            assert!((u_axis.length() - 1.0).abs() < 1e-12);
            assert!((v_axis.length() - 1.0).abs() < 1e-12);
        }
    }
    // Cylinder radius carries the scale; axis stays unit.
    let cyl = SurfaceGeometry::Cylinder {
        origin: crate::geometry::Point3::new(0.0, 0.0, 0.0),
        axis: crate::geometry::Vec3::new(0.0, 0.0, 1.0),
        radius: 0.08,
    };
    let scaled_cyl = super::brep_transforms::transform_surface(&cyl, &t, 0.6);
    match scaled_cyl {
        SurfaceGeometry::Cylinder { axis, radius, .. } => {
            assert!((axis.length() - 1.0).abs() < 1e-12);
            assert!((radius - 0.048).abs() < 1e-12);
        }
        other => panic!("cylinder must stay a cylinder, got {other:?}"),
    }
}

#[test]
fn brep_transform_declines_non_uniform_scale() {
    use crate::geometry::transform::Transform3;
    let brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let err = brep.transform(&Transform3::scale(1.0, 2.0, 1.0));
    assert!(
        matches!(err, Err(crate::error::CoreError::PreconditionFailed(_))),
        "non-uniform scale must decline with a typed error, got {err:?}"
    );
}

#[test]
fn similarity_scale_classifies_transforms() {
    use crate::geometry::transform::Transform3;
    let rigid = Transform3::rotation_z(0.7).compose(&Transform3::translation(1.0, 2.0, 3.0));
    assert!((rigid.similarity_scale().unwrap() - 1.0).abs() < 1e-12);
    let scaled = Transform3::scale_uniform(0.6).compose(&Transform3::rotation_x(0.3));
    assert!((scaled.similarity_scale().unwrap() - 0.6).abs() < 1e-12);
    assert!(
        Transform3::scale(1.0, 2.0, 1.0)
            .similarity_scale()
            .is_none()
    );
    let mirror = Transform3::reflection(
        crate::geometry::Point3::new(0.0, 0.0, 0.0),
        crate::geometry::Vec3::new(0.0, 0.0, 1.0),
    )
    .unwrap();
    assert!((mirror.similarity_scale().unwrap() - 1.0).abs() < 1e-12);
}

#[test]
fn face_outer_loop_positions_unit_cube() {
    let brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let face_handle = brep.faces.keys().next().unwrap();
    let positions = brep.face_outer_loop_positions(face_handle);
    assert_eq!(positions.len(), 4, "each cube face should have 4 vertices");
    for pos in &positions {
        assert!(
            pos.x >= -1e-10 && pos.x <= 1.0 + 1e-10,
            "x={} out of [0,1]",
            pos.x
        );
        assert!(
            pos.y >= -1e-10 && pos.y <= 1.0 + 1e-10,
            "y={} out of [0,1]",
            pos.y
        );
        assert!(
            pos.z >= -1e-10 && pos.z <= 1.0 + 1e-10,
            "z={} out of [0,1]",
            pos.z
        );
    }
}

#[test]
fn face_outer_loop_positions_invalid_face_returns_empty() {
    let brep = BRep::new(Tolerance::default());
    let bogus = FaceHandle::default();
    let positions = brep.face_outer_loop_positions(bogus);
    assert!(
        positions.is_empty(),
        "invalid handle should yield empty vec"
    );
}

#[test]
fn max_imprecision_default_from_tolerance() {
    let tol = Tolerance::default().with_linear(1e-4);
    let brep = BRep::new(tol);
    assert!((brep.max_imprecision - 1e-2).abs() < 1e-15);
}

#[test]
fn precision_summary_empty_brep() {
    let brep = BRep::new(Tolerance::default());
    let summary = brep.precision_summary();
    assert_eq!(summary.health, crate::precision::PrecisionHealth::Good);
    assert!((summary.health_score - 1.0).abs() < 1e-15);
    assert!(summary.worst_entity.is_none());
    assert_eq!(summary.degraded_count, 0);
    assert_eq!(summary.critical_count, 0);
}

#[test]
fn precision_summary_with_degraded_vertex() {
    let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let vh = brep.vertices.keys().next().unwrap();
    // Set tolerance to 60% of max_imprecision => degraded
    let val = brep.max_imprecision * 0.6;
    brep.set_entity_tolerance(EntityHandle::Vertex(vh), val);
    let summary = brep.precision_summary();
    assert_eq!(summary.health, crate::precision::PrecisionHealth::Degraded);
    assert_eq!(summary.degraded_count, 1);
    assert_eq!(summary.critical_count, 0);
    assert!(summary.worst_entity.is_some());
}

#[test]
fn precision_summary_with_critical_edge() {
    let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let eh = brep.edges.keys().next().unwrap();
    // Set tolerance to 95% of max_imprecision => critical
    let val = brep.max_imprecision * 0.95;
    brep.set_entity_tolerance(EntityHandle::Edge(eh), val);
    let summary = brep.precision_summary();
    assert_eq!(summary.health, crate::precision::PrecisionHealth::Critical);
    assert_eq!(summary.critical_count, 1);
}

#[test]
fn degraded_entities_sorted_worst_first() {
    let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let vhs: Vec<_> = brep.vertices.keys().collect();
    brep.set_entity_tolerance(EntityHandle::Vertex(vhs[0]), 1e-4);
    brep.set_entity_tolerance(EntityHandle::Vertex(vhs[1]), 5e-4);
    brep.set_entity_tolerance(EntityHandle::Vertex(vhs[2]), 2e-4);
    let entities = brep.degraded_entities();
    assert_eq!(entities.len(), 3);
    assert!(entities[0].1 >= entities[1].1);
    assert!(entities[1].1 >= entities[2].1);
}

#[test]
fn recover_entity_tolerance_tightens() {
    let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let vh = brep.vertices.keys().next().unwrap();
    brep.set_entity_tolerance(EntityHandle::Vertex(vh), 1e-3);
    assert!((brep.vertex_tolerance(vh) - 1e-3).abs() < 1e-15);
    // Recover to a tighter value
    brep.recover_entity_tolerance(EntityHandle::Vertex(vh), 1e-5);
    assert!((brep.vertex_tolerance(vh) - 1e-5).abs() < 1e-15);
}

#[test]
fn degradation_tracking_records() {
    use crate::precision::OperationPrecisionReport;
    let mut brep = BRep::new(Tolerance::default());
    brep.enable_degradation_tracking(10);
    let report = OperationPrecisionReport {
        entity_deltas: vec![],
        max_imprecision: 1e-6,
        exceeded_bound: vec![],
        diagnostics: vec![],
    };
    brep.record_operation("test_op", &report);
    // No panic, tracking is enabled
}

#[test]
fn set_entity_tolerance_face() {
    let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let fh = brep.faces.keys().next().unwrap();
    let val = brep.max_imprecision * 0.7;
    brep.set_entity_tolerance(EntityHandle::Face(fh), val);
    assert!((brep.face_tolerance(fh) - val).abs() < 1e-15);
}
