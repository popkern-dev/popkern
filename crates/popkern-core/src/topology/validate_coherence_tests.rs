use crate::geometry::Point3;
use crate::tolerance::Tolerance;
use crate::topology::BRep;
use crate::topology::coherence::{CoherenceValidationScope, PendingCoherenceMarker};

#[test]
fn unit_cube_has_no_l1_or_l3_violations() {
    let brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let report = brep.validate_coherence(&CoherenceValidationScope::AllSteadyState);
    assert!(
        report.l1_missing_edge_curves.is_empty(),
        "L1 on unit cube should be empty: {:?}",
        report.l1_missing_edge_curves
    );
    assert!(
        report.l3_missing_pcurves.is_empty(),
        "L3 on unit cube should be empty: {:?}",
        report.l3_missing_pcurves
    );
}

#[test]
fn empty_brep_has_no_violations() {
    let brep = BRep::new(Tolerance::default());
    let report = brep.validate_coherence(&CoherenceValidationScope::AllSteadyState);
    assert!(report.is_empty());
}

#[test]
fn mid_surgery_brep_reports_l1_for_curveless_edge() {
    // A mev edge has no curve and no surface attached — L1 fires.
    let mut brep = BRep::new(Tolerance::default());
    let mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
    let mev = brep
        .mev_at_vertex(mvfs.vertex, mvfs.loop_, Point3::new(1.0, 0.0, 0.0))
        .unwrap();

    let report = brep.validate_coherence(&CoherenceValidationScope::AllSteadyState);
    assert!(
        report.l1_missing_edge_curves.contains(&mev.edge),
        "L1 should report the curveless edge: {:?}",
        report.l1_missing_edge_curves
    );
}

#[test]
fn pending_marker_excuses_l1_violation() {
    // Record a pending marker for the curveless edge — L1 no
    // longer reports it as a violation.
    let mut brep = BRep::new(Tolerance::default());
    let mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
    let mev = brep
        .mev_at_vertex(mvfs.vertex, mvfs.loop_, Point3::new(1.0, 0.0, 0.0))
        .unwrap();
    brep.record_pending_authoring(PendingCoherenceMarker::EdgeCurvePending {
        edge: mev.edge,
        owner_operator: "test_op",
    });

    let report = brep.validate_coherence(&CoherenceValidationScope::AllSteadyState);
    assert!(
        !report.l1_missing_edge_curves.contains(&mev.edge),
        "L1 should NOT report the edge once a pending marker exists"
    );
}

#[test]
fn pending_only_scope_checks_marked_entities() {
    // PendingOnly walks the marker list, not the whole BRep.
    // Recording a marker for an edge that DOES lack a curve
    // surfaces the violation; an edge with a curve does not.
    let mut brep = BRep::new(Tolerance::default());
    let mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
    let mev = brep
        .mev_at_vertex(mvfs.vertex, mvfs.loop_, Point3::new(1.0, 0.0, 0.0))
        .unwrap();
    brep.record_pending_authoring(PendingCoherenceMarker::EdgeCurvePending {
        edge: mev.edge,
        owner_operator: "test_op",
    });

    // PendingOnly: the marker IS on a violating edge → report.
    let report = brep.validate_coherence(&CoherenceValidationScope::PendingOnly);
    assert!(report.l1_missing_edge_curves.contains(&mev.edge));
}

#[test]
fn pending_only_skips_resolved_edges() {
    // A marker on an edge that has since been given a curve does
    // not surface as a violation under PendingOnly.
    use crate::geometry::Vec3;
    use crate::geometry::curve::{CurveGeometry, CurveOnEdge};

    let mut brep = BRep::new(Tolerance::default());
    let mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
    let mev = brep
        .mev_at_vertex(mvfs.vertex, mvfs.loop_, Point3::new(1.0, 0.0, 0.0))
        .unwrap();
    brep.record_pending_authoring(PendingCoherenceMarker::EdgeCurvePending {
        edge: mev.edge,
        owner_operator: "test_op",
    });
    // Resolve the obligation: author the curve.
    brep.edges[mev.edge].curve = Some(CurveOnEdge {
        geometry: CurveGeometry::Line {
            origin: Point3::new(0.0, 0.0, 0.0),
            direction: Vec3::new(1.0, 0.0, 0.0),
        },
        t_start: 0.0,
        t_end: 1.0,
    });

    let report = brep.validate_coherence(&CoherenceValidationScope::PendingOnly);
    assert!(
        !report.l1_missing_edge_curves.contains(&mev.edge),
        "PendingOnly should NOT report an edge whose curve has been authored"
    );
}

#[test]
fn unit_cube_has_no_l2_violations() {
    // L2 holds on the cube: every cube face's planar surface
    // shares the same Z (or X / Y) plane equation with its
    // neighbour at the shared edge, so the surfaces agree on
    // the edge curve within tolerance.
    let brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let report = brep.validate_coherence(&CoherenceValidationScope::AllSteadyState);
    assert!(
        report.l2_face_boundary_disagreements.is_empty(),
        "L2 on unit cube should be empty: {:#?}",
        report.l2_face_boundary_disagreements
    );
}

#[test]
fn l2_fires_when_face_surface_is_displaced() {
    // Synthetic L2 violation: take the unit cube, displace one
    // face's plane origin by 0.1, and verify L2 surfaces at
    // least one disagreement.
    use crate::geometry::Vec3;
    use crate::geometry::surface::{SurfaceGeometry, SurfaceOnFace};

    let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let face = brep.faces.keys().next().unwrap();
    let surface = brep.faces[face].surface.as_ref().unwrap().clone();
    // Move the plane origin by 0.1 m along its normal — the
    // displaced surface no longer agrees with neighbouring faces
    // at the shared edges.
    let SurfaceGeometry::Plane {
        origin,
        u_axis,
        v_axis,
        normal,
    } = surface.geometry
    else {
        panic!("expected unit cube face surface to be Plane");
    };
    let shifted_origin = origin + Vec3::new(normal.x * 0.1, normal.y * 0.1, normal.z * 0.1);
    brep.faces[face].surface = Some(SurfaceOnFace {
        geometry: SurfaceGeometry::Plane {
            origin: shifted_origin,
            u_axis,
            v_axis,
            normal,
        },
        orientation: surface.orientation,
    });

    let report = brep.validate_coherence(&CoherenceValidationScope::AllSteadyState);
    assert!(
        !report.l2_face_boundary_disagreements.is_empty(),
        "L2 should report at least one disagreement after displacing a face surface by 0.1 m"
    );
    let worst_gap = report
        .l2_face_boundary_disagreements
        .iter()
        .map(|d| d.worst_gap)
        .fold(0.0_f64, f64::max);
    assert!(
        (worst_gap - 0.1).abs() < 1.0e-6,
        "worst L2 gap should be ≈ 0.1 m (the displacement); got {worst_gap}"
    );
}

#[test]
fn l3_only_fires_when_face_has_surface() {
    // A half-edge without a pcurve is NOT an L3 violation when
    // the owning face has no surface (no canonical UV space to
    // map into).
    let mut brep = BRep::new(Tolerance::default());
    let _mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
    // The mvfs face has no surface attached; the loop's half-edge
    // exists but has no pcurve — that's NOT an L3 violation.
    let report = brep.validate_coherence(&CoherenceValidationScope::AllSteadyState);
    assert!(
        report.l3_missing_pcurves.is_empty(),
        "L3 should be empty when faces have no surfaces: {:?}",
        report.l3_missing_pcurves
    );
}
