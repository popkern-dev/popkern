use super::boundary::PreparedExactBoundaryCurve;
use super::corner::{
    GeneratedInnerSpoke, PreparedExactBoundaryInnerSpokes, PreparedExactBoundarySectorBoundary,
};
use super::pcurves::{
    ExactBoundaryPcurveBudget, ExactBoundaryPcurveDeclineReason, prepare_exact_boundary_pcurves,
};
use super::ribbon::{
    ExactBoundaryRibbonSample, PreparedExactBoundaryRibbon, PreparedExactBoundaryRibbons,
};
use super::sector::ExactBoundaryOuterSide;
use super::solve::{ExactBoundarySolveBudget, solve_exact_boundary_surfaces};
use super::support::ExactBoundarySupportSurfaceFamily;
use popkern_core::geometry::pcurve::PcurveGeometry;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{NurbsCurve, Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle};

#[test]
fn pcurves_follow_locked_boundaries_and_shared_parameters() {
    let spokes = fixture_spokes();
    let ribbons = fixture_ribbons(spokes.sectors[0].fillet_boundary.edge, false);
    let solved =
        solve_exact_boundary_surfaces(&spokes, &ribbons, ExactBoundarySolveBudget::default())
            .unwrap();

    let prepared = prepare_exact_boundary_pcurves(
        &spokes,
        &ribbons,
        &solved,
        ExactBoundaryPcurveBudget::default(),
    )
    .unwrap();

    assert_eq!(prepared.pairs.len(), 2);
    for pair in &prepared.pairs {
        assert!(pair.max_cap_residual <= 1.0e-6);
        assert!(pair.max_support_residual <= 1.0e-6);
        assert!(pair.max_pair_residual <= 1.0e-6);
        assert_eq!(pair.sector_name, "direct_apex.cap_sector.0");
        assert_eq!(pair.sector_index, 0);
        assert!(matches!(
            pair.support_pcurve.geometry,
            PcurveGeometry::Nurbs(_)
        ));
    }
}

#[test]
fn cap_residual_declines_with_stable_diagnostic() {
    let spokes = fixture_spokes();
    let ribbons = fixture_ribbons(spokes.sectors[0].fillet_boundary.edge, false);
    let mut solved =
        solve_exact_boundary_surfaces(&spokes, &ribbons, ExactBoundarySolveBudget::default())
            .unwrap();
    solved.surfaces[0].surface.control_points[4].z += 0.25;

    let diagnostic = prepare_exact_boundary_pcurves(
        &spokes,
        &ribbons,
        &solved,
        ExactBoundaryPcurveBudget::default(),
    )
    .expect_err("mutated cap boundary must decline");

    assert_eq!(
        diagnostic.reason,
        ExactBoundaryPcurveDeclineReason::CapResidualExceeded
    );
    assert_eq!(
        diagnostic.code,
        "exact_boundary_cap.pcurve.cap_residual_exceeded"
    );
}

#[test]
fn support_residual_declines_with_stable_diagnostic() {
    let spokes = fixture_spokes();
    let ribbons = fixture_ribbons(spokes.sectors[0].fillet_boundary.edge, true);
    let solved =
        solve_exact_boundary_surfaces(&spokes, &ribbons, ExactBoundarySolveBudget::default())
            .unwrap();

    let diagnostic = prepare_exact_boundary_pcurves(
        &spokes,
        &ribbons,
        &solved,
        ExactBoundaryPcurveBudget::default(),
    )
    .expect_err("stale support samples must decline");

    assert_eq!(
        diagnostic.reason,
        ExactBoundaryPcurveDeclineReason::SupportResidualExceeded
    );
    assert_eq!(
        diagnostic.code,
        "exact_boundary_cap.pcurve.support_residual_exceeded"
    );
}

fn fixture_spokes() -> PreparedExactBoundaryInnerSpokes {
    let edge = dummy_edge();
    let shared = Point3::new(0.0, 0.0, 0.0);
    let fillet_anchor = Point3::new(1.0, 0.0, 0.0);
    let face_anchor = Point3::new(0.0, 1.0, 0.0);
    let pole = Point3::new(1.0, 1.0, 0.2);
    PreparedExactBoundaryInnerSpokes {
        sectors: vec![PreparedExactBoundarySectorBoundary {
            sector_index: 0,
            sector_name: "direct_apex.cap_sector.0".to_string(),
            fillet_boundary: prepared_boundary(
                ExactBoundaryOuterSide::FilletEnd,
                edge,
                cubic_line(shared, fillet_anchor),
            ),
            face_boundary: prepared_boundary(
                ExactBoundaryOuterSide::FaceEdge,
                edge,
                cubic_line(shared, face_anchor),
            ),
            fillet_inner_spoke: generated_spoke(
                ExactBoundaryOuterSide::FilletEnd,
                fillet_anchor,
                pole,
            ),
            face_inner_spoke: generated_spoke(ExactBoundaryOuterSide::FaceEdge, face_anchor, pole),
            outer_corner_gap: 0.0,
        }],
        pole,
    }
}

fn fixture_ribbons(edge: EdgeHandle, stale_support: bool) -> PreparedExactBoundaryRibbons {
    let support_half_edge = dummy_half_edge();
    let support_face = dummy_face();
    PreparedExactBoundaryRibbons {
        ribbons: vec![
            ribbon(
                ExactBoundaryOuterSide::FilletEnd,
                edge,
                support_half_edge,
                support_face,
                stale_support,
            ),
            ribbon(
                ExactBoundaryOuterSide::FaceEdge,
                edge,
                support_half_edge,
                support_face,
                false,
            ),
        ],
    }
}

fn ribbon(
    side: ExactBoundaryOuterSide,
    edge: EdgeHandle,
    support_half_edge: HalfEdgeHandle,
    support_face: FaceHandle,
    stale_support: bool,
) -> PreparedExactBoundaryRibbon {
    PreparedExactBoundaryRibbon {
        sector_index: 0,
        sector_name: "direct_apex.cap_sector.0".to_string(),
        side,
        edge,
        support_half_edge,
        support_face,
        support_surface: plane_surface(),
        support_family: ExactBoundarySupportSurfaceFamily::Plane,
        support_parameter_reversed: false,
        support_half_edge_curve_forward: true,
        max_support_residual: if stale_support { 0.25 } else { 0.0 },
        samples: (0..4)
            .map(|i| {
                let parameter = i as f64 / 3.0;
                let boundary_point = match side {
                    ExactBoundaryOuterSide::FilletEnd => Point3::new(parameter, 0.0, 0.0),
                    ExactBoundaryOuterSide::FaceEdge => Point3::new(0.0, parameter, 0.0),
                };
                let support_point = if stale_support {
                    boundary_point + Vec3::new(0.0, 0.0, 0.25)
                } else {
                    boundary_point
                };
                ExactBoundaryRibbonSample {
                    parameter,
                    support_uv: (boundary_point.x, boundary_point.y),
                    boundary_point,
                    support_point,
                    boundary_tangent: Vec3::new(1.0, 0.0, 0.0),
                    support_normal: Vec3::new(0.0, 0.0, 1.0),
                    cross_derivative: Vec3::new(0.0, 0.2, 0.1),
                }
            })
            .collect(),
    }
}

fn plane_surface() -> SurfaceGeometry {
    SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    }
}

fn prepared_boundary(
    side: ExactBoundaryOuterSide,
    edge: EdgeHandle,
    curve: NurbsCurve,
) -> PreparedExactBoundaryCurve {
    PreparedExactBoundaryCurve {
        sector_index: 0,
        sector_name: "direct_apex.cap_sector.0".to_string(),
        side,
        edge,
        curve,
    }
}

fn generated_spoke(
    source_side: ExactBoundaryOuterSide,
    anchor: Point3,
    pole: Point3,
) -> GeneratedInnerSpoke {
    GeneratedInnerSpoke {
        name: "direct_apex.cap_sector.0.inner_spoke".to_string(),
        source_side,
        anchor,
        pole,
        curve: NurbsCurve {
            degree: 1,
            knots: vec![0.0, 0.0, 1.0, 1.0],
            control_points: vec![anchor, pole],
            weights: vec![1.0, 1.0],
        },
    }
}

fn cubic_line(start: Point3, end: Point3) -> NurbsCurve {
    NurbsCurve {
        degree: 3,
        knots: vec![0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0],
        control_points: (0..4)
            .map(|i| start + (end - start) * (i as f64 / 3.0))
            .collect(),
        weights: vec![1.0; 4],
    }
}

fn dummy_edge() -> EdgeHandle {
    BRep::unit_cube(Tolerance::default())
        .unwrap()
        .edges
        .keys()
        .next()
        .unwrap()
}

fn dummy_half_edge() -> HalfEdgeHandle {
    BRep::unit_cube(Tolerance::default())
        .unwrap()
        .half_edges
        .keys()
        .next()
        .unwrap()
}

fn dummy_face() -> FaceHandle {
    BRep::unit_cube(Tolerance::default())
        .unwrap()
        .faces
        .keys()
        .next()
        .unwrap()
}
