//! GeometryCache pcurve resolution: foreign basis surfaces, 2-D conics,
//! rational weights.

use super::*;
use popkern_core::geometry::Point3;
use popkern_core::geometry::nurbs_curve::NurbsCurve;

use crate::step::document::{OpaqueBag, StepDocument, StepProductNode, UnitContext};
use crate::step::export::StepExportOptions;
use crate::step::export::{ExportContext, ref_, str_param};
use crate::step::parser::entities::parse_step_file;
use crate::step::parser::graph::EntityGraph;
use crate::step::serializer::serialize_part21;

/// A face plane (#10), an identical copy (#11) and a copy shifted
/// by 1 along its normal (#12). Surface curve #40's pcurve lives on
/// the copy, #41's on the shifted plane (OpenCASCADE writes each
/// PCURVE on its own basis-surface entity).
const FOREIGN_BASIS_STEP: &str = "ISO-10303-21;
HEADER;
FILE_DESCRIPTION((''),'2;1');
FILE_NAME('','',(''),(''),'','','');
FILE_SCHEMA(('AUTOMOTIVE_DESIGN'));
ENDSEC;
DATA;
#1=CARTESIAN_POINT('',(0.,0.,0.));
#2=DIRECTION('',(0.,0.,1.));
#3=DIRECTION('',(1.,0.,0.));
#4=AXIS2_PLACEMENT_3D('',#1,#2,#3);
#5=CARTESIAN_POINT('',(0.,0.,1.));
#6=AXIS2_PLACEMENT_3D('',#5,#2,#3);
#10=PLANE('',#4);
#11=PLANE('',#4);
#12=PLANE('',#6);
#20=CARTESIAN_POINT('',(0.,0.));
#21=CARTESIAN_POINT('',(1.,0.));
#22=B_SPLINE_CURVE_WITH_KNOTS('',1,(#20,#21),.UNSPECIFIED.,.F.,.F.,(2,2),(0.,1.),.UNSPECIFIED.);
#23=(GEOMETRIC_REPRESENTATION_CONTEXT(2) PARAMETRIC_REPRESENTATION_CONTEXT() REPRESENTATION_CONTEXT('2D SPACE',''));
#24=DEFINITIONAL_REPRESENTATION('',(#22),#23);
#25=PCURVE('',#11,#24);
#26=PCURVE('',#12,#24);
#30=VECTOR('',#3,1.);
#31=LINE('',#1,#30);
#40=SURFACE_CURVE('',#31,(#25),.PCURVE_S1.);
#41=SURFACE_CURVE('',#31,(#26),.PCURVE_S1.);
ENDSEC;
END-ISO-10303-21;
";

#[test]
fn a_pcurve_on_an_agreeing_copy_of_the_face_surface_resolves() {
    let parsed = parse_step_file(FOREIGN_BASIS_STEP.as_bytes()).expect("parses");
    let mut log = DiagnosticLog::new();
    let graph = EntityGraph::build(parsed.entities, &mut log);
    let mut geometry = GeometryCache::new(&graph);

    let on_copy = geometry.resolve_surface_curve_pcurves(40, 10, 1e-6, &mut log);
    assert_eq!(on_copy.len(), 1, "the copy agrees with the face plane");
    assert_eq!(on_copy[0].pcurve_id, 25);
    let uv = on_copy[0].pcurve.point_at(on_copy[0].pcurve.t_end);
    assert!((uv.x - 1.0).abs() < 1e-12 && uv.y.abs() < 1e-12, "{uv:?}");
    assert!(log.is_empty(), "{:?}", log.entries());

    // Negative control: the shifted plane disagrees by 1 along the
    // pcurve, so no pcurve of #41 fits the face; declined, reported.
    let on_shifted = geometry.resolve_surface_curve_pcurves(41, 10, 1e-6, &mut log);
    assert!(on_shifted.is_empty());
    let warning = log
        .warnings()
        .find(|w| w.message.contains("no PCURVE of #41 fits face surface #10"))
        .unwrap_or_else(|| panic!("decline reported: {:?}", log.entries()));
    assert!(warning.message.contains("1.000e0"), "{}", warning.message);
}

#[test]
fn a_2d_circle_pcurve_resolves_exactly_on_its_2d_placement() {
    // OpenCASCADE's form: a CIRCLE on an AXIS2_PLACEMENT_2D with a
    // two-coordinate location and direction.
    let step = FOREIGN_BASIS_STEP.replace(
        "#22=B_SPLINE_CURVE_WITH_KNOTS('',1,(#20,#21),.UNSPECIFIED.,.F.,.F.,(2,2),(0.,1.),.UNSPECIFIED.);",
        "#50=CARTESIAN_POINT('',(17.,-27.));\n#51=DIRECTION('',(0.,-1.));\n\
         #52=AXIS2_PLACEMENT_2D('',#50,#51);\n#22=CIRCLE('',#52,3.);",
    );
    let parsed = parse_step_file(step.as_bytes()).expect("parses");
    let mut log = DiagnosticLog::new();
    let graph = EntityGraph::build(parsed.entities, &mut log);
    let mut geometry = GeometryCache::new(&graph);
    let resolved = geometry.resolve_surface_curve_pcurves(40, 10, 1e-6, &mut log);
    assert_eq!(resolved.len(), 1, "{:?}", log.entries());
    let pcurve = &resolved[0].pcurve;
    assert!(matches!(pcurve.geometry, PcurveGeometry::Elliptical { .. }));
    // Parameterized by angle from the placement's x direction (0,-1).
    for (t, (u, v)) in [
        (0.0, (17.0, -30.0)),
        (std::f64::consts::FRAC_PI_2, (20.0, -27.0)),
    ] {
        let uv = pcurve.point_at(t);
        assert!(
            (uv.x - u).abs() < 1e-12 && (uv.y - v).abs() < 1e-12,
            "t={t}: {uv:?}"
        );
    }
    assert!(log.is_empty(), "{:?}", log.entries());
}

#[test]
fn serialized_rational_pcurve_roundtrip_preserves_weights() {
    let mut export = ExportContext::new();
    let surface_id = export.emit("PLANE", vec![str_param(""), ref_(999)]);
    let pcurve_id = export.export_pcurve(
        &NurbsCurve {
            degree: 2,
            knots: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            control_points: vec![
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(0.5, 0.5, 0.0),
                Point3::new(1.0, 0.0, 0.0),
            ],
            weights: vec![1.0, 0.5, 1.0],
        },
        surface_id,
    );
    let bytes = serialize_part21(
        &export.entities,
        &StepDocument {
            products: Vec::<StepProductNode>::new(),
            roots: Vec::new(),
            units: UnitContext::default(),
            opaque_bag: OpaqueBag::default(),
        },
        &StepExportOptions::default(),
    )
    .expect("serialized pcurve STEP should write");
    let parsed = parse_step_file(&bytes).expect("serialized pcurve STEP should parse");
    let mut log = DiagnosticLog::new();
    let graph = EntityGraph::build(parsed.entities, &mut log);
    let mut geometry = GeometryCache::new(&graph);

    let resolved = geometry
        .resolve_pcurve_on_surface(pcurve_id, surface_id, 1e-6, &mut log)
        .unwrap_or_else(|gap| panic!("pcurve should resolve on its exported surface: {gap:?}"));

    assert_eq!(resolved.pcurve_id, pcurve_id);
    assert_ne!(resolved.reference_curve_id, resolved.pcurve_id);

    match resolved.pcurve.geometry {
        popkern_core::geometry::pcurve::PcurveGeometry::Nurbs(nurbs) => {
            assert_eq!(nurbs.weights, vec![1.0, 0.5, 1.0]);
        }
        other => panic!("expected NURBS pcurve, got {other:?}"),
    }
}
