use super::*;
use crate::step::parser::types::{EntityInstance, Parameter};
use popkern_core::geometry::Vec3;
use popkern_core::topology::edge_segment;

fn make_entity(id: EntityId, type_name: &str, params: Vec<Parameter>) -> EntityInstance {
    EntityInstance {
        id,
        type_name: type_name.into(),
        params,
        complex_parts: vec![],
    }
}

fn p_str(s: &str) -> Parameter {
    Parameter::String(s.into())
}

fn p_ref(id: u64) -> Parameter {
    Parameter::EntityRef(id)
}

fn p_enum(s: &str) -> Parameter {
    Parameter::Enum(s.into())
}

fn p_real(v: f64) -> Parameter {
    Parameter::Real(v)
}

fn p_list(items: Vec<Parameter>) -> Parameter {
    Parameter::List(items)
}

#[test]
fn recovers_edge_curve_trim_from_topological_vertices() {
    let curve = CurveOnEdge {
        geometry: CurveGeometry::Circle {
            center: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            radius: 1.0,
        },
        t_start: 0.0,
        t_end: std::f64::consts::TAU,
    };
    let start = Point3::new(1.0, 0.0, 0.0);
    let end = Point3::new(0.0, 1.0, 0.0);
    let mut log = DiagnosticLog::new();

    let trimmed =
        recover_edge_curve_trim(99, curve, start, end, true, Tolerance::default(), &mut log);

    assert!(
        trimmed
            .geometry
            .point_at(trimmed.t_start)
            .distance_to(&start)
            <= Tolerance::default().linear()
    );
    assert!(
        trimmed.geometry.point_at(trimmed.t_end).distance_to(&end) <= Tolerance::default().linear()
    );
    assert!(
        (trimmed.t_end - trimmed.t_start).abs() < std::f64::consts::TAU,
        "trim recovery must not leave a full circle parameter range"
    );
    assert!(
        log.errors().next().is_none(),
        "trim recovery should not emit errors: {:?}",
        log.entries()
    );
}

#[test]
fn recovers_periodic_trim_across_zero_angle_with_same_sense() {
    let curve = CurveOnEdge {
        geometry: CurveGeometry::Circle {
            center: Point3::new(0.0, 0.0, 0.0),
            normal: Vec3::new(0.0, 0.0, 1.0),
            radius: 1.0,
        },
        t_start: 0.0,
        t_end: std::f64::consts::TAU,
    };
    let start_angle = 7.0 * std::f64::consts::TAU / 8.0;
    let end_angle = std::f64::consts::TAU / 8.0;
    let start = curve.geometry.point_at(start_angle);
    let end = curve.geometry.point_at(end_angle);
    let mut log = DiagnosticLog::new();

    let trimmed =
        recover_edge_curve_trim(99, curve, start, end, true, Tolerance::default(), &mut log);

    assert!(
        trimmed.t_end > trimmed.t_start,
        "same-sense periodic trims crossing zero must unwrap forward"
    );
    assert!(
        (trimmed.t_end - trimmed.t_start) < std::f64::consts::PI,
        "same-sense periodic trim should take the short forward span"
    );
    assert!(
        trimmed
            .geometry
            .point_at(trimmed.t_start)
            .distance_to(&start)
            <= Tolerance::default().linear()
    );
    assert!(
        trimmed.geometry.point_at(trimmed.t_end).distance_to(&end) <= Tolerance::default().linear()
    );
}

#[test]
fn keeps_bounded_nurbs_trim_when_domain_endpoints_match_topology() {
    let curve = CurveOnEdge {
        geometry: CurveGeometry::Nurbs(popkern_core::geometry::NurbsCurve {
            degree: 2,
            knots: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            control_points: vec![
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(0.5, 0.2, 0.0),
                Point3::new(1.0, 0.0, 0.0),
            ],
            weights: vec![1.0, 1.0, 1.0],
        }),
        t_start: 0.0,
        t_end: 1.0,
    };
    let start = curve.geometry.point_at(curve.t_start);
    let end = curve.geometry.point_at(curve.t_end);
    let mut log = DiagnosticLog::new();

    let trimmed =
        recover_edge_curve_trim(99, curve, start, end, true, Tolerance::default(), &mut log);

    assert_eq!(trimmed.t_start, 0.0);
    assert_eq!(trimmed.t_end, 1.0);
    assert!(
        log.entries().is_empty(),
        "trusted bounded NURBS domain should not emit trim recovery diagnostics: {:?}",
        log.entries()
    );
}

#[test]
fn reconciles_imported_edge_curve_from_valid_surface_trim_evidence() {
    let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let edge = brep
        .edges
        .keys()
        .find(|&edge| {
            brep.measure_edge_pcurve_deviation(edge)
                .is_some_and(|report| !report.half_edges.is_empty())
        })
        .expect("unit cube should have surfaced pcurve-backed edges");
    let (start, end) = edge_segment(&brep, edge).expect("edge has vertices");
    let offset = Vec3::new(0.0, 0.0, 0.05);
    brep.edges[edge].curve = Some(
        CurveOnEdge::line_segment(start + offset, end + offset, 1e-9)
            .expect("offset edge remains non-degenerate"),
    );
    let before = brep
        .measure_edge_pcurve_deviation(edge)
        .expect("corrupted edge remains measurable");
    assert!(before.max_deviation > 1e-3);
    let mut log = DiagnosticLog::new();

    reconcile_step_surface_curve_edges(&mut brep, Tolerance::default(), &mut log);

    let after = brep
        .measure_edge_pcurve_deviation(edge)
        .expect("reconciled edge remains measurable");
    assert!(after.max_deviation <= Tolerance::default().linear());
    assert!(
        log.entries()
            .iter()
            .any(|entry| entry.message.contains("normalized 1 STEP SURFACE_CURVE")),
        "accepted import normalization must be observable through structured diagnostics: {:?}",
        log.entries()
    );
}

#[test]
fn import_fidelity_diagnostics_report_step_entity_provenance() {
    let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
    let edge = brep
        .edges
        .keys()
        .find(|&edge| {
            brep.measure_edge_pcurve_deviation(edge)
                .is_some_and(|report| !report.half_edges.is_empty())
        })
        .expect("unit cube should have pcurve-backed edges");
    let initial = brep
        .measure_edge_pcurve_deviation(edge)
        .expect("selected edge should be measurable");
    let (start, end) = edge_segment(&brep, edge).expect("edge has vertices");
    brep.edges[edge].curve = Some(
        CurveOnEdge::line_segment(start + Vec3::new(0.0, 0.0, 0.05), end, 1e-9)
            .expect("corrupted edge remains non-degenerate"),
    );

    let mut provenance = StepTopologyImportProvenance::default();
    provenance.edge_curves.insert(edge, 300);
    for half_edge in initial.half_edges {
        provenance.half_edges.insert(
            half_edge.half_edge,
            StepHalfEdgeImportProvenance {
                edge_curve_id: 300,
                oriented_edge_id: 400,
                face_id: 500,
                surface_id: 600,
                pcurve_id: Some(700),
                reference_curve_id: Some(800),
                selected_pcurve_gap: Some(0.0),
            },
        );
    }
    let mut log = DiagnosticLog::new();

    report_step_import_fidelity_diagnostics(&brep, &provenance, Tolerance::default(), &mut log);

    assert!(
        log.entries().iter().any(|entry| {
            entry.entity_id == Some(300)
                && entry
                    .message
                    .contains("STEP import pcurve fidelity blocker")
                && entry.message.contains("edge_curve=#300")
                && entry.message.contains("oriented_edge=#400")
                && entry.message.contains("pcurve=#700")
                && entry.message.contains("reference_curve=#800")
        }),
        "pcurve fidelity diagnostics must carry STEP topology and pcurve provenance: {:?}",
        log.entries()
    );
}

/// Build a minimal triangle solid from STEP entities and verify topology.
#[test]
fn build_triangle_solid() {
    let entities = vec![
        make_entity(
            1,
            "CARTESIAN_POINT",
            vec![
                p_str(""),
                p_list(vec![p_real(0.0), p_real(0.0), p_real(0.0)]),
            ],
        ),
        make_entity(
            2,
            "CARTESIAN_POINT",
            vec![
                p_str(""),
                p_list(vec![p_real(1.0), p_real(0.0), p_real(0.0)]),
            ],
        ),
        make_entity(
            3,
            "CARTESIAN_POINT",
            vec![
                p_str(""),
                p_list(vec![p_real(0.0), p_real(1.0), p_real(0.0)]),
            ],
        ),
        make_entity(
            4,
            "DIRECTION",
            vec![
                p_str(""),
                p_list(vec![p_real(0.0), p_real(0.0), p_real(1.0)]),
            ],
        ),
        make_entity(
            5,
            "DIRECTION",
            vec![
                p_str(""),
                p_list(vec![p_real(1.0), p_real(0.0), p_real(0.0)]),
            ],
        ),
        make_entity(
            6,
            "AXIS2_PLACEMENT_3D",
            vec![p_str(""), p_ref(1), p_ref(4), p_ref(5)],
        ),
        make_entity(7, "PLANE", vec![p_str(""), p_ref(6)]),
        make_entity(
            8,
            "DIRECTION",
            vec![
                p_str(""),
                p_list(vec![p_real(1.0), p_real(0.0), p_real(0.0)]),
            ],
        ),
        make_entity(9, "VECTOR", vec![p_str(""), p_ref(8), p_real(1.0)]),
        make_entity(
            10,
            "DIRECTION",
            vec![
                p_str(""),
                p_list(vec![p_real(-1.0), p_real(1.0), p_real(0.0)]),
            ],
        ),
        make_entity(11, "VECTOR", vec![p_str(""), p_ref(10), p_real(1.0)]),
        make_entity(
            12,
            "DIRECTION",
            vec![
                p_str(""),
                p_list(vec![p_real(0.0), p_real(-1.0), p_real(0.0)]),
            ],
        ),
        make_entity(13, "VECTOR", vec![p_str(""), p_ref(12), p_real(1.0)]),
        make_entity(14, "LINE", vec![p_str(""), p_ref(1), p_ref(9)]),
        make_entity(15, "LINE", vec![p_str(""), p_ref(2), p_ref(11)]),
        make_entity(16, "LINE", vec![p_str(""), p_ref(3), p_ref(13)]),
        make_entity(20, "VERTEX_POINT", vec![p_str(""), p_ref(1)]),
        make_entity(21, "VERTEX_POINT", vec![p_str(""), p_ref(2)]),
        make_entity(22, "VERTEX_POINT", vec![p_str(""), p_ref(3)]),
        make_entity(
            30,
            "EDGE_CURVE",
            vec![p_str(""), p_ref(20), p_ref(21), p_ref(14), p_enum("T")],
        ),
        make_entity(
            31,
            "EDGE_CURVE",
            vec![p_str(""), p_ref(21), p_ref(22), p_ref(15), p_enum("T")],
        ),
        make_entity(
            32,
            "EDGE_CURVE",
            vec![p_str(""), p_ref(22), p_ref(20), p_ref(16), p_enum("T")],
        ),
        make_entity(
            40,
            "ORIENTED_EDGE",
            vec![
                p_str(""),
                Parameter::Derived,
                Parameter::Derived,
                p_ref(30),
                p_enum("T"),
            ],
        ),
        make_entity(
            41,
            "ORIENTED_EDGE",
            vec![
                p_str(""),
                Parameter::Derived,
                Parameter::Derived,
                p_ref(31),
                p_enum("T"),
            ],
        ),
        make_entity(
            42,
            "ORIENTED_EDGE",
            vec![
                p_str(""),
                Parameter::Derived,
                Parameter::Derived,
                p_ref(32),
                p_enum("T"),
            ],
        ),
        make_entity(
            50,
            "EDGE_LOOP",
            vec![p_str(""), p_list(vec![p_ref(40), p_ref(41), p_ref(42)])],
        ),
        make_entity(
            60,
            "FACE_OUTER_BOUND",
            vec![p_str(""), p_ref(50), p_enum("T")],
        ),
        make_entity(
            70,
            "ADVANCED_FACE",
            vec![p_str(""), p_list(vec![p_ref(60)]), p_ref(7), p_enum("T")],
        ),
        make_entity(80, "CLOSED_SHELL", vec![p_str(""), p_list(vec![p_ref(70)])]),
        make_entity(90, "MANIFOLD_SOLID_BREP", vec![p_str("body"), p_ref(80)]),
    ];

    let mut log = DiagnosticLog::new();
    let graph = EntityGraph::build(entities, &mut log);
    assert!(!log.has_errors(), "graph build errors: {:?}", log.entries());

    let mut geo_cache = GeometryCache::new(&graph);
    let tolerance = Tolerance::default();

    let solids = map_brep_solids(&graph, &mut geo_cache, &mut log, tolerance);
    assert_eq!(solids.len(), 1, "expected 1 solid, got {}", solids.len());

    let (entity_id, brep) = &solids[0];
    assert_eq!(*entity_id, 90);
    assert_eq!(brep.vertices.len(), 3);
    assert_eq!(brep.edges.len(), 3);
    assert_eq!(brep.faces.len(), 1);
    assert_eq!(brep.solids.len(), 1);
    assert_eq!(brep.shells.len(), 1);
    assert_eq!(
        brep.boundary_edges.len(),
        3,
        "single-face CLOSED_SHELL is invalid closed incidence and must import as open boundary topology"
    );
    assert!(
        log.warnings().any(|entry| entry
            .message
            .contains("downgrading imported shell to open boundary topology")),
        "closed-shell incidence downgrade should be diagnosed: {:?}",
        log.entries()
    );
}

/// A CLOSED rim (start vertex == end vertex) projects both trim
/// parameters to the same value; the periodic normalization must
/// expand the trim to the full period rather than leaving a
/// zero-length range (which collapsed every imported seam/rim
/// circle to a point).
#[test]
fn closed_rim_trim_recovery_spans_full_period() {
    use popkern_core::geometry::{Point3, Vec3};
    let mut start_t = 1.25;
    let mut end_t = 1.25;
    let circle = popkern_core::geometry::curve::CurveGeometry::Circle {
        center: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        radius: 0.43,
    };
    super::normalize_periodic_parameter_order(&circle, true, &mut start_t, &mut end_t);
    assert!(
        (end_t - start_t - std::f64::consts::TAU).abs() < 1e-9,
        "same_sense closed rim must cover a full period, got [{start_t}, {end_t}]"
    );
    let (mut s2, mut e2) = (2.0, 2.0);
    super::normalize_periodic_parameter_order(&circle, false, &mut s2, &mut e2);
    assert!(
        (s2 - e2 - std::f64::consts::TAU).abs() < 1e-9,
        "reversed closed rim must cover a full period, got [{s2}, {e2}]"
    );
}

/// Vertex projection can land turns away (Newton from spread starts):
/// the trim is the vertex angles folded into one turn, so an open arc
/// never spans more than a period (turbine.stp rims came out −15.6 →
/// 3.2, three turns, and their faces' pcurves livelocked the mesher).
#[test]
fn periodic_trims_fold_to_one_turn() {
    use popkern_core::geometry::{Point3, Vec3};
    let circle = popkern_core::geometry::curve::CurveGeometry::Circle {
        center: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(1.0, 0.0, 0.0),
        radius: 105.45,
    };
    let tau = std::f64::consts::TAU;
    for same_sense in [true, false] {
        let (mut s, mut e) = (-15.6483345620379, 3.239043580358187);
        super::normalize_periodic_parameter_order(&circle, same_sense, &mut s, &mut e);
        let span = (e - s).abs();
        assert!(
            span > 0.0 && span < tau,
            "same_sense {same_sense}: [{s}, {e}]"
        );
        assert!((0.0..2.0 * tau).contains(&s) && (0.0..2.0 * tau).contains(&e));
        // Same vertices as before, modulo the period.
        let same =
            |a: f64, b: f64| ((a - b).rem_euclid(tau)).min(tau - (a - b).rem_euclid(tau)) < 1e-9;
        assert!(same(s, -15.6483345620379) && same(e, 3.239043580358187));
    }
}
