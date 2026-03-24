//! End-to-end integration tests for the Part 21 parser pipeline.

use popkern_io::step::diagnostics::DiagnosticLog;
use popkern_io::step::parser::entities::parse_step_file;
use popkern_io::step::parser::graph::EntityGraph;
use popkern_io::step::parser::types::Parameter;

const MINIMAL_STEP: &[u8] = b"ISO-10303-21;
HEADER;
FILE_DESCRIPTION(('A minimal STEP file'), '2;1');
FILE_NAME('minimal.step', '2026-03-24', ('popkern'), (''), 'popkern-io', 'popkern', '');
FILE_SCHEMA(('AUTOMOTIVE_DESIGN'));
ENDSEC;
DATA;
#1 = CARTESIAN_POINT('origin', (0.0, 0.0, 0.0));
#2 = CARTESIAN_POINT('x', (1.0, 0.0, 0.0));
#3 = DIRECTION('z', (0.0, 0.0, 1.0));
#4 = DIRECTION('x', (1.0, 0.0, 0.0));
#5 = AXIS2_PLACEMENT_3D('placement', #1, #3, #4);
#6 = PLANE('', #5);
ENDSEC;
END-ISO-10303-21;
";

#[test]
fn parse_minimal_step_file() {
    let result = parse_step_file(MINIMAL_STEP).unwrap();
    assert_eq!(result.entities.len(), 6);
    assert_eq!(result.header.schemas, vec!["AUTOMOTIVE_DESIGN"]);
    assert!(result.parse_errors.is_empty());
    assert_eq!(result.header.name, "minimal.step");
    assert_eq!(result.header.implementation_level, "2;1");
}

#[test]
fn entity_graph_resolves_refs() {
    let result = parse_step_file(MINIMAL_STEP).unwrap();
    let mut log = DiagnosticLog::new();
    let graph = EntityGraph::build(result.entities, &mut log);
    assert!(
        !log.has_errors(),
        "unexpected errors: {:?}",
        log.errors().collect::<Vec<_>>()
    );
    assert_eq!(graph.len(), 6);

    // Verify AXIS2_PLACEMENT_3D (#5) references #1, #3, #4
    let axis = graph.get(5).unwrap();
    assert_eq!(axis.type_name, "AXIS2_PLACEMENT_3D");
    assert_eq!(axis.params[1], Parameter::EntityRef(1));
    assert_eq!(axis.params[2], Parameter::EntityRef(3));
    assert_eq!(axis.params[3], Parameter::EntityRef(4));

    // Verify PLANE (#6) references #5
    let plane = graph.get(6).unwrap();
    assert_eq!(plane.type_name, "PLANE");
    assert_eq!(plane.params[1], Parameter::EntityRef(5));
}

#[test]
fn parse_with_broken_ref_is_lenient() {
    let input = b"ISO-10303-21;
HEADER;
FILE_DESCRIPTION(('broken ref test'), '2;1');
FILE_NAME('broken.step', '', (''), (''), '', '', '');
FILE_SCHEMA(('AUTOMOTIVE_DESIGN'));
ENDSEC;
DATA;
#1 = CARTESIAN_POINT('origin', (0.0, 0.0, 0.0));
#2 = AXIS2_PLACEMENT_3D('', #1, #999, $);
ENDSEC;
END-ISO-10303-21;
";
    let result = parse_step_file(input).unwrap();
    assert!(result.parse_errors.is_empty());
    assert_eq!(result.entities.len(), 2);

    let mut log = DiagnosticLog::new();
    let graph = EntityGraph::build(result.entities, &mut log);
    assert!(log.has_errors());
    assert_eq!(log.errors().count(), 1);
    // Both entities are still present
    assert_eq!(graph.len(), 2);
    assert!(graph.get(1).is_some());
    assert!(graph.get(2).is_some());
}

#[test]
fn parse_empty_data_section() {
    let input = b"ISO-10303-21;
HEADER;
FILE_DESCRIPTION(('empty'), '2;1');
FILE_NAME('empty.step', '', (''), (''), '', '', '');
FILE_SCHEMA(('AUTOMOTIVE_DESIGN'));
ENDSEC;
DATA;
ENDSEC;
END-ISO-10303-21;
";
    let result = parse_step_file(input).unwrap();
    assert_eq!(result.entities.len(), 0);
    assert_eq!(result.header.schemas, vec!["AUTOMOTIVE_DESIGN"]);
    assert!(result.parse_errors.is_empty());
}
