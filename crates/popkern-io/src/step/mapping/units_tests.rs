use super::*;
use crate::step::parser::entities::parse_step_file;

fn units_of(data: &str) -> FileUnits {
    let text =
        format!("ISO-10303-21;\nHEADER;\nENDSEC;\nDATA;\n{data}\nENDSEC;\nEND-ISO-10303-21;\n");
    let parsed = parse_step_file(text.as_bytes()).expect("parses");
    let mut log = DiagnosticLog::new();
    let graph = EntityGraph::build(parsed.entities, &mut log);
    resolve_file_units(&graph, &mut log)
}

const DEGREES_MM: &str = "#1204=(
GEOMETRIC_REPRESENTATION_CONTEXT(3)
GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT((#1205))
GLOBAL_UNIT_ASSIGNED_CONTEXT((#1206,#1207,#1212))
REPRESENTATION_CONTEXT('ID1','3D')
);
#1205=UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE(1.E-06),#1212,'closure','');
#1206=(NAMED_UNIT(*) SI_UNIT($,.STERADIAN.) SOLID_ANGLE_UNIT());
#1207=(CONVERSION_BASED_UNIT('DEGREES',#1208) NAMED_UNIT(#1209) PLANE_ANGLE_UNIT());
#1208=PLANE_ANGLE_MEASURE_WITH_UNIT(PLANE_ANGLE_MEASURE(0.01745329252),#1210);
#1209=DIMENSIONAL_EXPONENTS(0.,0.,0.,0.,0.,0.,0.);
#1210=(NAMED_UNIT(*) PLANE_ANGLE_UNIT() SI_UNIT($,.RADIAN.));
#1212=(LENGTH_UNIT() NAMED_UNIT(*) SI_UNIT(.MILLI.,.METRE.));";

#[test]
fn a_degree_file_reads_angles_in_degrees() {
    let u = units_of(DEGREES_MM);
    assert_eq!(u.mm_per_length, 1.0);
    assert!((u.radians_per_angle - 0.01745329252).abs() < 1e-15, "{u:?}");
}

#[test]
fn an_inch_file_reads_lengths_in_inches() {
    let text = DEGREES_MM.replace(
        "#1212=(LENGTH_UNIT() NAMED_UNIT(*) SI_UNIT(.MILLI.,.METRE.));",
        "#1212=(CONVERSION_BASED_UNIT('INCH',#1213) LENGTH_UNIT() NAMED_UNIT(#1209));\n\
         #1213=LENGTH_MEASURE_WITH_UNIT(LENGTH_MEASURE(25.4),#1214);\n\
         #1214=(LENGTH_UNIT() NAMED_UNIT(*) SI_UNIT(.MILLI.,.METRE.));",
    );
    let u = units_of(&text);
    assert!((u.mm_per_length - 25.4).abs() < 1e-12, "{u:?}");
}

#[test]
fn si_prefixes_scale_to_millimetres() {
    let text = DEGREES_MM.replace("SI_UNIT(.MILLI.,.METRE.)", "SI_UNIT(.CENTI.,.METRE.)");
    assert!((units_of(&text).mm_per_length - 10.0).abs() < 1e-12);
    let text = DEGREES_MM.replace("SI_UNIT(.MILLI.,.METRE.)", "SI_UNIT($,.METRE.)");
    assert!((units_of(&text).mm_per_length - 1000.0).abs() < 1e-12);
}

#[test]
fn no_unit_context_is_millimetres_and_radians() {
    assert_eq!(
        units_of("#1=CARTESIAN_POINT('',(0.,0.,0.));"),
        FileUnits::default()
    );
}
