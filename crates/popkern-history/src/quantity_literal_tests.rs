//! Unit-suffixed literals carry their dimension (PK-754): `L - 180mm`
//! types as a length and evaluates in SI; a bare number summed with a
//! length is an error that names the fix; formatting writes the literal
//! back in the unit it was written in, so format -> parse is a fixed
//! point; documents serialize the new variant and keep reading old ones.

use crate::error::HistoryError;
use crate::expression::{Expression, dimension_kind, evaluate};
use crate::format::format_expression;
use crate::measure_query::MeasureContext;
use crate::parameter::ParameterTable;
use crate::parse::parse_expression;
use crate::unit::{DimensionKind, Unit};

fn table() -> ParameterTable {
    let mut t = ParameterTable::new();
    t.add_global("L", 0.425, Unit::Meters).expect("L");
    t.add_global("a", 0.5, Unit::Radians).expect("a");
    t
}

fn kind_and_value(src: &str) -> (DimensionKind, f64) {
    let t = table();
    let e = parse_expression(src, &t).expect("parse");
    let k = dimension_kind(&e, &t).expect("types");
    let v = evaluate(&e, &t, &MeasureContext::empty()).expect("evaluates");
    (k, v)
}

#[test]
fn suffixed_literals_type_and_evaluate_with_their_dimension() {
    let (k, v) = kind_and_value("L - 180mm");
    assert_eq!(k, DimensionKind::Length);
    assert!((v - 0.245).abs() < 1e-12);
    let (k, v) = kind_and_value("2*L + 5mm");
    assert_eq!(k, DimensionKind::Length);
    assert!((v - 0.855).abs() < 1e-12);
    let (k, v) = kind_and_value("a + 10deg");
    assert_eq!(k, DimensionKind::Angle);
    assert!((v - (0.5 + 10f64.to_radians())).abs() < 1e-12);
    // A bare factor still scales.
    let (k, _) = kind_and_value("L * 0.5");
    assert_eq!(k, DimensionKind::Length);
}

#[test]
fn a_bare_number_summed_with_a_length_names_the_suffix() {
    let t = table();
    let e = parse_expression("L - 180", &t).expect("parse");
    let err = dimension_kind(&e, &t).expect_err("bare literal in a sum");
    assert!(
        matches!(
            err,
            HistoryError::UnitlessLiteralInSum {
                dimension: DimensionKind::Length
            }
        ),
        "{err:?}"
    );
    assert!(err.to_string().contains("unit suffix"), "{err}");
}

#[test]
fn formatting_keeps_the_written_unit_and_round_trips() {
    let t = table();
    for src in [
        "L - 180mm",
        "2 * L + 5mm",
        "a + 10deg",
        "L - 2.5in",
        "L + 3cm",
    ] {
        let e = parse_expression(src, &t).expect("parse");
        let text = format_expression(&e, &t);
        let again = parse_expression(&text, &t).expect("reparse");
        assert_eq!(format_expression(&again, &t), text, "{src} -> {text}");
    }
    let e = parse_expression("L - 180mm", &t).expect("parse");
    assert!(format_expression(&e, &t).contains("180mm"));
}

#[test]
fn documents_serialize_quantities_and_read_old_literals() {
    let q = Expression::Quantity {
        value: 0.18,
        unit: Unit::Millimeters,
    };
    let json = serde_json::to_string(&q).expect("serialize");
    let back: Expression = serde_json::from_str(&json).expect("deserialize");
    assert!(
        matches!(back, Expression::Quantity { value, unit: Unit::Millimeters } if value == 0.18)
    );
    // A document written before PK-754 carries plain literals.
    let old: Expression = serde_json::from_str(r#"{"Literal":0.18}"#).expect("old literal");
    assert!(matches!(old, Expression::Literal(v) if v == 0.18));
}
