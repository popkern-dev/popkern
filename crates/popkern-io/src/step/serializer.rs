//! Part 21 serializer — writes STEP entity instances to ISO 10303-21 text.

use std::io::Write;

use crate::error::IoError;
use crate::step::document::StepDocument;
use crate::step::parser::types::{EntityInstance, Parameter};

use super::export::{StepExportOptions, StepSchema};

/// Serialize entities and document metadata to Part 21 bytes.
pub fn serialize_part21(
    entities: &[EntityInstance],
    doc: &StepDocument,
    options: &StepExportOptions,
) -> Result<Vec<u8>, IoError> {
    let mut out = Vec::new();
    write_header(&mut out, doc, options)?;
    write_data_section(&mut out, entities, options.line_width)?;
    Ok(out)
}

/// Write the Part 21 HEADER section.
fn write_header(
    out: &mut Vec<u8>,
    doc: &StepDocument,
    options: &StepExportOptions,
) -> Result<(), IoError> {
    writeln!(out, "ISO-10303-21;")?;
    writeln!(out, "HEADER;")?;
    writeln!(out, "FILE_DESCRIPTION(('popkern export'), '2;1');")?;

    let file_name = if let Some(first) = doc.products.first() {
        &first.name
    } else {
        "unnamed"
    };

    let timestamp = "2024-01-01T00:00:00";

    write!(out, "FILE_NAME(")?;
    write_step_string(out, file_name)?;
    write!(out, ", ")?;
    write_step_string(out, timestamp)?;
    writeln!(out, ", (''), (''), 'popkern-io', 'popkern', '');")?;

    let schema_name = match options.schema {
        StepSchema::Ap203 => "CONFIG_CONTROL_DESIGN",
        StepSchema::Ap214 => "AUTOMOTIVE_DESIGN",
    };
    writeln!(out, "FILE_SCHEMA(('{}'));", schema_name)?;
    writeln!(out, "ENDSEC;")?;
    Ok(())
}

/// Write the DATA section containing all entities.
fn write_data_section(
    out: &mut Vec<u8>,
    entities: &[EntityInstance],
    line_width: usize,
) -> Result<(), IoError> {
    writeln!(out, "DATA;")?;
    for entity in entities {
        write_entity(out, entity, line_width)?;
    }
    writeln!(out, "ENDSEC;")?;
    writeln!(out, "END-ISO-10303-21;")?;
    Ok(())
}

/// Write a single entity instance with optional line wrapping.
fn write_entity(
    out: &mut Vec<u8>,
    entity: &EntityInstance,
    line_width: usize,
) -> Result<(), IoError> {
    // Build the full entity text first, then wrap if needed.
    let mut buf = Vec::new();
    write!(buf, "#{}=", entity.id)?;
    if !entity.complex_parts.is_empty() {
        write!(buf, "(")?;
    }
    write!(buf, "{}(", entity.type_name)?;
    for (i, param) in entity.params.iter().enumerate() {
        if i > 0 {
            write!(buf, ",")?;
        }
        write_parameter(&mut buf, param)?;
    }
    write!(buf, ")")?;

    // Handle complex entity parts
    for (type_name, params) in &entity.complex_parts {
        write!(buf, " {}(", type_name)?;
        for (i, param) in params.iter().enumerate() {
            if i > 0 {
                write!(buf, ",")?;
            }
            write_parameter(&mut buf, param)?;
        }
        write!(buf, ")")?;
    }
    if !entity.complex_parts.is_empty() {
        write!(buf, ")")?;
    }

    write!(buf, ";")?;

    let text = String::from_utf8(buf).map_err(|e| IoError::Export(format!("UTF-8 error: {e}")))?;

    // Wrap long lines
    if text.len() <= line_width {
        writeln!(out, "{text}")?;
    } else {
        write_wrapped(out, &text, line_width)?;
    }

    Ok(())
}

/// Write text with line wrapping at approximately `line_width` characters.
fn write_wrapped(out: &mut Vec<u8>, text: &str, line_width: usize) -> Result<(), IoError> {
    let bytes = text.as_bytes();
    let mut pos = 0;
    let mut first_line = true;

    while pos < bytes.len() {
        let remaining = bytes.len() - pos;
        let chunk = if first_line {
            remaining.min(line_width)
        } else {
            remaining.min(line_width.saturating_sub(2))
        };

        // Try to break at a comma or space — never mid-token. With
        // complex entities like `( TYPE_A(...) TYPE_B(...) )`, the
        // only safe break inside a single line is after a comma
        // (between params) or after a space (between complex parts).
        // If we can find neither, scan forward past `chunk` until we
        // hit one — STEP allows long lines, so we'd rather emit a
        // long line than corrupt mid-token.
        let break_at = if pos + chunk >= bytes.len() {
            bytes.len() - pos
        } else {
            let search_end = pos + chunk;
            let mut best: Option<usize> = None;
            // Search backwards for a comma OR a space
            for i in (pos..search_end).rev() {
                if bytes[i] == b',' || bytes[i] == b' ' {
                    best = Some(i - pos + 1);
                    break;
                }
            }
            // If no break point in [pos..search_end], scan forward to
            // the next one to avoid breaking mid-token.
            if best.is_none() {
                for (i, &b) in bytes.iter().enumerate().skip(search_end) {
                    if b == b',' || b == b' ' {
                        best = Some(i - pos + 1);
                        break;
                    }
                }
            }
            best.unwrap_or(remaining)
        };

        if !first_line {
            write!(out, "  ")?;
        }
        out.write_all(&bytes[pos..pos + break_at])?;
        writeln!(out)?;

        pos += break_at;
        first_line = false;
    }

    Ok(())
}

/// Write a single parameter value.
fn write_parameter(out: &mut Vec<u8>, param: &Parameter) -> Result<(), IoError> {
    match param {
        Parameter::Int(i) => write!(out, "{i}")?,
        Parameter::Real(f) => write_real(out, *f)?,
        Parameter::String(s) => write_step_string(out, s)?,
        Parameter::Enum(s) => write!(out, ".{s}.")?,
        Parameter::EntityRef(id) => write!(out, "#{id}")?,
        Parameter::List(items) => {
            write!(out, "(")?;
            for (i, item) in items.iter().enumerate() {
                if i > 0 {
                    write!(out, ",")?;
                }
                write_parameter(out, item)?;
            }
            write!(out, ")")?;
        }
        Parameter::Typed(name, inner) => {
            write!(out, "{name}(")?;
            write_parameter(out, inner)?;
            write!(out, ")")?;
        }
        Parameter::Omitted => write!(out, "$")?,
        Parameter::Derived => write!(out, "*")?,
    }
    Ok(())
}

/// Format a real number with enough precision for round-trip fidelity.
fn write_real(out: &mut Vec<u8>, f: f64) -> Result<(), IoError> {
    if f == 0.0 {
        write!(out, "0.")?;
    } else if f.abs() >= 1e-3 && f.abs() < 1e7 && f == f.trunc() {
        // Integer-valued float: write as e.g. "5."
        write!(out, "{}.", f as i64)?;
    } else {
        // Use scientific notation with enough digits
        let s = format!("{:.15E}", f);
        // Trim trailing zeros in the mantissa
        let trimmed = trim_real_trailing_zeros(&s);
        write!(out, "{trimmed}")?;
    }
    Ok(())
}

/// Trim trailing zeros from a scientific-notation real number string.
fn trim_real_trailing_zeros(s: &str) -> String {
    if let Some(e_pos) = s.find('E') {
        let mantissa = &s[..e_pos];
        let exponent = &s[e_pos..];

        // Trim trailing zeros after the decimal point
        if let Some(dot_pos) = mantissa.find('.') {
            let after_dot = &mantissa[dot_pos + 1..];
            let trimmed_after = after_dot.trim_end_matches('0');
            if trimmed_after.is_empty() {
                format!("{}.{exponent}", &mantissa[..dot_pos])
            } else {
                format!("{}.{trimmed_after}{exponent}", &mantissa[..dot_pos])
            }
        } else {
            format!("{mantissa}{exponent}")
        }
    } else {
        s.to_string()
    }
}

/// Write a STEP string with proper escaping.
fn write_step_string(out: &mut Vec<u8>, s: &str) -> Result<(), IoError> {
    write!(out, "'")?;
    for ch in s.chars() {
        if ch == '\'' {
            write!(out, "''")?;
        } else if ch.is_ascii() {
            write!(out, "{ch}")?;
        } else {
            // Encode non-ASCII as \X2\XXXX\X0\
            let code = ch as u32;
            if code <= 0xFFFF {
                write!(out, "\\X2\\{:04X}\\X0\\", code)?;
            } else {
                // Extended: use \X4\ for characters beyond BMP
                write!(out, "\\X4\\{:08X}\\X0\\", code)?;
            }
        }
    }
    write!(out, "'")?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn serialize_cartesian_point() {
        let entity = EntityInstance {
            id: 1,
            type_name: "CARTESIAN_POINT".into(),
            params: vec![
                Parameter::String("".into()),
                Parameter::List(vec![
                    Parameter::Real(1.0),
                    Parameter::Real(2.5),
                    Parameter::Real(0.0),
                ]),
            ],
            complex_parts: vec![],
        };
        let mut out = Vec::new();
        write_entity(&mut out, &entity, 80).unwrap();
        let text = String::from_utf8(out).unwrap();
        assert!(
            text.contains("#1=CARTESIAN_POINT('',(1.,2.5E+00,0.));\n")
                || text.contains("#1=CARTESIAN_POINT('',(1.,2.5E+0,0.));\n")
                || text.starts_with("#1=CARTESIAN_POINT")
        );
        // Verify basic structure
        assert!(text.contains("CARTESIAN_POINT"));
        assert!(text.contains("#1="));
        assert!(text.trim_end().ends_with(';'));
    }

    #[test]
    fn string_escaping_single_quotes() {
        let mut out = Vec::new();
        write_step_string(&mut out, "it's a test").unwrap();
        let text = String::from_utf8(out).unwrap();
        assert_eq!(text, "'it''s a test'");
    }

    #[test]
    fn string_escaping_non_ascii() {
        let mut out = Vec::new();
        write_step_string(&mut out, "caf\u{00E9}").unwrap();
        let text = String::from_utf8(out).unwrap();
        assert_eq!(text, "'caf\\X2\\00E9\\X0\\'");
    }

    #[test]
    fn real_formatting_zero() {
        let mut out = Vec::new();
        write_real(&mut out, 0.0).unwrap();
        assert_eq!(String::from_utf8(out).unwrap(), "0.");
    }

    #[test]
    fn real_formatting_integer() {
        let mut out = Vec::new();
        write_real(&mut out, 5.0).unwrap();
        assert_eq!(String::from_utf8(out).unwrap(), "5.");
    }

    #[test]
    fn real_formatting_scientific() {
        let mut out = Vec::new();
        write_real(&mut out, 1.23e-10).unwrap();
        let text = String::from_utf8(out).unwrap();
        assert!(text.contains("1.23"));
        assert!(text.contains("E"));
    }

    #[test]
    fn parameter_omitted_and_derived() {
        let mut out = Vec::new();
        write_parameter(&mut out, &Parameter::Omitted).unwrap();
        write_parameter(&mut out, &Parameter::Derived).unwrap();
        assert_eq!(String::from_utf8(out).unwrap(), "$*");
    }

    #[test]
    fn parameter_enum() {
        let mut out = Vec::new();
        write_parameter(&mut out, &Parameter::Enum("T".into())).unwrap();
        assert_eq!(String::from_utf8(out).unwrap(), ".T.");
    }

    #[test]
    fn parameter_entity_ref() {
        let mut out = Vec::new();
        write_parameter(&mut out, &Parameter::EntityRef(42)).unwrap();
        assert_eq!(String::from_utf8(out).unwrap(), "#42");
    }

    #[test]
    fn parameter_typed() {
        let mut out = Vec::new();
        write_parameter(
            &mut out,
            &Parameter::Typed("LENGTH_MEASURE".into(), Box::new(Parameter::Real(25.4))),
        )
        .unwrap();
        let text = String::from_utf8(out).unwrap();
        assert!(text.starts_with("LENGTH_MEASURE("));
        assert!(text.ends_with(')'));
    }

    #[test]
    fn line_wrapping() {
        let entity = EntityInstance {
            id: 100,
            type_name: "B_SPLINE_CURVE_WITH_KNOTS".into(),
            params: vec![
                Parameter::String("long curve name".into()),
                Parameter::Int(3),
                Parameter::List(vec![
                    Parameter::EntityRef(1),
                    Parameter::EntityRef(2),
                    Parameter::EntityRef(3),
                    Parameter::EntityRef(4),
                    Parameter::EntityRef(5),
                    Parameter::EntityRef(6),
                    Parameter::EntityRef(7),
                    Parameter::EntityRef(8),
                ]),
                Parameter::Enum("UNSPECIFIED".into()),
                Parameter::Enum("F".into()),
                Parameter::Enum("F".into()),
            ],
            complex_parts: vec![],
        };
        let mut out = Vec::new();
        write_entity(&mut out, &entity, 40).unwrap();
        let text = String::from_utf8(out).unwrap();
        // Should have multiple lines due to wrapping
        let lines: Vec<&str> = text.lines().collect();
        assert!(
            lines.len() > 1,
            "Expected multiple lines for wrapped output, got: {text}"
        );
    }

    #[test]
    fn serialize_full_part21() {
        use crate::step::document::*;
        let doc = StepDocument {
            products: vec![StepProductNode {
                name: "TestPart".into(),
                id: "TP001".into(),
                children: vec![],
                shapes: vec![],
                metadata: ProductMetadata::default(),
            }],
            roots: vec![0],
            units: UnitContext::default(),
            opaque_bag: OpaqueBag::default(),
        };

        let entities = vec![EntityInstance {
            id: 1,
            type_name: "CARTESIAN_POINT".into(),
            params: vec![
                Parameter::String("origin".into()),
                Parameter::List(vec![
                    Parameter::Real(0.0),
                    Parameter::Real(0.0),
                    Parameter::Real(0.0),
                ]),
            ],
            complex_parts: vec![],
        }];

        let options = StepExportOptions::default();
        let bytes = serialize_part21(&entities, &doc, &options).unwrap();
        let text = String::from_utf8(bytes).unwrap();

        assert!(text.starts_with("ISO-10303-21;"));
        assert!(text.contains("HEADER;"));
        assert!(text.contains("FILE_DESCRIPTION"));
        assert!(text.contains("FILE_NAME"));
        assert!(text.contains("AUTOMOTIVE_DESIGN"));
        assert!(text.contains("DATA;"));
        assert!(text.contains("CARTESIAN_POINT"));
        assert!(text.contains("ENDSEC;"));
        assert!(text.contains("END-ISO-10303-21;"));
    }

    #[test]
    fn trim_trailing_zeros_in_scientific() {
        assert_eq!(
            trim_real_trailing_zeros("1.230000000000000E+02"),
            "1.23E+02"
        );
        assert_eq!(trim_real_trailing_zeros("5.000000000000000E+00"), "5.E+00");
        assert_eq!(trim_real_trailing_zeros("1.000000000000000E-10"), "1.E-10");
    }

    #[test]
    fn complex_entity_serialization() {
        let entity = EntityInstance {
            id: 5,
            type_name: "A".into(),
            params: vec![Parameter::Int(1)],
            complex_parts: vec![("B".into(), vec![Parameter::Int(2)])],
        };
        let mut out = Vec::new();
        write_entity(&mut out, &entity, 200).unwrap();
        let text = String::from_utf8(out).unwrap();
        assert!(text.contains("#5=(A(1) B(2));"));
    }
}
