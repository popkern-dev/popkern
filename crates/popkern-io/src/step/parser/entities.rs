//! Part 21 entity parser.
//!
//! Parses a tokenized Part 21 file into a [`ParseResult`] containing the
//! header, entity instances, and any parse errors encountered.

use super::tokenizer::{Token, tokenize};
use super::types::{EntityId, EntityInstance, Parameter, StepHeader};

/// A single parse error with its token offset.
#[derive(Debug, Clone)]
pub struct ParseError {
    /// Token offset where the error occurred.
    pub offset: usize,
    /// Human-readable description of the error.
    pub message: String,
}

/// Result of parsing a complete Part 21 file.
#[derive(Debug)]
pub struct ParseResult {
    /// Parsed HEADER section.
    pub header: StepHeader,
    /// Parsed entity instances from the DATA section.
    pub entities: Vec<EntityInstance>,
    /// Non-fatal parse errors encountered during parsing.
    pub parse_errors: Vec<ParseError>,
}

/// Parse a complete Part 21 file from bytes.
///
/// Tokenizes the input and then walks the token stream to extract the HEADER
/// and DATA sections.  Malformed entities are skipped and recorded in
/// [`ParseResult::parse_errors`].
///
/// # Errors
///
/// Returns an error string if the file structure is fundamentally broken
/// (e.g., missing ISO sentinel or sections).
pub fn parse_step_file(input: &[u8]) -> Result<ParseResult, String> {
    let tokens = tokenize(input)?;
    let mut pos = 0;
    let len = tokens.len();

    // Expect: ISO-10303-21 ;
    expect_keyword(&tokens, &mut pos, "ISO-10303-21")?;
    expect_token(&tokens, &mut pos, &Token::Semi, "';' after ISO sentinel")?;

    // Expect: HEADER ;
    expect_keyword(&tokens, &mut pos, "HEADER")?;
    expect_token(&tokens, &mut pos, &Token::Semi, "';' after HEADER")?;

    // Parse header entities until ENDSEC ;
    let header = parse_header(&tokens, &mut pos)?;

    // Expect: DATA ;
    expect_keyword(&tokens, &mut pos, "DATA")?;
    expect_token(&tokens, &mut pos, &Token::Semi, "';' after DATA")?;

    // Parse data entities until ENDSEC ;
    let mut entities = Vec::new();
    let mut parse_errors = Vec::new();

    while pos < len {
        if matches!(tokens.get(pos), Some(Token::Keyword(k)) if k == "ENDSEC") {
            pos += 1; // skip ENDSEC
            expect_token(&tokens, &mut pos, &Token::Semi, "';' after ENDSEC")?;
            break;
        }
        match parse_entity(&tokens, &mut pos) {
            Ok(entity) => entities.push(entity),
            Err(msg) => {
                parse_errors.push(ParseError {
                    offset: pos,
                    message: msg,
                });
                skip_to_semi(&tokens, &mut pos);
            }
        }
    }

    // Expect: END-ISO-10303-21 ;
    expect_keyword(&tokens, &mut pos, "END-ISO-10303-21")?;
    expect_token(&tokens, &mut pos, &Token::Semi, "';' after END sentinel")?;

    Ok(ParseResult {
        header,
        entities,
        parse_errors,
    })
}

/// Expect a specific keyword at the current position.
fn expect_keyword(tokens: &[Token], pos: &mut usize, name: &str) -> Result<(), String> {
    match tokens.get(*pos) {
        Some(Token::Keyword(k)) if k == name => {
            *pos += 1;
            Ok(())
        }
        other => Err(format!(
            "expected keyword '{name}' at token {pos}, got {other:?}",
            pos = *pos
        )),
    }
}

/// Expect a specific token at the current position.
fn expect_token(
    tokens: &[Token],
    pos: &mut usize,
    expected: &Token,
    ctx: &str,
) -> Result<(), String> {
    match tokens.get(*pos) {
        Some(t) if *t == *expected => {
            *pos += 1;
            Ok(())
        }
        other => Err(format!(
            "expected {ctx} at token {pos}, got {other:?}",
            pos = *pos
        )),
    }
}

/// Skip tokens until (and including) the next semicolon.
fn skip_to_semi(tokens: &[Token], pos: &mut usize) {
    while *pos < tokens.len() {
        if tokens[*pos] == Token::Semi {
            *pos += 1;
            return;
        }
        *pos += 1;
    }
}

/// Parse the HEADER section, returning a populated [`StepHeader`].
fn parse_header(tokens: &[Token], pos: &mut usize) -> Result<StepHeader, String> {
    let mut header = StepHeader::default();

    while *pos < tokens.len() {
        if matches!(tokens.get(*pos), Some(Token::Keyword(k)) if k == "ENDSEC") {
            *pos += 1; // skip ENDSEC
            expect_token(tokens, pos, &Token::Semi, "';' after ENDSEC")?;
            return Ok(header);
        }

        // Each header entity: Keyword LParen params RParen Semi
        let kw = match tokens.get(*pos) {
            Some(Token::Keyword(k)) => {
                let k = k.clone();
                *pos += 1;
                k
            }
            other => {
                return Err(format!(
                    "expected header keyword at token {pos}, got {other:?}",
                    pos = *pos
                ));
            }
        };

        expect_token(tokens, pos, &Token::LParen, "'(' after header keyword")?;
        let params = parse_param_list(tokens, pos)?;
        expect_token(tokens, pos, &Token::RParen, "')' after header params")?;
        expect_token(tokens, pos, &Token::Semi, "';' after header entity")?;

        match kw.as_str() {
            "FILE_DESCRIPTION" => {
                // FILE_DESCRIPTION((strings...), impl_level)
                if let Some(Parameter::List(strings)) = params.first() {
                    header.description = strings
                        .iter()
                        .filter_map(|p| {
                            if let Parameter::String(s) = p {
                                Some(s.clone())
                            } else {
                                None
                            }
                        })
                        .collect();
                }
                if let Some(Parameter::String(level)) = params.get(1) {
                    header.implementation_level = level.clone();
                }
            }
            "FILE_NAME" => {
                if let Some(Parameter::String(name)) = params.first() {
                    header.name = name.clone();
                }
            }
            "FILE_SCHEMA" => {
                // FILE_SCHEMA(('SCHEMA_NAME'))
                if let Some(Parameter::List(schemas)) = params.first() {
                    header.schemas = schemas
                        .iter()
                        .filter_map(|p| {
                            if let Parameter::String(s) = p {
                                Some(s.clone())
                            } else {
                                None
                            }
                        })
                        .collect();
                }
            }
            _ => { /* ignore unknown header entities */ }
        }
    }

    Err("unexpected end of tokens in HEADER section".into())
}

/// Parse a single entity instance from the DATA section.
fn parse_entity(tokens: &[Token], pos: &mut usize) -> Result<EntityInstance, String> {
    // EntityRef(n)
    let id = match tokens.get(*pos) {
        Some(Token::EntityRef(n)) => {
            let n = *n;
            *pos += 1;
            n
        }
        other => {
            return Err(format!(
                "expected entity ref at token {pos}, got {other:?}",
                pos = *pos
            ));
        }
    };

    // Eq
    expect_token(tokens, pos, &Token::Eq, "'=' after entity ref")?;

    // Check for complex entity: LParen (Keyword LParen params RParen)+ RParen Semi
    if matches!(tokens.get(*pos), Some(Token::LParen)) {
        return parse_complex_entity(tokens, pos, id);
    }

    // Simple entity: Keyword(type_name) LParen params RParen Semi
    let type_name = match tokens.get(*pos) {
        Some(Token::Keyword(k)) => {
            let k = k.clone();
            *pos += 1;
            k
        }
        other => {
            return Err(format!(
                "expected type keyword at token {pos}, got {other:?}",
                pos = *pos
            ));
        }
    };

    expect_token(tokens, pos, &Token::LParen, "'(' after type name")?;
    let params = parse_param_list(tokens, pos)?;
    expect_token(tokens, pos, &Token::RParen, "')' after entity params")?;
    expect_token(tokens, pos, &Token::Semi, "';' after entity")?;

    Ok(EntityInstance {
        id,
        type_name,
        params,
        complex_parts: vec![],
    })
}

/// Parse a complex (AND-connected) entity instance.
fn parse_complex_entity(
    tokens: &[Token],
    pos: &mut usize,
    id: EntityId,
) -> Result<EntityInstance, String> {
    // Skip the opening LParen
    expect_token(tokens, pos, &Token::LParen, "'(' for complex entity")?;

    let mut parts: Vec<(String, Vec<Parameter>)> = Vec::new();

    while *pos < tokens.len() {
        if matches!(tokens.get(*pos), Some(Token::RParen)) {
            *pos += 1;
            break;
        }

        let kw = match tokens.get(*pos) {
            Some(Token::Keyword(k)) => {
                let k = k.clone();
                *pos += 1;
                k
            }
            other => {
                return Err(format!(
                    "expected type keyword in complex entity at token {pos}, got {other:?}",
                    pos = *pos
                ));
            }
        };

        expect_token(tokens, pos, &Token::LParen, "'(' after complex type name")?;
        let params = parse_param_list(tokens, pos)?;
        expect_token(tokens, pos, &Token::RParen, "')' after complex type params")?;

        parts.push((kw, params));
    }

    expect_token(tokens, pos, &Token::Semi, "';' after complex entity")?;

    let (type_name, params) = if let Some(first) = parts.first() {
        (first.0.clone(), first.1.clone())
    } else {
        return Err("empty complex entity".into());
    };

    let complex_parts = if parts.len() > 1 {
        parts[1..].to_vec()
    } else {
        vec![]
    };

    Ok(EntityInstance {
        id,
        type_name,
        params,
        complex_parts,
    })
}

/// Parse a comma-separated list of parameters (does NOT consume the closing RParen).
fn parse_param_list(tokens: &[Token], pos: &mut usize) -> Result<Vec<Parameter>, String> {
    let mut params = Vec::new();

    // Empty parameter list
    if matches!(tokens.get(*pos), Some(Token::RParen)) {
        return Ok(params);
    }

    params.push(parse_parameter(tokens, pos)?);

    while matches!(tokens.get(*pos), Some(Token::Comma)) {
        *pos += 1; // skip comma
        params.push(parse_parameter(tokens, pos)?);
    }

    Ok(params)
}

/// Parse a single parameter value (recursive for lists and typed params).
fn parse_parameter(tokens: &[Token], pos: &mut usize) -> Result<Parameter, String> {
    match tokens.get(*pos) {
        Some(Token::Int(n)) => {
            let n = *n;
            *pos += 1;
            Ok(Parameter::Int(n))
        }
        Some(Token::Real(f)) => {
            let f = *f;
            *pos += 1;
            Ok(Parameter::Real(f))
        }
        Some(Token::String(s)) => {
            let s = s.clone();
            *pos += 1;
            Ok(Parameter::String(s))
        }
        Some(Token::Enum(s)) => {
            let s = s.clone();
            *pos += 1;
            Ok(Parameter::Enum(s))
        }
        Some(Token::EntityRef(id)) => {
            let id = *id;
            *pos += 1;
            Ok(Parameter::EntityRef(id))
        }
        Some(Token::Omitted) => {
            *pos += 1;
            Ok(Parameter::Omitted)
        }
        Some(Token::Derived) => {
            *pos += 1;
            Ok(Parameter::Derived)
        }
        Some(Token::LParen) => {
            *pos += 1;
            let items = parse_param_list(tokens, pos)?;
            expect_token(tokens, pos, &Token::RParen, "')' in list parameter")?;
            Ok(Parameter::List(items))
        }
        Some(Token::Keyword(name)) => {
            // Typed parameter: Keyword LParen params RParen
            let name = name.clone();
            *pos += 1;
            expect_token(tokens, pos, &Token::LParen, "'(' after typed param name")?;
            let items = parse_param_list(tokens, pos)?;
            expect_token(tokens, pos, &Token::RParen, "')' after typed param")?;
            Ok(Parameter::Typed(name, Box::new(Parameter::List(items))))
        }
        other => Err(format!(
            "unexpected token at position {pos}: {other:?}",
            pos = *pos
        )),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const SIMPLE_FILE: &[u8] = b"ISO-10303-21;
HEADER;
FILE_DESCRIPTION(('test'), '2;1');
FILE_NAME('test.step', '2026-03-24', (''), (''), '', '', '');
FILE_SCHEMA(('CONFIG_CONTROL_DESIGN'));
ENDSEC;
DATA;
#1 = CARTESIAN_POINT('origin', (0.0, 0.0, 0.0));
ENDSEC;
END-ISO-10303-21;
";

    #[test]
    fn parse_simple_entity() {
        let result = parse_step_file(SIMPLE_FILE).unwrap();
        assert_eq!(result.entities.len(), 1);
        let e = &result.entities[0];
        assert_eq!(e.id, 1);
        assert_eq!(e.type_name, "CARTESIAN_POINT");
        assert_eq!(e.params.len(), 2);
        assert_eq!(e.params[0], Parameter::String("origin".into()));
        if let Parameter::List(coords) = &e.params[1] {
            assert_eq!(coords.len(), 3);
        } else {
            panic!("expected list parameter for coordinates");
        }
    }

    #[test]
    fn parse_header_extracts_schema() {
        let result = parse_step_file(SIMPLE_FILE).unwrap();
        assert_eq!(result.header.schemas, vec!["CONFIG_CONTROL_DESIGN"]);
        assert_eq!(result.header.implementation_level, "2;1");
        assert_eq!(result.header.name, "test.step");
    }

    #[test]
    fn parse_multiple_entities() {
        let input = b"ISO-10303-21;
HEADER;
FILE_DESCRIPTION(('multi'), '2;1');
FILE_NAME('multi.step', '', (''), (''), '', '', '');
FILE_SCHEMA(('AUTOMOTIVE_DESIGN'));
ENDSEC;
DATA;
#1 = CARTESIAN_POINT('', (0.0, 0.0, 0.0));
#2 = DIRECTION('', (0.0, 0.0, 1.0));
#3 = AXIS2_PLACEMENT_3D('', #1, #2, $);
ENDSEC;
END-ISO-10303-21;
";
        let result = parse_step_file(input).unwrap();
        assert_eq!(result.entities.len(), 3);
        assert_eq!(result.entities[2].type_name, "AXIS2_PLACEMENT_3D");
        // Check entity ref and omitted params
        assert_eq!(result.entities[2].params[1], Parameter::EntityRef(1));
        assert_eq!(result.entities[2].params[2], Parameter::EntityRef(2));
        assert_eq!(result.entities[2].params[3], Parameter::Omitted);
    }

    #[test]
    fn parse_complex_entity() {
        let input = b"ISO-10303-21;
HEADER;
FILE_DESCRIPTION(('complex'), '2;1');
FILE_NAME('complex.step', '', (''), (''), '', '', '');
FILE_SCHEMA(('AUTOMOTIVE_DESIGN'));
ENDSEC;
DATA;
#1 = (GEOMETRIC_REPRESENTATION_CONTEXT(3) GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT((#2)) REPRESENTATION_CONTEXT('plan', ''));
ENDSEC;
END-ISO-10303-21;
";
        let result = parse_step_file(input).unwrap();
        assert_eq!(result.entities.len(), 1);
        let e = &result.entities[0];
        assert_eq!(e.id, 1);
        assert_eq!(e.type_name, "GEOMETRIC_REPRESENTATION_CONTEXT");
        assert_eq!(e.params, vec![Parameter::Int(3)]);
        assert_eq!(e.complex_parts.len(), 2);
        assert_eq!(e.complex_parts[0].0, "GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT");
        assert_eq!(e.complex_parts[1].0, "REPRESENTATION_CONTEXT");
    }
}
