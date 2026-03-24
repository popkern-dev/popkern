//! Part 21 tokenizer (lexer).
//!
//! Converts raw bytes into a flat sequence of [`Token`] values, stripping
//! whitespace and comments.  This is the first stage of the Part 21 parser
//! pipeline.

/// A single lexical token from a Part 21 file.
#[derive(Debug, Clone, PartialEq)]
pub enum Token {
    /// Entity reference (`#N`).
    EntityRef(u64),
    /// Integer literal.
    Int(i64),
    /// Real (floating-point) literal.
    Real(f64),
    /// String literal (decoded from Part 21 escaping).
    String(String),
    /// Enumeration value (e.g., `.BOTH.`).
    Enum(String),
    /// Keyword / type name (e.g., `CARTESIAN_POINT`).
    Keyword(String),
    /// Left parenthesis `(`.
    LParen,
    /// Right parenthesis `)`.
    RParen,
    /// Comma `,`.
    Comma,
    /// Semicolon `;`.
    Semi,
    /// Equals sign `=`.
    Eq,
    /// Omitted parameter `$`.
    Omitted,
    /// Derived parameter `*`.
    Derived,
}

/// Tokenize a Part 21 byte stream into a vector of tokens.
///
/// Comments (`/* ... */`) and whitespace are silently discarded.
///
/// # Errors
///
/// Returns a human-readable error string if the input contains bytes that
/// cannot be mapped to any token.
pub fn tokenize(input: &[u8]) -> Result<Vec<Token>, String> {
    let mut pos = 0;
    let mut tokens = Vec::new();

    while pos < input.len() {
        // Skip whitespace
        if input[pos].is_ascii_whitespace() {
            pos += 1;
            continue;
        }

        // Skip comments: /* ... */
        if pos + 1 < input.len() && input[pos] == b'/' && input[pos + 1] == b'*' {
            pos += 2;
            loop {
                if pos + 1 >= input.len() {
                    return Err("Unterminated comment".into());
                }
                if input[pos] == b'*' && input[pos + 1] == b'/' {
                    pos += 2;
                    break;
                }
                pos += 1;
            }
            continue;
        }

        // Entity ref: #N
        if input[pos] == b'#' {
            pos += 1;
            let start = pos;
            while pos < input.len() && input[pos].is_ascii_digit() {
                pos += 1;
            }
            if start == pos {
                return Err(format!("Expected digits after '#' at byte {start}"));
            }
            let s = std::str::from_utf8(&input[start..pos])
                .map_err(|e| format!("Invalid UTF-8 in entity ref: {e}"))?;
            let n: u64 = s
                .parse()
                .map_err(|e| format!("Invalid entity ref number: {e}"))?;
            tokens.push(Token::EntityRef(n));
            continue;
        }

        // Number (integer or real): optional '-', digits, optional '.' + digits, optional E/e
        if input[pos].is_ascii_digit() || (input[pos] == b'-' && is_digit_ahead(input, pos + 1)) {
            let start = pos;
            if input[pos] == b'-' {
                pos += 1;
            }
            while pos < input.len() && input[pos].is_ascii_digit() {
                pos += 1;
            }
            let mut is_real = false;
            if pos < input.len() && input[pos] == b'.' {
                is_real = true;
                pos += 1;
                while pos < input.len() && input[pos].is_ascii_digit() {
                    pos += 1;
                }
            }
            if pos < input.len() && (input[pos] == b'E' || input[pos] == b'e') {
                is_real = true;
                pos += 1;
                if pos < input.len() && (input[pos] == b'+' || input[pos] == b'-') {
                    pos += 1;
                }
                while pos < input.len() && input[pos].is_ascii_digit() {
                    pos += 1;
                }
            }
            let s = std::str::from_utf8(&input[start..pos])
                .map_err(|e| format!("Invalid UTF-8 in number: {e}"))?;
            if is_real {
                let v: f64 = s.parse().map_err(|e| format!("Invalid real: {e}"))?;
                tokens.push(Token::Real(v));
            } else {
                let v: i64 = s.parse().map_err(|e| format!("Invalid integer: {e}"))?;
                tokens.push(Token::Int(v));
            }
            continue;
        }

        // String: '...' with '' for escaped quote
        if input[pos] == b'\'' {
            pos += 1;
            let mut bytes = Vec::new();
            loop {
                if pos >= input.len() {
                    return Err("Unterminated string literal".into());
                }
                if input[pos] == b'\'' {
                    if pos + 1 < input.len() && input[pos + 1] == b'\'' {
                        bytes.push(b'\'');
                        pos += 2;
                    } else {
                        pos += 1;
                        break;
                    }
                } else {
                    bytes.push(input[pos]);
                    pos += 1;
                }
            }
            let raw =
                std::str::from_utf8(&bytes).map_err(|e| format!("Invalid UTF-8 in string: {e}"))?;
            let decoded = decode_step_string(raw)?;
            tokens.push(Token::String(decoded));
            continue;
        }

        // Enum: .NAME.
        if input[pos] == b'.' {
            pos += 1;
            let start = pos;
            while pos < input.len() && input[pos] != b'.' {
                pos += 1;
            }
            if pos >= input.len() {
                return Err("Unterminated enumeration literal".into());
            }
            let name = std::str::from_utf8(&input[start..pos])
                .map_err(|e| format!("Invalid UTF-8 in enum: {e}"))?;
            tokens.push(Token::Enum(name.to_string()));
            pos += 1; // skip closing '.'
            continue;
        }

        // Keyword: [A-Z][A-Z0-9_]*
        if input[pos].is_ascii_uppercase() {
            let start = pos;
            pos += 1;
            while pos < input.len()
                && (input[pos].is_ascii_uppercase()
                    || input[pos].is_ascii_digit()
                    || input[pos] == b'_')
            {
                pos += 1;
            }
            let s = std::str::from_utf8(&input[start..pos])
                .map_err(|e| format!("Invalid UTF-8 in keyword: {e}"))?;
            tokens.push(Token::Keyword(s.to_string()));
            continue;
        }

        // Punctuation
        match input[pos] {
            b'(' => {
                tokens.push(Token::LParen);
                pos += 1;
            }
            b')' => {
                tokens.push(Token::RParen);
                pos += 1;
            }
            b',' => {
                tokens.push(Token::Comma);
                pos += 1;
            }
            b';' => {
                tokens.push(Token::Semi);
                pos += 1;
            }
            b'=' => {
                tokens.push(Token::Eq);
                pos += 1;
            }
            b'$' => {
                tokens.push(Token::Omitted);
                pos += 1;
            }
            b'*' => {
                tokens.push(Token::Derived);
                pos += 1;
            }
            ch => {
                return Err(format!(
                    "Unexpected byte 0x{ch:02X} ('{}') at position {pos}",
                    ch as char
                ));
            }
        }
    }

    Ok(tokens)
}

/// Returns true if there is a digit at the given position.
fn is_digit_ahead(input: &[u8], pos: usize) -> bool {
    pos < input.len() && input[pos].is_ascii_digit()
}

/// Decode Part 21 `\X2\XXXX\X0\` Unicode escape sequences in a string.
fn decode_step_string(raw: &str) -> Result<String, String> {
    let mut result = String::with_capacity(raw.len());
    let mut chars = raw.chars().peekable();

    while let Some(ch) = chars.next() {
        if ch == '\\' {
            // Check for \X2\ or \X4\ Unicode escapes
            if chars.peek() == Some(&'X') {
                let mut escape = String::from('\\');
                // Collect until we have enough to identify the escape
                if let Some(x) = chars.next() {
                    escape.push(x);
                }
                if let Some(digit) = chars.next() {
                    escape.push(digit);
                }
                if let Some(backslash) = chars.next() {
                    escape.push(backslash);
                }

                if escape == "\\X2\\" {
                    // Read hex chars until \X0\
                    let mut hex_buf = String::new();
                    loop {
                        match chars.next() {
                            Some('\\') => {
                                // Check for \X0\
                                let mut end_seq = String::from('\\');
                                for _ in 0..3 {
                                    if let Some(c) = chars.next() {
                                        end_seq.push(c);
                                    }
                                }
                                if end_seq == "\\X0\\" {
                                    break;
                                }
                                // Not the terminator, treat as hex
                                hex_buf.push_str(&end_seq);
                            }
                            Some(c) => hex_buf.push(c),
                            None => {
                                return Err("Unterminated \\X2\\ Unicode escape in string".into());
                            }
                        }
                    }
                    // Each 4-hex-digit group is a Unicode code point
                    let mut i = 0;
                    while i + 4 <= hex_buf.len() {
                        let hex = &hex_buf[i..i + 4];
                        let code = u32::from_str_radix(hex, 16)
                            .map_err(|e| format!("Invalid hex in \\X2\\ escape '{hex}': {e}"))?;
                        let c = char::from_u32(code)
                            .ok_or_else(|| format!("Invalid Unicode code point U+{code:04X}"))?;
                        result.push(c);
                        i += 4;
                    }
                    if i != hex_buf.len() {
                        return Err(format!(
                            "\\X2\\ escape has {} leftover hex chars (not a multiple of 4)",
                            hex_buf.len() - i
                        ));
                    }
                } else {
                    // Unknown escape, pass through
                    result.push_str(&escape);
                }
            } else {
                // Unknown escape, pass through
                result.push(ch);
                if let Some(next) = chars.next() {
                    result.push(next);
                }
            }
        } else {
            result.push(ch);
        }
    }

    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tokenize_entity_ref() {
        let tokens = tokenize(b"#42").unwrap();
        assert_eq!(tokens, vec![Token::EntityRef(42)]);
    }

    #[test]
    fn tokenize_integer() {
        let tokens = tokenize(b"123").unwrap();
        assert_eq!(tokens, vec![Token::Int(123)]);
    }

    #[test]
    fn tokenize_negative_integer() {
        let tokens = tokenize(b"-5").unwrap();
        assert_eq!(tokens, vec![Token::Int(-5)]);
    }

    #[test]
    fn tokenize_real() {
        let tokens = tokenize(b"3.14").unwrap();
        assert_eq!(tokens, vec![Token::Real(3.14)]);
    }

    #[test]
    fn tokenize_real_scientific() {
        let tokens = tokenize(b"1.0E-3").unwrap();
        assert_eq!(tokens, vec![Token::Real(1.0e-3)]);
    }

    #[test]
    fn tokenize_string() {
        let tokens = tokenize(b"'hello'").unwrap();
        assert_eq!(tokens, vec![Token::String("hello".into())]);
    }

    #[test]
    fn tokenize_string_escaped_quote() {
        let tokens = tokenize(b"'it''s'").unwrap();
        assert_eq!(tokens, vec![Token::String("it's".into())]);
    }

    #[test]
    fn tokenize_enum() {
        let tokens = tokenize(b".BOTH.").unwrap();
        assert_eq!(tokens, vec![Token::Enum("BOTH".into())]);
    }

    #[test]
    fn tokenize_keyword() {
        let tokens = tokenize(b"CARTESIAN_POINT").unwrap();
        assert_eq!(tokens, vec![Token::Keyword("CARTESIAN_POINT".into())]);
    }

    #[test]
    fn tokenize_punctuation() {
        let tokens = tokenize(b"(,);=$*").unwrap();
        assert_eq!(
            tokens,
            vec![
                Token::LParen,
                Token::Comma,
                Token::RParen,
                Token::Semi,
                Token::Eq,
                Token::Omitted,
                Token::Derived,
            ]
        );
    }

    #[test]
    fn tokenize_full_entity_line() {
        let input = b"#1 = CARTESIAN_POINT('', (0.0, 0.0, 0.0));";
        let tokens = tokenize(input).unwrap();
        assert_eq!(
            tokens,
            vec![
                Token::EntityRef(1),
                Token::Eq,
                Token::Keyword("CARTESIAN_POINT".into()),
                Token::LParen,
                Token::String("".into()),
                Token::Comma,
                Token::LParen,
                Token::Real(0.0),
                Token::Comma,
                Token::Real(0.0),
                Token::Comma,
                Token::Real(0.0),
                Token::RParen,
                Token::RParen,
                Token::Semi,
            ]
        );
    }

    #[test]
    fn tokenize_skips_comments() {
        let tokens = tokenize(b"/* comment */ #1").unwrap();
        assert_eq!(tokens, vec![Token::EntityRef(1)]);
    }

    #[test]
    fn tokenize_unicode_escape() {
        let tokens = tokenize(b"'caf\\X2\\00E9\\X0\\'").unwrap();
        assert_eq!(tokens, vec![Token::String("caf\u{00E9}".into())]);
    }
}
