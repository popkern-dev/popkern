//! Tokenizer for the expression parser.
//!
//! Converts a string into a flat sequence of [`Token`]s, resolving unit
//! suffixes and building dotted identifiers.

use crate::parse::ParseError;
use crate::unit::Unit;

// ---------------------------------------------------------------------------
// Token types
// ---------------------------------------------------------------------------

/// Tokens produced by the expression lexer.
#[derive(Debug, Clone, PartialEq)]
pub(super) enum Token {
    /// A numeric literal in SI units (unit already applied).
    Number(f64),
    /// A unit-suffixed numeric literal: its SI value and its unit.
    Quantity(f64, Unit),
    /// An identifier, possibly dotted (e.g. `"extrude1.distance"`).
    Ident(String),
    /// `+`
    Plus,
    /// `-`
    Minus,
    /// `*`
    Star,
    /// `/`
    Slash,
    /// `>`
    Gt,
    /// `<`
    Lt,
    /// `>=`
    Ge,
    /// `<=`
    Le,
    /// `==`
    EqEq,
    /// `!=`
    BangEq,
    /// `&&`
    AndAnd,
    /// `||`
    OrOr,
    /// `!`
    Bang,
    /// `?`
    Question,
    /// `:`
    Colon,
    /// `(`
    LParen,
    /// `)`
    RParen,
    /// `,`
    Comma,
}

// ---------------------------------------------------------------------------
// Lexer
// ---------------------------------------------------------------------------

struct Lexer<'a> {
    chars: &'a [char],
    pos: usize,
}

impl<'a> Lexer<'a> {
    fn new(chars: &'a [char]) -> Self {
        Self { chars, pos: 0 }
    }

    fn peek(&self) -> Option<char> {
        self.chars.get(self.pos).copied()
    }

    fn peek2(&self) -> Option<char> {
        self.chars.get(self.pos + 1).copied()
    }

    fn advance(&mut self) -> Option<char> {
        let ch = self.chars.get(self.pos).copied();
        if ch.is_some() {
            self.pos += 1;
        }
        ch
    }

    fn skip_whitespace(&mut self) {
        while let Some(c) = self.peek() {
            if c.is_ascii_whitespace() {
                self.pos += 1;
            } else {
                break;
            }
        }
    }

    /// Lex the next token and return it along with its start column.
    fn next_token(&mut self) -> Result<Option<(Token, usize)>, ParseError> {
        self.skip_whitespace();
        let col = self.pos;

        let ch = match self.peek() {
            None => return Ok(None),
            Some(c) => c,
        };

        if ch.is_ascii_digit() {
            return self.lex_number(col).map(|t| Some((t, col)));
        }

        if ch.is_alphabetic() || ch == '_' {
            return self.lex_ident(col).map(|t| Some((t, col)));
        }

        self.advance();

        let tok = match ch {
            '+' => Token::Plus,
            '-' => Token::Minus,
            '*' => Token::Star,
            '/' => Token::Slash,
            '(' => Token::LParen,
            ')' => Token::RParen,
            ',' => Token::Comma,
            '?' => Token::Question,
            ':' => Token::Colon,
            '!' => {
                if self.peek() == Some('=') {
                    self.advance();
                    Token::BangEq
                } else {
                    Token::Bang
                }
            }
            '>' => {
                if self.peek() == Some('=') {
                    self.advance();
                    Token::Ge
                } else {
                    Token::Gt
                }
            }
            '<' => {
                if self.peek() == Some('=') {
                    self.advance();
                    Token::Le
                } else {
                    Token::Lt
                }
            }
            '=' => {
                if self.peek() == Some('=') {
                    self.advance();
                    Token::EqEq
                } else {
                    return Err(ParseError::new(col, "unexpected '='; did you mean '=='?"));
                }
            }
            '&' => {
                if self.peek() == Some('&') {
                    self.advance();
                    Token::AndAnd
                } else {
                    return Err(ParseError::new(col, "expected '&&'"));
                }
            }
            '|' => {
                if self.peek() == Some('|') {
                    self.advance();
                    Token::OrOr
                } else {
                    return Err(ParseError::new(col, "expected '||'"));
                }
            }
            other => {
                return Err(ParseError::new(
                    col,
                    format!("unexpected character '{other}'"),
                ));
            }
        };

        Ok(Some((tok, col)))
    }

    fn lex_number(&mut self, col: usize) -> Result<Token, ParseError> {
        let mut s = String::new();
        while let Some(c) = self.peek() {
            if c.is_ascii_digit() {
                s.push(c);
                self.advance();
            } else {
                break;
            }
        }
        if self.peek() == Some('.') && self.peek2().map(|c| c.is_ascii_digit()).unwrap_or(false) {
            s.push('.');
            self.advance();
            while let Some(c) = self.peek() {
                if c.is_ascii_digit() {
                    s.push(c);
                    self.advance();
                } else {
                    break;
                }
            }
        }
        let raw: f64 = s
            .parse()
            .map_err(|_| ParseError::new(col, "invalid number"))?;

        Ok(match self.try_lex_unit() {
            Some(u) => Token::Quantity(u.to_si(raw), u),
            None => Token::Number(raw),
        })
    }

    /// Try to consume a unit suffix immediately after a number (no whitespace).
    fn try_lex_unit(&mut self) -> Option<Unit> {
        let remaining: String = self.chars[self.pos..].iter().collect();

        // Longest match first.
        let units: &[(&str, Unit)] = &[
            ("mm", Unit::Millimeters),
            ("cm", Unit::Centimeters),
            ("deg", Unit::Degrees),
            ("rad", Unit::Radians),
            ("in", Unit::Inches),
            ("m", Unit::Meters),
        ];
        for (suffix, unit) in units {
            if remaining.starts_with(suffix) {
                let next_pos = self.pos + suffix.len();
                let next_char = self.chars.get(next_pos).copied();
                // Only consume if followed by a non-alphanumeric character or end.
                if next_char
                    .map(|c| !c.is_alphanumeric() && c != '_')
                    .unwrap_or(true)
                {
                    self.pos += suffix.len();
                    return Some(*unit);
                }
            }
        }
        None
    }

    fn lex_ident(&mut self, _col: usize) -> Result<Token, ParseError> {
        let mut s = String::new();
        while let Some(c) = self.peek() {
            if c.is_alphanumeric() || c == '_' {
                s.push(c);
                self.advance();
            } else {
                break;
            }
        }
        // Handle dotted names: "extrude1.distance" or "body.centroid.x"
        while self.peek() == Some('.')
            && self
                .peek2()
                .map(|c| c.is_alphabetic() || c == '_')
                .unwrap_or(false)
        {
            s.push('.');
            self.advance(); // consume '.'
            while let Some(c) = self.peek() {
                if c.is_alphanumeric() || c == '_' {
                    s.push(c);
                    self.advance();
                } else {
                    break;
                }
            }
        }
        Ok(Token::Ident(s))
    }
}

// ---------------------------------------------------------------------------
// Public (crate-internal) entry point
// ---------------------------------------------------------------------------

/// Tokenize the entire input string into a `Vec<(Token, column)>`.
pub(super) fn tokenize(input: &str) -> Result<Vec<(Token, usize)>, ParseError> {
    let chars: Vec<char> = input.chars().collect();
    let mut lexer = Lexer::new(&chars);
    let mut tokens = Vec::new();
    loop {
        match lexer.next_token()? {
            None => break,
            Some(tok) => tokens.push(tok),
        }
    }
    Ok(tokens)
}
