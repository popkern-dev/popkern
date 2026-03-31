//! Recursive-descent parser for parameter expressions.
//!
//! Converts a string like `"width * 2 + 3mm"` into an [`Expression`] AST,
//! resolving identifier names through a [`ParameterTable`].
//!
//! The module is split internally: [`parse_lexer`] handles tokenization, this
//! file contains the recursive-descent parser and the public API.

#[path = "parse_lexer.rs"]
mod parse_lexer;

use parse_lexer::{Token, tokenize};

use crate::expression::{BinOp, CmpOp, Expression, UnaryFunc};
use crate::parameter::ParameterTable;

/// A parse error with column information.
#[derive(Debug, thiserror::Error)]
#[error("parse error at column {column}: {message}")]
pub struct ParseError {
    /// Zero-based column offset where the error was detected.
    pub column: usize,
    /// Human-readable description of the error.
    pub message: String,
}

impl ParseError {
    pub(crate) fn new(column: usize, message: impl Into<String>) -> Self {
        Self {
            column,
            message: message.into(),
        }
    }
}

// ---------------------------------------------------------------------------
// Recursive-descent parser
// ---------------------------------------------------------------------------

struct Parser<'a> {
    tokens: Vec<(Token, usize)>,
    pos: usize,
    table: &'a ParameterTable,
    input_len: usize,
}

impl<'a> Parser<'a> {
    fn new(tokens: Vec<(Token, usize)>, table: &'a ParameterTable, input_len: usize) -> Self {
        Self {
            tokens,
            pos: 0,
            table,
            input_len,
        }
    }

    fn peek(&self) -> Option<&Token> {
        self.tokens.get(self.pos).map(|(t, _)| t)
    }

    fn current_col(&self) -> usize {
        self.tokens
            .get(self.pos)
            .map(|(_, c)| *c)
            .unwrap_or(self.input_len)
    }

    fn advance(&mut self) -> Option<&Token> {
        let tok = self.tokens.get(self.pos).map(|(t, _)| t);
        if tok.is_some() {
            self.pos += 1;
        }
        tok
    }

    fn consume_ident(&mut self) -> Result<(String, usize), ParseError> {
        let col = self.current_col();
        if let Some((Token::Ident(s), c)) = self.tokens.get(self.pos) {
            let name = s.clone();
            let col = *c;
            self.pos += 1;
            Ok((name, col))
        } else {
            Err(ParseError::new(col, "expected identifier"))
        }
    }

    // -- Grammar productions --------------------------------------------------

    fn parse_expr(&mut self) -> Result<Expression, ParseError> {
        self.parse_ternary()
    }

    fn parse_ternary(&mut self) -> Result<Expression, ParseError> {
        let cond = self.parse_logic()?;
        if self.peek() == Some(&Token::Question) {
            self.advance(); // consume '?'
            let then_expr = self.parse_expr()?;
            let colon_col = self.current_col();
            match self.peek() {
                Some(Token::Colon) => {
                    self.advance();
                }
                _ => return Err(ParseError::new(colon_col, "expected ':' in ternary")),
            }
            let else_expr = self.parse_expr()?;
            return Ok(Expression::Conditional(
                Box::new(cond),
                Box::new(then_expr),
                Box::new(else_expr),
            ));
        }
        Ok(cond)
    }

    fn parse_logic(&mut self) -> Result<Expression, ParseError> {
        let mut left = self.parse_compare()?;
        loop {
            let op = match self.peek() {
                Some(Token::AndAnd) => BinOp::And,
                Some(Token::OrOr) => BinOp::Or,
                _ => break,
            };
            self.advance();
            let right = self.parse_compare()?;
            left = Expression::BinOp(Box::new(left), op, Box::new(right));
        }
        Ok(left)
    }

    fn parse_compare(&mut self) -> Result<Expression, ParseError> {
        let left = self.parse_addition()?;
        let op = match self.peek() {
            Some(Token::Gt) => CmpOp::Gt,
            Some(Token::Lt) => CmpOp::Lt,
            Some(Token::Ge) => CmpOp::Ge,
            Some(Token::Le) => CmpOp::Le,
            Some(Token::EqEq) => CmpOp::Eq,
            Some(Token::BangEq) => CmpOp::Ne,
            _ => return Ok(left),
        };
        self.advance();
        let right = self.parse_addition()?;
        Ok(Expression::Compare(Box::new(left), op, Box::new(right)))
    }

    fn parse_addition(&mut self) -> Result<Expression, ParseError> {
        let mut left = self.parse_multiply()?;
        loop {
            let op = match self.peek() {
                Some(Token::Plus) => BinOp::Add,
                Some(Token::Minus) => BinOp::Sub,
                _ => break,
            };
            self.advance();
            let right = self.parse_multiply()?;
            left = Expression::BinOp(Box::new(left), op, Box::new(right));
        }
        Ok(left)
    }

    fn parse_multiply(&mut self) -> Result<Expression, ParseError> {
        let mut left = self.parse_unary()?;
        loop {
            let op = match self.peek() {
                Some(Token::Star) => BinOp::Mul,
                Some(Token::Slash) => BinOp::Div,
                _ => break,
            };
            self.advance();
            let right = self.parse_unary()?;
            left = Expression::BinOp(Box::new(left), op, Box::new(right));
        }
        Ok(left)
    }

    fn parse_unary(&mut self) -> Result<Expression, ParseError> {
        match self.peek() {
            Some(Token::Minus) => {
                self.advance();
                let inner = self.parse_unary()?;
                Ok(Expression::Func(UnaryFunc::Neg, Box::new(inner)))
            }
            Some(Token::Bang) => {
                self.advance();
                let inner = self.parse_unary()?;
                Ok(Expression::Func(UnaryFunc::Not, Box::new(inner)))
            }
            _ => self.parse_call(),
        }
    }

    fn parse_call(&mut self) -> Result<Expression, ParseError> {
        // If we see an IDENT followed by '(', it's a function call.
        if let Some(Token::Ident(_)) = self.peek()
            && self.tokens.get(self.pos + 1).map(|(t, _)| t) == Some(&Token::LParen)
        {
            let (name, col) = self.consume_ident()?;
            self.advance(); // consume '('
            return self.parse_call_args(&name, col);
        }
        self.parse_primary()
    }

    /// Parse arguments for a function call (opening '(' already consumed).
    fn parse_call_args(&mut self, name: &str, col: usize) -> Result<Expression, ParseError> {
        match name {
            "min" | "max" => {
                let a = self.parse_expr()?;
                let comma_col = self.current_col();
                match self.peek() {
                    Some(Token::Comma) => {
                        self.advance();
                    }
                    _ => {
                        return Err(ParseError::new(
                            comma_col,
                            format!("{name}() requires exactly 2 arguments"),
                        ));
                    }
                }
                let b = self.parse_expr()?;
                let rparen_col = self.current_col();
                match self.peek() {
                    Some(Token::RParen) => {
                        self.advance();
                    }
                    _ => {
                        return Err(ParseError::new(rparen_col, "expected ')'"));
                    }
                }
                let op = if name == "min" {
                    BinOp::Min
                } else {
                    BinOp::Max
                };
                Ok(Expression::BinOp(Box::new(a), op, Box::new(b)))
            }
            _ => self.parse_single_arg_call(name, col),
        }
    }

    /// Parse a single-argument function call (opening '(' already consumed).
    fn parse_single_arg_call(&mut self, name: &str, col: usize) -> Result<Expression, ParseError> {
        let func = match name {
            "sin" => UnaryFunc::Sin,
            "cos" => UnaryFunc::Cos,
            "tan" => UnaryFunc::Tan,
            "asin" => UnaryFunc::Asin,
            "acos" => UnaryFunc::Acos,
            "atan" => UnaryFunc::Atan,
            "sqrt" => UnaryFunc::Sqrt,
            "abs" => UnaryFunc::Abs,
            "floor" => UnaryFunc::Floor,
            "ceil" => UnaryFunc::Ceil,
            "round" => UnaryFunc::Round,
            "ln" => UnaryFunc::Ln,
            _ => {
                return Err(ParseError::new(col, format!("unknown function '{name}'")));
            }
        };
        let arg = self.parse_expr()?;
        let extra_col = self.current_col();
        if self.peek() == Some(&Token::Comma) {
            return Err(ParseError::new(
                extra_col,
                format!("{name}() takes exactly 1 argument"),
            ));
        }
        let rparen_col = self.current_col();
        match self.peek() {
            Some(Token::RParen) => {
                self.advance();
            }
            _ => {
                return Err(ParseError::new(rparen_col, "expected ')'"));
            }
        }
        Ok(Expression::Func(func, Box::new(arg)))
    }

    fn parse_primary(&mut self) -> Result<Expression, ParseError> {
        let col = self.current_col();
        match self.peek().cloned() {
            Some(Token::Number(v)) => {
                self.advance();
                Ok(Expression::Literal(v))
            }
            Some(Token::Ident(name)) => {
                self.advance();
                match self.table.resolve_name(&name) {
                    Some(id) => Ok(Expression::Param(id)),
                    None => Err(ParseError::new(col, format!("unknown parameter '{name}'"))),
                }
            }
            Some(Token::LParen) => {
                self.advance(); // consume '('
                let inner = self.parse_expr()?;
                let rparen_col = self.current_col();
                match self.peek() {
                    Some(Token::RParen) => {
                        self.advance();
                    }
                    _ => return Err(ParseError::new(rparen_col, "expected ')'")),
                }
                Ok(inner)
            }
            None => Err(ParseError::new(col, "unexpected end of input")),
            _ => Err(ParseError::new(col, "unexpected token")),
        }
    }
}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/// Parse an expression string into an [`Expression`] AST.
///
/// Identifiers are resolved via `table.resolve_name`; unknown names produce a
/// [`ParseError`].  Numbers with a unit suffix are converted to SI immediately:
/// `50mm` becomes `Literal(0.05)`.
///
/// # Errors
///
/// Returns [`ParseError`] on any syntax error or unknown parameter reference.
pub fn parse_expression(input: &str, table: &ParameterTable) -> Result<Expression, ParseError> {
    let trimmed = input.trim();
    if trimmed.is_empty() {
        return Err(ParseError::new(0, "empty expression"));
    }
    let tokens = tokenize(trimmed)?;
    if tokens.is_empty() {
        return Err(ParseError::new(0, "empty expression"));
    }
    let input_len = trimmed.len();
    let mut parser = Parser::new(tokens, table, input_len);
    let expr = parser.parse_expr()?;
    if parser.peek().is_some() {
        let col = parser.current_col();
        return Err(ParseError::new(col, "unexpected token after expression"));
    }
    Ok(expr)
}

#[cfg(test)]
#[path = "parse_tests.rs"]
mod tests;
