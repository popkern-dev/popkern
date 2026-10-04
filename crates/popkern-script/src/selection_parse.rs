//! The selection grammar (spec §6.1) as text: parse a canonical
//! `describe()` spelling back into a selector, and resolve a selector
//! against an arbitrary body. Skins that carry queries as text (MCP
//! tool arguments, the canonical emitter's verification pass) share
//! this one reader.

use popkern_core::topology::BRep;
use popkern_part::{Dir, SurfaceKind};

use crate::error::ScriptError;
use crate::model::Model;
use crate::selection::{EdgeSelector, EntityRef, FaceSelector, SelectionSnapshot};

/// Why a selection text did not parse.
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
#[error("selection text '{text}': {reason}")]
pub struct SelectionParseError {
    /// The text.
    pub text: String,
    /// What was wrong.
    pub reason: String,
}

struct Cursor<'a> {
    text: &'a str,
    pos: usize,
}

impl<'a> Cursor<'a> {
    fn fail(&self, reason: impl Into<String>) -> SelectionParseError {
        SelectionParseError {
            text: self.text.to_string(),
            reason: format!("{} at offset {}", reason.into(), self.pos),
        }
    }

    fn peek(&self) -> Option<char> {
        self.text[self.pos..].chars().next()
    }

    fn eat(&mut self, literal: &str) -> bool {
        if self.text[self.pos..].starts_with(literal) {
            self.pos += literal.len();
            true
        } else {
            false
        }
    }

    fn require(&mut self, literal: &str) -> Result<(), SelectionParseError> {
        if self.eat(literal) {
            Ok(())
        } else {
            Err(self.fail(format!("expected '{literal}'")))
        }
    }

    fn skip_ws(&mut self) {
        while self.peek().is_some_and(char::is_whitespace) {
            self.pos += 1;
        }
    }

    fn identifier(&mut self) -> Result<&'a str, SelectionParseError> {
        let start = self.pos;
        while self
            .peek()
            .is_some_and(|c| c.is_ascii_alphanumeric() || c == '_')
        {
            self.pos += 1;
        }
        if start == self.pos {
            return Err(self.fail("expected an identifier"));
        }
        Ok(&self.text[start..self.pos])
    }

    /// One argument: a quoted string, a number, or a bare word.
    fn argument(&mut self) -> Result<Arg, SelectionParseError> {
        self.skip_ws();
        if self.eat("\"") {
            let mut out = String::new();
            loop {
                match self.peek() {
                    None => return Err(self.fail("unterminated string")),
                    Some('"') => {
                        self.pos += 1;
                        return Ok(Arg::Text(out));
                    }
                    Some('\\') => {
                        self.pos += 1;
                        match self.peek() {
                            Some(c) => {
                                out.push(c);
                                self.pos += c.len_utf8();
                            }
                            None => return Err(self.fail("unterminated escape")),
                        }
                    }
                    Some(c) => {
                        out.push(c);
                        self.pos += c.len_utf8();
                    }
                }
            }
        }
        let start = self.pos;
        while self
            .peek()
            .is_some_and(|c| c.is_ascii_alphanumeric() || matches!(c, '_' | '.' | '-' | '+'))
        {
            self.pos += 1;
        }
        let word = &self.text[start..self.pos];
        if word.is_empty() {
            return Err(self.fail("expected an argument"));
        }
        if let Ok(v) = word.parse::<f64>() {
            return Ok(Arg::Number(v));
        }
        Ok(Arg::Word(word.to_string()))
    }

    fn arguments(&mut self) -> Result<Vec<Arg>, SelectionParseError> {
        self.require("(")?;
        let mut args = Vec::new();
        self.skip_ws();
        if self.eat(")") {
            return Ok(args);
        }
        loop {
            args.push(self.argument()?);
            self.skip_ws();
            if self.eat(")") {
                return Ok(args);
            }
            self.require(",")?;
        }
    }
}

enum Arg {
    Text(String),
    Number(f64),
    Word(String),
}

fn surface_kind(c: &Cursor<'_>, arg: &Arg) -> Result<SurfaceKind, SelectionParseError> {
    let word = match arg {
        Arg::Word(w) | Arg::Text(w) => w.as_str(),
        Arg::Number(_) => return Err(c.fail("expected a surface kind")),
    };
    Ok(match word {
        "plane" => SurfaceKind::Plane,
        "cylinder" => SurfaceKind::Cylinder,
        "cone" => SurfaceKind::Cone,
        "sphere" => SurfaceKind::Sphere,
        "torus" => SurfaceKind::Torus,
        "nurbs" => SurfaceKind::Nurbs,
        "gregory" => SurfaceKind::Gregory,
        other => return Err(c.fail(format!("unknown surface kind '{other}'"))),
    })
}

fn dir(c: &Cursor<'_>, arg: &Arg) -> Result<Dir, SelectionParseError> {
    let word = match arg {
        Arg::Word(w) | Arg::Text(w) => w.as_str(),
        Arg::Number(_) => return Err(c.fail("expected a direction")),
    };
    Ok(match word {
        "x" => Dir::X,
        "y" => Dir::Y,
        "z" => Dir::Z,
        other => return Err(c.fail(format!("unknown direction '{other}'"))),
    })
}

fn number(c: &Cursor<'_>, arg: &Arg) -> Result<f64, SelectionParseError> {
    match arg {
        Arg::Number(v) => Ok(*v),
        _ => Err(c.fail("expected a number")),
    }
}

fn text(c: &Cursor<'_>, arg: &Arg) -> Result<String, SelectionParseError> {
    match arg {
        Arg::Text(t) => Ok(t.clone()),
        _ => Err(c.fail("expected a quoted string")),
    }
}

fn arity(c: &Cursor<'_>, name: &str, args: &[Arg], n: usize) -> Result<(), SelectionParseError> {
    if args.len() == n {
        Ok(())
    } else {
        Err(c.fail(format!("{name} takes {n} argument(s), got {}", args.len())))
    }
}

fn edge_query(c: &mut Cursor<'_>) -> Result<EdgeSelector, SelectionParseError> {
    c.skip_ws();
    c.require("edges()")?;
    let mut sel = crate::selection::edges();
    loop {
        c.skip_ws();
        if !c.eat(".") {
            return Ok(sel);
        }
        let name = c.identifier()?;
        if name == "or" {
            c.require("(")?;
            let other = edge_query(c)?;
            c.skip_ws();
            c.require(")")?;
            sel = sel.or(other);
            continue;
        }
        let args = c.arguments()?;
        sel = match name {
            "circular" => sel.circular(),
            "straight" => sel.straight(),
            "conic" => sel.conic(),
            "nurbs" => sel.nurbs(),
            "not_seam" => sel.not_seam(),
            "allow_empty" => sel.allow_empty(),
            "between" => {
                arity(c, name, &args, 2)?;
                sel.between(surface_kind(c, &args[0])?, surface_kind(c, &args[1])?)
            }
            "at_height" => {
                arity(c, name, &args, 1)?;
                sel.at_height(number(c, &args[0])?)
            }
            "direction" => {
                arity(c, name, &args, 1)?;
                sel.direction(dir(c, &args[0])?)
            }
            "chord_length_between" => {
                arity(c, name, &args, 2)?;
                sel.chord_length_between(number(c, &args[0])?, number(c, &args[1])?)
            }
            "labeled" => {
                arity(c, name, &args, 1)?;
                sel.labeled(&text(c, &args[0])?)
            }
            "named" => {
                arity(c, name, &args, 1)?;
                sel.spelled(&text(c, &args[0])?)
            }
            other => return Err(c.fail(format!("unknown edge predicate '{other}'"))),
        };
    }
}

fn face_query(c: &mut Cursor<'_>) -> Result<FaceSelector, SelectionParseError> {
    c.skip_ws();
    c.require("faces()")?;
    let mut sel = crate::selection::faces();
    loop {
        c.skip_ws();
        if !c.eat(".") {
            return Ok(sel);
        }
        let name = c.identifier()?;
        let args = c.arguments()?;
        sel = match name {
            "allow_empty" => sel.allow_empty(),
            "kind" => {
                arity(c, name, &args, 1)?;
                sel.kind(surface_kind(c, &args[0])?)
            }
            "facing" => {
                arity(c, name, &args, 1)?;
                sel.facing(dir(c, &args[0])?)
            }
            "facing_neg" => {
                arity(c, name, &args, 1)?;
                sel.facing_neg(dir(c, &args[0])?)
            }
            "at_height" => {
                arity(c, name, &args, 1)?;
                sel.at_height(number(c, &args[0])?)
            }
            "labeled" => {
                arity(c, name, &args, 1)?;
                sel.labeled(&text(c, &args[0])?)
            }
            "named" => {
                arity(c, name, &args, 1)?;
                sel.spelled(&text(c, &args[0])?)
            }
            other => return Err(c.fail(format!("unknown face predicate '{other}'"))),
        };
    }
}

fn finish<T>(c: &Cursor<'_>, value: T) -> Result<T, SelectionParseError> {
    if c.pos == c.text.len() {
        Ok(value)
    } else {
        Err(c.fail("trailing text"))
    }
}

/// Parse an edge query in its canonical spelling.
pub fn parse_edges(text: &str) -> Result<EdgeSelector, SelectionParseError> {
    let mut c = Cursor { text, pos: 0 };
    let sel = edge_query(&mut c)?;
    c.skip_ws();
    finish(&c, sel)
}

/// Parse a face query in its canonical spelling.
pub fn parse_faces(text: &str) -> Result<FaceSelector, SelectionParseError> {
    let mut c = Cursor { text, pos: 0 };
    let sel = face_query(&mut c)?;
    c.skip_ws();
    finish(&c, sel)
}

impl Model {
    /// Resolve an edge query against an explicit body (labels and
    /// spelled names resolve through this document).
    pub fn resolve_edges_in(
        &self,
        brep: &BRep,
        selector: &EdgeSelector,
    ) -> Result<SelectionSnapshot, ScriptError> {
        let resolve = |r: EntityRef| self.resolve_entity_ref(r);
        crate::selection::resolve_edges(brep, selector, &resolve)
    }

    /// Resolve a face query against an explicit body.
    pub fn resolve_faces_in(
        &self,
        brep: &BRep,
        selector: &FaceSelector,
    ) -> Result<SelectionSnapshot, ScriptError> {
        let resolve = |r: EntityRef| self.resolve_entity_ref(r);
        crate::selection::resolve_faces(brep, selector, &resolve)
    }

    /// `body` as it stood before the feature at `index` ran, read from
    /// the step cache (steps the cache does not hold execute without
    /// touching the model's current bodies or report).
    pub fn body_before(
        &mut self,
        index: usize,
        body: popkern_history::node::BodyId,
    ) -> Result<Option<BRep>, ScriptError> {
        let state = self.tree.state_at(index, self.tolerance)?;
        Ok(state.body(body).cloned())
    }
}
