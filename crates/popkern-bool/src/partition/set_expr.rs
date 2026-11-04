//! Which regions of a partition a request keeps: a set expression over
//! the bodies, by input index (design 5.1, 6). Data, so a feature tree
//! records and replays it; a region id would change with any edit.

use serde::{Deserialize, Serialize};

use super::BodySet;

/// A set expression over a partition's bodies.
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum SetExpr {
    /// Inside body `i`.
    Body(usize),
    /// Inside any body (∪).
    Any,
    /// Inside every body (∩).
    All,
    /// Inside exactly one body.
    ExactlyOne,
    /// Inside an odd number of bodies (the N-ary ⊕; equal to
    /// `ExactlyOne` for two bodies).
    Odd,
    /// Inside at least this many bodies.
    AtLeast(usize),
    /// Both.
    And(Box<SetExpr>, Box<SetExpr>),
    /// Either.
    Or(Box<SetExpr>, Box<SetExpr>),
    /// The first and not the second.
    Minus(Box<SetExpr>, Box<SetExpr>),
}

impl SetExpr {
    /// Inside body `i`.
    pub fn body(i: usize) -> Self {
        SetExpr::Body(i)
    }

    /// Inside both `self` and `other`.
    pub fn and(self, other: SetExpr) -> Self {
        SetExpr::And(Box::new(self), Box::new(other))
    }

    /// Inside `self` or `other`.
    pub fn or(self, other: SetExpr) -> Self {
        SetExpr::Or(Box::new(self), Box::new(other))
    }

    /// Inside `self` and not `other`.
    pub fn minus(self, other: SetExpr) -> Self {
        SetExpr::Minus(Box::new(self), Box::new(other))
    }

    /// Whether a region inside exactly `bodies` (of `count`) is kept.
    pub fn holds(&self, bodies: &BodySet, count: usize) -> bool {
        match self {
            SetExpr::Body(i) => bodies.contains(*i),
            SetExpr::Any => !bodies.is_empty(),
            SetExpr::All => count > 0 && bodies.len() == count,
            SetExpr::ExactlyOne => bodies.len() == 1,
            SetExpr::Odd => bodies.len() % 2 == 1,
            SetExpr::AtLeast(k) => bodies.len() >= *k,
            SetExpr::And(a, b) => a.holds(bodies, count) && b.holds(bodies, count),
            SetExpr::Or(a, b) => a.holds(bodies, count) || b.holds(bodies, count),
            SetExpr::Minus(a, b) => a.holds(bodies, count) && !b.holds(bodies, count),
        }
    }

    /// The highest body index the expression names, if any.
    pub fn max_body(&self) -> Option<usize> {
        match self {
            SetExpr::Body(i) => Some(*i),
            SetExpr::And(a, b) | SetExpr::Or(a, b) | SetExpr::Minus(a, b) => {
                a.max_body().max(b.max_body())
            }
            _ => None,
        }
    }
}
