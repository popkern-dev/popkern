//! Edges chained end to end: a sweep path of several edges and a
//! fill's boundary are each given as edge names in any order and
//! direction; the chain orders them from the first, each piece run the
//! way the chain goes.

use popkern_core::geometry::NurbsCurve;
use popkern_core::geometry::curve::{CurveEval, CurveOnEdge};
use popkern_core::geometry::{Point3, Vec3};

/// One edge's curve in the chain, `reversed` when the chain runs
/// against the curve's own direction.
#[derive(Debug, Clone)]
pub(crate) struct ChainPiece {
    /// The edge's curve.
    pub curve: CurveOnEdge,
    /// The chain runs from the curve's end to its start.
    pub reversed: bool,
}

impl ChainPiece {
    /// The piece's start and end in chain order.
    pub(crate) fn ends(&self) -> (Point3, Point3) {
        let (a, b) = (self.at(0.0), self.at(1.0));
        (a, b)
    }

    /// The point a fraction `s` along the piece in chain order.
    pub(crate) fn at(&self, s: f64) -> Point3 {
        self.curve.geometry.point_at(self.param(s))
    }

    /// The unit tangent a fraction `s` along the piece in chain order.
    pub(crate) fn tangent(&self, s: f64) -> Vec3 {
        let t = self.curve.geometry.tangent_at(self.param(s));
        let t = if (self.curve.t_end >= self.curve.t_start) == self.reversed {
            -t
        } else {
            t
        };
        t.normalized().unwrap_or(t)
    }

    fn param(&self, s: f64) -> f64 {
        let s = if self.reversed { 1.0 - s } else { s };
        self.curve.t_start + (self.curve.t_end - self.curve.t_start) * s
    }

    /// The piece as a NURBS curve running in chain order.
    pub(crate) fn nurbs(&self) -> NurbsCurve {
        let mut n = self
            .curve
            .geometry
            .to_nurbs_in_range((self.curve.t_start, self.curve.t_end));
        let backwards = (self.curve.t_end < self.curve.t_start) != self.reversed;
        if backwards {
            n.control_points.reverse();
            n.weights.reverse();
            let (a, b) = (n.knots[0], n.knots[n.knots.len() - 1]);
            n.knots = n.knots.iter().rev().map(|k| a + b - k).collect();
        }
        n
    }
}

/// Why edges do not chain.
#[derive(Debug, Clone, PartialEq)]
pub(crate) enum ChainError {
    /// No edges.
    Empty,
    /// After `chained` pieces no remaining edge meets the chain's end.
    Broken {
        /// Pieces chained before the break.
        chained: usize,
    },
}

/// Chain `curves` end to end from the first, run forward; ends meet
/// within `join`.
pub(crate) fn chain(
    mut curves: Vec<CurveOnEdge>,
    join: f64,
) -> Result<Vec<ChainPiece>, ChainError> {
    if curves.is_empty() {
        return Err(ChainError::Empty);
    }
    let mut out = vec![ChainPiece {
        curve: curves.remove(0),
        reversed: false,
    }];
    while !curves.is_empty() {
        let (_, tail) = out[out.len() - 1].ends();
        let found = curves.iter().enumerate().find_map(|(k, c)| {
            let forward = ChainPiece {
                curve: c.clone(),
                reversed: false,
            };
            let (a, b) = forward.ends();
            if a.distance_to(&tail) <= join {
                Some((k, false))
            } else if b.distance_to(&tail) <= join {
                Some((k, true))
            } else {
                None
            }
        });
        let Some((k, reversed)) = found else {
            return Err(ChainError::Broken { chained: out.len() });
        };
        out.push(ChainPiece {
            curve: curves.remove(k),
            reversed,
        });
    }
    Ok(out)
}
