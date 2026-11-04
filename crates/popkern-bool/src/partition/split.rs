//! A body cut into its pieces by tools, every piece kept (design 6).

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use super::{BodySet, Interface, SetExpr, partition};
use crate::error::BoolError;

/// One piece of a split.
#[derive(Debug, Clone)]
pub struct SplitPiece {
    /// The piece, a solid.
    pub body: BRep,
    /// The tools (by index into the `tools` given) the piece lies inside.
    pub inside_tools: BodySet,
}

/// The pieces of a split and the faces neighbouring pieces share.
#[derive(Debug, Clone)]
pub struct SplitResult {
    /// Every piece of the target.
    pub pieces: Vec<SplitPiece>,
    /// Faces two pieces share (indices into `pieces`).
    pub interfaces: Vec<Interface>,
}

/// `target` cut by `tools` into its pieces: every connected region of
/// the target between the tools' boundaries, each its own solid, with
/// the tools it lies inside. A tool that misses the target, or only
/// touches it, leaves one piece: the target.
///
/// # Errors
///
/// The partition's and extraction's errors (see [`partition`]).
pub fn split(
    target: &BRep,
    tools: &[&BRep],
    tolerance: Tolerance,
) -> Result<SplitResult, BoolError> {
    let mut bodies = vec![target];
    bodies.extend_from_slice(tools);
    let p = partition(&bodies, tolerance)?;
    let cells = p.region_solids(&SetExpr::body(0))?;
    let pieces = cells
        .solids
        .into_iter()
        .map(|s| {
            let inside: Vec<usize> = s.bodies.iter().filter(|&i| i > 0).map(|i| i - 1).collect();
            SplitPiece {
                body: s.body,
                inside_tools: BodySet::of(&inside),
            }
        })
        .collect();
    Ok(SplitResult {
        pieces,
        interfaces: cells.interfaces,
    })
}
