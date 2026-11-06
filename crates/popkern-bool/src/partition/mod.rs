//! The regions N solids cut space into (design: operations on the cell
//! complex, section 7; PK-1046).
//!
//! One cell complex over every body: each body's topology in one
//! intersection graph under its own operand index, the records of every
//! body pair whose boxes overlap, the crossings of records of different
//! pairs on a face they share, one walk and one complex. Every region
//! of the complex knows which bodies contain it. A boolean of N bodies,
//! a split or a general fuse is then a question asked of the regions.

mod build;
mod crossings;
mod lamina_part;
mod nary;
mod region_nesting;
mod region_solids;
mod select;
mod set_expr;
mod split;

pub use nary::{FALLBACK_COUNTING_LIMIT, nary_boolean, nary_boolean_where};
pub use region_solids::{Interface, RegionSolid, RegionSolids};
pub use set_expr::SetExpr;
pub use split::{SplitPiece, SplitResult, split};

use popkern_core::geometry::Point3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, BodyKind};

use crate::error::BoolError;

/// The most bodies one partition holds in this stage.
pub const MAX_BODIES: usize = 64;

/// A set of a partition's bodies, by input index.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct BodySet {
    bits: u64,
}

impl BodySet {
    /// The set of the bodies at `indices` (indices at or past
    /// [`MAX_BODIES`] are left out).
    pub fn of(indices: &[usize]) -> Self {
        let mut bits = 0;
        for &i in indices {
            if i < MAX_BODIES {
                bits |= 1 << i;
            }
        }
        Self { bits }
    }

    /// Whether body `index` is in the set.
    pub fn contains(&self, index: usize) -> bool {
        index < MAX_BODIES && self.bits & (1 << index) != 0
    }

    /// How many bodies.
    pub fn len(&self) -> usize {
        self.bits.count_ones() as usize
    }

    /// Whether no body is in the set.
    pub fn is_empty(&self) -> bool {
        self.bits == 0
    }

    /// The bodies, ascending.
    pub fn iter(&self) -> impl Iterator<Item = usize> + '_ {
        (0..MAX_BODIES).filter(|&i| self.contains(i))
    }
}

/// One region: the bodies holding it, and a point strictly inside it.
#[derive(Debug, Clone, PartialEq)]
pub struct RegionInfo {
    /// The bodies containing the region.
    pub bodies: BodySet,
    /// A point strictly inside the region.
    pub interior_point: Point3,
}

/// The regions N solids cut space into, and what is needed to extract
/// them as bodies.
pub struct Partition {
    body_count: usize,
    bodies: Vec<BRep>,
    tolerance: Tolerance,
    built: build::Built,
}

impl std::fmt::Debug for Partition {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Partition")
            .field("body_count", &self.body_count)
            .field("regions", &self.built.regions)
            .finish()
    }
}

impl Partition {
    /// How many bodies the partition was built from.
    pub fn body_count(&self) -> usize {
        self.body_count
    }

    /// Every bounded region (in at least one body).
    pub fn regions(&self) -> &[RegionInfo] {
        &self.built.regions
    }
}

/// Partition `bodies` (solids, and sheets as laminae that cut regions
/// without bounding volume) into the regions they cut space into.
///
/// # Errors
///
/// [`BoolError::TooManyBodies`] past [`MAX_BODIES`];
/// [`BoolError::UnsupportedBodyKind`] for a body that is neither a
/// solid nor a sheet; [`BoolError::DegenerateInput`] for a face with
/// no surface;
/// [`BoolError::MergeFailed`] when the complex cannot be built (its
/// decline named in the reason).
pub fn partition(bodies: &[&BRep], tolerance: Tolerance) -> Result<Partition, BoolError> {
    if bodies.len() > MAX_BODIES {
        return Err(BoolError::TooManyBodies {
            given: bodies.len(),
            limit: MAX_BODIES,
        });
    }
    for (operand, body) in bodies.iter().enumerate() {
        let kind = body.body_kind();
        if kind != BodyKind::Solid && kind != BodyKind::Sheet {
            return Err(BoolError::UnsupportedBodyKind {
                operation: "partition",
                operand,
                kind,
                accepted: &[BodyKind::Solid, BodyKind::Sheet],
            });
        }
        if let Some((face, _)) = body.faces.iter().find(|(_, f)| f.surface.is_none()) {
            return Err(BoolError::DegenerateInput(format!(
                "body {operand} face {face:?} has no surface"
            )));
        }
    }
    // Sheets enter as laminae: their free-boundary mirror loops are no
    // holes (PK-978).
    let (views, lamina_edges): (Vec<BRep>, Vec<_>) = bodies
        .iter()
        .map(|b| {
            if b.body_kind() == BodyKind::Sheet {
                crate::sheet::lamina::pipeline_view(b)
            } else {
                ((*b).clone(), Default::default())
            }
        })
        .unzip();
    let view_refs: Vec<&BRep> = views.iter().collect();
    let built = build::build(&view_refs, &lamina_edges, &tolerance)?;
    let mut partition = Partition {
        body_count: bodies.len(),
        // The bodies as the pipeline read them (a sheet's view): the
        // assembly reads the same topology the graph was built from.
        bodies: views,
        tolerance,
        built,
    };
    partition.nest_regions()?;
    Ok(partition)
}

#[cfg(test)]
#[path = "partition_tests.rs"]
mod tests;

#[cfg(test)]
#[path = "sheet_tests.rs"]
mod sheet_tests;

#[cfg(test)]
#[path = "partition_star_tests.rs"]
mod star_tests;

#[cfg(test)]
#[path = "partition_corpus_tests.rs"]
mod corpus_tests;
