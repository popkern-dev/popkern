//! The chart service (F4, design `docs/architecture/chart-service.md`,
//! PK-1228): a surface's chart artifacts (the period cut and the singular
//! rows) described once, and charts free of them around any point. The
//! only place in the workspace that handles chart artifacts; consumers
//! ask here (the `structural_chart_service` ratchet counts the rest).
//!
//! Three meanings of "seam" stay apart: chart artifacts (this module),
//! topological seam edges (topology, asking here for chart questions) and
//! continuity seams between faces (blend certification, unrelated).

pub mod charted;
pub mod descriptor;
pub mod ops;
pub mod periods;
pub mod regular;
pub mod singularity;
pub mod wrap;

pub use charted::ChartedSurface;
pub use descriptor::{
    Axis, ChartDescriptor, ChartPeriod, Continuity, PeriodKind, SingularKind, SingularRow,
};
/// Named period arithmetic (`wrap_from`, `offset_mod`, `modulo`).
pub use popkern_math::periodic;
pub use regular::{ChartError, RegularChart, regular_chart};

#[cfg(test)]
mod tests;
