//! Shirman-Séquin vertex-enclosure solve — promoted to
//! `popkern_core::geometry::gregory::vertex_enclosure` (PK-248: the
//! core hole-healing polar fan needs it; the math is pure geometry).
//! This shim keeps the fillet-side paths stable.

#[allow(unused_imports)]
pub(super) use popkern_core::geometry::gregory::vertex_enclosure::{
    EnclosureSolve, shirman_sequin_vertex_enclosure_solve,
};
