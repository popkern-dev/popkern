//! Typed indices into the cell complex's arenas.
//!
//! The complex is built in batch and never deletes entities (later
//! stages mark, they do not remove), so dense `u32` indices are stable
//! for its whole life and order-preserving: iteration in index order is
//! deterministic.

macro_rules! complex_id {
    ($(#[$doc:meta])* $name:ident) => {
        $(#[$doc])*
        #[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
        pub(crate) struct $name(pub(crate) u32);

        impl $name {
            /// The arena slot this id names.
            pub(crate) fn index(self) -> usize {
                self.0 as usize
            }

            /// The id of arena slot `i`.
            pub(crate) fn from_index(i: usize) -> Self {
                Self(u32::try_from(i).unwrap_or(u32::MAX))
            }
        }
    };
}

complex_id!(
    /// A vertex of the complex.
    VertexId
);
complex_id!(
    /// An edge of the complex (a 1-cell; any number of fins).
    EdgeId
);
complex_id!(
    /// One face's use of one edge, inside one loop of that face.
    FinId
);
complex_id!(
    /// A boundary loop of a complex face.
    LoopId
);
complex_id!(
    /// A complex face (a 2-cell with two sides).
    FaceId
);
complex_id!(
    /// A connected set of face sides bounding one region locally.
    ShellId
);
