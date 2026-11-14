//! Contact analysis at tangent points: re-exported from popkern-core,
//! where the intersector that needs it lives (PK-1125).
//!
//! The classifier here before subtracted each surface's second
//! fundamental form in its own parameter basis, which compares
//! different bases (at Viviani's crossing it read a touch). Its
//! replacement, [`classify_contact`], expresses both in one frame of the
//! common tangent plane and judges the relative curvature against the
//! explicit tolerance and extent of a [`ContactScale`].

pub use popkern_core::geometry::intersection::contact_classify::{
    ContactClassification, ContactScale, ContactSide, SurfaceContactKind, classify_contact,
};
