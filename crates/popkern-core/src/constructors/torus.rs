//! Torus constructor (stub — deferred to v1b).

use crate::error::CoreError;
use crate::tolerance::Tolerance;
use crate::topology::BRep;

impl BRep {
    /// Build a torus — deferred to v1b.
    ///
    /// The torus CW-complex topology (V=1, E=2, F=1) requires edge
    /// identification operations that are not yet implemented. This
    /// constructor returns an error until those operations are available.
    ///
    /// # Errors
    ///
    /// Always returns [`CoreError::PreconditionFailed`].
    pub fn torus(
        _major_radius: f64,
        _minor_radius: f64,
        _tolerance: Tolerance,
    ) -> Result<Self, CoreError> {
        Err(CoreError::PreconditionFailed(
            "torus constructor deferred to v1b — requires edge identification operations".into(),
        ))
    }
}
