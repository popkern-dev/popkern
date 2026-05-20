//! Cap kind decision for Phase 6c-DirectApex apex construction.
//!
//! `CapKind` enumerates how the apex region's cap face is built.
//! For trihedral symmetric corners (uniform radii, coincident T's),
//! a single Sphere face suffices.  For variable-radius corners and
//! higher-valence apices, an N-sided Gregory cap is constructed.

use popkern_core::geometry::Point3;
use serde::{Deserialize, Serialize};

/// Cap face construction strategy decided by [`plan_apex`](super::plan_apex).
///
/// The choice between `Sphere` and `Gregory` is driven by whether
/// adjacent FaceEdge bridges between fillet terminations collapse
/// to (effectively) zero length and the apex is trihedral — which
/// indicates the corner is symmetric enough that a single
/// tangent-circle Sphere covers all fillets.  Variable-radius and
/// higher-valence corners need an N-sided Gregory cap.
///
/// See ADR §3.1 for full type design.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[non_exhaustive]
pub enum CapKind {
    /// Single sphere cap, tangent to each fillet's cylinder along a
    /// tangent circle.  Used for trihedral corners when adjacent
    /// fillet terminations collapse (uniform radii / symmetric
    /// corners).
    ///
    /// All fillets terminate at points on this sphere, on a single
    /// shared tangent circle.  The cap face has N edges (one per
    /// adjacent fillet pair).
    Sphere {
        /// Cap-sphere centre.
        center: Point3,
        /// Cap-sphere radius.  Equal to the (uniform) fillet radius.
        radius: f64,
        /// Endpoints of the tangent arc on each fillet's cylinder,
        /// in cap-loop walk order.  Length = number of fillets.
        tangent_arc_endpoints: Vec<Point3>,
    },

    /// N-sided Gregory cap with distinct per-fillet sphere centres
    /// and non-degenerate FaceEdge bridges.  Used when fillet radii
    /// vary across the corner.
    ///
    /// Each fillet has its own sphere (cap-tangent circle on its
    /// cylinder), distinct from neighbours'.  The cap face has 2N
    /// edges: N FilletEnd arcs (cap-vs-Cyl) and N FaceEdge bridges
    /// (cap-vs-Plane), alternating.
    Gregory {
        /// Per-fillet sphere-fit terminations [T_a, T_b], one per
        /// fillet.  T_a is on face_A side, T_b on face_B side.
        per_fillet_terminations: Vec<[Point3; 2]>,
        /// Per-fillet (sphere centre, sphere radius).  May vary
        /// across fillets.
        per_fillet_caps: Vec<(Point3, f64)>,
    },

    /// No cap is constructed — the planner's [Open
    /// rung](super::telemetry::ConstructionRung::Open) when both
    /// Direct and Setback failed.  Fillets terminate at their
    /// legacy rolling-ball footprints; the apex region is left as
    /// the pre-DirectApex substitute-then-patch overlay would have
    /// produced (or unfilleted, depending on caller).
    ///
    /// Carries the count of fillets that *would* have been at the
    /// apex, for telemetry / UI feedback only.
    Open {
        /// Number of fillets that meet at the apex (purely
        /// informational — no cap geometry exists).
        n_fillets_attempted: usize,
    },
}

impl CapKind {
    /// Number of fillets meeting at the apex (= N).  For
    /// [`CapKind::Open`] this is the count of fillets that *would*
    /// have been capped.
    pub fn n_fillets(&self) -> usize {
        match self {
            CapKind::Sphere {
                tangent_arc_endpoints,
                ..
            } => tangent_arc_endpoints.len(),
            CapKind::Gregory {
                per_fillet_terminations,
                ..
            } => per_fillet_terminations.len(),
            CapKind::Open {
                n_fillets_attempted,
            } => *n_fillets_attempted,
        }
    }

    /// True if this cap is the simpler Sphere kind.
    pub fn is_sphere(&self) -> bool {
        matches!(self, CapKind::Sphere { .. })
    }

    /// True if this cap is the N-sided Gregory kind.
    pub fn is_gregory(&self) -> bool {
        matches!(self, CapKind::Gregory { .. })
    }

    /// True if no cap is constructed (Open rung outcome).
    pub fn is_open(&self) -> bool {
        matches!(self, CapKind::Open { .. })
    }
}
