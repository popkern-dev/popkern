//! Cap kind decision for the sphere-fit edge-wiring workirectApex apex construction.
//!
//! `CapKind` enumerates how the apex region's cap face is built.
//! For trihedral symmetric corners (uniform radii, coincident T's),
//! a single Sphere face suffices.  For variable-radius corners and
//! higher-valence apices, an N-sided Gregory cap is constructed.

use popkern_core::geometry::Point3;
use serde::{Deserialize, Serialize};

/// Kind of the corner between adjacent fillet terminations on a
/// repeated apex (corner `k` sits between fillet `k`'s T_b and
/// fillet `k+1`'s T_a).
///
/// On a UNIFORM-radius apex the two rolling-ball contact curves on
/// each support face intersect at one point, so adjacent
/// terminations coincide and the FaceEdge bridge between them has
/// no length — the corner is a single welded vertex (the
/// shared-corner class, PK-216/PK-219). Mixed radii keep the
/// corners apart and the bridge is a real edge.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub enum ApexCornerKind {
    /// Non-degenerate FaceEdge bridge between the adjacent
    /// terminations (the historical 2N Gregory contract).
    Bridged,
    /// The adjacent terminations coincide within the weld gate:
    /// the corner is ONE vertex at `point` (the settled corner
    /// position, authoritative for both adjacent fillet-end arcs)
    /// and no FaceEdge exists.
    Welded {
        /// Settled corner position. Both adjacent fillet-end arc
        /// endpoints are snapped to this exact point so the cap
        /// boundary closes bitwise at the corner.
        point: Point3,
    },
}

impl ApexCornerKind {
    /// True for the welded (shared-corner) kind.
    pub fn is_welded(&self) -> bool {
        matches!(self, ApexCornerKind::Welded { .. })
    }
}

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
        /// For a corner `k` classified [`ApexCornerKind::Welded`],
        /// `per_fillet_terminations[k][1]` and
        /// `per_fillet_terminations[(k+1)%N][0]` are BOTH the weld
        /// point (bitwise equal).
        per_fillet_terminations: Vec<[Point3; 2]>,
        /// Per-fillet (sphere centre, sphere radius).  May vary
        /// across fillets.
        per_fillet_caps: Vec<(Point3, f64)>,
        /// Corner classification, length N: `corner_kinds[k]` is the
        /// corner between fillet `k` and fillet `k+1` (mod N).
        /// Empty (older serialized plans) means all corners
        /// Bridged — the pre-shared-corner contract.
        #[serde(default)]
        corner_kinds: Vec<ApexCornerKind>,
    },

    /// No cap is constructed — the planner's [Open
    /// rung](super::telemetry::ConstructionRung::Open) when both
    /// Direct and Setback failed.  Fillets terminate at their
    /// rolling-ball footprints; the apex region is left as
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

    /// Corner kinds for a Gregory cap (empty slice for Sphere/Open
    /// and for plans that predate corner classification —
    /// both mean "no welded corners").
    pub fn corner_kinds(&self) -> &[ApexCornerKind] {
        match self {
            CapKind::Gregory { corner_kinds, .. } => corner_kinds,
            _ => &[],
        }
    }

    /// True when corner `k` (between fillet `k` and `k+1`) is
    /// welded. Out-of-range or unclassified corners read Bridged.
    pub fn corner_is_welded(&self, k: usize) -> bool {
        self.corner_kinds()
            .get(k)
            .is_some_and(ApexCornerKind::is_welded)
    }

    /// Number of welded corners (0 for Sphere/Open/pre-classification plans).
    pub fn welded_corner_count(&self) -> usize {
        self.corner_kinds()
            .iter()
            .filter(|kind| kind.is_welded())
            .count()
    }

    /// True if no cap is constructed (Open rung outcome).
    pub fn is_open(&self) -> bool {
        matches!(self, CapKind::Open { .. })
    }
}
