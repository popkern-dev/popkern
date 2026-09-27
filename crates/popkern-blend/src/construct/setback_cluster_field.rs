//! Affine corner-cluster displacement fields (PK-575) shared by the rim
//! jet (G1) and curvature (G2) reconciles.
//!
//! Where two rim sides end at one non-dual corner cluster, both
//! reconciles used to move the three clamp-band cells of every incident
//! side by ONE shared vector. A translation keeps every seam relation
//! across the fan-internal seam there exact, but it cannot give the two
//! sides different slopes. So under mixed radii it was frozen at zero
//! (hexagon6_mixed kept a -0.69 curvature gap at every such junction).
//!
//! The cluster instead carries a displacement FIELD
//! `δ(x) = c + L·(x − x₀)/h`, evaluated at each band cell's own
//! position. Every seam continuity relation (C^k, or G^k with shape
//! factors) is linear in control points with coefficients summing to
//! zero, and such relations are invariant under an affine map applied
//! to all their participants. So the field is seam-safe for exactly the
//! reason the translation was, with twelve degrees of freedom instead
//! of three. Measured over the apex family (G2 request, whole-seam
//! exact-jet profile): worst rim g2 1.6e-2 → 4.8e-5 (flagship), 0.30 →
//! 8.1e-4 (pentagon5_blunt), 0.70 → 1.8e-4 (hexagon6_mixed), worst rim
//! g1 4.6e-6 → 1.6e-8.

use popkern_core::geometry::gregory::GregoryPatch;
use popkern_core::geometry::{Point3, Vec3};

/// Unknowns per cluster field: the constant `c` and the three columns
/// of `L`, each a `Vec3`.
pub(super) const CLUSTER_FIELD_UNKNOWNS: usize = 4;

/// Parametric-C1 defect below which a fan seam counts as born by
/// subdividing one macro-patch (measured ~1e-13 on such seams, 0.6–0.8
/// on spokes).
pub(super) const PARAMETRIC_SEAM_DEFECT: f64 = 1.0e-6;

/// The corner clusters of one reconcile pass: their fields'
/// normalization and the rim-side ends incident at each.
pub(super) struct ClusterFields {
    centers: Vec<Point3>,
    scale: Vec<f64>,
    sides: Vec<Vec<(usize, usize, bool)>>,
    spoke: Vec<bool>,
}

impl ClusterFields {
    /// Fields over `centers`, unnormalized until cells are included.
    pub(super) fn new(centers: Vec<Point3>) -> Self {
        let n = centers.len();
        ClusterFields {
            centers,
            scale: vec![0.0; n],
            sides: vec![Vec::new(); n],
            spoke: vec![false; n],
        }
    }

    /// Register a band cell of `cluster` at `x`: the field's length
    /// scale `h` is the farthest band cell, so `(x − x₀)/h` stays O(1)
    /// and the normal equations stay well scaled.
    pub(super) fn include(&mut self, cluster: usize, x: Point3) {
        self.scale[cluster] = self.scale[cluster].max(x.distance_to(&self.centers[cluster]));
    }

    /// Register a rim side `(piece, side)` ending at `cluster` (at its
    /// `k = 9` end when `at_k9`).
    pub(super) fn include_side(&mut self, cluster: usize, pi: usize, side: usize, at_k9: bool) {
        self.sides[cluster].push((pi, side, at_k9));
    }

    /// Classify each cluster: only a two-piece junction whose fan seam
    /// is VERIFIED parametric C1 (a subdivision seam) takes the affine
    /// field. Every other cluster keeps a translation: spokes (their
    /// per-side freedom comes from dual-jet claims), one-sided ends, and
    /// junctions whose seam rows do not align control point for control
    /// point. The saddle fan's hybrid-border pieces are that last class:
    /// an affine field there doubled its spoke G1 (7.8e-3 → 1.45e-2).
    pub(super) fn classify(&mut self, patches: &[GregoryPatch]) {
        for c in 0..self.centers.len() {
            let s = &self.sides[c];
            let parametric = s.len() == 2
                && s[0].0 != s[1].0
                && seam_c1_defect(patches, s[0], s[1]).is_some_and(|d| d < PARAMETRIC_SEAM_DEFECT);
            self.spoke[c] = !parametric;
        }
    }

    /// Basis values `[1, (x − x₀)/h]` of the field at `x`; a cluster
    /// not verified parametric keeps a translation (`[1, 0, 0, 0]`). At
    /// a spoke junction the per-side freedom comes from the dual-jet
    /// claims.
    pub(super) fn coords(&self, cluster: usize, x: Point3) -> [f64; CLUSTER_FIELD_UNKNOWNS] {
        if self.spoke[cluster] {
            return [1.0, 0.0, 0.0, 0.0];
        }
        let h = self.scale[cluster].max(1e-12);
        let o = self.centers[cluster];
        [1.0, (x.x - o.x) / h, (x.y - o.y) / h, (x.z - o.z) / h]
    }

    /// The solved field's displacement at `x` (unknowns start at
    /// `CLUSTER_FIELD_UNKNOWNS · cluster`).
    pub(super) fn delta(&self, solution: &[Vec3], cluster: usize, x: Point3) -> Vec3 {
        let q = self.coords(cluster, x);
        (0..CLUSTER_FIELD_UNKNOWNS).fold(Vec3::zero(), |acc, i| {
            acc + solution[CLUSTER_FIELD_UNKNOWNS * cluster + i] * q[i]
        })
    }
}

/// The fan-seam side of a patch at a rim side's end: the patch boundary
/// through that end, perpendicular to the rim side.
pub(super) fn fan_side(rim_side: usize, at_k9: bool) -> usize {
    match (rim_side < 2, at_k9) {
        (true, false) => 2,
        (true, true) => 3,
        (false, false) => 0,
        (false, true) => 1,
    }
}

/// Parametric-C1 defect of the fan seam two incident rim-side ends
/// share: `max_i |c_A,i + c_B,i| / |c_A,i|` over the seam's cross
/// vectors `c = row1 − row0`, with side B's rows aligned to side A's
/// direction. Zero for a seam born by subdividing one macro-patch.
pub(super) fn seam_c1_defect(
    patches: &[GregoryPatch],
    a: (usize, usize, bool),
    b: (usize, usize, bool),
) -> Option<f64> {
    use super::setback_corner_dual_jet::row_for_side;
    let sa = fan_side(a.1, a.2);
    let sb = fan_side(b.1, b.2);
    let a0 = row_for_side(patches, a.0, sa, 0)?;
    let a1 = row_for_side(patches, a.0, sa, 1)?;
    let b0 = row_for_side(patches, b.0, sb, 0)?;
    let b1 = row_for_side(patches, b.0, sb, 1)?;
    let reversed = a0[0].distance_to(&b0[9]) < a0[0].distance_to(&b0[0]);
    let mut worst = 0.0f64;
    for i in 0..10 {
        let j = if reversed { 9 - i } else { i };
        if a0[i].distance_to(&b0[j]) > 1e-6 {
            return None;
        }
        let ca = a1[i] - a0[i];
        let cb = b1[j] - b0[j];
        let len = ca.length().max(1e-300);
        worst = worst.max((ca + cb).length() / len);
    }
    Some(worst)
}

/// Whether a rim side's end is a dual end at a SPOKE junction (the
/// activated corner's other boundary is a fan seam, not a second rim
/// side of the same piece — that would be a setback corner).
pub(super) fn spoke_dual_end(
    rim_sides: &[(usize, usize)],
    pi: usize,
    side: usize,
    at_k9: bool,
) -> bool {
    !rim_sides.contains(&(pi, fan_side(side, at_k9)))
}
