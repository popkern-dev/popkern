//! API types for fillet and chamfer specifications.

use popkern_core::topology::{EdgeHandle, FaceHandle};

use crate::error::BlendWarning;

/// Result of a successful blend operation.
#[derive(Debug, Clone)]
pub struct BlendResult {
    /// The modified BRep with blends applied.
    pub brep: popkern_core::topology::BRep,
    /// Per-edge warnings emitted during blending.
    pub edge_warnings: Vec<(EdgeHandle, BlendWarning)>,
}

/// Specification for a single fillet edge.
#[derive(Debug, Clone)]
pub struct FilletSpec {
    /// The edge to fillet.
    pub edge: EdgeHandle,
    /// The fillet radius profile.
    pub radius: FilletRadius,
}

/// Fillet radius profile along an edge.
#[derive(Debug, Clone)]
pub enum FilletRadius {
    /// Constant radius along the entire edge.
    Constant(f64),
    /// Linearly varying radius from start to end of the edge.
    Variable {
        /// Radius at the start vertex.
        start: f64,
        /// Radius at the end vertex.
        end: f64,
    },
    /// Radius defined at specific parameter positions along the edge.
    ByPoints(Vec<(f64, f64)>),
}

impl From<f64> for FilletRadius {
    /// Convenience: a bare `f64` becomes a constant-radius profile.
    /// Lets `modify_fillet_radius(brep, face, 0.2, tol)` and
    /// `modify_fillet_radius(brep, face, FilletRadius::Variable { ... }, tol)`
    /// share the same entry point.
    fn from(r: f64) -> Self {
        FilletRadius::Constant(r)
    }
}

/// Options controlling fillet behaviour.
#[derive(Debug, Clone)]
pub struct FilletOptions {
    /// Whether to automatically propagate fillets along tangent-continuous edge chains.
    pub propagate_tangent: bool,
    /// How to handle vertex blends where three or more fillets meet.
    pub vertex_mode: VertexMode,
    /// Build the blend surface with **G2** (curvature-matching)
    /// cross-sections rather than the default G1 (circular-arc)
    /// sweep. Eliminates highlight-line breaks at the blend
    /// boundaries — required for §15 Class-A surfaces.
    ///
    /// **Default `false`.** Setting `true` routes the blend through
    /// the spine-walker contact solver and skins per-sample
    /// [`contact::ContactSolution`](crate::contact::ContactSolution)s
    /// into a degree-5 Bezier sweep via
    /// [`g2_cross_section::g2_sweep_from_contact_solutions`](crate::g2_cross_section::g2_sweep_from_contact_solutions).
    /// The G2 path applies to **any supported surface pair**
    /// (plane-plane, plane-cylinder, sphere-plane, cone-plane, …),
    /// not just cones.
    ///
    /// **Graceful degradation.** If the spine walker fails on a
    /// non-cone pair (e.g. degenerate normals, periodicity edge
    /// case), the routing silently falls back to the legacy G1
    /// arc sweep rather than failing the whole fillet operation.
    /// Callers who want G2 to be strict should audit the returned
    /// surface type (a Torus/analytic result for plane-plane
    /// indicates the fallback fired).
    ///
    /// **Cost.** G2 is materially more expensive than G1 — the
    /// walker performs ~20 contact-solver iterations per spine
    /// sample, where G1 does a single analytic offset intersect.
    /// Leave this `false` for interactive preview and enable it
    /// only for final / Class-A output.
    ///
    /// See §2.4 / §7 G2 in `docs/requirements/fillets_blends.md`.
    ///
    /// # Example
    ///
    /// ```no_run
    /// use popkern_fillet::{FilletOptions, FilletSpec, FilletRadius, fillet_edges};
    /// # use popkern_core::topology::BRep;
    /// # use popkern_core::tolerance::Tolerance;
    /// # let brep: BRep = unimplemented!();
    /// # let edge = unimplemented!();
    /// let specs = vec![FilletSpec {
    ///     edge,
    ///     radius: FilletRadius::Constant(0.1),
    /// }];
    /// let opts = FilletOptions {
    ///     g2_cross_section: true,
    ///     ..FilletOptions::default()
    /// };
    /// let result = fillet_edges(&brep, &specs, &opts, Tolerance::default(), None);
    /// ```
    pub g2_cross_section: bool,
}

impl Default for FilletOptions {
    fn default() -> Self {
        Self {
            propagate_tangent: true,
            vertex_mode: VertexMode::Smooth,
            g2_cross_section: false,
        }
    }
}

/// Vertex blend strategy at multi-fillet junctions.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum VertexMode {
    /// Smooth G1-continuous vertex blend.
    Smooth,
    /// Mitred (flat) vertex blend.
    Mitred,
}

/// Specification for a single chamfer edge.
#[derive(Debug, Clone)]
pub struct ChamferSpec {
    /// The edge to chamfer.
    pub edge: EdgeHandle,
    /// The chamfer geometry.
    pub kind: ChamferKind,
}

/// Chamfer geometry definition.
#[derive(Debug, Clone)]
pub enum ChamferKind {
    /// Equal-distance chamfer (45-degree bevel).
    Equal(f64),
    /// Two-distance chamfer with specified face orientation.
    TwoDistance {
        /// Reference face for the first distance.
        face_1: FaceHandle,
        /// Distance on the first face.
        dist_1: f64,
        /// Distance on the second (opposite) face.
        dist_2: f64,
    },
    /// Distance-and-angle chamfer.
    DistanceAngle {
        /// Reference face for the distance measurement.
        face: FaceHandle,
        /// Distance from the edge on the reference face.
        distance: f64,
        /// Angle in radians from the reference face.
        angle: f64,
    },
}

/// Options controlling chamfer behaviour.
#[derive(Debug, Clone)]
pub struct ChamferOptions {
    /// Whether to automatically propagate chamfers along tangent-continuous edge chains.
    pub propagate_tangent: bool,
}

impl Default for ChamferOptions {
    fn default() -> Self {
        Self {
            propagate_tangent: true,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fillet_options_default() {
        let opts = FilletOptions::default();
        assert!(opts.propagate_tangent);
        assert_eq!(opts.vertex_mode, VertexMode::Smooth);
    }

    #[test]
    fn chamfer_options_default() {
        let opts = ChamferOptions::default();
        assert!(opts.propagate_tangent);
    }

    #[test]
    fn fillet_radius_constant() {
        let r = FilletRadius::Constant(2.5);
        match r {
            FilletRadius::Constant(v) => assert!((v - 2.5).abs() < f64::EPSILON),
            _ => panic!("expected Constant"),
        }
    }
}
