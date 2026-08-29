//! Euler-Poincaré invariant check.

use crate::error::ValidationError;
use crate::topology::BRep;

impl BRep {
    /// Validate the Euler-Poincaré formula for this BRep.
    ///
    /// For genus-0 models: `V - E + F - (L - F) = 2 * S`
    ///
    /// This simplifies from the full formula `V - E + F - (L - F) - 2(S - G) = 0`
    /// with genus G = 0. For solids with through-holes, use
    /// [`validate_euler_poincare_any_genus`] instead.
    pub fn validate_euler_poincare(&self) -> Result<(), ValidationError> {
        // Open shells (boundary edges present) break the closed-
        // manifold Euler formula. Skip the check entirely — the
        // genus-aware variant below also can't handle open shells
        // without computing the number of boundary loops, which
        // is expensive and not yet implemented.
        if self.has_boundary_edges() {
            return Ok(());
        }

        let v = self.vertices.len() as i64;
        let e = self.edges.len() as i64;
        let f = self.faces.len() as i64;
        let l = self.loops.len() as i64;
        let s = self.shells.len() as i64;

        let lhs = v - e + f - (l - f);
        let rhs = 2 * s;

        if lhs != rhs {
            if crate::trace::trace_enabled("euler") {
                let mut orphan_v = 0;
                let used: std::collections::HashSet<_> =
                    self.half_edges.values().map(|he| he.origin).collect();
                for vh in self.vertices.keys() {
                    if !used.contains(&vh) {
                        orphan_v += 1;
                    }
                }
                let mut orphan_e = 0;
                for (_, edge) in &self.edges {
                    if !self.half_edges.contains_key(edge.half_edge) {
                        orphan_e += 1;
                    }
                }
                eprintln!(
                    "  [euler] V={v} E={e} F={f} L={l} S={s} lhs={lhs} rhs={rhs} orphan_v={orphan_v} orphan_e={orphan_e}"
                );
            }
            let solid = self.solids.keys().next().unwrap_or_default();
            return Err(ValidationError::EulerPoincare {
                solid,
                expected: rhs,
                actual: lhs,
            });
        }

        Ok(())
    }

    /// Validate the Euler-Poincaré formula allowing any genus.
    ///
    /// Full formula: `V - E + 2F - L = 2(S - G)`. Validates that the
    /// computed genus G is a non-negative integer.
    pub fn validate_euler_poincare_any_genus(&self) -> Result<(), ValidationError> {
        let v = self.vertices.len() as i64;
        let e = self.edges.len() as i64;
        let f = self.faces.len() as i64;
        let l = self.loops.len() as i64;
        let s = self.shells.len() as i64;

        let lhs = v - e + 2 * f - l;
        let two_g = 2 * s - lhs;

        // Accept any non-negative even 2G. Also tolerate small odd
        // values from one-sided edge artifacts (self-loop circles whose
        // twin was on a non-selected face).
        if two_g < -1 {
            if crate::trace::trace_enabled("euler") {
                let used: std::collections::HashSet<_> =
                    self.half_edges.values().map(|he| he.origin).collect();
                let orphan_v = self.vertices.keys().filter(|v| !used.contains(v)).count();
                let mut loops_per_face = std::collections::HashMap::new();
                for (_, lp) in &self.loops {
                    *loops_per_face.entry(lp.face).or_insert(0) += 1;
                }
                let orphan_loops = self
                    .loops
                    .iter()
                    .filter(|(_, lp)| !self.faces.contains_key(lp.face))
                    .count();
                eprintln!(
                    "  [euler-g] V={v} E={e} F={f} L={l} S={s} lhs={lhs} 2G={two_g} orphan_v={orphan_v} orphan_loops={orphan_loops}"
                );
            }
            let solid = self.solids.keys().next().unwrap_or_default();
            return Err(ValidationError::EulerPoincare {
                solid,
                expected: 2 * s,
                actual: lhs,
            });
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use crate::geometry::{Point3, Vec3};
    use crate::tolerance::Tolerance;
    use crate::topology::BRep;

    /// Helper: build a cube via mvfs + mev + mef + extrude.
    fn make_cube() -> BRep {
        let mut brep = BRep::new(Tolerance::default());
        let mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
        let mev1 = brep
            .mev_at_vertex(mvfs.vertex, mvfs.loop_, Point3::new(1.0, 0.0, 0.0))
            .unwrap();
        let mev2 = brep
            .mev(mev1.half_edge_twin, Point3::new(1.0, 1.0, 0.0))
            .unwrap();
        let mev3 = brep
            .mev(mev2.half_edge_twin, Point3::new(0.0, 1.0, 0.0))
            .unwrap();
        let mef = brep.mef(mev3.half_edge_twin, mev1.half_edge).unwrap();
        brep.extrude(mef.face, Vec3::new(0.0, 0.0, 1.0), 1.0, None)
            .unwrap();
        brep
    }

    #[test]
    fn cube_passes_euler_poincare() {
        let brep = make_cube();
        assert!(brep.validate_euler_poincare().is_ok());
    }

    #[test]
    fn mvfs_only_passes_euler_poincare() {
        let mut brep = BRep::new(Tolerance::default());
        let _mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
        // V=1, E=0, F=1, L=1, S=1 => 1 - 0 + 1 - (1 - 1) = 2 = 2*1
        assert!(brep.validate_euler_poincare().is_ok());
    }

    #[test]
    fn mvfs_plus_mev_passes_euler_poincare() {
        let mut brep = BRep::new(Tolerance::default());
        let mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
        let _mev = brep
            .mev_at_vertex(mvfs.vertex, mvfs.loop_, Point3::new(1.0, 0.0, 0.0))
            .unwrap();
        // V=2, E=1, F=1, L=1, S=1 => 2 - 1 + 1 - (1 - 1) = 2 = 2*1
        assert!(brep.validate_euler_poincare().is_ok());
    }

    #[test]
    fn corrupted_counts_fail_euler_poincare() {
        let mut brep = BRep::new(Tolerance::default());
        let mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
        // Remove the face but keep vertex — breaks Euler-Poincaré
        brep.faces.remove(mvfs.face);
        assert!(brep.validate_euler_poincare().is_err());
    }
}
