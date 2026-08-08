//! Geometric and topological validation for BRep models.

pub mod bvh;
mod closure;

pub use closure::NonManifoldEdgeReport;
mod curve_surface;
mod euler_poincare;
mod half_edge_checks;
mod manifold;
mod orientation;
mod precision;
mod self_intersection;
mod self_intersection_trim;
pub use self_intersection_trim::InteriorOverlapSample;
mod winding_honesty;

pub use winding_honesty::FaceWindingHonesty;

use crate::error::ValidationError;
use crate::topology::BRep;

impl BRep {
    /// Run all validation checks (structural + geometric). Returns on the first failure.
    pub fn validate(&self) -> Result<(), ValidationError> {
        // Structural checks (1-6). Euler–Poincaré is the GENUS-AWARE
        // form: solids with through-holes (a bored plate, a torus)
        // are legitimate kernel outputs — the strict genus-0 form
        // rejected the kernel's own torus constructor. Contexts that
        // require sphere-like topology call
        // `validate_euler_poincare()` directly.
        self.validate_euler_poincare_any_genus()?;
        self.validate_half_edge_pairing()?;
        self.validate_half_edge_cycles()?;
        self.validate_vertex_consistency()?;
        self.validate_manifold_edges()?;
        self.validate_vertex_valence()?;
        // Geometric checks (7-10). Check 9 (curve-surface) is NOT
        // in the fast path: pipelines carry transient curve lies
        // mid-rewrite that later stages repair (kve-chi-contract
        // analogue) and acceptance ladders gate on this method;
        // serving boundaries (Part ops, validate_all) arm it.
        self.validate_face_orientation()?;
        self.validate_loop_closure()?;
        self.validate_watertight()?;
        // Precision checks (11)
        self.validate_precision_fast()?;
        // Note: validate_self_intersection is intentionally excluded from
        // the fast validation path. It is O(n²) in face count, involves
        // surface-surface intersection solving, and should only be used
        // in thorough validation (validate_all). Builders and operators
        // use this method for quick correctness checks.
        Ok(())
    }

    /// Collect all errors from the FAST checks — everything in
    /// [`BRep::validate`] plus the Gregory fold check, with thorough
    /// precision, but WITHOUT the O(n²) surface-surface
    /// self-intersection sweep. This is the collecting variant for
    /// in-pipeline acceptance gates (blend/apex candidate screens),
    /// which inspect structural error classes and cannot afford the
    /// full sweep: on G² NURBS blend sweeps `validate_all` does not
    /// terminate in practical time (PK-094: >16 CPU-minutes on a
    /// 5-fillet pyramid apex).
    pub fn validate_collect_fast(&self) -> Vec<ValidationError> {
        let mut errors = Vec::new();
        if let Err(e) = self.validate_euler_poincare_any_genus() {
            errors.push(e);
        }
        if let Err(e) = self.validate_half_edge_pairing() {
            errors.push(e);
        }
        if let Err(e) = self.validate_half_edge_cycles() {
            errors.push(e);
        }
        if let Err(e) = self.validate_vertex_consistency() {
            errors.push(e);
        }
        if let Err(e) = self.validate_manifold_edges() {
            errors.push(e);
        }
        if let Err(e) = self.validate_vertex_valence() {
            errors.push(e);
        }
        if let Err(e) = self.validate_face_orientation() {
            errors.push(e);
        }
        if let Err(e) = self.validate_loop_closure() {
            errors.push(e);
        }
        if let Err(e) = self.validate_watertight() {
            errors.push(e);
        }
        errors.extend(self.validate_precision_thorough());
        if let Err(e) = self.validate_gregory_self_intersection() {
            errors.push(e);
        }
        errors
    }

    /// Referential integrity of the topology DAG: every loop's face
    /// must exist and own that loop (as outer or inner), and every
    /// half-edge's loop must exist. The strict Euler arithmetic used
    /// to catch dangling references only as a count imbalance; the
    /// genus-aware form (correctly) tolerates the counts, so the
    /// reference structure is validated directly.
    pub fn validate_topology_references(&self) -> Result<(), ValidationError> {
        // Anchored at HALF-EDGES: only REACHABLE topology is
        // validated. A half-edge whose loop is missing, or whose
        // loop's face is missing, is walkable corruption (the
        // mutilated-solid class: a face removed while its boundary
        // ring stays wired).
        //
        // Placement: `validate_all` ONLY — deliberately NOT in the
        // fast `validate()` / `validate_collect_fast` paths, which
        // construction operators use as MID-SURGERY checkpoints
        // where a face legitimately does not exist yet (measured:
        // running this at checkpoints flipped the pyramid3_mixed
        // exact-boundary route from Accepted to Declined). Final
        // results must pass it; intermediate states need not.
        for (heh, he) in self.half_edges.iter() {
            let Some(l) = self.loops.get(he.loop_) else {
                return Err(ValidationError::DanglingReference {
                    what: format!("half-edge {heh:?} references missing loop {:?}", he.loop_),
                });
            };
            if !self.faces.contains_key(l.face) {
                return Err(ValidationError::DanglingReference {
                    what: format!(
                        "half-edge {heh:?} walks loop {:?} whose face {:?} is missing",
                        he.loop_, l.face
                    ),
                });
            }
        }
        Ok(())
    }

    /// Run all validation checks and collect all errors.
    pub fn validate_all(&self) -> Vec<ValidationError> {
        let mut errors = Vec::new();

        // Structural checks (1-6; genus-aware Euler, see validate())
        if let Err(e) = self.validate_euler_poincare_any_genus() {
            errors.push(e);
        }
        if let Err(e) = self.validate_topology_references() {
            errors.push(e);
        }
        if let Err(e) = self.validate_half_edge_pairing() {
            errors.push(e);
        }
        if let Err(e) = self.validate_half_edge_cycles() {
            errors.push(e);
        }
        if let Err(e) = self.validate_vertex_consistency() {
            errors.push(e);
        }
        if let Err(e) = self.validate_manifold_edges() {
            errors.push(e);
        }
        if let Err(e) = self.validate_vertex_valence() {
            errors.push(e);
        }
        // Geometric checks (7-10)
        if let Err(e) = self.validate_face_orientation() {
            errors.push(e);
        }
        // Check 7's deferred half, armed at the serving boundary
        // (PK-245): the stored orientation flag must agree with the
        // measured outer-loop winding on every determinate face.
        // Mid-pipeline states may carry transient dishonesty that
        // operator tails reconcile, so this is validate_all-only.
        errors.extend(self.validate_winding_honesty_all());
        if let Err(e) = self.validate_loop_closure() {
            errors.push(e);
        }
        // Check 9 uses the collect-all form (PK-235): the early-exit
        // form reported one error carrying the FIRST deviation found,
        // so a model with a dozen broken edges looked like one mild
        // defect (saddle fixture: reported 4.685e-4, true worst
        // ~1.3e-1).
        errors.extend(self.validate_curve_surface_all());
        if let Err(e) = self.validate_watertight() {
            errors.push(e);
        }
        // Precision checks (11)
        errors.extend(self.validate_precision_thorough());
        // Geometric check 12: self-intersection (most expensive — run last,
        // and only when there are no structural errors to avoid false reports
        // from broken topology).
        if errors.is_empty()
            && let Err(e) = self.validate_self_intersection()
        {
            errors.push(e);
        }
        // §5.3 item #10 — within-face Gregory self-intersection check.
        // Fast (O(face count) with bounded per-face work); run
        // unconditionally on Gregory-carrying BReps even when other
        // errors are present, because a folded Gregory sub-patch often
        // CAUSES the downstream manifold / intersection errors.
        if let Err(e) = self.validate_gregory_self_intersection() {
            errors.push(e);
        }

        errors
    }
}

#[cfg(test)]
#[path = "mod_tests.rs"]
mod tests;
