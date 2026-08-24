//! Pcurve-presence audit (PK-331).
//!
//! A face that carries surface geometry but whose loop half-edges lack
//! pcurves is trim-unreadable: the winding-honesty audit cannot measure
//! it (and defers to its tangency witness), tessellation falls back to
//! untrimmed full-grid faces, and STEP export loses the trim. `thicken`
//! served exactly such bodies for an unknown period — every validator
//! stayed silent. This is a serving-boundary audit (`validate_all`
//! only): mid-pipeline states legitimately hold faces whose pcurves are
//! cleared for re-authoring.

use crate::error::ValidationError;
use crate::topology::BRep;

impl BRep {
    /// One [`ValidationError::MissingPcurves`] per face that has
    /// surface geometry but loop half-edges without pcurves.
    pub(crate) fn validate_pcurve_presence_all(&self) -> Vec<ValidationError> {
        let mut errors = Vec::new();
        for (fh, face) in &self.faces {
            if face.surface.is_none() {
                continue;
            }
            let mut missing = 0_usize;
            let mut total = 0_usize;
            let mut loops = vec![face.outer_loop];
            loops.extend(face.inner_loops.iter().copied());
            for lh in loops {
                let Some(loop_data) = self.loops.get(lh) else {
                    continue;
                };
                let start = loop_data.half_edge;
                if !self.half_edges.contains_key(start) {
                    continue;
                }
                let mut cur = start;
                for _ in 0..self.half_edges.len() + 1 {
                    total += 1;
                    if self.half_edges[cur].pcurve.is_none() {
                        missing += 1;
                    }
                    cur = self.half_edges[cur].next;
                    if cur == start || !self.half_edges.contains_key(cur) {
                        break;
                    }
                }
            }
            if missing > 0 {
                errors.push(ValidationError::MissingPcurves {
                    face: fh,
                    missing,
                    total,
                });
            }
        }
        errors
    }
}

#[cfg(test)]
mod tests {
    use crate::error::ValidationError;
    use crate::tolerance::Tolerance;
    use crate::topology::BRep;

    #[test]
    fn constructors_carry_full_pcurves() {
        for brep in [
            BRep::unit_cube(Tolerance::default()).unwrap(),
            BRep::cylinder(1.0, 2.0, Tolerance::default()).unwrap(),
            BRep::cone(1.0, 2.0, Tolerance::default()).unwrap(),
        ] {
            assert!(brep.validate_pcurve_presence_all().is_empty());
        }
    }

    #[test]
    fn stripped_pcurves_are_reported_per_face() {
        let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
        let hes: Vec<_> = brep.half_edges.keys().collect();
        for heh in hes {
            brep.half_edges[heh].pcurve = None;
        }
        let errors = brep.validate_pcurve_presence_all();
        assert_eq!(errors.len(), 6, "one report per cube face");
        assert!(errors.iter().all(|e| matches!(
            e,
            ValidationError::MissingPcurves {
                missing: 4,
                total: 4,
                ..
            }
        )));
        // And validate_all surfaces the class.
        assert!(
            brep.validate_all()
                .iter()
                .any(|e| matches!(e, ValidationError::MissingPcurves { .. })),
            "validate_all must include the pcurve-presence audit"
        );
    }
}
