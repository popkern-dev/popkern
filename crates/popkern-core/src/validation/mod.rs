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
        // Geometric checks (7-10). Check 9 (curve-surface
        // consistency) is intentionally NOT part of the fast path:
        // construction pipelines legitimately carry transient curve
        // lies mid-rewrite that later stages repair (the same
        // mid-pipeline exemption as the kve chi contract), and the
        // per-stage acceptance ladders gate on this method. Serving
        // boundaries (Part ops, export, validate_all) arm it.
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
        if let Err(e) = self.validate_loop_closure() {
            errors.push(e);
        }
        if let Err(e) = self.validate_curve_surface() {
            errors.push(e);
        }
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
mod tests {
    use crate::geometry::CurveGeometry;
    use crate::geometry::curve::CurveEval;
    use crate::geometry::{Point3, Vec3};
    use crate::tolerance::Tolerance;
    use crate::topology::BRep;

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
    fn cube_passes_full_validation() {
        let brep = make_cube();
        assert!(brep.validate().is_ok());
    }

    #[test]
    fn cube_validate_all_returns_no_errors() {
        let brep = make_cube();
        assert!(brep.validate_all().is_empty());
    }

    #[test]
    fn corrupted_brep_validate_all_collects_errors() {
        let mut brep = BRep::new(Tolerance::default());
        let mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
        let mev = brep
            .mev_at_vertex(mvfs.vertex, mvfs.loop_, Point3::new(1.0, 0.0, 0.0))
            .unwrap();
        // Corrupt the twin pointer.
        brep.half_edges[mev.half_edge].twin = Default::default();
        let errors = brep.validate_all();
        // Should have at least 1 error (broken twin).
        assert!(!errors.is_empty());
    }

    #[test]
    fn mvfs_only_passes_validation() {
        let mut brep = BRep::new(Tolerance::default());
        let _mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
        assert!(brep.validate().is_ok());
    }

    // --- Task 20: geometric checks ---

    #[test]
    fn geometric_checks_pass_on_topology_only_brep() {
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

        assert!(brep.validate_face_orientation().is_ok());
        assert!(brep.validate_loop_closure().is_ok());
        assert!(brep.validate_curve_surface().is_ok());
        assert!(brep.validate_watertight().is_ok());
    }

    #[test]
    fn full_validate_passes_for_cube() {
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

        assert!(brep.validate().is_ok());
    }

    #[test]
    fn validate_all_collects_multiple_errors() {
        let mut brep = BRep::new(Tolerance::default());
        let mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
        let mev = brep
            .mev_at_vertex(mvfs.vertex, mvfs.loop_, Point3::new(1.0, 0.0, 0.0))
            .unwrap();
        // Corrupt twin AND next to trigger multiple validation errors
        brep.half_edges[mev.half_edge].twin = Default::default();
        brep.half_edges[mev.half_edge_twin].twin = Default::default();
        let errors = brep.validate_all();
        // Should have multiple errors (broken twin, broken cycle, manifold, etc.)
        assert!(
            errors.len() >= 2,
            "Expected at least 2 errors, got {}",
            errors.len()
        );
    }

    #[test]
    fn validate_all_hits_all_check_branches() {
        use crate::geometry::surface::{SurfaceGeometry, SurfaceOnFace};

        // Build a cube with geometry, then corrupt it to fail many checks
        let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();

        // Remove a face to break Euler-Poincare (check 1)
        let face_h = brep.faces.keys().next().unwrap();
        brep.faces.remove(face_h);

        // Corrupt a twin pointer to break pairing (check 2)
        let he_h = brep.half_edges.keys().next().unwrap();
        brep.half_edges[he_h].twin = Default::default();

        // Attach degenerate surface to trigger orientation failure (check 7)
        let face2_h = brep.faces.keys().next().unwrap();
        brep.faces[face2_h].surface = Some(SurfaceOnFace {
            geometry: SurfaceGeometry::Sphere {
                center: Point3::new(0.0, 0.0, 0.0),
                radius: 0.0,
            },
            orientation: true,
        });

        let errors = brep.validate_all();
        // Should have errors from multiple check categories
        assert!(
            errors.len() >= 3,
            "Expected at least 3 errors from different checks, got {}",
            errors.len()
        );
    }

    #[test]
    fn validate_all_on_cube_with_surface_geometry() {
        let cube =
            crate::topology::BRep::unit_cube(crate::tolerance::Tolerance::default()).unwrap();
        let errors = cube.validate_all();
        assert!(
            errors.is_empty(),
            "Cube with geometry should pass all validation: {:?}",
            errors
        );
    }

    #[test]
    fn loop_closure_from_vertex_positions() {
        let mut brep = BRep::new(Tolerance::default());
        let mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
        let mev1 = brep
            .mev_at_vertex(mvfs.vertex, mvfs.loop_, Point3::new(1.0, 0.0, 0.0))
            .unwrap();
        let mev2 = brep
            .mev(mev1.half_edge_twin, Point3::new(0.5, 1.0, 0.0))
            .unwrap();
        let _mef = brep.mef(mev2.half_edge_twin, mev1.half_edge).unwrap();

        assert!(brep.validate_loop_closure().is_ok());
    }

    #[test]
    fn validate_all_hits_vertex_consistency_error() {
        let mut brep = BRep::new(Tolerance::default());
        let mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
        let mev = brep
            .mev_at_vertex(mvfs.vertex, mvfs.loop_, Point3::new(1.0, 0.0, 0.0))
            .unwrap();
        // Corrupt vertex outgoing to trigger vertex consistency error
        brep.vertices[mvfs.vertex].outgoing = mev.half_edge_twin;
        let errors = brep.validate_all();
        assert!(!errors.is_empty());
    }

    #[test]
    fn validate_all_hits_vertex_valence_error() {
        // Build a single-edge topology then corrupt one half-edge's loop
        // so the twins are on *different* faces, making the seam exemption
        // inapplicable and triggering LowValence.
        let mut brep = BRep::new(Tolerance::default());
        let mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
        let mev = brep
            .mev_at_vertex(mvfs.vertex, mvfs.loop_, Point3::new(1.0, 0.0, 0.0))
            .unwrap();
        // Create a second loop/face and move one half-edge into it so
        // the seam check (same-face) fails.
        let new_loop = brep.loops.insert(crate::topology::Loop {
            half_edge: mev.half_edge_twin,
            face: Default::default(),
        });
        let new_face = brep.faces.insert(crate::topology::Face {
            outer_loop: new_loop,
            inner_loops: Vec::new(),
            shell: brep.faces[mvfs.face].shell,
            surface: None,
            name: None,
            tolerance: None,
            role: Default::default(),
        });
        brep.loops[new_loop].face = new_face;
        brep.half_edges[mev.half_edge_twin].loop_ = new_loop;
        // Now the edge's two half-edges are on different faces => valence 1, no seam exemption
        let errors = brep.validate_all();
        assert!(
            errors
                .iter()
                .any(|e| matches!(e, crate::error::ValidationError::LowValence { .. }))
        );
    }

    #[test]
    fn validate_all_collects_loop_closure_error() {
        // Build a cube with geometry, then corrupt a vertex to cause loop closure failure
        let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
        let v_handle = brep.vertices.keys().next().unwrap();
        brep.vertices[v_handle].position = Point3::new(100.0, 100.0, 100.0);
        let errors = brep.validate_all();
        assert!(!errors.is_empty());
    }

    #[test]
    fn validate_early_return_on_first_error() {
        // validate() should return on the first failure (Euler-Poincare)
        let mut brep = BRep::new(Tolerance::default());
        let mvfs = brep.mvfs(Point3::new(0.0, 0.0, 0.0)).unwrap();
        let mev = brep
            .mev_at_vertex(mvfs.vertex, mvfs.loop_, Point3::new(1.0, 0.0, 0.0))
            .unwrap();
        // Corrupt twin to break pairing => validate fails on check 2
        brep.half_edges[mev.half_edge].twin = Default::default();
        let result = brep.validate();
        assert!(result.is_err());
    }

    #[test]
    fn validate_all_skips_self_intersection_when_structural_errors() {
        // validate_all should skip self-intersection check (line 73-76)
        // when there are structural errors
        let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();

        // Remove a face to break Euler formula — structural error
        let face_h = brep.faces.keys().next().unwrap();
        brep.faces.remove(face_h);

        let errors = brep.validate_all();
        // Should have structural errors, but NOT a self-intersection error
        assert!(!errors.is_empty());
        assert!(
            !errors
                .iter()
                .any(|e| matches!(e, crate::error::ValidationError::SelfIntersection { .. })),
            "Self-intersection check should be skipped when structural errors exist"
        );
    }

    #[test]
    fn validate_all_runs_self_intersection_when_no_structural_errors() {
        // A valid cube should have no self-intersection
        // This test verifies the self_intersection branch of validate_all is reached.
        let brep = BRep::unit_cube(Tolerance::default()).unwrap();
        let errors = brep.validate_all();
        // Valid cube: no errors (including self-intersection check ran successfully)
        assert!(
            errors.is_empty(),
            "Valid cube should pass all checks including self-intersection: {:?}",
            errors
        );
    }

    #[test]
    fn validate_all_collects_degenerate_surface_error() {
        use crate::geometry::surface::{SurfaceGeometry, SurfaceOnFace};

        // Build a cube with geometry, then set a face's surface to a degenerate
        // zero-radius sphere -- normal is zero everywhere, failing orientation.
        let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
        let face_h = brep.faces.keys().next().unwrap();

        brep.faces[face_h].surface = Some(SurfaceOnFace {
            geometry: SurfaceGeometry::Sphere {
                center: Point3::new(0.0, 0.0, 0.0),
                radius: 0.0,
            },
            orientation: true,
        });

        let errors = brep.validate_all();
        assert!(!errors.is_empty());
    }

    #[test]
    fn validate_all_collects_watertight_error() {
        // Build a cube with geometry, then corrupt an edge's vertex to trigger watertight failure
        let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
        // Find an edge with curve geometry
        let edge_h = brep
            .edges
            .keys()
            .find(|&eh| brep.edges[eh].curve.is_some())
            .unwrap();
        let he_h = brep.edges[edge_h].half_edge;
        let origin = brep.half_edges[he_h].origin;
        // Move the vertex far away
        brep.vertices[origin].position = Point3::new(999.0, 999.0, 999.0);

        let errors = brep.validate_all();
        assert!(!errors.is_empty());
    }

    #[test]
    fn validate_watertight_accepts_reversed_edge_curve_direction() {
        let mut brep = BRep::unit_cube(Tolerance::default()).unwrap();
        let edge = brep
            .edges
            .keys()
            .find(|&eh| brep.edges[eh].curve.is_some())
            .unwrap();
        let curve = brep.edges[edge].curve.as_mut().unwrap();
        let p_start = curve.geometry.point_at(curve.t_start);
        let p_end = curve.geometry.point_at(curve.t_end);
        curve.geometry = CurveGeometry::Line {
            origin: p_end,
            direction: p_start - p_end,
        };
        curve.t_start = 0.0;
        curve.t_end = 1.0;

        assert!(brep.validate_watertight().is_ok());
    }
}
