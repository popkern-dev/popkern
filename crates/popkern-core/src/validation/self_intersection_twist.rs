//! Per-face fold validators (extracted from `self_intersection.rs`
//! under the file-size ratchet): the Jacobian-sign twist detector and
//! the Gregory self-intersection report check.

use crate::error::ValidationError;
use crate::geometry::surface::SurfaceEval;
use crate::topology::BRep;

impl BRep {
    /// Check for twisted surfaces: faces where the normal reverses direction.
    ///
    /// Samples the surface normal on a grid and checks for sign changes in
    /// the determinant of the Jacobian (partial_u × partial_v). A sign
    /// change indicates the surface folds over itself within the face.
    pub fn validate_twisted_surfaces(&self) -> Result<(), ValidationError> {
        for (fh, face) in &self.faces {
            let surface = match face.surface.as_ref() {
                Some(s) => &s.geometry,
                None => continue,
            };

            // Get face parameter range
            let positions = self.face_outer_loop_positions(fh);
            if positions.is_empty() {
                continue;
            }
            let ((u_min, u_max), (v_min, v_max)) =
                crate::geometry::uv_projection::tight_bounded_parameter_range(&positions, surface);

            // Sample normals on a 4x4 grid
            let n = 4;
            let mut first_sign = None;
            let mut twisted = false;

            for i in 0..n {
                let u = u_min + (u_max - u_min) * (i as f64 + 0.5) / (n as f64);
                for j in 0..n {
                    let v = v_min + (v_max - v_min) * (j as f64 + 0.5) / (n as f64);
                    let normal = surface.normal_at(u, v);
                    let len = normal.length();
                    if len < 1e-12 {
                        continue; // degenerate normal, skip
                    }

                    // Use z-component sign as a proxy for orientation
                    // (for general surfaces, compare against first sample's direction)
                    let sign = match &first_sign {
                        None => {
                            first_sign = Some(normal);
                            true
                        }
                        Some(first) => {
                            let dot = normal.dot(first);
                            dot > 0.0
                        }
                    };

                    if first_sign.is_some() && !sign {
                        twisted = true;
                        break;
                    }
                }
                if twisted {
                    break;
                }
            }

            if twisted {
                return Err(ValidationError::TwistedSurface { face: fh });
            }
        }

        Ok(())
    }

    /// §5.3 item #10 — check every Gregory-surface face for
    /// within-face self-intersection (fold-over or curl-around).
    /// Delegates per-face to
    /// [`GregoryPatch::self_intersects_with_tolerance`].
    pub fn validate_gregory_self_intersection(&self) -> Result<(), ValidationError> {
        use crate::geometry::surface::SurfaceGeometry;

        let tol_linear = self.tolerance.linear().max(1e-6);
        for (fh, face) in &self.faces {
            let Some(sof) = face.surface.as_ref() else {
                continue;
            };
            let SurfaceGeometry::Gregory(ref patch) = sof.geometry else {
                continue;
            };
            if patch.self_intersects_with_tolerance(tol_linear) {
                return Err(ValidationError::TwistedSurface { face: fh });
            }
        }
        Ok(())
    }
}
