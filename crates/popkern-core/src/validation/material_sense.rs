//! Component material sense — which way the flag-folded normals face.
//!
//! No per-face orientation witness says which side of a face is
//! MATERIAL: `SurfaceOnFace::orientation` is honest to the loop
//! winding (see `winding_honesty`), but a component's loops may wind
//! either way globally, so the folded normal `flag ? n : -n` is
//! uniformly outward-of-material on some bodies and uniformly inward
//! on others (measured: every face of a served filleted cube folds
//! INWARD — planes with obvious material sides included). The mesher
//! already resolves this with the sign of the folded-normal wedge
//! volume as its final authority; this module makes the same
//! resolution available to any consumer that needs out-of-material
//! directions (the serving footprint repair pushed control points
//! INTO the material for lack of it).

use crate::geometry::surface::SurfaceEval;
use crate::topology::{BRep, FaceHandle};

/// Which way a component's flag-folded surface normals face relative
/// to its material.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MaterialSense {
    /// `flag ? n : -n` points OUT of the material.
    FoldedOutward,
    /// `flag ? n : -n` points INTO the material.
    FoldedInward,
}

impl MaterialSense {
    /// Multiplier converting a flag-folded direction into an
    /// out-of-material direction.
    pub fn outward_multiplier(self) -> f64 {
        match self {
            MaterialSense::FoldedOutward => 1.0,
            MaterialSense::FoldedInward => -1.0,
        }
    }
}

impl BRep {
    /// Resolve the material sense of the closed component containing
    /// `seed`: the sign of the folded-normal wedge integral
    /// (divergence theorem — positive volume means the folded normals
    /// face outward). Faces are sampled on a coarse grid over their
    /// exact bounded parameter windows with trim-interior
    /// classification, so the estimate honours the served trims; the
    /// magnitude is approximate but the SIGN is what the contract
    /// promises. Requires honest orientation flags (the armed
    /// winding-honesty checks guarantee them on served bodies).
    /// Returns `None` when the component's integral is too small to
    /// call — degenerate, open, or unsampleable geometry.
    pub fn component_material_sense(&self, seed: FaceHandle) -> Option<MaterialSense> {
        // Component = faces reachable over shared edges.
        let mut component = Vec::new();
        let mut visited = std::collections::BTreeSet::new();
        let mut queue = vec![seed];
        while let Some(face) = queue.pop() {
            if !visited.insert(face) {
                continue;
            }
            component.push(face);
            let Some(data) = self.faces.get(face) else {
                continue;
            };
            let Some(outer_loop) = self.loops.get(data.outer_loop) else {
                continue;
            };
            let start = outer_loop.half_edge;
            let mut he = start;
            // Open sheets have loopless (or absent) twins at free
            // rims — no neighbour across that edge.
            while let Some(he_data) = self.half_edges.get(he) {
                if let Some(twin_data) = he_data.mate().and_then(|m| self.half_edges.get(m))
                    && let Some(partner_loop) = self.loops.get(twin_data.loop_)
                {
                    let partner = partner_loop.face;
                    if !visited.contains(&partner) {
                        queue.push(partner);
                    }
                }
                he = he_data.next;
                if he == start {
                    break;
                }
            }
        }

        let wedge = self.folded_wedge_volume(&component)?;
        Some(if wedge > 0.0 {
            MaterialSense::FoldedOutward
        } else {
            MaterialSense::FoldedInward
        })
    }

    /// The folded-normal wedge integral (the enclosed volume, signed)
    /// of the closed face set `faces`: positive when their folded
    /// normals face out of the volume they bound, negative when in.
    /// Sampled as [`BRep::component_material_sense`] describes (a 6 × 6
    /// grid per face, samples within the trim classifier's boundary
    /// band left out), so the magnitude is approximate and the sign is
    /// the contract. `None` when the integral is too small to call (an
    /// open or degenerate set, unsampleable geometry).
    pub fn folded_wedge_volume(&self, faces: &[FaceHandle]) -> Option<f64> {
        let band = self.sense_parametric_band();
        self.wedge_integral(&unsigned(faces), 6, band)
    }

    /// [`BRep::folded_wedge_volume`] sampled finely enough to judge a
    /// body by: a 12 × 12 grid per face, each sample counted by strict
    /// trim containment. Leaving out the samples within the boundary
    /// band drops a margin of every face from the integral, a tenth of
    /// a unit cube's at a 1e-2 tolerance (measured: a cube less a
    /// sphere read −1.08 where it encloses +0.2).
    pub fn enclosed_folded_volume(&self, faces: &[FaceHandle]) -> Option<f64> {
        self.wedge_integral(&unsigned(faces), 12, 1.0e-12)
    }

    /// [`BRep::enclosed_folded_volume`] of face sides: each face counted
    /// with its sign, +1 where the region lies behind it, −1 where it
    /// lies in front (a closed set of sides bounding one region, read
    /// whichever way each face points).
    pub fn enclosed_signed_volume(&self, sides: &[(FaceHandle, f64)]) -> Option<f64> {
        self.wedge_integral(sides, 12, 1.0e-12)
    }

    /// The trim classifier's parametric boundary band
    /// (`face_trim_classify_uv`).
    fn sense_parametric_band(&self) -> f64 {
        self.tolerance
            .parametric_base()
            .max(self.tolerance.linear())
            .max(1.0e-9)
            * 10.0
    }

    /// The wedge integral over `faces` on a `grid` × `grid` sample grid
    /// per face, a sample counted when the trim classifier with
    /// parametric band `classify_band` reads it interior.
    fn wedge_integral(
        &self,
        faces: &[(FaceHandle, f64)],
        grid: usize,
        classify_band: f64,
    ) -> Option<f64> {
        let band = self.sense_parametric_band();
        // The wedge integral of a closed surface is the same about any
        // centre, but the coarse grid leaves a residual ∫n dA ≠ 0 that
        // the centre multiplies: about the world origin, a thin body
        // far from it reads the residual times its offset (measured: a
        // 0.04-radius bore at x = 0.61 read FoldedInward at two of
        // eight holes and the cell complex kept the plate's outside).
        // About the component's boundary-sample centroid the residual
        // scales with the body, not its placement.
        let centre = {
            let mut sum = crate::geometry::Vec3::new(0.0, 0.0, 0.0);
            let mut count = 0usize;
            for &(face, _) in faces {
                for p in self.face_boundary_sample_points(face, 2) {
                    sum = sum + crate::geometry::Vec3::new(p.x, p.y, p.z);
                    count += 1;
                }
            }
            if count == 0 {
                return None;
            }
            sum * (1.0 / count as f64)
        };
        let mut wedge = 0.0;
        let mut magnitude = 0.0;
        let mut sampled = 0usize;
        for &(face, sign) in faces {
            let Some(sof) = self.faces.get(face).and_then(|f| f.surface.as_ref()) else {
                continue;
            };
            // Boundary CURVE samples, not just vertices: a min-vertex
            // analytic body (cylinder: two seam vertices, caps: one rim
            // vertex each) collapses every vertex-derived window to zero
            // width, which used to make the whole component unsampleable
            // (None) — the sign contract needs the loop's real extent.
            let positions = self.face_boundary_sample_points(face, 8);
            let ((u0, mut u1), (v0, mut v1)) =
                crate::geometry::uv_projection::exact_bounded_parameter_range(
                    &positions,
                    &sof.geometry,
                );
            if !(u1 - u0).is_finite() || !(v1 - v0).is_finite() {
                continue;
            }
            // A seam-only boundary un-projects to one station on a
            // periodic axis, collapsing that window to zero width —
            // but a face whose boundary lies entirely on the seam IS
            // the full wrap (the 1-face sphere's meridian, a closed
            // torus chart), so the exact window is the full period.
            // "Zero width" is the trim classifier's boundary band: the
            // seam's projections scatter by ulps (measured: a pip
            // sphere's window 2.2e-15 wide, past an ε·period cutoff),
            // and a window inside the band has no interior to sample.
            let (period_u, period_v) = sof.geometry.parameter_periodicity();
            if let Some(p) = period_u
                && (u1 - u0).abs() < band
            {
                u1 = u0 + p;
            }
            if let Some(p) = period_v
                && (v1 - v0).abs() < band
            {
                v1 = v0 + p;
            }
            let du = (u1 - u0) / grid as f64;
            let dv = (v1 - v0) / grid as f64;
            if du.abs() < f64::EPSILON || dv.abs() < f64::EPSILON {
                continue;
            }
            for i in 0..grid {
                for j in 0..grid {
                    let u = u0 + du * (i as f64 + 0.5);
                    let v = v0 + dv * (j as f64 + 0.5);
                    let inside = matches!(
                        self.face_trim_classify_uv_band(
                            face,
                            &sof.geometry,
                            super::self_intersection_trim::UvPointForSense { u, v },
                            classify_band,
                        ),
                        Some(super::self_intersection_trim::TrimContainmentForSense::Interior)
                    );
                    if !inside {
                        continue;
                    }
                    let p = sof.geometry.point_at(u, v);
                    let su = sof.geometry.partial_u(u, v);
                    let sv = sof.geometry.partial_v(u, v);
                    let n = su.cross(&sv);
                    let folded = if sof.orientation { n } else { n * -1.0 };
                    // p · n dA with n dA = folded * du * dv (the volume
                    // is a third of its integral).
                    let (x, y, z) = (p.x - centre.x, p.y - centre.y, p.z - centre.z);
                    let contribution =
                        sign * (x * folded.x + y * folded.y + z * folded.z) * du * dv;
                    wedge += contribution;
                    magnitude += contribution.abs();
                    sampled += 1;
                }
            }
        }
        if sampled == 0 {
            return None;
        }
        // Indeterminate when the integral cancels to a sliver of its
        // own terms (an open sheet, a degenerate body), or encloses
        // less than a tolerance cube. The floor is relative: an
        // absolute (100 · linear)³ read a 0.34-volume cone as
        // indeterminate at a 1e-2 tolerance.
        // The volume is a third of the wedge integral (the divergence
        // of p is 3): read whole, every magnitude was three times the
        // volume (a unit cube 3.0), which the sign checks never saw.
        let volume = wedge / 3.0;
        let linear = self.tolerance.linear();
        if wedge.abs() < 1e-6 * magnitude || volume.abs() < linear * linear * linear {
            return None;
        }
        Some(volume)
    }
}

/// `faces`, each counted once forward.
fn unsigned(faces: &[FaceHandle]) -> Vec<(FaceHandle, f64)> {
    faces.iter().map(|&f| (f, 1.0)).collect()
}

#[cfg(test)]
#[path = "material_sense_tests.rs"]
mod tests;
