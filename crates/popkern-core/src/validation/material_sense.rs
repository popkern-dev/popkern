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
            let start = self.loops[data.outer_loop].half_edge;
            let mut he = start;
            loop {
                let twin = self.half_edges[he].twin;
                let partner = self.loops[self.half_edges[twin].loop_].face;
                if !visited.contains(&partner) {
                    queue.push(partner);
                }
                he = self.half_edges[he].next;
                if he == start {
                    break;
                }
            }
        }

        const GRID: usize = 6;
        let mut wedge = 0.0;
        let mut sampled = 0usize;
        for &face in &component {
            let Some(sof) = self.faces.get(face).and_then(|f| f.surface.as_ref()) else {
                continue;
            };
            // Boundary CURVE samples, not just vertices: a min-vertex
            // analytic body (cylinder: two seam vertices, caps: one rim
            // vertex each) collapses every vertex-derived window to zero
            // width, which used to make the whole component unsampleable
            // (None) — the sign contract needs the loop's real extent.
            let positions = self.face_boundary_sample_points(face, 8);
            let ((u0, u1), (v0, v1)) =
                crate::geometry::uv_projection::exact_bounded_parameter_range(
                    &positions,
                    &sof.geometry,
                );
            if !(u1 - u0).is_finite() || !(v1 - v0).is_finite() {
                continue;
            }
            let du = (u1 - u0) / GRID as f64;
            let dv = (v1 - v0) / GRID as f64;
            if du.abs() < f64::EPSILON || dv.abs() < f64::EPSILON {
                continue;
            }
            for i in 0..GRID {
                for j in 0..GRID {
                    let u = u0 + du * (i as f64 + 0.5);
                    let v = v0 + dv * (j as f64 + 0.5);
                    let inside = matches!(
                        self.face_trim_classify_uv(
                            face,
                            &sof.geometry,
                            super::self_intersection_trim::UvPointForSense { u, v },
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
                    // p · n dA with n dA = folded * du * dv.
                    let contribution = (p.x * folded.x + p.y * folded.y + p.z * folded.z) * du * dv;
                    wedge += contribution;
                    sampled += 1;
                }
            }
        }
        if sampled == 0 {
            return None;
        }
        // Indeterminate when the integral is numerically tiny for the
        // body's scale.
        let scale = self.tolerance.linear() * 100.0;
        if wedge.abs() < scale * scale * scale {
            return None;
        }
        Some(if wedge > 0.0 {
            MaterialSense::FoldedOutward
        } else {
            MaterialSense::FoldedInward
        })
    }
}

#[cfg(test)]
#[path = "material_sense_tests.rs"]
mod tests;
