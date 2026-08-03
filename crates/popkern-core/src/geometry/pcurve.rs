//! Parametric curves in surface UV space.
//!
//! `Pcurve` is the analogue of [`crate::geometry::curve::CurveOnEdge`]
//! for parametric (UV) curves on a surface. Stored either as an
//! analytic primitive (for trig surfaces where the canonical pcurve
//! is exact — e.g. a linear segment in cyl UV traces a circle on
//! the cylinder) or as a general NURBS for non-analytic
//! intersections. UV coordinates are encoded as `Point3 { x: u, y:
//! v, z: 0.0 }` by convention, matching popkern's existing pcurve
//! storage.
//!
//! See ADR `docs/decisions/2026-04-25-edge-triple-factory.md` for
//! the design rationale (B1: analytic pcurve types, commercial-style).

use serde::{Deserialize, Serialize};

use crate::geometry::nurbs_curve::NurbsCurve;
use crate::geometry::{Point3, Vec3};

#[path = "pcurve_nurbs_export.rs"]
mod pcurve_nurbs_export;
#[path = "pcurve_span.rs"]
mod pcurve_span;

/// A parametric curve binding to one side of an edge.
///
/// Holds the geometric form plus the parameter range
/// `[t_start, t_end]` over which the curve is defined.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Pcurve {
    /// The geometric form (analytic or NURBS).
    pub geometry: PcurveGeometry,
    /// Parameter at the start vertex.
    pub t_start: f64,
    /// Parameter at the end vertex.
    pub t_end: f64,
}

/// Geometric form of a pcurve. Mirrors
/// [`crate::geometry::curve::CurveGeometry`] but with 2D-native
/// analytic variants.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum PcurveGeometry {
    /// Linear segment in UV space: the canonical pcurve for
    /// horizontal circles on cylinders, isoparametric edges on any
    /// surface, and any constant-rate `(u, v)` traversal. Exact by
    /// construction.
    ///
    /// `point_at(t) = start + s · (end - start)` where
    /// `s = (t - t_start) / (t_end - t_start)`.
    Linear {
        /// Start `(u, v)`.
        start: (f64, f64),
        /// End `(u, v)`.
        end: (f64, f64),
    },
    /// Circular or elliptical arc in UV space. The canonical pcurve
    /// for circles on planes (after projection through the plane's
    /// UV axes) and for conic-section intersections.
    ///
    /// At parameter `t`, the angle is
    /// `α = start_angle + s · (end_angle - start_angle)` where
    /// `s = (t - t_start) / (t_end - t_start)`. The `(u, v)` point
    /// is `(u_radius·cos(α), v_radius·sin(α))` rotated by `rotation`
    /// and translated by `center`.
    Elliptical {
        /// Centre `(u, v)`.
        center: (f64, f64),
        /// Radius along the major axis.
        u_radius: f64,
        /// Radius along the minor axis.
        v_radius: f64,
        /// CCW rotation of the major axis from the U direction.
        rotation: f64,
        /// Sweep start angle (radians, in the rotated local frame).
        start_angle: f64,
        /// Sweep end angle (radians).
        end_angle: f64,
    },
    /// Sinusoidal pcurve: `u` is linear in `t`, `v` is
    /// `const + v_cos · cos(t) + v_sin · sin(t)`. Canonical for
    /// oblique plane×cylinder intersections, where the curve
    /// traces a sinusoid in cylinder UV (linear in the angular
    /// parameter, sinusoidal in the axial direction). Exact by
    /// construction.
    ///
    /// `point_at(t) = (u_offset + t, v_const + v_cos · cos(t) + v_sin · sin(t))`
    Sinusoidal {
        /// Constant offset on the linear `u` coordinate; full form
        /// is `u(t) = u_offset + t`.
        u_offset: f64,
        /// Constant term in the sinusoidal `v` coordinate.
        v_const: f64,
        /// Coefficient on `cos(t)` in the sinusoidal `v` term.
        v_cos: f64,
        /// Coefficient on `sin(t)` in the sinusoidal `v` term.
        v_sin: f64,
    },
    /// General NURBS pcurve. Control points stored as
    /// `Point3 { x: u, y: v, z: 0 }`.
    Nurbs(NurbsCurve),
    /// Symmetric rational-quadratic-Bezier circle arc projected
    /// onto cylinder UV.  Encodes the EXACT angular trajectory
    /// `θ(t)` of a `(1, cos α, 1)`-weighted rational quadratic
    /// circle arc, parameterised by Bezier t ∈ [0, 1] (or any
    /// sub-range via `Pcurve::t_start, t_end`).
    ///
    /// With α = (θ_b − θ_a)/2 and θ_m = (θ_a + θ_b)/2:
    /// ```text
    ///   θ(t) = θ_m + atan2( sin α · (2t − 1),
    ///                       cos α · (1 − 2t(1−t)) + 2t(1−t) )
    ///   u(t) = θ(t)
    ///   v(t) = v_const
    /// ```
    /// The trajectory matches `cylinder.point_at(θ(t), v_const)`
    /// to machine precision when paired with a `CurveGeometry::
    /// Nurbs` curve3d that IS the rational quadratic circle arc
    /// `arc_as_nurbs(centre, radius, u_axis, v_axis, θ_a, θ_b)`
    /// (as produced by `sphere_fit_termination` for sphere-fit
    /// FilletEnd substitutions where the cap centre lies on the
    /// cylinder axis).
    ///
    /// Used by Phase 6c-D Stage 4.7d on the cylinder-twin side
    /// of cap-vs-Cyl edges.  A `Linear` pcurve over the same
    /// parameter range straight-line-interpolates `(θ, v_const)`
    /// in UV — but the rational quad's `θ(t)` is non-linear in
    /// `t`, so the linear UV-interpolant lands at angles offset
    /// from the curve3d's actual angle by O(rad²·sin α) per cap-
    /// tangent angle, giving up to 1 cm 3-D dev on the larger-
    /// radius pyramid corner fillets.
    CylinderArcViaRationalQuad {
        /// Start angle of the full-arc rational quadratic
        /// (radians, in the cylinder's `orthonormal_basis` frame).
        theta_a: f64,
        /// End angle of the full-arc rational quadratic.  Must
        /// satisfy `|θ_b − θ_a| ≤ π` for the symmetric weights
        /// `(1, cos((θ_b − θ_a)/2), 1)` representation.
        theta_b: f64,
        /// Axial coordinate on the cylinder (constant for this
        /// arc — sphere-fit termination produces a circle in the
        /// plane perpendicular to the cylinder axis at this
        /// `v_const`).
        v_const: f64,
    },
}

#[path = "pcurve_chart_ops.rs"]
mod pcurve_chart_ops;

impl Pcurve {
    /// Wrap a `NurbsCurve` as a NURBS-form pcurve. Inherits the
    /// curve's intrinsic parameter range.
    pub fn from_nurbs(curve: NurbsCurve) -> Self {
        let (t_start, t_end) = nurbs_param_range(&curve);
        Self {
            geometry: PcurveGeometry::Nurbs(curve),
            t_start,
            t_end,
        }
    }

    /// Build a linear pcurve from `start` to `end` over the given
    /// parameter range `[t_start, t_end]`.
    pub fn linear(start: (f64, f64), end: (f64, f64), t_start: f64, t_end: f64) -> Self {
        Self {
            geometry: PcurveGeometry::Linear { start, end },
            t_start,
            t_end,
        }
    }

    /// Build an elliptical pcurve from explicit parameters.
    #[allow(clippy::too_many_arguments)]
    pub fn elliptical(
        center: (f64, f64),
        u_radius: f64,
        v_radius: f64,
        rotation: f64,
        start_angle: f64,
        end_angle: f64,
        t_start: f64,
        t_end: f64,
    ) -> Self {
        Self {
            geometry: PcurveGeometry::Elliptical {
                center,
                u_radius,
                v_radius,
                rotation,
                start_angle,
                end_angle,
            },
            t_start,
            t_end,
        }
    }

    /// Build a sinusoidal pcurve. See
    /// [`PcurveGeometry::Sinusoidal`] for the exact form.
    pub fn sinusoidal(
        u_offset: f64,
        v_const: f64,
        v_cos: f64,
        v_sin: f64,
        t_start: f64,
        t_end: f64,
    ) -> Self {
        Self {
            geometry: PcurveGeometry::Sinusoidal {
                u_offset,
                v_const,
                v_cos,
                v_sin,
            },
            t_start,
            t_end,
        }
    }

    /// Evaluate the pcurve at parameter `t`. Returns
    /// `Point3 { x: u, y: v, z: 0 }`.
    pub fn point_at(&self, t: f64) -> Point3 {
        match &self.geometry {
            PcurveGeometry::Linear { start, end } => {
                let s = self.normalised_t(t);
                let u = start.0 + s * (end.0 - start.0);
                let v = start.1 + s * (end.1 - start.1);
                Point3::new(u, v, 0.0)
            }
            PcurveGeometry::Elliptical {
                center,
                u_radius,
                v_radius,
                rotation,
                start_angle,
                end_angle,
            } => {
                let s = self.normalised_t(t);
                let alpha = start_angle + s * (end_angle - start_angle);
                let local_u = u_radius * alpha.cos();
                let local_v = v_radius * alpha.sin();
                let cr = rotation.cos();
                let sr = rotation.sin();
                let u = center.0 + local_u * cr - local_v * sr;
                let v = center.1 + local_u * sr + local_v * cr;
                Point3::new(u, v, 0.0)
            }
            PcurveGeometry::Sinusoidal {
                u_offset,
                v_const,
                v_cos,
                v_sin,
            } => {
                let u = u_offset + t;
                let v = v_const + v_cos * t.cos() + v_sin * t.sin();
                Point3::new(u, v, 0.0)
            }
            PcurveGeometry::Nurbs(n) => {
                use crate::geometry::curve::CurveEval;
                n.point_at(t)
            }
            PcurveGeometry::CylinderArcViaRationalQuad {
                theta_a,
                theta_b,
                v_const,
            } => {
                let s = self.normalised_t(t);
                let alpha = 0.5 * (theta_b - theta_a);
                let theta_m = 0.5 * (theta_a + theta_b);
                let s_alpha = alpha.sin();
                let c_alpha = alpha.cos();
                let denom_factor = 2.0 * s * (1.0 - s);
                let sin_off = s_alpha * (2.0 * s - 1.0);
                let cos_off = c_alpha * (1.0 - denom_factor) + denom_factor;
                let u = theta_m + sin_off.atan2(cos_off);
                Point3::new(u, *v_const, 0.0)
            }
        }
    }

    /// First-order tangent at parameter `t`. Returns the analytic
    /// derivative for `Linear`, `Elliptical`, and `Sinusoidal`;
    /// defers to NURBS de Boor otherwise.
    pub fn tangent_at(&self, t: f64) -> Vec3 {
        match &self.geometry {
            PcurveGeometry::Linear { start, end } => {
                let span = self.t_end - self.t_start;
                let scale = if span.abs() > 0.0 { 1.0 / span } else { 0.0 };
                Vec3::new((end.0 - start.0) * scale, (end.1 - start.1) * scale, 0.0)
            }
            PcurveGeometry::Elliptical {
                u_radius,
                v_radius,
                rotation,
                start_angle,
                end_angle,
                ..
            } => {
                let s = self.normalised_t(t);
                let alpha = start_angle + s * (end_angle - start_angle);
                let span = self.t_end - self.t_start;
                let dalpha_dt = if span.abs() > 0.0 {
                    (end_angle - start_angle) / span
                } else {
                    0.0
                };
                let dlocal_u = -u_radius * alpha.sin() * dalpha_dt;
                let dlocal_v = v_radius * alpha.cos() * dalpha_dt;
                let cr = rotation.cos();
                let sr = rotation.sin();
                Vec3::new(
                    dlocal_u * cr - dlocal_v * sr,
                    dlocal_u * sr + dlocal_v * cr,
                    0.0,
                )
            }
            PcurveGeometry::Sinusoidal { v_cos, v_sin, .. } => {
                let du_dt = 1.0;
                let dv_dt = -v_cos * t.sin() + v_sin * t.cos();
                Vec3::new(du_dt, dv_dt, 0.0)
            }
            PcurveGeometry::Nurbs(n) => {
                use crate::geometry::curve::CurveEval;
                n.tangent_at(t)
            }
            PcurveGeometry::CylinderArcViaRationalQuad {
                theta_a, theta_b, ..
            } => {
                // dθ/dt for the symmetric rational quadratic
                // closed-form, derived from the atan2 expression
                // and `arctan` derivative chain rule.
                //
                //   N(s) = sin α · (2s − 1)
                //   D(s) = cos α · (1 − 2s(1−s)) + 2s(1−s)
                //   θ(s) = θ_m + atan2(N, D)
                //   dθ/dt = (D·dN/dt − N·dD/dt) / (N² + D²)
                let span = self.t_end - self.t_start;
                let ds_dt = if span.abs() > 0.0 { 1.0 / span } else { 0.0 };
                let s = self.normalised_t(t);
                let alpha = 0.5 * (theta_b - theta_a);
                let s_alpha = alpha.sin();
                let c_alpha = alpha.cos();
                let n_t = s_alpha * (2.0 * s - 1.0);
                let denom_factor = 2.0 * s * (1.0 - s);
                let d_t = c_alpha * (1.0 - denom_factor) + denom_factor;
                let dn_dt = 2.0 * s_alpha;
                // d/ds[2s(1−s)] = 2(1−2s)
                let d_denom_factor_dt = 2.0 * (1.0 - 2.0 * s);
                let dd_dt = -c_alpha * d_denom_factor_dt + d_denom_factor_dt;
                let mag_sq = n_t * n_t + d_t * d_t;
                let dtheta_dt = if mag_sq > 0.0 {
                    ((d_t * dn_dt - n_t * dd_dt) / mag_sq) * ds_dt
                } else {
                    0.0
                };
                Vec3::new(dtheta_dt, 0.0, 0.0)
            }
        }
    }

    fn normalised_t(&self, t: f64) -> f64 {
        let span = self.t_end - self.t_start;
        if span.abs() > 0.0 {
            (t - self.t_start) / span
        } else {
            0.0
        }
    }

    /// Parameter range `[t_start, t_end]`.
    pub fn parameter_range(&self) -> (f64, f64) {
        (self.t_start, self.t_end)
    }

    /// Borrow the inner `NurbsCurve` if the pcurve is in NURBS
    /// form. Returns `None` for analytic variants — callers wanting
    /// uniform NURBS access should call [`Pcurve::to_nurbs`].
    pub fn as_nurbs(&self) -> Option<&NurbsCurve> {
        match &self.geometry {
            PcurveGeometry::Nurbs(n) => Some(n),
            _ => None,
        }
    }

    /// Mutable borrow of the inner `NurbsCurve` if the pcurve is
    /// in NURBS form. Used by band-aid passes
    /// (`gregory_pcurve_clamp`, `pcurve_refit`) that operate on CPs
    /// directly. Analytic variants return `None` because they're
    /// inherently within the surface UV box (G4) and don't need
    /// post-hoc clamping.
    pub fn as_nurbs_mut(&mut self) -> Option<&mut NurbsCurve> {
        match &mut self.geometry {
            PcurveGeometry::Nurbs(n) => Some(n),
            _ => None,
        }
    }

    /// Convert to a `NurbsCurve` in all cases. For analytic
    /// variants, builds the equivalent NURBS representation. Used
    /// by STEP export and any consumer that requires uniform NURBS
    /// data.
    pub fn to_nurbs(&self) -> NurbsCurve {
        match &self.geometry {
            PcurveGeometry::Nurbs(n) => {
                use crate::geometry::curve::CurveEval;
                let (full_start, full_end) = n.parameter_range();
                let (slice_start, slice_end, reversed) = if self.t_start <= self.t_end {
                    (self.t_start, self.t_end, false)
                } else {
                    (self.t_end, self.t_start, true)
                };
                let eps = 1.0e-12 * (full_end - full_start).abs().max(1.0);
                let restricted = if (slice_start - full_start).abs() <= eps
                    && (slice_end - full_end).abs() <= eps
                {
                    n.clone()
                } else {
                    pcurve_span::slice_nurbs_curve(n, slice_start, slice_end)
                        .unwrap_or_else(|| pcurve_nurbs_export::sampled_to_nurbs_approx(self))
                };
                if reversed {
                    pcurve_span::reverse_nurbs_curve(&restricted)
                } else {
                    restricted
                }
            }
            PcurveGeometry::Linear { start, end } => NurbsCurve {
                degree: 1,
                knots: vec![self.t_start, self.t_start, self.t_end, self.t_end],
                control_points: vec![
                    Point3::new(start.0, start.1, 0.0),
                    Point3::new(end.0, end.1, 0.0),
                ],
                weights: vec![1.0, 1.0],
            },
            PcurveGeometry::Elliptical { .. } => {
                // Sample-and-fit fallback for STEP export when an
                // analytic ellipse pcurve has to be serialised as
                // NURBS. A degree-2 rational arc is ~ulp-exact for
                // sweeps under π/2; for full circles, a 9-CP
                // standard rational construction. For the general
                // case we use a degree-3 LSQ fit over 17 samples
                // (residual < 1e-12 over a half-circle).
                pcurve_nurbs_export::sampled_to_nurbs_approx(self)
            }
            PcurveGeometry::Sinusoidal { .. } => pcurve_nurbs_export::sampled_to_nurbs_approx(self),
            PcurveGeometry::CylinderArcViaRationalQuad { .. } => {
                pcurve_nurbs_export::cylinder_arc_via_rational_quad_to_nurbs(self)
            }
        }
    }

    /// True if the pcurve has an analytic form (not a NURBS
    /// fallback). Diagnostics and provenance reporting read this.
    pub fn is_analytic(&self) -> bool {
        !matches!(&self.geometry, PcurveGeometry::Nurbs(_))
    }
}

impl From<NurbsCurve> for Pcurve {
    fn from(c: NurbsCurve) -> Self {
        Self::from_nurbs(c)
    }
}

fn nurbs_param_range(curve: &NurbsCurve) -> (f64, f64) {
    use crate::geometry::curve::CurveEval;
    curve.parameter_range()
}

#[cfg(test)]
#[path = "pcurve_tests.rs"]
mod tests;
