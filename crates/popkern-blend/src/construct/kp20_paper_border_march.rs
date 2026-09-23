//! The `OrientedSegment` ring-march arsenal: χ̃-calibrated march
//! vectors, the cross-field sampler with its pointwise-curvature
//! normal component (PK-515), and the per-surface ring-march
//! snapping. Split from `kp20_paper_border_extract` for the
//! file-size budget; the measurement semantics live on the items.

use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{Point3, Vec3};

use super::{
    ChiEndCalibration, ChiPieceCalibration, CornerJet, CrossFieldSamples, OrientedSegment,
};

impl OrientedSegment {
    /// Pointwise normal curvature of the owning surface at `p` in the
    /// tangent direction `dir` (second fundamental form over first),
    /// signed against the caller's `reference` normal (the chart's
    /// own orientation is arbitrary). Finite differences in the uv
    /// chart with a two-step Richardson guard; `None` where the chart
    /// is degenerate or the estimates disagree.
    pub(super) fn normal_curvature(&self, p: Point3, dir: Vec3, reference: Vec3) -> Option<f64> {
        use popkern_core::geometry::surface::SurfaceEval;
        use popkern_core::geometry::uv_projection::project_point_to_uv;
        let (u, v) = project_point_to_uv(&self.surface, &p)?;
        let s: &dyn SurfaceEval = &self.surface;
        let h = 1e-6;
        let su = (s.point_at(u + h, v) - s.point_at(u - h, v)) * (0.5 / h);
        let sv = (s.point_at(u, v + h) - s.point_at(u, v - h)) * (0.5 / h);
        let mut n = su.cross(&sv).normalized()?;
        if n.dot(&reference) < 0.0 {
            n = n * -1.0;
        }
        let (e, f_, g) = (su.dot(&su), su.dot(&sv), sv.dot(&sv));
        let (r1, r2) = (dir.dot(&su), dir.dot(&sv));
        let det = e * g - f_ * f_;
        if det.abs() < 1e-30 {
            return None;
        }
        let (a, b) = ((r1 * g - r2 * f_) / det, (r2 * e - r1 * f_) / det);
        let estimate = |step: f64| -> Option<f64> {
            let p0 = s.point_at(u, v);
            let p1 = s.point_at(u + a * step, v + b * step);
            let p2 = s.point_at(u - a * step, v - b * step);
            let second = (p1 - p0) + (p2 - p0);
            let first = (p1 - p2).length() * 0.5;
            (first >= 1e-30).then(|| second.dot(&n) / (first * first))
        };
        let k1 = estimate(1e-4)?;
        let k2 = estimate(5e-5)?;
        ((k1 - k2).abs() <= 0.1 * (1.0 + k2.abs())).then(|| (4.0 * k2 - k1) / 3.0)
    }

    /// Rim row: cubic through four curve samples at t = 0, 1/3, 2/3,
    /// 1 (endpoint-tangent Hermite is unreliable on the boundary
    /// NURBS, whose end derivatives are not chord-scaled).
    pub(super) fn rim_row(&self) -> [Point3; 4] {
        let p0 = self.point(0.0);
        let q1 = self.point(1.0 / 3.0);
        let q2 = self.point(2.0 / 3.0);
        let p3 = self.point(1.0);
        let interior = |a: f64, b: f64, c: f64, d: f64| {
            Point3::new(
                (a * p0.x + b * q1.x + c * q2.x + d * p3.x) / 6.0,
                (a * p0.y + b * q1.y + c * q2.y + d * p3.y) / 6.0,
                (a * p0.z + b * q1.z + c * q2.z + d * p3.z) / 6.0,
            )
        };
        [
            p0,
            interior(-5.0, 18.0, -9.0, 2.0),
            interior(2.0, -9.0, 18.0, -5.0),
            p3,
        ]
    }

    /// χ̃-calibrated march vector at one end of the rim row (`col` 0 or
    /// 6): the rim tangent rotated by the characteristic signed angle
    /// in this surface's tangent plane and scaled by the
    /// characteristic cross/rim ratio (measured on cubic χ̃ rows, so
    /// the rim leg is the CUBIC-EQUIVALENT `2·(b₁ − b₀)` of the
    /// degree-6 row).
    pub(super) fn calibrated_march(
        &self,
        row0: &[Point3; 7],
        col: usize,
        cal: &ChiEndCalibration,
        vertex: Point3,
    ) -> Vec3 {
        let rim = if col == 0 {
            (row0[1] - row0[0]) * 2.0
        } else {
            (row0[6] - row0[5]) * 2.0
        };
        let rim_len = rim.length().max(1.0e-14);
        let t_hat = rim * (1.0 / rim_len);
        let p = row0[col];
        let normal = self.surface_normal(p).unwrap_or(Vec3::new(0.0, 0.0, 1.0));
        // Orient the frame so `normal × t_hat` points away from the
        // blend vertex (into the ring), matching the χ̃ frame's
        // outward side.
        let mut side = normal.cross(&t_hat);
        if side.dot(&(p - vertex)) < 0.0 {
            side = -side;
        }
        let side = side.normalized().unwrap_or(side);
        (t_hat * cal.angle.cos() + side * cal.angle.sin()) * (cal.ratio * rim_len)
    }

    /// Sample the cap-side cross-jet and cross-curvature fields along
    /// the rim row (PK-398). The march vector blends `march_ends`
    /// linearly, kept in the tangent plane at each true curve point
    /// `S(t)`. Marching from the row's control points — the earlier
    /// construction — started every interior march off the surface by
    /// the control polygon's sagitta (1.4e-3 at r = 0.1), which the
    /// second difference read as cross-curvature (−1.4e-3/ℓ² ≈ −2.0
    /// against a straight sweep's 0) and served as the cap↔strip
    /// transverse G2 class (1.4 normalized on every rim+verticals seam).
    pub(super) fn march_cross_field(
        &self,
        row0: &[Point3; 7],
        march_ends: [Vec3; 2],
        first_order_marches: &mut usize,
    ) -> CrossFieldSamples {
        const FIELD_SAMPLES: usize = 25;
        let mut jets = Vec::with_capacity(FIELD_SAMPLES);
        let mut curvatures = Vec::with_capacity(FIELD_SAMPLES);
        for j in 0..FIELD_SAMPLES {
            let t = j as f64 / (FIELD_SAMPLES - 1) as f64;
            let rim_point = popkern_core::geometry::gregory::eval_sextic_row(row0, t);
            // March from the rim point's SURFACE foot: at a gapped
            // corner (trim curves that miss each other — 3.4e-3 on the
            // blunt pentagon G2 sweeps) the rim row interpolates a
            // corner 1.6e-3 off the strip, and a march started there
            // reads the offset as cross-curvature (offset/ℓ² ≈ 16)
            // exactly like the control-point march did. The fields are
            // relative to the row, so the rim reconcile's later landing
            // of the row carries them along.
            let q0 = self.surface_foot(rim_point).unwrap_or(rim_point);
            let n_hat = self
                .surface_normal(q0)
                .and_then(|n| n.normalized())
                .unwrap_or(Vec3::new(0.0, 0.0, 0.0));
            // Smootherstep blend: zero slope AND curvature at both
            // ends, so the corner block's third-order slots agree
            // between the two pieces by construction. The χ̃ scale
            // convention makes the cross leg ~2× longer at the spoke
            // point m than the corner-compatible leg at c, so a linear
            // blend leaves a length KINK at m (both halves rising
            // toward it) and a twist step at c; the field must be C²
            // across m and twist-free at c for the shared junction
            // pins to agree with the samples (measured: linear blend
            // forced the constrained fit 1.1e-2 off its samples).
            let w = t * t * t * (10.0 + t * (-15.0 + 6.0 * t));
            let blend = march_ends[0] * (1.0 - w) + march_ends[1] * w;
            let in_plane = blend - n_hat * blend.dot(&n_hat);
            let v = in_plane
                .normalized()
                .map(|u| u * blend.length())
                .unwrap_or(blend);
            let (q1, q2) = self.ring_march_vector(q0, v, first_order_marches);
            // Cap-side C² prolongation. The WP1b point reflection
            // r1 = 2q0 − q1 mirrors a SURFACE CHORD through q0: the
            // chord leaves the strip's tangent plane by ~κℓ²/2, so
            // the cap's cross-jet tilted out of the surface by
            // ~κℓ/2 radians — measured as a UNIFORM 2.6°–6.9° rim
            // tangency error along every setback arc (PK-212). G1
            // across the rim requires the jet IN the tangent plane:
            // project the chord's normal component out before
            // reflecting.
            let chord = q1 - q0;
            let jet = -(chord - n_hat * chord.dot(&n_hat));
            // Curvature semantics of a Bernstein depth row: derivatives
            // at the rim are 3·(r1 − r0) and 6·(r2 − 2r1 + r0). With
            // the leg equal to the sample spacing ℓ the row represents
            // the surface at depth parameter α = 3ℓ, so its second
            // difference is α²f''/6 = 1.5·ℓ²f'' — 1.5× the SAMPLE
            // second difference (≈ ℓ²f''). The in-plane part is
            // parameterization (row shaping) and keeps that model.
            let second = (q2 - q1) - (q1 - q0);
            let mut sample = second * 1.5;
            // The NORMAL part is geometry: the cap's border normal
            // curvature in the march direction is n·d2·30/|6·d1|², so
            // matching the strip pointwise needs n·d2 = 1.2·κ·|jet|².
            // The finite march's average (1.5·n·second ≈ 1.5·ℓ²·κ_avg
            // over [0, 2ℓ]) both overshoots the Bernstein factor by
            // 25% and, where κ varies at the leg scale — every G2
            // section's contact end, where κ climbs from 0 — reads
            // the wrong station entirely: measured 3.5e-4 at a corner
            // whose pointwise jet (and shared corner block) is 0, a
            // transverse rim curvature error of ~1 against the 5e-2
            // commercial gate (PK-515). Pointwise II is exact for
            // both. Where the curvature solve declines (degenerate
            // chart), the march average stays as the fallback.
            if let (Some(kappa), true) = (
                self.normal_curvature(q0, jet * -1.0, n_hat),
                n_hat.length() > 0.5,
            ) {
                let in_plane_part = sample - n_hat * sample.dot(&n_hat);
                sample = in_plane_part + n_hat * (1.2 * kappa * jet.dot(&jet));
            }
            jets.push((t, jet));
            curvatures.push((t, sample));
        }
        CrossFieldSamples { jets, curvatures }
    }

    /// Cross-field samples with χ̃-calibrated ends, or corner-compatible
    /// jets at `c` where `corner_jets` are given (the cap-side jet is
    /// the negated march vector).
    pub(super) fn border_piece_calibrated(
        &self,
        row0: &[Point3; 7],
        calibration: &ChiPieceCalibration,
        vertex: Point3,
        first_order_marches: &mut usize,
        corner_jets: [Option<CornerJet>; 2],
    ) -> CrossFieldSamples {
        let start = match corner_jets[0] {
            Some(corner) => -corner.d1,
            None => self.calibrated_march(row0, 0, &calibration.start, vertex),
        };
        let end = match corner_jets[1] {
            Some(corner) => -corner.d1,
            None => self.calibrated_march(row0, 6, &calibration.end, vertex),
        };
        self.march_cross_field(row0, [start, end], first_order_marches)
    }

    /// March into the ring by the leg vector `v` on the owning
    /// surface: returns surface points near `q0 + v` and `q0 + 2v`
    /// (exact in-plane; radially projected onto cylinders;
    /// closest-point-projected onto NURBS; tangent-plane first-order
    /// otherwise).
    pub(super) fn ring_march_vector(
        &self,
        q0: Point3,
        v: Vec3,
        first_order_marches: &mut usize,
    ) -> (Point3, Point3) {
        let raw1 = q0 + v;
        let raw2 = q0 + v * 2.0;
        match &self.surface {
            SurfaceGeometry::Plane { origin, normal, .. } => {
                let snap = |p: Point3| {
                    let d = (p - *origin).dot(normal);
                    p - *normal * d
                };
                (snap(raw1), snap(raw2))
            }
            SurfaceGeometry::Cylinder {
                origin,
                axis,
                radius,
                ..
            } => {
                let snap = |p: Point3| {
                    let along = (p - *origin).dot(axis);
                    let foot = *origin + *axis * along;
                    let radial = p - foot;
                    let unit = radial.normalized().unwrap_or(radial);
                    foot + unit * *radius
                };
                (snap(raw1), snap(raw2))
            }
            SurfaceGeometry::Nurbs(nurbs) => {
                // Snap the march onto the NURBS sweep via closest-point
                // projection (second-order accurate ring sampling).
                use popkern_core::geometry::surface::SurfaceEval;
                use popkern_core::geometry::uv_projection::project_point_to_uv;
                let snap = |p: Point3| -> Option<Point3> {
                    let (u, v) = project_point_to_uv(&self.surface, &p)?;
                    if popkern_core::trace::trace_enabled("kp-march") {
                        let out = nurbs.point_at(u, v);
                        eprintln!(
                            "      nurbs snap: raw=({:.5},{:.5},{:.5}) uv=({:.4},{:.4}) out=({:.5},{:.5},{:.5})",
                            p.x, p.y, p.z, u, v, out.x, out.y, out.z,
                        );
                    }
                    Some(nurbs.point_at(u, v))
                };
                match (snap(raw1), snap(raw2)) {
                    (Some(q1), Some(q2)) => (q1, q2),
                    _ => {
                        *first_order_marches += 1;
                        (raw1, raw2)
                    }
                }
            }
            _ => {
                *first_order_marches += 1;
                (raw1, raw2)
            }
        }
    }
}
