//! PK-090: corner-compatible two-piece KP20 border from a blend
//! boundary's 2n segments.
//!
//! The apex blend boundary alternates fillet-end arcs and face-edge
//! segments — two genuine surface pieces per KP20 sector, meeting at
//! the sector corner `c`. This extractor builds the paper cap's
//! [`Kp20PaperGuideBorder`] directly from those segments, with depth
//! rows sampled from the OWNING surfaces (exact ring marches for the
//! analytic cylinder/plane cases; first-order tangent-plane march
//! otherwise, recorded in the report). This replaces the invalid
//! smooth-cubic midpoint split, whose halves are not corner-compatible
//! at `c` and kink the assembled rim.
//!
//! Layout: sector `k` = (`FilletEnd[k]`, `FaceEdge[k]`); `c` is the
//! junction where a fillet end meets its trailing face segment, and
//! the spoke points `m` are the face→fillet junctions. Depth rows are
//! the cap-side C² prolongation of the ring (the WP1b convention):
//! `row1 = 2 q0 − q1`, `row2 = 4 q0 − 4 q1 + q2` with `q_d` true ring
//! samples at depth `d`.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::gregory::{Kp20PaperGuideBorder, Kp20PaperGuideBorderPiece};
use popkern_core::geometry::{Point3, SurfaceGeometry, Vec3};
use popkern_core::topology::BRep;

use crate::blend::{BlendBoundary, BoundarySegment};
use crate::fillet_surface::FilletSurfaceResult;

/// Extraction evidence for diagnostics and gating.
#[derive(Debug, Clone)]
pub struct Kp20PaperBorderExtractReport {
    /// Valence (number of sectors).
    pub n: usize,
    /// Uniform ring-march depth used for the cross rows.
    pub depth: f64,
    /// Worst gap between consecutive segment endpoints along the loop.
    pub max_loop_gap: f64,
    /// Segments whose owning surface had no exact ring march and fell
    /// back to a first-order tangent-plane march.
    pub first_order_marches: usize,
    /// `true` when the loop order was reversed to make the boundary
    /// wind counterclockwise about its mean normal.
    pub reversed_for_orientation: bool,
}

/// Build the paper cap's two-piece border from an apex blend
/// boundary.
pub fn kp20_paper_border_from_blend_boundary(
    boundary: &BlendBoundary,
    fillet_surfaces: &[FilletSurfaceResult],
    brep: &BRep,
) -> Result<(Kp20PaperGuideBorder, Kp20PaperBorderExtractReport), String> {
    let segment_count = boundary.segments.len();
    if segment_count < 6 || !segment_count.is_multiple_of(2) {
        return Err(format!(
            "paper border needs an alternating 2n-segment boundary with n >= 3, got {segment_count} segments"
        ));
    }
    let n = segment_count / 2;

    // Orient every segment along the loop: segment k must start where
    // segment k-1 ends.
    let mut curves: Vec<OrientedSegment> = Vec::with_capacity(segment_count);
    for segment in &boundary.segments {
        curves.push(OrientedSegment::new(segment, fillet_surfaces, brep)?);
    }
    let mut max_loop_gap = 0.0_f64;
    for k in 0..segment_count {
        let previous_end = curves[(k + segment_count - 1) % segment_count].end();
        let (gap_fwd, gap_rev) = {
            let this = &curves[k];
            (
                previous_end.distance_to(&this.start()),
                previous_end.distance_to(&this.end()),
            )
        };
        if gap_rev < gap_fwd {
            curves[k].reversed = !curves[k].reversed;
            max_loop_gap = max_loop_gap.max(gap_rev);
        } else {
            max_loop_gap = max_loop_gap.max(gap_fwd);
        }
    }

    // Normalize the phase so segment 0 is a FilletEnd (sector k =
    // FilletEnd[k] then FaceEdge[k], with c at their junction).
    if matches!(boundary.segments[0], BoundarySegment::FaceEdge { .. }) {
        curves.rotate_left(1);
    }

    // Wind counterclockwise about the boundary's mean normal so the
    // sector order matches the cap chart's handedness.
    let centroid = loop_centroid(&curves);
    let mean_normal = loop_mean_normal(&curves, centroid);
    let winding = loop_winding(&curves, centroid, mean_normal);
    let reversed_for_orientation = winding < 0.0;
    if reversed_for_orientation {
        curves.reverse();
        for segment in curves.iter_mut() {
            segment.reversed = !segment.reversed;
        }
        // Keep the FilletEnd-first phase after reversal.
        if !curves[0].is_fillet_end {
            curves.rotate_left(1);
        }
    }

    // Uniform ring-march depth from the mean segment chord.
    let mean_chord = curves
        .iter()
        .map(|segment| segment.start().distance_to(&segment.end()))
        .sum::<f64>()
        / segment_count as f64;
    let depth = 0.8 * mean_chord;

    // Boundary conditioning to the g∘χ̃ shape class via macro-edge
    // fitting. The characteristic layout corners ONLY at the n sector
    // corners c and is C²-smooth at the n spoke points m; a raw
    // fillet boundary corners at all 2n junctions, and a single bi-3
    // span lacks the freedom to be C² at m without disturbing its c
    // corner (measured: second-difference unification flattens c
    // from 119° to 28°). Each macro-edge (c_k → c_{k+1}, through m)
    // is therefore refit as a 2-span cubic spline with:
    // - both c endpoints and the m junction interpolated exactly,
    // - the TRUE segment tangent directions preserved at both c ends
    //   (the real ~119° corner, matching χ̃'s corner class),
    // - C¹/C² at m by construction (the constraints eliminate all
    //   but two scalar unknowns — the c-end tangent magnitudes),
    // - the remaining two scalars least-squares fit to dense samples
    //   of the true boundary segments.
    // Junction j sits at the start of segment j; with the
    // FilletEnd-first phase, ODD junctions are the sector corners c.
    let mut rims: Vec<[Point3; 4]> = curves.iter().map(|segment| segment.rim_row()).collect();
    for k in 0..n {
        let j = 2 * k + 1; // c junction between segment j-1 (arc) and j (bridge)
        let bridge_index = j % segment_count;
        let arc_index = (j + 1) % segment_count;
        let c_start = rims[bridge_index][0];
        let m_point = rims[arc_index][0];
        let c_end = rims[arc_index][3];
        // True corner directions from the raw interpolated rims.
        let t0 = (rims[bridge_index][1] - rims[bridge_index][0])
            .normalized()
            .unwrap_or(Vec3::new(1.0, 0.0, 0.0));
        let t1 = (rims[arc_index][3] - rims[arc_index][2])
            .normalized()
            .unwrap_or(Vec3::new(1.0, 0.0, 0.0));
        // Solve min over (alpha, beta) of the sample distance, with
        // P1 = c_start + alpha t0, Q2 = c_end - beta t1,
        // P2 = (c_start + 4m - c_end + alpha t0 + beta t1)/4,
        // Q1 = 2m - P2 (the C1/C2-at-m eliminations).
        const SAMPLES_PER_SEGMENT: usize = 9;
        let mut ata = [[0.0_f64; 2]; 2];
        let mut atb = [0.0_f64; 2];
        for (span, segment_index) in [(0_usize, bridge_index), (1, arc_index)] {
            for sample in 0..SAMPLES_PER_SEGMENT {
                let t = (sample + 1) as f64 / (SAMPLES_PER_SEGMENT + 1) as f64;
                let target = curves[segment_index].point(t);
                // Cubic Bernstein weights on this span.
                let sfp = 1.0 - t;
                let w = [
                    sfp * sfp * sfp,
                    3.0 * sfp * sfp * t,
                    3.0 * sfp * t * t,
                    t * t * t,
                ];
                // Express the span point as constant + alpha·a + beta·b.
                // Span 0 (P): [c_start, P1, P2, m]; span 1 (Q): [m, Q1, Q2, c_end].
                let combine = |points: [(Point3, f64); 3]| -> [f64; 3] {
                    let mut acc = [0.0_f64; 3];
                    for (point, weight) in points {
                        acc[0] += point.x * weight;
                        acc[1] += point.y * weight;
                        acc[2] += point.z * weight;
                    }
                    acc
                };
                let p2_const = Point3::new(
                    (c_start.x + 4.0 * m_point.x - c_end.x) / 4.0,
                    (c_start.y + 4.0 * m_point.y - c_end.y) / 4.0,
                    (c_start.z + 4.0 * m_point.z - c_end.z) / 4.0,
                );
                let (constant, a_vec, b_vec) = if span == 0 {
                    let constant =
                        combine([(c_start, w[0] + w[1]), (p2_const, w[2]), (m_point, w[3])]);
                    (constant, t0 * (w[1] + w[2] / 4.0), t1 * (w[2] / 4.0))
                } else {
                    let q1_const = Point3::new(
                        2.0 * m_point.x - p2_const.x,
                        2.0 * m_point.y - p2_const.y,
                        2.0 * m_point.z - p2_const.z,
                    );
                    let constant =
                        combine([(m_point, w[0]), (q1_const, w[1]), (c_end, w[2] + w[3])]);
                    (constant, t0 * (-w[1] / 4.0), t1 * (-w[1] / 4.0 - w[2]))
                };
                let r = Vec3::new(
                    target.x - constant[0],
                    target.y - constant[1],
                    target.z - constant[2],
                );
                ata[0][0] += a_vec.dot(&a_vec);
                ata[0][1] += a_vec.dot(&b_vec);
                ata[1][0] += a_vec.dot(&b_vec);
                ata[1][1] += b_vec.dot(&b_vec);
                atb[0] += a_vec.dot(&r);
                atb[1] += b_vec.dot(&r);
            }
        }
        let det = ata[0][0] * ata[1][1] - ata[0][1] * ata[1][0];
        if det.abs() < 1.0e-18 {
            continue;
        }
        let alpha = (atb[0] * ata[1][1] - atb[1] * ata[0][1]) / det;
        let beta = (ata[0][0] * atb[1] - ata[1][0] * atb[0]) / det;
        let p1 = c_start + t0 * alpha;
        let q2 = c_end - t1 * beta;
        let p2 = Point3::new(
            (c_start.x + 4.0 * m_point.x - c_end.x + alpha * t0.x + beta * t1.x) / 4.0,
            (c_start.y + 4.0 * m_point.y - c_end.y + alpha * t0.y + beta * t1.y) / 4.0,
            (c_start.z + 4.0 * m_point.z - c_end.z + alpha * t0.z + beta * t1.z) / 4.0,
        );
        let q1 = Point3::new(
            2.0 * m_point.x - p2.x,
            2.0 * m_point.y - p2.y,
            2.0 * m_point.z - p2.z,
        );
        rims[bridge_index] = [c_start, p1, p2, m_point];
        rims[arc_index] = [m_point, q1, q2, c_end];
    }

    // Pass 2: cross rows calibrated to the characteristic border's
    // own parameterization structure. For each of the four piece-end
    // sites (piece A/B × start/end), the χ̃ border for this valence
    // supplies the cross-leg's length ratio to the rim leg and its
    // SIGNED in-plane angle from the rim tangent (measured: ~90° at
    // m with ratio ~1.08; ~94°/86° mirror-skewed at c with ratio
    // 1.0 at n = 5). The real border's cross rows reproduce those
    // ratios and angles in each surface's tangent plane — the
    // ρ-reparameterization freedom applied so the input reaches the
    // assembly in exactly the parameterization the τ/b machinery was
    // built around. Marches are projected onto the owning surfaces.
    let calibration = chi_march_calibration(n)
        .ok_or_else(|| "characteristic march calibration failed".to_string())?;
    let mut first_order_marches = 0_usize;
    let mut sectors = Vec::with_capacity(n);
    for k in 0..n {
        let piece_a = curves[2 * k].border_piece_calibrated(
            &rims[2 * k],
            &calibration[0],
            boundary.vertex_position,
            &mut first_order_marches,
        );
        let piece_b = curves[2 * k + 1].border_piece_calibrated(
            &rims[2 * k + 1],
            &calibration[1],
            boundary.vertex_position,
            &mut first_order_marches,
        );
        // Sector rows run [m^k → c] then [c → m^{k+1}]: the FilletEnd
        // piece starts the sector at m and ends at c.
        sectors.push([
            Kp20PaperGuideBorderPiece { rows: piece_a },
            Kp20PaperGuideBorderPiece { rows: piece_b },
        ]);
    }

    Ok((
        Kp20PaperGuideBorder { sectors },
        Kp20PaperBorderExtractReport {
            n,
            depth,
            max_loop_gap,
            first_order_marches,
            reversed_for_orientation,
        },
    ))
}

/// Cross-leg calibration at one piece end: length ratio to the rim
/// leg and signed in-plane angle from the rim tangent.
#[derive(Debug, Clone, Copy)]
struct ChiEndCalibration {
    ratio: f64,
    angle: f64,
}

/// Per-piece calibration (start and end).
#[derive(Debug, Clone, Copy)]
struct ChiPieceCalibration {
    start: ChiEndCalibration,
    end: ChiEndCalibration,
}

/// Measure the characteristic border's cross-leg structure for this
/// valence: `[piece A ([m→c]-like near_c), piece B ([c→m]-like
/// near_m)]`, each with start/end calibrations.
fn chi_march_calibration(n: usize) -> Option<[ChiPieceCalibration; 2]> {
    use popkern_core::geometry::gregory::cc_lifted_ring_border_pieces;
    let pieces = cc_lifted_ring_border_pieces(n, |_, _| 0.0).ok()?;
    let [near_c, near_m] = &pieces[0];
    let plane_normal = Vec3::new(0.0, 0.0, 1.0);
    let measure = |rows: &[[Point3; 4]; 3], col: usize| -> Option<ChiEndCalibration> {
        let rim = if col == 0 {
            rows[0][1] - rows[0][0]
        } else {
            rows[0][3] - rows[0][2]
        };
        let cross = rows[1][col] - rows[0][col];
        let rim_len = rim.length();
        let cross_len = cross.length();
        if rim_len < 1.0e-14 || cross_len < 1.0e-14 {
            return None;
        }
        let t_hat = rim * (1.0 / rim_len);
        let side = plane_normal.cross(&t_hat);
        let angle = cross.dot(&side).atan2(cross.dot(&t_hat));
        Some(ChiEndCalibration {
            ratio: cross_len / rim_len,
            angle,
        })
    };
    Some([
        ChiPieceCalibration {
            start: measure(near_c, 0)?,
            end: measure(near_c, 3)?,
        },
        ChiPieceCalibration {
            start: measure(near_m, 0)?,
            end: measure(near_m, 3)?,
        },
    ])
}

/// One boundary segment with loop orientation and its owning surface.
struct OrientedSegment {
    curve: popkern_core::geometry::NurbsCurve,
    surface: SurfaceGeometry,
    is_fillet_end: bool,
    reversed: bool,
}

impl OrientedSegment {
    fn new(
        segment: &BoundarySegment,
        fillet_surfaces: &[FilletSurfaceResult],
        brep: &BRep,
    ) -> Result<Self, String> {
        match segment {
            BoundarySegment::FilletEnd {
                fillet_index,
                curve,
            } => {
                let surface = fillet_surfaces
                    .get(*fillet_index)
                    .ok_or_else(|| format!("fillet index {fillet_index} out of range"))?
                    .surface
                    .clone();
                Ok(Self {
                    curve: curve.clone(),
                    surface,
                    is_fillet_end: true,
                    reversed: false,
                })
            }
            BoundarySegment::FaceEdge { face, curve } => {
                let surface = brep.faces[*face]
                    .surface
                    .as_ref()
                    .ok_or_else(|| format!("face {face:?} has no surface"))?
                    .geometry
                    .clone();
                Ok(Self {
                    curve: curve.clone(),
                    surface,
                    is_fillet_end: false,
                    reversed: false,
                })
            }
        }
    }

    fn parameter(&self, t: f64) -> f64 {
        let (t0, t1) = self.curve.parameter_range();
        let s = if self.reversed { 1.0 - t } else { t };
        t0 + s * (t1 - t0)
    }

    fn point(&self, t: f64) -> Point3 {
        self.curve.point_at(self.parameter(t))
    }

    fn start(&self) -> Point3 {
        self.point(0.0)
    }

    fn end(&self) -> Point3 {
        self.point(1.0)
    }

    /// Surface normal near `p` on the owning surface (used to build
    /// the shared tangent plane at smooth spoke junctions).
    fn surface_normal(&self, p: Point3) -> Option<Vec3> {
        use popkern_core::geometry::surface::SurfaceEval;
        use popkern_core::geometry::uv_projection::project_point_to_uv;
        let uv = project_point_to_uv(&self.surface, &p)?;
        Some(self.surface.normal_at(uv.0, uv.1))
    }

    /// Rim row: cubic through four curve samples at t = 0, 1/3, 2/3,
    /// 1 (endpoint-tangent Hermite is unreliable on the boundary
    /// NURBS, whose end derivatives are not chord-scaled).
    fn rim_row(&self) -> [Point3; 4] {
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

    /// Bi-3 depth-2 border piece with χ̃-calibrated cross rows: at
    /// each end the march vector is the rim tangent rotated by the
    /// characteristic signed angle in this surface's tangent plane
    /// and scaled by the characteristic cross/rim ratio; interior
    /// columns blend the end vectors.
    fn border_piece_calibrated(
        &self,
        row0: &[Point3; 4],
        calibration: &ChiPieceCalibration,
        vertex: Point3,
        first_order_marches: &mut usize,
    ) -> [[Point3; 4]; 3] {
        let row0 = *row0;
        let end_vector = |col: usize, cal: &ChiEndCalibration| -> Vec3 {
            let rim = if col == 0 {
                row0[1] - row0[0]
            } else {
                row0[3] - row0[2]
            };
            let rim_len = rim.length().max(1.0e-14);
            let t_hat = rim * (1.0 / rim_len);
            let p = row0[col];
            let normal = self.surface_normal(p).unwrap_or(Vec3::new(0.0, 0.0, 1.0));
            // Orient the frame so `normal × t_hat` points away from
            // the blend vertex (into the ring), matching the χ̃
            // frame's outward side.
            let mut side = normal.cross(&t_hat);
            if side.dot(&(p - vertex)) < 0.0 {
                side = -side;
            }
            let side = side.normalized().unwrap_or(side);
            (t_hat * cal.angle.cos() + side * cal.angle.sin()) * (cal.ratio * rim_len)
        };
        let start_vec = end_vector(0, &calibration.start);
        let end_vec = end_vector(3, &calibration.end);
        if std::env::var_os("KP20_MARCH_TRACE").is_some() {
            eprintln!(
                "    march trace: start_vec=({:.5},{:.5},{:.5}) end_vec=({:.5},{:.5},{:.5}) fillet_end={}",
                start_vec.x,
                start_vec.y,
                start_vec.z,
                end_vec.x,
                end_vec.y,
                end_vec.z,
                self.is_fillet_end,
            );
        }
        let mut rows = [row0, row0, row0];
        for (index, &q0) in row0.iter().enumerate() {
            let blend = index as f64 / 3.0;
            let v = start_vec * (1.0 - blend) + end_vec * blend;
            let (q1, q2) = self.ring_march_vector(q0, v, first_order_marches);
            // Cap-side C² prolongation (WP1b convention):
            // r1 = 2 q0 − q1, r2 = 4 q0 − 4 q1 + q2.
            rows[1][index] = Point3::new(2.0 * q0.x - q1.x, 2.0 * q0.y - q1.y, 2.0 * q0.z - q1.z);
            rows[2][index] = Point3::new(
                4.0 * q0.x - 4.0 * q1.x + q2.x,
                4.0 * q0.y - 4.0 * q1.y + q2.y,
                4.0 * q0.z - 4.0 * q1.z + q2.z,
            );
        }
        rows
    }

    /// March into the ring by the leg vector `v` on the owning
    /// surface: returns surface points near `q0 + v` and `q0 + 2v`
    /// (exact in-plane; radially projected onto cylinders;
    /// closest-point-projected onto NURBS; tangent-plane first-order
    /// otherwise).
    fn ring_march_vector(
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
                    if std::env::var_os("KP20_MARCH_TRACE").is_some() {
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

fn loop_centroid(curves: &[OrientedSegment]) -> Point3 {
    let mut acc = Vec3::new(0.0, 0.0, 0.0);
    let mut count = 0.0;
    for segment in curves {
        let p = segment.start();
        acc = acc + Vec3::new(p.x, p.y, p.z);
        count += 1.0;
    }
    Point3::new(acc.x / count, acc.y / count, acc.z / count)
}

fn loop_mean_normal(curves: &[OrientedSegment], centroid: Point3) -> Vec3 {
    let mut acc = Vec3::new(0.0, 0.0, 0.0);
    for segment in curves {
        let a = segment.start() - centroid;
        let b = segment.end() - centroid;
        acc = acc + a.cross(&b);
    }
    acc.normalized().unwrap_or(Vec3::new(0.0, 0.0, 1.0))
}

fn loop_winding(curves: &[OrientedSegment], centroid: Point3, normal: Vec3) -> f64 {
    let mut acc = 0.0;
    for segment in curves {
        let a = segment.start() - centroid;
        let b = segment.end() - centroid;
        acc += a.cross(&b).dot(&normal);
    }
    acc
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{FilletRadius, FilletSpec, preview_fillet};
    use popkern_core::geometry::gregory::kp20_paper_border_corner_mismatch;
    use popkern_core::tolerance::Tolerance;

    /// End-to-end extraction gate on real rolling-ball data: the loop
    /// closes exactly, every ring march is exact (cylinder/plane),
    /// and the corner mismatch is finite and bounded — the measured
    /// G¹-input incompatibility at c, recorded rather than hidden.
    #[test]
    fn pyramid_apex_extraction_is_closed_and_exact() {
        let tol = Tolerance::default();
        let pyramid =
            popkern_core::topology::BRep::n_gonal_pyramid(5, 1.0, 1.2, tol).expect("pyramid");
        let apex = pyramid
            .vertices
            .iter()
            .find(|(_, v)| {
                v.position.x.abs() < 1e-9 && v.position.y.abs() < 1e-9 && v.position.z > 0.0
            })
            .map(|(h, _)| h)
            .expect("apex");
        let apex_edges: Vec<_> = pyramid
            .edges
            .iter()
            .filter_map(|(eh, _)| {
                let he = pyramid.edges[eh].half_edge;
                let origin = pyramid.half_edges[he].origin;
                let twin_origin = pyramid.half_edges[pyramid.half_edges[he].twin].origin;
                (origin == apex || twin_origin == apex).then_some(eh)
            })
            .collect();
        let specs: Vec<FilletSpec> = apex_edges
            .iter()
            .map(|&edge| FilletSpec {
                edge,
                radius: FilletRadius::Constant(0.12),
            })
            .collect();
        let preview = preview_fillet(&pyramid, &specs, tol).expect("preview");
        let fillet_surfaces: Vec<FilletSurfaceResult> = preview
            .per_edge
            .iter()
            .map(|eg| FilletSurfaceResult {
                surface: eg.blend_surface.clone(),
                trim_curve_1: eg.trim_curve_1.clone(),
                trim_curve_2: eg.trim_curve_2.clone(),
                continuity_downgrade: None,
            })
            .collect();
        let indices: Vec<usize> = (0..fillet_surfaces.len()).collect();
        let boundary =
            crate::blend::compute_blend_boundary(apex, &indices, &fillet_surfaces, &pyramid, &tol)
                .expect("boundary");
        let (border, report) =
            kp20_paper_border_from_blend_boundary(&boundary, &fillet_surfaces, &pyramid)
                .expect("extract");
        assert_eq!(report.n, 5);
        assert!(
            report.max_loop_gap < 1.0e-9,
            "loop gap {}",
            report.max_loop_gap
        );
        assert_eq!(
            report.first_order_marches, 0,
            "cylinder/plane marches must be exact"
        );
        let c = (core::f64::consts::TAU / 5.0).cos();
        let mut mismatch = 0.0_f64;
        for sector in &border.sectors {
            mismatch = mismatch.max(kp20_paper_border_corner_mismatch(
                5,
                &sector[0].rows,
                &sector[1].rows,
                2.0 * c / 3.0,
            ));
        }
        // Rolling-ball input is G1-class: the corner mismatch is real
        // and must be REPORTED, not zero. Bound it loosely so gross
        // extraction bugs (orders of magnitude off) still fail.
        assert!(
            mismatch > 1.0e-4 && mismatch < 1.0,
            "corner mismatch out of expected G1-class band: {mismatch:e}",
        );
    }
}

/// Corner-setback border: the vertex-blend boundary as n fillet-end
/// arcs set back to the TRIM-INTERSECTION corners — for each pair of
/// angularly adjacent fillets, the point where their contact curves
/// on the shared support face intersect. The boundary then lies
/// entirely on the fillet surfaces: each sector's two pieces are the
/// halves of ONE arc on ONE G² sweep (C² at the spoke point m by
/// construction), and the sector corners c sit where two sweeps meet
/// tangent to the same face (curvature-compatible corner data — the
/// property the G² sweep buys). This is the classical setback vertex
/// blend, shaped exactly like the characteristic layout.
pub fn kp20_paper_border_from_fillet_setback(
    fillet_surfaces: &[FilletSurfaceResult],
    vertex_position: Point3,
) -> Result<(Kp20PaperGuideBorder, Kp20PaperBorderExtractReport), String> {
    use popkern_core::geometry::surface::SurfaceEval;
    use popkern_core::geometry::uv_projection::project_point_to_uv;
    let n = fillet_surfaces.len();
    if n < 3 {
        return Err(format!("setback border needs >= 3 fillets, got {n}"));
    }

    // Angular order of the fillets about the vertex (best-fit frame
    // from the trim-curve midpoints).
    let rep_point = |fi: usize| -> Point3 {
        let c = &fillet_surfaces[fi].trim_curve_1;
        let (t0, t1) = c.parameter_range();
        c.point_at(0.5 * (t0 + t1))
    };
    let centroid = {
        let mut acc = [0.0_f64; 3];
        for fi in 0..n {
            let p = rep_point(fi);
            acc[0] += p.x;
            acc[1] += p.y;
            acc[2] += p.z;
        }
        Point3::new(acc[0] / n as f64, acc[1] / n as f64, acc[2] / n as f64)
    };
    let mut axis = Vec3::new(0.0, 0.0, 0.0);
    for fi in 0..n {
        let a = rep_point(fi) - centroid;
        let b = rep_point((fi + 1) % n) - centroid;
        axis = axis + a.cross(&b);
    }
    let axis = axis
        .normalized()
        .ok_or_else(|| "degenerate fillet arrangement".to_string())?;
    let seed = rep_point(0) - centroid;
    let u0 = (seed - axis * seed.dot(&axis))
        .normalized()
        .ok_or_else(|| "degenerate angular frame".to_string())?;
    let v0 = axis.cross(&u0);
    let angle_of = |p: Point3| -> f64 {
        let d = p - centroid;
        d.dot(&v0).atan2(d.dot(&u0))
    };
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&a, &b| {
        angle_of(rep_point(a))
            .partial_cmp(&angle_of(rep_point(b)))
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    // Trim-intersection corner between adjacent fillets: the closest
    // approach between any trim-curve pair (one from each side),
    // sample-refined.
    let closest_between = |ca: &popkern_core::geometry::NurbsCurve,
                           cb: &popkern_core::geometry::NurbsCurve|
     -> (f64, Point3) {
        const S: usize = 96;
        let (a0, a1) = ca.parameter_range();
        let (b0, b1) = cb.parameter_range();
        let mut best = (f64::INFINITY, Point3::new(0.0, 0.0, 0.0));
        for i in 0..=S {
            let ta = a0 + (a1 - a0) * i as f64 / S as f64;
            let pa = ca.point_at(ta);
            for j in 0..=S {
                let tb = b0 + (b1 - b0) * j as f64 / S as f64;
                let pb = cb.point_at(tb);
                let d = pa.distance_to(&pb);
                if d < best.0 {
                    best = (
                        d,
                        Point3::new(
                            0.5 * (pa.x + pb.x),
                            0.5 * (pa.y + pb.y),
                            0.5 * (pa.z + pb.z),
                        ),
                    );
                }
            }
        }
        best
    };
    let mut corners = Vec::with_capacity(n);
    for k in 0..n {
        let fa = &fillet_surfaces[order[k]];
        let fb = &fillet_surfaces[order[(k + 1) % n]];
        let mut best = (f64::INFINITY, Point3::new(0.0, 0.0, 0.0));
        for ca in [&fa.trim_curve_1, &fa.trim_curve_2] {
            for cb in [&fb.trim_curve_1, &fb.trim_curve_2] {
                let candidate = closest_between(ca, cb);
                if candidate.0 < best.0 {
                    best = candidate;
                }
            }
        }
        corners.push(best);
    }
    let max_corner_gap = corners.iter().map(|c| c.0).fold(0.0_f64, f64::max);

    // One arc per fillet between its two corners, sampled ON the
    // sweep along the uv segment; halves = the sector pieces.
    let mut sectors = Vec::with_capacity(n);
    let mut first_order_marches = 0_usize;
    let calibration = chi_march_calibration(n)
        .ok_or_else(|| "characteristic march calibration failed".to_string())?;
    for k in 0..n {
        let fillet = &fillet_surfaces[order[k]];
        let c_prev = corners[(k + n - 1) % n].1;
        let c_next = corners[k].1;
        let (ua, va) = project_point_to_uv(&fillet.surface, &c_prev)
            .ok_or_else(|| "corner projection failed".to_string())?;
        let (ub, vb) = project_point_to_uv(&fillet.surface, &c_next)
            .ok_or_else(|| "corner projection failed".to_string())?;
        let on_surface = |t: f64| -> Point3 {
            let u = ua + (ub - ua) * t;
            let v = va + (vb - va) * t;
            fillet.surface.point_at(u, v)
        };
        // Cubic through 4 samples, then de Casteljau halves — the
        // sector pieces [m→c] never cross a surface boundary.
        let p0 = on_surface(0.0);
        let q1 = on_surface(1.0 / 3.0);
        let q2 = on_surface(2.0 / 3.0);
        let p3 = on_surface(1.0);
        let interior = |a: f64, b: f64, c: f64, d: f64| {
            Point3::new(
                (a * p0.x + b * q1.x + c * q2.x + d * p3.x) / 6.0,
                (a * p0.y + b * q1.y + c * q2.y + d * p3.y) / 6.0,
                (a * p0.z + b * q1.z + c * q2.z + d * p3.z) / 6.0,
            )
        };
        let arc = [
            p0,
            interior(-5.0, 18.0, -9.0, 2.0),
            interior(2.0, -9.0, 18.0, -5.0),
            p3,
        ];
        let (first_half, second_half) = split_cubic_row(&arc);
        // The oriented segment shim reuses the calibrated march.
        let shim = OrientedSegment {
            curve: fillet.trim_curve_1.clone(),
            surface: fillet.surface.clone(),
            is_fillet_end: true,
            reversed: false,
        };
        let piece_a = shim.border_piece_calibrated(
            &first_half,
            &calibration[0],
            vertex_position,
            &mut first_order_marches,
        );
        let piece_b = shim.border_piece_calibrated(
            &second_half,
            &calibration[1],
            vertex_position,
            &mut first_order_marches,
        );
        sectors.push((order[k], piece_a, piece_b));
    }

    // Sector s spans [m_s → c_s → m_{s+1}]: piece A = second half of
    // arc s (m_s → c_s side? — the arc runs c_{s-1} → c_s, so its
    // halves are [c_{s-1} → m_s] and [m_s → c_s]). Regroup: sector s
    // = (arc_s second half, arc_{s+1} first half).
    type PiecePair = ([[Point3; 4]; 3], [[Point3; 4]; 3]);
    let pieces: Vec<PiecePair> = sectors.iter().map(|(_, a, b)| (*a, *b)).collect();
    let mut grouped = Vec::with_capacity(n);
    for s in 0..n {
        let piece_a = pieces[s].1; // [m_s → c_s]
        let piece_b = pieces[(s + 1) % n].0; // [c_s → m_{s+1}]
        grouped.push([
            Kp20PaperGuideBorderPiece { rows: piece_a },
            Kp20PaperGuideBorderPiece { rows: piece_b },
        ]);
    }

    // Corner face-tangency conditioning: at each corner c the two
    // G² sweeps are tangent to the shared support face, so the TRUE
    // ring data at c lies in that tangent plane and the cap must stay
    // on its outward side (the face wedge occupies the plane up to
    // c). Any below-plane component in the corner-adjacent cap-side
    // cross rows is extraction/mismatch noise that pushes the
    // assembled cap into the wedge (measured as wedge blades stabbing
    // the dome at larger radii); project it out.
    let loop_centroid = {
        let mut acc = [0.0_f64; 3];
        let mut count = 0.0;
        for sector in &grouped {
            for piece in sector {
                for cp in &piece.rows[0] {
                    acc[0] += cp.x;
                    acc[1] += cp.y;
                    acc[2] += cp.z;
                    count += 1.0;
                }
            }
        }
        Point3::new(acc[0] / count, acc[1] / count, acc[2] / count)
    };
    for s in 0..n {
        let c_point = grouped[s][0].rows[0][3];
        // Tangent plane at c from the owning sweep (= the shared
        // face plane, by G² tangency), oriented outward.
        let fillet = &fillet_surfaces[order[(s + 1) % n]];
        let Some((cu, cv)) = project_point_to_uv(&fillet.surface, &c_point) else {
            continue;
        };
        let mut normal = fillet.surface.normal_at(cu, cv);
        if normal.dot(&(c_point - loop_centroid)) < 0.0 {
            normal = -normal;
        }
        let clamp = |cp: &mut Point3| {
            let h = (*cp - c_point).dot(&normal);
            if h < 0.0 {
                *cp = *cp - normal * h;
            }
        };
        // Piece A of sector s ends at c (columns 3, 2); piece B of
        // sector s starts at c (columns 0, 1).
        for row in 1..3 {
            for col in [3_usize, 2] {
                clamp(&mut grouped[s][0].rows[row][col]);
            }
            for col in [0_usize, 1] {
                clamp(&mut grouped[s][1].rows[row][col]);
            }
        }
    }

    Ok((
        Kp20PaperGuideBorder { sectors: grouped },
        Kp20PaperBorderExtractReport {
            n,
            depth: 0.0,
            max_loop_gap: max_corner_gap,
            first_order_marches,
            reversed_for_orientation: false,
        },
    ))
}

/// De Casteljau halves of one cubic row.
fn split_cubic_row(row: &[Point3; 4]) -> ([Point3; 4], [Point3; 4]) {
    let mid =
        |a: Point3, b: Point3| Point3::new(0.5 * (a.x + b.x), 0.5 * (a.y + b.y), 0.5 * (a.z + b.z));
    let p01 = mid(row[0], row[1]);
    let p12 = mid(row[1], row[2]);
    let p23 = mid(row[2], row[3]);
    let p012 = mid(p01, p12);
    let p123 = mid(p12, p23);
    let center = mid(p012, p123);
    ([row[0], p01, p012, center], [center, p123, p23, row[3]])
}
