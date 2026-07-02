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

    // Pass 1: rim rows for every segment.
    let mut rims: Vec<[Point3; 4]> = curves.iter().map(|segment| segment.rim_row()).collect();

    // Boundary conditioning to the g∘χ̃ shape class. The
    // characteristic layout corners ONLY at the n sector corners c
    // and is SMOOTH at the n spoke points m, while a raw fillet apex
    // boundary corners at all 2n junctions. Junction j sits at the
    // start of segment j; with the FilletEnd-first phase, EVEN
    // junctions are the spoke points m (bridge → arc) and ODD
    // junctions are the sector corners c (arc → bridge). The real
    // corner angle at c already matches the characteristic class
    // (~119° vs ~115° at n = 5), so c is kept as-is; at every m the
    // two adjacent rim legs are rotated onto their common unified
    // tangent (endpoints fixed, inner control legs only — a re-trim
    // within the ring surfaces' tolerance band), making the rim G¹
    // at m exactly.
    let mut m_tangents: Vec<Option<Vec3>> = vec![None; segment_count];
    for j in (0..segment_count).step_by(2) {
        let p_index = (j + segment_count - 1) % segment_count;
        let q_index = j;
        let leg_in = rims[p_index][3] - rims[p_index][2];
        let leg_out = rims[q_index][1] - rims[q_index][0];
        let unified = (leg_in.normalized().unwrap_or(leg_in)
            + leg_out.normalized().unwrap_or(leg_out))
        .normalized();
        let Some(t_hat) = unified else { continue };
        // Unified length as well: the χ̃-calibrated cross magnitudes
        // derive from the rim legs, so cross-sector consistency at m
        // needs equal legs on both sides.
        let len = 0.5 * (leg_in.length() + leg_out.length());
        rims[p_index][2] = rims[p_index][3] - t_hat * len;
        rims[q_index][1] = rims[q_index][0] + t_hat * len;
        m_tangents[j] = Some(t_hat);
    }
    // Second pass at m: unify rim curvature (C² across the spoke, the
    // eq (13)/(14) prerequisite) by averaging the two sides' second
    // BB differences and pushing the correction into the third
    // control point from each side.
    for j in (0..segment_count).step_by(2) {
        let p_index = (j + segment_count - 1) % segment_count;
        let q_index = j;
        let pt = |p: Point3| [p.x, p.y, p.z];
        let p3 = pt(rims[p_index][3]);
        let p2 = pt(rims[p_index][2]);
        let p1 = pt(rims[p_index][1]);
        let q0 = pt(rims[q_index][0]);
        let q1 = pt(rims[q_index][1]);
        let q2 = pt(rims[q_index][2]);
        let mut new_p1 = [0.0_f64; 3];
        let mut new_q2 = [0.0_f64; 3];
        for axis in 0..3 {
            let d2_in = p1[axis] - 2.0 * p2[axis] + p3[axis];
            let d2_out = q2[axis] - 2.0 * q1[axis] + q0[axis];
            let mean = 0.5 * (d2_in + d2_out);
            new_p1[axis] = 2.0 * p2[axis] - p3[axis] + mean;
            new_q2[axis] = 2.0 * q1[axis] - q0[axis] + mean;
        }
        rims[p_index][1] = Point3::new(new_p1[0], new_p1[1], new_p1[2]);
        rims[q_index][2] = Point3::new(new_q2[0], new_q2[1], new_q2[2]);
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
