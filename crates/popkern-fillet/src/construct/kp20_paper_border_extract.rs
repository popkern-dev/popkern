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

    // Junction-harmonized transversal directions: the ring march
    // direction is a PARAMETERIZATION choice (the paper's ρ freedom).
    // Using each surface's natural direction leaves the two sides of
    // every junction with incompatible cross-rows on G¹ input; instead,
    // average the two sides' natural directions at every junction and
    // blend along each piece, projecting marches back onto the true
    // surfaces so the rows still sample real ring geometry.
    let mut junction_dirs: Vec<Vec3> = Vec::with_capacity(segment_count);
    for k in 0..segment_count {
        let previous = &curves[(k + segment_count - 1) % segment_count];
        let this = &curves[k];
        let p = this.start();
        let d_prev =
            previous.natural_march_direction(p, previous.tangent(1.0), boundary.vertex_position);
        let d_this = this.natural_march_direction(p, this.tangent(0.0), boundary.vertex_position);
        let sum = d_prev + d_this;
        junction_dirs.push(sum.normalized().unwrap_or(d_this));
    }

    let mut first_order_marches = 0_usize;
    let mut sectors = Vec::with_capacity(n);
    for k in 0..n {
        let piece_a = curves[2 * k].border_piece(
            depth,
            junction_dirs[2 * k],
            junction_dirs[(2 * k + 1) % segment_count],
            &mut first_order_marches,
        );
        let piece_b = curves[2 * k + 1].border_piece(
            depth,
            junction_dirs[(2 * k + 1) % segment_count],
            junction_dirs[(2 * k + 2) % segment_count],
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

    fn tangent(&self, t: f64) -> Vec3 {
        let raw = self.curve.tangent_at(self.parameter(t));
        if self.reversed { -raw } else { raw }
    }

    fn start(&self) -> Point3 {
        self.point(0.0)
    }

    fn end(&self) -> Point3 {
        self.point(1.0)
    }

    /// Unit in-surface transversal at `p` pointing away from the
    /// blend vertex (the surface's natural march direction).
    fn natural_march_direction(&self, p: Point3, boundary_tangent: Vec3, vertex: Point3) -> Vec3 {
        match &self.surface {
            SurfaceGeometry::Cylinder { axis, .. } => {
                let dir = if axis.dot(&(p - vertex)) >= 0.0 {
                    *axis
                } else {
                    -*axis
                };
                dir.normalized().unwrap_or(dir)
            }
            _ => {
                let normal =
                    surface_normal_near(&self.surface, p).unwrap_or(Vec3::new(0.0, 0.0, 1.0));
                in_surface_away_direction(p, boundary_tangent, normal, vertex)
            }
        }
    }

    /// Bi-3 depth-2 border piece: row 0 is a cubic Hermite fit of the
    /// segment; rows 1..2 are the cap-side C² prolongation of ring
    /// marches along junction-harmonized directions, projected onto
    /// the owning surface.
    fn border_piece(
        &self,
        depth: f64,
        start_dir: Vec3,
        end_dir: Vec3,
        first_order_marches: &mut usize,
    ) -> [[Point3; 4]; 3] {
        // Row 0: cubic Hermite from endpoints + endpoint tangents
        // (parameterization-honest: tangents scaled by the segment's
        // own parameter speed over [0, 1]).
        let p0 = self.point(0.0);
        let p3 = self.point(1.0);
        let m0 = self.tangent(0.0);
        let m1 = self.tangent(1.0);
        let (t0, t1) = self.curve.parameter_range();
        let speed = (t1 - t0).abs();
        let row0 = [p0, p0 + m0 * (speed / 3.0), p3 - m1 * (speed / 3.0), p3];

        // Ring marches at each row-0 control point along the blended
        // junction direction, projected onto the owning surface.
        let mut rows = [row0, row0, row0];
        for (index, &q0) in row0.iter().enumerate() {
            let blend = index as f64 / 3.0;
            let mixed = start_dir * (1.0 - blend) + end_dir * blend;
            let dir = mixed.normalized().unwrap_or(start_dir);
            let (q1, q2) = self.ring_march(q0, dir, depth, first_order_marches);
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

    /// March into the ring along `dir` on the owning surface:
    /// returns surface points near depth/3 and 2·depth/3 (exact
    /// in-plane; radially projected onto cylinders; tangent-plane
    /// first-order otherwise).
    fn ring_march(
        &self,
        q0: Point3,
        dir: Vec3,
        depth: f64,
        first_order_marches: &mut usize,
    ) -> (Point3, Point3) {
        let step = depth / 3.0;
        let raw1 = q0 + dir * step;
        let raw2 = q0 + dir * (2.0 * step);
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
            _ => {
                *first_order_marches += 1;
                (raw1, raw2)
            }
        }
    }
}

/// Unit direction in the tangent plane, perpendicular to the boundary
/// tangent, pointing away from the blend vertex.
fn in_surface_away_direction(
    q0: Point3,
    boundary_tangent: Vec3,
    normal: Vec3,
    vertex: Point3,
) -> Vec3 {
    let side = normal.cross(&boundary_tangent);
    let unit = side.normalized().unwrap_or(side);
    if unit.dot(&(q0 - vertex)) >= 0.0 {
        unit
    } else {
        -unit
    }
}

fn surface_normal_near(surface: &SurfaceGeometry, point: Point3) -> Option<Vec3> {
    use popkern_core::geometry::surface::SurfaceEval;
    use popkern_core::geometry::uv_projection::project_point_to_uv;
    let uv = project_point_to_uv(surface, &point)?;
    Some(surface.normal_at(uv.0, uv.1))
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
