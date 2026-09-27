//! Mesh-vs-surface normal disagreement of a zebra frame (PK-503):
//! the number that says how much of an image's stripe structure the
//! tessellation owns.

use popkern_core::topology::FaceHandle;

use super::{SurfaceHit, ZebraFrame, ZebraNormals, angle_between};

/// The worst pixel of a mode comparison.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ZebraWorstPixel {
    /// Pixel column.
    pub x: u32,
    /// Pixel row.
    pub y: u32,
    /// Owning face.
    pub face: FaceHandle,
    /// Native `(u, v)` on that face.
    pub uv: (f64, f64),
    /// Angle between the two normals (radians).
    pub angle_rad: f64,
}

/// Per-face share of a mode comparison.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ZebraFaceDisagreement {
    /// The face.
    pub face: FaceHandle,
    /// Evaluated pixels on the face.
    pub pixels: usize,
    /// Largest normal angle on the face (radians).
    pub max_angle_rad: f64,
    /// Mean normal angle on the face (radians).
    pub mean_angle_rad: f64,
}

/// How far the mesh-normal and surface-normal views of one frame
/// disagree.
#[derive(Debug, Clone, PartialEq)]
pub struct ZebraModeDisagreement {
    /// Pixels covered by the part.
    pub covered_pixels: usize,
    /// Covered pixels with an evaluated surface normal on the face the
    /// mesh shows.
    pub evaluated_pixels: usize,
    /// Evaluated pixels whose surface point is the ray hit (on the
    /// claiming face or, beside a seam, a neighbouring face), or the
    /// claiming face's domain point nearest the ray.
    pub ray_hit_pixels: usize,
    /// Covered pixels whose surface view shows a NEIGHBOURING face
    /// (between the mesh seam and the true seam). They are excluded
    /// from the angle statistics, which compare one face's normals.
    pub reassigned_pixels: usize,
    /// Evaluated pixels that fell back to the closest surface point.
    pub closest_point_pixels: usize,
    /// Evaluated pixels that fell back to the interpolated parameter.
    pub interpolated_pixels: usize,
    /// Largest angle between the two normals (radians).
    pub max_angle_rad: f64,
    /// Mean angle over evaluated pixels (radians).
    pub mean_angle_rad: f64,
    /// 99th-percentile angle over evaluated pixels (radians).
    pub p99_angle_rad: f64,
    /// Where the largest angle occurs.
    pub worst: Option<ZebraWorstPixel>,
    /// Evaluated pixels whose stripe band differs between the views
    /// for the compared stripe count.
    pub band_mismatch_pixels: usize,
    /// Per-face breakdown, largest `max_angle_rad` first.
    pub per_face: Vec<ZebraFaceDisagreement>,
}

impl ZebraFrame {
    /// Measure how far the mesh-normal and surface-normal views
    /// disagree, pixel by pixel, for `stripe_count` stripes.
    pub fn mode_disagreement(&self, stripe_count: u32) -> ZebraModeDisagreement {
        let mut covered = 0usize;
        let mut angles: Vec<f64> = Vec::new();
        let mut worst: Option<ZebraWorstPixel> = None;
        let mut band_mismatch = 0usize;
        let (mut ray_hits, mut closest, mut interpolated) = (0usize, 0usize, 0usize);
        let mut reassigned = 0usize;
        let mut per_face: Vec<(FaceHandle, usize, f64, f64)> = Vec::new();
        for (idx, sample) in self.samples.iter().enumerate() {
            let Some(s) = sample else { continue };
            covered += 1;
            if s.surface_face != s.face {
                // The views show different faces here (a pixel between
                // the mesh seam and the true seam): not a tessellation
                // error of one face's normals.
                reassigned += 1;
                continue;
            }
            let Some(sn) = s.surface_normal else { continue };
            let a = angle_between(s.mesh_normal, sn);
            angles.push(a);
            match s.hit {
                SurfaceHit::Ray | SurfaceHit::NeighbourRay | SurfaceHit::DomainBoundary => {
                    ray_hits += 1
                }
                SurfaceHit::ClosestPoint => closest += 1,
                SurfaceHit::Interpolated | SurfaceHit::NoSurface => interpolated += 1,
            }
            if self.band(s, ZebraNormals::Mesh, stripe_count)
                != self.band(s, ZebraNormals::Surface, stripe_count)
            {
                band_mismatch += 1;
            }
            if worst.is_none_or(|w| a > w.angle_rad) {
                worst = Some(ZebraWorstPixel {
                    x: (idx % self.width as usize) as u32,
                    y: (idx / self.width as usize) as u32,
                    face: s.face,
                    uv: s.uv,
                    angle_rad: a,
                });
            }
            match per_face.iter_mut().find(|e| e.0 == s.face) {
                Some(e) => {
                    e.1 += 1;
                    e.2 = e.2.max(a);
                    e.3 += a;
                }
                None => per_face.push((s.face, 1, a, a)),
            }
        }
        let evaluated = angles.len();
        let mean = if evaluated > 0 {
            angles.iter().sum::<f64>() / evaluated as f64
        } else {
            0.0
        };
        angles.sort_by(f64::total_cmp);
        let p99 = if evaluated > 0 {
            angles[((evaluated - 1) as f64 * 0.99).round() as usize]
        } else {
            0.0
        };
        let mut per_face: Vec<ZebraFaceDisagreement> = per_face
            .into_iter()
            .map(|(face, pixels, max, sum)| ZebraFaceDisagreement {
                face,
                pixels,
                max_angle_rad: max,
                mean_angle_rad: sum / pixels as f64,
            })
            .collect();
        per_face.sort_by(|a, b| b.max_angle_rad.total_cmp(&a.max_angle_rad));
        ZebraModeDisagreement {
            covered_pixels: covered,
            evaluated_pixels: evaluated,
            ray_hit_pixels: ray_hits,
            reassigned_pixels: reassigned,
            closest_point_pixels: closest,
            interpolated_pixels: interpolated,
            max_angle_rad: angles.last().copied().unwrap_or(0.0),
            mean_angle_rad: mean,
            p99_angle_rad: p99,
            worst,
            band_mismatch_pixels: band_mismatch,
            per_face,
        }
    }
}

/// One image-row crossing of a face seam in a zebra view: how the
/// reflection coordinate `s` behaves across it.
///
/// Along a row, `s(x)` is smooth inside a face. Across a seam, a G1
/// break shows as a VALUE jump (the stripes break), a G2 break as a
/// SLOPE jump (the stripes kink but stay connected). Both are
/// estimated from three pixels either side of the seam pixel pair:
/// `value_jump` is the gap between the two sides' quadratic
/// extrapolations to the seam midpoint, `slope_jump` the difference
/// of their extrapolated slopes there.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ZebraSeamCrossing {
    /// Face left of the seam.
    pub left: FaceHandle,
    /// Face right of the seam.
    pub right: FaceHandle,
    /// Last pixel column of the left face.
    pub x: u32,
    /// Pixel row.
    pub y: u32,
    /// Reflection-coordinate discontinuity at the seam.
    pub value_jump: f64,
    /// Reflection-coordinate slope discontinuity at the seam (per
    /// pixel).
    pub slope_jump: f64,
    /// The larger one-sided slope magnitude (per pixel): the scale a
    /// `slope_jump` is read against.
    pub slope_scale: f64,
    /// Angle (radians) between the two sides' normals extrapolated to
    /// the seam midpoint the same way: the G1 break in 3-D, blind to
    /// the stripe axis (the reflection coordinate sees one normal
    /// component only).
    pub normal_jump_rad: f64,
}

/// Columns either side within which a seam must reappear on the
/// neighbouring rows: a seam crossing rows at less than
/// `atan(1/2)` ≈ 27° runs too nearly along the row for a two-pixel
/// one-sided slope to straddle it honestly.
const SEAM_STEEPNESS_COLUMNS: usize = 2;

impl ZebraFrame {
    /// Whether the `left → right` transition at `(x, y)` also occurs
    /// on the rows above and below within [`SEAM_STEEPNESS_COLUMNS`].
    fn seam_crosses_rows(
        &self,
        x: usize,
        y: usize,
        left: FaceHandle,
        right: FaceHandle,
        normals: ZebraNormals,
    ) -> bool {
        let w = self.width as usize;
        let transition_on = |yy: usize| {
            let lo = x.saturating_sub(SEAM_STEEPNESS_COLUMNS);
            let hi = (x + SEAM_STEEPNESS_COLUMNS).min(w.saturating_sub(2));
            (lo..=hi).any(|xx| {
                let a = self.samples[yy * w + xx];
                let b = self.samples[yy * w + xx + 1];
                matches!((a, b), (Some(a), Some(b))
                    if a.face_in(normals) == left && b.face_in(normals) == right)
            })
        };
        y > 0 && y + 1 < self.height as usize && transition_on(y - 1) && transition_on(y + 1)
    }

    /// Every row crossing of a face seam with two covered pixels on
    /// each side, and steep enough to cross the neighbouring rows too,
    /// measured in the given view.
    pub fn seam_crossings(&self, normals: ZebraNormals) -> Vec<ZebraSeamCrossing> {
        let w = self.width as usize;
        let mut out = Vec::new();
        for y in 0..self.height as usize {
            let row = &self.samples[y * w..(y + 1) * w];
            for x in 2..w.saturating_sub(3) {
                let (Some(a2), Some(a1), Some(a0), Some(b0), Some(b1), Some(b2)) = (
                    row[x - 2],
                    row[x - 1],
                    row[x],
                    row[x + 1],
                    row[x + 2],
                    row[x + 3],
                ) else {
                    continue;
                };
                let f = |k: &super::ZebraSample| k.face_in(normals);
                let (fa, fb) = (f(&a0), f(&b0));
                if fa == fb
                    || f(&a2) != fa
                    || f(&a1) != fa
                    || f(&b1) != fb
                    || f(&b2) != fb
                    || !self.seam_crosses_rows(x, y, fa, fb, normals)
                {
                    continue;
                }
                let s = |k: &super::ZebraSample| self.reflection_of(k, normals);
                let (sa2, sa1, sa0) = (s(&a2), s(&a1), s(&a0));
                let (sb0, sb1, sb2) = (s(&b0), s(&b1), s(&b2));
                // Quadratic (Newton-form) extrapolation of each side's
                // three pixels to the seam midpoint, half a pixel past
                // its nearest sample: value and slope there are exact
                // up to third-order terms, so a smooth side whose
                // curvature varies fast (a G2 round's rim) reads no
                // spurious slope jump.
                let (d1, d2) = (sa0 - sa1, sa0 - 2.0 * sa1 + sa2);
                let (e1, e2) = (sb1 - sb0, sb2 - 2.0 * sb1 + sb0);
                let left_slope = d1 + d2;
                let right_slope = e1 - e2;
                let left_at_mid = sa0 + 0.5 * d1 + 0.375 * d2;
                let right_at_mid = sb0 - 0.5 * e1 + 0.375 * e2;
                let n = |k: &super::ZebraSample| k.normal(normals);
                let (na2, na1, na0) = (n(&a2), n(&a1), n(&a0));
                let (nb0, nb1, nb2) = (n(&b0), n(&b1), n(&b2));
                let left_n: [f64; 3] = std::array::from_fn(|i| {
                    let (d1, d2) = (na0[i] - na1[i], na0[i] - 2.0 * na1[i] + na2[i]);
                    na0[i] + 0.5 * d1 + 0.375 * d2
                });
                let right_n: [f64; 3] = std::array::from_fn(|i| {
                    let (e1, e2) = (nb1[i] - nb0[i], nb2[i] - 2.0 * nb1[i] + nb0[i]);
                    nb0[i] - 0.5 * e1 + 0.375 * e2
                });
                out.push(ZebraSeamCrossing {
                    left: fa,
                    right: fb,
                    x: x as u32,
                    y: y as u32,
                    value_jump: (right_at_mid - left_at_mid).abs(),
                    slope_jump: (right_slope - left_slope).abs(),
                    slope_scale: left_slope.abs().max(right_slope.abs()),
                    normal_jump_rad: angle_between(left_n, right_n),
                });
            }
        }
        out
    }
}
