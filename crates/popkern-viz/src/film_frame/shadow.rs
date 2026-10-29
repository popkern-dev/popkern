//! The contact shadow (`FilmLook.placeGround`): a key light above and a
//! little in front of the model (`c + (0.6r, −0.9r, 4r)`) casts the
//! model onto a ground plane just under its bounds; the ground shows
//! only the shadow, darkening the background by up to 16 %.
//!
//! The light is orthographic over ±2.5 r (as the page's shadow camera),
//! so the ground's darkness is a texture in the light's coordinates: a
//! depth map of every drawn triangle, compared with the ground plane's
//! own depth (occluded = 1), then blurred for the soft edge the page's
//! PCF gives. A background pixel looks its ground point up bilinearly:
//! one lookup instead of a filter per pixel.

use glam::{DVec3, Vec3};

use super::{FrameScene, Tri};

/// Texels per side of the darkness texture.
const SIZE: usize = 512;
/// Blur radius in texels (the page's `shadow.radius` is 6 at 2048).
const BLUR: usize = 4;

/// The ground plane's shadow, ready to look up.
pub(crate) struct ShadowMap {
    /// Height of the ground plane.
    pub ground_z: f32,
    centre: Vec3,
    right: Vec3,
    up: Vec3,
    half: f32,
    darkness: Vec<f32>,
}

impl ShadowMap {
    /// The shadow of the drawn triangles on the ground under `bounds`.
    pub(crate) fn build(
        scene: &FrameScene<'_>,
        opaque: &[Tri],
        translucent: &[Tri],
        bounds: [DVec3; 2],
    ) -> ShadowMap {
        let [lo, hi] = bounds;
        let c = ((lo + hi) * 0.5).as_vec3();
        let r = (((hi - lo).length() * 0.5).max(1e-9)) as f32;
        let light = c + Vec3::new(0.6 * r, -0.9 * r, 4.0 * r);
        let dir = (c - light).normalize();
        let right = dir.cross(Vec3::Y).normalize_or(Vec3::X);
        let up = right.cross(dir);
        let half = 2.5 * r;
        let ground_z = (lo.z as f32) - 1e-3 * r;
        let texel = |p: Vec3| {
            let d = p - c;
            let u = (d.dot(right) / half * 0.5 + 0.5) * SIZE as f32;
            let v = (d.dot(up) / half * 0.5 + 0.5) * SIZE as f32;
            (u, v, (p - light).dot(dir))
        };
        let mut depth = vec![f32::INFINITY; SIZE * SIZE];
        let _ = scene;
        for tri in opaque.iter().chain(translucent) {
            let s = tri.world.map(texel);
            let area =
                (s[1].0 - s[0].0) * (s[2].1 - s[0].1) - (s[2].0 - s[0].0) * (s[1].1 - s[0].1);
            if !area.is_finite() || area.abs() < 1e-12 {
                continue;
            }
            let min_x = s
                .iter()
                .map(|p| p.0)
                .fold(f32::INFINITY, f32::min)
                .floor()
                .max(0.0) as usize;
            let max_x = s
                .iter()
                .map(|p| p.0)
                .fold(f32::NEG_INFINITY, f32::max)
                .ceil()
                .min(SIZE as f32 - 1.0);
            let min_y = s
                .iter()
                .map(|p| p.1)
                .fold(f32::INFINITY, f32::min)
                .floor()
                .max(0.0) as usize;
            let max_y = s
                .iter()
                .map(|p| p.1)
                .fold(f32::NEG_INFINITY, f32::max)
                .ceil()
                .min(SIZE as f32 - 1.0);
            if max_x < 0.0 || max_y < 0.0 {
                continue;
            }
            for y in min_y..=max_y as usize {
                for x in min_x..=max_x as usize {
                    let (px, py) = (x as f32 + 0.5, y as f32 + 0.5);
                    let w0 = ((s[1].0 - px) * (s[2].1 - py) - (s[2].0 - px) * (s[1].1 - py)) / area;
                    let w1 = ((s[2].0 - px) * (s[0].1 - py) - (s[0].0 - px) * (s[2].1 - py)) / area;
                    let w2 = 1.0 - w0 - w1;
                    if w0 < 0.0 || w1 < 0.0 || w2 < 0.0 {
                        continue;
                    }
                    let d = w0 * s[0].2 + w1 * s[1].2 + w2 * s[2].2;
                    let i = y * SIZE + x;
                    depth[i] = depth[i].min(d);
                }
            }
        }
        // Occluded where something lies nearer the light than the ground.
        let bias = 0.01 * r;
        let mut hard = vec![0.0f32; SIZE * SIZE];
        for y in 0..SIZE {
            for x in 0..SIZE {
                let u = ((x as f32 + 0.5) / SIZE as f32 * 2.0 - 1.0) * half;
                let v = ((y as f32 + 0.5) / SIZE as f32 * 2.0 - 1.0) * half;
                // The ground point under this texel along the light.
                let p = c + right * u + up * v;
                if dir.z.abs() < 1e-9 {
                    continue;
                }
                let s = (ground_z - p.z) / dir.z;
                let g = p + dir * s;
                let ground_depth = (g - light).dot(dir);
                if depth[y * SIZE + x] < ground_depth - bias {
                    hard[y * SIZE + x] = 1.0;
                }
            }
        }
        ShadowMap {
            ground_z,
            centre: c,
            right,
            up,
            half,
            darkness: blur(&blur(&hard, true), false),
        }
    }

    /// How dark the ground is at ground point `g` (0..1).
    pub(crate) fn darkness(&self, g: Vec3) -> f32 {
        let d = g - self.centre;
        // The texel whose light ray lands on `g`: project along the
        // light onto the texture plane through the centre.
        let dir = self.right.cross(self.up);
        let back = d - dir * d.dot(dir);
        let u = (back.dot(self.right) / self.half * 0.5 + 0.5) * SIZE as f32 - 0.5;
        let v = (back.dot(self.up) / self.half * 0.5 + 0.5) * SIZE as f32 - 0.5;
        if !(0.0..(SIZE - 1) as f32).contains(&u) || !(0.0..(SIZE - 1) as f32).contains(&v) {
            return 0.0;
        }
        let (x, y) = (u.floor() as usize, v.floor() as usize);
        let (fx, fy) = (u - x as f32, v - y as f32);
        let at = |x: usize, y: usize| self.darkness[y * SIZE + x];
        let top = at(x, y) * (1.0 - fx) + at(x + 1, y) * fx;
        let bottom = at(x, y + 1) * (1.0 - fx) + at(x + 1, y + 1) * fx;
        top * (1.0 - fy) + bottom * fy
    }
}

/// One pass of a box blur, along x or y.
fn blur(src: &[f32], along_x: bool) -> Vec<f32> {
    let mut out = vec![0.0; src.len()];
    let r = BLUR as isize;
    for y in 0..SIZE as isize {
        for x in 0..SIZE as isize {
            let mut sum = 0.0;
            for k in -r..=r {
                let (sx, sy) = if along_x { (x + k, y) } else { (x, y + k) };
                if (0..SIZE as isize).contains(&sx) && (0..SIZE as isize).contains(&sy) {
                    sum += src[sy as usize * SIZE + sx as usize];
                }
            }
            out[y as usize * SIZE + x as usize] = sum / (2 * r + 1) as f32;
        }
    }
    out
}
