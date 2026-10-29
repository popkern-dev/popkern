//! The looks: per-pixel shading of a fragment in its part's colour.
//!
//! Lights follow popkern live's viewer: a hemisphere light (white sky
//! over a slate ground, about +Z) and a key light from `(3, −4, 6)`,
//! fixed in the world, so a turntable turns the part under the light.
//! `studio`/`cad` add a soft specular; `toon` quantises the key light's
//! diffuse term to the page's three-tone ramp (90 / 175 / 255).

use glam::Vec3;
use popkern_core::geometry::PortableMath;

/// A film look (`storyboard(#{ look })`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Look {
    /// The CAD shading on a plain background.
    Studio,
    /// Flat three-tone shading, black outlines.
    Toon,
    /// The viewer's CAD look.
    Cad,
}

impl Look {
    /// The look named `name` (`"toon"`, `"studio"`, `"cad"`; anything
    /// else is the studio look, as the page treats it).
    pub fn named(name: &str) -> Look {
        match name {
            "toon" => Look::Toon,
            "cad" => Look::Cad,
            _ => Look::Studio,
        }
    }
}

/// sRGB → linear.
pub fn srgb_to_linear(c: f32) -> f32 {
    if c <= 0.04045 {
        c / 12.92
    } else {
        ((c + 0.055) / 1.055).ppowf(2.4)
    }
}

/// Linear → sRGB.
pub fn linear_to_srgb(c: f32) -> f32 {
    let c = c.clamp(0.0, 1.0);
    if c <= 0.0031308 {
        c * 12.92
    } else {
        1.055 * c.ppowf(1.0 / 2.4) - 0.055
    }
}

/// An sRGB triple in linear light.
pub fn linear(rgb: [f32; 3]) -> Vec3 {
    Vec3::new(
        srgb_to_linear(rgb[0]),
        srgb_to_linear(rgb[1]),
        srgb_to_linear(rgb[2]),
    )
}

/// The key light's direction (towards the light).
pub fn key_direction() -> Vec3 {
    Vec3::new(3.0, -4.0, 6.0).normalize()
}

/// The hover blue (`HOVER` 0x3b82f6), linear.
fn highlight_color() -> Vec3 {
    linear([
        0x3b as f32 / 255.0,
        0x82 as f32 / 255.0,
        0xf6 as f32 / 255.0,
    ])
}

/// The toon ramp: three tones, nearest (`toonRamp`).
fn toon_tone(x: f32) -> f32 {
    let i = ((x.clamp(0.0, 1.0) * 3.0) as usize).min(2);
    [90.0, 175.0, 255.0][i] / 255.0
}

/// A lit fragment: `albedo` (linear) with unit normal `n` facing the
/// viewer, `view` the unit direction from the fragment to the eye.
pub fn lit(look: Look, albedo: Vec3, n: Vec3, view: Vec3, highlight: f32) -> Vec3 {
    let sky = Vec3::ONE;
    let ground = linear([
        0x33 as f32 / 255.0,
        0x44 as f32 / 255.0,
        0x55 as f32 / 255.0,
    ]);
    let hemi = ground.lerp(sky, 0.5 * n.z + 0.5);
    let l = key_direction();
    let ndl = n.dot(l);
    let color = match look {
        Look::Toon => {
            let tone = toon_tone(0.5 * ndl + 0.5);
            albedo * (0.25 * hemi + 0.8 * tone)
        }
        Look::Studio | Look::Cad => {
            let diffuse = albedo * (0.5 * hemi + 0.7 * ndl.max(0.0));
            let h = (l + view).normalize_or_zero();
            let spec = 0.12 * n.dot(h).max(0.0).ppowf(48.0);
            diffuse + Vec3::splat(spec)
        }
    };
    color + highlight_color() * (0.45 * highlight)
}

/// The flat tone of a section cut seen through the plane (a back face
/// behind the cut): the part's colour, dimmed.
pub fn cap(albedo: Vec3) -> Vec3 {
    albedo * 0.55
}

/// The colour edges are drawn in (`CAD_EDGE`).
pub const EDGE: Vec3 = Vec3::ZERO;
