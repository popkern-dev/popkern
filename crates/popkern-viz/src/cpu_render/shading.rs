//! Lambertian shading with sRGB gamma correction.
//!
//! Matches the GPU mesh.wgsl shader: directional light at
//! `[0.577, 0.577, 0.577]` with 30% ambient. The GPU writes
//! linear values to `Rgba8UnormSrgb` (hardware gamma); the CPU
//! path applies sRGB conversion explicitly.

/// Normalized light direction (approximately `[1,1,1] / sqrt(3)`).
const LIGHT_DIR: [f32; 3] = [0.577, 0.577, 0.577];

/// Ambient lighting factor (matches GPU shader).
const AMBIENT: f32 = 0.3;

/// Convert a linear-space value to sRGB gamma (IEC 61966-2-1).
///
/// Values below the linear threshold use a simple scale factor;
/// values above use the standard power curve with offset.
pub fn linear_to_srgb(linear: f32) -> f32 {
    if linear <= 0.0031308 {
        linear * 12.92
    } else {
        1.055 * linear.powf(1.0 / 2.4) - 0.055
    }
}

/// Compute shaded color for a vertex.
///
/// Applies Lambertian diffuse + ambient lighting matching the GPU
/// `mesh.wgsl` fragment shader, then converts each RGB channel from
/// linear to sRGB gamma. Alpha is passed through unchanged.
pub fn shade_vertex(normal: [f32; 3], color: [f32; 4]) -> [f32; 4] {
    // Normalize the normal vector.
    let len = (normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]).sqrt();
    let n = if len > 1e-12 {
        [normal[0] / len, normal[1] / len, normal[2] / len]
    } else {
        [0.0, 0.0, 0.0]
    };

    // Dot product with light direction, clamped to [0, 1].
    let dot = (n[0] * LIGHT_DIR[0] + n[1] * LIGHT_DIR[1] + n[2] * LIGHT_DIR[2]).max(0.0);

    // Lighting = ambient + (1 - ambient) * diffuse.
    let lighting = AMBIENT + (1.0 - AMBIENT) * dot;

    // Apply lighting in linear space, then convert to sRGB.
    [
        linear_to_srgb((color[0] * lighting).clamp(0.0, 1.0)),
        linear_to_srgb((color[1] * lighting).clamp(0.0, 1.0)),
        linear_to_srgb((color[2] * lighting).clamp(0.0, 1.0)),
        color[3],
    ]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn front_facing_fully_lit() {
        // Normal pointing toward the light → high brightness.
        let result = shade_vertex([0.577, 0.577, 0.577], [1.0, 1.0, 1.0, 1.0]);
        // With full diffuse, lighting = 0.3 + 0.7 * 1.0 = 1.0,
        // so linear 1.0 → sRGB 1.0.
        for i in 0..3 {
            assert!(
                result[i] > 0.95,
                "channel {i} should be near 1.0, got {}",
                result[i]
            );
        }
        assert!((result[3] - 1.0).abs() < f32::EPSILON, "alpha unchanged");
    }

    #[test]
    fn back_facing_ambient_only() {
        // Normal pointing away from the light → ambient only.
        let result = shade_vertex([-0.577, -0.577, -0.577], [1.0, 1.0, 1.0, 1.0]);
        // lighting = 0.3 + 0.7 * 0.0 = 0.3
        // sRGB(0.3) ≈ 0.584
        let expected = linear_to_srgb(0.3);
        for i in 0..3 {
            assert!(
                (result[i] - expected).abs() < 0.01,
                "channel {i} should be ~{expected}, got {}",
                result[i]
            );
        }
    }

    #[test]
    fn srgb_gamma_applied() {
        // linear 0.5 → sRGB ≈ 0.735
        let srgb = linear_to_srgb(0.5);
        assert!(
            (srgb - 0.735).abs() < 0.01,
            "sRGB(0.5) should be ~0.735, got {srgb}"
        );
    }

    #[test]
    fn srgb_zero_and_one() {
        assert!((linear_to_srgb(0.0)).abs() < f32::EPSILON, "sRGB(0) = 0");
        assert!((linear_to_srgb(1.0) - 1.0).abs() < 0.001, "sRGB(1) = 1");
    }
}
