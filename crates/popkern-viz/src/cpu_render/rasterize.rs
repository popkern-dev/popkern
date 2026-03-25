//! Scanline triangle rasterization with z-buffer and back-face culling.
//!
//! Provides a software framebuffer and a triangle rasterizer that matches
//! the GPU pipeline's winding-order convention (CCW = front-facing).

// ---------------------------------------------------------------------------
// Screen-space vertex
// ---------------------------------------------------------------------------

/// Screen-space vertex after projection.
///
/// Coordinates are in pixel space (fractional). Color is already-shaded
/// RGBA in `[u8; 4]` form (sRGB gamma applied).
#[derive(Clone, Debug)]
pub struct ScreenVertex {
    /// Pixel X coordinate (fractional).
    pub x: f32,
    /// Pixel Y coordinate (fractional).
    pub y: f32,
    /// Depth: 0 = near, 1 = far.
    pub z: f32,
    /// Already-shaded RGBA bytes.
    pub color: [u8; 4],
}

// ---------------------------------------------------------------------------
// Framebuffer
// ---------------------------------------------------------------------------

/// RGBA + depth framebuffer for CPU rendering.
///
/// Color is stored as per-pixel `[u8; 4]` (RGBA, sRGB gamma).
/// Depth is per-pixel `f32` in `[0, 1]` where 0 is near, 1 is far.
pub struct Framebuffer {
    /// Width in pixels.
    pub width: u32,
    /// Height in pixels.
    pub height: u32,
    /// Per-pixel RGBA color, row-major. Initialized to `[0, 0, 0, 0]`.
    pub color: Vec<[u8; 4]>,
    /// Per-pixel depth, row-major. Initialized to `1.0` (far plane).
    pub depth: Vec<f32>,
}

impl Framebuffer {
    /// Create a new framebuffer initialized to transparent black + far depth.
    pub fn new(width: u32, height: u32) -> Self {
        let len = (width as usize) * (height as usize);
        Self {
            width,
            height,
            color: vec![[0, 0, 0, 0]; len],
            depth: vec![1.0; len],
        }
    }

    /// Clear the framebuffer with a background color and reset depth.
    pub fn clear(&mut self, bg: [u8; 4]) {
        self.color.fill(bg);
        self.depth.fill(1.0);
    }

    /// Read the color of a pixel at `(x, y)`.
    ///
    /// Returns `[0, 0, 0, 0]` if coordinates are out of bounds.
    pub fn pixel(&self, x: u32, y: u32) -> [u8; 4] {
        if x >= self.width || y >= self.height {
            return [0, 0, 0, 0];
        }
        self.color[(y as usize) * (self.width as usize) + (x as usize)]
    }

    /// Flatten the framebuffer to a contiguous RGBA byte vector.
    ///
    /// Pixels are in row-major order: `[R, G, B, A, R, G, B, A, ...]`.
    pub fn to_rgba(&self) -> Vec<u8> {
        let mut out = Vec::with_capacity(self.color.len() * 4);
        for px in &self.color {
            out.extend_from_slice(px);
        }
        out
    }
}

// ---------------------------------------------------------------------------
// Triangle rasterization
// ---------------------------------------------------------------------------

/// Rasterize a single triangle into the framebuffer.
///
/// Uses the standard edge-function approach with back-face culling:
/// - CCW winding (positive signed area) = front-facing, drawn.
/// - CW or degenerate (signed area <= 0) = culled.
///
/// Depth test: smaller z wins (closer to camera).
pub fn rasterize_triangle(
    fb: &mut Framebuffer,
    v0: &ScreenVertex,
    v1: &ScreenVertex,
    v2: &ScreenVertex,
) {
    // Signed area * 2 (cross product of screen-space edges).
    let area = (v1.x - v0.x) * (v2.y - v0.y) - (v2.x - v0.x) * (v1.y - v0.y);
    if area <= 0.0 {
        // Back-facing or degenerate — cull.
        return;
    }
    let inv_area = 1.0 / area;

    // Bounding box, clamped to framebuffer.
    let min_x = v0.x.min(v1.x).min(v2.x).max(0.0).floor() as i32;
    let max_x = v0.x.max(v1.x).max(v2.x).min(fb.width as f32 - 1.0).ceil() as i32;
    let min_y = v0.y.min(v1.y).min(v2.y).max(0.0).floor() as i32;
    let max_y = v0.y.max(v1.y).max(v2.y).min(fb.height as f32 - 1.0).ceil() as i32;

    if min_x > max_x || min_y > max_y {
        return;
    }

    for py in min_y..=max_y {
        for px in min_x..=max_x {
            let px_f = px as f32 + 0.5;
            let py_f = py as f32 + 0.5;

            // Barycentric coordinates via edge functions.
            let w0 = (v1.x - v0.x) * (py_f - v0.y) - (v1.y - v0.y) * (px_f - v0.x);
            let w1 = (v2.x - v1.x) * (py_f - v1.y) - (v2.y - v1.y) * (px_f - v1.x);
            let w2 = (v0.x - v2.x) * (py_f - v2.y) - (v0.y - v2.y) * (px_f - v2.x);

            // Standard edge function: all three must be >= 0 for the
            // pixel center to be inside the CCW triangle.
            if w0 < 0.0 || w1 < 0.0 || w2 < 0.0 {
                continue;
            }

            // Normalize barycentric coords.
            let b0 = w1 * inv_area; // weight for v0
            let b1 = w2 * inv_area; // weight for v1
            let b2 = w0 * inv_area; // weight for v2

            // Interpolate depth.
            let z = b0 * v0.z + b1 * v1.z + b2 * v2.z;

            // Depth test (less = closer = wins).
            let idx = (py as usize) * (fb.width as usize) + (px as usize);
            if z >= fb.depth[idx] {
                continue;
            }

            // Interpolate color.
            let color = [
                (b0 * v0.color[0] as f32 + b1 * v1.color[0] as f32 + b2 * v2.color[0] as f32) as u8,
                (b0 * v0.color[1] as f32 + b1 * v1.color[1] as f32 + b2 * v2.color[1] as f32) as u8,
                (b0 * v0.color[2] as f32 + b1 * v1.color[2] as f32 + b2 * v2.color[2] as f32) as u8,
                (b0 * v0.color[3] as f32 + b1 * v1.color[3] as f32 + b2 * v2.color[3] as f32) as u8,
            ];

            fb.color[idx] = color;
            fb.depth[idx] = z;
        }
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn single_triangle_covers_pixels() {
        // A large CCW triangle covering the center of a 10x10 framebuffer.
        let mut fb = Framebuffer::new(10, 10);
        let green = [0, 255, 0, 255];
        let v0 = ScreenVertex {
            x: 1.0,
            y: 1.0,
            z: 0.5,
            color: green,
        };
        let v1 = ScreenVertex {
            x: 9.0,
            y: 1.0,
            z: 0.5,
            color: green,
        };
        let v2 = ScreenVertex {
            x: 5.0,
            y: 9.0,
            z: 0.5,
            color: green,
        };
        rasterize_triangle(&mut fb, &v0, &v1, &v2);

        // Center pixel (5, 5) should be inside the triangle.
        let px = fb.pixel(5, 5);
        assert_eq!(px, green, "center pixel should be green");

        // Corner pixel (0, 0) should still be transparent.
        let corner = fb.pixel(0, 0);
        assert_eq!(corner, [0, 0, 0, 0], "corner should be untouched");
    }

    #[test]
    fn z_buffer_closer_wins() {
        let mut fb = Framebuffer::new(10, 10);
        let red = [255, 0, 0, 255];
        let blue = [0, 0, 255, 255];

        // Far red triangle (z=0.8).
        let v0 = ScreenVertex {
            x: 0.0,
            y: 0.0,
            z: 0.8,
            color: red,
        };
        let v1 = ScreenVertex {
            x: 9.0,
            y: 0.0,
            z: 0.8,
            color: red,
        };
        let v2 = ScreenVertex {
            x: 5.0,
            y: 9.0,
            z: 0.8,
            color: red,
        };
        rasterize_triangle(&mut fb, &v0, &v1, &v2);

        // Close blue triangle (z=0.2) overlapping same area.
        let v0 = ScreenVertex {
            x: 0.0,
            y: 0.0,
            z: 0.2,
            color: blue,
        };
        let v1 = ScreenVertex {
            x: 9.0,
            y: 0.0,
            z: 0.2,
            color: blue,
        };
        let v2 = ScreenVertex {
            x: 5.0,
            y: 9.0,
            z: 0.2,
            color: blue,
        };
        rasterize_triangle(&mut fb, &v0, &v1, &v2);

        // Blue should win at center.
        let px = fb.pixel(5, 3);
        assert_eq!(px, blue, "closer blue triangle should win");
    }

    #[test]
    fn backface_culled() {
        // CW triangle (negative signed area) → should be culled entirely.
        let mut fb = Framebuffer::new(10, 10);
        let red = [255, 0, 0, 255];
        let v0 = ScreenVertex {
            x: 1.0,
            y: 1.0,
            z: 0.5,
            color: red,
        };
        let v1 = ScreenVertex {
            x: 5.0,
            y: 9.0,
            z: 0.5,
            color: red,
        };
        let v2 = ScreenVertex {
            x: 9.0,
            y: 1.0,
            z: 0.5,
            color: red,
        };
        rasterize_triangle(&mut fb, &v0, &v1, &v2);

        // All pixels should remain transparent.
        for y in 0..10 {
            for x in 0..10 {
                assert_eq!(
                    fb.pixel(x, y),
                    [0, 0, 0, 0],
                    "pixel ({x},{y}) should be untouched after backface cull"
                );
            }
        }
    }

    #[test]
    fn framebuffer_to_rgba() {
        let mut fb = Framebuffer::new(2, 2);
        fb.color[0] = [255, 0, 0, 255]; // top-left red
        fb.color[1] = [0, 255, 0, 255]; // top-right green
        fb.color[2] = [0, 0, 255, 255]; // bottom-left blue
        fb.color[3] = [255, 255, 0, 255]; // bottom-right yellow

        let bytes = fb.to_rgba();
        assert_eq!(bytes.len(), 16, "2x2 RGBA = 16 bytes");
        assert_eq!(&bytes[0..4], &[255, 0, 0, 255]);
        assert_eq!(&bytes[4..8], &[0, 255, 0, 255]);
        assert_eq!(&bytes[8..12], &[0, 0, 255, 255]);
        assert_eq!(&bytes[12..16], &[255, 255, 0, 255]);
    }
}
