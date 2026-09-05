//! CPU wireframe rendering via Bresenham line drawing.
//!
//! Draws edge overlays into an existing [`Framebuffer`] using integer
//! Bresenham's line algorithm. Pixels outside the framebuffer bounds
//! are silently clipped.

use crate::cpu_render::rasterize::Framebuffer;

/// Wireframe line color: black (visible against dark gray background).
pub const WIREFRAME_COLOR: [u8; 4] = [0, 0, 0, 255];

/// Draw a line between two screen-space points using Bresenham's algorithm.
///
/// Coordinates are in pixel space (fractional, rounded to nearest integer).
/// Lines that extend outside the framebuffer are clipped per-pixel — no
/// panic on out-of-bounds coordinates.
pub fn draw_line(fb: &mut Framebuffer, x0: f32, y0: f32, x1: f32, y1: f32, color: [u8; 4]) {
    draw_line_depth(
        fb,
        x0,
        y0,
        f32::NEG_INFINITY,
        x1,
        y1,
        f32::NEG_INFINITY,
        0.0,
        color,
    );
}

/// [`draw_line`] with a per-endpoint depth and a z-buffer test, so a
/// wireframe overlay only marks edges the camera can see. Undepthed
/// callers pass `-inf` and always win. `bias` pulls the interpolated
/// depth toward the camera: an edge lies exactly ON its faces'
/// surface, and with zero bias the faces' own rasterized depth
/// (offset by facet sagitta) stitches through the line. The caller
/// scales the bias to the scene's rasterized depth range — a constant
/// mis-scales one way or the other (measured on a propeller: hidden
/// blade outlines drawn straight through the hub at 2e-3, visible
/// silhouette rims eaten at -5e-3).
#[allow(clippy::too_many_arguments)]
pub fn draw_line_depth(
    fb: &mut Framebuffer,
    x0: f32,
    y0: f32,
    z0: f32,
    x1: f32,
    y1: f32,
    z1: f32,
    bias: f32,
    color: [u8; 4],
) {
    let mut ix0 = x0.round() as i32;
    let mut iy0 = y0.round() as i32;
    let ix1 = x1.round() as i32;
    let iy1 = y1.round() as i32;

    let dx = (ix1 - ix0).abs();
    let dy = -(iy1 - iy0).abs();
    let sx: i32 = if ix0 < ix1 { 1 } else { -1 };
    let sy: i32 = if iy0 < iy1 { 1 } else { -1 };
    let mut err = dx + dy;
    let steps = dx.max(-dy).max(1) as f32;
    let mut step = 0.0f32;

    loop {
        // Write pixel if within bounds and in front of the z-buffer.
        if ix0 >= 0 && iy0 >= 0 && (ix0 as u32) < fb.width && (iy0 as u32) < fb.height {
            let idx = (iy0 as usize) * (fb.width as usize) + (ix0 as usize);
            // Undepthed endpoints are -inf; (z1 - z0) would be NaN
            // there, so interpolate only over a finite span.
            let dz = z1 - z0;
            let z = if dz.is_finite() {
                z0 + dz * (step / steps) - bias
            } else {
                z0.min(z1)
            };
            if z <= fb.depth[idx] {
                fb.color[idx] = color;
            }
        }

        if ix0 == ix1 && iy0 == iy1 {
            break;
        }
        step += 1.0;

        let e2 = 2 * err;
        if e2 >= dy {
            err += dy;
            ix0 += sx;
        }
        if e2 <= dx {
            err += dx;
            iy0 += sy;
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
    fn depth_tested_line_hides_behind_nearer_surface() {
        let mut fb = Framebuffer::new(10, 10);
        // A "face" already rasterized at depth 0.6 over row 5.
        for x in 0..10 {
            let idx = 5 * 10 + x;
            fb.depth[idx] = 0.6;
        }
        let white = [255, 255, 255, 255];
        // A line BEHIND the face (0.9) stays hidden.
        draw_line_depth(&mut fb, 1.0, 5.0, 0.9, 8.0, 5.0, 0.9, 1.0e-3, white);
        for x in 1..=8 {
            assert_eq!(fb.pixel(x, 5), [0, 0, 0, 0], "hidden at ({x},5)");
        }
        // A line IN FRONT (0.3) draws.
        draw_line_depth(&mut fb, 1.0, 5.0, 0.3, 8.0, 5.0, 0.3, 1.0e-3, white);
        for x in 1..=8 {
            assert_eq!(fb.pixel(x, 5), white, "visible at ({x},5)");
        }
    }

    #[test]
    fn horizontal_line() {
        let mut fb = Framebuffer::new(10, 10);
        let white = [255, 255, 255, 255];
        draw_line(&mut fb, 1.0, 5.0, 8.0, 5.0, white);

        // All pixels from x=1..=8 at y=5 should be white.
        for x in 1..=8 {
            assert_eq!(fb.pixel(x, 5), white, "pixel ({x}, 5) should be white");
        }
        // Pixel outside the line should be untouched.
        assert_eq!(fb.pixel(0, 5), [0, 0, 0, 0]);
        assert_eq!(fb.pixel(9, 5), [0, 0, 0, 0]);
    }

    #[test]
    fn vertical_line() {
        let mut fb = Framebuffer::new(10, 10);
        let white = [255, 255, 255, 255];
        draw_line(&mut fb, 5.0, 1.0, 5.0, 8.0, white);

        // All pixels from y=1..=8 at x=5 should be white.
        for y in 1..=8 {
            assert_eq!(fb.pixel(5, y), white, "pixel (5, {y}) should be white");
        }
        // Pixel outside the line should be untouched.
        assert_eq!(fb.pixel(5, 0), [0, 0, 0, 0]);
        assert_eq!(fb.pixel(5, 9), [0, 0, 0, 0]);
    }

    #[test]
    fn line_clipped_to_bounds() {
        // Line extends well outside the 10x10 framebuffer — must not panic.
        let mut fb = Framebuffer::new(10, 10);
        let red = [255, 0, 0, 255];
        draw_line(&mut fb, -5.0, 5.0, 15.0, 5.0, red);

        // Visible portion should be drawn.
        assert_eq!(fb.pixel(0, 5), red);
        assert_eq!(fb.pixel(9, 5), red);
    }
}
